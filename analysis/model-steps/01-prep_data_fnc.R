
prep_data <- function(param_list) {
  library(tidyverse)
  library(lubridate)
  library(rstan)
  library(foreach)
  library(doParallel)
  library(covidmodeldata)
  
 # define params in environment --------------------------------------------
  DATE <- param_list$DATE
  NYT_FILE <- param_list$NYT_FILE
  ACS_FILE <- param_list$ACS_FILE
  DATE_0 <- param_list$DATE_0
  DATE_N <- param_list$DATE_N
  SAMPLES_ROOT <- param_list$SAMPLES_ROOT
  NWORKERS <- param_list$NWORKERS
  NNODES <- param_list$NNODES
  RESULTS_FILE <- param_list$RESULTS_FILE
  UTILS_FUNC_FILE <- param_list$UTILS_FUNC_FILE
  CLEAN_DIR <- param_list$CLEAN_DIR
  SSH_KEY <- param_list$SSH_KEY
  cades_workers <- param_list$cades_workers
  NITER <- param_list$NITER
  NTHIN <- param_list$NTHIN
  NCHAINS <- param_list$NCHAINS
  WARMUP <- param_list$WARMUP
  DATA_DIR <- param_list$DATA_DIR
  DIAG_DF_LOC <- param_list$DIAG_DF_LOC
  ZERO_PAD <- param_list$ZERO_PAD
  TPRED <- param_list$TPRED
  SPL_K <- param_list$SPL_K
  
  source(UTILS_FUNC_FILE)
  
  if(CLEAN_DIR) {
    unlink(DATA_DIR, recursive=TRUE)
    print(glue::glue("{DATA_DIR} removed before rewriting"))
  }
  
  if (is.null(NYT_FILE)) {
    nyt_data <- get_nyt() %>%
      format_nyt(distribute_unknowns = FALSE) # Assign KC, but that's it.
    
    # on 7/14 data for Northern Mariana Islands started being reported.
    # we dont have population numbers aligned for those currently.
    # so for now remove any geoids that are not in our population data
    nyt_data <- nyt_data %>%
      filter(
        geoid %in% covidmodeldata::landscan_usa$geoid
      )
    
  } else nyt_data <- read_csv(NYT_FILE)
  
  nyt_data <- mutate(nyt_data, date = lubridate::as_date(date))
  
  # Create tibble with cols
  #  geoid
  #  date
  #  total_cases
  #  Y
  #  t
  covid_df <- proc_covid(
    raw_data = nyt_data %>% select(geoid, date, total_cases=total_cases_mdl),
    nyt_data = nyt_data,
    DATE_0 = DATE_0,
    ZERO_PAD = ZERO_PAD
  ) 
 
  covid_df <- covid_df %>%
    left_join(nyt_data %>% select(geoid, date, new_cases)) %>%
    filter(date <= min(as.Date(DATE_N), as.Date(DATE)))

  if (is.null(ACS_FILE)) {
    acs_data <- acs_data # from the covidmodeldata package 
  } else {acs_data <- readr::read_rds(ACS_FILE)}
  
  # Create tibble with
  # geoid, mystate, i, group1, group2, pop
  county_df <- proc_county(
    raw_data = acs_data,
    geoid.list = unique(covid_df$geoid)
  )
  
  covid_df <- remove_prisons(covid_df, county_df) %>%
    filter(!is.na(Y))
  
  date_df <- tibble(
    date = seq.Date(from=as.Date(DATE_0),
                    t = max(covid_df$date)+TPRED,
                    by=1)
  ) %>% 
    mutate(t = as.integer(date-as.Date(DATE_0)+1))
  
  
  ###################################################
  # Create the time splines:
  X_dow <- date_df %>%
    mutate(dow = factor(weekdays(date), 
                        levels= c("Monday",  "Tuesday", "Wednesday", 
                                  "Thursday", "Friday", "Saturday","Sunday" ))) %>%
    arrange(t) %>%
    model.matrix(~dow-1, data=.)
  
  ######################################
  # Create spline basis for the time series
  # With one knot per date, the spline is B * alpha = I * alpha
  # The coefficients alpha have a random walk
  # 1. Set up a random walk smoothing matrix S=DD'
  # i.e. alpha ~ N(0,DD')
  # 2. Reparameterize so that Ba = Za', a' ~ N(0,I)
  
  ####################################
  # RE method 1
  # Create the intrinsic RW penalty on coef of a Identity matrix spline
  
  # KEY OUTPUTS OF THIS SECTION
  #  Z, the spline
  #  krig_wt,
  #  krig_var_chol
  T <- as.integer(max(covid_df$date)-as_date(DATE_0)+1)
  
  L <- pracma::tril(toeplitz(c(1,-2,1, rep(0, T+TPRED-4 ))))
  
  # Create a spline over observable time
  svd <- svd(L[1:(T-1),1:(T-1)])
  Z <- svd$v %*% diag(1/svd$d)
  
  # Flip order so most important is first:
  Z <- Z[,rev(1:ncol(Z))]
  
  ZZ <- Z
  # Now, take a low-rank-approximation
  Z <- Z[,1:SPL_K]
  # Now add the first time to it
  Z <- rbind(0, Z)
  
  # Now calculate the kriging weight and kriging variance
  Cov <- tcrossprod(solve(L))
  C11 <- Cov[1:(T-1),1:(T-1)]
  C22 <- Cov[T:(T+TPRED-1),T:(T+TPRED-1)]
  C12 <- Cov[1:(T-1),T:(T+TPRED-1) ]
  krig_wt <- zapsmall(solve(C11, C12))
  krig_var_chol <- t(zapsmall(chol(C22 - t(C12) %*% krig_wt)))
  krig_wt <- rbind(0, krig_wt)
  
  ## PRIOR on variance os spline
  ## Expect the value to grow by 3 orders of magnitude, over 90 days
  ## The variance seems to follow a quadratic equation with time, but I
  ## don't know what the parameters are.
  ## Trial and error suggests .02 is a good prior for the variance.
  
  # variance on last day is:
  max_var <- Cov[T,T]
  # sd that should give us log(10^4) change  
  sd_scale <- log(10^2) / sqrt(max_var)
  
  
  if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive=TRUE)
  save(covid_df, county_df, X_dow, Z, date_df, 
       file=file.path(DATA_DIR, 'data_frames.Rdata'))
  
  
  # We need a list of states to loop through.
  # NO LONGER NECESSARY Check the coviddf to make sure that there are enough cases statewide to model.
  # Except for DC (11), which we'll process with Maryland (24)
  
  # Step 1. get a list of states from the NYT file, 
  # and order from most counties to least.
  state_list <- county_df %>%
    group_by(mystate) %>%
    distinct(county_name) %>%
    summarize(n=n())  %>%
    arrange(desc(n)) %>%
    pull(mystate)
  
  #####################################################
  # CREATE A LIST WITH ONE ENTRY PER CHAIN (50 * NITER)
  # 
  stan_fit_list <- vector(mode='list', length=length(state_list)*NCHAINS)
  
  slot = 1 # A counter
  for(i in 1:length(state_list)){
    
    STATE_FIPS = state_list[i]
    
    # Create sample directory
    STATE_SAMPLES_DIR <- file.path(SAMPLES_ROOT, DATE, STATE_FIPS)
    if(!dir.exists(STATE_SAMPLES_DIR)) dir.create(STATE_SAMPLES_DIR, recursive=TRUE)
    
    
    ####################
    # Pull out geo data for this state
    this_county_df <- county_df %>%
      filter(mystate == STATE_FIPS) %>%
      arrange(i)
    
    # Pull out covid data for this state and attribute with i, j, and t
    this_covid_df <- covid_df %>%
      filter(geoid %in% this_county_df$geoid) %>%
      left_join(this_county_df %>% select(geoid, i, group1, group2),
                by = 'geoid')
    
    
    # Create the county-level regression data matrix
    Xdf <- this_county_df %>%
      arrange(i) %>%
      select(
        i, 
        group1,
        group2,
        pop
      ) %>%
      mutate(
        pop_s = as.vector(scale(pop)),
        lpop_s = as.vector(scale(log(pop)))
      )
    
    group1 = if(max(this_county_df$group1)>1) Xdf %>% pull(group1) else  rep(0, nrow(Xdf))
    group1_bin <- which(group1>0)
    group1_bin_id <- group1[group1_bin]
    stan_data <- list(Pop    = Xdf %>% pull(pop),
                      Y      = this_covid_df %>% pull(Y), # using modeled cases
                      Y_i    = this_covid_df %>% pull(i),
                      Y_t    = this_covid_df %>% pull(t),
                      Y_lin_idx = col_major_sub2ind(i=this_covid_df %>% pull(i),
                                                    j=this_covid_df %>% pull(t),
                                                    I=max(this_covid_df$i),
                                                    J=max(this_covid_df$t)),  
                      I      = max(this_covid_df$i),
                      T      = max(covid_df$t),
                      TPRED  = TPRED,
                      N      = nrow(this_covid_df),
                      Ki      = 1,
                      Kt      = 7, 
                      Kit     = 0,
                      Xi      = Xdf %>% select(lpop_s),
                      Xt      = X_dow,
                      Xit     = array(0,dim=c(nrow(this_covid_df)*(max(covid_df$t)+7),0)),
                      krig_wt = krig_wt,
                      krig_var_chol = krig_var_chol,
                      spl_K  = SPL_K,
                      Z_spl  = Z,
                      J1     = if(max(this_county_df$group1)>1) max(Xdf$group1) else 0,
                      J2     = max(Xdf$group2),
                      group1 = group1,
                      N1     = length(group1_bin),
                      group1_bin = group1_bin,
                      group1_bin_id = group1_bin_id,
                      group2 =  Xdf %>% pull(group2),
                      taub0_scale = .5*sd_scale,
                      taub1_scale = .5*sd_scale,
                      taub2_scale = sd_scale,
                      sample_flag = TRUE,
                      lppd_flag = FALSE,
                      post_pred = FALSE
    )
    
    for(j in 1:NCHAINS){
      
      sample_file_name <- file.path(STATE_SAMPLES_DIR, paste0('samples_grw_',j,'.csv'))
      diagnostic_file_name <- file.path(STATE_SAMPLES_DIR, paste0('diagnostic_grw_',j,'.csv'))
      
      stan_fit_list[[slot]] <- list(stan_data   = stan_data,
                                    sample_file = sample_file_name,
                                    diagnostic_file = diagnostic_file_name,
                                    Xdf = Xdf,
                                    covid_df = this_covid_df)
      slot = slot + 1
    }
    stan_dat <- stan_fit_list[[slot-1]]
    save(stan_dat, 
         file = file.path(dirname(stan_fit_list[[slot-1]]$sample_file), 'standata.RData'))
  }
  
  save(stan_fit_list, file=file.path(DATA_DIR, 'stan_fit_list.Rdata'))
  
  
  
  
} # end prep_data
