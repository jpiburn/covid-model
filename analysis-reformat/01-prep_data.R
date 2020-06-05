#library(sf)
library(tidyverse)
library(lubridate)
library(rstan)
#library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
################################################################################

source('analysis-reformat/00-PARAMS.R')
source('analysis-reformat/00-functions.R')
source('analysis-reformat/00-proc_covid.R')
source('analysis-reformat/00-proc_county.R')

#if(CLEAN_DIR) unlink(DATA_DIR, recursive=TRUE)

if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(distribute_unknowns = FALSE) # Assign KC, but that's it.
  
  #    format_nyt(skip_assignment = c("09","44")) # don't assign Rhode Island cases 
  
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
  DATE_0 = DATE_0,
  ZERO_PAD = ZERO_PAD
) %>%
  left_join(nyt_data %>% select(geoid, date, new_cases))


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

# D11 <- (Dbig%*%t(Dbig))[1:(NTIME-1-TPRED),1:(NTIME-1-TPRED)]
# D12 <- (Dbig%*%t(Dbig))[(NTIME-1-TPRED+1):(NTIME-1),1:(NTIME-1-TPRED)]
# D22 = (Dbig%*%t(Dbig))[(NTIME-1-TPRED+1):(NTIME-1), (NTIME-1-TPRED+1):(NTIME-1)]
# 
# pred_cov <- 
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
sd_scale <- log(10^4) / sqrt(max_var)


# if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive=TRUE)
save(covid_df, county_df, X_dow, Z, date_df,
     file=file.path(DATA_DIR, 'data_frames.Rdata'))
file.copy(from = 'analysis-reformat/00-PARAMS.R',
          to = file.path(DATA_DIR, '00-PARAMS.R'))




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

#state_list <- state_list[c(20:25)]

#state_list <- c('15','44','02','56','30')
# state_list <- '09'

if(TN_ONLY) state_list <- '47'


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
    stan_fit_list[[slot]] <- list(stan_data=stan_data,
                                  sample_file = file.path(STATE_SAMPLES_DIR, paste0('samples_grw_',j,'.csv')),
                                  diagnostic_file = file.path(STATE_SAMPLES_DIR, paste0('diagnostic_grw_',j,'.csv')),
                                  Xdf = Xdf,
                                  covid_df = this_covid_df)
    slot = slot + 1
  }
  stan_dat <- stan_fit_list[[slot-1]]
  save(stan_dat, file = file.path(dirname(stan_fit_list[[slot-1]]$sample_file), 'standata.RData'))
}



