library(sf)
library(tidyverse)
library(lubridate)
library(rstan)
library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
################################################################################

source('analysis/00-PARAMS.R')

if(CLEAN_DIR) unlink(DATA_DIR, recursive=TRUE)

if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(skip_assignment = c("44")) # don't assign Rhode Island cases 
  
} else nyt_data <- read_csv(NYT_FILE)
nyt_data <- mutate(nyt_data, date=as_date(date))


if (is.null(ACS_FILE)) {
  acs_data <- acs_data # from the covidmodeldata package 
  } else {acs_data <- readr::read_rds(ACS_FILE)}
  
geodf <- acs_data


#####################################################
# Calculate the first day of record for each county
first_day_df <- nyt_data %>%
  group_by(geoid) %>%
  summarize(first_date = min(date))

#####################################################
# Create date, geoid, and geoid/date templates
last_data_day <- max(nyt_data$date)
date_df <- tibble(date = as_date(seq.Date(from = as_date(DATE_0),
                                          to = last_data_day + TPRED,
                                          by = 1)))
county_df <- acs_data %>%
  st_drop_geometry() %>%
  select(geoid) %>%
  arrange(geoid)

geo_date_template <- crossing(county_df, date_df)

#####################################################
# Clean covid data and add zero records for one week prior to first observation
covid_df <- geo_date_template %>%
  left_join(first_day_df, by='geoid') %>%
  left_join(nyt_data, by=c('geoid', 'date') ) %>%
  filter(date >= (first_date - 7)) %>%
  mutate(new_cases_mdl = ifelse(date < first_date, 0, new_cases_mdl)) %>%
  select(geoid, date, total_cases:new_deaths_mdl)
  
################################################
# There are some negatives in the cases.
# These appear to be corrections to previous records. For now, just delete them. It will create some bias
covid_df <- covid_df %>%
  mutate(new_cases_mdl = pmax(new_cases_mdl, 0))

################################################
# Add DC to Maryland
covid_df <- covid_df %>%
  mutate(mygeoid = ifelse(geoid=='11001', '24XDC', geoid))


################################################
# Create a sequential county id:
county_i_df <- covid_df %>%
  select(mygeoid, geoid) %>%
  mutate(mystate=substring(mygeoid,1,2)) %>%
  group_by(mygeoid) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  group_by(mystate) %>%
  mutate(i=row_number())
  
county_df <- 
  county_df %>%
  left_join(
    county_i_df,
    by='geoid')
  


##################################################
# Create a time index
date_df <- date_df %>%
  mutate(t = as.integer(1 + date - as_date(DATE_0)))


##################################################
# Create a metro code within each state:
# First, recode metros t otry to delete one county groups
metro_recode_df <-
  county_df %>% select(geoid, mygeoid, mystate) %>%
  na.omit() %>%
  left_join(
    geodf %>% 
      st_drop_geometry()  %>%
      select(geoid, state_fips, csa_code, csa_title),
    by = c('geoid')
  ) %>%
  # Replace missings with non-cbsa
  mutate(
    csa_code = ifelse(is.na(csa_code), '999', csa_code),
    csa_title = ifelse(is.na(csa_title), 'Non-csa', csa_title))   %>%
  # Calculate number of counties in each 
  group_by(
    mystate, 
    csa_code, 
    csa_title) %>%
  mutate(
    num_counties = n()) %>%
  ungroup() %>%
  # Recode any csas with only one county
  mutate(
    csa_code = ifelse(num_counties < 2, '999', csa_code),
    csa_title = ifelse(num_counties < 2, 'Non-csa', csa_title))

metro_j_df <- metro_recode_df %>%
  # Summarize csas
  group_by(
    mystate, 
    csa_title, 
    csa_code) %>%
  summarize(
    n=n()) %>%
  # give index values
  arrange(mystate, csa_code) %>%
  group_by(
    mystate) %>%
  mutate(
    j = row_number()
  )
  
# Add j back to the county database
metro_recode_df <- metro_recode_df %>%
  left_join(
    metro_j_df,
    by = c('mystate', 'csa_code', 'csa_title')) %>%
  select(
    geoid, 
    mygeoid, 
    mystate,
    csa_code,
    csa_title,
    j
  )
  

county_df <-
  county_df %>%
  left_join(
    metro_recode_df,
    by=c('geoid','mygeoid','mystate'))

  
#####################################################
# Add population and income back to the county 
county_df <- county_df %>%
  left_join(acs_data %>% 
              st_drop_geometry() %>%
              select(
                geoid, 
                state_name, 
                county_name, 
                pop = acs_total_pop_e, 
                inc = acs_median_income_e,
                age = acs_median_age_e),
            by = 'geoid')


# Rearrange for convenience
covid_df <- covid_df %>%
  select(geoid, mygeoid, date, everything())




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
# This will require specifying slope separately
D <- diff(diag(max(date_df$t)),differences = 2)
#D[T-2,] <- D[T-2,]*100 # really penalize curves at the end of time
D.svd <- svd(t(D))
Z <- D.svd$u %*% diag(1/D.svd$d)
# Now, take a low-rank-approximation
Z <- Z[,(ncol(Z)+1-SPL_K):ncol(Z)]
####################################
## RW method 2
## intercept + tp*beta, beta ~ N(0,sigma)
## truncated polynomial basis
## This requires that there by NO separate slope
#tp <- spam::toeplitz.spam(y=seq(1:(T+TPRED))-1, x=rep(0,(T+TPRED)))
#tp <- tp[,1:(T+TPRED-1)]
## qr decomposition
#tp.qr <- qr(tp, LAPACK=FALSE)
#Q <- qr.Q(tp.qr)
#R <- qr.R(tp.qr)
## svd decomp of R
#R.svd <- svd(R)
#Z <- Q %*% R.svd$u %*% diag(R.svd$d)
#Z <- Z[1:spl_K]
## Now, we can do intercept + Za, a~N(0,sigma)
## compared to tp, Z has orthogonal columns
## PRIOR on variance os spline
## Expect the value to grow by 3 orders of magnitude over 30 days
## The variance seems to follow a quadratic equation with time, but I
## don't know what the parameters are.
## Trial and error suggests .02 is a good prior for the variance.


if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive=TRUE)
save(covid_df, county_df, X_dow, Z, date_df, 
     file=file.path(DATA_DIR, 'data_frames.Rdata'))




# We need a list of states to loop through.
# NO LONGER NECESSARY Check the coviddf to make sure that there are enough cases statewide to model.
# Except for DC (11), which we'll process with Maryland (24)

# Step 1. get a list of states from the NYT file, 
# and order from most counties to least.
state_list <- nyt_data %>%
  group_by(state_fips) %>%
  distinct(county_name) %>%
  summarize(n=n()) %>%
  filter(!(state_fips %in%  c('11'))) %>% # DC 
  filter(!is.na(state_fips)) %>%
  arrange(desc(n)) %>%
  pull(state_fips)

#state_list <- state_list[c(20:25)]

#state_list <- c('15','44','02','56','30')



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
    left_join(this_county_df %>% select(geoid, i, j),
              by = 'geoid') %>%
    left_join(date_df %>% select(date, t), 
              by='date')
  
  
  # Create the county-level regression data matrix
  Xdf <- this_county_df %>%
    arrange(i) %>%
    select(
      i, 
      j,
      pop, 
      inc
      ) %>%
    mutate(
      pop_s = as.vector(scale(pop)),
      inc_s = as.vector(scale(inc)),
      lpop_s = as.vector(scale(log(pop)))
      )
  
  
  stan_data <- list(Pop    = Xdf %>% pull(pop),
                    Y      = this_covid_df %>% filter(new_cases_mdl>0) %>% pull(new_cases_mdl), # using modeled cases
                    Y_i    = this_covid_df %>% filter(new_cases_mdl>0) %>% pull(i),
                    Y_t    = this_covid_df %>% filter(new_cases_mdl>0) %>% pull(t),
                    Y0_i   = this_covid_df %>% filter(new_cases_mdl==0) %>% pull(i),
                    Y0_t   = this_covid_df %>% filter(new_cases_mdl==0) %>% pull(t),
                    I      = max(this_covid_df$i),
                    T      = max(date_df$t),
                    N      = nrow(this_covid_df %>% filter(new_cases_mdl>0) ),
                    N0     = nrow(this_covid_df %>% filter(new_cases_mdl==0) ),
                    K      = 1,
                    X      = Xdf %>% select(lpop_s),
                    X_dow  = X_dow,
                    spl_K  = SPL_K,
                    Z_spl  = Z,
                    J1     = if(max(this_county_df$j)>1) max(Xdf$j) else 0,
                    group1 = if(max(this_county_df$j)>1) Xdf %>% pull(j) else  rep(0, nrow(Xdf)),
                    taub_scale = .05,
                    sample_flag = TRUE
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
  save(stan_dat, 
       file = file.path(dirname(stan_fit_list[[slot-1]]$sample_file), 'standata.RData'))
}


######################################################
# COMPILE MODEL
# Negative Binomial with Gaussian Random Walk on slope

stan_mod <- stan_model(file='analysis/nb1_spline.stan', model_name='nb1_spline')





#################################################################
# Fit model here
######################################################
# Set up cluster here
cl <- makeCluster(NNODES)
registerDoParallel(cl)

##############################################
# Initialize function
init_fun <- function(chain_id) list( a = runif(1,-13,-9),
                                     b = runif(1,0,.4),
                                     tau_a0 = runif(1,0,2),
                                     tau_a1 = runif(1,0,2),
                                     raw_tau_b0 = runif(1,0,2), 
                                     raw_tau_splb = runif(1,0,2),
                                     phi = runif(1,.1,.9))


r <- foreach(i = 1:length(stan_fit_list),
             .verbose=TRUE, .packages='rstan') %dopar% {
  stan_fit <- try(sampling(object =          stan_mod,
                           data =            stan_fit_list[[i]]$stan_data,
                           iter =            NITER, 
                           thin =            NTHIN,
                           chains =          1,
                           cores =           1,
                           warmup =          WARMUP*NTHIN,
                           sample_file =     stan_fit_list[[i]]$sample_file,
                           diagnostic_file = stan_fit_list[[i]]$diagnostic_file,
                           append_samples =  FALSE,
                           init =            init_fun,
                           control =         list(adapt_delta = 0.99, max_treedepth=10)) ) 
  #return(1)
}
stopCluster(cl)

#stan_opt <- try(optimizing(object =          stan_mod,
#                         data =            stan_fit_list[[i]]$stan_data,
#                         init =            init_fun(),
#                         init_alpha =      1e-8,
#                         tol_rel_grad =    1e4,
#                         algorithm =       "BFGS",
#                         verbose=          TRUE)) 

#stan_pars <- tibble(value = stan_opt$par, variable = names(stan_opt$par))

