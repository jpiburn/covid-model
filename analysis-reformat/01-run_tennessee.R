library(sf)
library(tidyverse)
library(lubridate)
library(rstan)
library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
################################################################################

source('analysis-reformat/00-PARAMS.R')
NNODES <- NCHAINS

if(CLEAN_DIR) unlink(DATA_DIR, recursive=TRUE)

if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(skip_assignment = c("09","44")) # don't assign Rhode Island cases 
  
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
  filter(date >= (first_date - ZERO_PAD)) %>%
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
# Our metro hierarchy will be a metro/nonmetro, then by csa within metros
# Csa do include micropolitans, but we'll exclude those from metro
# First, recode metros to try to delete one county groups
metro_recode_df <-
  county_df %>% 
  select(geoid, mygeoid, mystate) %>%
  na.omit() %>%
  left_join(
    geodf %>% 
      st_drop_geometry()  %>%
      select(geoid, state_fips, csa_code, csa_title, metro = metropolitan_micropolitan_statistical_area),
    by = c('geoid')
  ) %>%
  mutate(
    metro = forcats::fct_explicit_na(metro) == 'Metropolitan Statistical Area'
    ) %>%
  # Change Connecticut, whch has one non-metro county (Litchfield) to metro.
  mutate(metro = if_else(mygeoid == '09005', TRUE, metro)) %>%
  mutate(
    my_metro_code= as.character(forcats::fct_explicit_na(csa_code, na_level='998'))
    )  %>%
  mutate(
    my_metro_title = if_else(is.na(csa_title), 'not_csa', csa_title)
    ) %>% 
  mutate(
    my_metro_code = ifelse(metro, my_metro_code, '999')
    ) %>%
  mutate(
    my_metro_title = if_else(metro, my_metro_title, 'not_metro')
  ) %>%
  # Calculate number of counties in each group
  group_by(
    mystate, 
    my_metro_code, 
    my_metro_title) %>%
  mutate(
    num_counties = n()) %>%
  ungroup()   %>%
  # Recode any csas with only one or two counties
  mutate(
    my_metro_code = ifelse(metro & num_counties <= 2, '998', my_metro_code),
    my_metro_title = ifelse(metro & num_counties <= 2, 'not_csa', my_metro_title),
    ) %>%
  # ReCalculate number of counties in each group
  group_by(
    mystate, 
    my_metro_code, 
    my_metro_title) %>%
  mutate(
    num_counties = n()) 



# Create a dataset with one row per metro area, and a unique id 'j'
metro_j_df <- metro_recode_df %>%
  # Summarize csas
  group_by(
    mystate, 
    my_metro_title, 
    my_metro_code) %>%
  summarize(
    n=n()) %>%
  # give index values
  arrange(mystate, my_metro_code) %>%
  group_by(
    mystate) %>%
  mutate(
    j = row_number()
  ) 

# Add j back to the county database
metro_recode_df <- metro_recode_df %>%
  left_join(
    metro_j_df,
    by = c('mystate', 'my_metro_code', 'my_metro_title')) %>%
  select(
    geoid, 
    mygeoid, 
    mystate,
    metro,
    my_metro_code,
    my_metro_title,
    j
  ) %>%
  # set j=0 if non-metro
  mutate(group1 = ifelse(metro, j, 0)) %>%
  mutate(group2 = 2-as.numeric(metro))

  

county_df <-
  county_df %>%
  left_join(
    metro_recode_df,
    by=c('geoid','mygeoid','mystate'))

  
#####################################################
# Add population and income metro status back to the county 
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
NTIME <- max(date_df$t)
D <- pracma::tril(toeplitz(c(1,-2,1, rep(0, NTIME-4 ))))
svd <- svd((D))
Z <- svd$v %*% diag(1/svd$d)

# Flip order so most important is first:
Z <- Z[,rev(1:ncol(Z))]

ZZ <- Z
# Now, take a low-rank-approximation
Z <- Z[,1:SPL_K]
# Now add the first time to it
Z <- rbind(0, Z)


## PRIOR on variance os spline
## Expect the value to grow by 3 orders of magnitude, over 90 days
## The variance seems to follow a quadratic equation with time, but I
## don't know what the parameters are.
## Trial and error suggests .02 is a good prior for the variance.

# variance on last day is:
max_var <- max(diag(chol2inv(chol(crossprod(D)))))
# sd that should give us log(10^4) change  
sd_scale <- log(10^4) / sqrt(max_var)
  
  
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

state_list = '47'


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
      group1,
      group2,
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
                    J1     = if(max(this_county_df$group1)>1) max(Xdf$group1) else 0,
                    J2     = max(Xdf$group2),
                    group1 = if(max(this_county_df$group1)>1) Xdf %>% pull(group1) else  rep(0, nrow(Xdf)),
                    group2 = Xdf %>% pull(group2),
                    taub0_scale = .5*sd_scale,
                    taub1_scale = .5*sd_scale,
                    taub2_scale = sd_scale,
                    sample_flag = TRUE,
                    lppd_flag = TRUE,
                    post_pred = TRUE
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

stan_mod <- stan_model(file='analysis-reformat/nb1_spline.stan', model_name='nb1_spline')





#################################################################
# Fit model here
######################################################
# Set up cluster here
cl <- makeCluster(NNODES)
registerDoParallel(cl)

##############################################
# Initialize function
init_fun <- function(chain_id) list( a = runif(1,-13,-9),
                                     tau_a0 = runif(1,0,.02),
                                     tau_a1 = runif(1,0,.02),
                                     tau_a2 = runif(1,0,2),
                                     raw_tau_splb0 = runif(1,0,.02),
                                     raw_tau_splb1 = runif(1,0,.02),
                                     raw_tau_splb2 = runif(1,0,2),
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
                           append_samples =  FALSE,
                           init =            init_fun,
                           control =         list(adapt_delta = 0.99, max_treedepth=12)) ) 
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

