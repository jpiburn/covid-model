format_predictions <- function(param_list) {
  # load libraries ----------------------------------------------------------
  require(tidyverse)
  require(lubridate)
  require(rstan)
  require(foreach)
  require(doParallel)
  require(covidmodeldata)
  require(future)
  require(furrr)
  require(sf)
  
  # define params in environment --------------------------------------------
  DATE <- param_list$DATE
  NYT_FILE <- param_list$NYT_FILE
  ACS_FILE <- param_list$ACS_FILE
  DATE_0 <- param_list$DATE_0
  DATE_N <- param_list$DATE_N
  SAMPLES_ROOT <- param_list$SAMPLES_ROOT
  NWORKERS <- param_list$NWORKERS
  NNODES <- param_list$NNODES
  RESULTS_DIR <- param_list$RESULTS_DIR
  RESULTS_FILE <- param_list$RESULTS_FILE
  UTILS_FUNC_FILE <- param_list$UTILS_FUNC_FILE
  STAN_MODEL_FILE <- param_list$STAN_MODEL_FILE
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
  load(RESULTS_FILE)
  
  df_export <- 
    data_out %>%
    transmute(
      geoid,
      state_name,
      county_name,
      date,
      obs_type = if_else(date <= as.Date(DATE), 'observed', 'predicted'),
      observed_new_cases = new_cases,
      estimated_trend_q90_lower  = lambda_q05 * pop,
      estimated_trend_q70_lower  = lambda_q15 * pop,
      estimated_trend_q50_lower  = lambda_q25 * pop,
      estimated_trend_median     = lambda_q50 * pop,
      estimated_trend_q50_higher = lambda_q75 * pop,
      estimated_trend_q70_higher = lambda_q85 * pop,
      estimated_trend_q90_higher = lambda_q95 * pop,
      estimated_new_cases_q90_lower  = Ysim_q05,
      estimated_new_cases_q70_lower  = Ysim_q15,
      estimated_new_cases_q50_lower  = Ysim_q25,
      estimated_new_cases_median     = Ysim_q50,
      estimated_new_cases_q50_higher = Ysim_q75,
      estimated_new_cases_q70_higher = Ysim_q85,
      estimated_new_cases_q90_higher = Ysim_q95
    )
  
  dir.create(glue::glue('{RESULTS_DIR}'), showWarnings = FALSE)
  readr::write_csv(df_export, glue::glue('{RESULTS_DIR}/{DATE}_output_predictions.csv.gz'))
}


