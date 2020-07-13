extract_summary <- function(cl, param_list) { 
  # load libraries ----------------------------------------------------------
  require(tidyverse)
  require(lubridate)
  require(rstan)
  require(foreach)
  require(doParallel)
  require(covidmodeldata)
  require(future)
  require(furrr)
  
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
  
  FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'
  load(file=file.path(DATA_DIR, 'data_frames.Rdata'))
  load(DIAG_DF_LOC)
  plan(list(tweak(cluster, workers = cl), multiprocess))
  

  # group split for cluster -------------------------------------------------
  summary_df <- diagnostic_df %>% select(State, good_files)
  
  summary_list <- 
    summary_df %>%
    mutate(
      split_id = rep(1:(length(cl) - 3), length = n()),
      split_id = as.numeric(split_id),
      split_id = if_else(State == "48", length(cl) + 1, split_id), # texas
      split_id = if_else(State == "51", length(cl) + 2, split_id), # virginia
      split_id = if_else(State == "13", length(cl) + 3, split_id)  # georgia 
    ) %>% 
    group_by(split_id) %>%
    group_split()
  
  
  # summarize lambda (predicted new cases per capita) -----------------------
  summary_lambda <- future_map(
    # cluster of workers
    .x = summary_list, 
    
    .f = ~ {
      
      state_df <- .x
      
      summary_lambda_i <- 
        state_df %>%
        mutate(
          lambda_summary = future_map(good_files, function(x) { 
            if(length(x)==0){
              return(tibble(lambda_mean = NA,
                            lambda_q05 = NA,
                            lambda_q15 = NA,
                            lambda_q25 = NA,
                            lambda_q50 = NA,
                            lambda_q75 = NA,
                            lambda_q85 = NA,
                            lambda_q95 = NA))
            }
            
            samples <- read_stan_draws(x, par_select=c(starts_with("log_lambda")), warmup = WARMUP)
            
            summary <- samples %>% 
              group_by(.variable) %>%
              summarize(lambda_mean = mean(exp(.value)),
                        lambda_q05 = exp(quantile(.value, .05)),
                        lambda_q15 = exp(quantile(.value, .15)),
                        lambda_q25 = exp(quantile(.value, .25)),
                        lambda_q50 = exp(quantile(.value, .5)),
                        lambda_q75 = exp(quantile(.value, .75)),
                        lambda_q85 = exp(quantile(.value, .85)),
                        lambda_q95 = exp(quantile(.value, .95))
              ) 
            
            summary <- summary %>%
              separate(col = .variable, into=c('.variable','i','t'), sep = '\\.') %>%
              mutate(
                i = as.integer(i), 
                t=as.integer(t)) %>%
              mutate(
                .variable = 'lambda'
              )
            
            return(summary)
          })
        ) %>%
        select(
          State, 
          lambda_summary
        ) %>%
        unnest(lambda_summary) %>%
        select(-.variable) 
      
      summary_lambda_i
    }
  )
  
  summary_lambda <- bind_rows(summary_lambda)
  
  
  # summarize Ysim (predicted new cases) ------------------------------------
  summary_Ysim <- future_map(
    # cluster of workers
    .x = summary_list, 
    
    .f = ~ {
      
      state_df <- .x
      
      summary_Ysim_i <- 
        state_df %>%
        mutate(
          Ysim_summary = future_map(good_files, function(x){
            if(length(x)==0){
              return(tibble(Ysim_mean = NA,
                            Ysim_q05 = NA,
                            Ysim_q15 = NA,
                            Ysim_q25 = NA,
                            Ysim_q50 = NA,
                            Ysim_q75 = NA,
                            Ysim_q85 = NA,
                            Ysim_q95 = NA))
            }
            samples <- read_stan_draws(x, par_select=c(starts_with("Y_sim")), warmup = WARMUP)
            summary <- samples %>% 
              group_by(.variable) %>%
              summarize(Ysim_mean = mean(exp(.value)),
                        Ysim_q05 = as.integer(quantile(.value, .05)),
                        Ysim_q15 = as.integer(quantile(.value, .15)),
                        Ysim_q25 = as.integer(quantile(.value, .25)),
                        Ysim_q50 = as.integer(quantile(.value, .5)),
                        Ysim_q75 = as.integer(quantile(.value, .75)),
                        Ysim_q85 = as.integer(quantile(.value, .85)),
                        Ysim_q95 = as.integer(quantile(.value, .95))) 
            summary <- summary %>%
              separate(col = .variable, into=c('.variable','i','t'), sep = '\\.') %>%
              mutate(
                i = as.integer(i), 
                t=as.integer(t)) %>%
              mutate(
                .variable = 'Ysim')
            return(summary)
          }
          )
        ) %>%
        select(
          State, 
          Ysim_summary) %>%
        unnest(
          Ysim_summary) %>%
        select(-.variable)
      
      summary_Ysim_i
    }
  )
  
  summary_Ysim <- bind_rows(summary_Ysim)
  

  # create data_out (output predictions) ------------------------------------
  data_out <- 
    crossing(
      county_df %>% select(geoid, mygeoid, i, mystate, state_name, county_name, pop), 
      date_df) %>%
    left_join(
      summary_lambda,
      by = c('mystate'='State', 'i', 't')
    ) %>%
    left_join(
      summary_Ysim,
      by = c('mystate'='State', 'i', 't')
    ) %>%
    left_join(
      covid_df %>% select(-pop),
      by = c('geoid', 'date','t')
    ) %>%
    select(-mygeoid, -mystate, -i, -t)
  
  
  # save output -------------------------------------------------------------
  message('saving output predictions')
  dir.create(RESULTS_DIR)
  save(data_out, file = RESULTS_FILE)
}
