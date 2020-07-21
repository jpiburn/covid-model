run_diagnostics <- function(cl, param_list) {
  
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
  plan(list(tweak(cluster, workers = cl), multiprocess))
 
  message('splitting up groups for cluster')
  # group split for cluster -------------------------------------------------
  diagnostic_list <- 
    tibble(State = strsplit(FIPS, split=' ')[[1]]) %>%
    mutate(
      files = future_map(.x=State, ~get_sample_paths(.x, DATA_DIR))
      ) %>%
    filter(
      map_lgl(files, function(x) (length(x)>0))
    ) %>%
    mutate(
      split_id = rep(1:(length(cl) - 3), length = n()),
      split_id = as.numeric(split_id),
      split_id = if_else(State == "48", length(cl) + 1, split_id), # texas
      split_id = if_else(State == "51", length(cl) + 2, split_id), # virginia
      split_id = if_else(State == "13", length(cl) + 3, split_id)  # georgia 
    ) %>% 
    group_by(split_id) %>%
    group_split()
  
  message('starting diagnostics...')
  # running diagnostics -----------------------------------------------------
  diagnostic_df <- future_map(
    # cluster of workers
    .x = diagnostic_list, 
    
    .f = ~ {
      
      state_df <- .x
      
      state_df <- state_df %>%
        mutate(
          chain_info = future_map(files, get_sampling_params),
          diag       = future_map(
            .x=files,
            ~read_and_diagnose(.x, warmup=WARMUP/NTHIN,
                               par_select=c(starts_with('phi'),
                                            starts_with('tau'),
                                            starts_with('log_lambda'))
            )
          )
        )
      
    state_df
    }
  )
  
  diagnostic_df <- bind_rows(diagnostic_df)
  
  message('looking for bad chains...')
  # remove any bad chains ---------------------------------------------------
  bad_chains <- 
  diagnostic_df %>%
    select(-files, -chain_info) %>%
    unnest(cols=diag) %>% 
    filter(Rhat > 1.02) %>%
    arrange(desc(Rhat))
  
  chains_to_drop <- 
    bad_chains %>% 
    select(starts_with("Rhat_drop_")) %>% 
    apply(1, function(x) which(x == min(x)))
  
  bad_chains <- 
    bad_chains %>%
    mutate(
      drop_chain_num = chains_to_drop
    ) %>%
    distinct(State, drop_chain_num)
  
  diagnostic_df <- mutate(diagnostic_df, good_files = files)
  
  if (nrow(bad_chains) > 0) {
    for (i in 1:nrow(bad_chains)) {
      
      state_i <- pull(bad_chains[i, "State"])
      drop_chain_i <- pull(bad_chains[i, "drop_chain_num"])
      
      good_files <- diagnostic_df$good_files[[which(diagnostic_df$State == state_i)]] 
      gooder_files <- good_files[setdiff(1:length(good_files), drop_chain_i)]
      
      diagnostic_df$good_files[[which(diagnostic_df$State == state_i)]] <-  gooder_files
      
      message(glue::glue('chain {drop_chain_i} was removed for State {state_i}'))
    }
  }
  

 # save output -------------------------------------------------------------
  message('saving diagnostic file')
  save(diagnostic_df, file=DIAG_DF_LOC)
}
