
run_model <- function(cl, param_list) {
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
  
  # load data and parameters ------------------------------------------------
  load(file.path(DATA_DIR, "stan_fit_list.Rdata"))
  
  
  # tweak cluster plan ------------------------------------------------------
  plan(list(tweak(cluster, workers = cl), multiprocess))
  
  
  # split stan fit list -----------------------------------------------------
  stan_fit_list_split <- clusterSplit(cl, stan_fit_list)
  
  message('compiling model...')
  # compile stan model ------------------------------------------------------
  stan_mod <- stan_model(file = STAN_MODEL_FILE, model_name = 'nb1_spline')
  
  
  # initialize function -----------------------------------------------------
  init_fun <- function(chain_id) list( a = runif(1,-13,-9),
                                       tau_a0 = runif(1,0,.02),
                                       tau_a1 = runif(1,0,.02),
                                       tau_a2 = runif(1,0,2),
                                       raw_tau_splb0 = runif(1,0,.02),
                                       raw_tau_splb1 = runif(1,0,.02),
                                       raw_tau_splb2 = runif(1,0,2),
                                       phi = runif(1,.1,.9))
  
  message('running model...')
  # running model -----------------------------------------------------------
  r <- future_map(
    # cluster of workers
    .x = stan_fit_list_split, 
    
    .f = ~ {
      
      outer_idx <- .x
      
      future_map(
        
        # multiprocess within each worker
        .x = 1:length(outer_idx), 
        
        .f = ~ {
          inner_idx <- .x
          stan_fit <- try(sampling(object =          stan_mod,
                                   data =            outer_idx[[inner_idx]]$stan_data,
                                   iter =            NITER,
                                   thin =            NTHIN,
                                   chains =          1,
                                   cores =           1,
                                   warmup =          WARMUP,
                                   sample_file =     outer_idx[[inner_idx]]$sample_file,
                                   append_samples =  FALSE,
                                   init =            init_fun,
                                   control =         list(adapt_delta = 0.90, max_treedepth = 12))
          )
          
          readr::write_lines(x = outer_idx[[inner_idx]]$sample_file, file.path(DATA_DIR, "chain_checkout.txt"), append = TRUE)
          
        }
      )
      
    }
  )
  
} # end run_model()





