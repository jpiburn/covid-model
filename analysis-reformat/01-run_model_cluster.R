# load libraries ----------------------------------------------------------
library(tidyverse)
library(lubridate)
library(rstan)
library(foreach)
library(doParallel)
library(covidmodeldata)
library(future)
library(furrr)


# load data and parameters ------------------------------------------------
source('analysis-reformat/00-PARAMS.R')
load(file.path(DATA_DIR, "stan_fit_list.Rdata"))


# create cluster connections ----------------------------------------------

worker_ips <- cades_workers$ip_address # from 00-PARAMS.R
worker_ips <- setdiff(worker_ips, c('172.22.3.211','172.22.3.204', '172.22.3.203'))

# Connect and create a cluster
cl <- makeClusterPSOCK(
  worker_ips,
  user = "root",
  rshopts = c(
    "-o", "StrictHostKeyChecking=no",
    "-o", "IdentitiesOnly=yes",
    "-i", SSH_KEY
  ),
  # Command to run on each remote machine
  # --net=host allows it to communicate back to this computer
  rscript = c("sudo", "docker", "run", "--net=host", "-v", "/covidmodeldata:/covidmodeldata", 
              "rstan/covid", "Rscript"),
  dryrun = FALSE,
  verbose = TRUE
)

plan(list(tweak(cluster, workers = cl), multiprocess))


# split stan fit list -----------------------------------------------------
stan_fit_list_split <- clusterSplit(cl, stan_fit_list)


# compile stan model ------------------------------------------------------
stan_mod <- stan_model(file='analysis-reformat/nb1_spline.stan', model_name='nb1_spline')


# initialize function -----------------------------------------------------
init_fun <- function(chain_id) list( a = runif(1,-13,-9),
                                     tau_a0 = runif(1,0,.02),
                                     tau_a1 = runif(1,0,.02),
                                     tau_a2 = runif(1,0,2),
                                     raw_tau_splb0 = runif(1,0,.02),
                                     raw_tau_splb1 = runif(1,0,.02),
                                     raw_tau_splb2 = runif(1,0,2),
                                     phi = runif(1,.1,.9))


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
                                 control =         list(adapt_delta = 0.90, max_treedepth=12))
                        )

        readr::write_lines(x = outer_idx[[inner_idx]]$sample_file, file.path(DATA_DIR, "chain_checkout.txt"), append = TRUE)
        
        }
    )
    
  }
)


# closing cluster ---------------------------------------------------------
plan(sequential)




