
#library(sf)
library(tidyverse)
library(lubridate)
library(rstan)
#library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
library(future)
library(furrr)
##########################
######################################################
# COMPILE MODEL
# Negative Binomial with Gaussian Random Walk on slope

stan_mod <- stan_model(file='covid-model-master/analysis-reformat/nb1_spline.stan', model_name='nb1_spline')





#################################################################
# Fit model here
######################################################
# Set up cluster here
# cl <- makeCluster(NNODES)
# registerDoParallel(cl)

# Public IP for droplet(s); this can also be a vector of IP addresses
ip <- c('172.22.3.223', '172.22.3.218', '172.22.3.222', '172.22.3.219', '172.22.3.212', 
        '172.22.3.217', '172.22.3.220', '172.22.3.209', '172.22.3.221')

# '172.22.3.211'

# Path to private SSH key that matches key uploaded to DigitalOcean
ssh_private_key_file <- "covid-model/id_rsa"

# Connect and create a cluster
cl <- makeClusterPSOCK(
  ip,
  user = "root",
  rshopts = c(
    "-o", "StrictHostKeyChecking=no",
    "-o", "IdentitiesOnly=yes",
    "-i", ssh_private_key_file
  ),
  # Command to run on each remote machine
  # --net=host allows it to communicate back to this computer
  rscript = c("sudo", "docker", "run", "--net=host", "-v", "/covidmodeldata:/covidmodeldata", 
              "rstan/covid", "Rscript"),
  dryrun = FALSE,
  verbose = TRUE
)




plan(list(tweak(cluster, workers = cl), multiprocess))

stan_fit_list_split <- clusterSplit(cl, stan_fit_list)



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



start_time <- system.time()
r <- future_map(
  
  # Map over the 10 instances
  .x = stan_fit_list_split, 
  
  .f = ~ {
    
    outer_idx <- .x
    
    future_map(
      
      # Each instance has 32 cores we can utilize
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
                                 control =         list(adapt_delta = 0.90, max_treedepth=12)) ) 
      }
    )
    
  }
)
stop_time <- system.time()

plan(sequential)
