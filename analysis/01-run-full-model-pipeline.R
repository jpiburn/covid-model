# load libraries ----------------------------------------------------------
library(tidyverse)
library(lubridate)
library(rstan)
library(foreach)
library(doParallel)
library(covidmodeldata)
library(future)
library(furrr)

# I'm not really satistified with how we set the date right now.
# I know there is a smarter way, but we'll get there soon. fine for now
DATE <- Sys.Date() - 1
# DATE <- '2020-07-08'
SAMPLES_ROOT = '/covidmodeldata'
DATA_DIR = file.path(SAMPLES_ROOT, DATE)
RESULTS_DIR <- file.path(SAMPLES_ROOT, 'results', DATE)

DATE_0 = '2020-03-01' # First date to use.
# This is something we'll need to test going forward.
# the model takes longer and longer to run b/c the # of days always increases. 
# we should experiment will moving this day and see how results change.

param_list <- list(
  DATE =  DATE, # last observered day
  NYT_FILE = NULL,
  ACS_FILE = NULL,
  DATE_0 = '2020-03-01', # First date to use
  DATE_N = DATE, # Last date to use 
  SAMPLES_ROOT = '/covidmodeldata',
  NWORKERS = 0,
  NNODES = 32,
  RESULTS_DIR = RESULTS_DIR,
  RESULTS_FILE = file.path(RESULTS_DIR, paste0('results_', DATE,'.RData')),
  UTILS_FUNC_FILE = '/home/cades/covid-model/analysis/model-steps/00-functions.R',
  STAN_MODEL_FILE = '/home/cades/covid-model/analysis/nb1_spline.stan',
  CLEAN_DIR = TRUE,
  SSH_KEY = "/home/cades/covid-model/.ssh/id_rsa",
  cades_workers = 
    tibble::tribble(
      ~instance_name,    ~ip_address,
      "covid-worker-20", "172.22.3.223",
      "covid-worker-19", "172.22.3.218",
      "covid-worker-18", "172.22.3.222",
      "covid-worker-17", "172.22.3.219",
      "covid-worker-16", "172.22.3.212",
      "covid-worker-15", "172.22.3.217",
      "covid-worker-14", "172.22.3.220",
      "covid-worker-13", "172.22.3.211",
      "covid-worker-12", "172.22.3.209",
      "covid-worker-11", "172.22.3.221",
      "covid-worker-10", "172.22.3.216",
      "covid-worker-9", "172.22.3.214",
      "covid-worker-8", "172.22.3.208",
      "covid-worker-7", "172.22.3.210",
      "covid-worker-6", "172.22.3.215",
      "covid-worker-5", "172.22.3.207",
      "covid-worker-4", "172.22.3.206",
      "covid-worker-3", "172.22.3.205",
      "covid-worker-2", "172.22.3.204",
      "covid-worker-1", "172.22.3.203"
    ),
  NITER = 1500, # Number of  iterations (before thinning)
  NTHIN = 2,
  NCHAINS = 8,
  WARMUP = 500, # WARMUP is number of retained (before thinning) samples
  DATA_DIR = file.path(SAMPLES_ROOT, DATE),
  DIAG_DF_LOC = file.path(DATA_DIR, 'diagnostic.Rdata'),
  ZERO_PAD = 7, # Number of days of zeros to add before first obs for each county
  TPRED = 7, # Number of days forward to predict (from last data day, not from Sys.Date())
  SPL_K = floor(as.numeric(as.Date(DATE)-as.Date(DATE_0))/14) # My logic here is we want fewer than one peak per week.
)


# load needed functions ---------------------------------------------------
fnc_scripts <- list.files('analysis/model-steps/', full.names = TRUE)
lapply(fnc_scripts, source)

# create cluster connections ----------------------------------------------
worker_ips <- param_list$cades_workers$ip_address # from 00-PARAMS.R

# Connect and create a cluster
cl <- makeClusterPSOCK(
  worker_ips,
  user = "root",
  rshopts = c(
    "-o", "StrictHostKeyChecking=no",
    "-o", "IdentitiesOnly=yes",
    "-i", param_list$SSH_KEY
  ),
  # Command to run on each remote machine
  # --net=host allows it to communicate back to this computer
  rscript = c("sudo", "docker", "run", "--net=host", "-v", "/covidmodeldata:/covidmodeldata", 
              "rstan/covid", "Rscript"),
  dryrun = FALSE,
  verbose = FALSE 
)



# run each model step -----------------------------------------------------
prep_data(param_list)

run_model(cl, param_list)

run_diagnostics(cl, param_list)

extract_summary(cl, param_list)

gen_report_maps(param_list)

format_predictions(param_list)

# closing cluster ---------------------------------------------------------
plan(sequential)