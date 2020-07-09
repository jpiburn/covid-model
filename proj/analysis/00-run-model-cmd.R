
# This is used b/c of weird permission issues when running the sricpt through Rstudio vs R cmd
# probaly won't need it

library(tidyverse)
library(lubridate)
library(rstan)
library(foreach)
library(doParallel)
library(covidmodeldata)
library(future)
library(furrr)

source('analysis-reformat/00-PARAMS.R')

system("sudo Rscript analysis-reformat/01-prep_data_only.R")

# run model
system("sudo Rscript covid-model-master/analysis-reformat/02-run_diagnostics.R")
system("sudo Rscript covid-model-master/analysis-reformat/03-extract_summary.R")
