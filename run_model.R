#!/usr/bin/env Rscript
#' @title     this is the command line script to create the models
#' 
# DATE <- '2020-06-24' # Date of run -- replaced with below
# DATE <- format(Sys.Date(), format = '%Y-%m-%d')
library(data.table)
library(snakecase)
library(tidycensus)
library(tidyverse)
library(lubridate)
library(rstan)
library(foreach)
library(doParallel)
library(covidmodeldata)
library(future)
library(furrr)

# This is used b/c of weird permission issues when running the sricpt through 
# Rstudio vs R cmd probaly won't need it

system("sudo Rscript analysis-reformat/01-prep_data_only.R")

# run model
system("sudo Rscript covid-model-master/analysis-reformat/02-run_diagnostics.R")
system("sudo Rscript covid-model-master/analysis-reformat/03-extract_summary.R")


#' Get Latest COVID-19 Data from New York Times Database
#'
#' @param level character. county level data or state level
#'
#' @return df
#' @export
get_nyt <- function(level = "county") {
  
  if (level == "county") {
    url <- "https://github.com/nytimes/covid-19-data/raw/master/us-counties.csv"
  }
  
  if (level == "state") {
    url <- "https://github.com/nytimes/covid-19-data/raw/master/us-states.csv"
  }
  
  df <- readr::read_csv(url)
  df <- janitor::clean_names(df)
  
  df
}

# improve the retrieval of the new york times data
# > format(object.size(dat), 'MB')
# [1] "10.6 Mb"
# > format(object.size(df), 'MB')
# [1] "14.2 Mb"
# data.table preferred
retrieve_acs_data <- function() {
  geoid <- affgeoid <- statefp <- countyfp <- name_x <- name_y <- 
    fips_state_code <- fips_county_code <- county_county_equivalent <- NULL
    
  options(tigris_use_cache = TRUE)
  api_key <- Sys.getenv("CENSUS_API_KEY")
  tidycensus::census_api_key(api_key)
  msa_url <- paste0("https://www2.census.gov/programs-surveys/",
                    'metro-micro/geographies/reference-files/2018/',
                    "delineation-files/list1_Sep_2018.xls")
  msa_file <- "data-raw/msa-list.xls"
  download.file(msa_url, msa_file, mode = "wb")
  # download acs ------------------------------------------------------------
  acs_vars <- c(
    acs_total_pop = "B25026_001",
    acs_median_income = "B07011_001",
    acs_median_age = "B01002_001"
  )
  
  acsdf <- get_acs(
    geography = "county", 
    variables = acs_vars,
    state = NULL, # get everystate
    geometry = TRUE, 
    keep_geo_vars = TRUE,
    output = "wide"
  ) %>%
    janitor::clean_names() %>%
    tidy::select(
      geoid = geoid,
      geoid_full = affgeoid,
      state_fips = statefp,
      county_fips = countyfp,
      county_name = name_x,
      county_name_full = name_y,
      starts_with("acs_")
    )
  # join counties to msa data -----------------------------------------------
  msadf <- readxl::read_xls(msa_file, skip = 2) %>%
    janitor::clean_names() %>%
    tidy::mutate(
      geoid = paste0(fips_state_code, fips_county_code)
    )
  
  acsdf <- tidy::left_join(acsdf, msadf, by = "geoid") %>%
    tidy::select(
      -(county_county_equivalent:fips_county_code))
  
  acsdf
}
# Ingest the data ---------------------------------------------------------
retrieve_nyt_data <- function(level = c('county', 'state')) {
  if (level == 'county') url <-
      "https://github.com/nytimes/covid-19-data/raw/master/us-counties.csv"
  else if (level == "state") url <- 
      "https://github.com/nytimes/covid-19-data/raw/master/us-states.csv"
  dat <- data.table::fread(url)
  data.table::setnames(dat, names(dat), snakecase::to_snake_case(names(dat)))
  dat
}

nyt_county <- retrieve_nyt_data('county')
nyt_state <- retrieve_nyt_data('state')
msa_list <- read_excel("data-raw/msa-list.xls")

nyt_county_datemax <- max(nyt_county[['date']])
nyt_state_datemax <- max(nyt_state[['date']])
acs_data <- retrieve_acs_data()

# create the model session information ------------------------------------
model_session_info <- list(
  DATE = nyt_county_datemax,
  TN_ONLY = FALSE,
  NYT_FILE = NULL,
  ACS_FILE = NULL,
  DATE_0 = '2020-03-01', # First date to use
  DATE_N = nyt_county_datemax, # Last date to use 
  SAMPLES_ROOT = '/covidmodeldata',
  NWORKERS = 0,
  NNODES = 32,
  RESULTS_DIR = file.path('/covidmodeldata/results'),
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
  # NTHIN = 2,
  # NCHAINS = 8,
  # WARMUP = (2/8)/2,
  WARMUP = 500, # WARMUP is number of retained (before thinning) samples
  DATA_DIR = file.path('/covidmodeldata', nyt_county_datemax),
  DIAG_DF_LOC = file.path('/covidmodeldata', 'diagnostic.Rdata'),
  ZERO_PAD = 7, # Number of days of zeros to add before first obs for each county
  # Number of days forward to predict (from last data day, not from Sys.Date())
  TPRED = 7,
  # My logic here is we want fewer than one peak per week.
  SPL_K = floor(as.numeric(as.Date(nyt_county_datemax) - 
                             as.Date('2020-03-01'))/14) 
  
)

# data prep ---------------------------------------------------------------
library(covidmodeldata)
# folder with date, directory by fips code (names must stay the same)
# 

