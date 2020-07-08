#!/usr/bin/env Rscript
#' @title     this is the command line script to create the models
#' 
# DATE <- '2020-06-24' # Date of run -- replaced with below
# DATE <- format(Sys.Date(), format = '%Y-%m-%d')
library(data.table)
library(tidyr)
library(dplyr)
library(snakecase)
library(tidycensus)
library(tidyselect)
library(tidyverse)
library(lubridate)
library(rstan)
library(parallel)
library(doParallel)
library(foreach)
library(covidmodeldata)
library(future)
library(furrr)
library(rlang)
library(ggplot2)
library(hms)
library(sf)
# original notes --------------------------------------------------------- {{{1
# This is used b/c of weird permission issues when running the sricpt through 
# Rstudio vs R cmd probaly won't need it

#system("sudo Rscript analysis-reformat/01-prep_data_only.R")

# run model
#system("sudo Rscript covid-model-master/analysis-reformat/02-run_diagnostics.R")
#system("sudo Rscript covid-model-master/analysis-reformat/03-extract_summary.R")
# 1}}} ------------------------------------------------------------------------
# define functions ------------------------------------------------------- {{{1
# get_nyt {{{2
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
} # 2}}}
# dt justification {{{2
# improve the retrieval of the new york times data
# > format(object.size(dat), 'MB')
# [1] "10.6 Mb"
# > format(object.size(df), 'MB')
# [1] "14.2 Mb"
# data.table preferred
# 2}}}
retrieve_acs_data <- function() { # {{{2
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
} # 2}}}
get_sample_paths <- function(STATE, ROOT) { # {{{2
  # Function to take a state and root directory, and return a data_frame
  #  with one column containing all sample csv paths
  SAMP_DIR <- file.path(ROOT, STATE)
  csv.files <- list.files(SAMP_DIR, pattern='samples.*[0-9]\\.csv$',
                          full.names=TRUE)
    
  return(csv.files)
} # 2}}}
read_stan_draws <- function(files, par_select, warmup = 500) { # {{{2
  # Function to create tidy object from list of files:
  #Input: files: a dataframe
  #csv_files_col_name: character: names of column with paths to csvs
  # Output: a long dataframe with samples
  
  par_enquo <- rlang::enquo(par_select)

  samples <- tibble::tibble(files=files) %>%
    dplyr::mutate(.chain=1:n()) %>%
    dplyr::mutate(samples = map(.x = files, 
                                ~fread_stan(.x, col_select = !!par_enquo))) %>%
    dplyr::mutate(samples = map(.x=samples, ~mutate(.x,.iteration=1:n()))) %>%
    tidyr::unnest(cols=samples) %>%
    dplyr::mutate(.draw = 1:n()) %>%
    dplyr::filter(.iteration > warmup) %>%
    dplyr::select(-files) %>%
    dplyr::select(.draw, .chain, .iteration, tidyselect::everything(),) %>%
    dplyr::ungroup() %>% 
    tidyr::pivot_longer(cols = !c(.draw, .chain, .iteration),
                        names_to = '.variable', 
                        values_to = '.value')
              
  return(samples)
} # 2}}}
calculate_rhat <- function(samples, warmup = 0, par) { # {{{2
  # Function to take a dataframe of samples, and a paramater string and 
  #  calculate rhat for all matching pars (start_with(par))
  temp <- samples %>% dplyr::filter(.iteration > warmup) %>%
    dplyr::select(.chain, .iteration, tidyselect::starts_with(par)) %>%
  tidyr::pivot_longer(cols = starts_with(par), names_to = '.variable', 
                      values_to = '.value')

  NCHAINS = max(temp$.chain)
  NITER = max(temp$.iteration)

  grand_mean <- temp %>% dplyr::group_by(.variable) %>% 
    dplyr::summarize(gmean = mean(.value))

  chain_summary <- temp %>% dplyr::group_by(.chain, .variable) %>% 
    dplyr::summarize(wmean = mean(.value), wvar = var(.value)) %>%
    dplyr::ungroup() %>% dplyr::left_join(grand_mean, by = '.variable') %>%
    dplyr::group_by(.variable) %>%
    dplyr::summarize(W = mean(wvar), 
                     B = NITER /(NCHAINS-1) * sum((wmean-gmean)^2)) %>%
    dplyr::mutate(V = (1-1/NITER) * W + 1/NITER * B, rhat = sqrt(V/W)) %>%
    dplyr::summarize(max_rhat = max(rhat), 
                     which_max = .variable[which.max(rhat)],
                     q99_rhat = quantile(rhat, .99, na.rm=TRUE),
                     frac_gt_101 = mean(rhat > 1.01))
  return(chain_summary)
} # 2}}}
read_and_summarize <- function(files, par_select, warmup = 0, # {{{2
                               k = NULL) { 
  # Function to take file names and read csv and calculate rhat summary
  # INPUTS:
  #  files: a character vector of stan sample files to process
  #  warmup: # number of warmup samples
  #  par_select: a tidy select of parameters to select 
  #              (e.g. c(starts_with('b0_raw')))
  #  k: will drop chain k from calculation r-hat (default NULL)
  par_expr <- rlang::expr(par_select)
  par_enquo <- rlang::enquo(par_select)
    
  if(length(files) < 3) {
    return(tibble::tibble(max_rhat = NA*0, which_max = as(NA,'character'),
                     q99_rhat = NA*0, frac_gt_101 = NA*0))
  }

  samples <- tibble::tibble(files = files) %>% 
    dplyr::mutate(.chain = 1:n()) %>% 
    dplyr::mutate(samples = map(.x = files, 
                                ~fread_stan(.x, col_select = !!par_enquo))) %>%
    dplyr::mutate(samples = map(.x = samples, 
                                ~mutate(.x, .iteration = 1:n()))) %>%
    tidyr::unnest(cols = samples) %>%
    dplyr::mutate(.draw = 1:n()) %>%
    dplyr::select(.draw, .chain, .iteration, tidyselect::everything()) %>%
    dplyr::ungroup()

  temp <- samples %>% dplyr::filter(.iteration > warmup) %>%
    dplyr::select(.chain, .iteration, !!par_enquo) %>%
    tidyr::pivot_longer(cols = !!par_enquo, names_to = '.variable', 
                        values_to = '.value')

  NCHAINS = max(temp$.chain)
  NITER = max(temp$.iteration)

  if(!is.null(k)) {
    # Delete a chain
    temp <- temp %>% filter(.chain != k)
    NCHAINS = NCHAINS - 1
  }

  grand_mean <- temp %>% dplyr::group_by(.variable) %>% 
    dplyr::summarize(gmean = mean(.value))

  chain_summary <- temp %>% dplyr::group_by(.chain, .variable) %>% 
    dplyr::summarize(wmean = mean(.value), wvar = var(.value)) %>%
    dplyr::ungroup() %>% dplyr::left_join(grand_mean, by = '.variable') %>%
    dplyr::group_by(.variable) %>%
    dplyr::summarize(W = mean(wvar), 
                     B = NITER /(NCHAINS-1) * sum((wmean-gmean)^2)) %>%
    dplyr::mutate(V = (1-1/NITER) * W + 1/NITER * B, rhat = sqrt(V/W)) %>%
    dplyr::summarize(max_rhat = max(rhat),
                     which_max = .variable[which.max(rhat)],
                     q99_rhat = quantile(rhat, .99, na.rm=TRUE),
                     frac_gt_101 = mean(rhat > 1.01))

  return(chain_summary)
} # 2}}}
vroom_stan <- function(file, ...) { # {{{2
  # Stan sample files have commented out lines in the start, middle and end of the file
  # vroom can figure out the start, but not the middle and end
  # Use grep to delete those
  
  tfile <- paste0(file, '.tmp')
  grepcmd <- paste0("grep -vh '^#' ", file, " > ", tfile)
  
  system(grepcmd)
  
  out <- vroom::vroom(tfile, delim = ',', num_threads = 1, 
                      col_types = cols(.default = col_double()), ...)
    
  unlink(tfile)

  return(out)
} # 2}}}
fread_stan <- function(file, col_select, ...) { # {{{2
  col_enquo <- rlang::expr(col_select)
    
  # was having trouble getting vroom_stan() to work on windows machine 
  # (what the current VM is)
  # same idea as vroom_stan, but avoids making the tmp file
  # fread is pretty fast so it shouldn't be much of a slow down
  # might be faster since no tmp files. convert to tibble for compatibility
    
  grepcmd <- paste0("grep -vh '^#' ", file)
    
  col_names <- data.table::fread(cmd = grepcmd, sep = ",", nThread = 1, nrows = 0)
  #col_select <- tidyselect::eval_select(col_expr, col_names)
  col_nums <- tidyselect::eval_select(col_enquo, col_names)
  names(col_nums) <- NULL
    
  # set nthread to 1 below,  because this may be called by future_map
  out <- data.table::fread(cmd = grepcmd, sep = ",", nThread = 1, 
                           colClasses = "numeric", select = col_nums, ...)

  out <- tibble::as_tibble(out)
    
  return(out)
} # 2}}}
diagnose_var <- function(df){ # {{{2
  # take a data.frame with a single variable, but all .chains and .iterations
  # Calculate Effective Sample Size, Rhat,
  # and Rhats dropping one chain at a time
  MN = nrow(df)
  M = max(df$.chain)
  N <- MN/M
  sample_mat <- df %>% dplyr::arrange(.chain, .iteration) %>% 
    dplyr::pull(.value) %>% matrix(data=., nrow=N, ncol=M)
    drop_k_rhat <- lapply(1:M, function(x) Rhat(sample_mat[,-x]))
    names(drop_k_rhat) <- paste0('Rhat_drop_',1:M)
    return(tibble::tibble(ess_bulk = ess_bulk(sample_mat),
                          ess_tail = ess_tail(sample_mat),
                          Rhat = Rhat(sample_mat)) %>%
      cbind(., as_tibble(drop_k_rhat)))
} # 2}}}
read_and_diagnose <- function(files, par_select, warmup = 0) { # {{{2
  par_expr <- rlang::expr(par_select)
  par_enquo <- rlang::enquo(par_select)
  samples <- tibble(files = files) %>% dplyr::mutate(.chain = 1:n()) %>%
    dplyr::mutate(
      samples = future_map(.x = files, 
                           ~fread_stan(.x, col_select = all_of(!!par_enquo)))
      ) %>%
    dplyr::mutate(samples = map(.x=samples, 
                                ~mutate(.x, .iteration=row_number()) %>%
      dplyr::filter(.iteration > warmup))) %>%
    dplyr::select(.chain, samples) %>% tidyr::unnest(cols=samples) %>%
    dplyr::mutate(.draw = row_number() ) %>%
    tidyr::pivot_longer(cols = !!par_enquo, names_to = '.variable',
                        values_to = '.value' ) %>% 
    dplyr::group_by(.variable) %>% tidyr::nest() %>%
    dplyr::mutate(diag = purrr::map(.x=data, .f=~diagnose_var(.x))) 

  samples <- samples %>% dplyr::select(.variable, diag) %>%
    tidyr::unnest(cols=diag) %>% dplyr::ungroup()
  return(samples)
} # 2}}}
col_major_sub2ind <- function(i,j,I,J){ # {{{2
  indx <- (j-1)*I + i
  return(indx)
} # 2}}}
get_sampling_params <- function(files) { # {{{2
  if (length(files) == 0) return(tibble::tibble())

  params_df <- lapply(files, parse_samples_file)
  params_df <- dplyr::bind_rows(params_df)
  params_df <- janitor::clean_names(params_df)
    
  params_df <- params_df %>%
    dplyr::mutate(chain_id = readr::parse_number(gsub(".*samples_grw_", "", 
                                                      sample_file))) %>%
    dplyr::mutate_if(is.character, readr::parse_guess)
      
  params_df
} # 2}}}
parse_samples_file <- function(samples_file) { # {{{2
  grep_cmd <- glue::glue("grep '^#' {samples_file}")
  params_raw <- as.data.frame(data.table::fread(cmd = grep_cmd, sep = ",", 
                                                sep2 = " ", fill = TRUE))
  params_vec <- gsub("#|# ", "", params_raw[1:nrow(params_raw), 1])
  params_vec <- stringr::str_trim(params_vec)
  params_equal <- params_vec[grep("=", params_vec)]
  params_split <- stringr::str_split(params_equal, "=")

  params_list <- lapply(params_split, function(x){
                          out_df <- tibble::tibble(x[[2]])
                          names(out_df) <- x[[1]]
                          out_df})

  params_df <- dplyr::bind_cols(params_list) %>%
  dplyr::mutate(
    time_warm_up = readr::parse_number(params_vec[grep("seconds \\(Warm-up\\)", 
                                                       params_vec)]),
    time_sampling = readr::parse_number(params_vec[grep("seconds \\(Sampling\\)",
                                                        params_vec)]),
    time_total = readr::parse_number(params_vec[grep("seconds \\(Total\\)",
                                                     params_vec)])
  )

  params_df
} # 2}}}
plot_run_time <- function(diagnostic_df) { # {{{2
  chain_df <- diagnostic_df %>% dplyr::select(-diag, -files) %>%
    tidyr::unnest(cols = chain_info) %>%
    dplyr::left_join(dplyr::distinct(covidmodeldata::acs_names, state_fips, 
                                     state_name), 
                     by = c("State" = "state_fips"))

    chain_plot <- chain_df %>%
      dplyr::mutate(
        state_name = forcats::fct_reorder(state_name, time_total)) %>%
      ggplot2::ggplot(ggplot2::aes(x = time_total, y = state_name)) +
        ggplot2::geom_boxplot() + ggplot2::theme_minimal() + 
        ggplot2::scale_x_time() + 
        ggplot2::labs(title = NULL, 
                      subtitle = glue::glue(
    "Shortest Running Chain: {hms::hms(min(chain_df$time_total))}",
    "Longest Running Chain:  {hms::hms(max(chain_df$time_total))}",
    "Warm Up: {unique(chain_df$warmup)}\tIterations: {unique(chain_df$iter)}",
    .sep = "\n"
    ),
    x = "Duration",
    y = NULL
    )

  chain_plot
} # 2}}}
proc_county <- function(raw_data, geoid.list){ # {{{2
  #############################################################################
  # Define the geographic hierarch and pull out descriptive data for counties.
  # Create a county_df data_frame with columns
  #   geoid
  #   mygeoid
  #   pop
  #   county_name
  #   state_name
  #   i: a sequenctial county code
  #   j1: a sequential level 1 code
  #   j2: a sequential level 2 code
  if(any(class(raw_data)=='sf')) raw_data <- raw_data %>% 
    sf::st_drop_geometry()

  expected_cols <- c('geoid', 'state_name', 'csa_code', 'csa_title', 
                     'metropolitan_micropolitan_statistical_area',
                     'acs_total_pop_e')
  for(i in expected_cols){
    if (!(i %in% colnames(raw_data))){
      error(sprintf('expected column named %s in raw_data', i))
    }
  }
  # filter out missing counties and create sequential id
  county_df <- raw_data %>% dplyr::filter(geoid %in% geoid.list) %>%
    dplyr::mutate(mygeoid = ifelse(geoid=='11001', 
                                   '24XDC', geoid)) %>% # add DC to Maryland
    dplyr::mutate(mystate = substring(mygeoid,1,2)) %>%
    dplyr::group_by(mystate) %>%
    dplyr::mutate(i = row_number())

  ##################################################
  # Create a metro code within each state:
  # Our metro hierarchy will be a metro/nonmetro, then by csa within metros
  # Csa do include micropolitans, but we'll exclude those from metro
  # First, recode metros to try to delete one county groups

  metro_recode_df <- county_df %>%
    dplyr::select(geoid, mygeoid, mystate, 
                  metropolitan_micropolitan_statistical_area, csa_code, 
                  csa_title) %>%
    dplyr::mutate(
      metro = forcats::fct_explicit_na(metropolitan_micropolitan_statistical_area) == 'Metropolitan Statistical Area') %>%
      # Change Connecticut, whch has one non-metro county (Litchfield) to metro.
      dplyr::mutate(metro = if_else(mygeoid == '09005', TRUE, metro)) %>%
      dplyr::mutate(my_metro_code= as.character(
        forcats::fct_explicit_na(csa_code, na_level='998')))  %>%
      dplyr::mutate(my_metro_title = if_else(is.na(csa_title), 
                                             'not_csa', csa_title)) %>% 
      dplyr::mutate(my_metro_code = ifelse(metro, my_metro_code, '999')) %>%
      dplyr::mutate(
        my_metro_title = if_else(metro, my_metro_title, 'not_metro')) %>%
      # Calculate number of counties in each group
      dplyr::group_by(mystate, my_metro_code, my_metro_title) %>%
      dplyr::mutate(num_counties =  dplyr::n()) %>%
      dplyr::ungroup()   %>%
      # Recode any csas with only one or two counties
      dplyr::mutate(
        my_metro_code = ifelse(metro & num_counties <= 2, '998', 
                               my_metro_code),
        my_metro_title = ifelse(metro & num_counties <= 2, 'not_csa', 
                                my_metro_title),) %>%
      # ReCalculate number of counties in each group
      dplyr::group_by(mystate, my_metro_code, my_metro_title) %>%
      dplyr::mutate(num_counties =  dplyr::n()) 

  # Create a dataset with one row per metro area, and a unique id 'j'
  metro_j_df <- metro_recode_df %>%
    # Summarize csas
    dplyr::group_by(mystate, my_metro_title, my_metro_code) %>%
    dplyr::summarize(n=dplyr::n()) %>%
    # give index values
    dplyr::arrange(mystate, my_metro_code) %>%
    dplyr::group_by(mystate) %>%
    dplyr::mutate(j = dplyr::row_number()) 

  # Add j back to the county database
  metro_recode_df <- metro_recode_df %>%
    dplyr::left_join(metro_j_df, 
                     by = c('mystate', 'my_metro_code', 'my_metro_title')) %>%
    dplyr::select(geoid, mygeoid, mystate, metro, my_metro_code,
                  my_metro_title, j) %>%
    # set j=0 if non-metro
    dplyr::mutate(group1 = ifelse(metro, j, 0),
  group2 = ifelse(metro, 1, 2),
  group2_name = ifelse(metro, 'Metropolitan', 'Non-Metropolitan'),
  group1_name = my_metro_title)

  county_df <- county_df %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metro_recode_df, by=c('geoid','mygeoid','mystate')) %>%
    dplyr::select(geoid, mygeoid, mystate, state_name, county_name, i, group1, 
                  group2, group1_name, group2_name, pop=acs_total_pop_e)

  return(county_df)
} # 2}}}
proc_covid <- function(raw_data, DATE_0, ZERO_PAD){ # {{{2
  #library(lubridate)
  #library(covidmodeldata)
  #############################################################################
  #############################################################################
  # Create a covid_df data_frame with columns
  #   geoid
  #   date
  #   Y (new_cases)
  #   new_cases
  # The data should extend from ZERO_PAD (7?) days before the first record to 
  # the current time.
  
  if(!('date' %in% colnames(raw_data))) stop(
    'Expected column called date in raw_date')
    if(!('geoid' %in% colnames(raw_data))) stop(
      'Expected column called geoid in raw_date')
    if(!('total_cases' %in% colnames(raw_data))) stop(
      'Expected column called total_cases in raw_date')

  raw_data <- raw_data %>% dplyr::filter(!is.na(geoid))

  LAST_DAY <- max(raw_data$date)
  GEOIDS <- unique(raw_data$geoid)

  # Create a data_frame with every day, county, and cases
  covid_df <- tidyr::crossing(geoid = unique(raw_data$geoid),
                              date  = seq.Date(from=as.Date(DATE_0), 
                                               to=max(raw_data$date), 
                                               by=1)) %>%
    dplyr::left_join(nyt_data %>% dplyr::select(geoid, date, total_cases))


  # Filter to only include days after ZERO_PAD
  covid_df <- covid_df %>% dplyr::group_by(geoid) %>%
    dplyr::mutate(first_day = min(date[!is.na(total_cases)])) %>%
    dplyr::filter(date >= first_day - ZERO_PAD) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(geoid, date) 

  #########################################################
  # Remove the negative cases;
  # the way I will do this is by making sure  that the confirmed cases is 
  # monotonic, then recomputing new cases as difference

  covid_df <- covid_df %>% dplyr::group_by(geoid) %>%
    dplyr::arrange(geoid, dplyr::desc(date)) %>%
    dplyr::mutate(future_min = cummin(total_cases)) %>%
    dplyr::arrange(geoid, date) %>%
    dplyr::mutate(
                  future_min = ifelse(date < first_day, 0, future_min),
                  Y = future_min - lag(future_min, n = 1L, default = 0)) %>%
    dplyr::ungroup() %>% dplyr::arrange(geoid, date) %>%
    dplyr::select(-first_day, -future_min) %>%
    dplyr::mutate(t = as.integer(date-as.Date(DATE_0) + 1))
  return(covid_df)
} # 2}}}
remove_prisons <- function(covid_df, county_df){ # {{{2
  expected_cols <- c('geoid', 'date', 'Y')
  for(i in expected_cols){
    if (!(i %in% colnames(covid_df))) stop(sprintf(
      'expected column named %s in covid_df', i))
    }
    
  expected_cols <- c('geoid', 'pop')
  # FIXME: this logic should be all(cols %in% colnames(county_df))
  for(i in expected_cols){
    if (!(i %in% colnames(county_df))) stop(sprintf(
      'expected column named %s in county_df', i))
  }

  # Filter for huge jumps that are likely prisons
  covid_df <- covid_df %>% 
    dplyr::left_join(county_df %>% dplyr::select(geoid, pop), by='geoid') %>%
    dplyr::arrange(geoid, date) %>%
    dplyr::mutate(Yold = Y)

  MAX <- 0
  for(i in 2:nrow(covid_df)){
    if(covid_df$geoid[i]!=covid_df$geoid[i-1]) MAX = 0
    if(is.na(covid_df$Y[i])) {next}
    if(covid_df$pop[i] > 30000) {next}
    if (covid_df$Y[i] > 20 & 
      sum(covid_df$Yold[(i-1):(i+1)], na.rm=TRUE) > 50 & 
      sum(covid_df$Yold[(i-1):(i+1)], na.rm=TRUE) > 8*MAX
    ) covid_df$Y[i] <- NA
    if (!is.na(covid_df$Y[i]) & covid_df$Y[i]>MAX) MAX = covid_df$Y[i]
  }
  return(covid_df %>% select(-Yold))
} # 2}}}
retrieve_nyt_data <- function(level = c('county', 'state')) { # {{{2
  if (level == 'county') url <-
      "https://github.com/nytimes/covid-19-data/raw/master/us-counties.csv"
  else if (level == "state") url <- 
      "https://github.com/nytimes/covid-19-data/raw/master/us-states.csv"
  dat <- data.table::fread(url)
  data.table::setnames(dat, names(dat), snakecase::to_snake_case(names(dat)))
  dat
} # 2}}}
# 1}}}
# Ingest the data -------------------------------------------------------- {{{1
nyt_county <- retrieve_nyt_data('county')
nyt_state <- retrieve_nyt_data('state')
msa_list <- read_excel("data-raw/msa-list.xls")

nyt_county_datemax <- max(nyt_county[['date']])
nyt_state_datemax <- max(nyt_state[['date']])
acs_data <- retrieve_acs_data()
# 1}}}
# create the model session information ----------------------------------- {{{1
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
  NTHIN = 2,
  NCHAINS = 8,
  #WARMUP = (2/8)/2,
  WARMUP = 500, # WARMUP is number of retained (before thinning) samples
  DATA_DIR = file.path('/covidmodeldata', nyt_county_datemax),
  DIAG_DF_LOC = file.path('/covidmodeldata', 'diagnostic.Rdata'),
  ZERO_PAD = 7, # Number of days of zeros to add before first obs for each county
  # Number of days forward to predict (from last data day, not from Sys.Date())
  TPRED = 7,
  # My logic here is we want fewer than one peak per week.
  SPL_K = floor(as.numeric(as.Date(nyt_county_datemax) - 
                             as.Date('2020-03-01')) / 14) 
  
) # 1}}}
# data prep -------------------------------------------------------------- {{{1
library(covidmodeldata)
# folder with date, directory by fips code (names must stay the same)
# 
if (model_session_info$CLEAN_DIR) {
  unlink(model_session_info$DATA_DIR, recursive=TRUE)
  print(glue::glue("{DATA_DIR} removed before rewriting"))
}

if (is.null(model_session_info$NYT_FILE)) {
  nyt_data <- retrieve_nyt_data('county') %>% 
    format_nyt(distribute_unknowns = FALSE)  # Assign KC, but that's it
} else nyt_data <- data.table::fread(model_session_info$NYT_FILE)

# tibble with cols: geoid, date, total_cases, Y, t
covid_df <- proc_covid(raw_data = nyt_data %>% 
                       dplyr::select(geoid, date, total_cases=total_cases_mdl),
                       DATE_0 = model_session_info$DATE_0,
                       ZERO_PAD = model_session_info$ZERO_PAD) %>%
  dplyr::left_join(nyt_data %>% dplyr::select(geoid, date, new_cases)) %>%
  dplyr::filter(date <= min(lubridate::as_date(model_session_info$DATE_N), 
                            lubridate::as_date(model_session_info$DATE)))

if (is.null(model_session_info$ACS_FILE)) {
  acs_data <- covidmodeldata::acs_data
} else acs_data <- data.table::fread(model_session_info$ACS_FILE)

# tibble with cols: geoid, mystate, i, group1, group2, pop
county_df <- proc_county(raw_data = acs_data,
                                         geoid.list = unique(covid_df$geoid))
# remove the prisons
covid_df <- remove_prisons(covid_df, county_df) %>%
  dplyr::filter(!is.na(Y))

# create the date tibble
date_df <- tibble::tibble(
  date = seq.Date(from = as.Date(model_session_info$DATE_0),
                  t = max(covid_df$date) + model_session_info$TPRED,
                  by = 1)) %>%
  dplyr::mutate(t = as.integer(date - as.Date(model_session_info$DATE_0) + 1))

# create the time splines
X_dow <- date_df %>% dplyr::mutate(weekdays(date), 
                                   levels = c('Monday', 'Tuesday', 'Wednesday',
                                              'Thursday', 'Friday', 'Saturday',
                                              'Sunday')) %>%
  dplyr::arrange(t) %>% model.matrix(~dow - 1, data = .)

# create the spline basis for the timeseries
# with one knot per date, spline is B * alpha = I * alpha
# the coefficients alpha have a random walk
# KEY OUTPUTS:
#   Z, the spline
#   krig_wt
#   krig_var_chal
T <- as.integer(max(covid_df$date) - as.date(model_session_info$DATE_0) + 1)
L <- pracma::tril(toeplitz(c(1, -2, 1, 
                             rep(0, T + model_session_info$TPRED - 4))))
# create a spline over observable time
svd <- svd(L[1:(T - 1), 1:(T - 1)])
Z <- svd$v %*% diag(1 / svd$d)
# Flip order so most important is first:
Z <- Z[,rev(1:ncol(Z))]
ZZ <- Z  # TODO is this needed?
# add a low-rank-approximation
Z <- z[,1:model_session_info$SPL_K]
# now add the first time to it
Z <- rbind(0, Z)
# calculate the kriging weight and kriging variance
Cov <- tcrossprod(solve(L))
C11 <- Cov[1:(T - 1), 1:(T-1)]
C22 <- Cov[T:(T + model_session_info$TPRED - 1), 
           T:(T + model_session_info$TPRED - 1)]
C12 <- Cov[1:(T - 1), T:(T + model_session_info$TPRED - 1)]
krig_wt <- zapsmall(solve(C11, C12))
krig_var_chol <- t(zapsmall(chol(C22 - t(C12) %*% krig_wt)))
krig_wt <- rbind(0, krig_wt)

# PRIOR on variabce as spline
# - expect the value to grow by 3 orders of magnitude, over 90 days
max_var <- Cov[T, T]  # variance on last day
# sd that should give us log(10^4) change
sd_scale <- log(10^2) / sqrt(max_var)

# save out the data created thus far
if (!dir.exists(model_session_info$DATA_DIR)) {
  dir.create(model_session_info$DATA_DIR, recursive = TRUE)
}
save(covid_df, county_df, X_dow, Z, date_df,
     file = file.path(model_session_info$DATA_DIR, 'data_frames.Rdata'))

# get a list of states to loop through
state_list <- county_df %>% dplyr::group_by(mystate) %>% 
  dplyr::distinct(county_name) %>% 
  dplyr::summarize(n = dplyr::n()) %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::pull(mystate)

# create a list with one entry per chain (50 * NITER)
stan_fit_list <- vector(mode = 'list', 
                        length = length(state_list) *
                          model_session_info$NCHAINS)
slot = 1  # a counter
for (i in 1:length(state_list)) {
  STATE_FIPS = state_list[i]
  # create the sample directory
  STATE_SAMPLES_DIR <- file.path(model_session_info$SAMPLES_ROOT,
                                 model_session_info$DATE,
                                 STATE_FIPS)
  if (!dir.exists(STATE_SAMPLES_DIR)) dir.create(STATE_SAMPLES_DIR,
                                                 recursive = TRUE)
  # pull out geo data for this state
  this_county_df <- county_df %>% dplyr::filter(mystate == STATE_FIPS) %>%
    dplyr::arrange(i)

  # pull out covid data for this state
  this_covid_df <- covid_df %>%
    dplyr::filter(geoid %in% this_county_df$geoid) %>%
    dplyr::left_join(this_county_df %>% 
                     dplyr::select(geoid, i, group1, group2), by = 'geoid')

  # create the county-level regression data matrix
  Xdf <- this_county_df %>% dplyr::arrange(i) %>%
    dplyr::select(i, group1, group2, pop) %>%
    dplyr::mutate(pop_s = as.vector(scale(pop)),
                  lpop_s = as.vector(scale(log(pop))))
  group1 <- if (max(this_county_df$group1) > 1) Xdf else rep(0, nrow(Xdf))
  group1_bin <- which(group1 > 0)
  group1_bin_id <- group1[group1_bin]
  stan_data <- list(Pop = Xdf %>% dplyr::pull(pop),
                    Y = this_covid_df %>% dplyr::pull(Y), # use modeled cases
                    Y_i = this_covid_df %>% dplyr::pull(i),
                    Y_t = this_covid_df %>% dplyr::pull(t),
                    Y_lin_idx = col_major_sub2ind(i = this_covid_df %>%
                                                    dplyr::pull(i),
                                                  j = this_covid_df %>%
                                                    dplyr::pull(t),
                                                  I = max(this_covid_df$i),
                                                  J = max(this_covid_df$t)),
                    I = max(this_covid_df$i),
                    T = max(covid_df$t),
                    TPRED = model_session_info$TPRED,
                    N = nrow(this_covid_df),
                    Ki = 1,
                    Kt = 7,
                    Kit = 0,
                    Xi = Xdf %>% dplyr::select(lpop_s),
                    Xt = X_dow,
                    Xit = array(0, dim = c(
                      nrow(this_covid_df) * (max(covid_df$t) + 7), 0)),
                    krig_wt = krig_wt,
                    krig_var_chol = krig_var_chol,
                    spl_K = model_session_info$SPL_K,
                    Z_spl = Z,
                    J1 = if (max(this_county_df$group1) > 1) max(Xdf$group1) else 0,
                    J2 = max(Xdf$group2),
                    group1 = group1,
                    N1 = length(group1_bin),
                    group1_bin = group1_bin,
                    group1_bin_id = group1_bin_id,
                    group2 = Xdf %>% dplyr::pull(group2),
                    taub0_scale = .5 * sd_scale,
                    taub1_scale = .5 * sd_scale,
                    taub2_scale = sd_scale,
                    sample_flag = TRUE,
                    lppd_flag = FALSE,
                    post_pred = FALSE)

  for (j in 1:model_session_info$NCHAINS) {
    sample_file_name <- file.path(STATE_SAMPLES_DIR, paste0('samples_grw_', j,
                                                            '.csv'))
    diagnostic_file_name <- file.path(STATE_SAMPLES_DIR, 
                                      paste0('diagnostic_grw_', j, '.csv'))
    stan_fit_list[[slot]] <- list(stan_data = stan_data,
                                  sample_file = sample_file_name,
                                  diagnostic_file = diagnostic_file_name,
                                  Xdf = Xdf,
                                  covid_df = this_covid_df)
    slot = slot + 1
  }
  stan_dat <- stan_fit_list[[slot - 1]]
  save(stan_dat,
       file = file.path(dirname(stan_fit_list[[slot - 1]]$sample_file),
                        'standata.RData'))
}

save(stan_fit_list, file=file.path(model_session_info$DATA_DIR,
                                   'stan_fit_list.Rdata'))
# 1}}}
# run model ------------------------------------------------------------ # {{{1
worker_ips <- model_session_data$cades_workers$ip_address
# Jesse says these ips are broken or wouldn't connect
broke_ips <- c('172.22.3.211', '172.22.3.204', '172.22.3.203')
worker_ips <- setdiff(worker_ips, broke_ips)
# connect and create a cluster
cl <- future::makeClusterPSOCK(
  workers = worker_ips,
  user = 'root',
  rshopts = c('-o', 'StrictHostKeyChecking=no', '-o', 'IdentitiesOnly=yes',
              '-i', model_session_info$SSH_KEY),
  rscript = c('sudo', 'docker', 'run', '--net=host', '-v', 
              '/covidmodeldata:/covidmodeldata', 'rstan/covid', 'Rscript'),
  dryrun = FALSE,
  verbose = TRUE
)
future::plan(list(future::tweak(cluster, workers = cl), future::multiprocess))
# split stan fit list
stan_fit_list_split <- parallel::clusterSplit(cl, stan_fit_list)
# compile the stan model
stan_mod <- stan_model(file='analysis/nb1_spline.stan')
# initialization function
init_fun <- function(chain_id) list(a = runif(1, -13, 9),
                                    tau_a0 = runif(1, 0, .02),
                                    tau_a1 = runif(1, 0, .02),
                                    tau_a2 = runif(1, 0, 2),
                                    raw_tau_splb0 = runif(1, 0, .02),
                                    raw_tau_splb1 = runif(1, 0, .02),
                                    raw_tau_splb2 = runif(1, 0, 2),
                                    phi = runif(1, .1, .9))
# run the model
r <- future::future_map(
  # cluster of workers
  .x = stan_fit_list_split,
  .f = ~ {
    outer_idx <- .x
    future::future_map(
      # multiprocess within each worker
      .x = 1:length(outer_idx),
      .f = ~ {
        inner_idx <- .x
        stan_fit <- try(rstan::sampling(
          object = stan_mod,
          data = outer_idx[[inner_idx]]$stan_data,
          iter = model_session_info$NITER,
          thin = model_session_info$NTHIN,
          chains = 1,
          cores = 1,
          warmup = model_session_info$WARMUP,
          sample_file = outer_idx[[inner_idx]]$sample_file,
          append_samples = FALSE,
          init = init_fun,
          control = list(adapt_delta = 0.90, max_treedepth = 12)))
        readr::write_lines(x = outer_idx[[inner_idx]]$sample_file,
                           file.path(model_session_info$DATA_DIR,
                                     'chain_checkout.txt'), append = TRUE)
      }
    )
  }
)
# close the cluster
future::plan(future::sequential)
# 1}}}
