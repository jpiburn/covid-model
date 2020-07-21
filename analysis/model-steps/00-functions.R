###############################################################
# FUNCTIONS

# Function to take a state and root directory, and return a data_frame
#  with one column containing all sample csv paths
get_sample_paths <- function(STATE, ROOT) {
  
  SAMP_DIR <- file.path(ROOT, STATE)
  csv.files <- list.files(SAMP_DIR,
                          pattern='samples.*[0-9]\\.csv$',
                          full.names=TRUE)
  
  return(csv.files)
}



# Function to create tidy object from list of files:
read_stan_draws <- function(files, par_select, warmup = 500) {
  
  #Input: files: a dataframe
  #csv_files_col_name: character: names of column with paths to csvs
  # Output: a long dataframe with samples
  
  par_enquo <- rlang::enquo(par_select)
  
  samples <- tibble(files=files) %>%
    mutate(.chain=1:n()) %>%
    mutate(
     # samples = map(.x = files, ~vroom_stan(.x, col_select=!!par_enquo))
      samples = map(.x = files, ~fread_stan(.x, col_select = !!par_enquo))
      #samples = map(.x = files, ~fread_stan(.x, col_select = par_select))
      ) %>%
    mutate(
      samples = map(.x=samples, ~mutate(.x,.iteration=1:n()))
      ) %>%
    unnest(cols=samples) %>%
    mutate(.draw = 1:n()) %>%
    filter(.iteration > warmup) %>%
    select(-files) %>%
    select(
      .draw, 
      .chain, 
      .iteration, 
      everything(),
      ) %>%
    ungroup()  %>% 
    pivot_longer(
      cols = !c(.draw, .chain, .iteration),
      names_to = '.variable', 
      values_to = '.value'
    )
  
  return(samples)
}


# Function to take a dataframe of samples, and a paramater string and 
#  calculate rhat for all matching pars (start_with(par))

calculate_rhat <- function(samples, warmup = 0, par) {
  
  temp <- samples %>%
    filter(.iteration > warmup) %>%
    select(
      .chain, 
      .iteration, 
      starts_with(par)
      ) %>%
    pivot_longer(
      cols = starts_with(par),
      names_to = '.variable', 
      values_to = '.value'
      )
  
  NCHAINS = max(temp$.chain)
  NITER = max(temp$.iteration)
  
  grand_mean <- temp %>% 
    group_by(.variable) %>% 
    summarize(
      gmean = mean(.value)
    )
  
  chain_summary <- temp %>% 
    group_by(
      .chain, 
      .variable
      ) %>% 
    summarize(
      wmean = mean(.value), 
      wvar = var(.value)
      ) %>%
    ungroup() %>%
    left_join(
      grand_mean, by = '.variable'
      ) %>%
    group_by(.variable) %>%
    summarize(
      W = mean(wvar), 
      B = NITER /(NCHAINS-1) * sum((wmean-gmean)^2)
      ) %>%
    mutate(
      V = (1-1/NITER) * W + 1/NITER * B,
      rhat = sqrt(V/W)
      ) %>%
    summarize(
      max_rhat = max(rhat),
      which_max = .variable[which.max(rhat)],
      q99_rhat = quantile(rhat, .99, na.rm=TRUE),
      frac_gt_101 = mean(rhat > 1.01)
      )
  
  
  return(chain_summary)
}


# Function to take file names and read csv and calculate rhat summary
# INPUTS:
#  files: a character vector of stan sample files to process
#  warmup: # number of warmup samples
#  par_select: a tidy select of parameters to select (e.g. c(starts_with('b0_raw')))
#  k: will drop chain k from calculation r-hat (default NULL)
read_and_summarize <- function(files, par_select, warmup = 0,  k = NULL) {
  
  par_expr <- rlang::expr(par_select)
  par_enquo <- rlang::enquo(par_select)
  
  if(length(files) < 3) {
    return(
      tibble(
        max_rhat = NA*0,
        which_max = as(NA,'character'),
        q99_rhat = NA*0,
        frac_gt_101 = NA*0)
      )
  }
  
  samples <- tibble(files = files) %>%
    mutate(
      .chain = 1:n()
      ) %>%
    mutate(
      # samples = map(.x = files, ~vroom_stan(.x, col_select = !!par_enquo))
      samples = map(.x = files, ~fread_stan(.x, col_select = !!par_enquo))
      ) %>%
    mutate(
      samples = map(.x = samples, ~mutate(.x, .iteration = 1:n()))
      ) %>%
    unnest(cols = samples) %>%
    mutate(
      .draw = 1:n()
      ) %>%
    select(
      .draw,
      .chain,
      .iteration,
      everything()
      ) %>%
    ungroup()
  
  
  temp <- samples %>%
    filter(.iteration > warmup) %>%
    select(
      .chain,
      .iteration,
      !!par_enquo
      ) %>%
    pivot_longer(
      cols = !!par_enquo,
      names_to = '.variable',
      values_to = '.value'
      )
  
  NCHAINS = max(temp$.chain)
  NITER = max(temp$.iteration)
  
  if(!is.null(k)) {
    # Delete a chain
    temp <- temp %>% filter(.chain != k)
    NCHAINS = NCHAINS - 1
  }
  
  grand_mean <- temp %>% 
    group_by(.variable) %>% 
    summarize(
      gmean = mean(.value)
      )
  
  chain_summary <- temp %>% 
    group_by(
      .chain,
      .variable
      ) %>% 
    summarize(
      wmean = mean(.value),
      wvar = var(.value)
      ) %>%
    ungroup() %>%
    left_join(grand_mean, by = '.variable') %>%
    group_by(.variable) %>%
    summarize(
      W = mean(wvar),
      B = NITER /(NCHAINS-1) * sum((wmean-gmean)^2)
      ) %>%
    mutate(
      V = (1-1/NITER) * W + 1/NITER * B,
      rhat = sqrt(V/W)
      ) %>%
    summarize(
      max_rhat = max(rhat),
      which_max = .variable[which.max(rhat)],
      q99_rhat = quantile(rhat, .99, na.rm=TRUE),
      frac_gt_101 = mean(rhat > 1.01)
      )
  
  
  return(chain_summary)
}



vroom_stan <- function(file, ...) {
  
  # Stan sample files have commented out lines in the start, middle and end of the file
  # vroom can figure out the start, but not the middle and end
  # Use grep to delete those
  
  tfile <- paste0(file, '.tmp')
  grepcmd <- paste0("grep -vh '^#' ", file, " > ", tfile)
  
  system(grepcmd)
  
  out <- vroom::vroom(tfile, delim = ',', num_threads = 1, 
                      col_types = cols(.default = col_double()), 
                      ...)
  
  unlink(tfile)
  

  return(out)
}


fread_stan <- function(file, col_select, ...) {
  require(tidyselect)
  #col_expr <- rlang::expr(all_of(col_select))
  col_enquo <- rlang::expr(col_select)
  
  # was having trouble getting vroom_stan() to work on windows machine (what the current VM is)
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
                           colClasses = "numeric", select = col_nums, 
                           ...)

  out <- tibble::as_tibble(out)
  
  
  return(out)
}

diagnose_var <- function(df){
  require(rstan)
  # take a data.frame with a single variable, but all .chains and .iterations
  # Calculate Effective Sample Size, Rhat,
  # and Rhats dropping one chain at a time
  MN = nrow(df)
  M = max(df$.chain)
  N <- MN/M
  sample_mat <- df %>% 
    arrange(.chain, .iteration) %>%
    pull(.value) %>%
    matrix(data=., nrow=N, ncol=M)
  drop_k_rhat <- lapply(1:M, 
                        function(x) Rhat(sample_mat[,-x]))
  names(drop_k_rhat) <- paste0('Rhat_drop_',1:M)
  return(
    tibble(
      ess_bulk = ess_bulk(sample_mat),
      ess_tail = ess_tail(sample_mat),
      Rhat = Rhat(sample_mat)) %>%
      cbind(., as_tibble(drop_k_rhat))
  )
}

read_and_diagnose <- function(files, par_select, warmup = 0) {
  require(tidyselect)
  par_expr <- rlang::expr(par_select)
  par_enquo <- rlang::enquo(par_select)
  samples <- tibble(files = files) %>%
    mutate(
      .chain = 1:n()
    ) %>%
    mutate(
      samples = future_map(.x = files, ~fread_stan(.x, col_select = all_of(!!par_enquo)))
    ) %>%
    mutate(
      samples = map(.x=samples, 
                    ~mutate(.x, .iteration=row_number()) %>%
                      filter(.iteration > warmup))
    ) %>%
    select(
      .chain, 
      samples
    ) %>%
    unnest(
      cols=samples) %>%
    mutate(
      .draw = row_number() ) %>%
    pivot_longer(
      cols = !!par_enquo,
      names_to = '.variable',
      values_to = '.value' ) %>% 
    group_by(
      .variable) %>%
    nest()   %>%
    mutate(
      diag = purrr::map(.x=data, 
                        .f=~diagnose_var(.x))) 
  
  samples <- samples %>%
    select(
      .variable, 
      diag) %>%
    unnest(
      cols=diag) %>%
    ungroup()
  return(samples)
}

col_major_sub2ind <- function(i,j,I,J){
  indx <- (j-1)*I + i
  return(indx)
}

get_sampling_params <- function(files) {
  
  if (length(files) == 0) return(tibble::tibble())
  
  params_df <- lapply(files, parse_samples_file)
  params_df <- dplyr::bind_rows(params_df)
  params_df <- janitor::clean_names(params_df)
  
  params_df <- params_df %>%
    dplyr::mutate(
      chain_id = readr::parse_number(gsub(".*samples_grw_", "", sample_file))
    ) %>%
    dplyr::mutate_if(is.character, readr::parse_guess)
  
  params_df
}

parse_samples_file <- function(samples_file) {
  
  grep_cmd <- glue::glue("grep '^#' {samples_file}")
  params_raw <- as.data.frame(data.table::fread(cmd = grep_cmd, sep = ",", sep2 = " ", fill = TRUE))
  
  params_vec <- gsub("#|# ", "", params_raw[1:nrow(params_raw), 1])
  params_vec <- stringr::str_trim(params_vec)
  params_equal <- params_vec[grep("=", params_vec)]
  params_split <- stringr::str_split(params_equal, "=")
  
  params_list <- lapply(params_split, function(x){
    out_df <- tibble::tibble(x[[2]])
    names(out_df) <- x[[1]]
    
    out_df
  })
  
  params_df <- dplyr::bind_cols(params_list) %>%
    dplyr::mutate(
      time_warm_up = readr::parse_number(params_vec[grep("seconds \\(Warm-up\\)", params_vec)]),
      time_sampling = readr::parse_number(params_vec[grep("seconds \\(Sampling\\)", params_vec)]),
      time_total = readr::parse_number(params_vec[grep("seconds \\(Total\\)", params_vec)])
    )
  
  params_df
}

plot_run_time <- function(diagnostic_df) {
  require(ggplot2)
  require(hms)
  
  chain_df <- 
    diagnostic_df %>%
    select(-diag, -files) %>%
    unnest(cols = chain_info) %>%
    left_join(
      distinct(covidmodeldata::acs_names, state_fips, state_name),
      by = c("State" = "state_fips")
    )
  
  chain_plot <- 
    chain_df %>%
    mutate(
      state_name = forcats::fct_reorder(state_name, time_total)
    ) %>%
    ggplot(
      aes(x = time_total, y = state_name)
    ) +
    geom_boxplot() +
    theme_minimal() +
    scale_x_time() +
    labs(
      title = NULL,
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
}



################################################################################
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

proc_county <- function(raw_data, geoid.list){
  if(any(class(raw_data)=='sf')) raw_data <- raw_data %>% sf::st_drop_geometry()
  
  expected_cols <- c('geoid', 'state_name', 'csa_code', 'csa_title', 'metropolitan_micropolitan_statistical_area',
                     'acs_total_pop_e')
  for( i in expected_cols){
    if ( !(i %in% colnames(raw_data))) error(sprintf('expected column named %s in raw_data', i))
  }
  
  # filter out missing counties and create sequential id
  county_df <- 
    raw_data %>%
    filter(geoid %in% geoid.list) %>%
    mutate(mygeoid = ifelse(geoid=='11001', '24XDC', geoid)) %>% # add DC to Maryland
    mutate(mystate = substring(mygeoid,1,2)) %>%
    group_by(mystate) %>%
    mutate(i = row_number())
  
  ##################################################
  # Create a metro code within each state:
  # Our metro hierarchy will be a metro/nonmetro, then by csa within metros
  # Csa do include micropolitans, but we'll exclude those from metro
  # First, recode metros to try to delete one county groups
  
  metro_recode_df <-
    county_df %>%
    select(geoid, mygeoid, mystate, metropolitan_micropolitan_statistical_area, csa_code, csa_title) %>%
    mutate(
      metro = forcats::fct_explicit_na(metropolitan_micropolitan_statistical_area) == 'Metropolitan Statistical Area'
    ) %>%
    # Change Connecticut, whch has one non-metro county (Litchfield) to metro.
    mutate(metro = if_else(mygeoid == '09005', TRUE, metro)) %>%
    mutate(
      my_metro_code= as.character(forcats::fct_explicit_na(csa_code, na_level='998'))
    )  %>%
    mutate(
      my_metro_title = if_else(is.na(csa_title), 'not_csa', csa_title)
    ) %>% 
    mutate(
      my_metro_code = ifelse(metro, my_metro_code, '999')
    ) %>%
    mutate(
      my_metro_title = if_else(metro, my_metro_title, 'not_metro')
    ) %>%
    # Calculate number of counties in each group
    group_by(
      mystate, 
      my_metro_code, 
      my_metro_title) %>%
    mutate(
      num_counties = n()) %>%
    ungroup()   %>%
    # Recode any csas with only one or two counties
    mutate(
      my_metro_code = ifelse(metro & num_counties <= 2, '998', my_metro_code),
      my_metro_title = ifelse(metro & num_counties <= 2, 'not_csa', my_metro_title),
    ) %>%
    # ReCalculate number of counties in each group
    group_by(
      mystate, 
      my_metro_code, 
      my_metro_title) %>%
    mutate(
      num_counties = n()) 
  
  # Create a dataset with one row per metro area, and a unique id 'j'
  metro_j_df <- metro_recode_df %>%
    # Summarize csas
    group_by(
      mystate, 
      my_metro_title, 
      my_metro_code) %>%
    summarize(
      n=n()) %>%
    # give index values
    arrange(mystate, my_metro_code) %>%
    group_by(
      mystate) %>%
    mutate(
      j = row_number()
    ) 
  
  # Add j back to the county database
  metro_recode_df <- metro_recode_df %>%
    left_join(
      metro_j_df,
      by = c('mystate', 'my_metro_code', 'my_metro_title')) %>%
    select(
      geoid, 
      mygeoid, 
      mystate,
      metro,
      my_metro_code,
      my_metro_title,
      j
    ) %>%
    # set j=0 if non-metro
    mutate(group1 = ifelse(metro, j, 0),
           group2 = ifelse(metro, 1, 2),
           group2_name = ifelse(metro, 'Metropolitan', 'Non-Metropolitan'),
           group1_name = my_metro_title)
  
  county_df <-
    county_df %>%
    ungroup() %>%
    left_join(
      metro_recode_df,
      by=c('geoid','mygeoid','mystate')) %>%
    select(geoid, mygeoid, mystate, state_name, county_name, i, group1, group2, group1_name, group2_name, pop=acs_total_pop_e)
  
  return(county_df)
  
}



################################################################################
# Create a covid_df data_frame with columns
#   geoid
#   date
#   Y (new_cases)
#   new_cases

# The data should extend from ZERO_PAD (7?) days before the first record to the current time.


proc_covid <- function(raw_data, nyt_data, DATE_0, ZERO_PAD) {
  require(tidyverse)
  #library(lubridate)
  #library(covidmodeldata)
  ################################################################################
  
  if(! ('date' %in% colnames(raw_data))) stop('Expected column called date in raw_date')
  if(! ('geoid' %in% colnames(raw_data))) stop('Expected column called geoid in raw_date')
  if(! ('total_cases' %in% colnames(raw_data))) stop('Expected column called total_cases in raw_date')
  
  raw_data <- raw_data %>%
    filter(!is.na(geoid))
  
  LAST_DAY = max(raw_data$date)
  GEOIDS <- unique(raw_data$geoid)
  
  # Create a data_frame with every day, county, and cases
  covid_df <- 
    crossing(
      geoid = unique(raw_data$geoid),
      date  = seq.Date(from=as.Date(DATE_0), to=max(raw_data$date), by=1) ) %>%
    left_join(
      nyt_data %>%
        select(geoid,
               date,
               total_cases)
    )
  
  
  # Filter to only include days after ZERO_PAD
  covid_df <- 
    covid_df %>%
    group_by(
      geoid) %>%
    mutate(
      first_day = min(date[!is.na(total_cases)])
    ) %>%
    filter(
      date >= first_day - ZERO_PAD
    ) %>%
    ungroup() %>%
    arrange(
      geoid, 
      date
    ) 
  
  
  #########################################################
  # Remove the negative cases;
  # the way I will do this is by making sure  that the confirmed cases is monotonic, then recomputing new cases as difference
  
  covid_df <- 
    covid_df %>%
    group_by(
      geoid) %>%
    arrange(
      geoid,
      desc(date)) %>%
    mutate(
      future_min = cummin(total_cases)) %>%
    arrange(
      geoid,
      date) %>%
    mutate(
      future_min = ifelse(date < first_day, 0, future_min),
      Y = future_min - lag(future_min, n = 1L, default = 0)
    ) %>%
    ungroup() %>%
    arrange(
      geoid,
      date
    ) %>%
    select(
      -first_day,
      -future_min
    ) %>%
    mutate(t = as.integer(date-as.Date(DATE_0) + 1))
  return(covid_df)
}

remove_prisons <- function(covid_df, county_df) {
  expected_cols <- c('geoid', 'date', 'Y')
  for( i in expected_cols){
    if ( !(i %in% colnames(covid_df))) stop(sprintf('expected column named %s in covid_df', i))
  }
  
  expected_cols <- c('geoid', 'pop')
  for( i in expected_cols){
    if ( !(i %in% colnames(county_df))) stop(sprintf('expected column named %s in county_df', i))
  }
  
  # Filter for huge jumps that are likely prisons
  covid_df <- covid_df %>%
    left_join(
      county_df %>%
        select(geoid, pop),
      by='geoid') %>%
    arrange(geoid, date) %>%
    mutate(Yold = Y)
  
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
  return(
    covid_df %>% select(-Yold)
  )
}



