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
      samples = map(.x = files, ~fread_stan(.x, col_select = all_of(!!par_enquo)))
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
                        .f=~diagnose_var(.x))) %>%
    select(
      .variable, 
      diag) %>%
    unnest(
      cols=diag) %>%
    ungroup()
}

