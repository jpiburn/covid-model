
gen_report_maps <- function(param_list) {
  # load libraries ----------------------------------------------------------
  require(tidyverse)
  require(lubridate)
  require(rstan)
  require(foreach)
  require(doParallel)
  require(covidmodeldata)
  require(future)
  require(furrr)
  require(sf)
  
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
  RESULTS_FILE <- param_list$RESULTS_FILE
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
  load(RESULTS_FILE)
  

  # format data -------------------------------------------------------------
  df_map <-
    data_out %>%
    filter(
      date > (as.Date(DATE) - TPRED),
      date <= (as.Date(DATE) + TPRED)
    ) %>%
    mutate(
      lambda_pop = lambda_q50 * pop,
      lambda_pop = if_else(lambda_pop <= 0, 0, lambda_pop),
      new_cases = if_else(new_cases <= 0, 0, new_cases),
      week = if_else(date <= as.Date(DATE), glue::glue('previous {TPRED} days'), glue::glue('next {TPRED} days'))
    ) %>%
    group_by(geoid, state_name, county_name, week) %>%
    summarise(
      total_new_cases_lambda_pop = round(sum(lambda_pop, na.rm = TRUE)),
      total_new_cases_ysim_q50 = round(sum(Ysim_q50, na.rm = TRUE)),
      total_new_cases_observed = round(sum(new_cases, na.rm = TRUE))
    ) %>%
    left_join(landscan_usa) %>%
    mutate(
      total_new_cases_lambda_pop_per_100k = total_new_cases_lambda_pop / night_pop_100k,
      total_new_cases_ysim_q50_per_100k = total_new_cases_ysim_q50 / night_pop_100k,
      total_new_cases_observed_per_100k = total_new_cases_observed / night_pop_100k,
      total_new_cases_lambda_pop_per_100k_bin = case_when(
        total_new_cases_lambda_pop_per_100k >= 200 ~ "6",
        total_new_cases_lambda_pop_per_100k >= 100 ~ "5",
        total_new_cases_lambda_pop_per_100k >= 50 ~ "4",
        total_new_cases_lambda_pop_per_100k >= 20 ~ "3",
        total_new_cases_lambda_pop_per_100k >= 10 ~ "2",
        total_new_cases_lambda_pop_per_100k >= 1 ~ "1",
        TRUE ~ "0"
      ),
      total_new_cases_observed_per_100k_bin = case_when(
        total_new_cases_observed_per_100k >= 200 ~ "6",
        total_new_cases_observed_per_100k >= 100 ~ "5",
        total_new_cases_observed_per_100k >= 50 ~ "4",
        total_new_cases_observed_per_100k >= 20 ~ "3",
        total_new_cases_observed_per_100k >= 10 ~ "2",
        total_new_cases_observed_per_100k >= 1 ~ "1",
        TRUE ~ "0"
      )
    ) %>%
    ungroup()
  

  # load data need for charts -----------------------------------------------
  geom_df <- covidmodeldata::acs_data
  colors <- c("white", RColorBrewer::brewer.pal(6, "Reds"))
  names(colors)  <- 0:6
  
  state_geoms <- geom_df %>%
    dplyr::mutate(
      state_fips = substr(geoid, 1, 2)
    ) %>%
    dplyr::group_by(
      state_fips
    ) %>%
    dplyr::summarise()
  

  # observed new cases of the previous period ---------------------------------
  obs_plot <-
    df_map %>%
      filter(
        week == glue::glue('previous {TPRED} days')
      ) %>%
      dplyr::right_join(
        geom_df,
        by = "geoid"
      ) %>%
      sf::st_as_sf() %>%
      ggplot() +
      geom_sf(
        aes(fill = total_new_cases_observed_per_100k_bin),
        lwd = 0.05
      ) + 
      geom_sf(
        data = state_geoms,
        lwd = 0.2,
        color = "black",
        fill = NA
      ) + scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank()
      )
  

  # predicted new cases for the next period ---------------------------------
  predict_plot <-
  df_map %>%
    filter(
      week == glue::glue('next {TPRED} days')
    ) %>%
    dplyr::right_join(
      geom_df,
      by = "geoid"
    ) %>%
    sf::st_as_sf() %>%
    ggplot() +
    geom_sf(
      aes(fill = total_new_cases_lambda_pop_per_100k_bin),
      lwd = 0.05
    ) + 
    geom_sf(
      data = state_geoms,
      lwd = 0.2,
      color = "black",
      fill = NA
    ) + scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_blank()
    )
  

  # last periods prediction accuracy ----------------------------------------
  
  acc_first_day <- as.Date(DATE) - (TPRED - 1)
  acc_last_day <- as.Date(DATE)
  
  load(glue::glue('{SAMPLES_ROOT}/results/{acc_last_day - TPRED}/results_{acc_last_day - TPRED}.RData'))
  
  df_pred <- data_out %>%
    filter(
      date >= acc_first_day & date <= acc_last_day 
    ) %>%
    select(
      -total_cases, -new_cases
    ) 
  
  
  df_obs_orig <- covidmodeldata::format_nyt(covidmodeldata::get_nyt(), distribute_unknowns = FALSE)
  
  df_obs <- df_obs_orig %>%
    select(-state_name, -county_name)
  
  
  df <- left_join(df_pred, df_obs) %>%
    ungroup() 
  
  df_acc_daily <-
    df %>%
    mutate(
      new_cases = if_else(new_cases < 0| is.na(new_cases), 0, new_cases),
      estimate = lambda_q50*pop
    ) %>%
    group_by(geoid, state_name, county_name, date) %>%
    # summarise_if(is.numeric, sum) %>%
    mutate(
      ysim_diff = Ysim_q50 - new_cases,
      lambda_diff = lambda_q50 - new_cases/pop,
      within_90 = if_else(new_cases <= Ysim_q95 & new_cases >= Ysim_q05, TRUE, FALSE),
      within_70 = if_else(new_cases <= Ysim_q85 & new_cases >= Ysim_q15, TRUE, FALSE),
      within_50 = if_else(new_cases <= Ysim_q75 & new_cases >= Ysim_q25, TRUE, FALSE),
      under_95 = if_else(new_cases <= Ysim_q95 , TRUE, FALSE),
      under_85 = if_else(new_cases <= Ysim_q85 , TRUE, FALSE),
      under_75 = if_else(new_cases <= Ysim_q75 , TRUE, FALSE),
      prediction_cat = case_when(
        within_50 == TRUE ~ "within_50",
        within_70 == TRUE ~ "within_70",
        within_90 == TRUE ~ "within_90",
        TRUE ~ "outside_90"
      )
    ) %>%
    ungroup() %>%
    left_join(covidmodeldata::acs_data)
  
  
  plot_daily_acc <- 
  df_acc_daily %>%
    sf::st_as_sf() %>%
    ggplot() +
    geom_sf(
      data = covidmodeldata::acs_data,
      fill = "light grey",
      lwd = 0.02
    ) +
    geom_sf(
      aes(fill = prediction_cat),
      lwd = 0.02
    ) +
    scale_fill_manual(
      values = c("within_50" = "#2166ac", "within_70" = "#92c5de", "within_90" = "#d1e5f0", "outside_90" = "#b2182b")
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA)
    ) +
    geom_sf(
      data = state_geoms,
      lwd = 0.2,
      color = "black",
      fill = NA
    ) +
    facet_wrap(~date, nrow = 2)
  
  dir.create(glue::glue('{RESULTS_DIR}/report-plots/'), recursive = TRUE)
  
  ggsave(glue::glue('{RESULTS_DIR}/report-plots/{DATE}-observed-per-capita-caes-previous-{TPRED}-days.png'), plot = obs_plot, height = 8, width = 12, units = "in", dpi = 600)
  ggsave(glue::glue('{RESULTS_DIR}/report-plots/{DATE}-observed-per-capita-caes-previous-{TPRED}-days.png'), plot = predict_plot, height = 8, width = 12, units = "in", dpi = 600)
  ggsave(glue::glue('{RESULTS_DIR}/report-plots/{DATE}-previous-{TPRED}-days-accuracy-map.png'), plot = plot_daily_acc, height = 8, width = 16, units = "in", dpi = 600)

}

