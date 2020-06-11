library(tidyverse)
library(covidmodeldata)
library(yardstick)


source('/covidmodeldata/2020-06-03/00-PARAMS.R')
source('analysis-reformat/00-functions.R')
load("results_2020-06-03.RData")

multi_metric <- metric_set(rmse, mape, mae, huber_loss_pseudo)
multi_metric <- metric_set(rmse)


df_pred <- data_out %>%
  select(
    -total_cases, -new_cases
  ) %>%
  filter(
    date >= as.Date(DATE_0)
  ) %>%
  mutate(
    obs_type = if_else(date <= as.Date(DATE_N), "observed", "predicted")
  )

df_obs <- format_nyt(get_nyt(), distribute_unknowns = FALSE) %>%
  select(-state_fips, -state_name, -county_name)


df <- left_join(df_pred, df_obs)



df_acc <-
  df %>%
  filter(
    !is.na(new_cases),
    obs_type == "predicted"
  ) %>%
  left_join(sf::st_drop_geometry(acs_data)) %>%
  mutate(
    low_pop = if_else(pop <= 100000, TRUE, FALSE),
    metro   = if_else(is.na(cbsa_title), FALSE, TRUE),
    estimate = lambda_q50*pop
  )



df_acc %>%
  group_by(date, state_name) %>%
  summarise(
    new_cases = sum(new_cases, na.rm = TRUE)
  ) %>%
 # multi_metric(new_cases, estimate) %>%
  mutate(
    state_name = forcats::fct_reorder(state_name, new_cases)
  ) %>%
  ggplot(
    aes(x = new_cases, y = state_name)
  ) +
  geom_boxplot() +
  theme_minimal()




df_acc %>%
  # filter(
  #   new_cases >= 0
  # ) %>%
  group_by(metro, state_name, date) %>%
  #multi_metric(new_cases, estimate) %>%
  multi_metric((new_cases/pop) , lambda_q50) %>%
  ungroup() %>%
  mutate(
    state_name = forcats::fct_reorder(state_name, .estimate, .desc = TRUE)
  ) %>%
  ggplot(
    aes(y = .estimate, x = date,color = metro, group = metro)
  ) +
  geom_line(colour = "grey") +
  geom_point() +
  theme_minimal() +
  scale_y_log10() +
  labs(
    title = "Performance Evalutaion of NOWcast Forecasts from June 3, 2020",
    subtitle = "States ordered from Highest to Lowest RMSE",
    y = "RMSE New Case Counts Per Person"
  ) +
  facet_wrap(~state_name)
