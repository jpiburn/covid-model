
library(tidyverse)
library(sf)
library(covidmodeldata)

load("/home/cades/old-stuff/covid-model-master/results_2020-05-28.RData")


df_map <- 
  data_out %>% 
  filter(
    date >= as.Date("2020-05-21")
  ) %>%
  mutate(
    lambda_pop = lambda_q50 * pop,
    lambda_pop = if_else(lambda_pop <= 0, 0, lambda_pop),
    new_cases = if_else(new_cases <= 0, 0, new_cases),
    week = if_else(date <= as.Date("2020-05-27"), "this week", "next week")
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
  ) 

geom_df <- acs_data

colors <- c("white", RColorBrewer::brewer.pal(6, "Reds"))
names(colors)  <- 0:6


#map_plot <-
df_map %>%
  filter(
    week == "next week"
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
  ) + scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    # legend.position="top",
    # legend.direction="horizontal",
    # legend.key.width=unit(2,"cm"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank()
  ) 
ggsave("covid-model-master/next-week.svg", height = 8, width = 12, units = "in", dpi = 800, bg = "transparent")




df_map %>%
  filter(
    week == "this week"
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
  ) + scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    #legend.position="top",
    #legend.direction="horizontal",
    #legend.key.width=unit(2,"cm"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank()
  ) 
ggsave("covid-model-master/this-week.svg", height = 8, width = 12, units = "in", dpi = 800, bg = "transparent")










# increase/decrease map ---------------------------------------------------
df_map %>%
  pivot_wider(
    names_from = week, 
    values_from = total_new_cases_lambda_pop:total_new_cases_observed_per_100k_bin
    ) %>%
  janitor::clean_names() %>%
  mutate(
    increase_cases_per_100k = total_new_cases_lambda_pop_per_100k_next_week - total_new_cases_observed_per_100k_this_week,
    increase_cases = total_new_cases_lambda_pop_next_week - total_new_cases_observed_this_week
  ) %>%
  filter(
    abs(increase_cases) < 1000
  ) %>%
  dplyr::right_join(
    geom_df,
    by = "geoid"
  ) %>%
  sf::st_as_sf() %>%
  ggplot() +
  geom_sf(
    aes(fill = increase_cases),
    lwd = 0.05
  )  + scale_fill_viridis_c()

  theme_minimal() 
  scale_fill_manual(values = colors) +
  theme(
    #legend.position="top",
    #legend.direction="horizontal",
    #legend.key.width=unit(2,"cm"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank()
  ) 



