
library(tidyverse)
library(sf)
library(covidmodeldata)

load("/home/cades/old-stuff/covid-model-master/results_2020-05-28.RData")

last_obs_date <- as.Date("2020-05-27")

df_pred <- data_out %>%
  select(
    -total_cases, -new_cases
  ) %>%
  mutate(
    obs_type = if_else(date <= last_obs_date, "observed", "predicted")
  )

df_obs <- format_nyt(get_nyt(), distribute_unknowns = FALSE) %>%
  select(-state_fips, -state_name, -county_name)


df <- left_join(df_pred, df_obs)

df_min_cases <- df %>%
  filter(
    obs_type == "observed"
  ) %>%
  filter( 
    date == max(date)
  )

counties_max_1 <- df_min_cases[df_min_cases$total_cases <= 1, "geoid"] %>% pull()
counties_max_10 <- df_min_cases[df_min_cases$total_cases <= 10, "geoid"] %>% pull()
counties_max_20 <- df_min_cases[df_min_cases$total_cases <= 20, "geoid"] %>% pull()
counties_max_50 <- df_min_cases[df_min_cases$total_cases <= 50, "geoid"] %>% pull()
counties_max_100 <- df_min_cases[df_min_cases$total_cases <= 100, "geoid"] %>% pull()

counties_min_100 <- df_min_cases[df_min_cases$total_cases >= 100, "geoid"] %>% pull()
counties_min_200 <- df_min_cases[df_min_cases$total_cases >= 200, "geoid"] %>% pull()
counties_min_500 <- df_min_cases[df_min_cases$total_cases >= 500, "geoid"] %>% pull()
counties_min_1000 <- df_min_cases[df_min_cases$total_cases >= 1000, "geoid"] %>% pull()



df_forecast <- df %>%
  left_join(sf::st_drop_geometry(acs_data)) %>%
  mutate(
    low_pop = if_else(pop <= 100000, TRUE, FALSE),
    metro   = if_else(is.na(cbsa_title), FALSE, TRUE),
    counties_max_1 = if_else(geoid %in% counties_max_1, TRUE, FALSE),
    counties_max_10 = if_else(geoid %in% counties_max_10, TRUE, FALSE),
    counties_max_20 = if_else(geoid %in% counties_max_20, TRUE, FALSE),
    counties_max_50 = if_else(geoid %in% counties_max_50, TRUE, FALSE),
    counties_max_100 = if_else(geoid %in% counties_max_100, TRUE, FALSE),
    counties_min_100 = if_else(geoid %in% counties_min_100, TRUE, FALSE),
    counties_min_200 = if_else(geoid %in% counties_min_200, TRUE, FALSE),
    counties_min_500 = if_else(geoid %in% counties_min_500, TRUE, FALSE),
    counties_min_1000 = if_else(geoid %in% counties_min_1000, TRUE, FALSE)
  ) %>%
  filter(!is.na(new_cases))

df_forecast %>% 
  #group_by(metro) %>% 
  filter(
    obs_type == "predicted"
  ) %>%
  summarise(
    int_90 = mean(new_cases >= Ysim_q05 & new_cases <= Ysim_q95), 
    int_70 = mean(new_cases >= Ysim_q15 & new_cases <= Ysim_q85), 
    int_50 = mean(new_cases >= Ysim_q25 & new_cases <= Ysim_q75)
    ) %>% View()



df_forecast %>%
  filter(
    #state_fips == "47",
    date >= as.Date("2020-05-14"),
   pop > 50000,
   state_name == "North Carolina"
  ) %>%
  ggplot(
    aes(x = date, y= new_cases, group = county_name)
  ) +
  geom_vline(xintercept = as.Date("2020-05-27"), size = .8, linetype=2, color = "dark grey") +
  geom_ribbon(aes(ymin = Ysim_q05, ymax = Ysim_q95), alpha = 0.1, fill = "blue") +
  geom_ribbon(aes(ymin = Ysim_q25, ymax = Ysim_q75), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = pop*lambda_q50), color = "blue", size = 1.2) +
  geom_point(size = 1, color = "red") +
  # geom_point(aes(y = Ysim_q50), size = 1, color = "black") +
  theme_minimal() +
  scale_x_date(
    labels = scales::label_date_short(),
    breaks = scales::breaks_width("3 days"),
    expand = c(0,0)
  ) +
  facet_wrap(~county_name, scales = "free")







df_forecast_map <- df_forecast %>% 
  group_by(geoid, state_name, county_name) %>% 
  summarise(
    int_90 = mean(new_cases >= Ysim_q05 & new_cases <= Ysim_q95), 
    int_70 = mean(new_cases >= Ysim_q15 & new_cases <= Ysim_q85), 
    int_50 = mean(new_cases >= Ysim_q25 & new_cases <= Ysim_q75)
  ) %>%

  
  
  # df %>%
  #   ggplot(
  #     aes(x = date, y= new_cases_mdl, group = county_name)
  #   ) +
  #   geom_ribbon(aes(ymin = Ysim_q05, ymax = Ysim_q95), alpha = 0.1, fill = "blue") +
  #   geom_ribbon(aes(ymin = Ysim_q25, ymax = Ysim_q75), alpha = 0.3, fill = "blue") +
  #   geom_line(aes(y = Ysim_q50), color = "blue", size = 1.2) +
  #   geom_point(size = 1, color = "grey") +
#   geom_point(aes(y = observed_new_cases_mdl), size = 2, color = "red") +
#   theme_minimal() +
#   facet_wrap(~county_name, scales = "free")





pdf_map <- 
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
    abs(increase_cases_per_100k) < 150
  ) %>%
  dplyr::right_join(
    geom_df,
    by = "geoid"
  ) %>%
  sf::st_as_sf() %>%
  ggplot() +
  geom_sf(
    aes(fill = increase_cases_per_100k),
    lwd = 0.05
  )  + 
  scale_fill_gradient2() +
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



