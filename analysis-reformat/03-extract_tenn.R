library(rstan)
library(tidyverse)
library(shinystan)
DATA_DIR <- '/data/covid/tmp/TN/2020-05-31/'

source('analysis-reformat/00-functions.R')
DATE <- '2020-05-31'
source(sprintf('%s/00-PARAMS.R', DATA_DIR))




# Load in the standata
load(file.path(DATA_DIR, '47', 'standata.RData'))
# Load in ancillary data (shapefiles and what not)
load(file=file.path(DATA_DIR, 'data_frames.Rdata'))
acs_data <- covidmodeldata::acs_data %>% sf::st_drop_geometry()


# Load draws
chain_files <- list.files(file.path(DATA_DIR, '47'), pattern='*.csv',full.names = TRUE)
draws <- read_stan_draws(chain_files,  par_select=c(everything()), warmup=WARMUP/NTHIN)


############################################################
# Create county summaries
draws_lambda <- 
  draws %>%
  filter(
    str_detect(.variable, 'log_lambda')) %>%
  separate(
    .variable, 
    sep='\\.', 
    into=c('variable', 'i','t')) %>%
  mutate(
    state_name='Tennessee',
    i = as.integer(i), 
    t=as.integer(t)
  )   %>%
  left_join(
    county_df %>%
      select(state_name, i, county_name, geoid, pop, group1,group2, group1_name, group2_name),
    by=c('state_name','i')
  ) %>%
  left_join(
    date_df %>%
      select(
        t,
        date
      )
  ) %>%
  left_join(
    covid_df %>% select(geoid, date, new_cases),
    by=c('geoid','date')
  )


draws_Y <- 
  draws %>%
  filter(
    str_detect(.variable, 'Y_sim')) %>%
  separate(
    .variable, 
    sep='\\.', 
    into=c('variable', 'i','t')) %>%
  mutate(
    state_name='Tennessee',
    i = as.integer(i), 
    t=as.integer(t)
  )   %>%
  left_join(
    county_df %>%
      select(state_name, i, county_name, geoid, pop, group1,group2, group1_name, group2_name),
    by=c('state_name','i')
  ) %>%
  left_join(
    date_df %>%
      select(
        t,
        date
      )
  ) %>%
  left_join(
    covid_df %>% select(geoid, date, new_cases),
    by=c('geoid','date')
  )

county_lambda <- 
  draws_lambda %>%
  group_by(
    state_name,
    variable, 
    county_name,
    geoid,
    date
    ) %>%
  summarize(
    pop = pop[1],
    q50 = exp(quantile(.value, .5)),
    q05 = exp(quantile(.value, .05)),
    q15 = exp(quantile(.value, .15)),
    q25 = exp(quantile(.value, .25)),
    q35 = exp(quantile(.value, .35)),
    q45 = exp(quantile(.value, .45)),
    q55 = exp(quantile(.value, .55)),
    q65 = exp(quantile(.value, .65)),
    q75 = exp(quantile(.value, .75)),
    q85 = exp(quantile(.value, .85)),
    q95 = exp(quantile(.value, .95)),
    new_cases = new_cases[1]) %>%
  ungroup() %>%
  mutate(q50 = q50*pop,
         q05 = q05*pop,
         q15 = q15*pop,
         q25 = q25*pop,
         q35 = q35*pop,
         q45 = q45*pop,
         q55 = q55*pop,
         q65 = q65*pop,
         q75 = q75*pop,
         q85 = q85*pop,
         q95 = q95*pop)

# Merge with case data and set to NA before first period
county_first_case <- county_lambda %>%
  filter(!is.na(new_cases)) %>%
  group_by(geoid) %>%
  summarize(
    first_case = min(date)
  ) %>%
  ungroup()

county_out <- 
  county_lambda %>%
  left_join(county_first_case) %>%
  filter(date >= first_case) %>%
  select(-first_case)

alpha = .1
plot_00 <- county_out %>%
  ggplot(aes(x=date)) +
  geom_ribbon(aes(ymin = q05, ymax=q95), alpha=alpha) +
  geom_ribbon(aes(ymin = q15, ymax=q85), alpha=alpha) +
  geom_ribbon(aes(ymin = q25, ymax=q75), alpha=alpha) +
  geom_ribbon(aes(ymin = q35, ymax=q65), alpha=alpha) +
  geom_ribbon(aes(ymin = q45, ymax=q55), alpha=3*alpha) +
  geom_point(aes(y=new_cases), color='red', alpha=.3) +
  geom_line(aes(y=q50)) +
  labs(
    y='Cases', 
    title = 'Actual and Predicted New Cases in Knox County') +
  facet_wrap(~county_name, scales = 'free_y')
  


# Merge with geoid, population and count and set lambda to NA if before the first time period.
#county_out <- 
#  county_lambda %>%
#  mutate(
#    state_name == 'Tennessee') %>%
#  left_join(
#    county_df %>%
#      select(geoid, state_name, county_name, pop)
#  )
  


# Calculate probability of decline over 14 days
trend_14 <- draws_lambda %>%
  filter(county_name=='Knox') %>%
  group_by(county_name,.draw) %>%
  arrange(date) %>%
  mutate(lag14 = pop*(exp(.value) - exp(lag(.value,14)))) %>%
  ungroup() %>%
  filter(date==as.Date(DATE)-1)

ggplot(data=trend_14,aes(x=lag14)) + 
  geom_histogram(aes(y=..density..)) + 
  labs(x = 'Change in Number of Cases over 14 days',
    title='Knox County: Change in expected cases over last 14 days',
        subtitle=sprintf('There is a %s chance that the rate has been declining over this period',
                         as.character(round(mean(trend_14$lag14<0),3))))

write_csv(county_out %>% select(-variable),
            path=file.path('results/test',sprintf('results_TN_county_%s.csv', DATE)))



############################################################
# Summarize by metro
metro_out <- 
  draws_lambda %>%
  group_by(.draw, variable, date, group1_name) %>% # aggregate each draw by group1/date
  mutate(estimate = sum(exp(.value)*pop),
            pop = sum(pop),
            new_cases = sum(new_cases, na.rm=TRUE) ) %>%
  ungroup() %>%
  group_by(group1_name, date) %>%
  summarize(
    pop = pop[1],
    new_cases = new_cases[1],
    q50 = (quantile(estimate, .5)),
    q05 = (quantile(estimate, .05)),
    q15 = (quantile(estimate, .15)),
    q25 = (quantile(estimate, .25)),
    q35 = (quantile(estimate, .35)),
    q45 = (quantile(estimate, .45)),
    q55 = (quantile(estimate, .55)),
    q65 = (quantile(estimate, .65)),
    q75 = (quantile(estimate, .75)),
    q85 = (quantile(estimate, .85)),
    q95 = (quantile(estimate, .95))) %>%
  ungroup()

ggplot(metro_out,
       mapping = aes(x=date)) +
  geom_line(aes(y=q50)) +
  geom_point(aes(y=new_cases), alpha=.3, color='red') +
  facet_wrap(~group1_name, scales='free_y')

write_csv(metro_out ,
          path=file.path('results',sprintf('results_TN_metro_%s.csv', DATE)))



############################################################
# Create county summaries
county_lambda <- 
  draws %>%
  group_by(
    .variable) %>%
  summarize(
    q50 = exp(quantile(.value, .5)),
    q05 = exp(quantile(.value, .05)),
    q95 = exp(quantile(.value, .95))) %>%
  separate(
    .variable, 
    sep='\\.', 
    into=c('variable', 'i','t')) %>%
  filter(
    variable == 'log_lambda')  %>%
  mutate(
    i = as.integer(i), 
    t=as.integer(t),
    variable = 'lambda')






