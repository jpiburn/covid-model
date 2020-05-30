
################################################################################
# Create a covid_df data_frame with columns
#   geoid
#   date
#   Y (new_cases)
#   new_cases

# The data should extend from ZERO_PAD (7?) days before the first record to the current time.


proc_covid <- function(raw_data, DATE_0, ZERO_PAD){
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

remove_prisons <- function(covid_df, county_df){
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


