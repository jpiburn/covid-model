
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