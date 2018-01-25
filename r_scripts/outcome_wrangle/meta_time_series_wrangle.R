# ------------------------------------------------------------------------------
# Title: Wrangling meta time-series from state specific time-series dataframes
# Author: Ryan Gan
# Date Created: 5/8/2017
# R Version 3.3.3
# ------------------------------------------------------------------------------

# load libraries ----
library(tidyverse)

# load datasets ----

# oregon health outcomes time series 2013
or_health <- read_csv("./data/or_2013_county_er_time_series.csv") 
# oregon smk time series 
or_smk <- read_csv("./data/or_2013_county_pop_wt_pm.csv") %>% 
  mutate(fips = as.character(fips))
# oregon county pop
or_county_pop <- read_csv("./data/or_census_pop_est_county.csv") %>% 
  # cleaning up data
  # seperate county from state
  separate(GEO.display.label, c("county", "state"), sep = ",") %>% 
  # remove the word 'county' from the county string (mind the space)
  mutate(county = gsub(county, pattern = " County", replacement = ""),
         GEO.id2 = as.character(GEO.id2), state = "OR") %>% 
  # rename variables
  rename(fips = GEO.id2, pop_2010 = respop72010, pop_2011 = respop72011,
         pop_2012 = respop72012, pop_2013 = respop72013, pop_2014 = respop72014,
         pop_2015 = respop72015, pop_2016 = respop72016) %>% 
  # limit to certain variables
  select(state, county, fips, pop_2013) %>% 
  rename(population = pop_2013)


# washington health outcomes time series 2012
wa_health <- read_csv("./data/wa_2012_county_er_time_series.csv")
# washington smk time series
wa_smk <- read_csv("./data/wa_2012_county_pop_wt_pm.csv")
# county population estimates
wa_county_pop <- read_csv("./data/wa_census_pop_est_county.csv") %>% 
    # cleaning up data
  # seperate county from state
  separate(GEO.display.label, c("county", "state"), sep = ",") %>% 
  # remove the word 'county' from the county string (mind the space)
  mutate(county = gsub(county, pattern = " County", replacement = ""),
         GEO.id2 = as.character(GEO.id2), state = "WA") %>% 
  # rename variables
  rename(fips = GEO.id2, pop_2010 = respop72010, pop_2011 = respop72011,
         pop_2012 = respop72012, pop_2013 = respop72013, pop_2014 = respop72014,
         pop_2015 = respop72015, pop_2016 = respop72016) %>% 
  # limit to certain variables
  select(state, county, fips, pop_2012) %>% 
  rename(population = pop_2012)

# data wrangling ----

# check datasets for variable names
glimpse(or_health)
# check county names
unique(or_health$county)
glimpse(or_county_pop)

glimpse(or_smk)
# check county names
unique(or_smk$county)
# check dates
summary(or_smk)

# join oregon smk to oregon health outcomes by date and county
or_ts <- or_health %>% right_join(or_smk, by = c("county", "date")) %>% 
  # some missing values indicating no ER claims on given data for given county
  # setting missing to 0
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .))) %>% 
  # reorganize variables
  select(state, county, fips, date, 3:15, 18:26) %>% 
  # turn . in to a space for county name
  mutate(county = gsub(county, pattern = "\\.", replacement = " ")) %>% 
  # join county population estimate in
  right_join(or_county_pop, by = c("county", "state", "fips")) %>% 
  select(state, county, fips, population, date, 5:26)

# check oregon timeseries 
glimpse(or_ts)
summary(or_ts)
unique(or_ts$county)

# check washington
glimpse(wa_health)
# check smoke
glimpse(wa_smk)
# join washington smk to washington health outcomes by date and county
wa_ts <- wa_health %>% right_join(wa_smk, by = c("county", "date")) %>% 
  # create state variable and add state code to fips code
  mutate(state = "WA", fips = paste0("53", fips)) %>% 
  # some missing values indicating no ER claims on given data for given county
  # setting missing to 0
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .))) %>% 
  # select variables of use
  select(state, county, fips, date, 3:8, 10:16, 19:22, 24:26, 28:29) %>% 
  # join county population estimate in
  right_join(wa_county_pop, by = c("county", "state", "fips")) %>% 
  select(state, county, fips, population, date, 5:26)

# check wash timeseries
glimpse(wa_ts)
summary(wa_ts)

# creating meta time series ----
# variables look the same in both oregon and washington time series dataframes
# bind two datasets together
meta_ts <- bind_rows(wa_ts, or_ts)

summary(meta_ts)

# write permanent dataframe
write_csv(meta_ts, "./data/meta_timeseries.csv")
