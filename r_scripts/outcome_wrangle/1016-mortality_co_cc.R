# ------------------------------------------------------------------------------
# Title: Creation of mortaility case-crossover dataframe for Colorado 2010-2015
# Author: Ryan Gan
# Date Created: 2018-03-12
# ------------------------------------------------------------------------------

# This script combines estimated population-weighted PM2.5, temperature, HMS,
# and air quality index data to create a time series from 2010 to 2015.
# I am considering extending this to 2006 to 2010; I would need to process
# temperature data.

# libraries -----
library(tidyverse)
library(lubridate)
library(case.crossover)

# read mortality data in
read_path <- paste0("../colorado_wildfire/data/health/co_death_2016.csv") 

co_mortality_2010_2016 <- read_csv(read_path) %>% 
  mutate(date_of_death = as.Date(dod, format = "%m%d%Y"),
         # race categories
         race_cat = as.factor(
           case_when(race == 0 ~ "other_asian_pacific",
                     race == 1 ~ "white",
                     race == 2 ~ "black", 
                     race == 3 ~ "american_indian",
                     race == 4 ~ "chinese",
                     race == 5 ~ "japanese",
                     race == 6 ~ "hawaiian",
                     race == 7 ~ "other_nonwhite",
                     race == 8 ~ "filipino",
                     race == 9 ~ "unknown")),
         hisp_cat = as.factor(
           case_when(origin == 0 ~ "non_hispanic",
                     origin == 1 ~ "mexican",
                     origin == 2 ~ "puerto_rican",
                     origin == 3 ~ "cuban",
                     origin == 4 ~ "central_south_american",
                     origin == 5 ~ "other_spanish",
                     origin == 9 ~ "unknown")),
         fips = paste0("08", stringr::str_sub(paste0("00",coor), start = -3)),
         # find cause of death respiratory
         resp_ind = case_when(str_sub(ucod, 1, 1) == "J" |
            str_sub(acme1, 1, 1) == "J" | str_sub(acme3, 1, 1) == "J" | 
            str_sub(acme4, 1, 1) == "J" | str_sub(acme4, 1, 1) == "J" |
            str_sub(acme5, 1, 1) == "J" | str_sub(acme6, 1, 1) == "J" | 
            str_sub(acme7, 1, 1) == "J" | str_sub(acme8, 1, 1) == "J" |
            str_sub(acme9, 1, 1) == "J" | str_sub(acme10, 1, 1) == "J" |
            str_sub(acme11, 1, 1) == "J" ~ 1),
         cod_resp = ifelse(is.na(resp_ind), 0, 1),
         ucod_rep = ifelse(str_sub(ucod, 1, 1) == "J", 1, 0),
         # find cause of death cvd
         cvd_ind = case_when(str_sub(ucod, 1, 1) == "I" |
            str_sub(acme1, 1, 1) == "I" | str_sub(acme3, 1, 1) == "I" | 
            str_sub(acme4, 1, 1) == "I" | str_sub(acme4, 1, 1) == "I" |
            str_sub(acme5, 1, 1) == "I" | str_sub(acme6, 1, 1) == "I" | 
            str_sub(acme7, 1, 1) == "I" | str_sub(acme8, 1, 1) == "I" |
            str_sub(acme9, 1, 1) == "I" | str_sub(acme10, 1, 1) == "I" |
            str_sub(acme11, 1, 1) == "I" ~ 1),
         cod_cvd = ifelse(is.na(cvd_ind), 0, 1),
         ucod_cvd = ifelse(str_sub(ucod, 1, 1) == "I", 1, 0))

# I'm only going to create dataframes for respiratory and cardiovascular major
# underlying cause of death ICD10 codes since I think Rish presents data this way
# Note, I could look at subsets of outcomes
copd <- co_mortality_2010_2016 %>% 
  filter(str_sub(ucod, 1, 3)== "J44")
# time-stratified case-crossover ----
start_time <- Sys.time()
casecross_list <- c("J", "I") %>% 
  # create lists of outcome dataframes
  map(~filter(co_mortality_2010_2016, str_sub(ucod, 1, 1) %in% .)) %>% 
  # apply ts casecross over function and join with pm data
  map(~casecross(data = ., id = "vsid", date = "date_of_death", period = "month", 
    covariate = c("ucod", "age", "sex", "race_cat", "hisp_cat", "zip", 
                  "fips", "wrfgrid_id"))) 
stop_time <- Sys.time()
time <- stop_time - start_time
print(time)
