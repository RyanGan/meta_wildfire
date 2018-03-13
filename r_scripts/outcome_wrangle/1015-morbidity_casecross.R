# ------------------------------------------------------------------------------
# Title: Creation of morbidity case-crossover dataframe for CO, WA, OR 2010-2015
# Author: Ryan Gan
# Date Created: 2018-03-12
# ------------------------------------------------------------------------------

# This script combines estimated population-weighted PM2.5, temperature, HMS,
# and air quality index data to create a time series from 2010 to 2015.
# I am considering extending this to 2006 to 2010; I would need to process
# temperature data.

# library ----
library(tidyverse)
library(lubridate)
library(case.crossover)

# read ER hosp
er_hosp <- read_csv("./data/health/1015-er_hosp.csv")

# apply time stratified case-crossover function ----
# time-stratified case-crossover ----
start_time <- Sys.time()
casecross_list <- icd9_outcomes %>% 
  # create lists of outcome dataframes
  map(~filter(er_hosp, dx1 %in% .)) %>% 
  # apply ts casecross over function and join with pm data
  map(~casecross(data = ., id = "id", date = "date", period = "month", 
    covariate = c("stay_length", "dx1", "age", "age_cat", "sex", "zip", 
                  "fips", "state"))) 
stop_time <- Sys.time()
time <- stop_time - start_time
print(time)

# head casecross list
head(casecross_list[1])

# saving casecross over list as Rdata
save(casecross_list, file = paste0("./data/health/",
  "1015-morbidity_casecross_list.RData"))

