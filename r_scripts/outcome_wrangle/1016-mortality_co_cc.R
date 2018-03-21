# ------------------------------------------------------------------------------
# Title: Creation of mortaility case-crossover dataframe for Colorado 2010-2015
# Author: Ryan Gan
# Date Created: 2018-03-12
# ------------------------------------------------------------------------------

# This script creates monthly time-stratified case-crossover dataframes for 
# colorado mortality outcomes from 2010 to 2016.
# Note: I may move this to the Colorado wildfire repo.

# libraries -----
library(tidyverse)
library(lubridate)
library(case.crossover)
library(parallel)

# load Rdata vector of outcomes ----
load("./data/health/icd10_outcome.RData")

# read mortality data -----
# changed this path to read on the server
read_path <- paste0("./data/health/co_death_1016.csv") 

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
         ucod_resp = ifelse(str_sub(ucod, 1, 1) == "J", 1, 0),
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

# subsetting copd underlying cause of death to evaulate while process is running
# commenting the following steps out since i don't need to run it each time
# copd <- co_mortality_2010_2016 %>% 
#   filter(str_sub(ucod, 1, 3)== "J44")
# 
# copd_casecross <- casecross(data=copd, id = "vsid", date = "date_of_death", 
#   period = "month", covariate = c("ucod", "age", "sex", "race_cat", 
#                                   "hisp_cat", "zip", "fips", "wrfgrid_id"))
# 
# # write copd casecrossover
# write_csv(copd_casecross, "./data/health/1015-copd_mortality_cc.csv")

# time-stratified case-crossover ----
# parallel distributed lag computing ----
# set up cluster of 8 cores to parallelize models
cores <- parallel::detectCores()
# check to see if cores detected
print(cores)
cl <- makeCluster(cores)
# check cluster
print(cl)

# load packages on each processor of the node/cluster
clusterCall(cl, function() c(library(tidyverse), library(lubridate),
                             library(case.crossover)))

# export new set of objects to global objects to cluster
clusterExport(cl, c("co_mortality_2010_2016", "icd10_outcomes"), 
              envir = .GlobalEnv)

# Distributed lag association for smoke ----------------------------------------
# start time
start <- Sys.time()

# distributed lag function
casecross_list <- parLapply(cl, icd10_outcomes, function(x){
  # filter co_mortality by underlying cause of death in string of icd10 codes
  df <- filter(co_mortality_2010_2016, str_sub(ucod, 1, 3) %in% x)
  # apply custom case.cross function
  df_cc <- casecross(data = df, id = "vsid", date = "date_of_death", 
    period = "month", covariate = c("ucod", "age", "sex", "race_cat", 
                                    "hisp_cat", "zip", "fips", "wrfgrid_id"))
  })

stop <- Sys.time()
time <- stop - start
print(time)

# check first observation of asthma
print(head(casecross_list[[2]]))

# save casecross list----
save(casecross_list, file = "./data/health/co_morality_cc_list.RData")

# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)

# stop cluster
stopCluster(cl)

