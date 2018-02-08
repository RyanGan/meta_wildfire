# ------------------------------------------------------------------------------
# Title: Time-stratified case-crossover dataframes for inpatient emergency 
#        department hospitalizations and results
# Author: Ryan Gan
# Date Created: 2018-02-06
# ------------------------------------------------------------------------------

# I think I'll run this on the server; maybe in parallel

# libraries ----
library(tidyverse)
library(haven)
library(survival)
library(case.crossover)

# creation of county key -----
COUNTYRES <- c(paste0("0", 1:9), as.character(10:39))
# county name
county <- c("Adams", "Asotin", "Benton", "Chelan", "Clallam", "Clark", 
            "Columbia", "Cowlitz", "Douglas", "Ferry", "Franklin", 
            "Garfield", "Grant", "Grays Harbor", "Island", "Jefferson",
            "King", "Kitsap", "Kittitas", "Klickitat", "Lewis", "Lincoln", 
            "Mason", "Okanogan", "Pacific", "Pend Oreille", "Pierce", "San Juan",
            "Skagit", "Skamania", "Snohomish", "Spokane", "Stevens", "Thurston",
            "Wahkiakum", "Walla Walla", "Whatcom", "Whitman", "Yakima")
# fips code
fips <- c(paste0("5300", seq(1, 9, by = 2)), paste0("530", seq(11, 77, by =2)))

# county key
county_key <- as_data_frame(cbind(COUNTYRES, fips, county))

# load icd9 outcomes data ----
load("./data/health/icd9_outcome_vectors.RData")

# read in pm 2.5 data
county_pm <- read_csv("./data/smoke/2015-smoke_wus_county_popwt.csv")

# may be worth creating in a function one day
gwr_lag <- county_pm %>% 
  select(fips, date, gwr_smk) %>% 
  arrange(fips, date) %>% 
  group_by(fips) %>% 
  mutate(gwr_smk10 = gwr_smk/10,
         gwr_smk_lag1 = lag(gwr_smk10, 1, order_by = fips),
         gwr_smk_lag2 = lag(gwr_smk10, 2, order_by = fips),
         gwr_smk_lag3 = lag(gwr_smk10, 3, order_by = fips),
         gwr_smk_lag4 = lag(gwr_smk10, 4, order_by = fips),         
         gwr_smk_lag5 = lag(gwr_smk10, 5, order_by = fips),
         gwr_smk_lag6 = lag(gwr_smk10, 6, order_by = fips)) %>% 
  select(-gwr_smk)

# read 2015 Washington CHARS data ----
# i decided against reading chars observation data as I think it mostly
# duplicates inpatient observations; inpatient is most comparable to colorado too

# read chars 2015 inpatient data
chars_2015_inpat <- haven::read_sas(paste0("./data/health/washington_chars/",
  "chr2015inpat.sas7bdat")) %>% 
  # select certain variables
  select(SEQ_NO_ENC, HOSPITAL, ZIPCODE, STATERES, COUNTYRES, AGE, SEX, 
         ADM_DATE, DIS_DATE, LENSTAYD, ADM_TYPE, DIAG1) %>% 
  # filter to ED admission
  filter(ADM_TYPE == "1") %>% 
  # filter to primary diagnosis of cardiopulmonary
  filter(DIAG1 %in% c(icd9_outcomes$resp, icd9_outcomes$cvd)) %>% 
  # fitler to admits in year
  mutate(year = as.character(lubridate::year(ADM_DATE))) %>% 
  filter(year == "2015") %>% 
  # add outcome indicators, age indicators
  mutate(date = ADM_DATE,
         resp = if_else(DIAG1 %in% icd9_outcomes$resp, 1, 0),
         asthma = if_else(DIAG1 %in% icd9_outcomes$asthma, 1, 0),
         copd = if_else(DIAG1 %in% icd9_outcomes$copd, 1, 0),
         acute_bronch = if_else(DIAG1 %in% icd9_outcomes$acute_bronch, 1, 0),
         pneum = if_else(DIAG1 %in% icd9_outcomes$pneumonia, 1, 0),
         cvd = if_else(DIAG1 %in% icd9_outcomes$cvd, 1, 0),
         arrhythmia = if_else(DIAG1 %in% icd9_outcomes$arrhythmia, 1, 0),
         cereb_vas = if_else(DIAG1 %in% icd9_outcomes$cereb, 1, 0),
         hf = if_else(DIAG1 %in% icd9_outcomes$hf, 1, 0), 
         ihd = if_else(DIAG1 %in% icd9_outcomes$ihd, 1, 0),
         mi = if_else(DIAG1 %in% icd9_outcomes$mi, 1, 0),
         age_cat = case_when(AGE < 15 ~ "age_under_15",
                             AGE >= 15 & AGE < 65 ~ "age_15_to_65",
                             AGE >= 65 ~ "age_over_65"),
         COUNTYRES = as.character(COUNTYRES),
         age = AGE, 
         sex = SEX) %>% 
  # join county key
  filter(STATERES == "WA") %>% 
  # filter to specific dates
  filter(date >= "2015-06-01" & date <= "2015-09-30") %>% 
  left_join(county_key, by = "COUNTYRES") 

# list of case-crossover dataframes joined to lagged pm values
# vector of outcomes
outcome_vector <- names(icd9_outcomes)


ts_casecross_list <- icd9_outcomes %>% 
  # create lists of outcome dataframes
  map(~filter(chars_2015_inpat, DIAG1 %in% .)) %>% 
  # filter to date ranges of pm data; function won't work if dates are different
  map(~filter(., date >= "2015-06-01" & date <= "2015-09-30")) %>% 
  # apply ts casecross over function and join with pm data
  map(~time_stratified(data=., id="SEQ_NO_ENC", 
        covariate = c("DIAG1", "ZIPCODE", "fips", "age", "age_cat", "sex"),
        admit_date = "date", start_date = "2015-06-01", end_date = "2015-09-30", 
        interval = 7) %>% 
      mutate(date = as.Date(date, format = "%Y-%m-%d"),
         outcome = as.numeric(as.character(outcome))) %>% 
      left_join(gwr_lag, by = c("fips", "date")))
  
head(ts_casecross_list[[2]])

# asthma casecross over test ----
copd <- chars_2015_inpat %>% 
  filter(copd == 1) %>% 
  filter(date >= "2015-06-01" & date <= "2015-09-30") %>% 
  filter(STATERES == "WA") %>% 
  filter(age_cat == "age_over_65")

copd_cc <- time_stratified(data=copd, id="SEQ_NO_ENC", 
  covariate = c("ZIPCODE", "fips", "AGE", "age_cat", "SEX"),
  admit_date = "date", start_date = "2015-06-01", end_date = "2015-09-30", 
  interval = 7) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         outcome = as.numeric(as.character(outcome))) %>% 
  left_join(county_pm, by = c("fips", "date")) %>% 
  mutate(gwr_smk10 = gwr_smk/10)




mod <- clogit(outcome ~ gwr_smk10 + strata(identifier), data = copd_cc)

summary(mod)

