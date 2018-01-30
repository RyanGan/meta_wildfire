# ------------------------------------------------------------------------------
# Title: Data wrangling and creation of county time series counts of 
#   cardiopulmonary emergency department visits in Colorado state 2010-2015
# Author: Ryan Gan
# Date Created: 2018-01-25
# ------------------------------------------------------------------------------

# libraries ----
# loading tidyverse
library(tidyverse)

# load Rdata vector of outcomes ----
load("./data/health/icd9_outcome_vectors.RData")

icd9_outcomes
# load hospital morbidity data ----
co_path <- paste0("../colorado_wildfire/data/health/co_hosp_1015.csv")
co_hosp <- read_csv(co_path, col_types = cols(.default = "c")) 

# check variables
glimpse(co_hosp)
# check patient type; it looks like Colorado ED visits are inpatient only
# this is similar to CHARS
summary(as.factor(co_hosp$PATTYPE))

summary(as.factor(co_hosp$PLACSERV))

xtabs(~PATTYPE + PLACSERV, data = co_hosp)

xtabs(~ADMSNO, data = co_hosp)
names(icd9_outcomes)

# I think I'll summarise up inpatients admitted from ER; should be similar to 
# what I did for CHARS
co_ed <- co_hosp %>%
  filter(ADMTNO == 1) %>% # subset to ADMTNO == 1 Emergency
  mutate(date = as.Date(admit, format = "%m/%d/%Y"),
    resp = if_else(dx1 %in% icd9_outcomes$resp, 1, 0),
    asthma = if_else(dx1 %in% icd9_outcomes$asthma, 1, 0),
    copd = if_else(dx1 %in% icd9_outcomes$copd, 1, 0),
    acute_bronch = if_else(dx1 %in% icd9_outcomes$acute_bronch, 1, 0),
    pneum = if_else(dx1 %in% icd9_outcomes$pneumonia, 1, 0),
    cvd = if_else(dx1 %in% icd9_outcomes$cvd, 1, 0),
    arrhythmia = if_else(dx1 %in% icd9_outcomes$arrhythmia, 1, 0),
    cereb_vas = if_else(dx1 %in% icd9_outcomes$cereb_vas, 1, 0),
    hf = if_else(dx1 %in% icd9_outcomes$hf, 1, 0), 
    ihd = if_else(dx1 %in% icd9_outcomes$ihd, 1, 0),
    mi = if_else(dx1 %in% icd9_outcomes$mi, 1, 0),
    sex_factor = case_when(SEX == "1" ~ "M",
                           SEX == "2" ~ "F"),
    AGEYRS = as.numeric(AGEYRS),
    age_cat = case_when(AGEYRS < 15 ~ "age_under_15",
                        AGEYRS >= 15 & AGEYRS < 65 ~ "age_15_to_65",
                        AGEYRS >= 65 ~ "age_over_65"),
    fips = paste0("08", county_final))

# time series ----
      # total counts
      co_total_ts <- co_ed %>% 
        group_by(date, fips) %>% 
        summarise_at(vars(resp:mi), funs(sum(.))) %>% 
        mutate(cardiopulm_n = resp+cvd,
               strata = "total") %>% 
        select(date, strata, fips, resp:cardiopulm_n) %>% 
        arrange(date, fips)

      # age strata
      co_age_ts <- co_ed  %>% 
        group_by(date, fips, age_cat) %>% 
        summarise_at(vars(resp:mi), funs(sum(.))) %>% 
        mutate(cardiopulm_n = resp+cvd) %>% 
        rename(strata = age_cat) %>% 
        arrange(date, fips)
      
      # sex strata
      co_sex_ts <- co_ed %>% 
        group_by(date, fips, sex_factor) %>% 
        summarise_at(vars(resp:mi), funs(sum(.))) %>% 
        mutate(cardiopulm_n = resp+cvd) %>% 
        rename(strata = sex_factor) %>% 
        arrange(date, fips)
      
      # rowbind all timeseries together
      co_ts <- bind_rows(co_total_ts, co_age_ts, co_sex_ts)

# may be some missing dates for some counties for strata
strata_seq <- c("total", "age_15_to_65", "age_over_65", "age_under_15",
                "M", "F")
      
date_seq <- as.character(seq.Date(as.Date("2010-01-01"), as.Date("2015-09-27"), 
                                        by = "day"))

fips <- unique(co_ts$fips)[!is.na(unique(co_ts$fips))]

# expand grid
date_fips <- expand.grid(date_seq, fips, strata_seq) %>% 
        rename(date = Var1, fips = Var2, strata = Var3) %>% 
        mutate(date = as.Date(date), fips = as.character(fips), 
               strata = as.character(strata))

# create time series dataframe -----
  # fead in years list to purrr map function
  colorado_timeseries <- co_ts %>% 
    # join to counties and dates missing counts
    full_join(date_fips, by = c("date", "fips", "strata")) %>% 
    # set missing values to 0
    mutate_at(vars(resp:cardiopulm_n), funs(if_else(is.na(.),0, .))) %>% 
    filter(!is.na(strata))

summary(colorado_timeseries)  
summary(as.factor(colorado_timeseries$strata))
summary(as.factor(colorado_timeseries$fips))
# write timeseries csv file ----
# write path
write_path <- paste0("./data/health/2010-2015_morbidity_co_ts.csv")
write_csv(colorado_timeseries, write_path)
