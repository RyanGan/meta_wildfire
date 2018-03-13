# ------------------------------------------------------------------------------
# Title: Creation of ER inpatient outcome dataframe for CO, WA, OR 2010-2015
# Author: Ryan Gan
# Date Created: 2018-03-13
# ------------------------------------------------------------------------------

# this script aggregates cardiopulmonary diagnosis observations for Colorado, 
# Washington, and Oregon.

# library ----
library(tidyverse)
library(lubridate)

# load Rdata vector of outcomes ----
load("./data/health/icd9_outcome_vectors.RData")

# load hospital morbidity data ----
# define path where hosp data is
co_path <- paste0("../colorado_wildfire/data/health/co_hosp_1015.csv")
# import path
co_hosp <- read_csv(co_path, col_types = cols(.default = "c")) %>% 
  # filter to ER admission
  filter(ADMTNO == 1) %>% 
  # filter to study date
  mutate(date = as.Date(admit, format = "%m/%d/%Y")) %>% 
  # code outcomes
  mutate(resp = if_else(dx1 %in% icd9_outcomes$resp, 1, 0),
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
         sex  = case_when(SEX == "1" ~ "M",
                          SEX == "2" ~ "F"),
         age = as.numeric(AGEYRS),
         age_cat = case_when(AGEYRS < 15 ~ "age_under_15",
                             AGEYRS >= 15 & AGEYRS < 65 ~ "age_15_to_65",
                             AGEYRS >= 65 ~ "age_over_65"),
         zip = as.character(ZIP),
         id = paste0(cdpheid,"c"),
         state = "colorado",
         fips = paste0("08", county_final),
         stay_length = as.numeric(LOS)) %>% 
  select(id, date, stay_length, dx1, age, age_cat, sex, zip, fips, state)

# check variables
glimpse(co_hosp)

# washington -------------------------------------------------------------------
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

# create list of years to read in
years_list <- list("2011", "2012", "2013", "2014", "2015")

# function to read in health data
wa_hosp <-  years_list %>% 
  map(., function(x){
    # read inpatient and observation
    inpat_path <- paste0("./data/health/washington_chars/chr",x,"inpat.sas7bdat")
    
    # read inpatient SAS files using haven package
    chars_inpat <- haven::read_sas(inpat_path)  %>% 
      mutate(COUNTYRES = as.character(COUNTYRES)) %>% 
      # filter to ER admission
      filter(ADM_TYPE == "1") %>% 
      # filter to primary diagnosis of cardiopulmonary
      filter(DIAG1 %in% c(icd9_outcomes$resp, icd9_outcomes$cvd)) %>% 
      # fitler to admits in year
      mutate(year = as.character(lubridate::year(ADM_DATE))) %>% 
      filter(year == x) %>% 
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
             # shared covariates
             age = AGE,
             age_cat = case_when(AGE < 15 ~ "age_under_15",
                                 AGE >= 15 & AGE < 65 ~ "age_15_to_65",
                                 AGE >= 65 ~ "age_over_65"),
             zip = as.character(ZIPCODE),
             id = paste0(SEQ_NO_ENC,"w"),
             state = "washington",
             sex = SEX,
             stay_length = as.numeric(LENSTAYD),
             dx1 = DIAG1) %>% 
      # join county key
      filter(STATERES == "WA") %>% 
      left_join(county_key, by = "COUNTYRES") %>% 
      select(id, date, stay_length, dx1, age, age_cat, sex, zip, fips, state)
  } # end function
  ) %>%  # end map
  # bind rows
  map_dfr(bind_rows)


# oregon ------
or_path <- paste0("./data/health/2013-oregon_er_inpatient.csv")

# i need jingyang's zip to county key
or_hosp <- read_csv(or_path, col_types = cols(.default = "c")) %>% 
  rename(id = personkey, sex = gender, 
         zip = ZIP, state = STATE) %>% 
  mutate(id = as.character(paste0(id,"o")),
         date = as.Date(fromdate),
         stay_length = as.numeric(los),
         fips = as.character(NA),
         state = "oregon", 
         age = as.numeric(age),
         age_cat = case_when(age < 15 ~ "age_under_15",
                             age >= 15 & age < 65 ~ "age_15_to_65",
                             age >= 65 ~ "age_over_65")) %>% 
  # filter to place of service inpatient
  filter(pos == 21) %>% 
  select(id, date, stay_length, dx1, age, age_cat, sex, zip, fips, state)

glimpse(or_hosp)
# combine ER inpatient datasets -----
er_hosp <- bind_rows(co_hosp, wa_hosp, or_hosp)

# saving er_hosp 
write_csv(er_hosp, "./data/health/1015-er_hosp.csv")