# ------------------------------------------------------------------------------
# Title: Time-stratified case-crossover dataframes for inpatient emergency 
#        department hospitalizations 
# Author: Ryan Gan
# Date Created: 2018-02-06
# ------------------------------------------------------------------------------

# I think I'll run this on the server; maybe in parallel

# libraries ----
library(tidyverse)
library(haven)

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
         COUNTYRES = as.character(COUNTYRES)) %>% 
  # join county key
  filter(STATERES == "WA") %>% 
  left_join(county_key, by = "COUNTYRES")

# casecross over for each outcome -----
# custom function
time.stratified <- function(data, id, covariate=F, admit_date,
                            start_date, end_date, interval){
  # if id value is given
  if(is.character(id)){
    # vector of ids
    id_vec <- data[,id]
  } else { # else if no value provided, create an id variable
    id_vec <- seq(1, nrow(data), by=1)
  }
  # vector of admit dates joined to the id vector
  admit_date <- data[,admit_date]
  id_date <- data.frame(id_vec, admit_date)

  # create list of vectors of referent dates and admission date
  referent_date_list <- apply(id_date, 1, function(x){
    date_prior <- as.character(seq(as.Date(x[2]),as.Date(start_date),
                                   by = (-interval)))
    date_posterior <- as.character(seq(as.Date(x[2]),as.Date(end_date),
                                       by = (interval)))
    date <- sort(unique(c(date_prior, date_posterior)))
    identifier <- rep(x[1], length(date))
    # outcome
    outcome <- ifelse(date %in% x[2], 1, 0)
    id_date_vec <- cbind(identifier, date, outcome)
    return(id_date_vec)
  }) # end apply
  # timestrat dataframe
  ts_data <- do.call(rbind, referent_date_list)
  rownames(ts_data) <- NULL
  ts_data <- as.data.frame(ts_data)
  # remove white space
  ts_data$identifier <- gsub("[[:space:]]", "", 
                             as.character(ts_data$identifier))
  # if covariates provided, join to the dataframe
  if(is.character(covariate)){
    cov_data <- as.data.frame(cbind(id_vec, data[,covariate]))
    # names of cov data
    colnames(cov_data) <- c("identifier", covariate)
    # conver identifier to character to merge
    cov_data$identifier <- as.character(cov_data$identifier)
    # merge with ts_data
    ts_data <- left_join(ts_data, cov_data, by = "identifier")
  }
  # return dataframe
  return(ts_data)
}

# I need to figure out how to update my case.crossover package

# asthma subset check
asthma <- chars_2015_inpat %>% 
  filter(asthma == 1) %>% 
  filter(date >= "2015-06-01" & date <= "2015-09-30") %>% 
  filter(STATERES == "WA")


asthma_cc <- time.stratified(data=asthma, id="SEQ_NO_ENC", 
                             covariate = c("ZIPCODE", "fips", "AGE", "SEX"),
                             admit_date = "date", start_date = "2015-06-01",
                             end_date = "2015-09-30", interval = 7)


glimpse(asthma)


