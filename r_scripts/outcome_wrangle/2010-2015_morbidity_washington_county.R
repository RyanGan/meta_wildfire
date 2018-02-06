# ------------------------------------------------------------------------------
# Title: Data wrangling and creation of county time series counts of 
#   cardiopulmonary emergency department visits in Washington state 2010-2015
# Author: Ryan Gan
# Date Created: 2018-01-25
# ------------------------------------------------------------------------------

# libraries ----
# loading tidyverse
library(tidyverse)

# load icd9 diagnosis codes ----
# read icd9 files
icd9_path <- paste0("./data/health/washington_chars/MedicareDiagCodeV32.xlsx")
icd9 <- readxl::read_xlsx(icd9_path) %>% 
  rename(icd9 = 'DIAGNOSIS CODE', long_desc = 'LONG DESCRIPTION',
         short_desc = 'SHORT DESCRIPTION')

# vector of cardio pulmonary outcomes
resp_vector <- icd9 %>% 
  slice(which(icd9 == '460'):which(icd9 == '5199')) %>% 
  select(icd9) %>% 
  as_vector()

# asthma icd9
asthma_vector <- icd9 %>% 
  slice(which(icd9 == '49300'):which(icd9 == '49392')) %>% 
  select(icd9) %>% 
  as_vector()

# copd icd9
copd_vector <- c('490', '4910','4911','49120','49121','49122','4918','4919', 
                 '4920','4928', '4940', '4941', '496') 

# pneumonia
pneumonia_vector <- icd9 %>% 
  slice(which(icd9 == '4800'):which(icd9 == '486')) %>% 
  select(icd9) %>% 
  as_vector()

# acute bronch
acute_bronch_vector <- icd9 %>% 
  slice(which(icd9 == '4660'):which(icd9 == '46619')) %>% 
  select(icd9) %>% 
  as_vector()

# cvd vector
cvd_vector <- icd9 %>% 
  slice(which(icd9 == '390'):which(icd9 == '4599')) %>% 
  select(icd9) %>% 
  as_vector()

# ihd
ihd_vector <- icd9 %>% 
  slice(which(icd9 == '41000'):which(icd9 == '4139')) %>% 
  select(icd9) %>% 
  as_vector()
  
# arrhythmia
arrhythmia_vector <- icd9 %>% 
  slice(which(icd9 == '4270'):which(icd9 == '4279')) %>% 
  select(icd9) %>% 
  as_vector()

# heart failur
hf_vector <- icd9 %>% 
  slice(which(icd9 == '4280'):which(icd9 == '4289')) %>% 
  select(icd9) %>% 
  as_vector()

# cerebrovas
cereb_vector <- icd9 %>% 
  slice(which(icd9 == '430'):which(icd9 == '4389')) %>% 
  select(icd9) %>% 
  as_vector()

# myocardial infarction
mi_vector <- c("41000", "41001", "41002", "41010", "41011", "41012", "41020", 
               "41021", "41022", "41030", "41031", "41032", "41040", "41041",
               "41042", "41050", "41051", "41052", "41060", "41061", "41062",
               "41070", "41071", "41072", "41080", "41081", "41082", "41090",
               "41091", "41092")

# bind outcomes in a list and save as an R object for other states

icd9_outcomes <- list(resp_vector, asthma_vector, copd_vector, 
  acute_bronch_vector, pneumonia_vector, cvd_vector, arrhythmia_vector,
  ihd_vector, mi_vector, hf_vector, cereb_vector)

# assign names
names(icd9_outcomes) <- c("resp", "asthma", "copd", "acute_bronch", "pneumonia",
                          "cvd", "arrhythmia", "ihd", "mi", "hf", "cereb_vas")

save(icd9_outcomes, file = "./data/health/icd9_outcome_vectors.RData")

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

# load CHARS data ----
# read revisit with PATIENTID
id_path <- paste0("./data/health/washington_chars/",
                  "confrevisit2011_2015.sas7bdat")
id <- haven::read_sas(id_path)

# create function purrr; need to be aware that it seems like after 2011
# obs file labels AGE_O

# define function to read in inpatient and observation datasets for each year
# and estimate time sereis of emergency room inpatient/hospitalizations
# for washington
timeseries <- function(x){
        # read inpatient and observation
        inpat_path <- paste0("./data/health/washington_chars/chr",x,"inpat.sas7bdat")
        obs_path <- paste0("./data/health/washington_chars/chr",x,"obs.sas7bdat")
        
        # read SAS files using haven package
        # import inpatient
        chars_inpat <- haven::read_sas(inpat_path)  %>% 
          # select specific variables
          select(SEQ_NO_ENC, ADM_TYPE, ADM_DATE, AGE, SEX, 
                 COUNTYRES, STATERES, COUNTYRES, DIAG1)
# Note 2018-02-06: I commented out observation datset; identical info -----
        # # import observation
        # chars_obs <- haven::read_sas(obs_path) %>% 
        #   # select specific variables
        #   select(SEQ_NO_ENC, ADM_TYPE, ADM_DATE, AGE, SEX, 
        #          COUNTYRES, STATERES, COUNTYRES, DIAG1)
        # 
        # row bind chars inpatient and obs together
        #chars <- bind_rows(chars_inpat, chars_obs) %>% 
        
        
        chars <- chars_inpat %>% 
          # join by seq_no_enc
          left_join(id, by = "SEQ_NO_ENC") %>% 
          # filter to ER visits
          filter(ADM_TYPE == "1") %>% 
          # filter to primary diagnosis of cardiopulmonary
          filter(DIAG1 %in% c(resp_vector, cvd_vector)) %>% 
          # fitler to admits in year
          mutate(year = as.character(lubridate::year(ADM_DATE))) %>% 
          filter(year == x) %>% 
          # add outcome indicators, age indicators
          mutate(date = ADM_DATE,
                 resp = if_else(DIAG1 %in% resp_vector, 1, 0),
                 asthma = if_else(DIAG1 %in% asthma_vector, 1, 0),
                 copd = if_else(DIAG1 %in% copd_vector, 1, 0),
                 acute_bronch = if_else(DIAG1 %in% acute_bronch_vector, 1, 0),
                 pneum = if_else(DIAG1 %in% pneumonia_vector, 1, 0),
                 cvd = if_else(DIAG1 %in% cvd_vector, 1, 0),
                 arrhythmia = if_else(DIAG1 %in% arrhythmia_vector, 1, 0),
                 cereb_vas = if_else(DIAG1 %in% cereb_vector, 1, 0),
                 hf = if_else(DIAG1 %in% hf_vector, 1, 0), 
                 ihd = if_else(DIAG1 %in% ihd_vector, 1, 0),
                 mi = if_else(DIAG1 %in% mi_vector, 1, 0),
                 age_cat = case_when(AGE < 15 ~ "age_under_15",
                                     AGE >= 15 & AGE < 65 ~ "age_15_to_65",
                                     AGE >= 65 ~ "age_over_65")) %>% 
          # join county key
          filter(STATERES == "WA") %>% 
          left_join(county_key, by = "COUNTYRES")
        
        # total counts
        wa_total_ts <- chars %>% 
          group_by(date, fips) %>% 
          summarise_at(vars(resp:mi), funs(sum(.))) %>% 
          mutate(cardiopulm_n = resp+cvd,
                 strata = "total") %>% 
          select(date, fips, strata, resp:cardiopulm_n) %>% 
          arrange(date, fips)
        
        # age strata
        wa_age_ts <- chars %>% 
          group_by(date, fips, age_cat) %>% 
          summarise_at(vars(resp:mi), funs(sum(.))) %>% 
          mutate(cardiopulm_n = resp+cvd) %>% 
          rename(strata = age_cat) %>% 
          arrange(date, fips)
        
        # sex strata
        wa_sex_ts <- chars %>% 
          group_by(date, fips, SEX) %>% 
          summarise_at(vars(resp:mi), funs(sum(.))) %>% 
          mutate(cardiopulm_n = resp+cvd) %>% 
          rename(strata = SEX) %>% 
          arrange(date, fips)
        
        # rowbind all timeseries together
        wa_ts <- bind_rows(wa_total_ts, wa_age_ts, wa_sex_ts)
        return(wa_ts)
      } # end function

# create list of years to read in
years_list <- list("2011", "2012", "2013", "2014", "2015")

xtabs(~strata, washington_timeseries)

# may be some missing dates for some counties for strata
strata_seq <- c("total", "age_15_to_65", "age_over_65", "age_under_15",
                  "M", "F")

date_seq <- as.character(seq.Date(as.Date("2011-01-01"), as.Date("2015-09-30"), 
  by = "day"))

date_fips <- expand.grid(date_seq, fips, strata_seq) %>% 
  rename(date = Var1, fips = Var2, strata = Var3) %>% 
  mutate(date = as.Date(date), fips = as.character(fips), 
         strata = as.character(strata))

# create time series dataframe -----
# fead in years list to purrr map function
washington_timeseries <- years_list %>% 
  map(., timeseries) %>% # apply timeseries function9
  map_dfr(bind_rows) %>%  # bind rows
  # get rid of missing strata
  filter(strata != "U" & !is.na(fips)) %>%
  # join to counties and dates missing counts
  full_join(date_fips, by = c("date", "fips", "strata")) %>% 
  # set missing values to 0
  mutate_at(vars(resp:cardiopulm_n), funs(if_else(is.na(.),0, .)))

summary(washington_timeseries)
summary(as.factor(washington_timeseries$fips))
summary(as.factor(washington_timeseries$strata))

# write timeseries csv file ----
# write path
write_path <- paste0("./data/health/2011-2015_morbidity_wa_ts.csv")
write_csv(washington_timeseries, write_path)
