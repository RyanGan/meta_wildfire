# ------------------------------------------------------------------------------
# Title: Creation of morbidity time-series dataframe for Western US 2010-2015
# Author: Ryan Gan
# Date Created: 2018-02-27
# ------------------------------------------------------------------------------

# library ----
library(tidyverse)

# import envrionmental/exposure data ----
# pm file path
pm_file_path <- paste0("./data/smoke/krig_county_popwt/")

# importing pm2.5 data 2010 to 2015
# krig pm  dataframe
krig_pm <- (list.files(pm_file_path, pattern="krig_pm")[5:10]) %>%  
  # read in each csv file
  map(~read_csv(paste0(pm_file_path, .))) %>% 
  # transform data wide to long 
  map_dfr(~ gather(., date, pm_krig, -fips) %>% 
            mutate(date = as.Date(str_sub(date, start=2), format="%Y%m%d"))) 

# background pm
bg_pm <- (list.files(pm_file_path, pattern="background_pm")[5:10]) %>% 
  # read in each csv file
  map(~read_csv(paste0(pm_file_path, .))) %>% 
  # transform data wide to long and retain only colorado and washington fips
  # bind them in a list
  map_dfr(~ gather(., date, bg_pm, -fips) %>% 
            mutate(date = as.Date(str_sub(date, start=2), format="%Y%m%d"))) 

# hms county proportion
hms <- (list.files(pm_file_path, pattern="hms")[5:10]) %>% 
  # read in each csv file
  map(~read_csv(paste0(pm_file_path, .))) %>% 
  # transform data wide to long and retain only colorado and washington fips
  # bind them in a list
  map_dfr(~ gather(., date, hms, -fips) %>% 
            mutate(date = as.Date(str_sub(date, start=2), format="%Y%m%d")))

# read in temp data
temp <- list.files(paste0("./data/smoke/county_air_temp/")) %>% 
  map(~read_csv(paste0("./data/smoke/county_air_temp/", .))) %>% 
  # transform data wide to long and retain only colorado and washington fips
  # bind them in a list
  map_dfr(~ gather(., date, temp_k, -fips) %>% 
            mutate(date = as.Date(str_sub(date, start=2), format="%Y%m%d"),
                   temp_f = round(temp_k * (9/5) - 459.67, 2)))

# state fp code to join
state_fp <- sort(unique(str_sub(krig_pm$fips, start = 1, end =2)))
# state character name
state <- c("Arizona", "California", "Colorado", "Idaho", "Montana",
           "Nevada", "New Mexico", "Oregon", "Utah", "Washington", "Wyoming")
# create key
state_key <- data.frame(cbind(state_fp, state)) %>% 
  mutate_if(is.factor, as.character)

# join all exposure values by fips and date
county_pm_ts <- list(krig_pm, bg_pm, hms, temp) %>% 
  # bind list together
  reduce(function(df1,df2)left_join(df1,df2, by = c("fips", "date"))) %>% 
  mutate(month = lubridate::month(date),
         pm_krig = ifelse(pm_krig < 0, 0, pm_krig), # set values lower than 0 to 0
         # create smoke pm variable
         pm_smk = ifelse(pm_krig - bg_pm > 0, pm_krig - bg_pm, 0),
         # create binary smoke variables
         smoke0 = ifelse(pm_smk > 0 & month %in% c(5:9), 1, 0),
         smoke5 = ifelse(pm_smk > 5 & month %in% c(5:9), 1, 0), 
         smoke10 = ifelse(pm_smk > 10 & month %in% c(5:9), 1, 0),
         smoke15 = ifelse(pm_smk > 15 & month %in% c(5:9), 1, 0),
         # create state variables based on fips codes
         state_fp = as.character(str_sub(fips, start = 1, end = 2))) %>% 
  # join state names in
  left_join(state_key, by = "state_fp") %>% 
  # create smoke wave (98th perctile for 2 days)
  group_by(fips) %>% 
  arrange(fips, date) %>% 
  mutate(high_pm_day = ifelse(pm_krig >= 22.31, 1, 0),
         smoke_wave = ifelse(lag(high_pm_day, order_by = fips)==1 &
                               high_pm_day==1 & month %in% c(5:9), 1, 0)) %>% 
  select(state, state_fp, fips, date, month, pm_krig, bg_pm, pm_smk, hms, temp_k,
         temp_f, smoke0, smoke5, smoke10, smoke15, high_pm_day, smoke_wave)

# write csv file for rish
write_csv(county_pm_ts, "./data/smoke/1015-county_popwt_pm.csv")

# removing files to save envrionment memory
# rm(bg_pm, hms, krig_pm, temp)

# county population denom ----
# colorado census
co_denom <- read_csv("./data/census/2016-colorado_population.csv") %>% 
  select(GEO.id2, 'GEO.display-label', respop72010:respop72015) %>% 
  rename(fips = GEO.id2, county = 'GEO.display-label') %>% 
  rename_at(vars(respop72010:respop72015), 
            funs(paste0("pop",str_sub(., start=8)))) %>% 
  # take out state name and "county" from county name
  mutate(county = str_split_fixed(county, " County", n=2)[,1])
# sum up counties to get a state total
co_denom_total <- co_denom %>% 
  gather(key="year", value = "population", pop2010:pop2015) %>% 
  mutate(year = as.numeric(str_sub(year, start = 4)),
         fips = as.character(fips)) 

# washington cesus
wa_denom <- read_csv("./data/census/2016-washington_population.csv") %>% 
  select(GEO.id2, 'GEO.display-label', respop72010:respop72015) %>% 
  rename(fips = GEO.id2, county = 'GEO.display-label') %>% 
  rename_at(vars(respop72010:respop72015), 
            funs(paste0("pop",str_sub(., start=8)))) %>% 
  # take out state name and "county" from county name
  mutate(county = str_split_fixed(county, " County", n=2)[,1])
# sum up counties to get a state total
wa_denom_total <- wa_denom %>% 
  gather(key="year", value = "population", pop2010:pop2015) %>% 
  mutate(year = as.numeric(str_sub(year, start = 4)),
         fips = as.character(fips))

# colorado/washington county denom
pop_denom <- bind_rows(co_denom_total, wa_denom_total)

# import morbidity count data ----
# washington
# define relative path
wa_path <- paste0("./data/health/2011-2015_morbidity_wa_ts.csv")
# read in dataframe
wa_ts <- read_csv(wa_path, progress = F) %>% 
  mutate(fips = as.character(fips)) %>% 
  # limiting to total strata 
  filter(strata == "total")

# colorado
# define relative path to data folder
co_path <- paste0("./data/health/2010-2015_morbidity_co_ts.csv")
# read in dataframe using read_csv function from tidyverse package
co_ts <- read_csv(co_path, progress = F) %>% 
  mutate(fips = as.character(fips)) %>% 
  # limiting to total strata for now
  filter(strata == "total")

# oregon data


# final morbidity time series ----
# bind with county counts 
morbidity_ts <- bind_rows(co_ts, wa_ts) %>% 
  mutate(year = lubridate::year(date)) %>% 
  # join population
  left_join(pop_denom, by = c("fips", "year")) %>% 
  # join pm
  left_join(county_pm_ts, by = c("fips", "date")) %>% 
  select(state, state_fp, fips, county, date, year, population,
         resp:cardiopulm_n, pm_krig:smoke_wave)

# note Washington missing values for washington 2010 dates that I don't have
# health outcomes data for
glimpse(morbidity_ts)
# write file
write_csv(morbidity_ts, "./data/health/1015-morbidity_pm_ts.csv")

