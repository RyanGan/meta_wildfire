# ------------------------------------------------------------------------------
# Title: Creation of mortaility underlying cause of death by ICD10 list
# Author: Ryan Gan
# Date Created: 2018-03-20
# ------------------------------------------------------------------------------

# library
library(tidyverse)
library(stringr)

# creating a list of vectors of casue of death
outcome <- c("resp", "asthma", "copd", "cvd", "hf", "card_arrest",
             "ihd", "mi", "cereb_vas")

# vector of icd10 codes of first letter and first 2 digits
resp <- paste0("J", c(paste0("0", c(0:9)),c(10:98)))
asthma <- paste0("J", c(45:46))
copd <- "J44"
cvd <- paste0("I", c(paste0("0", 0:9), c(10:78)))
hf <- "I50"
card_arrest <- "I46"
ihd <- paste0("I", c(20:25))
mi <- paste0("I", c(21:22))
cereb_vas <- paste0("I", c(60:69))

# create list of first letter and first 2 digits of ICD10 code
icd10_outcomes <- list(resp, asthma, copd, cvd, hf, card_arrest, ihd, mi, 
                       cereb_vas)
# assign names
names(icd10_outcomes) <- outcome
# save list as R file
save(icd10_outcomes, file = "./data/health/icd10_outcome.RData")

