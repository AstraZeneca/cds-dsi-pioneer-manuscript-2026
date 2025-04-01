if(!require(here))install.packages("here")
library(here)

i_am("historical_data_integration/main_pipeline_data.R")


source(here("historical_data_integration", "get_data.R"))
source(here("historical_data_integration","wrangle_data.R"))
source(here("historical_data_integration","other_data.R"))

if(!require(schoolmath))install.packages("schoolmath")
library(schoolmath) # Custom function needs it

if(!require(dplyr))install.packages("dplyr")
library(dplyr) # wrangle_data needs it

if(!require(solvebio))install.packages("solvebio")
library(solvebio) # Only if you need datasets not already in Domino

if(!require(readxl))install.packages("readxl")
library(readxl) # needed for dictionary

if(!require(tidyr))install.packages("tidyr")
library(tidyr)

if(!require(stringr))install.packages("stringr")
library(stringr)

if(!require(data.table))install.packages("data.table")
library(data.table) # Only for _speed functions. NOT used

if(!require(lubridate))install.packages("lubridate")
library(lubridate) # Only for _speed functions. NOT used

login() # Needs: Access to solvebio, PAT created in QuartzBio, PAT from QB added
# into domino as Env variable. Not needed if the data is already in Domino.

help_get_dataset(information="general")
help_get_dataset(information="dataset_type")


## ADAM data

# Clinical data
c_db <- download_dataset_API(dataset="clinical", download=FALSE)

# Clinical longitudinal data
cl_db_adam <- download_dataset_API(dataset="clinical_longitudinal", 
                                   download=FALSE, data_source="adam")

# Response Longitudinal data
rl_db_adam <- download_dataset_API(dataset="response_longitudinal",
                                   download=FALSE, data_source="adam")

# Measurements longitudinal data
ml_db_adam <- download_dataset_API(dataset="measurements_longitudinal",
                                   download=FALSE, data_source="adam")

## Treat data

clinical_db <- transform_numbers_debug(c_db)
clinical_db <- transform_dates(clinical_db)

clinical_longitudinal <- transform_numbers_debug(cl_db_adam)
clinical_longitudinal <- transform_dates(clinical_longitudinal)
clinical_longitudinal <- transform_names(clinical_longitudinal,
                                         data_source="adam") # human names


response_longitudinal <- transform_numbers_debug(rl_db_adam)
response_longitudinal <- transform_dates(response_longitudinal)
response_longitudinal <- transform_names(response_longitudinal,
                                         data_source="adam") # human names


measurements_longitudinal <- transform_numbers_debug(ml_db_adam)
measurements_longitudinal <- transform_dates(measurements_longitudinal)
measurements_longitudinal <- transform_names(measurements_longitudinal,
                                             data_source="adam") # human names
# QC

# TREATMENT ARM
ntable(clinical_db$actual_arm)

# Patient nbr and discordance between datasets
len(unique(clinical_db$USUBJID)) # 1840
len(unique(clinical_longitudinal$unique_subject_identifier)) # 1537 all are included in clinical
len(unique(response_longitudinal$unique_subject_identifier)) # 1099 all are included in clinical all in clinical_longitudinal and in measurements
len(unique(measurements_longitudinal$unique_subject_identifier)) # 1115 all are included in clinical; E6210008 not in clinical_longitudinal


# Merge data into schema #1

# Identify SF
sfs_ids_list <- clinical_db$USUBJID[which(clinical_db$actual_arm %in%
                                            c("SCREEN FAILURE", "NOT ASSIGNED"))]
clinical_db_no_SF <- subset(clinical_db, !clinical_db$USUBJID %in% sfs_ids_list)

ntable(clinical_db_no_SF$actual_arm) # Still some "not treated and SOC"

sfs_ids_list <- clinical_db$USUBJID[which(clinical_db$actual_arm %in%
                                            c("SCREEN FAILURE", "NOT ASSIGNED",
                                              "NOT TREATED", "SOC TREATMENT ARM"
                                              ))]


clinical_longitudinal_no_SF <- subset(
  clinical_longitudinal,
  !clinical_longitudinal$unique_subject_identifier %in% sfs_ids_list)

len(unique(clinical_longitudinal_no_SF$unique_subject_identifier))


measurements_longitudinal_no_SF <- subset(
  measurements_longitudinal,
  !measurements_longitudinal$unique_subject_identifier %in% sfs_ids_list)

len(unique(measurements_longitudinal_no_SF$unique_subject_identifier))


response_longitudinal_no_SF <- subset(
  response_longitudinal,
  !response_longitudinal$unique_subject_identifier %in% sfs_ids_list)

len(unique(response_longitudinal_no_SF$unique_subject_identifier))


#=================================
#     ASSESSMENT VISIT DATA      #
#=================================

# loop_assessment_visit_date_no_SF <- generate_assessment_visit_date_dataset(measurements_longitudinal_dataset = measurements_longitudinal_no_SF , response_longitudinal_dataset = response_longitudinal_no_SF)
response_longitudinal_no_SF$analysis_value_c <- 
  response_longitudinal_no_SF$`analysis_value_(c)`

claude_assessment_visit_date_no_SF <- 
  claude_opt_generate_assessment_visit_date_dplyr(
    measurements_longitudinal_dataset = measurements_longitudinal_no_SF,
    response_longitudinal_dataset = response_longitudinal_no_SF
    )

claude_assessment_visit_date_no_SF$response <- ifelse(
  claude_assessment_visit_date_no_SF$response=="",
  NA,
  claude_assessment_visit_date_no_SF$response
  )


# QC (2)
colSums(is.na(claude_assessment_visit_date_no_SF))
weirdos <- claude_assessment_visit_date_no_SF %>%
  filter(is.na(response)) %>%
  filter(visitnum != "1.00" & visitnum != "1")


assessment_visit_data <- claude_assessment_visit_date_no_SF 

assessment_visit_data <- calc_visit_date(assessment_visit_data) # READY!

#=================================
#           PATIENT DATA         #
#=================================

patient_data <- generate_patient_data_dataset(
  clinical_dataset=clinical_db_no_SF,
  assessment_visit_data = assessment_visit_data,
  measurements_longitudinal_data = measurements_longitudinal_no_SF,
  clinical_longitudinal_data = clinical_longitudinal_no_SF
) # Ready!

pacient_data <- remove_patient_data_dups(patient_data) # One patient is repeated 2 times!
patient_data <- remove_patient_data_dups(pacient_data) # READY!
# ============================
#            SPARE           #
# ============================

# ## STDM data
# 
# # c_db <- get_dataset(dataset="clinical")
# # c_db <- download_dataset_API(dataset="clinical", download=TRUE)
# 
# 
# #cl_db <- get_dataset(dataset="clinical_longitudinal")
# #cl_db <- download_dataset_API(dataset="clinical_longitudinal", download=TRUE)
# cl_db <- download_dataset_API(dataset="clinical_longitudinal", download=FALSE)
# 
# 
# # rl_db <- get_dataset(dataset="response_longitudinal")
# # rl_db <- download_dataset_API(dataset="response_longitudinal", download=TRUE)
# rl_db <- download_dataset_API(dataset="response_longitudinal", download=FALSE)
# 
# # ml_db <- get_dataset(dataset="measurements_longitudinal")
# # ml_db <- download_dataset_API(dataset="measurements_longitudinal", download=TRUE)
# ml_db <- download_dataset_API(dataset="measurements_longitudinal", download=FALSE)
# 
# 
# 
# 
# 
# response_longitudinal <- transform_numbers(rl_db)
# response_longitudinal <- transform_dates(response_longitudinal)
# response_longitudinal <- transform_names(response_longitudinal) # human names
# 
# # loop_assessment_visit_date <- generate_assessment_visit_date_dataset()
# # dplyr_assessment_visit_date <- generate_assessment_visit_date_dataset_dplyr()
# 
# 
# 
# 
# which(!paste0("D419AC00001/",unique(loop_assessment_visit_date$usubjid)) %in% unique(clinical_db$USUBJID))
# sum(unique(clinical_longitudinal$unique_subject_identifier) %in% unique(clinical_db$USUBJID))
# sum(unique(response_longitudinal$unique_subject_identifier) %in% unique(clinical_db$USUBJID))
# sum(unique(measurements_longitudinal$unique_subject_identifier) %in% unique(clinical_db$USUBJID))
# 
# 
# which(!unique(measurements_longitudinal$unique_subject_identifier) %in% unique(clinical_longitudinal$unique_subject_identifier))
# unique(measurements_longitudinal$unique_subject_identifier)[781]
# 
# sum(unique(response_longitudinal$unique_subject_identifier) %in% unique(clinical_longitudinal$unique_subject_identifier))
# sum(unique(response_longitudinal$unique_subject_identifier) %in% unique(measurements_longitudinal$unique_subject_identifier))
# 
# clinical_db$actual_arm[clinical_db$USUBJID=="E0302003"]
# 
# which(is.na(measurements_longitudinal$unique_subject_identifier))
# 
# # Resum:
# # Tenim 1840 pacients, 1 dels quals no te ID. Tenim 1537 pacients amb dades longitudinals, 1099 amb informacio longitudinal de resposta i 1115 amb mesures del tac. Dels 1115 amb responstes, 1 no 
# # te info longitudinal clinica.
# 
# # Q:
# # 1840-1099 no tenim resposta? O van ser PD primer Ct, van sortir i per tant no tenim info longitudinal? Son SF
# 
# ids_no_respons_long <- unique(clinical_db$USUBJID)[which(!unique(clinical_db$USUBJID) %in% unique(response_longitudinal$unique_subject_identifier))]
# clinical_db_weirdos <- subset(clinical_db, clinical_db$USUBJID %in% ids_no_respons_long)
# unique(clinical_db_weirdos$actual_arm)
# 
# 
# clinical_db_no_DF <- subset(clinical_db, !clinical_db$USUBJID %in% ids_no_respons_long)
# 
# # Q:
# # 1099 (clinical no SF) - 1115 measurements es perque son PD o AE sense ctSCAN? Measurements de SF
# 
# ids_no_measus <- unique(measurements_longitudinal$unique_subject_identifier)[which(!unique(measurements_longitudinal$unique_subject_identifier) %in% unique(clinical_db_no_DF$USUBJID))]
# measus_db_weirdos <- subset(measurements_longitudinal, measurements_longitudinal$unique_subject_identifier %in% ids_no_measus)
# 
# ids_no_measus %in% clinical_db_weirdos$USUBJID
# 
# # Q:
# # 1537 long clin - 1099 clin no SF, SF que tenim longit? YES
# 
# ids_clin_long <- unique(clinical_longitudinal$unique_subject_identifier)[which(!unique(clinical_longitudinal$unique_subject_identifier) %in% unique(clinical_db_no_DF$USUBJID))]
# clong_db_weirdos <- subset(clinical_longitudinal, clinical_longitudinal$unique_subject_identifier %in% ids_clin_long)
# 
# clong_db_weirdos_cdb <- clinical_db[clinical_db$USUBJID %in% ids_clin_long,]
# unique(clong_db_weirdos_cdb$actual_arm)
# 
# sfs_ids_list <- clinical_db$USUBJID[which(clinical_db$actual_arm %in% c("SCREEN FAILURE", "NOT ASSIGNED"))]
# clinical_longitudinal_no_SF <- subset(clinical_longitudinal, !clinical_longitudinal$unique_subject_identifier %in% sfs_ids_list)
# len(unique(clinical_longitudinal_no_SF$unique_subject_identifier))
# 
# 
# measurements_longitudinal_no_SF <- subset(measurements_longitudinal, !measurements_longitudinal$unique_subject_identifier %in% sfs_ids_list)
# len(unique(measurements_longitudinal_no_SF$unique_subject_identifier))
# 
# 
# response_longitudinal_no_SF <- subset(response_longitudinal, !response_longitudinal$unique_subject_identifier %in% sfs_ids_list)
# len(unique(response_longitudinal_no_SF$unique_subject_identifier))
# 
# 
# loop_assessment_visit_date_no_SF <- generate_assessment_visit_date_dataset(measurements_longitudinal_dataset = measurements_longitudinal_no_SF , response_longitudinal_dataset = response_longitudinal_no_SF)
# #########################
# 
# clinical_db$best_overall_response
# clinical_longitudinal$best
# 
# pfs(start_date=clinical_db$treatment_start_date,
#     pd_date = clinical_db$date,
#     fu_date = clinical_db$death_date,
#     eot=clinical_db$treatment_end_date,
#     cens_os=clinical_db$overall_survival_censor)
# pfs_time
# pfs_cens
