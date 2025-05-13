if(!require(here))install.packages("here")
library(here)

i_am("historical_data_integration/main_pipeline_data.R")


source(here("historical_data_integration", "get_data.R"))
source(here("historical_data_integration","wrangle_data.R"))
source(here("historical_data_integration","other_data.R"))

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
len(unique(measurements_longitudinal$unique_subject_identifier)) # 1115 all are included in clinical; 1 not in clinical_longitudinal

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

clinical_db_no_SF <- subset(clinical_db, !clinical_db$USUBJID %in% sfs_ids_list)

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

assessment_visit_date_no_SF <- 
  generate_assessment_visit_date_dplyr(
    measurements_longitudinal_dataset = measurements_longitudinal_no_SF,
    response_longitudinal_dataset = response_longitudinal_no_SF
    )

assessment_visit_date_no_SF$response <- ifelse(
  assessment_visit_date_no_SF$response=="",
  NA,
  assessment_visit_date_no_SF$response
  )


# QC (2)
colSums(is.na(assessment_visit_date_no_SF))
weirdos <- assessment_visit_date_no_SF %>%
  filter(is.na(response)) %>%
  filter(visitnum != "1.00" & visitnum != "1")


assessment_visit_data <- 
  assessment_visit_date_no_SF 

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


#==========================#
#     RESOLVING ISSUES     #
#==========================#


historical_patient_data <- patient_data
historical_visit_data <- assessment_visit_data
historical_patient_data |> filter(is.na(patient_max_t)) %>% semi_join(historical_visit_data, ., by = "usubjid")
# we check which patients have patient_max_t as NA
# and we look at their assessment_data

ids_with_max_t_NA <-historical_patient_data$usubjid[is.na(historical_patient_data$patient_max_t)] # guilty people
len(unique(ids_with_max_t_NA)) # 93 people

visits_not_BL <-historical_visit_data[historical_visit_data$visitnum!="1",] # we get the visits that are not the BL
visits_not_BL[visits_not_BL$usubjid %in% ids_with_max_t_NA,] 
# we get the people that are both NA in patient_max_t, and that are in the "visits_not_BL"
# So, we are digging people out with NA as patient_max_t and that have more than just a baseline visit.
# 1 comes out. With 3 visits not BL. Still, we dont have ady or week for this so...

# Since ady2 = date_time_of_tumor_measurement - baseline_visit 

weird_issue2_pt_id <- as.character(unique(visits_not_BL$usubjid[visits_not_BL$usubjid %in% ids_with_max_t_NA]))
# we get the weird pt ID

historical_patient_data <- historical_patient_data[!historical_patient_data$usubjid %in% ids_with_max_t_NA,]
historical_visit_data <- historical_visit_data[!historical_visit_data$usubjid %in% ids_with_max_t_NA,]

historical_patient_data$pfs <- ifelse((!is.na(historical_patient_data$progress_week) & is.na(historical_patient_data$pfs)), 1, historical_patient_data$pfs)
historical_patient_data$progress_week <- ifelse((is.na(historical_patient_data$progress_week) & !is.na(historical_patient_data$pfs)), historical_patient_data$patient_max_t, historical_patient_data$progress_week)

historical_patient_data$pfs <- ifelse(is.na(historical_patient_data$pfs) & is.na(historical_patient_data$progress_week), 1, historical_patient_data$pfs)
historical_patient_data$progress_week <- ifelse(is.na(historical_patient_data$pfs) & is.na(historical_patient_data$progress_week), 1, historical_patient_data$progress_week)

historical_patient_data$interval_censored <- ifelse(is.na(historical_patient_data$interval_censored), (historical_patient_data$progress_week-historical_patient_data$pfs), historical_patient_data$interval_censored)

historical_visit_data <- historical_visit_data |>
  filter(!is.na(historical_visit_data$visitnum))

#=======================#
#   DATA EXPORT         #
#=======================#

# write.csv(historical_visit_data, "../../mnt/data/PIONEER_2025_Historical_data/assessment_visit_data_290425.csv")
# write.csv(historical_patient_data, "../../mnt/data/PIONEER_2025_Historical_data/cooked_patient_data_290425.csv")
