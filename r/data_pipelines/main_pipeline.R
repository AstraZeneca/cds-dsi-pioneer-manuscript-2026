
renv::restore()
library(here)
library(dplyr) # wrangle_data needs it
library(tidyr)
library(stringr)
library(testthat)
library(survival)
library(survminer)
library(mice)
library(rlang)

# conflicts_prefer(dplyr::filter)
# conflicts_prefer(base::match)

# i_am("r/data_preparation_pipeline/sclc/main_pipeline_data_sclc.R")

# source(here("r/data_preparation_pipeline/sclc","wrangle_data_sclc.R"))
# source(here("r/data_preparation_pipeline/generic_functions","02_data_cleaning_funcs.R"))
# source(here("r/data_preparation_pipeline/generic_functions","04_analysis_helper_funcs.R"))

directory_path <- "..//..//mnt/imported/data/flatiron_test"




#=======================================================#
# Load raw data and generate combined sclc dataset  #
#=======================================================#

# List all files with the .RDS extension
rds_files <- list.files(
  path = directory_path,
  pattern = "\\.RData$",
  full.names = TRUE
)

# Load each RDS file and assign it to a variable named after the file
for (file_path in rds_files) {
  load(file_path)
}

brad_directory <- "/mnt/data/Pioneer_data/"
brad_data <- read.csv(paste0(brad_directory, "pluvcohort.csv")) 

# =================#
#    BRAD data    #
# =================#

# Patient data

# studyid - 
# usubjid - patientid
# trtedt - 
# treatment_end_week +
# treatment_end_day
# calendar_week *
# calendar_day *
# patient_min_t*
# patient_max_t*
# patient_first_visit

colnames(brad_data)[grep("id", colnames(brad_data), ignore.case = TRUE)]



# ================#
#    RAW Data    #
# ================#


# Patient data

# studyid - NONE
# usubjid - patientid (all tables)
# trtsdt -  treatment_visit_starts (generated below) Notice we have no indicator this is actually a start date. It's simply the first visit we have
# trtedt - treatment visit ends (generated below), As above, we are using the last visit (treatment) we have. 
# treatment_end_week +
# treatment_end_day
# calendar_week *
# calendar_day *
# patient_min_t*
# patient_max_t*
# patient_first_visit

# Search for columns matching a pattern across all loaded data frames
search_pattern <- "treatment"  # Change this to search for different patterns

for(file_path in rds_files) {
    obj_name <- tools::file_path_sans_ext(basename(file_path))
    obj <- get(obj_name)
    matching_cols <- grep(search_pattern, names(obj), ignore.case = TRUE, value = TRUE)
    cat("\n[", obj_name, "]\n")
    print(matching_cols)
}

treatment_visit_starts <- visit %>%
  group_by(patientid) %>%
  filter(istreatmentvisit==T) %>%
  arrange(visitdate) %>%
  slice(1) %>%
  select(patientid, visitdate)

treatment_visit_ends <- visit %>%
  group_by(patientid) %>%
  filter(istreatmentvisit == T) %>%
  arrange(visitdate) %>%
  slice(n()) %>%
  select(patientid, visitdate)
