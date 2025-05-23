# The following is just a script detailing the issues found during the data validation phase,
# the purpose of this is to be able to document and track patients excluded in the future
# This script can be safely ignored. A the bottom there's a comment with the actions that were taken to solve
# this. Comments throughout the script should guide the rationale behind what's going on.

# THIS IS NOT A DOCUMENT DETAILING THE PATIENTS TO BE INCLUDED/EXCLUDED, ITS JUST A TOOL
# IN PLACE OF SAID DOCUMENTATION WHILE THAT IS NOT READY.


## Issue 1: calendar_week and calendar_day are 0

weird_user1 <- historical_patient_data |>
  filter(calendar_week <= 0) |>
  select(usubjid)

fpsd <- sort(clinical_db_no_SF$treatment_start_date)[1] # First patient start date
first_patient_id <- clinical_db_no_SF$subjid[clinical_db_no_SF$treatment_start_date==fpsd]
first_patient_id == weird_user1$usubjid 

## Formula for calendar week is: (treatment_start_date - first_patient_sd)-1)%/%7+1

(as.numeric(fpsd - fpsd)-1)%/%7+1 # This is 0!

as.numeric(fpsd - fpsd) # In days is also 0

# We will use this format from now on for calendar_week and calendar_day so we dont have 0s

calendar_to_study_date <- function(date, treatment_start_date, unit = "weeks") {
  floor(time_length(date - treatment_start_date, unit = unit)) + 1
}


## Issue 2: patient_max_t is na sometimes # triple-check these have only bl, filter them. 

historical_patient_data <- patient_data
historical_visit_data <- assessment_visit_data
historical_patient_data |> filter(is.na(patient_max_t)) %>% semi_join(historical_visit_data, ., by = "usubjid") # we check which patients have patient_max_t as NA
# and we look at their assessment_data

ids_with_max_t_NA <-historical_patient_data$usubjid[is.na(historical_patient_data$patient_max_t)] # max_t_NA people
len(unique(ids_with_max_t_NA)) # 93 people

visits_not_BL <-historical_visit_data[historical_visit_data$visitnum!="1",] # we get the visits that are not the BL
visits_not_BL[visits_not_BL$usubjid %in% ids_with_max_t_NA,] # we get the people that are both NA in patient_max_t, and that are in the "visits_not_BL"
# So, we are digging people out with NA as patient_max_t and that have more than just a baseline visit.
# 1 comes out. With 3 visits not BL. Still, we dont have ady or week for this so...

# Since ady2 = date_time_of_tumor_measurement - baseline_visit 

weird_issue2_pt_id <- as.character(unique(visits_not_BL$usubjid[visits_not_BL$usubjid %in% ids_with_max_t_NA])) # we get the weird pt ID


View(measurements_longitudinal[measurements_longitudinal$subject_identifier_for_the_study==weird_issue2_pt_id,c("date_time_of_tumor_measurement", "baseline_record_flag", "visit_number")])
# So, this one has visit_number 1, however it is not flagged as Baseline. As such we filter it out for now


colSums(is.na(historical_visit_data))

clinical_db_no_SF[clinical_db_no_SF$subjid %in% ids_with_R_NA,c("treatment_start_date", "treatment_end_date")]

clinical_db_no_SF %>%
  filter(subjid %in% ids_with_R_NA) %>%
  select(treatment_start_date, treatment_end_date)


ids_with_R_NA %>%
  filter(is.na(usubjid))

colSums(is.na(ids_with_R_NA))

measurements_longitudinal_no_SF %>%
  filter(subject_identifier_for_the_study %in% ids_with_R_NA) %>%
  distinct(subject_identifier_for_the_study, visit_number)

measurements_longitudinal_no_SF$subject_identifier_for_the_study

ids_with_max_t_NA %in% ids_with_R_NA
#It happens only for those patients without a response.

# Filter those out.

### Sub-issue ADY and ADY2 ## Not baseline, use treatment_start_date
## Snippet from the "claude_opt_generate_assessment_visit_date_dplyr" function
## in wrangle data
measurements_longitudinal_dataset <- measurements_longitudinal_no_SF

# We get the baseline visit for every patient ## Not baseline, use treatment_start_date

baseline_ids <- measurements_longitudinal_dataset %>%
  group_by(subject_identifier_for_the_study) %>%
  filter(baseline_record_flag=="Y") %>%
  select(date_time_of_tumor_measurement)

baseline_visits_by_id <- baseline_ids %>%
  group_by(subject_identifier_for_the_study) %>%
  summarise(
    baseline_visit = first(sort(date_time_of_tumor_measurement))
  )

# Since patients will be repeated in the OG df, we will join them so that
# we can vectorize the calculation

mldj <- full_join(baseline_visits_by_id, measurements_longitudinal_dataset,
                  by="subject_identifier_for_the_study")

pre_avd <- mldj %>%
  mutate(
    usubjid = str_split_fixed(unique_subject_identifier, "/", 2)[,2], # I'm dumb. This is just the patient ID.
    ady2 = as.numeric(date_time_of_tumor_measurement - baseline_visit),  ## Not baseline, use treatment_start_date
    week = as.numeric(date_time_of_tumor_measurement - baseline_visit-1) %/% 7 +1
  ) %>%
  select(
    studyid = study_identifier,
    usubjid,
    visitnum = visit_number,
    ady = study_day_of_tumor_measurement,
    ady2 = ady2,
    #study_day_of_tumor_measurement,
    week,
    mmsumdiam = numeric_result_finding_in_standard_units,
    tatn = tumor_assessment_test_name,
    accepted_record_flag = accepted_record_flag
  ) %>%
  filter(tatn == "Sum of Diameter", accepted_record_flag == "Y")
# Notice ady and ady2 are not the same as ady is not taking baseline into account. However, I'm using this later to 
# merge datasets and have a reference point for the measurements (hence why I keep it)

# Going back the the issue. What does this have to do with the is.na(patient_max_t)?

# The data generated on the previous step is used in the next one. So the week calculation is
# the one we used as.numeric(date_time_of_tumor_measurement - baseline_visit-1) %/% 7 +1
clinical_dataset <- clinical_db_no_SF
assessment_visit_data <- historical_visit_data
assessment_visit_data <- calc_visit_date(assessment_visit_data) # READY!

first_patient_sd <- sort(clinical_dataset$treatment_start_date)[1] # we get the first visit of the first patient

week_filtered2 <- assessment_visit_data %>% 
  group_by(usubjid) %>%
  filter(week>1) %>%
  summarise(
    lweek = last(week), # idem with the last one
    .groups = "drop"
  )

# lweek is the patient_max_t

sum(is.na(week_filtered2$lweek))

# Issue 3. Interval censored is NA (and not 0) # Default to 1 if there's a ctSCAN or modify the code (ideally)

historical_patient_data |> filter(is.na(interval_censored))

# Let's create a dataset that will contain baseline date and all visits date in a "wide" format, one patient per row
# This way we can see what can be the issue. interval_censored is the span between
# the lower and upper bounds of the interval... Well, what if this patients lack one of those bounds?

interval_weirdos <- historical_patient_data |> filter(is.na(interval_censored))

response_longitudinal_no_SF |>
  filter(subject_identifier_for_the_study %in% interval_weirdos$usubjid) |>
  # filter(baseline_flag!="Y") |>
  select(date_of_visit, subject_identifier_for_the_study)

measurements_unique <- measurements_longitudinal_no_SF %>%
  group_by(subject_identifier_for_the_study, visit_number) %>%
  summarise(date_time_of_tumor_measurement = first(date_time_of_tumor_measurement), .groups = 'drop')

# Ensure unique combinations in `assessment_visit_data`
assessment_unique <- assessment_visit_data %>%
  group_by(usubjid, visitnum) %>%
  summarise(across(everything(), first), .groups = 'drop')

# Join the datasets on unique keys
complete_data <- assessment_unique %>%
  left_join(
    measurements_unique,
    by = c("usubjid" = "subject_identifier_for_the_study", "visitnum" = "visit_number")
  )

# Pivot wider only using the actual dates
pivoted_data <- complete_data %>%
  group_by(usubjid) %>%
  arrange(date_time_of_tumor_measurement) %>%
  filter(!is.na(date_time_of_tumor_measurement)) %>%
  mutate(visit_number = row_number()) %>% # Re-determine visit order based on date
  pivot_wider(
    names_from = visit_number,
    values_from = date_time_of_tumor_measurement,
    names_prefix = "visit_"
  ) %>%
  mutate(baseline_date = visit_1) %>% # Assume first visit as baseline
  select(usubjid, baseline_date, starts_with("visit_")) %>%
  ungroup()

View(pivoted_data[pivoted_data$usubjid %in% interval_weirdos$usubjid,])

# Check for each condition
results <- pivoted_data %>%
  mutate(
    # Does the patient have a baseline?
    has_baseline = !is.na(baseline_date),
    
    # Is baseline the only visit they have?
    only_baseline_visit = rowSums(!is.na(select(., starts_with("visit_")))) == 1, 
    
    # Fetch patient information like right_censored status
    right_censored_status = map_dbl(usubjid, ~ patient_data$right_censored[patient_data$usubjid == .]),
    
    # Check if there is a progression before the cut-off
    pd_before_cutoff = !right_censored_status,
    
    # Determine if a PD exists without a baseline
    pd_no_baseline = pd_before_cutoff & !has_baseline
  )

# Select columns to review patients for each condition
review_results <- results %>%
  select(usubjid, has_baseline, only_baseline_visit, pd_before_cutoff, pd_no_baseline)

aggregated_results_by_usubjid <- review_results |>
  group_by(usubjid) |>
  count(has_baseline & pd_before_cutoff)

with(aggregated_results_by_usubjid,len(unique(usubjid[`has_baseline & pd_before_cutoff`==T]))) # Usable patients


aggregated_results_by_usubjid2 <- review_results |>
  group_by(usubjid) |>
  count(pd_no_baseline)

with(aggregated_results_by_usubjid2,len(unique(usubjid[pd_no_baseline==T]))) # Not usable patients

sum(unique(aggregated_results_by_usubjid$usubjid[aggregated_results_by_usubjid$`has_baseline & pd_before_cutoff`==T]) %in%
      unique(aggregated_results_by_usubjid2$usubjid[aggregated_results_by_usubjid2$pd_no_baseline==T]))

ntable2(review_results$has_baseline, ifelse(review_results$pd_no_baseline==TRUE, "true", "false"))
# So, 294 visits have no Baseline, and no PD (since they can't have 2 FALSES and have no baseline but PD without baseline)
# 715 visits have baseline
# 2618 visits have pd without a baseline
# Last one is a QC, it has to be 0 as again, we cant have both being true at the same time by definition.

inner_join(clinical_db_no_SF, measurements_longitudinal_no_SF, by=c("subjid" = "subject_identifier_for_the_study")) |>
  group_by(subjid) |>
  filter(min_rank(as.Date(date_of_visit))==n()) |>
  ungroup()|>
  select(treatment_end_date, date_of_visit, subjid) |>
  arrange(as.Date(date_of_visit))

inner_join(clinical_db_no_SF, measurements_longitudinal_no_SF, by=c("subjid" = "subject_identifier_for_the_study")) |>
  filter(!is.na(date_of_visit)) |>
  group_by(subjid) |>
  arrange(as.Date(date_of_visit)) |>
  summarise(
    dov = last(date_of_visit),
    treatment_end_date = treatment_end_date[1],
    .groups = "drop"
  ) |>
  select(treatment_end_date, dov, subjid)


measurements_longitudinal_no_SF$date_of_visit[measurements_longitudinal_no_SF$subject_identifier_for_the_study=="E0302003"]


# Ensure 'assessment_data' has 'usubjid', 'ady', 'response' and necessary date columns correctly initialized

# Convert dates in assessment data if required
assessment_unique <- assessment_visit_data %>%
  mutate(visit_date = as.Date(ady, origin = "1970-01-01")) %>%  # Assuming `ady` needs conversion
  select(usubjid, visit_date, response) %>%
  group_by(usubjid, response) %>%
  summarise(visit_date = first(visit_date), .groups = 'drop')

# Merge response data into the processed pivoted data
complete_data <- assessment_unique %>%
  left_join(
    pivoted_data %>%
      mutate(across(starts_with("visit_"), as.Date)),  # Ensure all visit columns are Date type
    by = "usubjid"
  )

# Implement the logic to check for PD and other conditions
results <- complete_data %>%
  rowwise() %>%
  mutate(
    # Does the patient have a baseline?
    has_baseline = !is.na(baseline_date),
    
    # Is baseline the only visit they have?
    only_baseline_visit = sum(!is.na(c_across(starts_with("visit_")))) == 1,
    
    # Find the earliest PD date leveraging the merged responses
    pd_date = min(visit_date[response == "PD"], na.rm = TRUE),
    
    # Assuming `patient_data` still holds right_censored info
    right_censored_status = patient_data %>%
      filter(usubjid == cur_data()$usubjid) %>%
      pull(right_censored),
    
    # Ensure these checks accurately consider PD and status
    pd_before_cutoff = !right_censored_status & !is.na(pd_date),
    
    # Determine if a PD exists without a baseline
    pd_no_baseline = !has_baseline & !is.na(pd_date)
  ) %>%
  ungroup()

# Select columns to review patients for each condition
review_results <- results %>%
  select(usubjid, has_baseline, only_baseline_visit, pd_date, pd_before_cutoff, pd_no_baseline)

# Output the results
print(review_results)



historical_visit_data |> count(usubjid, visitnum) |> filter(n > 1) %>% semi_join(historical_visit_data, ., by = "usubjid")



# Actions DONE:

# Issue 1: We will use this format from now on for calendar_week and calendar_day so we dont have 0s

# calendar_to_study_date <- function(date, treatment_start_date, unit = "weeks") {
#   floor(time_length(date - treatment_start_date, unit = unit)) + 1
# }

## DONE wrangle_data > generate_patient_data_dataset

# Issue 2: patient_max_t is na sometimes # triple-check these have only bl, filter them.

## DONE, these patients do not have other visits beyond baseline, except for one, that has other visits but not a flagged baseline one.
# filter them out

# Issue 3. Interval censored is NA (and not 0) # Default to 1 if there's a ctSCAN or modify the code (ideally)
# same for PFS

# These patients have either a PD after treatment ends so no Upper bound, or no LB (because they progress too early)
# Just impute a 1 to the LB

# Drop NAs caused by the full_join. Add the actual dates of the visitnum.