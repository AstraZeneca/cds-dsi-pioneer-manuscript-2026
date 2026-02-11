# =============================================================================
# FlatIron to HISTORICAL Schema Transformation Pipeline
# =============================================================================
# 
# Purpose: Transform FlatIron RWD (pluvcohort) into HISTORICAL trial schema format
#          Produces: flatiron_patient_data and flatiron_visit_data
#
# Input:   pluvcohort.csv (from Brad's pipeline)
#          Raw FlatIron tables: visit, lab, enhanced_metpc_progression
#
# Output:  Patient-level data matching HISTORICAL schema
#          Visit-level (longitudinal) data matching HISTORICAL schema
#
# Author:  Data Science Team
# Date:    2026-02-11
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

# Source utility functions (use dirname of this script)
SCRIPT_DIR <- if (interactive()) {
  "r/data_pipelines"
} else {
  dirname(sys.frame(1)$ofile)
}
source(file.path(SCRIPT_DIR, "historical_utils.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input paths
PLUVCOHORT_PATH <- "/mnt/data/Pioneer_data/pluvcohort.csv"
LONGLABS_PATH <- "/mnt/data/Pioneer_data/pluvictolonglabs.csv"

# Study identifier for this dataset
STUDY_ID <- "FLATIRON_PLUVICTO"

# Masking window for progression attribution (FlatIron convention)
PD_MASKING_DAYS <- 14

# Data cutoff date
DATA_CUTOFF_DATE <- as.Date("2025-07-15")


# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading pluvcohort data...\n")
pluvcohort <- read.csv(PLUVCOHORT_PATH, stringsAsFactors = FALSE) %>%
  mutate(
    across(c(indexstartdt, indexenddt, deathdt, lkfdt, diagnosisdate, 
             mcrpcdt, firstdt_pluv, enddt_pluv, nextlotstartdt,
             priorlotstartdt, priorlotenddt, endobsdt), 
           ~as.Date(.x))
  )

cat(sprintf("Loaded %d patients from pluvcohort\n", nrow(pluvcohort)))


# =============================================================================
# PART 1: PATIENT-LEVEL DATA (flatiron_patient_data)
# =============================================================================
# Matches HISTORICAL historical_patient_data schema

generate_patient_data <- function(pluvcohort) {
  
  cat("Generating patient-level data...\n")
  
  # Calculate first patient start date for calendar week/day
  first_patient_start <- min(pluvcohort$indexstartdt, na.rm = TRUE)
  
  patient_data <- pluvcohort %>%
    mutate(
      # =====================================================================
      # Core Identification & Treatment Timing
      # =====================================================================
      studyid = STUDY_ID,
      usubjid = patientid,
      trtsdt = indexstartdt,
      trtedt = indexenddt,
      
      # Treatment timing calculations
      treatment_end_week = ifelse(
        !is.na(trtedt) & !is.na(trtsdt),
        as.integer((as.numeric(trtedt - trtsdt)) %/% 7 + 1),
        NA_integer_
      ),
      treatment_end_day = ifelse(
        !is.na(trtedt) & !is.na(trtsdt),
        as.integer(as.numeric(trtedt - trtsdt) + 1),
        NA_integer_
      ),
      
      # Calendar week/day relative to first patient
      calendar_week = as.integer(as.numeric(trtsdt - first_patient_start) %/% 7 + 1),
      calendar_day = as.integer(as.numeric(trtsdt - first_patient_start) + 1),
      
      # =====================================================================
      # Visit Timing
      # =====================================================================
      # patient_min_t: First post-baseline assessment week
      # For prostate cancer: based on first PSA/visit after treatment start
      patient_first_visit = trtsdt,  # Approximation: treatment start
      patient_last_visit = lkfdt,    # Last known follow-up
      
      # patient_min_t and patient_max_t will be calculated from visit data
      # For now, use approximations based on available data
      patient_min_t = 1L,  # First assessment after baseline
      patient_max_t = ifelse(
        !is.na(lkfdt) & !is.na(trtsdt),
        as.integer((as.numeric(lkfdt - trtsdt)) %/% 7 + 1),
        NA_integer_
      ),
      patient_t_width = ifelse(
        !is.na(patient_max_t) & !is.na(patient_min_t),
        patient_max_t - patient_min_t + 1L,
        NA_integer_
      ),
      
      # =====================================================================
      # Survival Endpoints
      # =====================================================================
      death = ifelse(!is.na(deathdt), TRUE, FALSE),
      death_week = ifelse(
        !is.na(deathdt) & !is.na(trtsdt),
        as.integer((as.numeric(deathdt - trtsdt)) %/% 7 + 1),
        NA_integer_
      ),
      
      # OS time in weeks (from risktime which is in months)
      os_time = round(risktime * 4.33, 1),  # Convert months to weeks
      os_cnsr = ifelse(death, 1L, 0L),  # 1=event, 0=censored
      
      # PFS endpoints (placeholder - need progression data for full calculation)
      # These will be populated from visit data or separate progression table
      pfs_cnsr = NA_integer_,
      progressor = NA,
      progress_week = NA_integer_,
      pfs = NA_integer_,
      interval_censored = NA_integer_,
      right_censored = !death,
      
      # =====================================================================
      # Demographics & Baseline Characteristics
      # =====================================================================
      age = ageatindex,
      age_group = case_when(
        age < 18 ~ "<18",
        age >= 18 & age < 40 ~ "18-40",
        age >= 40 & age < 65 ~ "40-65",
        age >= 65 & age < 75 ~ "65-75",
        age >= 75 ~ ">75",
        TRUE ~ NA_character_
      ),
      sex = "M",  # Prostate cancer cohort is male-only
      race = as.character(racecat),
      
      # =====================================================================
      # Clinical Characteristics
      # =====================================================================
      ecogbl = as.integer(ecogvalue_index),
      histology = histology,
      stage = denovo,  # De novo vs Recurrent as staging proxy
      arm = indexregcat,
      
      # =====================================================================
      # Metastatic Sites
      # =====================================================================
      liver_mets = as.character(met_liver),
      brain_mets = as.character(met_brain),
      bone_mets = as.character(met_bone),
      lung_mets = as.character(met_lung),
      lymph_mets = as.character(met_lymph),
      visceral_mets = as.character(visceralmet),
      
      # =====================================================================
      # Treatment History (Prostate Cancer Specific)
      # =====================================================================
      prev_lines = indexlinenum - 1L,
      chemo_flag = ifelse(numpriortaxaneregs > 0, 1L, 0L),
      chemo_naive = ifelse(numpriortaxaneregs == 0, 1L, 0L),
      first_line = ifelse(indexlinenum == 1, 1L, 0L),
      prior_arpi = ifelse(numpriorarpiregs > 0, 1L, 0L),
      prior_taxane_count = numpriortaxaneregs,
      prior_arpi_count = numpriorarpiregs,
      
      # =====================================================================
      # Baseline Laboratory Values (from pluvcohort)
      # Using FlatIron labs, not forcing HISTORICAL labs
      # =====================================================================
      baseline_albumin = blrslt_albumin,
      baseline_alp = blrslt_alp,
      baseline_alt = blrslt_alanine,
      baseline_ast = blrslt_aspar,
      baseline_bilirubin = blrslt_bili,
      baseline_gfr = blrslt_gfr,
      baseline_hemoglobin = blrslt_heme,
      baseline_ldh = blrslt_ldh,
      baseline_leukocytes = blrslt_leuk,
      baseline_neutrophils = blrslt_neut,
      baseline_platelets = blrslt_plat,
      baseline_psa = blrslt_psa,
      
      # NLR approximation (using leukocytes as proxy for lymphocytes ratio)
      baseline_nlr = ifelse(
        !is.na(blrslt_neut) & !is.na(blrslt_leuk) & blrslt_leuk > 0,
        blrslt_neut / (blrslt_leuk - blrslt_neut),  # Approximate
        NA_real_
      ),
      
      # =====================================================================
      # QC Flags
      # =====================================================================
      imputed = 0L  # No imputation in this transformation
      
    ) %>%
    
    # Select and order columns to match HISTORICAL schema
    select(
      # Core identification
      studyid, usubjid, trtsdt, trtedt,
      treatment_end_week, treatment_end_day,
      calendar_week, calendar_day,
      
      # Visit timing
      patient_min_t, patient_max_t, patient_t_width,
      patient_first_visit, patient_last_visit,
      
      # Survival endpoints
      death, death_week, os_time, os_cnsr,
      pfs_cnsr, progressor, progress_week, pfs,
      interval_censored, right_censored,
      
      # Demographics
      age, age_group, sex, race,
      
      # Clinical characteristics
      ecogbl, histology, stage, arm,
      
      # Metastatic sites
      liver_mets, brain_mets, bone_mets, lung_mets, lymph_mets, visceral_mets,
      
      # Treatment history
      prev_lines, chemo_flag, chemo_naive, first_line,
      prior_arpi, prior_taxane_count, prior_arpi_count,
      
      # Labs
      baseline_albumin, baseline_alp, baseline_alt, baseline_ast,
      baseline_bilirubin, baseline_gfr, baseline_hemoglobin, baseline_ldh,
      baseline_leukocytes, baseline_neutrophils, baseline_platelets,
      baseline_psa, baseline_nlr,
      
      # QC
      imputed
    )
  
  cat(sprintf("Generated patient data for %d patients\n", nrow(patient_data)))
  
  return(patient_data)
}


# =============================================================================
# PART 2: VISIT-LEVEL (LONGITUDINAL) DATA (flatiron_visit_data)
# =============================================================================
# Matches HISTORICAL historical_visit_data schema
# Uses PSA assessments as the "tumor assessment" equivalent for prostate cancer

generate_visit_data <- function(pluvcohort, psa_data = NULL, progression_data = NULL) {
  
  cat("Generating visit-level data...\n")
  
  # If PSA data not provided, create placeholder from pluvcohort
  if (is.null(psa_data)) {
    cat("Note: No longitudinal PSA data provided. Creating minimal visit data from pluvcohort.\n")
    cat("To get full longitudinal data, provide psa_data from lab table.\n")
    
    # Create minimal visit data with baseline and available info
    visit_data <- pluvcohort %>%
      select(
        patientid, indexstartdt, lkfdt, endobsdt,
        blrslt_psa
      ) %>%
      # Create baseline visit
      mutate(
        studyid = STUDY_ID,
        usubjid = patientid,
        visitnum = 0L,
        visit_date = indexstartdt,
        ady = 0L,
        day = 0L,
        week = 0L,
        psa_value = blrslt_psa,
        psa_pct_change = 0,  # Baseline is 0% change
        response = "BL",  # Baseline
        is_baseline = TRUE
      ) %>%
      select(
        studyid, usubjid, visitnum, visit_date, ady, day, week,
        psa_value, psa_pct_change, response, is_baseline
      )
    
  } else {
    # Full longitudinal PSA data processing
    cat("Processing longitudinal PSA data...\n")
    
    # Ensure date columns
    psa_data <- psa_data %>%
      mutate(testdate = as.Date(testdate))
    
    # Join with pluvcohort to get index dates and baseline PSA
    visit_data <- psa_data %>%
      inner_join(
        pluvcohort %>% 
          select(patientid, indexstartdt, endobsdt, blrslt_psa),
        by = "patientid"
      ) %>%
      # Filter to post-index PSA values only (within observation period)
      filter(testdate >= indexstartdt & testdate <= endobsdt) %>%
      # Calculate study day, week, and percent change
      mutate(
        studyid = STUDY_ID,
        usubjid = patientid,
        visit_date = testdate,
        ady = as.integer(as.numeric(testdate - indexstartdt) + 1),
        day = ady,
        week = as.integer((ady - 1) %/% 7 + 1),
        psa_value = testresultcleaned,
        psa_pct_change = ifelse(
          !is.na(blrslt_psa) & blrslt_psa > 0,
          100 * (psa_value - blrslt_psa) / blrslt_psa,
          NA_real_
        ),
        is_baseline = (ady == 1)
      ) %>%
      # Assign visit numbers
      arrange(usubjid, visit_date) %>%
      group_by(usubjid) %>%
      mutate(visitnum = row_number() - 1L) %>%  # 0-indexed for baseline
      ungroup() %>%
      # Derive PSA response category
      mutate(
        response = case_when(
          is_baseline ~ "BL",
          psa_pct_change <= -90 ~ "PSA90",  # >=90% decline
          psa_pct_change <= -50 ~ "PSA50",  # >=50% decline
          psa_pct_change < 0 ~ "DECLINE",   # Any decline
          psa_pct_change >= 25 ~ "RISE",    # >=25% rise (PSA progression criteria)
          TRUE ~ "STABLE"
        )
      ) %>%
      select(
        studyid, usubjid, visitnum, visit_date, ady, day, week,
        psa_value, psa_pct_change, response, is_baseline
      )
  }
  
  # Add progression events if provided
  if (!is.null(progression_data)) {
    cat("Adding progression events to visit data...\n")
    
    progression_data <- progression_data %>%
      mutate(progressiondate = as.Date(progressiondate))
    
    pd_visits <- progression_data %>%
      inner_join(
        pluvcohort %>% select(patientid, indexstartdt),
        by = "patientid"
      ) %>%
      # Apply +14 day masking
      filter(progressiondate >= indexstartdt + PD_MASKING_DAYS) %>%
      # Take first PD per patient
      arrange(patientid, progressiondate) %>%
      group_by(patientid) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(
        studyid = STUDY_ID,
        usubjid = patientid,
        visitnum = 999L,  # Special visit number for PD
        visit_date = progressiondate,
        ady = as.integer(as.numeric(progressiondate - indexstartdt) + 1),
        day = ady,
        week = as.integer((ady - 1) %/% 7 + 1),
        psa_value = NA_real_,
        psa_pct_change = NA_real_,
        response = "PD",
        is_baseline = FALSE
      ) %>%
      select(
        studyid, usubjid, visitnum, visit_date, ady, day, week,
        psa_value, psa_pct_change, response, is_baseline
      )
    
    # Combine PSA visits with PD events
    visit_data <- bind_rows(visit_data, pd_visits) %>%
      arrange(usubjid, ady)
  }
  
  cat(sprintf("Generated visit data: %d visits for %d patients\n", 
              nrow(visit_data), n_distinct(visit_data$usubjid)))
  
  return(visit_data)
}


# =============================================================================
# PART 2b: VISIT-LEVEL DATA FROM LONGITUDINAL LABS (pluvictolonglabs.csv)
# =============================================================================
# Uses the actual longitudinal lab data including PSA measurements

generate_visit_data_from_longlabs <- function(pluvcohort, psa_data, longlabs) {
  
  cat("Generating visit-level data from longitudinal labs...\n")
  
  # Get index dates for each patient (ensure unique)
  patient_index <- pluvcohort %>%
    select(patientid, indexstartdt, endobsdt, lkfdt) %>%
    distinct(patientid, .keep_all = TRUE)
  
  # Process PSA data (already extracted in main function)
  visit_data <- psa_data %>%
    inner_join(patient_index, by = "patientid", relationship = "many-to-one") %>%
    mutate(
      studyid = STUDY_ID,
      usubjid = patientid,
      visit_date = testdate,
      ady = as.integer(as.numeric(testdate - indexstartdt) + 1),
      day = ady,
      # Use the week from longlabs (already calculated as weektime)
      is_baseline = (week == 0),
      
      # Derive PSA response category based on percent change
      response = case_when(
        is_baseline ~ "BL",
        psa_pct_change <= -90 ~ "PSA90",  # >=90% decline (deep response)
        psa_pct_change <= -50 ~ "PSA50",  # >=50% decline (confirmed response)
        psa_pct_change < 0 ~ "DECLINE",    # Any decline
        psa_pct_change >= 25 ~ "RISE",     # >=25% rise (PSA progression)
        TRUE ~ "STABLE"
      )
    ) %>%
    # Assign visit numbers per patient
    arrange(usubjid, visit_date) %>%
    group_by(usubjid) %>%
    mutate(visitnum = row_number() - 1L) %>%  # 0-indexed for baseline
    ungroup() %>%
    # Select final columns
    select(
      studyid, usubjid, visitnum, visit_date, ady, day, week,
      psa_value, psa_pct_change, response, is_baseline,
      baseline_psa, bestchngcat, worstchngcat
    )
  
  # Add patients without PSA data (baseline only from pluvcohort)
  patients_with_psa <- unique(visit_data$usubjid)
  patients_without_psa <- pluvcohort %>%
    filter(!patientid %in% patients_with_psa) %>%
    select(patientid, indexstartdt, blrslt_psa) %>%
    mutate(
      studyid = STUDY_ID,
      usubjid = patientid,
      visitnum = 0L,
      visit_date = indexstartdt,
      ady = 0L,
      day = 0L,
      week = 0L,
      psa_value = blrslt_psa,
      psa_pct_change = 0,
      response = "BL",
      is_baseline = TRUE,
      baseline_psa = blrslt_psa,
      bestchngcat = NA_character_,
      worstchngcat = NA_character_
    ) %>%
    select(
      studyid, usubjid, visitnum, visit_date, ady, day, week,
      psa_value, psa_pct_change, response, is_baseline,
      baseline_psa, bestchngcat, worstchngcat
    )
  
  visit_data <- bind_rows(visit_data, patients_without_psa) %>%
    arrange(usubjid, ady)
  
  cat(sprintf("Generated visit data: %d visits for %d patients\n", 
              nrow(visit_data), n_distinct(visit_data$usubjid)))
  cat(sprintf("  Patients with longitudinal PSA: %d\n", length(patients_with_psa)))
  cat(sprintf("  Patients with baseline only: %d\n", nrow(patients_without_psa)))
  
  # Response summary
  cat("\nPSA Response distribution:\n")
  print(table(visit_data$response, useNA = "ifany"))
  
  return(visit_data)
}


# =============================================================================
# PART 3: UPDATE PATIENT DATA WITH PFS FROM VISIT DATA
# =============================================================================

update_patient_pfs <- function(patient_data, visit_data) {
  
  cat("Updating patient PFS from visit data...\n")
  
  # For prostate cancer: RISE (>=25% PSA increase) = PSA progression
  # This is analogous to PD in RECIST for solid tumors
  
  # Calculate PFS metrics from visit data
  pfs_summary <- visit_data %>%
    group_by(usubjid) %>%
    summarize(
      # Progression info - include both "PD" (if present) and "RISE" (PSA progression)
      has_pd = any(response %in% c("PD", "RISE"), na.rm = TRUE),
      progress_week_new = ifelse(
        has_pd, 
        min(week[response %in% c("PD", "RISE")], na.rm = TRUE), 
        NA_integer_
      ),
      
      # Last assessment before PD (lower bound for interval censoring)
      last_clean_week = ifelse(
        has_pd,
        max(week[!response %in% c("PD", "RISE") & week < progress_week_new], na.rm = TRUE),
        max(week, na.rm = TRUE)  # For censored patients, use last assessment
      ),
      
      # Interval censored gap
      interval_censored_new = ifelse(
        has_pd & !is.na(progress_week_new) & !is.infinite(last_clean_week),
        progress_week_new - last_clean_week - 1L,
        NA_integer_
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      # Handle infinite values
      last_clean_week = ifelse(is.infinite(last_clean_week), NA_integer_, as.integer(last_clean_week))
    )
  
  # Remove old PFS columns from patient_data before joining
  cols_to_remove <- intersect(
    c("progressor", "progress_week", "pfs", "interval_censored", "pfs_cnsr", "right_censored"),
    names(patient_data)
  )
  
  # Update patient data with PFS info
  patient_data_updated <- patient_data %>%
    select(-any_of(cols_to_remove)) %>%
    left_join(pfs_summary, by = "usubjid") %>%
    mutate(
      progressor = coalesce(has_pd, FALSE),
      progress_week = progress_week_new,
      pfs = last_clean_week,  # Lower bound
      interval_censored = interval_censored_new,
      pfs_cnsr = ifelse(progressor, 1L, 0L),  # 1=event, 0=censored
      right_censored = !progressor & !death
    ) %>%
    select(-c(has_pd, progress_week_new, last_clean_week, interval_censored_new))
  
  cat("Patient PFS data updated\n")
  
  return(patient_data_updated)
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

run_flatiron_to_historical_pipeline <- function(
    pluvcohort_path = PLUVCOHORT_PATH,
    longlabs_path = LONGLABS_PATH,
    output_dir = "/mnt/data/processed/"
) {
  
  cat("=============================================================\n")
  cat("FlatIron to HISTORICAL Schema Transformation Pipeline\n")
  cat("=============================================================\n\n")
  
  # Load pluvcohort
  cat("Step 1: Loading pluvcohort data...\n")
  pluvcohort <- read.csv(pluvcohort_path, stringsAsFactors = FALSE) %>%
    mutate(
      across(c(indexstartdt, indexenddt, deathdt, lkfdt, diagnosisdate, 
               mcrpcdt, firstdt_pluv, enddt_pluv, nextlotstartdt,
               priorlotstartdt, priorlotenddt, endobsdt), 
             ~as.Date(.x))
    )
  cat(sprintf("  Loaded %d patients\n", nrow(pluvcohort)))
  
  # Load longitudinal labs
  cat("\nStep 1b: Loading longitudinal labs data...\n")
  longlabs <- read.csv(longlabs_path, stringsAsFactors = FALSE) %>%
    mutate(
      indexdate = as.Date(indexdate),
      labdt = as.Date(labdt)
    )
  cat(sprintf("  Loaded %d lab records for %d patients\n", 
              nrow(longlabs), n_distinct(longlabs$patientid)))
  
  # Extract PSA data from longlabs
  psa_data <- longlabs %>%
    filter(labtype == "Prostate Specific Antigen (ng/mL)") %>%
    select(patientid, labdt, weektime, blrslt, labrslt, pctchng, bestchngcat, worstchngcat) %>%
    rename(
      testdate = labdt,
      week = weektime,
      baseline_psa = blrslt,
      psa_value = labrslt,
      psa_pct_change = pctchng
    )
  cat(sprintf("  Extracted %d PSA measurements for %d patients\n", 
              nrow(psa_data), n_distinct(psa_data$patientid)))
  
  # Generate patient-level data
  cat("\nStep 2: Generating patient-level data...\n")
  flatiron_patient_data <- generate_patient_data(pluvcohort)
  
  # Generate visit-level data using PSA longitudinal data
  cat("\nStep 3: Generating visit-level data...\n")
  flatiron_visit_data <- generate_visit_data_from_longlabs(pluvcohort, psa_data, longlabs)
  
  # Update patient PFS from visit data
  cat("\nStep 4: Updating patient PFS metrics...\n")
  flatiron_patient_data <- update_patient_pfs(flatiron_patient_data, flatiron_visit_data)
  
  # Save outputs
  cat("\nStep 5: Saving outputs...\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(flatiron_patient_data, 
            file.path(output_dir, "flatiron_patient_data.csv"),
            row.names = FALSE)
  
  write.csv(flatiron_visit_data,
            file.path(output_dir, "flatiron_visit_data.csv"),
            row.names = FALSE)
  
  # Summary
  cat("\n=============================================================\n")
  cat("Pipeline Complete!\n")
  cat("=============================================================\n")
  cat(sprintf("Patient-level data: %d patients\n", nrow(flatiron_patient_data)))
  cat(sprintf("Visit-level data:   %d visits, %d patients\n", 
              nrow(flatiron_visit_data), n_distinct(flatiron_visit_data$usubjid)))
  cat(sprintf("\nOutputs saved to: %s\n", output_dir))
  cat("  - flatiron_patient_data.csv\n")
  cat("  - flatiron_visit_data.csv\n")
  
  # Return as list
  return(list(
    patient_data = flatiron_patient_data,
    visit_data = flatiron_visit_data
  ))
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# Basic usage (with just pluvcohort):
# results <- run_flatiron_to_historical_pipeline()

# With full longitudinal data (if you have access to raw FlatIron tables):
# 
# # Load PSA data from lab table
# psa_data <- lab %>%
#   filter(testbasename == 'prostate specific antigen' & !is.na(testdate)) %>%
#   select(patientid, testdate, testresult, testresultcleaned)
#
# # Load progression data
# progression_data <- pfsdat %>%
#   mutate(progressiondate = case_when(
#     progressiondategranularity == 'Month' ~ 
#       as.Date(paste(year(progressiondate), month(progressiondate), 15, sep="-")),
#     TRUE ~ progressiondate
#   )) %>%
#   select(patientid, progressiondate)
#
# # Run pipeline with full data
# results <- run_flatiron_to_historical_pipeline(
#   psa_data = psa_data,
#   progression_data = progression_data
# )

# Access results:
# patient_data <- results$patient_data
# visit_data <- results$visit_data
