# =============================================================================
# FlatIron to HISTORICAL Schema Transformation Pipeline
# =============================================================================
# 
# Purpose: Transform FlatIron RWD into HISTORICAL trial schema format
#          Uses RAW FlatIron RData files directly
#
# Input:   Raw FlatIron .RData files from /mnt/imported/data/flatiron_test/
#          - enhanced_metpc_alphabetaemitters.RData (Pluvicto cohort)
#          - enhanced_metpc_progression.RData (PD events)
#          - lab.RData (PSA and other labs)
#          - visit.RData (visit dates)
#          - demographics.RData
#          - enhanced_metprostate.RData (disease info)
#          - lineoftherapy.RData (treatment history)
#          - enhanced_mortality_v2.RData (death data)
#
# Output:  Patient-level data matching HISTORICAL schema
#          Visit-level (longitudinal) data matching HISTORICAL schema
#
# Author:  Data Science Innovation Team
# Date:    2026-02-12
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

# Raw data directory
RAW_DATA_DIR <- "/mnt/imported/data/flatiron_test/"

# Study identifier for this dataset
STUDY_ID <- "FLATIRON_PLUVICTO"

# Pluvicto drug name in raw data
PLUVICTO_DRUGNAME <- "vipivotide tetraxetan"

# Masking window for progression attribution (FlatIron convention)
# PD events within this window of treatment start are attributed to prior line
PD_MASKING_DAYS <- 14

# Data cutoff date
DATA_CUTOFF_DATE <- as.Date("2025-07-15")


# =============================================================================
# HELPER: Load RData file into fresh environment
# =============================================================================

load_rdata <- function(filename) {
  filepath <- file.path(RAW_DATA_DIR, filename)
  env <- new.env()
  load(filepath, envir = env)
  # Return the first (and usually only) object
  obj_name <- ls(envir = env)[1]
  return(env[[obj_name]])
}


# =============================================================================
# HELPER: Apply Brad's Cohort Criteria
# =============================================================================
# Implements inclusion/exclusion criteria from Brad's d300_cohortbuilder.Rmd:
#
# For Pluvicto Initiator Population:
#   1. Has mCRPC date
#   2. Line of therapy must be specific Pluvicto regimens:
#      - Vipivotide Tetraxetan (monotherapy)
#      - Abiraterone,Vipivotide Tetraxetan
#      - Apalutamide,Vipivotide Tetraxetan
#      - Darolutamide,Vipivotide Tetraxetan
#      - Enzalutamide,Vipivotide Tetraxetan
#   3. mCRPC line number 1-6 (high lines filtered due to misclassification)
#   4. Index date year >= 2022
#   5. Consecutive Pluvicto lines collapsed together
#
# Parameters:
#   pluvcohort        - Patient cohort with Pluvicto data
#   alpha_emitters    - Raw drug administration data
#   progression_data  - Raw progression events
#   master            - enhanced_metprostate with mCRPC date
#   lot               - Raw line of therapy data (lineoftherapy.RData)
#   drugepisode       - Raw drug episode data (drugepisode.RData) - optional
#
# Returns:
#   Filtered cohort matching Brad's criteria with indexregcat column

apply_brads_cohort_criteria <- function(pluvcohort, alpha_emitters, progression_data, 
                                        master, lot, drugepisode = NULL) {
  
  cat("Applying Brad's cohort criteria (Pluvicto Initiator Population)...\n")
  n_initial <- nrow(pluvcohort)
  
  # -------------------------------------------------------------------------
  # Step 1: Get mCRPC date from master
  # -------------------------------------------------------------------------
  mcrpc_dates <- master %>%
    select(patientid, metdiagnosisdate, crpcdate) %>%
    mutate(
      metdiagnosisdate = as.Date(metdiagnosisdate),
      crpcdate = as.Date(crpcdate),
      mcrpcdt = pmax(metdiagnosisdate, crpcdate, na.rm = FALSE)
    ) %>%
    select(patientid, mcrpcdt)
  
  # -------------------------------------------------------------------------
  # Step 2: Build clean line of therapy in mCRPC setting
  # -------------------------------------------------------------------------
  cat("  Building mCRPC line of therapy data...\n")
  
  # Filter to mCRPC lines and join with mCRPC dates
  cleanlot_mcrpc <- lot %>%
    filter(linesetting == "mCRPC") %>%
    filter(!is.na(linenumber)) %>%
    arrange(patientid, startdate) %>%
    select(patientid, linenumber, linename, startdate, enddate) %>%
    mutate(
      startdate = as.Date(startdate),
      enddate = as.Date(enddate)
    )
  
  cat(sprintf("  mCRPC lines: %d records for %d patients\n", 
              nrow(cleanlot_mcrpc), n_distinct(cleanlot_mcrpc$patientid)))
  
  # -------------------------------------------------------------------------
  # Step 3: Collapse consecutive Pluvicto-containing lines
  # -------------------------------------------------------------------------
  cat("  Collapsing consecutive Pluvicto lines...\n")
  
  # Function to collapse consecutive Pluvicto lines (iterative like Brad's code)
  collapse_consecutive_pluv <- function(lot_data) {
    for (iter in 1:3) {  # Brad does 3 iterations
      lot_data <- lot_data %>%
        arrange(patientid, linenumber) %>%
        group_by(patientid) %>%
        mutate(
          prev_has_pluv = lag(grepl("vipivotide|tetraxetan", linename, ignore.case = TRUE), default = FALSE),
          curr_has_pluv = grepl("vipivotide|tetraxetan", linename, ignore.case = TRUE),
          consecutive_pluv = curr_has_pluv & prev_has_pluv
        ) %>%
        mutate(
          # If consecutive Pluvicto, merge into prior line
          linenumber = ifelse(consecutive_pluv, linenumber - 1, linenumber)
        ) %>%
        ungroup() %>%
        # Re-aggregate by new line number
        mutate(drug = strsplit(as.character(linename), ",")) %>%
        unnest(drug) %>%
        group_by(patientid, linenumber) %>%
        summarize(
          linename = paste(unique(drug), collapse = ","),
          startdate = min(startdate, na.rm = TRUE),
          enddate = max(enddate, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        # Renumber lines to be contiguous
        arrange(patientid, linenumber) %>%
        group_by(patientid) %>%
        mutate(linenumber = row_number()) %>%
        ungroup()
    }
    return(lot_data)
  }
  
  cleanlot_pluv <- collapse_consecutive_pluv(cleanlot_mcrpc)
  
  # -------------------------------------------------------------------------
  # Step 4: Filter to valid Pluvicto regimens
  # -------------------------------------------------------------------------
  cat("  Filtering to valid Pluvicto regimens...\n")
  
  # Valid Pluvicto regimens (exact line names from Brad's code)
  valid_regimens <- c(
    "Vipivotide Tetraxetan",
    "Abiraterone,Vipivotide Tetraxetan",
    "Apalutamide,Vipivotide Tetraxetan",
    "Darolutamide,Vipivotide Tetraxetan",
    "Enzalutamide,Vipivotide Tetraxetan"
  )
  
  # Also handle reversed order (Vipivotide first)
  valid_regimens_rev <- c(
    valid_regimens,
    "Vipivotide Tetraxetan,Abiraterone",
    "Vipivotide Tetraxetan,Apalutamide",
    "Vipivotide Tetraxetan,Darolutamide",
    "Vipivotide Tetraxetan,Enzalutamide"
  )
  
  pluv_lines <- cleanlot_pluv %>%
    filter(linename %in% valid_regimens_rev) %>%
    mutate(
      indexregcat = case_when(
        grepl("rone|amide", linename, ignore.case = TRUE) ~ "Pluvicto + ARPI combo",
        TRUE ~ "Pluvicto monotherapy"
      )
    )
  
  cat(sprintf("  Valid Pluvicto regimen lines: %d records for %d patients\n", 
              nrow(pluv_lines), n_distinct(pluv_lines$patientid)))
  
  # -------------------------------------------------------------------------
  # Step 5: Apply additional filters
  # -------------------------------------------------------------------------
  
  # Filter: Line number 1-6 in mCRPC
  pluv_lines <- pluv_lines %>%
    filter(between(linenumber, 1, 6))
  
  n_after_linenum <- n_distinct(pluv_lines$patientid)
  cat(sprintf("  After line 1-6 filter: %d patients\n", n_after_linenum))
  
  # Filter: Index date year >= 2022
  pluv_lines <- pluv_lines %>%
    filter(year(startdate) >= 2022)
  
  n_after_year <- n_distinct(pluv_lines$patientid)
  cat(sprintf("  After year >= 2022 filter: %d patients\n", n_after_year))
  
  # -------------------------------------------------------------------------
  # Step 6: Join with cohort data and apply mCRPC filter
  # -------------------------------------------------------------------------
  
  # Join pluvcohort with the valid Pluvicto lines
  pluvcohort_filtered <- pluvcohort %>%
    inner_join(
      pluv_lines %>%
        select(patientid, linenumber, linename, indexregcat,
               indexlinestartdt = startdate, indexlineenddt = enddate),
      by = "patientid"
    ) %>%
    left_join(mcrpc_dates, by = "patientid") %>%
    # Must have mCRPC date
    filter(!is.na(mcrpcdt))
  
  n_after_mcrpc <- n_distinct(pluvcohort_filtered$patientid)
  cat(sprintf("  After mCRPC date filter: %d patients\n", n_after_mcrpc))
  
  # Update indexstartdt to first Pluvicto dose within the line
  # (Brad uses startdt_pluv from pluvdatfixed, we'll use alpha_emitters)
  first_pluv_in_line <- alpha_emitters %>%
    filter(drugname == PLUVICTO_DRUGNAME) %>%
    filter(!is.na(administrationdate)) %>%
    select(patientid, pluv_dose_dt = administrationdate)
  
  pluvcohort_filtered <- pluvcohort_filtered %>%
    left_join(
      first_pluv_in_line %>%
        inner_join(
          pluvcohort_filtered %>% select(patientid, indexlinestartdt, indexlineenddt),
          by = "patientid"
        ) %>%
        filter(pluv_dose_dt >= indexlinestartdt & pluv_dose_dt <= indexlineenddt) %>%
        group_by(patientid, indexlinestartdt) %>%
        summarize(
          startdt_pluv = min(pluv_dose_dt),
          enddt_pluv = max(pluv_dose_dt),
          .groups = "drop"
        ),
      by = c("patientid", "indexlinestartdt")
    ) %>%
    # Update index start to first Pluvicto dose
    mutate(
      indexstartdt = coalesce(startdt_pluv, indexstartdt),
      enddt_pluv = coalesce(enddt_pluv, indexenddt)
    )
  
  # -------------------------------------------------------------------------
  # Step 7: Keep only first Pluvicto line per patient (like Brad)
  # -------------------------------------------------------------------------
  pluvcohort_filtered <- pluvcohort_filtered %>%
    arrange(patientid, linenumber) %>%
    group_by(patientid) %>%
    slice(1) %>%
    ungroup()
  
  # -------------------------------------------------------------------------
  # Step 8: Summary
  # -------------------------------------------------------------------------
  cat(sprintf("\n  Final cohort: %d patients (%.1f%% of initial %d)\n", 
              n_distinct(pluvcohort_filtered$patientid),
              100 * n_distinct(pluvcohort_filtered$patientid) / n_initial,
              n_initial))
  
  cat("\n  Regimen distribution:\n")
  regimen_dist <- table(pluvcohort_filtered$indexregcat)
  for (reg in names(regimen_dist)) {
    cat(sprintf("    %s: %d\n", reg, regimen_dist[reg]))
  }
  
  return(pluvcohort_filtered)
}


# =============================================================================
# PART 1: BUILD PLUVICTO COHORT FROM RAW FILES
# =============================================================================

build_pluvicto_cohort_from_raw <- function() {
  
  cat("Building Pluvicto cohort from raw files...\n")
  
  # Load raw data files
  cat("  Loading raw data files...\n")
  alpha_emitters <- load_rdata("enhanced_metpc_alphabetaemitters.RData")
  demographics <- load_rdata("demographics.RData")
  master <- load_rdata("enhanced_metprostate.RData")
  mortality <- load_rdata("enhanced_mortality_v2.RData")
  lot <- load_rdata("lineoftherapy.RData")
  
  # Identify Pluvicto patients (vipivotide tetraxetan = Lu-177 PSMA)
  pluvicto_admin <- alpha_emitters %>%
    filter(drugname == PLUVICTO_DRUGNAME) %>%
    mutate(administrationdate = as.Date(administrationdate)) %>%
    group_by(patientid) %>%
    summarize(
      indexstartdt = min(administrationdate, na.rm = TRUE),  # First Pluvicto dose
      indexenddt = max(administrationdate, na.rm = TRUE),    # Last Pluvicto dose
      n_doses = n(),
      .groups = "drop"
    )
  
  cat(sprintf("  Found %d Pluvicto patients\n", nrow(pluvicto_admin)))
  
  # Join with demographics
  cohort <- pluvicto_admin %>%
    left_join(demographics, by = "patientid") %>%
    left_join(master, by = "patientid") %>%
    left_join(
      mortality %>% 
        select(patientid, dateofdeath) %>%
        mutate(dateofdeath = as.Date(dateofdeath)),
      by = "patientid"
    ) %>%
    mutate(
      # Calculate derived fields
      deceased = !is.na(dateofdeath),
      deathdt = dateofdeath,
      ageatindex = as.integer(difftime(indexstartdt, birthdate, units = "days") / 365.25),
      
      # Last known follow-up (use death date or data cutoff)
      lkfdt = pmin(
        coalesce(dateofdeath, DATA_CUTOFF_DATE),
        DATA_CUTOFF_DATE,
        na.rm = TRUE
      ),
      
      # End of observation
      endobsdt = lkfdt
    )
  
  # Get prior treatment info from line of therapy
  prior_lot <- lot %>%
    mutate(
      startdate = as.Date(startdate),
      enddate = as.Date(enddate)
    )
  
  # For each Pluvicto patient, count prior lines
  cohort <- cohort %>%
    left_join(
      prior_lot %>%
        inner_join(pluvicto_admin %>% select(patientid, indexstartdt), by = "patientid") %>%
        filter(startdate < indexstartdt) %>%  # Lines before Pluvicto
        group_by(patientid) %>%
        summarize(
          priorlotnum = n(),
          .groups = "drop"
        ),
      by = "patientid"
    ) %>%
    mutate(priorlotnum = coalesce(priorlotnum, 0L))
  
  cat(sprintf("  Cohort built: %d patients\n", nrow(cohort)))
  
  return(cohort)
}


# =============================================================================
# PART 2: GENERATE PATIENT-LEVEL DATA (RAW DATA VERSION)
# =============================================================================
# Simplified version that works with raw data fields

generate_patient_data_from_raw <- function(pluvcohort) {
  
  cat("Generating patient-level data from raw data...\n")
  
  first_patient_start <- min(pluvcohort$indexstartdt, na.rm = TRUE)
  
  patient_data <- pluvcohort %>%
    mutate(
      # Core Identification
      studyid = STUDY_ID,
      usubjid = patientid,
      trtsdt = indexstartdt,
      trtedt = indexenddt,
      
      # Follow-up timing
      lkfdt = endobsdt,
      
      # Treatment timing
      treatment_end_week = ifelse(
        !is.na(trtedt) & !is.na(trtsdt),
        as.integer((as.numeric(trtedt - trtsdt)) %/% 7 + 1),
        NA_integer_
      ),
      
      # Calendar timing
      calendar_week = as.integer(as.numeric(trtsdt - first_patient_start) %/% 7 + 1),
      
      # Follow-up window
      patient_max_t = ifelse(
        !is.na(lkfdt) & !is.na(trtsdt),
        as.integer((as.numeric(lkfdt - trtsdt)) %/% 7 + 1),
        NA_integer_
      ),
      
      # Death endpoints
      death = deceased,
      death_week = ifelse(
        deceased & !is.na(deathdt) & !is.na(trtsdt),
        as.integer((as.numeric(deathdt - trtsdt)) %/% 7 + 1),
        NA_integer_
      ),
      
      # OS time in weeks
      os_time = ifelse(
        !is.na(lkfdt) & !is.na(trtsdt),
        round(as.numeric(lkfdt - trtsdt) / 7, 1),
        NA_real_
      ),
      os_cnsr = ifelse(death, 1L, 0L),
      
      # PFS placeholders (will be filled by add_real_pd)
      pfs_cnsr = NA_integer_,
      progressor = FALSE,
      progress_week = NA_integer_,
      pfs = NA_integer_,
      interval_censored = NA_integer_,
      right_censored = !death,
      
      # Demographics
      age = ageatindex,
      age_group = case_when(
        age < 65 ~ "<65",
        age >= 65 & age < 75 ~ "65-74",
        age >= 75 ~ ">=75",
        TRUE ~ NA_character_
      ),
      sex = "M",
      race = coalesce(race, "Unknown"),
      
      # Baseline PSA
      baseline_psa = blrslt_psa,
      
      # Treatment history
      prev_lines = priorlotnum,
      
      # QC flag
      imputed = FALSE,
      
      # Regimen category (from Brad's criteria)
      indexregcat = if ("indexregcat" %in% names(pluvcohort)) indexregcat else NA_character_,
      indexlinereg = if ("linename" %in% names(pluvcohort)) linename else NA_character_
    ) %>%
    select(
      studyid, usubjid, trtsdt, trtedt, lkfdt,
      treatment_end_week, calendar_week, patient_max_t,
      death, death_week, deathdt,
      os_time, os_cnsr,
      pfs_cnsr, progressor, progress_week, pfs, interval_censored, right_censored,
      age, age_group, sex, race,
      baseline_psa, prev_lines, imputed,
      indexregcat, indexlinereg, indexstartdt
    )
  
  cat(sprintf("Generated patient data for %d patients\n", nrow(patient_data)))
  
  return(patient_data)
}


# =============================================================================
# PART 2B: GENERATE PATIENT-LEVEL DATA (PROCESSED FILE VERSION)
# =============================================================================

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
# PART 4: ADD REAL PROGRESSION DATA AND CALCULATE PFS BOUNDS
# =============================================================================
# Uses enhanced_metpc_progression.RData for actual PD events
# Calculates interval-censored PFS:
#   - Upper bound (progress_week): Week of real PD from raw files
#   - Lower bound (pfs): Last assessment week BEFORE real PD date

add_real_pd <- function(patient_data, pluvcohort, progression_data, visit_data = NULL, raw_labs = NULL) {
  
  cat("Adding real progression data and calculating PFS bounds...\n")
  
  # Get index dates for each patient
  patient_index <- pluvcohort %>%
    select(patientid, indexstartdt, endobsdt) %>%
    distinct(patientid, .keep_all = TRUE)
  
  # Process progression data - filter to events after index date with +14 day masking
  pd_events <- progression_data %>%
    inner_join(patient_index, by = "patientid") %>%
    # Apply +14 day masking rule: PD must be >= 14 days after treatment start
    filter(progressiondate >= (indexstartdt + PD_MASKING_DAYS)) %>%
    # Filter to within observation period
    filter(progressiondate <= endobsdt) %>%
    # Calculate week of progression
    mutate(
      real_pd_week = as.integer((as.numeric(progressiondate - indexstartdt)) %/% 7 + 1)
    ) %>%
    # Get first PD event per patient
    group_by(patientid) %>%
    arrange(progressiondate) %>%
    slice(1) %>%
    ungroup() %>%
    select(
      patientid,
      indexstartdt,  # Keep for downstream calculations
      real_pd_date = progressiondate,
      real_pd_week,
      real_pd_granularity = progressiondategranularity,
      real_pd_radiographic = isradiographicevidence,
      real_pd_pathologic = ispathologicevidence,
      real_pd_tumor_marker = istumormarkerevidence,
      real_pd_clinical = isclinicalassessmentonly,
      real_pd_mixed = ismixedresponse
    )
  
  cat(sprintf("  Found %d patients with real PD events (after +14 day masking)\n", 
              nrow(pd_events)))
  
  # =========================================================================
  # Calculate lower bound: Last assessment BEFORE real PD date
  # =========================================================================
  
  # Load raw lab data for assessments if not provided
  if (is.null(raw_labs)) {
    cat("  Loading raw lab data for assessment dates...\n")
    raw_labs <- load_rdata("lab.RData")
  }
  
  # Load raw visit data
  cat("  Loading raw visit data...\n")
  raw_visits <- load_rdata("visit.RData")
  
  # Get all assessment dates (labs + visits) for each patient
  # Focus on PSA/tumor marker labs as they are the relevant assessments
  psa_labs <- raw_labs %>%
    filter(testbasename == "prostate specific ag") %>%
    select(patientid, assessment_date = testdate) %>%
    mutate(assessment_date = as.Date(assessment_date), assessment_type = "PSA")
  
  visits <- raw_visits %>%
    select(patientid, assessment_date = visitdate) %>%
    mutate(assessment_type = "Visit")
  
  # Combine all assessments
  all_assessments <- bind_rows(psa_labs, visits) %>%
    distinct(patientid, assessment_date, .keep_all = TRUE)
  
  cat(sprintf("  Total assessments: %d\n", nrow(all_assessments)))
  
  # For each patient with real PD, find the last assessment BEFORE PD date
  last_clean_assessments <- pd_events %>%
    select(patientid, real_pd_date, real_pd_week, indexstartdt) %>%
    left_join(all_assessments, by = "patientid", relationship = "many-to-many") %>%
    # Filter to assessments BEFORE PD date and AFTER treatment start
    filter(
      assessment_date >= indexstartdt,
      assessment_date < real_pd_date
    ) %>%
    # Get the last assessment before PD
    group_by(patientid) %>%
    arrange(desc(assessment_date)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      last_clean_date = assessment_date,
      last_clean_week = as.integer((as.numeric(last_clean_date - indexstartdt)) %/% 7 + 1),
      # Interval censored gap (weeks between last clean and PD)
      interval_gap = real_pd_week - last_clean_week - 1L
    ) %>%
    select(patientid, last_clean_date, last_clean_week, interval_gap, last_assessment_type = assessment_type)
  
  cat(sprintf("  Patients with last clean assessment before PD: %d\n", nrow(last_clean_assessments)))
  
  # Update pd_events with lower bound info
  pd_events_with_bounds <- pd_events %>%
    left_join(last_clean_assessments, by = "patientid")
  
  # =========================================================================
  # Join to patient data
  # =========================================================================
  
  patient_data_updated <- patient_data %>%
    left_join(pd_events_with_bounds, by = c("usubjid" = "patientid")) %>%
    mutate(
      # Real PD flag
      real_pd = !is.na(real_pd_date),
      # Composite progressor: either PSA rise OR real PD
      progressor_any = progressor | real_pd,
      
      # Interval-censored PFS (based on REAL PD, not PSA rise)
      # Upper bound: real PD week
      pfs_upper = real_pd_week,
      # Lower bound: last clean assessment week
      pfs_lower = last_clean_week,
      # Interval gap
      pfs_interval_gap = interval_gap
    )
  
  # Summary
  n_psa_only <- sum(patient_data_updated$progressor & !patient_data_updated$real_pd, na.rm = TRUE)
  n_real_only <- sum(!patient_data_updated$progressor & patient_data_updated$real_pd, na.rm = TRUE)
  n_both <- sum(patient_data_updated$progressor & patient_data_updated$real_pd, na.rm = TRUE)
  
  cat(sprintf("  PSA progression only: %d patients\n", n_psa_only))
  cat(sprintf("  Real PD only: %d patients\n", n_real_only))
  cat(sprintf("  Both PSA + Real PD: %d patients\n", n_both))
  cat(sprintf("  Any progression: %d patients (%.1f%%)\n", 
              sum(patient_data_updated$progressor_any, na.rm = TRUE),
              100 * mean(patient_data_updated$progressor_any, na.rm = TRUE)))
  cat(sprintf("  Patients with valid PFS bounds (lower < upper): %d\n",
              sum(!is.na(patient_data_updated$pfs_lower) & 
                  patient_data_updated$pfs_lower < patient_data_updated$pfs_upper, na.rm = TRUE)))
  
  return(patient_data_updated)
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

run_flatiron_to_historical_pipeline <- function(
    use_raw_files = TRUE,  # Use raw RData files directly
    apply_brads_criteria = TRUE,  # Apply Brad's cohort criteria (mCRPC, no PD before 3rd dose)
    pluvcohort_path = "/mnt/data/Pioneer_data/pluvcohort.csv",  # Fallback if not using raw
    output_dir = "/mnt/data/processed/"
) {
  
  cat("=============================================================\n")
  cat("FlatIron to HISTORICAL Schema Transformation Pipeline\n")
  cat("=============================================================\n\n")
  
  if (apply_brads_criteria && !use_raw_files) {
    warning("apply_brads_criteria=TRUE requires use_raw_files=TRUE. Setting apply_brads_criteria=FALSE.")
    apply_brads_criteria <- FALSE
  }
  
  if (use_raw_files) {
    cat("Using RAW FlatIron RData files from:", RAW_DATA_DIR, "\n\n")
    
    # =========================================================================
    # Step 1: Load all raw data files
    # =========================================================================
    cat("Step 1: Loading raw data files...\n")
    
    cat("  Loading enhanced_metpc_alphabetaemitters (Pluvicto)...\n")
    alpha_emitters <- load_rdata("enhanced_metpc_alphabetaemitters.RData")
    
    cat("  Loading demographics...\n")
    demographics <- load_rdata("demographics.RData")
    
    cat("  Loading enhanced_metprostate (disease info)...\n")
    master <- load_rdata("enhanced_metprostate.RData")
    
    cat("  Loading enhanced_mortality_v2 (death data)...\n")
    mortality <- load_rdata("enhanced_mortality_v2.RData")
    
    cat("  Loading lineoftherapy...\n")
    lot <- load_rdata("lineoftherapy.RData")
    
    cat("  Loading enhanced_metpc_progression (PD events)...\n")
    progression_data <- load_rdata("enhanced_metpc_progression.RData") %>%
      mutate(progressiondate = as.Date(progressiondate))
    
    cat("  Loading lab (for PSA)...\n")
    raw_labs <- load_rdata("lab.RData")
    
    cat("  Loading visit...\n")
    raw_visits <- load_rdata("visit.RData")
    
    # =========================================================================
    # Step 2: Build Pluvicto cohort from raw files
    # =========================================================================
    cat("\nStep 2: Building Pluvicto cohort...\n")
    
    # Identify Pluvicto patients (vipivotide tetraxetan = Lu-177 PSMA)
    pluvicto_admin <- alpha_emitters %>%
      filter(drugname == PLUVICTO_DRUGNAME) %>%
      filter(!is.na(administrationdate)) %>%
      group_by(patientid) %>%
      summarize(
        indexstartdt = min(administrationdate),  # First Pluvicto dose
        indexenddt = max(administrationdate),    # Last Pluvicto dose
        n_doses = n(),
        .groups = "drop"
      )
    
    cat(sprintf("  Found %d Pluvicto patients\n", nrow(pluvicto_admin)))
    
    # Helper function to parse month-level dates (YYYY-MM format to day 15)
    parse_month_date <- function(x) {
      ifelse(is.na(x) | x == "", NA_character_,
             paste0(x, "-15")) %>%
        as.Date()
    }
    
    # Build cohort with demographics, disease info, mortality
    pluvcohort <- pluvicto_admin %>%
      left_join(demographics, by = "patientid") %>%
      left_join(master, by = "patientid") %>%
      left_join(
        mortality %>% 
          select(patientid, dateofdeath) %>%
          mutate(dateofdeath = parse_month_date(dateofdeath)),  # YYYY-MM format
        by = "patientid"
      ) %>%
      mutate(
        deceased = !is.na(dateofdeath),
        deathdt = dateofdeath,
        ageatindex = as.integer(year(indexstartdt) - birthyear),  # Approximate age from birth year
        lkfdt = pmin(coalesce(dateofdeath, DATA_CUTOFF_DATE), DATA_CUTOFF_DATE, na.rm = TRUE),
        endobsdt = lkfdt
      )
    
    # Add prior treatment info
    prior_lot <- lot %>% mutate(startdate = as.Date(startdate), enddate = as.Date(enddate))
    
    pluvcohort <- pluvcohort %>%
      left_join(
        prior_lot %>%
          inner_join(pluvicto_admin %>% select(patientid, indexstartdt), by = "patientid") %>%
          filter(startdate < indexstartdt) %>%
          group_by(patientid) %>%
          summarize(priorlotnum = n(), .groups = "drop"),
        by = "patientid"
      ) %>%
      mutate(priorlotnum = coalesce(priorlotnum, 0L))
    
    cat(sprintf("  Cohort built: %d patients\n", nrow(pluvcohort)))
    
    # =========================================================================
    # Step 2b: Apply Brad's cohort criteria (if enabled)
    # =========================================================================
    if (apply_brads_criteria) {
      cat("\nStep 2b: Applying Brad's cohort criteria...\n")
      pluvcohort <- apply_brads_cohort_criteria(
        pluvcohort, 
        alpha_emitters, 
        progression_data, 
        master,
        lot  # Pass line of therapy data for regimen filtering
      )
    }
    
    # =========================================================================
    # Step 3: Get baseline PSA from raw labs
    # =========================================================================
    cat("\nStep 3: Extracting baseline labs...\n")
    
    # Get PSA labs (test name in FlatIron is "prostate specific ag")
    psa_labs <- raw_labs %>%
      filter(testbasename == "prostate specific ag") %>%
      filter(!is.na(testdate))
    
    cat(sprintf("  Total PSA records: %d\n", nrow(psa_labs)))
    
    # Get baseline PSA (closest to or before index date)
    # Use pluvcohort's indexstartdt which may have been updated by Brad's criteria
    baseline_psa <- psa_labs %>%
      inner_join(pluvcohort %>% select(patientid, indexstartdt), by = "patientid") %>%
      filter(testdate <= indexstartdt + 7) %>%  # Within 7 days after index
      group_by(patientid) %>%
      arrange(abs(as.numeric(testdate - indexstartdt))) %>%
      slice(1) %>%
      ungroup() %>%
      select(patientid, blrslt_psa = testresultcleaned)
    
    pluvcohort <- pluvcohort %>%
      left_join(baseline_psa, by = "patientid")
    
    cat(sprintf("  Patients with baseline PSA: %d\n", sum(!is.na(pluvcohort$blrslt_psa))))
    
  } else {
    # Fallback: Load processed pluvcohort.csv
    cat("Using processed pluvcohort.csv file\n\n")
    cat("Step 1: Loading pluvcohort data...\n")
    pluvcohort <- read.csv(pluvcohort_path, stringsAsFactors = FALSE) %>%
      mutate(
        across(c(indexstartdt, indexenddt, deathdt, lkfdt, diagnosisdate, 
                 mcrpcdt, firstdt_pluv, enddt_pluv, nextlotstartdt,
                 priorlotstartdt, priorlotenddt, endobsdt), 
               ~as.Date(.x))
      )
    cat(sprintf("  Loaded %d patients\n", nrow(pluvcohort)))
    
    # Load progression data
    cat("\nStep 2: Loading progression data...\n")
    progression_data <- load_rdata("enhanced_metpc_progression.RData") %>%
      mutate(progressiondate = as.Date(progressiondate))
    cat(sprintf("  Loaded %d progression events\n", nrow(progression_data)))
    
    raw_labs <- NULL  # Will be loaded in add_real_pd
  }
  
  # =========================================================================
  # Step 4: Generate patient-level data
  # =========================================================================
  cat("\nStep 4: Generating patient-level data...\n")
  if (use_raw_files) {
    flatiron_patient_data <- generate_patient_data_from_raw(pluvcohort)
  } else {
    flatiron_patient_data <- generate_patient_data(pluvcohort)
  }
  
  # =========================================================================
  # Step 5: Generate visit-level data from raw PSA labs
  # =========================================================================
  cat("\nStep 5: Generating visit-level data from raw labs...\n")
  
  if (use_raw_files) {
    # Use raw PSA labs directly
    flatiron_visit_data <- psa_labs %>%
      inner_join(
        pluvcohort %>% select(patientid, indexstartdt, endobsdt, blrslt_psa),
        by = "patientid"
      ) %>%
      filter(testdate >= indexstartdt & testdate <= endobsdt) %>%
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
        is_baseline = (week <= 1)
      ) %>%
      arrange(usubjid, visit_date) %>%
      group_by(usubjid) %>%
      mutate(visitnum = row_number() - 1L) %>%
      ungroup() %>%
      mutate(
        response = case_when(
          is_baseline ~ "BL",
          psa_pct_change <= -90 ~ "PSA90",
          psa_pct_change <= -50 ~ "PSA50",
          psa_pct_change < 0 ~ "DECLINE",
          psa_pct_change >= 25 ~ "RISE",
          TRUE ~ "STABLE"
        )
      ) %>%
      select(
        studyid, usubjid, visitnum, visit_date, ady, day, week,
        psa_value, psa_pct_change, response, is_baseline, blrslt_psa
      )
    
    cat(sprintf("  Generated %d PSA visits for %d patients\n", 
                nrow(flatiron_visit_data), n_distinct(flatiron_visit_data$usubjid)))
    
  } else {
    # Use processed longlabs.csv
    longlabs <- read.csv("/mnt/data/Pioneer_data/pluvictolonglabs.csv", stringsAsFactors = FALSE) %>%
      mutate(indexdate = as.Date(indexdate), labdt = as.Date(labdt))
    
    psa_data <- longlabs %>%
      filter(labtype == "Prostate Specific Antigen (ng/mL)") %>%
      select(patientid, labdt, weektime, blrslt, labrslt, pctchng, bestchngcat, worstchngcat) %>%
      rename(testdate = labdt, week = weektime, baseline_psa = blrslt, psa_value = labrslt, psa_pct_change = pctchng)
    
    flatiron_visit_data <- generate_visit_data_from_longlabs(pluvcohort, psa_data, longlabs)
  }
  
  # =========================================================================
  # Step 6: Update patient PFS from visit data (PSA-based progression)
  # =========================================================================
  cat("\nStep 6: Updating patient PFS metrics (PSA-based)...\n")
  flatiron_patient_data <- update_patient_pfs(flatiron_patient_data, flatiron_visit_data)
  
  # =========================================================================
  # Step 7: Add REAL PD from raw progression files and calculate PFS bounds
  # =========================================================================
  cat("\nStep 7: Adding real progression data and calculating interval-censored PFS...\n")
  flatiron_patient_data <- add_real_pd(
    flatiron_patient_data, 
    pluvcohort, 
    progression_data, 
    flatiron_visit_data,
    if (use_raw_files) raw_labs else NULL
  )
  
  # =========================================================================
  # Step 8: Save outputs
  # =========================================================================
  cat("\nStep 8: Saving outputs...\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(flatiron_patient_data, 
            file.path(output_dir, "flatiron_patient_data.csv"),
            row.names = FALSE)
  
  write.csv(flatiron_visit_data,
            file.path(output_dir, "flatiron_visit_data.csv"),
            row.names = FALSE)
  
  # =========================================================================
  # Summary
  # =========================================================================
  cat("\n=============================================================\n")
  cat("Pipeline Complete!\n")
  cat("=============================================================\n")
  cat(sprintf("Patient-level data: %d patients\n", nrow(flatiron_patient_data)))
  cat(sprintf("Visit-level data:   %d visits, %d patients\n", 
              nrow(flatiron_visit_data), n_distinct(flatiron_visit_data$usubjid)))
  
  # PFS summary
  n_real_pd <- sum(flatiron_patient_data$real_pd, na.rm = TRUE)
  n_valid_bounds <- sum(!is.na(flatiron_patient_data$pfs_lower) & 
                        flatiron_patient_data$pfs_lower < flatiron_patient_data$pfs_upper, na.rm = TRUE)
  cat(sprintf("\nProgression summary:\n"))
  cat(sprintf("  Patients with real PD: %d (%.1f%%)\n", n_real_pd, 100 * n_real_pd / nrow(flatiron_patient_data)))
  cat(sprintf("  Patients with valid PFS interval bounds: %d\n", n_valid_bounds))
  
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

# Using raw files (recommended):
# results <- run_flatiron_to_historical_pipeline(use_raw_files = TRUE)

# Using processed files (fallback):
# results <- run_flatiron_to_historical_pipeline(use_raw_files = FALSE)

# Access results:
# patient_data <- results$patient_data
# visit_data <- results$visit_data
