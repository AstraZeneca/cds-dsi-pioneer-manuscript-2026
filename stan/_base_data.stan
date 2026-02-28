/*
 * Data Structure for Proportional Hazard Survival Model
 *
 * This model uses a hierarchical structure to represent patient, tumor, and measurement data
 * across multiple trials. The data is organized in ragged arrays to efficiently handle
 * varying numbers of tumors per patient and measurements per tumor.
 *
 * Study and Calendar Dates: There are two types of ways to handle time in the trial. First, most
 * commonly used, is the study time: the offset from the day/week of treatment, with all <=0 are
 * baseline periods and >0 are post-treatment intervals. Second, calendar date are with respect to
 * a single point of time, typically the earliest treatment date in the data. This is usually used
 * for managing data cuts.
 *
 * Key Dimensions:
 * - n_trials: Number of clinical trials
 * - n_patients: Total number of patients across all trials
 * - n_tumor_locations: Number of possible tumor locations
 *
 * Data Hierarchy:
 * Trials > Patients > Tumors > Measurements
 */

#include "_base_hierarchy_data.stan"

array[n_patients] int<lower = 1> n_patient_visits;

array[sum(n_patient_visits)] int t_patient_visits;

// Calendar Information. These are the days/weeks each patient started treatment relative
// to all the patients in the trials modeled.
array[n_patients] int<lower = 1> calendar_week;
array[n_patients] int<lower = 1> calendar_day;

// Time Horizon Extension. Sometimes we want to extrapolate beyond the latest visit observed
// in the data, we extend it by this number of weeks.
int<lower = 1> extend_max_all_t;

/*
 * Note on Ragged Arrays:
 *
 * This data structure uses ragged arrays to efficiently represent varying numbers
 * of measurements per patient across time. The arrays are "flattened" into 1D arrays,
 * with the hierarchical structure maintained through careful indexing.
 *
 * Longitudinal measurements (tumor SLD, PSA, etc.) are defined in separate modular
 * files (tumor/data.stan, psa/data.stan) but share this common visit schedule structure.
 */
