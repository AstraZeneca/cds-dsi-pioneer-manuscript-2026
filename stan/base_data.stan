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

// Trial and Patient Level Data
int<lower = 1> n_trials;
int<lower = 0> n_patients;

// Patient to Trial Mapping
array[n_patients] int<lower = 1, upper = n_trials> patient_trial;

array[n_patients] int<lower = 1> n_patient_visits;


array[sum(n_patient_visits)] int t_patient_visits;
array[sum(n_patient_visits)] int t_patient_visits_day; // Days

/*
 * Diagram for t_measure and t_day_measure:
 * 
 * [-1, 0, 6, 8, -1, 0, 6, 10, -2, 1, ...]
 *  ^         ^   ^             ^  ^
 *  |         |   |             |  |
 *  Tumor 1   |   Tumor 2       |  Tumor 4 (Patient 2)
 *  (Patient 1)   (Patient 1)   Tumor 3 (Patient 1)
 */

vector<lower = 0>[sum(n_patient_visits)] sum_tumor_size; // cm

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
 * of tumors per patient and measurements per tumor. The arrays are "flattened"
 * into 1D arrays, with the hierarchical structure maintained through careful indexing.
 *
 * For example, to access the tumor sizes for the j-th tumor of the i-th patient:
 * 1. Calculate the start index for the i-th patient's tumors
 * 2. Calculate the start index for the j-th tumor's measurements
 * 3. Use n_measures to determine how many measurements to read
 */