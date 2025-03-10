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
int<lower = 1> n_tumor_locations;

// Patient to Trial Mapping
array[n_patients] int<lower = 1, upper = n_trials> patient_trial;

// Number of Tumors per Patient
array[n_patients] int<lower = 1> n_patient_tumors;

/*
 * Diagram for n_patient_tumors:
 * 
 * [3, 2, 4, 1, ...]
 *  ^  ^  ^  ^
 *  |  |  |  |
 *  Patient 1 has 3 tumors
 *     Patient 2 has 2 tumors
 *        Patient 3 has 4 tumors
 *           Patient 4 has 1 tumor
 */

// Number of Measurements per Tumor
array[sum(n_patient_tumors)] int<lower = 1> n_measures;

/*
 * Diagram for n_measures:
 * 
 * [4, 3, 5, 2, 3, 4, 3, ...]
 *  ^     ^     ^  ^
 *  |     |     |  |
 *  Tumor 1     |  Tumor 4 (Patient 2)
 *    Tumor 2   Tumor 3
 *    (Patient 1)
 */

// Tumor Location
array[sum(n_patient_tumors)] int<lower = 1, upper = n_tumor_locations> tumor_location;

// Measurement Times
array[sum(n_measures)] int t_measure; // Periods
array[sum(n_measures)] int t_day_measure; // Days

/*
 * Diagram for t_measure and t_day_measure:
 * 
 * [-1, 0, 6, 8, -1, 0, 6, 10, -2, 1, ...]
 *  ^         ^   ^             ^  ^
 *  |         |   |             |  |
 *  Tumor 1   |   Tumor 2       |  Tumor 4 (Patient 2)
 *  (Patient 1)   (Patient 1)   Tumor 3 (Patient 1)
 */

// Tumor Sizes
vector<lower = 0>[sum(n_measures)] tumor_size; // cm

/*
 * Diagram for tumor_size:
 * 
 * [2.1, 2.3, 2.0, 1.8, 3.2, 3.0, 2.8, 1.5, 1.3, ...]
 *  ^         ^    ^              ^    ^
 *  |         |    |              |    |
 *  Tumor 1   |    Tumor 2        |    Tumor 4 (Patient 2)
 *  (Patient 1)    (Patient 1)    Tumor 3 (Patient 1)
 */

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