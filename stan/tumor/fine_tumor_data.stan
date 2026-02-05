int<lower = 1> n_tumor_locations;

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

