/*
 * Data for Full Biomarker Models (tumor/PSA state-space).
 *
 * Requires in scope (include before this file):
 *   _hierarchy_data.stan  — n_trials, n_patients, patient_trial, hierarchy
 *   _visit_data.stan      — n_patient_visits, t_patient_visits
 *
 * Study and Calendar Dates: There are two types of ways to handle time in the
 * trial. First, most commonly used, is the study time: the offset from the
 * day/week of treatment, with all <=0 are baseline periods and >0 are
 * post-treatment intervals. Second, calendar date are with respect to a single
 * point of time, typically the earliest treatment date in the data. This is
 * usually used for managing data cuts.
 *
 * Longitudinal measurements (tumor SLD, PSA, etc.) are defined in separate
 * modular files (tumor/data.stan, psa/data.stan) but share this common visit
 * schedule structure.
 */

// HMC/Laplace routing (full models only; standalones define n_forecast_patients as constant)
int<lower=0> n_forecast_patients;
int<lower=0> forecast_split_level;
int<lower=0> forecast_group;

array[sum(n_patient_visits)] int t_patient_visits_day; // Study days (for Day 82 gate and LFO calendar-day cutoffs)

// Calendar Information. These are the days/weeks each patient started treatment
// relative to all the patients in the trials modeled.
array[n_patients] int<lower = 1> calendar_week;
array[n_patients] int<lower = 1> calendar_day;

// Time Horizon Extension. Sometimes we want to extrapolate beyond the latest
// visit observed in the data, we extend it by this number of weeks.
int<lower = 1> extend_max_all_t;

// Measurement noise degrees of freedom (Student-t). Shared across all biomarkers.
// Use ~5 for robust noise, positive_infinity() → Gaussian.
real<lower=2> measure_nu;
