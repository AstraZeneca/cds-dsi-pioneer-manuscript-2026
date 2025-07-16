# Stan Function Testing Pattern: All-in-One R/Stan Tests

This repository uses a robust, unified pattern for testing Stan functions (e.g., `cutoff_visits`, `fine_cutoff_visits`, `get_oos_patients_idx`). All test cases for a function are run in a single Stan call, with all data prepared in R and all outputs checked in R.

## How to Add a New Stan Test

### 1. Stan Model File

- Place your Stan test model in `tests/testthat/stan/` (e.g., `test_myfunction_all.stan`).
- Include all required Stan files at the top (`#include "../lfo.stan"`, etc).
- The model should loop over all test cases (vectorized if possible).
- Output all results as rectangular arrays (pad with zeros if needed).
- Use clear, unique output variable names.

### 2. R Test File

- Place your R test script in `tests/testthat/` (e.g., `test-stan-myfunction.R`).
- Inline all test cases as a single R list, or source them from a separate file if very large.
- Each test case should be a named list with all required inputs and expected outputs.
- Prepare a single rectangular data list for Stan, padding arrays as needed.
- Use the provided helper (`test_stan_function`) to run the Stan model (this always recompiles and surfaces Stan errors).
- Extract all outputs for each case and check them with `testthat::expect_equal()` or similar.
- Always check edge cases and boundary conditions.

### 3. Naming and Structure

- Stan test models: `stan/test_<function>_all.stan`
- R test scripts: `test-stan-<function>.R`
- Inline all test cases in the R file unless they are very large (then source them).
- Do not use legacy per-case files (e.g., `cutoff_visits_cases.R`); all cases should be inlined or sourced as a single list.
- Only keep one test file per function; delete or archive old/fragmented test files.

### 4. Adding Test Cases

- Add new cases directly to the unified list in the R test file.
- Each case should include all necessary input data and expected outputs.
- For ragged arrays, use position arrays and pad as needed.
- Always provide at least one element for each array (Stan does not allow zero-length arrays).
- For new edge cases, add them to the list and update the expected outputs.

### 5. Debugging and Best Practices

- If a test fails, check the Stan output, the R data prep, and the expected values.
- Use diagnostic `cat()` prints in R to compare expected and actual outputs.
- If you hit Stan errors, run the model standalone with the generated JSON to surface the true error.
- Always review previous all-in-one test files for working patterns.

### 6. Example Files

- Stan: `tests/testthat/stan/test_lfo_all.stan`
- R: `tests/testthat/test-stan-lfo.R`
- Helper: `tests/testthat/helper-stan.R`

### 7. What to Avoid

- Do not use per-case R files or Stan files.
- Do not keep unused JSON or helper files.
- Do not skip zero-patient cases by logic; instead, always provide at least one patient.

---

**Summary Table**

| What to do                        | How/Where                                      |
|------------------------------------|------------------------------------------------|
| Add Stan test model                | `tests/testthat/stan/test_<function>_all.stan`  |
| Add R test script                  | `tests/testthat/test-stan-<function>.R`        |
| Inline all test cases              | As a single list in the R file                  |
| Prepare rectangular Stan data      | In the R test script                            |
| Use helper for Stan runs           | `test_stan_function()` in R                     |
| Check all outputs                  | `expect_equal()` in R                           |
| Add new cases                      | Directly to the unified list                    |
| Remove unused/legacy files         | Clean up as you go                              |

---

**Replicate this pattern for all new Stan function tests to ensure robust, maintainable, and debuggable validation.**

