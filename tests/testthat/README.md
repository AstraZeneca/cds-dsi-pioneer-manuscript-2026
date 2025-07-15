# Stan Function Testing Pattern for Survival Analysis (Kaplan-Meier Example)

This document describes the robust testing pattern used for validating Stan survival analysis functions (e.g., `estimate_kaplan_meier`) in this repository. Follow this pattern to add new tests for other Stan functions.

## Design Pattern Overview

- **All test cases for a function are run from a single Stan model file.**
  - This ensures that any change to the Stan function or its includes will force recompilation and all scenarios are tested together.
- **Test data is defined in R as a list of cases.**
  - Each test case is a list with all required inputs (e.g., event times, censoring, max_t, offsets).
  - Zero-patient cases are excluded (all cases must have at least one patient).
- **R test script prepares a single data list for Stan.**
  - Data is packed into rectangular arrays/matrices for Stan input.
  - Indices are tracked so each test case's results can be extracted from the Stan output.
- **A single Stan test model runs all test cases in a vectorized way.**
  - The Stan model loops over all cases, calling the function under test for each, and outputs results as rectangular arrays.
- **R test script extracts and checks results for each case.**
  - Output arrays are parsed per-case.
  - All relevant properties (length, monotonicity, at-risk counts, event/censor counts, etc.) are checked for each case.
  - Output is compared to R reference (e.g., `survfit2`) where possible.
- **No skip logic for zero-patient cases.**
  - All test cases are valid and robust.
- **Test helper always forces Stan recompilation and prints compiler errors.**


## What You Need to Add a New Stan Function Test

Before adding a new Stan function test and creating a new all-in-one Stan testing model, gather the following information:

- **Stan function signature**: Name, argument types, and expected output type/shape (scalar, vector, matrix, etc.).
- **Test cases**: For each, the full set of input values and the expected output (including edge cases).
- **R reference implementation**: If available, for expected results and debugging.
- **Stan includes and dependencies**: List of all required Stan files to include (e.g., utility, math, or domain-specific functions). 
- **Output conventions**: How the Stan model should output results (rectangular arrays, padding, sentinel values, etc.).
- **Naming/location conventions**: Where to place the Stan test model and R test file (see previous files for pattern).
- **Edge case handling**: Any special handling for zero-length, NA, or boundary cases.

### Checklist for Creating a New Test

1. **Review previous all-in-one Stan test model files** (e.g., `stan/test_estimate_kaplan_meier_all.stan`, `stan/test_calculate_target_recist_all.stan`) for conventions:
    - How includes are handled (relative paths, required dependencies)
    - How data is declared and passed in
    - How outputs are structured (rectangular arrays, padding)
    - Any required blocks (functions, data, generated quantities)
2. **Create a new all-in-one Stan test model** for your function, following the conventions above.
3. **Add a new R test file** in this directory, using the pattern in the existing all-in-one test files. Define a list of test cases, each with all required inputs and expected outputs.
4. **Update helper-stan.R** if new helper logic is needed for your function.
5. **Run the test file** with `testthat::test_file()` and debug as needed. If you hit errors, check the Stan model includes, output structure, and compare to previous working models.

**Important:**
- Always review previous all-in-one Stan test models before starting a new one. Many issues (include paths, output shape, block structure) are solved by copying the working pattern.
- If you hit errors, compare your new model to a previous one line by line. Most problems are due to missing includes, non-rectangular output, or block structure mismatches.

| What you need                | Why it's needed                                    |
|-----------------------------|----------------------------------------------------|
| Stan function signature     | To know what to test and how to call it            |
| Test cases (inputs/outputs) | To check correctness and edge case handling        |
| R reference implementation  | For expected results and debugging                 |
| Stan includes/dependencies  | To avoid missing function errors                   |
| Output conventions          | For robust R-side extraction and comparison        |
| Naming/location conventions | To keep the test suite organized                   |
| Edge case handling          | To ensure robustness and avoid silent failures     |
| Review of previous models   | To avoid common pitfalls and follow conventions    |

---

1. **Create a new Stan test model** in `tests/testthat/stan/` (e.g., `test_myfunction_all.stan`).
   - Include all necessary Stan files.
   - Loop over all test cases and output results as rectangular arrays.
2. **Write an R test script** in `tests/testthat/` (e.g., `test-stan-myfunction-<signature>.R`).
   - Define a list of test cases (no zero-patient cases).
   - Prepare a single data list for Stan.
   - Use the test helper to run the Stan model (forces recompilation, prints errors).
   - Extract and check all relevant outputs for each case.
3. **Follow the conventions in the Kaplan-Meier test** for naming, structure, and checks.

## Example Files

- Stan model: `tests/testthat/stan/test_estimate_kaplan_meier_all.stan`
- R test: `tests/testthat/test-stan-estimate_kaplan_meier-array_int-array_int-int-int.R`
- Helper: `tests/testthat/helper-stan.R`

## Key Conventions

- All test cases in one Stan model per function.
- Rectangular output arrays for robust R-side extraction.
- No zero-patient cases or unnecessary skip logic.
- All relevant properties and edge cases are checked.
- Output matches R reference where possible.

---

**Replicate this pattern for all new Stan function tests to ensure robust, maintainable, and debuggable validation.**
