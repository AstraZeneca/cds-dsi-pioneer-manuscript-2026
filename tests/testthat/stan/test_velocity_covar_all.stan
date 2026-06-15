// Tests for standardize_velocity() and central_difference_row() (_burden.stanfunctions)
//
// standardize_velocity(v, median, iqr) = (v - median) / iqr
// central_difference_row(g): forward at col 1, central interior, backward at last col.
//   v[1]      = g[2] - g[1]
//   v[k]      = (g[k+1] - g[k-1]) / 2   for 1 < k < n
//   v[n]      = g[n] - g[n-1]
//
// Output: n_failures == 0. R test asserts n_failures == 0.

functions {
  #include "_burden.stanfunctions"
}

data {
  int<lower=0> dummy;
}

generated quantities {
  int n_failures = 0;
  real tol = 1e-10;

  // --- standardize_velocity ---
  if (abs(standardize_velocity(3.0, 1.0, 2.0) - 1.0) > tol) {
    print("FAIL standardize_velocity basic"); n_failures += 1;
  }
  if (abs(standardize_velocity(1.0, 1.0, 2.0) - 0.0) > tol) {
    print("FAIL standardize_velocity centered-at-median"); n_failures += 1;
  }

  // --- central_difference_row on a known parabola g(w) = w^2, w = 1..5 ---
  // g = [1, 4, 9, 16, 25]
  // v[1] = 4-1 = 3 (forward)
  // v[2] = (9-1)/2 = 4
  // v[3] = (16-4)/2 = 6
  // v[4] = (25-9)/2 = 8
  // v[5] = 25-16 = 9 (backward)
  {
    row_vector[5] g = [1, 4, 9, 16, 25];
    row_vector[5] v = central_difference_row(g);
    row_vector[5] expected = [3, 4, 6, 8, 9];
    for (k in 1:5) {
      if (abs(v[k] - expected[k]) > tol) {
        print("FAIL central_difference_row[", k, "]=", v[k], " expected ", expected[k]);
        n_failures += 1;
      }
    }
  }

  // --- central_difference equals analytic derivative for a quadratic (interior) ---
  // For g(w) = b0 + b1*w + b2*w^2, central diff at interior = b1 + 2*b2*w exactly.
  // b0=0, b1=1, b2=0.5 → g(w) = w + 0.5 w^2; g'(w) = 1 + w.
  // w=1..5: g = [1.5, 4, 7.5, 12, 17.5]; interior v should equal 1+w = [_,3,4,5,_].
  {
    row_vector[5] g = [1.5, 4, 7.5, 12, 17.5];
    row_vector[5] v = central_difference_row(g);
    if (abs(v[2] - 3.0) > tol) { print("FAIL quad interior v[2]=", v[2]); n_failures += 1; }
    if (abs(v[3] - 4.0) > tol) { print("FAIL quad interior v[3]=", v[3]); n_failures += 1; }
    if (abs(v[4] - 5.0) > tol) { print("FAIL quad interior v[4]=", v[4]); n_failures += 1; }
  }

  // --- single-column edge case: n=1 → velocity 0 (no neighbour) ---
  {
    row_vector[1] g = [7];
    row_vector[1] v = central_difference_row(g);
    if (abs(v[1] - 0.0) > tol) { print("FAIL n=1 v[1]=", v[1]); n_failures += 1; }
  }
}
