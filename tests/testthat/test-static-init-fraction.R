library(testthat)

# Mirror the Stan nested-logit math in R to lock the contract.
nested_split <- function(logit_loc, logit_static = NULL) {
  pi_decrease <- plogis(logit_loc)
  rest <- 1 - pi_decrease
  if (is.null(logit_static)) {
    return(c(decrease = pi_decrease, static = 0, growth = rest))
  }
  pi_static <- rest * plogis(logit_static)
  pi_growth <- rest * (1 - plogis(logit_static))
  c(decrease = pi_decrease, static = pi_static, growth = pi_growth)
}

test_that("2-way split (flag off) puts all non-decrease mass in growth", {
  p <- nested_split(0.3)
  expect_equal(unname(p["static"]), 0)
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that("3-way split is a valid simplex and recovers 2-way as static logit -> -Inf", {
  p3 <- nested_split(0.3, logit_static = 0.5)
  expect_equal(sum(p3), 1, tolerance = 1e-12)
  expect_true(all(p3 >= 0))
  p_recover <- nested_split(0.3, logit_static = -Inf)
  expect_equal(unname(p_recover["static"]), 0, tolerance = 1e-12)
  expect_equal(unname(p_recover["growth"]), unname(nested_split(0.3)["growth"]), tolerance = 1e-12)
})
