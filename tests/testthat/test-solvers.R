# R Solver mathematical test suite using testthat
library(testthat)

# When running standalone, source the solver files
# (adjust paths if sourced from elsewhere)
if (!exists("solve_topsis")) {
  # We assume we are running in tests/testthat/
  if (file.exists("../../R/solvers.R")) {
    source("../../R/solvers.R")
    source("../../R/monte_carlo.R")
  } else if (file.exists("R/solvers.R")) {
    source("R/solvers.R")
    source("R/monte_carlo.R")
  }
}

test_that("SMARTER ROC weights calculation", {
  weights <- calculate_roc_weights(3)
  expect_equal(length(weights), 3)
  # w_1 = 11/18 = 0.6111
  # w_2 = 5/18 = 0.2777
  # w_3 = 2/18 = 0.1111
  expect_true(abs(weights[1] - 11/18) < 1e-4)
  expect_true(abs(weights[2] - 5/18) < 1e-4)
  expect_true(abs(weights[3] - 2/18) < 1e-4)
})

test_that("Linear utility normalization", {
  mat <- matrix(c(
    10.0, 100.0,
    20.0, 50.0,
    30.0, 0.0
  ), nrow=3, ncol=2, byrow=TRUE)
  types <- c("benefit", "cost")
  norm <- normalize_linear(mat, types)
  
  expect_equal(norm[1, 1], 0.0)
  expect_equal(norm[2, 1], 0.5)
  expect_equal(norm[3, 1], 1.0)
  expect_equal(norm[1, 2], 0.0)
  expect_equal(norm[2, 2], 0.5)
  expect_equal(norm[3, 2], 1.0)
})

test_that("AHP weights using row geometric mean", {
  mat <- matrix(c(
    1.0, 3.0, 9.0,
    1/3, 1.0, 3.0,
    1/9, 1/3, 1.0
  ), nrow=3, ncol=3, byrow=TRUE)
  
  res <- calculate_geometric_mean_weights(mat)
  expect_true(abs(res$weights[1] - 9/13) < 1e-4)
  expect_true(abs(res$weights[2] - 3/13) < 1e-4)
  expect_true(abs(res$weights[3] - 1/13) < 1e-4)
  expect_true(res$cr < 0.01)
})

test_that("BWM weights optimization via lpSolve", {
  best_to_others <- c(1.0, 3.0, 9.0)
  others_to_worst <- c(9.0, 3.0, 1.0)
  res <- solve_bwm_weights_lp(3, 1, 3, best_to_others, others_to_worst)
  
  expect_true(res$success)
  expect_true(abs(res$weights[1] - 9/13) < 1e-3)
  expect_true(abs(res$weights[2] - 3/13) < 1e-3)
  expect_true(abs(res$weights[3] - 1/13) < 1e-3)
  expect_true(res$xi < 1e-4)
})

test_that("BWT piecewise linear bisection interpolation", {
  expect_equal(interpolate_bisection(10.0, 10.0, 30.0, 15.0, TRUE), 0.0)
  expect_equal(interpolate_bisection(12.5, 10.0, 30.0, 15.0, TRUE), 0.25)
  expect_equal(interpolate_bisection(15.0, 10.0, 30.0, 15.0, TRUE), 0.5)
  expect_equal(interpolate_bisection(22.5, 10.0, 30.0, 15.0, TRUE), 0.75)
  expect_equal(interpolate_bisection(30.0, 10.0, 30.0, 15.0, TRUE), 1.0)
})

test_that("MACBETH scale weights optimization via lpSolve", {
  mat <- matrix(c(
    0.0, 3.0, 5.0,
    -3.0, 0.0, 2.0,
    -5.0, -2.0, 0.0
  ), nrow=3, ncol=3, byrow=TRUE)
  res <- solve_macbeth_lp(mat)
  expect_true(res$success)
  expect_true(res$weights[1] > res$weights[2])
  expect_true(res$weights[2] > res$weights[3])
})

test_that("TOPSIS compromise resolution", {
  matrix_data <- list(
    list(10.0, 100.0),
    list(20.0, 50.0),
    list(30.0, 0.0)
  )
  types <- c("benefit", "cost")
  pref <- list(weights = list("1" = 10.0, "2" = 10.0))
  res <- solve_topsis(matrix_data, types, pref, c(1, 2), c(1, 2, 3))
  
  expect_equal(res$ranks[3], 1)
  expect_equal(res$ranks[1], 3)
})

test_that("VIKOR compromise resolution", {
  matrix_data <- list(
    list(10.0, 100.0),
    list(20.0, 50.0),
    list(30.0, 0.0)
  )
  types <- c("benefit", "cost")
  pref <- list(weights = list("1" = 10.0, "2" = 10.0), v = 0.5)
  res <- solve_vikor(matrix_data, types, pref, c(1, 2), c(1, 2, 3))
  
  expect_equal(res$ranks[3], 1)
  expect_equal(res$ranks[1], 3)
})

test_that("ELECTRE overranking kernel resolution", {
  matrix_data <- list(
    list(10.0, 100.0),
    list(20.0, 50.0),
    list(30.0, 0.0)
  )
  types <- c("benefit", "cost")
  pref <- list(weights = list("1" = 10.0, "2" = 10.0), electre_version = "I", concordance_threshold = 0.5, discordance_threshold = 0.5)
  res <- solve_electre(matrix_data, types, pref, c(1, 2), c(1, 2, 3))
  
  expect_true("kernel" %in% names(res$extra))
})

test_that("PROMETHEE net flows resolution", {
  matrix_data <- list(
    list(10.0, 100.0),
    list(20.0, 50.0),
    list(30.0, 0.0)
  )
  types = c("benefit", "cost")
  pref = list(weights = list("1" = 10.0, "2" = 10.0), promethee_version = "II")
  res <- solve_promethee(matrix_data, types, pref, c(1, 2), c(1, 2, 3))
  
  expect_equal(res$ranks[3], 1)
  expect_equal(res$ranks[1], 3)
})

test_that("Monte Carlo simulation running in R", {
  matrix_data <- list(
    list(10.0, 100.0),
    list(20.0, 50.0),
    list(30.0, 0.0)
  )
  types <- c("benefit", "cost")
  weights <- c(0.5, 0.5)
  variations_pct <- c(10.0, 10.0)
  pref <- list(weights = list("1" = 50.0, "2" = 50.0))
  
  sim <- run_monte_carlo(matrix_data, types, weights, variations_pct, 100, "topsis", pref, c(1, 2))
  expect_equal(length(sim$first_place_probabilities), 3)
  expect_equal(length(sim$average_ranks), 3)
})

print("All R testthat solver mathematical tests defined successfully!")
