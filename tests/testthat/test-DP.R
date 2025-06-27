# --- Comprehensive Unit Tests for the dp() Function ---

set.seed(42)
test_data <- simulate_dp_data(n_per_group = 100, L = 15, K = 3, interval = 5)
test_dat <- test_data$dat
test_datA <- test_data$datA


# --- Unit Tests for All Combinations of dp() ---

testthat::test_that("dp runs for link='linear' with stabilized weights", {
   result <- dp(
      dat = test_dat, datA = test_datA, L = 15, K = 3, interval = 5,
      method = "link", link = "linear", weights = "stabilized"
   )
   expect_true(is.list(result)); expect_named(result, c("betahat", "se"))
   expect_true(is.numeric(result$betahat)); expect_equal(length(result$betahat), length(result$se))
})

testthat::test_that("dp runs for link='linear' with unstabilized weights", {
   result <- dp(
      dat = test_dat, datA = test_datA, L = 15, K = 3, interval = 5,
      method = "link", link = "linear", weights = "unstabilized"
   )
   expect_true(is.list(result)); expect_named(result, c("betahat", "se"))
})

testthat::test_that("dp runs for link='log' with stabilized weights", {
   result <- dp(
      dat = test_dat, datA = test_datA, L = 15, K = 3, interval = 5,
      method = "link", link = "log", weights = "stabilized"
   )
   expect_true(is.list(result)); expect_named(result, c("betahat", "se"))
})

testthat::test_that("dp runs for link='log' with unstabilized weights", {
   result <- dp(
      dat = test_dat, datA = test_datA, L = 15, K = 3, interval = 5,
      method = "link", link = "log", weights = "unstabilized"
   )
   expect_true(is.list(result)); expect_named(result, c("betahat", "se"))
})

testthat::test_that("dp runs for stratified='add' with stabilized weights", {
   # This test assumes the stratified method is implemented in the dp function
   result <- dp(
      dat = test_dat, datA = test_datA, L = 15, K = 3, interval = 5,
      method = "stratified", stratified = "add", weights = "stabilized"
   )
   expect_true(is.list(result))
   expect_true(is.na(result$betahat) || is.numeric(result$betahat))
})

testthat::test_that("dp runs for stratified='multi' with unstabilized weights", {
   # This test assumes the stratified method is implemented in the dp function
   result <- dp(
      dat = test_dat, datA = test_datA, L = 15, K = 3, interval = 5,
      method = "stratified", stratified = "multi", weights = "unstabilized"
   )
   expect_true(is.list(result))
   expect_true(is.na(result$betahat) || is.numeric(result$betahat))
})

test_that("dp returns NA when no subjects are left after stacking", {
   bad_datA <- test_datA
   bad_datA$A <- 0
   expect_warning(
      result <- dp(
         dat = test_dat, datA = bad_datA, L = 15, K = 3, interval = 5,
         method = "link", link = "linear"
      ),
      "No subjects at risk"
   )
   expect_true(is.na(result$betahat))
})

testthat::test_that("dp returns NA when no treatment events exist", {
   bad_dat <- test_dat
   bad_dat$deltaT <- 0
   expect_warning(
      result <- dp(
         dat = bad_dat, datA = test_datA, L = 15, K = 3, interval = 5,
         method = "link", link = "linear"
      ),
      "No treatment events"
   )
   expect_true(is.na(result$betahat))
})

