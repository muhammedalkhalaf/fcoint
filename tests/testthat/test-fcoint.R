test_that("fcoint runs FADL and returns expected structure", {
  set.seed(1)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.3)
  res <- fcoint(y, x, test = "fadl", max_freq = 2)
  expect_s3_class(res, "fcoint")
  expect_equal(res$test, "fadl")
  r <- res$results[["fadl"]]
  expect_true(is.numeric(r$tstat))
  expect_true(is.numeric(r$cv5))
  expect_true(r$nobs > 0)
})

test_that("fcoint FEG returns numeric tstat", {
  set.seed(2)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.4 * x + rnorm(n, sd = 0.5)
  res <- fcoint(y, x, test = "feg", max_freq = 2)
  r <- res$results[["feg"]]
  expect_true(is.numeric(r$tstat))
})

test_that("fcoint FEG2 includes rho2 component", {
  set.seed(3)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.6 * x + rnorm(n, sd = 0.4)
  res <- fcoint(y, x, test = "feg2", max_freq = 2)
  r <- res$results[["feg2"]]
  expect_true(is.numeric(r$rho2))
})

test_that("fcoint Tsong returns CI_stat and F_stat", {
  set.seed(4)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.7 * x + rnorm(n, sd = 0.3)
  res <- fcoint(y, x, test = "tsong", max_freq = 2)
  r <- res$results[["tsong"]]
  expect_true(is.numeric(r$CI_stat))
  expect_true(is.numeric(r$F_stat))
})

test_that("fcoint all runs without error", {
  set.seed(5)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.3)
  res <- fcoint(y, x, test = "all", max_freq = 2)
  expect_equal(length(res$results), 4)
})

test_that("fcoint input validation works", {
  expect_error(fcoint(1:10, matrix(1:10), test = "fadl"))
  y <- rnorm(80); x <- matrix(rnorm(80), ncol = 1)
  expect_error(fcoint(y, x, max_freq = 0))
})

test_that("print.fcoint produces output", {
  set.seed(6)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n)
  res <- fcoint(y, x, test = "fadl", max_freq = 2)
  out <- capture.output(print(res))
  expect_true(length(out) > 0)
})
