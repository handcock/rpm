library(testthat)
library(rpm)
data(fauxmatching)
context("rpm.R")
test_that("rpm-fit",{
  fit <- rpm(~match("edu") + WtoM_diff("edu",3),
          Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
          X_w="X_w", Z_w="Z_w",
          pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
          sampled="sampled")
  expect_true(all(!is.na(fit$coef)))
})
