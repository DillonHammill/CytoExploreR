context("spillover-edit")
# This file is for testing the applications in the apps/ directory.

library(shinytest)

test_that("spillover_edit", {
  skip_on_travis()
  shinytest::expect_pass(testApp("apps/spillover_edit/", compareImages = FALSE))
})