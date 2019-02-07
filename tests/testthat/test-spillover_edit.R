context("spillover-edit")
# This file is for testing the applications in the apps/ directory.

library(shinytest)

# spillover_edit ---------------------------------------------------------------

test_that("spillover_edit", {
  skip("Skip spillover_edit test")
  shinytest::expect_pass(testApp("apps/spillover_edit/", compareImages = FALSE))
})