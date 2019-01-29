context("zzz")

# zzz ---------------------------------------------------------------

test_that("CytoRSuite loading", {
  CytoRSuite:::.onLoad()
  expect_true(getOption("CytoRSuite_interact") == interactive())
  expect_true(getOption("CytoRSuite_wd_check") == TRUE)
  expect_true(openCyto:::.isRegistered("gate_manual"))
  expect_true(openCyto:::.isRegistered("pp_gate_draw"))
  expect_true(openCyto:::.isRegistered("gate_draw"))
})