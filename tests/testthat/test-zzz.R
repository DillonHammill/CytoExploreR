context("zzz")

# zzz ---------------------------------------------------------------

test_that("CytoExploreR loading", {
  CytoExploreR:::.onLoad()
  expect_false(getOption("locatorBell"))
  expect_equal(getOption("CytoExploreR_gatingTemplate"), NULL)
  expect_true(getOption("CytoExploreR_wd_check") == TRUE)
  expect_equal(getOption("cyto_plot_call"), NULL)
  expect_equal(getOption("cyto_plot_match"), NULL)
  expect_equal(getOption("cyto_plot_theme"), NULL)
  expect_false(getOption("cyto_plot_save"))
  expect_equal(getOption("cyto_plot_method"), NULL)
  expect_false(getOption("cyto_plot_custom"))
  expect_false(getOption("cyto_plot_grid"))
  expect_equal(getOption("cyto_plot_label_coords"), NULL)  
  expect_true(openCyto:::.isRegistered("cyto_gate_manual"))
  expect_true(openCyto:::.isRegistered("pp_cyto_gate_draw"))
  expect_true(openCyto:::.isRegistered("cyto_gate_draw"))
})