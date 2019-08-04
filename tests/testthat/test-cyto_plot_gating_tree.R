context("cyto_plot_gating_tree")

# CYTO_PLOT_GATING_TREE --------------------------------------------------------

test_that("cyto_plot_gating_tree", {
  
  # gatingTemplate -------------------------------------------------------------
  
  expect_error(cyto_plot_gating_tree(gt), NA)
  
  # GatingSet ------------------------------------------------------------------
  
  expect_error(cyto_plot_gating_tree(gs), NA)
  
  # GatingHierarchy ------------------------------------------------------------
  
  expect_error(cyto_plot_gating_tree(gs[[32]], stat = "percent"), NA)
  expect_error(cyto_plot_gating_tree(gs[[32]], stat = "count"), NA)
  
})
