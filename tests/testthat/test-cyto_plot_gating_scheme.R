context("cyto_plot_gating_scheme")

# cyto_plot_gating_scheme GatingHierarchy method ----------------------------------------------

test_that("cyto_plot_gating_scheme GatingHierarchy method", {
  
  # gatingTemplate not in working directory
  expect_error(cyto_plot_gating_scheme(gs[[4]], gatingTemplate = "Test.csv"), "Supplied gatingTemplate file does not exist in current working directory.")
  
  p <- function() cyto_plot_gating_scheme(gs[[4]])
  expect_doppelganger("cyto_plot_gating_scheme-gh1", p)
  
  p <- function() cyto_plot_gating_scheme(gs[[4]],
                                          back_gate = TRUE, 
                                          gate_track = TRUE,
                                          gatingTemplate = "Activation-gatingTemplate.csv")
  expect_doppelganger("cyto_plot_gating_scheme-gh2", p)
  
})

# cyto_plot_gating_scheme GatingSet method ---------------------------------------------------

test_that("cyto_plot_gating_scheme GatingSet method", {
  
  # gatingTemplate not in working directory
  expect_error(cyto_plot_gating_scheme(gs, gatingTemplate = "Test.csv"), "Supplied gatingTemplate file does not exist in current working directory.")

  p <- function() cyto_plot_gating_scheme(gs[1], gatingTemplate = "Activation-gatingTemplate.csv")
  expect_doppelganger("cyto_plot_gating_scheme-gs1", p)
  
  p <- function() cyto_plot_gating_scheme(gs[1], back_gate = TRUE, gate_track = TRUE, group_by = "name", show_all = TRUE)
  expect_doppelganger("cyto_plot_gating_scheme-gs2", p)
  
})