context("cyto_plot_profile")

# cyto_plot_profile flowFrame method -------------------------------------------
test_that("cyto_plot_profile flowFrame method", {
  
  p <- function() {
    cyto_plot_profile(getData(gs, "T Cells")[[4]],
                      axes_trans = trans)
  }
  expect_doppelganger("cyto_plot_profile-fr1", p)
  
})

# cyto_plot_profile flowSet method ---------------------------------------------
test_that("cyto_plot_profile flowSet method", {
  
  p <- function() {
    cyto_plot_profile(getData(gs, "T Cells"),
                      group_by = "all",
                      axes_trans = trans)
  }
  expect_doppelganger("cyto_plot_profile-fs1", p)
  
})

# cyto_plot_profile GatingSet method -------------------------------------------

test_that("cyto_plot_profile GatingSet method", {
  
  # Parent missing
  expect_error(cyto_plot_profile(gs), 
               "Please supply the name of the parent population to plot.")
  
  p <- function() {
    cyto_plot_profile(gs,
                      parent = "T Cells")
  }
  expect_doppelganger("cyto_plot_profile-gs1", p)
  
})