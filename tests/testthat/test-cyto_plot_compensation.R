context("cyto_plot_compensation")

# Make sure working directory checks are turned off
options("CytoRSuite_wd_check" = FALSE)

# cyto_plot_compensation flowFrame method --------------------------------------

test_that("cyto_plot_compensation flowFrame method", {
  
  p <- function() {
    cyto_plot_compensation(getData(gsc, "Single Cells")[[1]],
                           channel = "7-AAD-A",
                           compensate = FALSE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-fr1", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc, "Single Cells")[[5]],
                           channel = channel_match_file,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv",
                           display = 500,
                           contour_lines = 10)
  }
  
  expect_doppelganger("cyto_plot_compensation-fr2", p)
  
})

# cyto_plot_compensation flowSet method ----------------------------------------

test_that("cyto_plot_compensation flowSet method", {
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(3,7)], "Single Cells"),
                           channel_match = channel_match_file,
                           compensate = FALSE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-fs1", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(3,7)], "Single Cells"),
                           channel_match = channel_match_file,
                           compensate = TRUE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-fs2", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(7,5)], "Single Cells"),
                           channel_match = channel_match_file,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv",
                           display = 500,
                           contour_lines = 10)
  }
  expect_doppelganger("cyto_plot_compensation-fs3", p)
  
})

# cyto_plot_compensation GatingHierarchy method --------------------------------

test_that("cyto_plot_compensation GatingHierarchy method", {
  
  p <- function() {
    cyto_plot_compensation(gsc[[1]],
                           channel_match = channel_match_file,
                           compensate = FALSE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-gh1", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[[1]],
                           channel_match = channel_match_file,
                           compensate = TRUE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-gh2", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[[2]],
                           channel_match = channel_match_file,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv",
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-gh3", p)
  
})

# cyto_plot_compensation GatingSet method --------------------------------------

test_that("cyto_plot_compensation GatingSet method", {
  
  p <- function() {
    cyto_plot_compensation(gsc[c(1,7)],
                           channel_match = channel_match_file,
                           compensate = FALSE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-gs1", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[c(1,7)],
                           channel_match = channel_match_file,
                           compensate = TRUE,
                           display = 500)
  }
  expect_doppelganger("cyto_plot_compensation-gs2", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[c(7,5)],
                           channel_match = channel_match_file,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv",
                           diaplsy = 500)
  }
  expect_doppelganger("cyto_plot_compensation-gs3", p)
  
})

# Switch working directory checks back on
options("CytoRSuite_wd_check" = TRUE)