context("cyto_plot_compensation")

# Transformations -
fr <- as(getData(gsc, "Single Cells"), "flowFrame")
trns <- estimateLogicle(fr, cyto_fluor_channels(fr))

# cyto_plot_compensation flowFrame method --------------------------------------

test_that("cyto_plot_compensation flowFrame method", {
  
  p <- function() {
    cyto_plot_compensation(getData(gsc, "Single Cells")[[1]],
                           channel = "7-AAD-A",
                           axes_trans = trns,
                           compensate = TRUE, 
                           overlay = transform(
                             getData(gsc, "Single Cells")[[7]], trns)
                           )
  }
    
  expect_doppelganger("cyto_plot_compensation-fr1", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc, "Single Cells")[[5]],
                           channel = "Alexa Fluor 488-A",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv",
                           overlay = transform(
                             getData(gsc, "Single Cells")[[7]], trns)
                           )
  }
  
  expect_doppelganger("cyto_plot_compensation-fr2", p)
  
  expect_error(
    cyto_plot_compensation(getData(gsc, "Single Cells")[[1]],
                           channel = "7-AAD-A",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Notinwd.csv"),
  "Notinwd.csv does not exist in this working directory.")
  
})

# cyto_plot_compensation flowSet method ----------------------------------------

test_that("cyto_plot_compensation flowSet method", {
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(3,7)], "Single Cells"),
                           channel_match = "Compensation-Channels.csv",
                           axes_trans = trns,
                           compensate = TRUE)
  }
  
  expect_doppelganger("cyto_plot_compensation-fs1", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(7,5)], "Single Cells"),
                           channel_match = "Compensation-Channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-fs2", p)
  
  expect_error(
    cyto_plot_compensation(getData(gsc[c(7,5)], "Single Cells"),
                           channel_match = "Compensation-Channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Notinwd.csv"),
    "Notinwd.csv does not exist in this working directory.")
  
  expect_doppelganger("cyto_plot_compensation-fs3", p)
  
})

# cyto_plot_compensation GatingSet method --------------------------------------

test_that("cyto_plot_compensation GatingSet method", {
  
  p <- function() {
    cyto_plot_compensation(gsc[c(1,7)],
                           channel_match = "Compensation-Channels.csv",
                           axes_trans = trns,
                           compensate = TRUE)
  }
  
  expect_doppelganger("cyto_plot_compensation-gs1", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[c(7,5)],
                           channel_match = "Compensation-Channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-gs2", p)
  
  expect_error(
    cyto_plot_compensation(gsc[c(7,5)],
                           channel_match = "Compensation-Channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Notinwd.csv"),
    "Notinwd.csv does not exist in this working directory.")
  
})