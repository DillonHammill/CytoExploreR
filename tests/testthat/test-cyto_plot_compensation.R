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
                             getData(gsc, "Single Cells")[[4]], trns)
                           )
  }
    
  expect_doppelganger("cyto_plot_compensation-fr1", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc, "Single Cells")[[5]],
                           channel = "PE-Cy7-A",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv",
                           overlay = transform(
                             getData(gsc, "Single Cells")[[4]], trns)
                           )
  }
  
  expect_doppelganger("cyto_plot_compensation-fr2", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc, "Single Cells")[[1]],
                           channel = "7-AAD-A",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Notinwd.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-fr3", p)
  
})


# cyto_plot_compensation flowSet method ----------------------------------------

test_that("cyto_plot_compensation flowSet method", {
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(1,4)], "Single Cells"),
                           channel_match = "Compensation-channels.csv",
                           axes_trans = trns,
                           compensate = TRUE)
  }
  
  expect_doppelganger("cyto_plot_compensation-fs1", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(4,5)], "Single Cells"),
                           channel_match = "Compensation-channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-fs2", p)
  
  p <- function() {
    cyto_plot_compensation(getData(gsc[c(4,5)], "Single Cells"),
                           channel_match = "Compensation-channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Notinwd.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-fs3", p)
  
})

# cyto_plot_compensation GatingSet method --------------------------------------

test_that("cyto_plot_compensation GatingSet method", {
  
  p <- function() {
    cyto_plot_compensation(gsc[c(1,4)],
                           channel_match = "Compensation-channels.csv",
                           axes_trans = trns,
                           compensate = TRUE)
  }
  
  expect_doppelganger("cyto_plot_compensation-gs1", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[c(4,5)],
                           channel_match = "Compensation-channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Ref-Spillover-Matrix.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-gs2", p)
  
  p <- function() {
    cyto_plot_compensation(gsc[c(4,5)],
                           channel_match = "Compensation-channels.csv",
                           axes_trans = trns,
                           compensate = TRUE,
                           spillover = "Notinwd.csv")
  }
  
  expect_doppelganger("cyto_plot_compensation-gs3", p)
  
})