context("cyto_spillover_compute")

# CYTO_SPILLOVER_COMPUTE -------------------------------------------------------

# UNIVERSAL UNSTAINED CONTROL
test_that("cyto_spillover_compute universal reference", {
  
  # GATES
  mock_locator <- mock(list("x" = c(1740, 4265),
                            "y" = c(50,50)),
                       list("x" = c(1800, 3100),
                            "y" = c(50,50)),
                       list("x" = c(1710, 2950),
                            "y" = c(50,50)),
                       list("x" = c(2720, 4115),
                            "y" = c(50,50)),
                       list("x" = c(2310, 3780),
                            "y" = c(50,50)),
                       list("x" = c(2350, 4120),
                            "y" = c(50, 50)),
                       cycle = TRUE)
  
  # SPILL
  SPILL <- read.csv("Reference-Universal-Spillover-Matrix.csv",
                    header = TRUE,
                    row.names = 1,
                    stringsAsFactors = FALSE)
  SPILL <- as.matrix(SPILL)
  colnames(SPILL) <- cyto_fluor_channels(gs_comp)
  rownames(SPILL) <- cyto_fluor_channels(gs_comp)
  
  # GATINGSET
  testthat::with_mock(locator = mock_locator,{
                        expect_equal(
                          cyto_spillover_compute(gs_comp,
                                                 parent = "Single Cells",
                                                 channel_match = "Reference-Compensation-Channels.csv",
                                                 spillover = "Spillover-Matrix.csv"),
                          SPILL)
                      })

  expect_true(file.exists("Spillover-Matrix.csv"))
  
  # SAVED MATRIX
  spill <- read.csv("Spillover-Matrix.csv",
                    header = TRUE,
                    row.names = 1,
                    stringsAsFactors = FALSE)
  colnames(spill) <- rownames(spill)
  spill <- as.matrix(spill)
  expect_equal(spill, SPILL)
  
})

# INTERNAL UNSTAINED REFERENCE
test_that("cyto_spillover_compute internal reference", {
  
  # CHANNELS
  mock_menu <- mock(4, 10, 11, 9, 1, 2, cycle = TRUE)
  
  # GATES
  mock_locator <- mock(list("x" = c(-850, 525),
                            "y" = c(50,50)),
                       list("x" = c(1770, 3950),
                            "y" = c(50,50)),
                       list("x" = c(125, 1080),
                            "y" = c(50,50)),
                       list("x" = c(1930, 2920),
                            "y" = c(50,50)),
                       list("x" = c(-30, 1030),
                            "y" = c(50,50)),
                       list("x" = c(1710, 2900),
                            "y" = c(50, 50)),
                       list("x" = c(185, 1700),
                            "y" = c(50,50)),
                       list("x" = c(2800, 4180),
                            "y" = c(50,50)),
                       list("x" = c(-235, 1325),
                            "y" = c(50,50)),
                       list("x" = c(2400, 3700),
                            "y" = c(50,50)),
                       list("x" = c(-125, 1495),
                            "y" = c(50,50)),
                       list("x" = c(2600, 3900),
                            "y" = c(50, 50)),
                       cycle = TRUE)
  
  # SPILL
  SPILL <- read.csv("Reference-Internal-Spillover-Matrix.csv",
                    header = TRUE,
                    row.names = 1,
                    stringsAsFactors = FALSE)
  SPILL <- as.matrix(SPILL)
  colnames(SPILL) <- cyto_fluor_channels(gs_comp)
  rownames(SPILL) <- cyto_fluor_channels(gs_comp)

  # TRANSFORMED GATINGSET
  testthat::with_mock(menu = mock_menu,
                      locator = mock_locator,{
                        expect_equal(
                          cyto_spillover_compute(gs_comp_trans[seq_len(length(gs_comp)-1)],
                                                 parent = "Single Cells",
                                                 spillover = "Spillover-Matrix.csv"),
                          SPILL)
                      })
  
  expect_true(file.exists("Spillover-Matrix.csv"))
  
  # SAVED MATRIX
  spill <- read.csv("Spillover-Matrix.csv",
                    header = TRUE,
                    row.names = 1,
                    stringsAsFactors = FALSE)
  colnames(spill) <- rownames(spill)
  spill <- as.matrix(spill)
  expect_equal(spill, SPILL)

})

# DELETE FILES
base::unlink("Spillover-Matrix.csv")
base::unlink(paste0(format(Sys.Date(), "%d%m%y"),"-Compensation-Channels.csv"))
