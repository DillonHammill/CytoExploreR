context("cyto_spillover_compute")

# CYTO_SPILLOVER_COMPUTE -------------------------------------------------------

# UNIVERSAL UNSTAINED CONTROL
test_that("cyto_spillover_compute universal reference", {
  
  # GATES
  mock_locator <- mock(list("x" = 1740,
                            "y" = 50),
                       list("x" = 4265,
                            "y" = 50),
                       list("x" = 1800,
                            "y" = 50),
                       list("x" = 3100,
                            "y" = 50),
                       list("x" = 1710,
                            "y" = 50),
                       list("x" = 2950,
                            "y" = 50),
                       list("x" = 2720,
                            "y" = 50),
                       list("x" = 4115,
                            "y" = 50),
                       list("x" = 2310,
                            "y" = 50),
                       list("x" = 3780,
                           "y" = 50),
                       list("x" = 2350,
                            "y" = 50),
                       list("x" = 4120,
                            "y" = 50),
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
  mock_locator <- mock(list("x" = -850,
                            "y" = 50),
                       list("x" = 525,
                            "y" = 50),
                       list("x" = 1770,
                            "y" = 50),
                       list("x" = 3950,
                            "y" = 50),
                       list("x" = 125,
                            "y" = 50),
                       list("x" = 1080,
                            "y" = 50),
                       list("x" = 1930,
                            "y" = 50),
                       list("x" = 2920,
                            "y" = 50),
                       list("x" = -30,
                            "y" = 50),
                       list("x" = 1030,
                            "y" = 50),
                       list("x" = 1710,
                            "y" = 50),
                       list("x" = 2900,
                            "y" = 50),
                       list("x" = 185,
                            "y" = 50),
                       list("x" = 1700,
                            "y" = 50),
                       list("x" = 2800,
                            "y" = 50),
                       list("x" = 4180,
                            "y" = 50),
                       list("x" = -235,
                            "y" = 50),
                       list("x" = 1325,
                            "y" = 50),
                       list("x" = 2400,
                            "y" = 50),
                       list("x" = 3700,
                            "y" = 50),
                       list("x" = -125,
                            "y" = 50),
                       list("x" = 1495,
                            "y" = 50),
                       list("x" = 2600,
                            "y" = 50),
                       list("x" = 3900,
                            "y" = 50),
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
                      locator = mock_locator, {
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
