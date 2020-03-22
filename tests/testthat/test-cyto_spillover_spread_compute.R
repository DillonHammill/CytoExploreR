context("cyto_spillover_spread_compute")

# CYTO_SPILLOVER_SPREAD_COMPUTE ------------------------------------------------

# UNIVERSAL UNSTAINED CONTROL
test_that("cyto_spillover_spread_compute universal reference", {

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

  # SPILLOVER SPREAD
  SPREAD <- read.csv("Reference-Universal-Spillover-Spread-Matrix.csv",
                     header = TRUE,
                     row.names = 1,
                     stringsAsFactors = FALSE)
  SPREAD <- as.matrix(SPREAD)
  colnames(SPREAD) <- cyto_fluor_channels(gs_comp)
  rownames(SPREAD) <- c("Alexa Fluor 488-A",
                        "PE-A",
                        "7-AAD-A",
                        "Alexa Fluor 647-A",
                        "Alexa Fluor 700-A",
                        "APC-Cy7-A")

  # GATINGSET
  testthat::with_mock(locator = mock_locator,{
    expect_equal(
      cyto_spillover_spread_compute(gs_comp,
                                    parent = "Single Cells",
                                    channel_match = "Reference-Compensation-Channels.csv",
                                    compensated = FALSE,
                                    spillover = "Reference-Universal-Spillover-Matrix.csv",
                                    spillover_spread = "Spillover-Spread-Matrix.csv"),
      SPREAD, tolerance = 0.01)
  })

  expect_true(file.exists("Spillover-Spread-Matrix.csv"))

  # SAVED MATRIX
  spread <- read.csv("Spillover-Spread-Matrix.csv",
                    header = TRUE,
                    row.names = 1,
                    stringsAsFactors = FALSE)
  spread <- as.matrix(spread)
  colnames(spread) <- cyto_fluor_channels(gs_comp)
  rownames(spread) <- c("Alexa Fluor 488-A",
                        "PE-A",
                        "7-AAD-A",
                        "Alexa Fluor 647-A",
                        "Alexa Fluor 700-A",
                        "APC-Cy7-A")

})

# INTERNAL UNSTAINED CONTROL
test_that("cyto_spillover_spread_compute internal reference", {

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

  # SPILLOVER SPREAD
  SPREAD <- read.csv("Reference-Internal-Spillover-Spread-Matrix.csv",
                     header = TRUE,
                     row.names = 1,
                     stringsAsFactors = FALSE)
  SPREAD <- as.matrix(SPREAD)
  colnames(SPREAD) <- cyto_fluor_channels(gs_comp)
  rownames(SPREAD) <- c("Alexa Fluor 488-A",
                        "PE-A",
                        "7-AAD-A",
                        "Alexa Fluor 647-A",
                        "Alexa Fluor 700-A",
                        "APC-Cy7-A")

  # GATINGSET
  testthat::with_mock(menu = mock_menu,
                      locator = mock_locator,{
    expect_equal(
      cyto_spillover_spread_compute(gs_comp[seq_len(length(gs_comp)-1)],
                                    parent = "Single Cells",
                                    compensated = FALSE,
                                    spillover_spread = "Spillover-Spread-Matrix.csv"),
      SPREAD, tolerance = 0.01)

  })

  expect_true(file.exists("Spillover-Spread-Matrix.csv"))

  # SAVED MATRIX
  spread <- read.csv("Spillover-Spread-Matrix.csv",
                     header = TRUE,
                     row.names = 1,
                     stringsAsFactors = FALSE)
  spread <- as.matrix(spread)
  colnames(spread) <- cyto_fluor_channels(gs_comp)
  rownames(spread) <- c("Alexa Fluor 488-A",
                        "PE-A",
                        "7-AAD-A",
                        "Alexa Fluor 647-A",
                        "Alexa Fluor 700-A",
                        "APC-Cy7-A")
})

# DELETE FILES
base::unlink("Spillover-Spread-Matrix.csv")
base::unlink(paste0(format(Sys.Date(), "%d%m%y"),"-Compensation-Channels.csv"))