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
  mock_data_editor <- mock(data.frame("name" = cyto_names(gs_comp[-length(gs_comp)]),
                                      "channel" = c("7-AAD-A",
                                                    "Alexa Fluor 700-A",
                                                    "APC-Cy7-A",
                                                    "Alexa Fluor 647-A",
                                                    "Alexa Fluor 488-A",
                                                    "PE-A"),
                                      stringsAsFactors = FALSE))
  # GATES
  mock_locator <- mock(list("x" = -0.109,
                            "y" = 50),
                       list("x" = 1.234,
                            "y" = 50),
                       list("x" = 1.927,
                            "y" = 50),
                       list("x" = 3.685,
                            "y" = 50),
                       list("x" = -0.204,
                            "y" = 50),
                       list("x" = 1.637,
                            "y" = 50),
                       list("x" = 2.012,
                            "y" = 50),
                       list("x" = 4.003,
                            "y" = 50),
                       list("x" = -0.126,
                            "y" = 50),
                       list("x" = 1.465,
                            "y" = 50),
                       list("x" = 1.83,
                            "y" = 50),
                       list("x" = 3.808,
                            "y" = 50),
                       list("x" = -0.838,
                            "y" = 50),
                       list("x" = 1.39,
                            "y" = 50),
                       list("x" = 2.602,
                            "y" = 50),
                       list("x" = 4.622,
                            "y" = 50),
                       list("x" = 0.306,
                            "y" = 50),
                       list("x" = 1.702,
                            "y" = 50),
                       list("x" = 2.587,
                            "y" = 50),
                       list("x" = 4.421,
                            "y" = 50),
                       list("x" = 0.435,
                            "y" = 50),
                       list("x" = 2.14,
                            "y" = 50),
                       list("x" = 2.891,
                            "y" = 50),
                       list("x" = 4.305,
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
  testthat::with_mock(data_editor = mock_data_editor,
                      locator = mock_locator,{
    expect_equal(
      cyto_spillover_spread_compute(gs_comp_trans[seq_len(length(gs_comp)-1)],
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