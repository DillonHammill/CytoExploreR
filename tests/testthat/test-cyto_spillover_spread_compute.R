# context("cyto_spillover_spread_compute")
# 
# # CYTO_SPILLOVER_SPREAD_COMPUTE ------------------------------------------------
# 
# # UNIVERSAL UNSTAINED CONTROL
# test_that("cyto_spillover_spread_compute universal reference", {
#   
#   # GATES
#   mock_locator <- mock(list("x" = c(1740, 4265),
#                             "y" = c(50,50)),
#                        list("x" = c(1800, 3100),
#                             "y" = c(50,50)),
#                        list("x" = c(1710, 2950),
#                             "y" = c(50,50)),
#                        list("x" = c(2720, 4115),
#                             "y" = c(50,50)),
#                        list("x" = c(2310, 3780),
#                             "y" = c(50,50)),
#                        list("x" = c(2350, 4120),
#                             "y" = c(50, 50)),
#                        cycle = TRUE)
#   
#   # SPILLOVER SPREAD
#   SPREAD <- read.csv("Reference-Universal-Spillover-Spread-Matrix.csv",
#                      header = TRUE,
#                      row.names = 1,
#                      stringsAsFactors = FALSE)
#   SPREAD <- as.matrix(SPREAD)
#   colnames(SPREAD) <- cyto_fluor_channels(gs_comp)
#   rownames(SPREAD) <- c("Alexa Fluor 488-A",
#                         "PE-A",
#                         "7-AAD-A",
#                         "Alexa Fluor 647-A",
#                         "Alexa Fluor 700-A",
#                         "APC-Cy7-A")
# 
#   # GATINGSET
#   testthat::with_mock(locator = mock_locator,{
#     expect_equal(
#       cyto_spillover_spread_compute(gs_comp,
#                                     parent = "Single Cells",
#                                     channel_match = "Reference-Compensation-Channels.csv",
#                                     compensated = FALSE,
#                                     spillover = "Reference-Universal-Spillover-Matrix.csv",
#                                     spillover_spread = "Spillover-Spread-Matrix.csv"),
#       SPREAD, tolerance = 0.01)
#   })
#   
#   expect_true(file.exists("Spillover-Spread-Matrix.csv"))
#   
#   # SAVED MATRIX
#   spread <- read.csv("Spillover-Spread-Matrix.csv",
#                     header = TRUE,
#                     row.names = 1,
#                     stringsAsFactors = FALSE)
#   spread <- as.matrix(spread)  
#   colnames(spread) <- cyto_fluor_channels(gs_comp)
#   rownames(spread) <- c("Alexa Fluor 488-A",
#                         "PE-A",
#                         "7-AAD-A",
#                         "Alexa Fluor 647-A",
#                         "Alexa Fluor 700-A",
#                         "APC-Cy7-A")
# 
# })
# 
# # INTERNAL UNSTAINED CONTROL
# test_that("cyto_spillover_spread_compute internal reference", {
# 
#   # CHANNELS
#   mock_menu <- mock(4, 10, 11, 9, 1, 2, cycle = TRUE)
# 
#   # GATES
#   mock_locator <- mock(list("x" = c(-850, 525),
#                             "y" = c(50,50)),
#                        list("x" = c(1770, 3950),
#                             "y" = c(50,50)),
#                        list("x" = c(125, 1080),
#                             "y" = c(50,50)),
#                        list("x" = c(1930, 2920),
#                             "y" = c(50,50)),
#                        list("x" = c(-30, 1030),
#                             "y" = c(50,50)),
#                        list("x" = c(1710, 2900),
#                             "y" = c(50, 50)),
#                        list("x" = c(185, 1700),
#                             "y" = c(50,50)),
#                        list("x" = c(2800, 4180),
#                             "y" = c(50,50)),
#                        list("x" = c(-235, 1325),
#                             "y" = c(50,50)),
#                        list("x" = c(2400, 3700),
#                             "y" = c(50,50)),
#                        list("x" = c(-125, 1495),
#                             "y" = c(50,50)),
#                        list("x" = c(2600, 3900),
#                             "y" = c(50, 50)),
#                        cycle = TRUE)
# 
#   # SPILLOVER SPREAD
#   SPREAD <- read.csv("Reference-Internal-Spillover-Spread-Matrix.csv",
#                      header = TRUE,
#                      row.names = 1,
#                      stringsAsFactors = FALSE)
#   SPREAD <- as.matrix(SPREAD)
#   colnames(SPREAD) <- cyto_fluor_channels(gs_comp)
#   rownames(SPREAD) <- c("Alexa Fluor 488-A",
#                         "PE-A",
#                         "7-AAD-A",
#                         "Alexa Fluor 647-A",
#                         "Alexa Fluor 700-A",
#                         "APC-Cy7-A")
# 
#   # GATINGSET
#   testthat::with_mock(locator = mock_locator,{
#     expect_equal(
#       cyto_spillover_spread_compute(gs_comp[seq_len(length(gs_comp)-1)],
#                                     parent = "Single Cells",
#                                     channel_match = "Reference-Compensation-Channels.csv",
#                                     compensated = FALSE,
#                                     spillover_spread = "Spillover-Spread-Matrix.csv"),
#       SPREAD, tolerance = 0.01)
# 
#   })
# 
#   expect_true(file.exists("Spillover-Spread-Matrix.csv"))
# 
#   # SAVED MATRIX
#   spread <- read.csv("Spillover-Spread-Matrix.csv",
#                      header = TRUE,
#                      row.names = 1,
#                      stringsAsFactors = FALSE)
#   spread <- as.matrix(spread)
#   colnames(spread) <- cyto_fluor_channels(gs_comp)
#   rownames(spread) <- c("Alexa Fluor 488-A",
#                         "PE-A",
#                         "7-AAD-A",
#                         "Alexa Fluor 647-A",
#                         "Alexa Fluor 700-A",
#                         "APC-Cy7-A")
# 
# })
# 
# # DELETE FILES
# base::unlink("Spillover-Spread-Matrix.csv")
