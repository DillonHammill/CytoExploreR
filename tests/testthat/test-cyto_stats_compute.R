# CYTO_STATS_COMPUTE -----------------------------------------------------------

test_that("cyto_stats_compute", {
  
  # ERROR
  expect_error(cyto_stats_compute(gs_sub,
                                  channels = c("CD4", "CD8")),
               "Supply the name of the statistical function to 'stat'.") 
  
  # MEDIAN
  expect_snapshot_output(cyto_stats_compute(gs_sub,
                                            alias = c("CD4 T Cells",
                                                      "CD8 T Cells"),
                                            stat = "median",
                                            select = 29:33))
  
  # FREQUENCY
  expect_snapshot_output(cyto_stats_compute(gs_sub,
                                            alias = c("CD4 T Cells",
                                                      "CD8 T Cells"),
                                            parent = c("Live Cells",
                                                       "T Cells"),
                                            stat = "freq",
                                            select = 29:33,
                                            format = "long"))
  
  # GEOMETRIC MEAN
  expect_snapshot_output(cyto_stats_compute(gs_sub,
                                            alias = "CD8 T Cells",
                                            stat = "geomean",
                                            select = 29:33,
                                            channels = c("FSC-A",
                                                         "PE-A"),
                                            details = FALSE,
                                            markers = FALSE))
  
  # MODE
  expect_snapshot_output(cyto_stats_compute(gs_sub,
                                            alias = "CD4 T Cells",
                                            stat = "mode",
                                            select = 29:33,
                                            channels = c("FSC-A",
                                                         "PE-A"),
                                            tibble = TRUE))
  
  # # AUC
  # expect_snapshot_output(cyto_stats_compute(gs_sub,
  #                                           alias = "CD4 T Cells",
  #                                           stat = "AUC",
  #                                           select = 29:33,
  #                                           channels = c("CD4",
  #                                                        "PE-A"),
  #                                           details = FALSE,
  #                                           markers = FALSE))
  
  # CUSTOM FUNCTION
  mean_and_count <- function(x) {
    c("mean" = sum(x)/length(x),
      "N" = length(x))
  }
  expect_snapshot_output(cyto_stats_compute(gs_sub,
                                            alias = "CD4 T Cells",
                                            stat = mean_and_count,
                                            select = 29:33,
                                            channels = c("CD4",
                                                         "PE-A"),
                                            details = FALSE,
                                            markers = FALSE,
                                            format = "long",
                                            input = "column",
                                            save_as = paste0(temp_dir, 
                                                             "stats.csv")))
  expect_true(file.exists(paste0(temp_dir, "stats.csv")))
  unlink(paste0(temp_dir, "stats.csv"))
  
})