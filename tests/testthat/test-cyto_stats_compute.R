context("cyto_stats_compute")

# NOTE: expect_equal does not work with tibbles so we need to coerce to
# data.frames for testing purposes.

# cyto_stats_compute flowFrame method ------------------------------------------

test_that("cyto_stats_compute flowFrame method", {
  
  # count -
  ref <- tibble("name" = "Activation_32.fcs", 
                "count" = 250)
  exp <- cyto_stats_compute(fr_test, stat = "count")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref))
  
  # mean -
  expect_message(cyto_stats_compute(fr_test, stat = "mean"),
                 paste("'trans' missing - statistics will be returned on the",
                       "current scale."))
  
  ref <- tibble("name" = rep("Activation_32.fcs",3), 
                "Marker" = c("FSC-A","CD8","Va2"),
                "MFI" = c(70372.987,7226.068,18939.233))  
  exp <- cyto_stats_compute(fr_test, 
                            c("FSC-A","CD8","Va2"),
                            trans,
                            stat = "mean")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
  # geometric mean
  expect_error(cyto_stats_compute(fr_test, stat = "geo mean"),
           "Supply transformList/transformerList to calculate geometric mean.")
  
  ref <- tibble("name" = rep("Activation_32.fcs",3), 
                "Marker" = c("FSC-A","CD8","Va2"),
                "GMFI" = c(67637.015,1878.195,17170.085))  
  exp <- cyto_stats_compute(fr_test, 
                            c("FSC-A","CD8","Va2"),
                            trans,
                            stat = "geo mean")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
  # median
  expect_message(cyto_stats_compute(fr_test, stat = "median"),
                 paste("'trans' missing - statistics will be returned on the",
                       "current scale."))
  
  ref <- tibble("name" = rep("Activation_32.fcs",3), 
                "Marker" = c("FSC-A","CD8","Va2"),
                "MedFI" = c(64774.15,6065.696,18448.630))  
  exp <- cyto_stats_compute(fr_test, 
                            c("FSC-A","CD8","Va2"),
                            trans,
                            stat = "median")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
  # mode
  expect_message(cyto_stats_compute(fr_test, stat = "mode"),
                 paste("'trans' missing - statistics will be returned on the",
                       "current scale."))
  
  ref <- tibble("name" = rep("Activation_32.fcs",3), 
                "Marker" = c("FSC-A","CD8","Va2"),
                "ModFI" = c(62786.47,1197.93,17915.35))  
  exp <- cyto_stats_compute(fr_test, 
                            c("FSC-A","CD8","Va2"),
                            trans,
                            stat = "mode")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
  # CV
  expect_message(cyto_stats_compute(fr_test, stat = "CV"),
                 paste("'trans' missing - statistics will be returned on the",
                       "current scale."))
  
  ref <- tibble("name" = rep("Activation_32.fcs",3), 
                "Marker" = c("FSC-A","CD8","Va2"),
                "CV" = c(16.64379,135.05146,39.77155))  
  exp <- cyto_stats_compute(fr_test, 
                            c("FSC-A","CD8","Va2"),
                            trans,
                            stat = "CV")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# cyto_stats_compute flowSet method --------------------------------------------

test_that("cyto_stats_compute flowSet method", {
  
  # NOTE: For some reason returned tibble has class <I(chr)> for "name" column,
  # not sure what this means. Expect_equal with data.frame coercion does not
  # work so instaed coerce to simpler matrix format. These test do not look into
  # factorisation or rownames but instead check that the core computation and
  # format remains unchanged.
  
  # count
  ref <- tibble("name" = sampleNames(fs_test),
                "OVAConc" = c(0,0.005,0.05,0.5),
                "Treatment" = rep("Stim-D", 4),
                "count" = c(202,262,247,250))
  exp <- cyto_stats_compute(fs_test, stat = "count")
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.matrix(exp), as.matrix(ref))
  
  # mean
  ref <- tibble("name" = rep(sampleNames(fs_test),3),
                "OVAConc" = rep(c(0,0.005,0.05,0.5),3),
                "Treatment" = rep("Stim-D", 12),
                "Marker" = c(rep("FSC-A",4),
                             rep("CD8", 4),
                             rep("Va2", 4)),
                "MFI" = c(63428,63267,64370,70373,
                          5182,4541,5862,7226,
                          19145,16613,17668,18939)) %>%
    mutate(MFI = round(MFI, 0))
  exp <- cyto_stats_compute(fs_test, 
                            channels = c("FSC-A","CD8","Va2"),
                            stat = "mean",
                            trans = trans) %>%
    mutate(MFI = round(MFI, 0))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))
  
  # geometric mean
  ref <- tibble("name" = rep(sampleNames(fs_test),3),
                "OVAConc" = rep(c(0,0.005,0.05,0.5),3),
                "Treatment" = rep("Stim-D", 12),
                "Marker" = c(rep("FSC-A",4),
                             rep("CD8", 4),
                             rep("Va2", 4)),
                "GMFI" = c(62724,62489,63157,67637,
                           1409,1281,1660,1878,
                           17565,15011,16057,17170)) %>%
    mutate(GMFI = round(GMFI, 0))
  exp <- cyto_stats_compute(fs_test, 
                            channels = c("FSC-A","CD8","Va2"),
                            stat = "geo mean",
                            trans = trans) %>%
    mutate(GMFI = round(GMFI, 0))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))
  # median
  ref <- tibble("name" = rep(sampleNames(fs_test),3),
                "OVAConc" = rep(c(0,0.005,0.05,0.5),3),
                "Treatment" = rep("Stim-D", 12),
                "Marker" = c(rep("FSC-A",4),
                             rep("CD8", 4),
                             rep("Va2", 4)),
                "MedFI" = c(63642,63219,62049,64774,
                           5184,4760,6097,6066,
                           18587,15385,16821,18449)) %>%
    mutate(MedFI = round(MedFI, 0))
  exp <- cyto_stats_compute(fs_test, 
                            channels = c("FSC-A","CD8","Va2"),
                            stat = "median",
                            trans = trans) %>%
    mutate(MedFI = round(MedFI, 0))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))
  
  # mode
  ref <- tibble("name" = rep(sampleNames(fs_test),3),
                "OVAConc" = rep(c(0,0.005,0.05,0.5),3),
                "Treatment" = rep("Stim-D", 12),
                "Marker" = c(rep("FSC-A",4),
                             rep("CD8", 4),
                             rep("Va2", 4)),
                "ModFI" = c(65097,64373,60554,62786,
                            382,398,920,1198,
                            17838,13554,14768,17915)) %>%
    mutate(ModFI = round(ModFI, 0))
  exp <- cyto_stats_compute(fs_test, 
                            channels = c("FSC-A","CD8","Va2"),
                            stat = "mode",
                            trans = trans) %>%
    mutate(ModFI = round(ModFI, 0))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))
  
  # CV
  ref <- tibble("name" = rep(sampleNames(fs_test),3),
                "OVAConc" = rep(c(0,0.005,0.05,0.5),3),
                "Treatment" = rep("Stim-D", 12),
                "Marker" = c(rep("FSC-A",4),
                             rep("CD8", 4),
                             rep("Va2", 4)),
                "CV" = c(12,13,13,17,
                         140,139,116,135,
                         37,44,40,40)) %>%
    mutate(CV = round(CV, 0))
  exp <- cyto_stats_compute(fs_test, 
                            channels = c("FSC-A","CD8","Va2"),
                            stat = "CV",
                            trans = trans) %>%
    mutate(CV = round(CV, 0))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))
  
})

# cyto_stats_compute GatingSet method --------------------------------------------

test_that("cyto_stats_compute GatingSet method", {
  
  # median
  exp <- cyto_stats_compute(gs_test,
                            alias = c("CD4 T Cells","CD8 T Cells"),
                            channels = c("CD44","CD69","CD8","CD4"),
                            stat = "median",
                            format = "long",
                            save_as = "Test-Median") %>%
    mutate(MedFI = round(MedFI, 2))
  
  # Check file has been saved
  expect_true(.file_wd_check("Test-Median.csv"))
  
  # Remove file
  base::unlink("Test-Median.csv")
  
  ref <- tibble("name" = c(rep("Activation_26.fcs",8),
                           rep("Activation_28.fcs",8),
                           rep("Activation_30.fcs",8),
                           rep("Activation_32.fcs", 8)),
                "OVAConc" = c(rep(0,8),
                              rep(0.005,8),
                              rep(0.05,8),
                              rep(0.5, 8)),
                "Treatment" = rep("Stim-D",32),
                "Population" = rep(c("CD4 T Cells","CD8 T Cells"), 16),
                "Marker" = rep(c("CD44","CD44","CD69","CD69",
                                 "CD8","CD8","CD4","CD4"), 4),
                "MedFI" = c(623.27, 408.52, 655.19, 250.52, 212.19, 9107.77, 
                            2719.59, 1.95, 762.06, 248.34, 416.65, 184.85, 
                            137.46, 7360.77, 2537.33, -19.37, 798.78, 369.31, 
                            623.58, 193.02, 178.90, 8353.14, 2515.42, -8.16, 
                            866.42, 615.80, 1003.24, 239.68, 214.7937, 9619.73, 
                            2614.18, 1.90)) %>%
    mutate(MedFI = round(MedFI, 2))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))

  # freq
  exp <- cyto_stats_compute(gs_test,
                            parent = c("Live Cells", "T Cells"),
                            alias = c("CD4 T Cells","CD8 T Cells"),
                            stat = "freq",
                            format = "long") %>%
    mutate(Frequency = round(Frequency, 2))
  
  ref <- tibble("name" = c(rep("Activation_26.fcs",4),
                           rep("Activation_28.fcs",4),
                           rep("Activation_30.fcs",4),
                           rep("Activation_32.fcs", 4)),
                "OVAConc" = c(rep(0,4),
                              rep(0.005,4),
                              rep(0.05,4),
                              rep(0.5,4)),
                "Treatment" = rep("Stim-D",16),
                "Population" = rep(c("CD4 T Cells","CD8 T Cells"), 8),
                "Parent" = rep(c("Live Cells","Live Cells",
                                 "T Cells","T Cells"), 4),
                "Frequency" = c(12.63, 18.77, 36.63, 54.46, 14.81, 24.15, 
                                35.11, 57.25, 13.71, 24.53, 32.79, 58.70, 
                                11.59, 22.41, 30.4, 58.8)) %>%
    mutate(Frequency = round(Frequency, 2))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  
  # matrix conversion requires trim = TRUE not supported by as.matrix
  expect_equal(vapply(exp, 
                      format, 
                      FUN.VALUE = character(nrow(exp)), 
                      trim = TRUE), 
               vapply(ref, 
                      format, 
                      FUN.VALUE = character(nrow(ref)), 
                      trim = TRUE))
  
})
