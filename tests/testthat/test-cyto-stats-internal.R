context("cyto-stats")

# NOTE: expect_equal does not work with tibbles so we need to coerce to
# data.frames for testing purposes.

# .cyto_count ------------------------------------------------------------------

test_that(".cyto_count", {
  
  ref <- tibble("count" = BiocGenerics::nrow(Va2[[4]]))
  expect_equal(.cyto_count(Va2[[4]])[,"count"],
               ref)
  
})

# .cyto_mean -------------------------------------------------------------------

test_that(".cyto_mean", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_mean(Va2[[4]], c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error - flowFrame objects only
  expect_error(.cyto_mean(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 70187.580779,
                "Va2" = 3.351131,
                "CD4" = 1.388182)
  exp <- .cyto_mean(Va2[[4]], c("FSC-A","Va2","CD4"))
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_geometric_mean ---------------------------------------------------------

test_that(".cyto_geometric_mean", {
  
  # Throw error if you try to get geometric mean of transformed channels
  # without supplying transformation object
  expect_error(.cyto_geometric_mean(Va2[[4]],
                                    c("FSC-A","Va2","CD4","Hoechst-430")),
           "Supply transformList/transformerList to calculate geometric mean.")
  
  # Error - flowFrame objects only
  expect_error(.cyto_geometric_mean(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 67578.12698,
                "Va2" = 18969.77070,
                "CD4" = 334.66689,
                "Hoechst-430" = 76.23341)
  exp <- .cyto_geometric_mean(Va2[[4]], 
                             c("FSC-A","Va2","CD4", "Hoechst-430"), 
                             trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_median -----------------------------------------------------------------

test_that(".cyto_median", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_median(Va2[[4]], c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error flowFrame objects only
  expect_error(.cyto_median(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 64473.5,
                "Va2" = 20394.90959,
                "CD4" = 89.22609)
  exp <- .cyto_median(Va2[[4]], c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_mode -------------------------------------------------------------------

test_that(".cyto_mode", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_mode(Va2[[4]], c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error - flowFrame objects only
  expect_error(.cyto_mode(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 62791.87068,
                "Va2" = 19952.17720,
                "CD4" = 14.76467)
  exp <- .cyto_mode(Va2[[4]], c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_CV ---------------------------------------------------------------------

test_that(".cyto_CV", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_CV(Va2[[4]], c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error - flowFrame objects only
  expect_error(.cyto_CV(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 14.382532,
                "Va2" = 34.00737,
                "CD4" = 418.05588)
  exp <- .cyto_CV(Va2[[4]], c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_density ----------------------------------------------------------------

test_that(".cyto_density", {
  
  # Error - flowFrame objects only
  expect_error(.cyto_density(Va2, channel = "CD4"),
               "'x' should be a flowFrame object.")
  
  exp <- .cyto_density(Va2[[4]], channel = "CD4")
  expect_s3_class(exp, "density")
  
})
