context("cyto-stats")

# NOTE: expect_equal does not work with tibbles so we need to coerce to
# data.frames for testing purposes.

# .cyto_count ------------------------------------------------------------------

test_that(".cyto_count", {
  
  ref <- tibble("count" = BiocGenerics::nrow(fr_test))
  exp <- .cyto_count(fr_test)[,"count"]
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref))
  
})

# .cyto_mean -------------------------------------------------------------------

test_that(".cyto_mean", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_mean(fr_test, c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error - flowFrame objects only
  expect_error(.cyto_mean(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 70372.9868,
                "Va2" = 18939.2329,
                "CD4" = 982.8482)
  exp <- .cyto_mean(fr_test, 
                    c("FSC-A","Va2","CD4"), 
                    trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_geometric_mean ---------------------------------------------------------

test_that(".cyto_geometric_mean", {
  
  # Error - flowFrame objects only
  expect_error(.cyto_geometric_mean(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 67637.01451,
                "Va2" = 17170.08487,
                "CD4" = 280.31106,
                "Hoechst-430" = 82.27335)
  exp <- .cyto_geometric_mean(fr_test, 
                             c("FSC-A","Va2","CD4", "Hoechst-430"), 
                             trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_median -----------------------------------------------------------------

test_that(".cyto_median", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_median(fr_test, c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error flowFrame objects only
  expect_error(.cyto_median(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 64774.15039,
                "Va2" = 18448.62970,
                "CD4" = 63.38221)
  exp <- .cyto_median(fr_test, c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_mode -------------------------------------------------------------------

test_that(".cyto_mode", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_mode(fr_test, c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error - flowFrame objects only
  expect_error(.cyto_mode(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 62786.47125,
                "Va2" = 17915.35399,
                "CD4" = 18.10466)
  exp <- .cyto_mode(fr_test, c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_CV ---------------------------------------------------------------------

test_that(".cyto_CV", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_CV(fr_test, c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error - flowFrame objects only
  expect_error(.cyto_CV(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  ref <- tibble("FSC-A" = 16.64379,
                "Va2" = 39.77155,
                "CD4" = 432.67621)
  exp <- .cyto_CV(fr_test, c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, c("tbl_df","tbl","data.frame"))
  expect_equal(as.data.frame(exp), as.data.frame(ref), tolerance = 0.001)
  
})

# .cyto_density ----------------------------------------------------------------

test_that(".cyto_density", {
  
  # Error - flowFrame objects only
  expect_error(.cyto_density(Va2, channel = "CD4"),
               "'x' should be a flowFrame object.")
  
  exp <- .cyto_density(fr_test, channel = "CD4")
  expect_s3_class(exp, "density")
  
})
