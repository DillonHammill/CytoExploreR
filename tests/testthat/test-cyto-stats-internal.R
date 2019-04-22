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
  
  ref <- tibble("FSC-A" = 63807.668,
                "Va2" = 21467.458,
                "CD4" = 1179.643)
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
  
  ref <- tibble("FSC-A" = 62652.05858,
                "Va2" = 19684.97676,
                "CD4" = 406.33749,
                "Hoechst-430" = 64.71426)
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
  
  ref <- tibble("FSC-A" = 61485.9004,
                "Va2" = 20932.5897,
                "CD4" = 144.1649)
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
  
  ref <- tibble("FSC-A" = 60567.6369,
                "Va2" = 20560.1825,
                "CD4" = 38.8483)
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
  
  ref <- tibble("FSC-A" = 12.57232,
                "Va2" = 35.77359,
                "CD4" = 262.80242)
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
