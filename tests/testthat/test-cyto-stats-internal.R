context("cyto-stats")

# .cyto_count ------------------------------------------------------------------

test_that(".cyto_count", {
  
  expect_equal(.cyto_count(Va2[[4]])[,"Count"],
               234)
  
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
  
  exp <- .cyto_mean(Va2[[4]], c("FSC-A","Va2","CD4"))
  expect_s3_class(exp, "data.frame")
  expect_equal(rownames(exp), "Activation4.fcs")
  expect_equal(colnames(exp), c("FSC-A","Va2","CD4"))
  expect_equal(exp[,"FSC-A"], 70187.58, tolerance = 0.001)
  expect_equal(exp[,"Va2"], 3.351131, tolerance = 0.001)
  expect_equal(exp[,"CD4"], 1.388182, tolerance = 0.001)
  
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
  
  exp <- .cyto_geometric_mean(Va2[[4]], 
                             c("FSC-A","Va2","CD4", "Hoechst-430"), 
                             trans)
  
  expect_s3_class(exp, "data.frame")
  expect_equal(rownames(exp), "Activation4.fcs")
  expect_equal(colnames(exp), c("FSC-A","Va2","CD4", "Hoechst-430"))
  expect_equal(exp[,"FSC-A"], 67578.13, tolerance = 0.001)
  expect_equal(exp[,"Va2"], 18969.77, tolerance = 0.001)
  expect_equal(exp[,"CD4"], 334.6669, tolerance = 0.001)
  expect_equal(exp[,"Hoechst-430"], 76.23341, tolerance = 0.001)
  
})

# .cyto_median -----------------------------------------------------------------

test_that(".cyto_median", {
  
  # Message - trans object required for transformed channels
  expect_message(.cyto_median(Va2[[4]], c("FSC-A","Va2","CD4")),
            paste(
              "'trans' requires a transformList/transformerList to calculate",
              "statistics on a linear scale for transformed channels."
            ))  
  
  # Error flowFrame onjects only
  expect_error(.cyto_median(Va2, c("FSC-A","Va2","CD4")),
               "'x' should be a flowFrame object.")
  
  exp <- .cyto_median(Va2[[4]], c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, "data.frame")
  expect_equal(rownames(exp), "Activation4.fcs")
  expect_equal(colnames(exp), c("FSC-A","Va2","CD4"))
  expect_equal(exp[,"FSC-A"], 64473.5, tolerance = 0.001)
  expect_equal(exp[,"Va2"], 20394.91, tolerance = 0.001)
  expect_equal(exp[,"CD4"], 89.22609, tolerance = 0.001)
  
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
  
  exp <- .cyto_mode(Va2[[4]], c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, "data.frame")
  expect_equal(rownames(exp), "Activation4.fcs")
  expect_equal(colnames(exp), c("FSC-A","Va2","CD4"))
  expect_equal(exp[,"FSC-A"], 62791.87, tolerance = 0.001)
  expect_equal(exp[,"Va2"], 19952.18, tolerance = 0.001)
  expect_equal(exp[,"CD4"], 14.76467, tolerance = 0.001)
  
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
  
  exp <- .cyto_CV(Va2[[4]], c("FSC-A","Va2","CD4"), trans)
  expect_s3_class(exp, "data.frame")
  expect_equal(rownames(exp), "Activation4.fcs")
  expect_equal(colnames(exp), c("FSC-A","Va2","CD4"))
  expect_equal(exp[,"FSC-A"], 14.38253, tolerance = 0.001)
  expect_equal(exp[,"Va2"], 34.00737, tolerance = 0.001)
  expect_equal(exp[,"CD4"], 418.0559, tolerance = 0.001)
  
})

# .cyto_density ----------------------------------------------------------------

test_that(".cyto_density", {
  
  # Error - flowFrame objects only
  expect_error(.cyto_density(Va2, channel = "CD4"),
               "'x' should be a flowFrame object.")
  
  exp <- .cyto_density(Va2[[4]], channel = "CD4")
  expect_s3_class(exp, "density")
  
})
