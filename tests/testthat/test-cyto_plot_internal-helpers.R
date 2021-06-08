# .CYTO_PLOT_DATA --------------------------------------------------------------

test_that(".cyto_plot_data", {
  
  cyto_data <- .cyto_plot_data(gs,
                               parent = "T Cells",
                               overlay = c("CD4 T Cells"),
                               merge_by = "Treatment",
                               select = 1:32,
                               display = 500)
  
  groups <- cyto_groups(gs[1:32], 
                        group_by = "Treatment")
  
  expect_equal(names(cyto_data), groups)
  expect_equal(names(cyto_data[[1]]), c("T Cells", "CD4 T Cells"))
  
})

# .CYTO_PLOT_GATES -------------------------------------------------------------

test_that(".cyto_plot_gates", {
  
  gates <- .cyto_plot_gates(gs,
                            parent = "T Cells",
                            alias = "",
                            channels = c("Alexa Fluor 700-A",
                                         "Alexa Fluor 488-A"),
                            merge_by = "Treatment",
                            select = 1:32,
                            negate = TRUE)

  groups <- cyto_groups(gs[1:32], 
                        group_by = "Treatment")
  
  expect_equal(names(gates), groups)
  expect_setequal(names(gates[[1]]), c("CD8 T Cells",
                                       "CD4 T Cells",
                                       "negate"))
  expect_setequal(sapply(gates[[1]], cyto_class), c("rectangleGate",
                                                    "rectangleGate",
                                                    "complementFilter"))
  expect_setequal(sapply(gates[[1]], function(z){
    z@filterId
  }), c("CD8 T Cells", 
        "CD4 T Cells", 
        "not CD4 T Cells or CD8 T Cells"))
  
})