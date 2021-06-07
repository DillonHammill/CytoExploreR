# .CYTO_PLOT_DATA --------------------------------------------------------------

test_that(".cyto_plot_data", {
  
  cyto_data <- .cyto_plot_data(gs,
                               parent = "T Cells",
                               overlay = c("CD4 T Cells"),
                               merge_by = "Treatment",
                               select = 1:32,
                               display = 500)
  
  groups <- cyto_groups(gs[1:32], group_by = "Treatment")
  
  expect_equal(names(cyto_data), groups)
  expect_equal(names(cyto_data[[1]]), c("T Cells", "CD4 T Cells"))
  
})