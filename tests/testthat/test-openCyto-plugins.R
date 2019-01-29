context("openCyto Plugins")

# .gate_manual -------------------------------------------------------------------------------

test_that("gate_manual", {
  
  expect_equal(unname(.gate_manual(fs[[1]], channels = c("FSC-A","SSC-A"), alias = "Cells", type = "r")), filters(list(rg)))
  
})

# .pp_gate_draw ------------------------------------------------------------------------------

test_that("pp_gate_draw", {
  
  expect_message(.pp_gate_draw(fs = fs, gs = gs, groupBy = 1), "Numeric groupBy is not supported, use pData variables instead. All samples will be grouped together.")
  expect_equal(.pp_gate_draw(fs = fs, gs = gs, groupBy = NA), 1)
  
})