context("openCyto Plugins")

# .cyto_gate_manual -----------------------------------------------------------------

test_that("cyto_gate_manual", {
  
  mock_locator <- mock(list("x" = c(0, 50000),
                            "y" = c(0, 50000)))
  testthat::with_mock(
    locator = mock_locator,
      expect_equal(.cyto_gate_manual(fs[[1]], 
                                   channels = c("FSC-A","SSC-A"), 
                                   alias = "A", 
                                   type = "r",
                                   plot = TRUE,
                                   display = 100), 
               list("A" = rg2))
  )

})

# .pp_cyto_gate_draw ----------------------------------------------------------------

test_that("pp_cyto_gate_draw", {
  
  expect_message(.pp_cyto_gate_draw(fs = fs, 
                               gs = gs,
                               groupBy = 1), 
                 "Numeric groupBy is not supported.")
  expect_equal(.pp_cyto_gate_draw(fs = fs, 
                             gs = gs, 
                             groupBy = NA), 1)
  
})