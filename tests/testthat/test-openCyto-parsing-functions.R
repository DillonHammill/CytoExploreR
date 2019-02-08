context("openCyto parsing functions")

# .argDeparser -----------------------------------------------------------------

test_that(".argDeparser", {
  
  expect_equal(.argDeparser(
    list(gate = filters(list(rg)))), 
    paste("gate = new('filters', .Data = list(new('rectangleGate',",
          "min = c('FSC-A' = 25000, 'SSC-A' = 5000),",
          "max = c('FSC-A' = 150000,", 
          "'SSC-A' = 150000),     parameters = new('parameters',",
          ".Data = list(new('unitytransform',        ",
          ".Data = function ()         NULL, parameters = 'FSC-A',",
          "transformationId = 'defaultUnityTransform'),        ",
          "new('unitytransform', .Data = function ()         NULL,",
          "parameters = 'SSC-A', transformationId =",
          "'defaultUnityTransform'))),     filterId = 'Cells')))"))
  
  expect_equal(.argDeparser(
    list(gate = filters(list(rg))), 
    split = FALSE), 
    paste("new(\"rectangleGate\", min = c(`FSC-A` = 25000, `SSC-A` = 5000),",
          "max = c(`FSC-A` = 150000, `SSC-A` = 150000),",
          "parameters = new(\"parameters\",",
          ".Data = list(new(\"unitytransform\", .Data = function () \nNULL,",
          "parameters = \"FSC-A\",",
          "transformationId = \"defaultUnityTransform\"),", 
          "new(\"unitytransform\", .Data = function () \nNULL,",
          "parameters = \"SSC-A\",",
          "transformationId = \"defaultUnityTransform\"))),",
          "filterId = \"Cells\")"))
  
})