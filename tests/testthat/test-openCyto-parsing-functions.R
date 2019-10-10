context("openCyto parsing functions")

# .argDeparser -----------------------------------------------------------------

test_that(".argDeparser", {
  
  expect_equal(.argDeparser(
    list(gate = filters(list(rg2)))), 
    paste0("gate = new('filters', .Data = list(new('rectangleGate'",
    ", min = c('FSC-A' = 0, 'SSC-A' = 0), max = c('FSC-A' = 50000",
    ", 'SSC-A' = 50000), parameters = new('parameters',     .Data = list",
    "(new('unitytransform', .Data = function ()     NULL, parameters = 'FSC-A'",
    ", transformationId = 'defaultUnityTransform'),         new('unitytran",
    "sform', .Data = function ()         NULL, parameters = 'SSC-A', transf",
    "ormationId = 'defaultUnityTransform'))),     filterId = 'Cells')))"))
  
  expect_equal(.argDeparser(
    list(gate = filters(list(rg2))), 
    split = FALSE), 
    paste("new(\"rectangleGate\", min = c(`FSC-A` = 0, `SSC-A` = 0),",
          "max = c(`FSC-A` = 50000, `SSC-A` = 50000),",
          "parameters = new(\"parameters\",",
          ".Data = list(new(\"unitytransform\", .Data = function () \nNULL,",
          "parameters = \"FSC-A\",",
          "transformationId = \"defaultUnityTransform\"),", 
          "new(\"unitytransform\", .Data = function () \nNULL,",
          "parameters = \"SSC-A\",",
          "transformationId = \"defaultUnityTransform\"))),",
          "filterId = \"Cells\")"))
  
})