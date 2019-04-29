context("Validation Internal Functions")

## .cyto_gate_type_check -------------------------------------------------------

test_that(".cyto_gate_type_check", {
  expect_equal(.cyto_gate_type_check(type = "r", 
                                     alias = c("A", "B", "C")), 
               c("rectangle", "rectangle", "rectangle"))
  expect_equal(.cyto_gate_type_check(type = "p", 
                                     alias = c("A", "B", "C")), 
               c("polygon", "polygon", "polygon"))
  expect_equal(.cyto_gate_type_check(type = "q", 
                                     alias = c("A", "B", "C", "D")), 
               "quadrant")
  expect_error(.cyto_gate_type_check(type = "q", 
                                     alias = c("A", "B", "C")), 
               "Supply the names of 4 poulations to alias for quadrant gates.", 
               fixed = TRUE)
  expect_error(.cyto_gate_type_check(type = c("v", "j"), 
                                     alias = c("A", "B", "C")), 
               "v & j are not valid gate types for gate_draw!", 
               fixed = TRUE)
  expect_error(.cyto_gate_type_check(type = "z",
                                     alias = c("A", "B", "C")), 
               "z is not a valid gate type for gate_draw!", 
               fixed = TRUE)
})

## .cyto_alias_check -----------------------------------------------------------

test_that(".cyto_alias_check stops gating process if alias is incorrect", {
  expect_error(.cyto_alias_check(type = "web"), 
               "Supply names of populations to 'alias' for checking.", 
               fixed = TRUE)
  expect_error(.cyto_alias_check(type = "quadrant", 
                                 alias = "Cells"), 
          "Supply 4 population names to 'alias' to construct quadrant gates.", 
               fixed = TRUE)
  expect_error(.cyto_alias_check(type = c("rectangle", "polygon"),
                                 alias = "Cells"), 
          "Length of alias must be the same length as type for multi-gates.",
               fixed = TRUE)
})

## .cyto_gatingTemplate_check --------------------------------------------------

test_that(".cyto_gatingTemplate_check", {
  
  expect_error(
    .cyto_gatingTemplate_check(parent = "Cells", 
                               alias = "Single Cells", 
                               gatingTemplate = "Compensation-gatingTemplate.csv"), 
    "Supply another gatingTemplate or edit gate(s) using gate_edit.",
    fixed = TRUE)
  
  expect_error(
    .cyto_gatingTemplate_check(parent = "Cells",
                               alias = "Single Cells",
                               gatingTemplate = gt),
    "'gatingTemplate' should be the name of the gatingTemplate csv file."
  )
  
})

## .cyto_stat_check ------------------------------------------------------------

test_that(".cyto_stat_check", {
  
  expect_error(.cyto_stat_check("average"), 
               "Supplied statistic not supported.")
  expect_equal(.cyto_stat_check("Mean"), "mean")
  expect_equal(.cyto_stat_check("Median"), "median")
  expect_equal(.cyto_stat_check("Mode"), "mode")
  expect_equal(.cyto_stat_check("Count"), "count")
  expect_equal(.cyto_stat_check("Freq"), "freq")
  expect_equal(.cyto_stat_check("Geo mean"), "geo mean")
  expect_equal(.cyto_stat_check("cv"), "CV")
  
})