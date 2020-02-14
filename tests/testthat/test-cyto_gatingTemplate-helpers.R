context("cyto_gatingTemplate-helpers")

# MISSING TESTS
# CYTO_GATINGTEMPLATE_EDIT
# CYTO_GATINGTEMPLATE_APPLY

# cyto_gatingTemplate_select/active --------------------------------------------

test_that("cyto_gatingTemplate_select", {
  
  cyto_gatingTemplate_select("Activation")
  expect_equal(getOption("CytoExploreR_gatingTemplate"),
               "Activation.csv")
  
  cyto_gatingTemplate_select("Activation_gatingTemplate.csv")
  expect_equal(getOption("CytoExploreR_gatingTemplate"),
               "Activation_gatingTemplate.csv")
  
  expect_equal(cyto_gatingTemplate_active(),
               "Activation_gatingTemplate.csv")
  
})

# cyto_gatingTemplate_create ---------------------------------------------------

test_that("cyto_gatingTemplate_create", {
  
  cyto_gatingTemplate_create("New-gatingTemplate.csv")
  new_gt <- read.csv("New-gatingTemplate.csv", header = TRUE)
  expect_equal(colnames(new_gt), c("alias",
                                   "pop",
                                   "parent",
                                   "dims",
                                   "gating_method",
                                   "gating_args",
                                   "collapseDataForGating",
                                   "groupBy",
                                   "preprocessing_method",
                                   "preprocessing_args"))
  
})

# cyto_gatingTemplate_update --------------------------------------------------

test_that("cyto_gatingTemplate_update", {
  
  # UPDATE CYTORSUITE GATINGTEMPLATE
  cyto_gatingTemplate_update("Ref-old-gatingTemplate.csv")
  expect_equal(read.csv("Updated-gatingTemplate.csv", header =TRUE),
               read.csv("Ref-updated-gatingTemplate.csv", header = TRUE))
  
})

# .cyto_gatingTemplate_check ---------------------------------------------------

test_that(".cyto_gatingTemplate_check", {
  
  # BYPASS CHECKS
  expect_equal(.cyto_gatingTemplate_check(parent = "T Cells",
                              alias = "CD4 T Cells",
                              "NOTINWD.csv"),
               NULL)
  
  # ENTRY EXISTS ALREADY 
  expect_error(.cyto_gatingTemplate_check(parent = "T Cells",
                                          alias = "CD4 T Cells",
                                         "Reference-Activation-gatingTemplate.csv"),
      "Supply another gatingTemplate or edit gate(s) using cyto_gate_edit.",
      fixed = TRUE)
  
})

unlink("Updated-gatingTemplate.csv")
unlink("New-gatingTemplate.csv")