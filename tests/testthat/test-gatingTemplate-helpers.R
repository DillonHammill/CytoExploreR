context("gatingTemplate-helpers")

# gatingTemplate_select --------------------------------------------------------

test_that("gatingTemplate_select", {
  
  gatingTemplate_select("Activation")
  expect_equal(getOption("CytoRSuite_gatingTemplate"),
               "Activation.csv")
  
  gatingTemplate_select("Activation_gatingTemplate.csv")
  expect_equal(getOption("CytoRSuite_gatingTemplate"),
               "Activation_gatingTemplate.csv")
  
})

# gatingTemplate_edit ----------------------------------------------------------

# gatingTemplate_convert -------------------------------------------------------

test_that("gatingTemplate_convert", {
  
  old.gt <- read.csv("Ref-old-gatingTemplate.csv", header = TRUE)
  
  gatingTemplate_convert(gs, "Ref-old-gatingTemplate.csv")
  
  gt_convert_ref <- read.csv("Converted-gatingTemplate.csv", header = TRUE)
  gt_convert <- read.csv("Ref-old-gatingTemplate.csv", header = TRUE)
  
  expect_equal(gt_convert, gt_convert_ref)
  
  write.csv(old.gt, "Ref-old-gatingTemplate.csv", row.names = FALSE)
})

write.csv(gtf, "Activation-gatingTemplate.csv", row.names = FALSE)
