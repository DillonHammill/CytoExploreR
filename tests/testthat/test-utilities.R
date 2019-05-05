context("utilities")

# .empty -----------------------------------------------------------------------

test_that(".empty", {
  
  expect_true(.empty(""))
  expect_true(.empty(c("","")))
  
})