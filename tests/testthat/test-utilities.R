context("utilities")

# .empty -----------------------------------------------------------------------

test_that(".empty", {
  
  expect_true(.empty(""))
  expect_false(.empty("test"))
  expect_true(.empty(c("","")))
  
})

# .all_na ----------------------------------------------------------------------

test_that(".all_na", {
  
  expect_true(.all_na(NA))
  expect_false(.all_na(""))
  expect_false(.all_na(c(NA,"test")))
  expect_true(.all_na(list(NA,NA,NA)))
  
})

# args_list --------------------------------------------------------------------

test_that("args_list", {
  
  f <- function(x,
                y = 1,
                z = 3){
    
    .args_list()
    
  }
  ref <- alist("x" = "", "y" = 1, "z" = 3)
  expect_equal(f()[["x"]], "")
  expect_equal(f()[["y"]], 1)
  expect_equal(f()[["z"]], 3)
  
})

# file_wd_check ----------------------------------------------------------------

test_that("file_wd_check", {
  
  expect_true(file_wd_check("Activation-gatingTemplate.csv"))
  expect_false(file_wd_check("Test.csv"))
  
})

# .ROUND -----------------------------------------------------------------------

test_that(".round", {
  expect_equal(.round(34.5678), 
               "34.57")
})

# LAPPLY -----------------------------------------------------------------------

test_that("LAPPLY", {
  expect_equal(LAPPLY(c("A","B","C"), function(x){x}),
               c("A","B","C"))
})

# .PAR -------------------------------------------------------------------------

test_that(".par", {
  expect_equal(.par("mfrow"), 
               list("mfrow" = c(1,1)))
  expect_equal(.par(c("mfrow","mfcol")), 
               list("mfrow" = c(1,1),
                    "mfcol" = c(1,1)))
})
