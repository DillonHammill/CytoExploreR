context("spillover_compute")

# .getCompleteTransList --------------------------------------------------------

test_that(".getCompleteTransList", {
  
  # transList of incorrect class
  expect_error(.getCompleteTransList(fs[[4]], "Test"), 
               "'trans' should be a transformList or transformerList object.")
  
  # flowFrame/flowSet transformed but no transList
  expect_error(.getCompleteTransList(getData(gs, "T Cells")[[1]]), 
               "Looks like the data is already transformed. 
 Please supply the transformList/transformerList used.")
  expect_error(.getCompleteTransList(getData(gs, "T Cells")), 
               "Looks like the data is already transformed. 
 Please supply the transformList/transformerList used.")
  
  fst <- transform(fs, cyto_trans_check(trans, inverse = FALSE))
  gst <- GatingSet(fst)
  
  expect_error(.getCompleteTransList(gst),
               "Looks like the data is already transformed. 
 Please supply the transformList/transformerList used.")
  
  # NULL transList -------------------------------------------------------------
  
  # flowFrame no transformations - no transList
  tr <- .getCompleteTransList(fs[[4]])
  
  expect_s4_class(tr, "transformList")
  expect_equal(names(tr@transforms), cyto_fluor_channels(fs[[4]]))
  
  # flowSet no transformations -  no transList
  tr <- .getCompleteTransList(fs)
  
  expect_s4_class(tr, "transformList")
  expect_equal(names(tr@transforms), cyto_fluor_channels(fs[[4]]))
  
  # GatingSet no transformations
  gst <- GatingSet(fs)
  tr <- .getCompleteTransList(gst)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # GatingSet transformed
  tr <- .getCompleteTransList(gs)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gs))
  
  # GatingSet with some transformations
  trns <- estimateLogicle(gst[[4]], c("PE-A", "Alexa Fluor 700-A"))
  gst <- transform(gst, trns)
  
  tr <- .getCompleteTransList(gst)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # Complete transList ---------------------------------------------------------
  
  # flowFrame complete transList
  tr <- .getCompleteTransList(fs[[4]], trans)
  
  expect_s4_class(tr, "transformList")
  expect_equal(names(tr@transforms), cyto_fluor_channels(fs[[4]]))
  
  # flowSet complete transList
  tr <- .getCompleteTransList(fs[[4]], trans)
  
  expect_s4_class(tr, "transformList")
  expect_equal(names(tr@transforms), cyto_fluor_channels(fs[[4]]))
  
  # GatingSet complete transList
  gst <- GatingSet(fs)
  tr <- .getCompleteTransList(gst, trans)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # GatingSet complete transformList
  tr <- .getCompleteTransList(gst, cyto_trans_check(trns, inverse = FALSE))
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # Transformed GatingSet - complete transList
  tr <- .getCompleteTransList(gs, trans)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gs))
  
  # Incomplete transList -------------------------------------------------------
  
  # flowFrame with some transformations
  trns <- estimateLogicle(fs[[4]], c("PE-A", "Alexa Fluor 488-A"))
  fst <- transform(fs, trns)
  
  tr <- .getCompleteTransList(fst[[4]], trns)
  
  expect_s4_class(tr, "transformList")
  expect_setequal(names(tr@transforms), cyto_fluor_channels(fs[[4]]))
  
  # flowSet with some transformations
  trns <- estimateLogicle(fs[[4]], c("PE-A", "Alexa Fluor 488-A"))
  fst <- transform(fs, trns)
  
  tr <- .getCompleteTransList(fst, trns)
  
  expect_s4_class(tr, "transformList")
  expect_setequal(names(tr@transforms), cyto_fluor_channels(fs[[4]]))
  
  # GatingSet no transforms - incomplete transList
  gst <- GatingSet(fs)
  trns <- estimateLogicle(gst[[4]], c("PE-A", "Alexa Fluor 488-A"))
  
  tr <- .getCompleteTransList(gst, trans)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # GatingSet some transforms - incomplete transList
  gst <- GatingSet(fs)
  trns <- estimateLogicle(gst[[4]], c("PE-A", "Alexa Fluor 488-A"))
  gst <- transform(gst, trns)
  
  tr <- .getCompleteTransList(gst, trns)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # GatingSet some transforms - incomplete transList
  gst <- GatingSet(fs)
  trns <- estimateLogicle(gst[[4]], c("PE-A", "Alexa Fluor 488-A"))
  gst <- transform(gst, trns)
  
  trns <- estimateLogicle(gst[[4]], c("PE-A", 
                                      "Alexa Fluor 488-A", 
                                      "Alexa Fluor 700-A"))
  
  tr <- .getCompleteTransList(gst, trns)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
  # GatingSet some transforms - complete transList
  tr <- .getCompleteTransList(gst, trans)
  
  expect_is(tr, "transformerList")
  expect_setequal(names(tr), cyto_fluor_channels(gst))
  
})

# .getTransformedData ----------------------------------------------------------

test_that(".getTransformedData", {
  expect_error(.getTransformedData(gs[[1]]), 
               "'x' must be either a flowFrame, flowSet or GatingSet.")
  
  # flowFrame raw
  fr <- .getTransformedData(fs[[1]])
  fst <- transform(fs, cyto_trans_check(trans, inverse = FALSE))
  
  expect_equal(pData(parameters(fr))[, "maxRange"], 
               pData(parameters(fst[[1]]))[, "maxRange"])
  
  # flowSet raw
  fst <- .getTransformedData(fs)
  fst2 <- transform(fs, cyto_trans_check(trans, inverse = FALSE))
  
  expect_equal(pData(parameters(fst[[1]]))[, "maxRange"], 
               pData(parameters(fst2[[1]]))[, "maxRange"])
  
  # GatingSet raw
  gst <- .getTransformedData(GatingSet(fs, trans))
  
  expect_equal(pData(parameters(getData(gst, "root")[[1]]))[, "maxRange"], 
               pData(parameters(getData(gs, "root")[[1]]))[, "maxRange"])
  
  # flowSet transformed
  fst <- transform(fs, cyto_trans_check(trans, inverse = FALSE))
  fst2 <- .getTransformedData(fst, trans)
  
  expect_equal(pData(parameters(fst2[[1]]))[, "maxRange"], 
               pData(parameters(fst[[1]]))[, "maxRange"])
  
  # flowSet some transformations
  trns <- estimateLogicle(fs[[4]], c("PE-A", "Alexa Fluor 700-A"))
  fst <- transform(fs, trns)
  fst2 <- .getTransformedData(fst, trns)
  fst3 <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
  
  expect_equal(pData(parameters(fst2[[1]]))[, "maxRange"], 
               pData(parameters(fst3[[1]]))[, "maxRange"])
  
  # GatingSet some transformations
  gst <- GatingSet(fs)
  trns <- estimateLogicle(gst[[4]], c("PE-A", "Alexa Fluor 700-A"))
  gst <- transform(gst, trns)
  gst2 <- .getTransformedData(gst, trns)
  gst3 <- transform(gst, estimateLogicle(gst[[4]], cyto_fluor_channels(fs)))
  
  expect_equal(pData(parameters(getData(gst2, "root")[[1]]))[, "maxRange"], 
               pData(parameters(getData(gst3, "root")[[1]]))[, "maxRange"])
})

# .getRawData ------------------------------------------------------------------

test_that(".getRawData", {
  
  # Transformed data without transList
  expect_error(.getRawData(getData(gs, "T Cells")), 
               "Supply a transform object to inverse transformations.")
  
  # Object of incorrect class
  expect_error(.getRawData(gs[[1]]), 
               "'x' must be either a flowFrame, flowSet or GatingSet.")
  
  # Untransformed --------------------------------------------------------------
  
  # flowFrame
  expect_equal(.getRawData(fs[[1]]), fs[[1]])
  
  # flowSet
  expect_equal(.getRawData(fs), fs)
  
  # GatingSet
  gst <- GatingSet(fs)
  expect_equal(pData(parameters(.getRawData(gst)[[1]])), 
               pData(parameters(fs[[1]])))
  
  # Transformed ----------------------------------------------------------------
  
  # flowFrame
  fr <- getData(gs, "T Cells")[[1]]
  inv <- cyto_trans_check(trans, inverse = TRUE)
  fr <- transform(fr, inv)
  
  expect_equal(.getRawData(getData(gs, "T Cells")[[1]], 
                           trans), 
               fr)
  expect_equal(.getRawData(getData(gs, "T Cells")[[1]], 
                           cyto_trans_check(trans)),
               fr)
  
  # flowSet
  fst <- getData(gs, "T Cells")
  inv <- cyto_trans_check(trans, inverse = TRUE)
  fst <- transform(fst, inv)
  
  expect_equal(.getRawData(getData(gs, "T Cells"), 
                           trans), 
               fst)
  expect_equal(.getRawData(getData(gs, "T Cells"), 
                           cyto_trans_check(trans)), 
               fst)
  
  # GatingSet
  fst <- getData(gs, "T Cells")
  inv <- cyto_trans_check(trans, inverse = TRUE)
  fst <- transform(fst, inv)
  
  expect_equal(pData(parameters(.getRawData(gs)[[1]])), 
               pData(parameters(fst[[1]])))
  
  # Some transformations - complete transList ----------------------------------
  fr <- fs[[1]]
  trns <- estimateLogicle(fs[[1]], c("PE-A", "Alexa Fluor 700-A"))
  fr <- transform(fr, trns)
  
  trns <- estimateLogicle(fs[[1]], cyto_fluor_channels(fs))
  
  expect_equal(pData(parameters(.getRawData(fr, trns))), 
               pData(parameters(fs[[1]])))
})

# spillover_compute ------------------------------------------------------------

test_that("spillover_compute", {
  
  # GatingSet method -----------------------------------------------------------
  spillover_compute(gsc, parent = "Single Cells")
  
  expect_true(.file_wd_check("Spillover-Matrix.csv"))
  
  sp <- read.csv("Spillover-Matrix.csv", header = TRUE, row.names = 1)
  colnames(sp) <- rownames(sp)
  sp <- round(sp, 3)
  
  expect_equal(sp, spill, tolerance = 0.001)
  
  # GatingSet - no parent - cmfile
  spillover_compute(gsc, channel_match = "Compensation-channels.csv")
  
  expect_true(.file_wd_check("Spillover-Matrix.csv"))
  
  sp <- read.csv("Spillover-Matrix.csv", header = TRUE, row.names = 1)
  colnames(sp) <- rownames(sp)
  sp <- round(sp, 3)
  
  expect_equal(sp, spill, tolerance = 0.001)
})

unlink("Spillover-Matrix.csv")