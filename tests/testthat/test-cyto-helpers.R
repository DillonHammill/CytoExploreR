context("cyto-helpers")

# cyto_check -------------------------------------------------------------------

test_that("cyto_check",{
  
  # flowFrame
  expect_true(cyto_check(fs[[1]]))
  
  # flowSet
  expect_true(cyto_check(fs))
  
  # ncdfFlowSet
  expect_true(suppressMessages(cyto_check(ncdfFlowSet(fs))))
  
  # GatingHierarchy
  expect_true(cyto_check(gs[[1]]))
  
  # GatingSet
  expect_true(cyto_check(gs))
  
})

# cyto_extract ----------------------------------------------------------------

test_that("cyto_extract", {
  
  exp <- getData(gs, "root")
  
  # flowFrame ------------------------------------------------------------------
  expect_equal(cyto_extract(fs[[1]]), fs[[1]])
  
  # flowSet --------------------------------------------------------------------
  expect_equal(cyto_extract(fs), fs)
  
  # GatingHierachy -------------------------------------------------------------
  expect_equal(cyto_extract(gs[[1]]), exp[[1]])
  
  # GatingSet ------------------------------------------------------------------
  expect_equal(cyto_extract(gs), exp)
  
})

# cyto_convert ----------------------------------------------------------------

test_that("cyto_convert", {
  
  # Note: no coercion for list objects (e.g. list(fs) <- list(fs[[1]],fs[[2]]))
  # flowFrame ------------------------------------------------------------------
  expect_equal(cyto_convert(fs[[1]], "flowFrame"), fs[[1]])
  expect_equal(cyto_convert(fs[[1]], "flowFrame list"), list(fs[[1]]))
  expect_equal(cyto_convert(fs[[1]], "flowSet"), flowSet(fs[[1]]))
  expect_equal(cyto_convert(fs[[1]], "flowSet list"), list(flowSet(fs[[1]])))
  
  # flowSet --------------------------------------------------------------------
  fr_exp <- as(fs,"flowFrame")
  if ("Original" %in% BiocGenerics::colnames(fr_exp)) {
    fr_exp <- suppressWarnings(
      fr_exp[, -match("Original", BiocGenerics::colnames(fr_exp))]
    )
  }
  fr_list_exp <- lapply(seq_len(length(fs)), function(x){fs[[x]]})
  expect_equal(cyto_convert(fs, "flowFrame"), fr_exp)
  expect_equal(cyto_convert(fs, "flowSet"), fs)
  expect_equal(cyto_convert(fs, "flowFrame list"), list(fr_exp))
  expect_equal(cyto_convert(fs, "list of flowFrames"), fr_list_exp)
  expect_equal(cyto_convert(fs, "flowSet list"), list(fs))
  
  # GatingHierarchy ------------------------------------------------------------
  expect_equal(cyto_convert(gs[[1]],  "root", "flowFrame"),
               getData(gs,"root")[[1]])
  expect_equal(cyto_convert(gs[[1]],"root" ,"flowFrame list"),
               list(getData(gs,"root")[[1]]))
  expect_equal(cyto_convert(gs[[1]], "root", "flowSet"),
               flowSet(getData(gs,"root")[[1]]))
  expect_equal(cyto_convert(gs[[1]],"root", "flowSet list"),
               list(flowSet(getData(gs,"root")[[1]])))
  
  # GatingSet ------------------------------------------------------------------
  fr_exp <- as(getData(gs,"root"),"flowFrame")
  if ("Original" %in% BiocGenerics::colnames(fr_exp)) {
    fr_exp <- suppressWarnings(
      fr_exp[, -match("Original", BiocGenerics::colnames(fr_exp))]
    )
  }
  fr_list_exp <- lapply(seq_len(length(gs)), 
                        function(y){getData(gs,"root")[[y]]})
  expect_equal(cyto_convert(gs, "root", "flowFrame"),
               fr_exp)
  expect_equal(cyto_convert(gs, "root", "list of flowFrames"),
               fr_list_exp)
  expect_equal(cyto_convert(gs, "root", "flowFrame list"),
               list(fr_exp))
  expect_equal(cyto_convert(gs, "root", "flowSet"),
               getData(gs,"root"))
  expect_equal(cyto_convert(gs, "root", "flowSet list"),
               list(getData(gs,"root")))
  
})

# cyto_filter ------------------------------------------------------------------

test_that("cyto_filter", {
  
  expect_error(cyto_filter(list(fs), Treatment = "Stim-C"),
               "'x' should be an object of class flowSet or GatingSet.")
  
  expect_equal(cyto_filter(fs, Treatment == "Stim-C"), 
               fs[c(17,18,19,20,21,22,23,24)])
  
  expect_equal(cyto_filter(fs, Treatment == "Stim-A", OVAConc %in% c(0,0.5)),
               fs[c(1,2,7,8)])
  
  # Filtered GatingSet will have different guid slot
  expect_equivalent(cyto_filter(gs, Treatment == "Stim-C"), 
                                gs[c(17, 18, 19, 20, 21, 22, 23, 24)])
  
  expect_equivalent(cyto_filter(gs, 
                                Treatment == "Stim-A", 
                                OVAConc %in% c(0,0.5)),
                                gs[c(1, 2, 7, 8)])
  
})

# cyto_select ------------------------------------------------------------------

test_that("cyto_select", {
  
  # Must be flowSet or GatingSet object
  expect_error(cyto_select(list(fs), Treatment = "Stim-C"),
               "'x' should be an object of class flowSet or GatingSet.")
  
  # Invalid variable names
  expect_error(cyto_select(fs, Treatment = "Stim-A", OvaConc = 0),
              "OvaConc is not a valid variable in pData(x).", fixed = TRUE)
  
  # Invalid levels for variables
  expect_error(cyto_select(fs, Treatment = "Stim-A", "OVAConc" = 10),
              "10 is not a valid level for OVAConc!", fixed = TRUE)
  
  expect_equal(cyto_select(fs, Treatment = "Stim-C"), 
               fs[c(17,18,19,20,21,22,23,24)])
  
  expect_equal(cyto_select(fs, Treatment = "Stim-A", OVAConc = c(0,0.5)),
               fs[c(1,2,7,8)])
  
  expect_equal(cyto_select(fs, 
                           list("Treatment" = "Stim-A", "OVAConc" = c(0,0.5))),
               fs[c(1,2,7,8)])
  
  # Filtered GatingSet will have different guid slot
  expect_equivalent(cyto_select(gs, Treatment = "Stim-C"), 
                    gs[c(17, 18, 19, 20, 21, 22, 23, 24)])
  
  expect_equivalent(cyto_select(gs, 
                                list("Treatment" = "Stim-A", 
                                     "OVAConc" = c(0,0.5))),
                    gs[c(1, 2, 7, 8)])
  
})

# cyto_group_by -----------------------------------------------------------------

test_that("cyto_group_by", {
  
  expect_error(cyto_group_by(list(fs)),
               "'x' should be an object of class flowSet or GatingSet.")
  
  expect_error(cyto_group_by(fs, c("Treatment","OVA")),
               "OVA is not a valid variable for this flowSet.")
  
  # flowSet
  expect_equal(cyto_group_by(fs,"Treatment"),
               list(fs[1:8],fs[9:16],fs[17:24],fs[25:32],fs[33]))
  
  # GatingSet - different guid slot for subset
  expect_equivalent(cyto_group_by(gs,"Treatment"),
                    list(gs[1:8],gs[9:16],gs[17:24],gs[25:32],gs[33]))
  
})

# cyto_sample ------------------------------------------------------------------

test_that("cyto_sample", {
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 1))), 2000)
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 0.5))), 1000)
})

# cyto_markers -----------------------------------------------------------------

# cyto_annotate ----------------------------------------------------------------

# cyto_compensate --------------------------------------------------------------

test_that("cyto_compensate", {
  
  # Write fs[[1]]@description$SPILL to csv file for testing
  spill <- fs[[1]]@description$SPILL
  rownames(spill) <- colnames(spill)
  write.csv(spill, "Test-Spillover-Matrix.csv")
  
  # flowFrame
  expect_equivalent(cyto_compensate(fs[[1]]),
                    getData(gs, "root")[[1]])
  expect_equivalent(cyto_compensate(fs[[1]], "Test-Spillover-Matrix.csv"),
                    getData(gs, "root")[[1]])
  
  # flowSet
  expect_equivalent(cyto_compensate(fs),
                    getData(gs, "root"))
  expect_equivalent(cyto_compensate(fs, "Test-Spillover-Matrix.csv"),
                    getData(gs, "root"))
  
  # GatingSet
  expect_equivalent(cyto_compensate(GatingSet(fs)),
                    gs)
  expect_equivalent(cyto_compensate(GatingSet(fs), "Test-Spillover-Matrix.csv"),
                    gs)
  
})

# cyto_transform ---------------------------------------------------------------

# cyto_inverse_transform -------------------------------------------------------


# Remove generated spillover matrix
base::unlink("Test-Spillover-Matrix.csv")