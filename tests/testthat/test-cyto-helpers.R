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
  expect_equal(cyto_extract(gs, getData(gs, "root")), exp)
  
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
  expect_equal(cyto_convert(fs, "flowFrame list"), fr_list_exp)
  expect_equal(cyto_convert(fs, "flowSet list"), list(fs))
  
  # GatingHierarchy ------------------------------------------------------------
  expect_equal(cyto_convert(gs[[1]], "flowFrame", "root"),
               getData(gs,"root")[[1]])
  expect_equal(cyto_convert(gs[[1]], "flowFrame list", "root"),
               list(getData(gs,"root")[[1]]))
  expect_equal(cyto_convert(gs[[1]], "flowSet", "root"),
               flowSet(getData(gs,"root")[[1]]))
  expect_equal(cyto_convert(gs[[1]], "flowSet list", "root"),
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
  expect_equal(cyto_convert(gs, "flowFrame", "root"),
               fr_exp)
  expect_equal(cyto_convert(gs, "flowFrame list", "root"),
               fr_list_exp)
  expect_equal(cyto_convert(gs, "flowSet", "root"),
               getData(gs,"root"))
  expect_equal(cyto_convert(gs, "flowSet list", "root"),
               list(getData(gs,"root")))
  
})

# cyto_sample ------------------------------------------------------------------

test_that("cyto_sample", {
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 1))), 2000)
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 0.5))), 1000)
})


