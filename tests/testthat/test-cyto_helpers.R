context("cyto-helpers")

# CYTO_LOAD --------------------------------------------------------------------

test_that("cyto_load", {
  
  # NCDFFLOWSET
  test_fs <- cyto_load(path = paste0(datadir, "/Activation")) 
  expect_true(is(test_fs, "ncdfFlowSet"))
  
})

# CYTO_SETUP -------------------------------------------------------------------

# CYTO_DETAILS -----------------------------------------------------------------

test_that("cyto_details", {
  
  pd <- data.frame("name" = paste0("Activation", "_", seq_len(33), ".fcs"),
                   "OVAConc" = c(rep(c(0,0,0.005,0.005,0.05,0.05,0.5,0.5), 4), 0),
                   "Treatment" = c(rep("Stim-A", 8),
                                   rep("Stim-B", 8),
                                   rep("Stim-C", 8),
                                   rep("Stim-D", 8),
                                   "NA"),
                   stringsAsFactors = FALSE)
  rownames(pd) <- paste0("Activation", "_", seq_len(33), ".fcs")
  pd$Treatment <- factor(pd$Treatment,
                         c(paste0("Stim", "-", c("A","B","C","D")), "NA"))
  pd$name <- factor(pd$name, levels = pd$name)
  expect_equal(cyto_details(gs),
               pd)
  
})

# CYTO_NAMES -------------------------------------------------------------------

test_that("cyto_names", {
  
  # flowFrame
  expect_equal(cyto_names(fs[[1]]), "Activation_1.fcs")
  
  # flowSet
  expect_equal(cyto_names(fs), paste0("Activation", "_", seq_len(33), ".fcs"))
  
  # GatingHierarchy
  expect_equal(cyto_names(gs[[1]]), "Activation_1.fcs")
  
  # GatingSet
  expect_equal(cyto_names(gs), paste0("Activation", "_", seq_len(33), ".fcs"))
  
})

# CYTO_CHECK -------------------------------------------------------------------

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

# CYTO_TRANSFORM ---------------------------------------------------------------

# CYTO_TRANSFORM_EXTRACT -------------------------------------------------------

test_that("cyto_transform_extract", {
  
  # NULL
  expect_equal(cyto_transform_extract(NULL), NULL)
  
  # transformList
  trns <- transformList(names(trans), lapply(trans, `[[`, "transform"))
  expect_equal(cyto_transform_extract(trns), trns)
  
  # transformerList
  expect_equal(cyto_transform_extract(trans),
               transformList(names(trans), lapply(trans, `[[`, "transform")))
  
  # transformerList - inverse
  expect_equal(cyto_transform_extract(trans, inverse = TRUE),
               transformList(names(trans), lapply(trans, `[[`, "inverse")))
  
})

# CYTO_EXTRACT -----------------------------------------------------------------

test_that("cyto_extract", {
  
  exp <- gs_pop_get_data(gs, "root")
  
  # flowFrame ------------------------------------------------------------------
  expect_equal(cyto_extract(fs[[1]]), fs[[1]])
  
  # flowSet --------------------------------------------------------------------
  expect_equal(cyto_extract(fs), fs)
  
  # GatingHierachy -------------------------------------------------------------
  expect_equal(cyto_extract(gs[[1]]), exp[[1]])
  
  # GatingSet ------------------------------------------------------------------
  expect_equal(cyto_extract(gs), exp)
  
})

# CYTO_CONVERT ----------------------------------------------------------------

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
               gs_pop_get_data(gs,"root")[[1]])
  expect_equal(cyto_convert(gs[[1]],"root" ,"flowFrame list"),
               list(gs_pop_get_data(gs,"root")[[1]]))
  expect_equal(cyto_convert(gs[[1]], "root", "flowSet"),
               flowSet(gs_pop_get_data(gs,"root")[[1]]))
  expect_equal(cyto_convert(gs[[1]],"root", "flowSet list"),
               list(flowSet(gs_pop_get_data(gs,"root")[[1]])))
  
  # GatingSet ------------------------------------------------------------------
  fr_exp <- as(gs_pop_get_data(gs,"root"),"flowFrame")
  if ("Original" %in% BiocGenerics::colnames(fr_exp)) {
    fr_exp <- suppressWarnings(
      fr_exp[, -match("Original", BiocGenerics::colnames(fr_exp))]
    )
  }
  fr_list_exp <- lapply(seq_len(length(gs)), 
                        function(y){gs_pop_get_data(gs,"root")[[y]]})
  expect_equal(cyto_convert(gs, "root", "flowFrame"),
               fr_exp)
  expect_equal(cyto_convert(gs, "root", "list of flowFrames"),
               fr_list_exp)
  expect_equal(cyto_convert(gs, "root", "flowFrame list"),
               list(fr_exp))
  expect_equal(cyto_convert(gs, "root", "flowSet"),
               gs_pop_get_data(gs,"root"))
  expect_equal(cyto_convert(gs, "root", "flowSet list"),
               list(gs_pop_get_data(gs,"root")))
  
})

# CYTO_FILTER ------------------------------------------------------------------

# CYTO_SELECT ------------------------------------------------------------------

test_that("cyto_select", {
  
  # Must be flowSet or GatingSet object
  expect_error(cyto_select(list(fs), Treatment = "Stim-C"),
               "'x' should be an object of class flowSet or GatingSet.")
  
  # Invalid variable names
  expect_error(cyto_select(fs, Treatment = "Stim-A", OvaConc = 0),
               "OvaConc is not a valid variable in cyto_details(x).", 
               fixed = TRUE)
  
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

# CYTO_GROUP_BY ----------------------------------------------------------------

test_that("cyto_group_by", {
  
  expect_error(cyto_group_by(list(fs)),
               "'x' should be an object of class flowSet or GatingSet.")
  
  expect_error(cyto_group_by(fs, c("Treatment","OVA")),
               "OVA is not a valid variable for this flowSet.")
  
  # flowSet
  exp <- list(fs[1:8],fs[9:16],fs[17:24],fs[25:32],fs[33])
  names(exp) <- c(paste0("Stim-",c("A","B","C","D")), "NA")
  expect_equal(cyto_group_by(fs,"Treatment"),
               exp)
  
  # GatingSet - different guid slot for subset
  exp <- list(gs[1:8],gs[9:16],gs[17:24],gs[25:32],gs[33])
  names(exp) <- c(paste0("Stim-",c("A","B","C","D")), "NA")
  expect_equivalent(cyto_group_by(gs,"Treatment"),
                    exp)
  
})

# CYTO_MERGE_BY ----------------------------------------------------------------

test_that("cyto_merge_by", {
  
})

# CYTO_SAMPLE ------------------------------------------------------------------

# CYTO_BARCODE -----------------------------------------------------------------

# CYTO_MARKERS_EDIT ------------------------------------------------------------

# CYTO_DETAILS_EDIT ------------------------------------------------------------

# CYTO_COMPENSATE --------------------------------------------------------------

# CYTO_NODES -------------------------------------------------------------------

# CYTO_CHANNEL_MATCH -----------------------------------------------------------

