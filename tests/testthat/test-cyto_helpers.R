context("cyto-helpers")

# CYTO_LOAD & CYTO_SETUP -------------------------------------------------------

# USED TO LOAD IN TEST DATASETS

# CYTO_DETAILS -----------------------------------------------------------------

test_that("cyto_details", {
  
  pd <- data.frame("name" = paste0("Activation", "_", seq_len(33), ".fcs"),
                   "OVAConc" = c(rep(c(0,0,5,5,50,50,500,500), 4), 0),
                   "Treatment" = c(rep("Stim-A", 8),
                                   rep("Stim-B", 8),
                                   rep("Stim-C", 8),
                                   rep("Stim-D", 8),
                                   "NA"),
                   stringsAsFactors = FALSE)
  rownames(pd) <- paste0("Activation", "_", seq_len(33), ".fcs")
  pd$OVAConc <- as.character(pd$OVAConc)
  
  # mac columns different order - test individually
  expect_true(all(colnames(cyto_details(gs)) %in% colnames(pd)))
  expect_equal(cyto_details(gs)$name, pd$name)
  expect_equal(cyto_details(gs)$OVAConc, pd$OVAConc)
  expect_equal(cyto_details(gs)$Treatment, pd$Treatment)

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
  expect_true(suppressMessages(cyto_check(flowSet_to_cytoset(fs))))
  
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
  expect_equivalent(cyto_convert(fs[[1]], "flowSet"), fs[1])
  expect_equivalent(cyto_convert(fs[[1]], "flowSet list"), list(fs[1]))
  
  # flowSet --------------------------------------------------------------------
  fr_exp <- as(fs,"flowFrame")
  if ("Original" %in% BiocGenerics::colnames(fr_exp)) {
    fr_exp <- suppressWarnings(
      fr_exp[, -match("Original", BiocGenerics::colnames(fr_exp))]
    )
  }
  fr_list_exp <- lapply(seq_len(length(fs)), function(x){fs[[x]]})
  names(fr_list_exp) <- cyto_names(fs)
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
  expect_equivalent(cyto_convert(gs[[1]], "root", "flowSet"),
               flowSet(gs_pop_get_data(gs,"root")[[1]]))
  expect_equivalent(cyto_convert(gs[[1]],"root", "flowSet list"),
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
  names(fr_list_exp) <- cyto_names(fs)
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

test_that("cyto_filter", {
  
  expect_error(cyto_filter(list(fs), Treatment = "Stim-C"),
               "'x' should be an object of class flowSet or GatingSet.")
  
  expect_equal(cyto_filter(fs, Treatment == "Stim-C"), 
               fs[c(17,18,19,20,21,22,23,24)])
  
  expect_equal(cyto_filter(fs, Treatment == "Stim-A", OVAConc %in% c(0,500)),
               fs[c(1,2,7,8)])
  
  # Filtered GatingSet will have different guid slot
  expect_equivalent(cyto_filter(gs, Treatment == "Stim-C"), 
                    gs[c(17, 18, 19, 20, 21, 22, 23, 24)])
  
  expect_equivalent(cyto_filter(gs, 
                                Treatment == "Stim-A", 
                                OVAConc %in% c(0,500)),
                    gs[c(1, 2, 7, 8)])
  
})

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
  
  expect_equal(cyto_select(fs, Treatment = "Stim-A", OVAConc = c(0,500)),
               fs[c(1,2,7,8)])
  
  expect_equal(cyto_select(fs, 
                           list("Treatment" = "Stim-A", "OVAConc" = c(0,500))),
               fs[c(1,2,7,8)])
  
  # Filtered GatingSet will have different guid slot
  expect_equivalent(cyto_select(gs, Treatment = "Stim-C"), 
                    gs[c(17, 18, 19, 20, 21, 22, 23, 24)])
  
  expect_equivalent(cyto_select(gs, 
                                list("Treatment" = "Stim-A", 
                                     "OVAConc" = c(0,500))),
                    gs[c(1, 2, 7, 8)])
  
})

# CYTO_GROUP_BY ----------------------------------------------------------------

test_that("cyto_group_by", {
  
  expect_error(cyto_group_by(list(fs)),
               "'x' should be an object of class flowSet or GatingSet.")
  
  expect_error(cyto_group_by(fs, c("Treatment","OVA")),
               "OVA is not a valid variable for this flowSet.")
  
  # flowSet
  exp <- list( fs[33], fs[1:8], fs[9:16], fs[17:24], fs[25:32])
  names(exp) <- c("NA", paste0("Stim-",c("A","B","C","D")))
  expect_equal(cyto_group_by(fs,"Treatment"), exp)
  
  # GatingSet - different guid slot for subset
  exp <- list(gs[33], gs[1:8], gs[9:16], gs[17:24], gs[25:32])
  names(exp) <- c("NA", paste0("Stim-",c("A","B","C","D")))
  expect_equivalent(cyto_group_by(gs,"Treatment"), exp)
  
})

# CYTO_MERGE_BY & CYTO_SPLIT ---------------------------------------------------

test_that("cyto_merge_by", {
  
  # CYTO_MERGE_BY
  fs <- cyto_extract(gs[c(1,2,9,10,17,18,25,26)], "root")
  fs <- cyto_barcode(fs)
  fs_list <- cyto_group_by(fs, "Treatment")
  fr_list <- lapply(seq_len(length(fs_list)), function(z){
    x <- cyto_convert(fs_list[[z]], "flowFrame")
    identifier(x) <- paste0("Stim-", c("A","B","C","D"))[z]
    return(x)
  })
  names(fr_list) <- paste0("Stim-", c("A","B","C","D"))
  expect_equal(cyto_merge_by(gs[c(1,2,9,10,17,18,25,26)], 
                             merge_by = "Treatment"),
               fr_list)
  
  # CYTO_SPLIT
  expect_equivalent(cyto_split(fr_list[[1]],
                          names = c("Activation_1.fcs",
                                    "Activation_2.fcs")),
                    list("Activation_1.fcs" = fs[[1]], 
                         "Activation_2.fcs" = fs[[2]]))
  
})

# CYTO_SAVE --------------------------------------------------------------------

test_that("cyto_save", {
  
  # GATINGSET
  cyto_save(gs[1], parent = "T Cells", save_as = "Saved-Samples")
  expect_true(dir.exists("Saved-Samples"))
  expect_true(grepl("Activation_1.fcs", list.files("Saved-Samples")))
  
  # GATINGHIERARCHY
  cyto_save(gs[[2]], parent = "T Cells", save_as = "Saved-Samples")
  expect_true(dir.exists("Saved-Samples"))
  expect_true(any(grepl("Activation_2.fcs", list.files("Saved-Samples"))))
  
})

# CYTO_SAMPLE ------------------------------------------------------------------

test_that("cyto_sample", {
  
  # flowFrame
  expect_equal(nrow(cyto_sample(fs[[1]], 0.1)), 200)
  expect_equal(nrow(cyto_sample(fs[[1]], 200)), 200)
  
  # flowSet
  fs_sample <- cyto_sample(fs, 0.1)
  fs_sample <- LAPPLY(seq_len(length(fs)), function(x){
    nrow(fs_sample[[x]])
  })
  expect_equal(fs_sample, rep(200, 33))
  
  fs_sample <- cyto_sample(fs, 200)
  fs_sample <- LAPPLY(seq_len(length(fs)), function(x){
    nrow(fs_sample[[x]])
  })
  expect_equal(fs_sample, rep(200, 33))
  
  # list (unexported - used in cyto_plot only)
  
})

# CYTO_BARCODE -----------------------------------------------------------------

test_that("cyto_barcode", {
  
  # SAMPLES
  fs <- fs[c(1,2)]
  barcode <- lapply(seq_len(2), function(z){
    matrix(rep(z, nrow(fs[[z]])),
           ncol = 1,
           dimnames = list(NULL, "Sample ID"))
  })
  fs_barcode <- fsApply(fs, function(x){
    fr_append_cols(x, barcode[[match(cyto_names(x), cyto_names(fs))]])
  })
  expect_equal(cyto_barcode(fs[1:2]),
               fs_barcode)
  
  # EVENTS
  fs <- fs[c(1,2)]
  barcode <- split(seq_len(4000), c(rep(1, 2000),
                                    rep(2, 2000)))
  barcode <- lapply(barcode, function(z){
    matrix(z,
           ncol = 1,
           dimnames = list(NULL, "Event ID"))
  })
  fs_barcode <- fsApply(fs, function(x){
    fr_append_cols(x, barcode[[match(cyto_names(x), cyto_names(fs))]])
  })
  expect_equal(cyto_barcode(fs[1:2],
                            type = "events"),
               fs_barcode)
  
})

# CYTO_MARKERS_EDIT ------------------------------------------------------------

test_that("cyto_markers_edit", {
  
  mock_data_editor <- mock(data.frame("channel" = cyto_channels(fs),
                                      "marker" = pData(parameters(fs[[1]]))$desc,
                                      stringsAsFactors = FALSE,
                                      row.names = NULL))
  testthat::with_mock(data_editor = mock_data_editor,
                      expect_equal(cyto_markers_edit(fs),
                                   fs))
  
})

# CYTO_DETAILS_EDIT ------------------------------------------------------------

test_that("cyto_details_edit", {
  
  pd <- cyto_details(fs)
  rownames(pd) <- NULL
  mock_data_editor <- mock(cyto_details(fs))
  testthat::with_mock(data_editor = mock_data_editor,
                      expect_equivalent(cyto_details_edit(fs),
                                   fs))
  
})

# CYTO_COMPENSATE --------------------------------------------------------------

test_that("cyto_compensate", {
  
  # Write fs[[1]]@description$SPILL to csv file for testing
  spill <- fs[[1]]@description$SPILL
  rownames(spill) <- colnames(spill)
  write.csv(spill, "Test-Spillover-Matrix.csv")
  
  # flowFrame
  expect_equivalent(cyto_compensate(fs[[1]]),
                    gs_pop_get_data(gs, "root")[[1]])
  expect_equivalent(cyto_compensate(fs[[1]], "Test-Spillover-Matrix.csv"),
                    gs_pop_get_data(gs, "root")[[1]])
  
  # flowSet
  expect_equivalent(cyto_compensate(fs),
                    gs_pop_get_data(gs, "root"))
  expect_equivalent(cyto_compensate(fs, "Test-Spillover-Matrix.csv"),
                    gs_pop_get_data(gs, "root"))
  
  # GatingSet
  expect_equivalent(cyto_compensate(GatingSet(fs)),
                    gs)
  expect_equivalent(cyto_compensate(GatingSet(fs), "Test-Spillover-Matrix.csv"),
                    gs)
  
})

# CYTO_NODES -------------------------------------------------------------------

test_that("cyto_nodes", {
  
  # GATINGSET
  expect_equal(cyto_nodes(gs), gs_get_pop_paths(gs))
  #GATINHIERARCHY
  expect_equal(cyto_nodes(gs[[1]]), gh_get_pop_paths(gs[[1]]))
  
})

# CYTO_CHANNEL_MATCH -----------------------------------------------------------

test_that("cyto_channel_match", {
  
  channel_match <- data.frame("name" = c("Compensation-7AAD.fcs",
                                          "Compensation-AF700.fcs",
                                          "Compensation-APC-Cy7.fcs",
                                          "Compensation-APC.fcs",
                                          "Compensation-FITC.fcs",
                                          "Compensation-PE.fcs",
                                          "Compensation-Unstained.fcs"),
                               "channel" = c("7-AAD-A",
                                             "Alexa Fluor 700-A",
                                             "APC-Cy7-A",
                                             "Alexa Fluor 647-A",
                                             "Alexa Fluor 488-A",
                                             "PE-A",
                                             "Unstained"))
  mock_edit <- mockery::mock(channel_match)
  mock_channel_match <- testthat::with_mock(
    edit = mock_edit,
    cyto_channel_match(fs_comp, 
                       menu = FALSE,
                       save_as = "Channel-Match"))
  expect_equal(read.csv("Channel-Match.csv"), channel_match)
  expect_equal(mock_channel_match, channel_match)
  
})

# CYTO_EMPTY -------------------------------------------------------------------

test_that("cyto_empty", {
  
  chans <- cyto_channels(fs)
  empty_flowFrame <- matrix(0,
                            ncol = length(chans),
                            nrow = 1,
                            byrow = TRUE)
  colnames(empty_flowFrame) <- chans
  empty_flowFrame <- flowFrame(empty_flowFrame)
  empty_flowFrame <- empty_flowFrame[-1, ]
  identifier(empty_flowFrame) <- "Activation_1.fcs"
  expect_error(cyto_empty("Activation_1.fcs"),
               "Supply the names of the channels to include in the flowFrame.")
  expect_equal(cyto_empty("Activation_1.fcs",
                          channels = cyto_channels(fs)),
               empty_flowFrame)
  
})

# REMOVE GENERATED FILES -------------------------------------------------------

# Spillover matrix
base::unlink("Test-Spillover-Matrix.csv")

# Channel match 
base::unlink("Channel-Match.csv")

# Experiment markers
base::unlink(paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Markers.csv"))

# Experiment details
base::unlink(paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Details.csv"))

# Saved samples
base::unlink("Saved-Samples", recursive = TRUE)
