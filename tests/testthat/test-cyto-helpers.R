# CYTO_IMPORT ------------------------------------------------------------------

# CYTO_EXPORT ------------------------------------------------------------------

# CYTO_LOAD --------------------------------------------------------------------

# CYTO_CLEAN -------------------------------------------------------------------

# CYTO_SETUP -------------------------------------------------------------------

# CYTO_CLASS -------------------------------------------------------------------

test_that("cyto_class", {
  # GATINGSET
  expect_equal(cyto_class(gs), "GatingSet", ignore_attr = TRUE)
  expect_equal(cyto_class(gs, class = FALSE), "GatingSet")
  expect_false(cyto_class(gs, "GatingHierarchy"))
  expect_true(cyto_class(gs, "GatingSet"))
  expect_true(cyto_class(gs, c("GatingHierarchy", "GatingSet")))
  # GATINGHIERARCHY
  expect_equal(cyto_class(gs[[1]]), "GatingHierarchy", ignore_attr = TRUE)
  expect_equal(
    cyto_class(gs[[1]], class = FALSE),
    c("GatingHierarchy", "GatingSet")
  )
  expect_true(cyto_class(gs[[1]], "GatingHierarchy"))
  expect_true(cyto_class(gs[[1]], "GatingSet"))
  expect_true(cyto_class(gs[[1]], c("GatingHierarchy", "GatingSet")))
  # CYTOSET
  expect_equal(cyto_class(cs), "cytoset", ignore_attr = TRUE)
  expect_equal(cyto_class(cs, class = FALSE), c("cytoset", "flowSet"))
  expect_false(cyto_class(cs, c("flowFrame", "GatingHierarchy", "GatingSet")))
  expect_true(cyto_class(cs, "cytoset"))
  expect_true(cyto_class(cs, "flowSet"))
  expect_false(cyto_class(cs, "flowSet", TRUE))
  # CYTOFRAME
  expect_equal(cyto_class(cs[[1]]), "cytoframe", ignore_attr = TRUE)
  expect_equal(cyto_class(cs[[1]], class = FALSE), c("cytoframe", "flowFrame"))
  expect_false(cyto_class(cs[[1]], c("flowSet", "GatingHierarchy", "GatingSet")))
  expect_true(cyto_class(cs[[1]], "cytoframe"))
  expect_true(cyto_class(cs[[1]], "flowFrame"))
  expect_false(cyto_class(cs[[1]], "flowFrame", class = TRUE))
})

# CYTO_DETAILS -----------------------------------------------------------------

test_that("cyto_details", {
  # DETAILS
  pd <- data.frame(
    "name" = paste0("Activation_", 1:33, ".fcs"),
    "Treatment" = c(
      rep("Stim-A", 8),
      rep("Stim-B", 8),
      rep("Stim-C", 8),
      rep("Stim-D", 8),
      "NA"
    ),
    "OVAConc" = as.character(
      c(rep(c(0, 0, 5, 5, 50, 50, 500, 500), 4), 0)
    ),
    row.names = paste0("Activation_", 1:33, ".fcs"),
    stringsAsFactors = FALSE
  )
  # GATINGSET
  expect_equal(cyto_details(gs)[, colnames(pd)], pd)
  # CYTOSET
  expect_equal(cyto_details(cs)[, colnames(pd)], pd)
  # GATINGHIERARCHY
  expect_equal(
    cyto_details(gs[[1]]),
    pd[1, ]
  )
  # CYTOFRAME
  expect_error(cyto_details(cs[[1]]),
               "'cyto_details' cannot be extracted from cytoframe objects!")
  # CONVERT TO FACTOR & DROP ROW NAMES
  pd$name <- factor(pd$name)
  pd$Treatment <- factor(pd$Treatment,
    levels = c(
      NA,
      "Stim-A",
      "Stim-B",
      "Stim-C",
      "Stim-D"
    )
  )
  pd$OVAConc <- as.numeric(pd$OVAConc)
  rownames(pd) <- NULL
  expect_equal(
    cyto_details(gs,
      drop = TRUE,
      convert = TRUE,
      factor = TRUE
    )[, colnames(pd)],
    pd
  )
})

# CYTO_NAMES -------------------------------------------------------------------

test_that("cyto_names", {
  # GATINGSET
  expect_equal(cyto_names(gs), paste0("Activation_", 1:33, ".fcs"))
  # GATINHIERRACHY
  expect_equal(cyto_names(gs[[1]]), "Activation_1.fcs")
  # CYTOSET
  expect_equal(cyto_names(cs), paste0("Activation_", 1:33, ".fcs"))
  # CYTOFRAME
  expect_error(cyto_names(cs[[1]]), 
               "'cyto_names' cannot be extracted from cytoframe objects!")
  # REPLACE - GATINGSET/CYTOSET
  cyto_names(gs)[1] <- "Activation_1.fcs"
  expect_equal(cyto_names(gs), paste0("Activation_", 1:33, ".fcs"))
  cyto_names(cs)[1] <- "Activation_1.fcs"
  expect_equal(cyto_names(cs), paste0("Activation_", 1:33, ".fcs"))
  # REPLACE - CYTOFRAME
  expect_error(cyto_names(cs[[1]]) <- "Test.fcs",
               "'cyto_names' cannot be replaced for objects of class cytoframe!")
})

# CYTO_NAMES_PARSE -------------------------------------------------------------

test_that("cyto_names_parse", {
  pd <- cyto_details(gs)
  pd <- cbind(
    "experiment" = "Activation",
    "sample" = as.character(1:33),
    pd
  )
  gs_clone <- cyto_copy(gs)
  cyto_names_parse(gs_clone, c("experiment", "sample"))
  cyto_details <- cyto_details(gs_clone)
  expect_equal(
    cyto_details[, colnames(pd)],
    pd
  )
})

# CYTO_TRANSFORM ---------------------------------------------------------------

# .CYTO_TRANSFORM --------------------------------------------------------------

# CYTO_TRANSFORM_EXTRACT -------------------------------------------------------

test_that("cyto_transform_extract", {
  transformer_list <- cyto_transformer_extract(gs)
  transform_list <- cyto_transform_extract(transformer_list)
  expect_equal(class(transform_list), "transformList", ignore_attr = TRUE)
  transform_list <- cyto_transform_extract(transformer_list, inverse = TRUE)
  expect_equal(class(transform_list), "transformList", ignore_attr = TRUE)
})

# CYTO_DATA_EXTRACT ------------------------------------------------------------

test_that("cyto_data_extract", {

  # GATINGSET -> CYTOSET
  cs <- cyto_data_extract(gs_sub,
    parent = "root",
    copy = TRUE,
    select = list(Treatment = "Stim-D"),
    channels = c("CD4", "CD8")
  )[["root"]]
  expect_true(cyto_class(cs, "cytoset"))
  expect_equal(cyto_channels(cs), c(
    "Alexa Fluor 700-A",
    "Alexa Fluor 488-A"
  ))

  # GATINGSET -> LIST OF MATRICES
  cs_raw <- cyto_data_extract(gs_sub,
    parent = "root",
    copy = TRUE,
    select = list(Treatment = "Stim-D"),
    channels = c("CD4", "CD8"),
    format = "matrix",
    markers = TRUE
  )[["root"]]
  expect_true(cyto_class(cs_raw, "list"))
  expect_true(cyto_class(cs_raw[[1]], "matrix"))
  expect_equal(colnames(cs_raw[[1]]), c("CD4", "CD8"))
})

# CYTO_FILTER ------------------------------------------------------------------

test_that("cyto_filter", {

  # CLASS CHECK
  expect_error(
    cyto_filter(cs_sub[[1]],
      Treatment = "Stim-C"
    ),
    "'x' should be an object of class cytoset or GatingSet."
  )

  # CYTOSET
  res <- cyto_filter(
    cs_sub,
    Treatment == "Stim-C"
  )
  expect_equal(length(res), 8)
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", 17:24, ".fcs")
  )

  # GATINGSET
  res <- cyto_filter(
    gs_sub,
    Treatment == "Stim-A",
    OVAConc %in% c(0, 500)
  )
  expect_equal(length(res), 4)
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", c(1, 2, 7, 8), ".fcs")
  )
  expect_equal(
    cyto_details(res)$OVAConc,
    as.character(c(0, 0, 500, 500))
  )

  # ERROR
  expect_message(
    cyto_filter(
      cs_sub,
      Treatment == "Stim-E"
    ),
    "No samples match the filtering criteria. Returning all samples."
  )
})


# CYTO_SELECT ------------------------------------------------------------------

test_that("cyto_select", {

  # CLASS CHECK
  expect_error(
    cyto_select(cs_sub[[1]],
      Treatment = "Stim-C"
    ),
    "'x' should be an object of class cytoset or GatingSet."
  )

  # INVALID VARIABLE
  expect_error(cyto_select(cs_sub,
    Treatment = "Stim-A",
    OvaConc = 0
  ),
  "OvaConc is not a valid variable in cyto_details(x).",
  fixed = TRUE
  )

  # INVALID LEVEL
  expect_error(cyto_select(cs_sub,
    Treatment = "Stim-A",
    "OVAConc" = 10
  ),
  "10 is not a valid level for OVAConc!",
  fixed = TRUE
  )

  # INDEX
  res <- cyto_select(cs_sub, 1:8)
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", 1:8, ".fcs")
  )

  # EXCLUDE
  res <- cyto_select(cs_sub, 1:8, exclude = TRUE)
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", 9:33, ".fcs")
  )

  # VARIABLES
  res <- cyto_select(cs_sub,
    Treatment = "Stim-A",
    OVAConc = c(0, 500)
  )
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", c(1, 2, 7, 8), ".fcs")
  )

  # LIST
  res <- cyto_select(
    cs_sub,
    list(
      Treatment = "Stim-A",
      OVAConc = c(0, 500)
    )
  )
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", c(1, 2, 7, 8), ".fcs")
  )
})

# CYTO_GROUPS ------------------------------------------------------------------

test_that("cyto_groups", {

  # CLASS CHECK
  expect_error(
    cyto_groups(cs_sub[[1]]),
    "'x' should be an object of class cytoset or GatingSet."
  )

  # NAME
  expect_equal(
    cyto_groups(cs_sub, "name"),
    paste0("Activation_", 1:33, ".fcs")
  )

  # TREATMENT
  expect_equal(
    cyto_groups(cs_sub, "Treatment"),
    c("NA", "Stim-A", "Stim-B", "Stim-C", "Stim-D")
  )

  # DETAILS
  pd <- cyto_details(cs_sub)
  expect_equal(
    cyto_groups(cs_sub,
      "Treatment",
      details = TRUE
    ),
    list(
      "NA" = pd[33, ],
      `Stim-A` = pd[1:8, ],
      `Stim-B` = pd[9:16, ],
      `Stim-C` = pd[17:24, ],
      `Stim-D` = pd[25:32, ]
    )
  )

  # ORDERED GROUPS
  expect_equal(
    cyto_groups(
      cs_sub,
      list(
        Treatment = c(
          "Stim-A",
          "Stim-B",
          "Stim-C",
          "Stim-D",
          "NA"
        ),
        OVAConc = c(
          500,
          50,
          5,
          0
        )
      )
    ),
    c(
      paste(c("Stim-A"), c(500, 50, 5, 0)),
      paste(c("Stim-B"), c(500, 50, 5, 0)),
      paste(c("Stim-C"), c(500, 50, 5, 0)),
      paste(c("Stim-D"), c(500, 50, 5, 0)),
      "NA 0"
    )
  )
})

# CYTO_SORT_BY -----------------------------------------------------------------

test_that("cyto_sort_by", {

  # CYTOSET
  res <- cyto_sort_by(
    cs_sub,
    list(Treatment = c(
      "Stim-A",
      "Stim-C",
      "Stim-B",
      "Stim-D",
      "NA"
    ))
  )
  expect_equal(
    cyto_details(res)$name,
    paste0("Activation_", c(1:8, 17:24, 9:16, 25:33), ".fcs")
  )
})

# CYTO_GROUP_BY ----------------------------------------------------------------

test_that("cyto_group_by", {

  # CYTOSET
  res <- cyto_group_by(
    cs_sub,
    list(Treatment = c(
      "Stim-A",
      "Stim-C",
      "Stim-B",
      "Stim-D",
      "NA"
    ))
  )
  expect_true(cyto_class(res, "list"))
  expect_true(all(unlist(lapply(res, "cyto_class", "flowSet"))))
  expect_equal(names(res), c(
    "Stim-A",
    "Stim-C",
    "Stim-B",
    "Stim-D",
    "NA"
  ))
})

# CYTO_MERGE_BY ----------------------------------------------------------------

test_that("cyto_merge_by", {

  # CYTOSET
  res <- cyto_merge_by(cs_sub,
    merge_by = list(Treatment = c(
      "Stim-A",
      "Stim-C",
      "Stim-B",
      "Stim-D",
      "NA"
    )),
    select = list(OVAConc = 0)
  )
  expect_true(cyto_class(res, "list"))
  expect_equal(names(res), c(
    "Stim-A",
    "Stim-C",
    "Stim-B",
    "Stim-D",
    "NA"
  ))
  expect_true(all(unlist(lapply(res, "cyto_class", "flowSet", TRUE))))
  expect_equal(unlist(lapply(res, "cyto_names")),
               c(
                 "Stim-A",
                 "Stim-C",
                 "Stim-B",
                 "Stim-D",
                 "NA"
               ))
})

# CYTO_SPLIT -------------------------------------------------------------------

test_that("cyto_split", {

  # MERGE
  cf_list <- cyto_merge_by(cs_sub[1:2])

  # CLASS CHECK
  expect_error(cyto_split(cs_sub),
    "cyto_split() requires a cytoframe object.",
    fixed = TRUE
  )

  # MERGE CHECK
  expect_error(cyto_split(cs_sub[[1]]),
    "Merged samples must be barcoded in cyto_merge().",
    fixed = TRUE
  )

  # NAMES CHECK
  expect_error(
    cyto_split(cf_list[[1]],
      names = "Test1.fcs"
    ),
    "Supply a name for each file."
  )

  # SPLIT
  res <- cyto_split(cf_list[[1]], names = c(
    "Test1.fcs",
    "Test2.fcs"
  ))
  expect_equal(cyto_names(res), c(
    "Test1.fcs" = "Test1.fcs",
    "Test2.fcs" = "Test2.fcs"
  ))
})

# CYTO_SAVE --------------------------------------------------------------------

# CYTO_SAMPLE ------------------------------------------------------------------

test_that("cyto_sample", {

  skip("cyto_sample")
  
  # SAMPLED IN HELPER-LIB
  expect_equal(nrow(gs_sub), structure(rep(list(2000), 33),
    names = paste0("Activation_", 1:33, ".fcs")
  ))
})

# CYTO_SAMPLE_TO_NODE ----------------------------------------------------------

test_that("cyto_sample_to_node", {

  skip("cyto_sample_to_node")
  
  # COPY
  gs_sub_copy <- cyto_copy(gs_sub[1:33])

  # ERROR - NO EVENTS
  expect_error(cyto_sample_to_node(gs_sub_copy,
    node = "T Cells",
    count = 100
  ),
  paste0("The following samples do not contain any events in the specified node:",
  paste0("\n", "Activation_33.fcs")),
  fixed = TRUE
  )

  # MIN COUNT
  gs_sub_copy <- cyto_sample_to_node(gs_sub_copy,
    node = "T Cells"
  )
  expect_equal(mean(unlist(nrow(cyto_extract(gs_sub_copy, "T Cells")))),
    176,
    tolerance = 2
  )
})

# CYTO_BARCODE -----------------------------------------------------------------

test_that("cyto_barcode", {
  
  # CYTOSET
  cs_sub_copy <- cyto_copy(cs_sub[1:4])
  res <- cyto_barcode(cs_sub_copy,
                      type = "both")
  expect_true("Sample ID" %in% cyto_channels(res))
  expect_true("Event ID" %in% cyto_channels(res))
  
})

# CYTO_MARKERS_EDIT ------------------------------------------------------------

test_that("cyto_markers_edit", {
  # MOCK DATA_EDIT
  stub(cyto_markers_edit,
       "data_edit",
       data.frame("channel" = cyto_channels(cs),
                    "marker" = pData(parameters(cs[[1]]))$desc,
                    stringsAsFactors = FALSE,
                    row.names = NULL))
  cyto_markers_edit(cs)
  expect_equal(cyto_markers(cs),
               c("Alexa Fluor 488-A" = "CD8",
                 "PE-A" = "Va2",
                 "7-AAD-A" = "CD69",
                 "Alexa Fluor 405-A" = "Hoechst-405",
                 "Alexa Fluor 430-A" = "Hoechst-430",
                 "Alexa Fluor 647-A" = "CD44",
                 "Alexa Fluor 700-A" = "CD4",
                 "APC-Cy7-A" = "CD11c"))
})

# CYTO_DETAILS_EDIT ------------------------------------------------------------

test_that("cyto_details_edit", {
  # MOCK DATA_EDIT
  pd <- cyto_details(cs)
  stub(cyto_details_edit,
       "data_edit",
       pd)
  cyto_details_edit(cs)
  expect_equal(cyto_details(cs)[, colnames(pd)], pd)
  
})

# CYTO_COMPENSATE --------------------------------------------------------------

# CYTO_NODES -------------------------------------------------------------------

test_that("cyto_nodes", {
  
  expect_snapshot_output(cyto_nodes(gs))
  expect_snapshot_output(cyto_nodes(gs, path = "auto"))
  expect_snapshot_output(cyto_nodes(gs[[1]]))
  expect_snapshot_output(cyto_nodes(gs[[1]], path = "auto"))
  
})

# CYTO_NODES_CONVERT -----------------------------------------------------------

test_that("cyto_nodes_convert", {
  
  expect_equal(cyto_nodes_convert(gs,
                                  nodes = "CD69+ CD4 T Cells",
                                  anchor = "T Cells"),
               "CD69+ CD4 T Cells")
  
})

# CYTO_NODES_ANCESTOR ----------------------------------------------------------

test_that("cyto_nodes_ancestor", {
  
  expect_equal(cyto_nodes_ancestor(gs,
                                   nodes = c("CD4 T Cells",
                                             "CD8 T Cells")),
               "T Cells")
  
  expect_equal(cyto_nodes_ancestor(gs,
                                   nodes = c("CD4 T Cells",
                                             "Dendritic Cells"),
                                   path = "full"),
               "/Cells/Single Cells/Live Cells")
  
})

# CYTO_EMPTY -------------------------------------------------------------------

# CYTO_SPILLOVER_EXTRACT -------------------------------------------------------

test_that("cyto_spillover_extract", {
  
  expect_snapshot_output(cyto_spillover_extract(gs[1]))
  expect_snapshot_output(cyto_spillover_extract(gs[[1]]))
  expect_snapshot_output(cyto_spillover_extract(cs[1]))
  expect_snapshot_output(cyto_spillover_extract(cs[[1]]))
  
})

# CYTO_CALIBRATE ---------------------------------------------------------------

test_that("cyto_calibrate", {
  
  # QUANTILE CALIBRATION
  cyto_calibrate(gs_sub[29:33],
                 parent = "root")
  expect_true(
    file.exists(
      paste0(temp_dir, .Platform$file.sep, "cyto_calibrate.rds")
      )
    )
  # RECALL
  expect_snapshot_output(.cyto_calibrate_recall())
  
  # RANGE CALIBRATION
  cyto_calibrate(gs_sub[29:33],
                 parent = "root",
                 type = "range")
  expect_true(
    file.exists(
      paste0(temp_dir, .Platform$file.sep, "cyto_calibrate.rds")
    )
  )
  # RECALL
  expect_snapshot_output(.cyto_calibrate_recall())
  
  # RESET
  cyto_calibrate_reset()
  expect_false(
    file.exists(
      paste0(temp_dir, .Platform$file.sep, "cyto_calibrate.rds")
    )
  )
  
})

# CYTO_APPLY -------------------------------------------------------------------

test_that("cyto_apply", {
  
  # CYTOSET
  expect_snapshot_output(
    cyto_apply(gs_sub[29:32],
               "cyto_names",
               input = "cytoset")
  )
  
  # CYTOFRAME
  expect_snapshot_output(
    cyto_apply(gs_sub[29:32],
               function(cf){
                 colMeans(exprs(cf))
               },
               input = "cytoframe",
               inverse = TRUE,
               channels = c("CD44", "CD69"))
  )
  
  # MATRIX
  expect_snapshot_output(
    cyto_apply(gs_sub[29:32],
               "nrow",
               input = "matrix",
               parent = "T Cells")
  )
  
  # COLUMN
  expect_snapshot_output(
    cyto_apply(gs_sub[29:32],
               function(z,
                        probs,
                        na.rm){
                 round(
                   quantile(z, 
                            probs = probs,
                            na.rm = na.rm),
                   2)
               },
               probs = c(0.5, 0.75),
               na.rm = TRUE,
               input = "column", 
               parent = c("CD4 T Cells", "CD8 T Cells"),
               channels = c("CD44", "CD69"))
  )
  
})

# CYTO_CBIND -------------------------------------------------------------------

test_that("cyto_cbind", {
  
  cs_copy <- cyto_copy(cs_sub[1:4])
  expect_error(cyto_cbind(cs_copy[[1]]),
               "'cols' must be a matrix!", fixed = TRUE)
  new_col <- matrix(1:8000,
                     ncol = 1,
                     dimnames = list(NULL, "Test"))
  cs_copy <- cyto_cbind(cs_copy, new_col)
  expect_true("Test" %in% cyto_channels(cs_copy))
  
})

# CLEAN UP ---------------------------------------------------------------------

# EXPERIMENT MARKERS
base::unlink(paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Markers.csv"))

# EXPERIMENT DETAILS
base::unlink(paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Details.csv"))
