
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
  expect_equal(cyto_class(gs[[1]], class = FALSE), 
               c("GatingHierarchy", "GatingSet"))
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
    "Treatment" = c(rep("Stim-A", 8),
                    rep("Stim-B", 8),
                    rep("Stim-C", 8),
                    rep("Stim-D", 8),
                    "NA"),
    "OVAConc" = as.character(
      c(rep(c(0,0,5,5,50,50,500,500), 4), 0)
      ),
    row.names = paste0("Activation_", 1:33, ".fcs"),
    stringsAsFactors = FALSE
  )
  # GATINGSET
  expect_equal(cyto_details(gs), pd)
  # CYTOSET
  expect_equal(cyto_details(cs), pd)
  # GATINGHIERARCHY
  expect_equal(cyto_details(gs[[1]]), 
               pd[1, ])
  # CYTOFRAME
  expect_equal(cyto_details(cs[[1]]), 
               data.frame("name" = "Activation_1.fcs",
                          row.names = "Activation_1.fcs",
                          stringsAsFactors = FALSE))
  # CONVERT TO FACTOR & DROP ROW NAMES
  pd$name <- factor(pd$name)
  pd$Treatment <- factor(pd$Treatment, 
                         levels = c(NA, 
                                    "Stim-A", 
                                    "Stim-B",
                                    "Stim-C",
                                    "Stim-D"))
  pd$OVAConc <- as.numeric(pd$OVAConc)
  rownames(pd) <- NULL
  expect_equal(cyto_details(gs, 
                            drop = TRUE,
                            convert = TRUE,
                            factor = TRUE),
               pd)
})

# CYTO_NAMES -------------------------------------------------------------------

# CYTO_NAMES_PARSE -------------------------------------------------------------

test_that("cyto_names_parse", {
  
  pd <- cyto_details(gs)
  pd <- cbind("experiment" = "Activation",
              "sample" = as.character(1:33),
              pd)
  gs_clone <- cyto_copy(gs)
  cyto_names_parse(gs_clone, c("experiment", "sample"))
  cyto_details <- cyto_details(gs_clone)
  expect_equal(cyto_details, 
               pd[, match(colnames(cyto_details), colnames(pd))])
  
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

# CYTO_EXTRACT -----------------------------------------------------------------


