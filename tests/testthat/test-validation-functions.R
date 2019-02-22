context("Validation Functions")

## .file_wd_check --------------------------------------------------------------

test_that(".file_wd_check", {
  expect_false(.file_wd_check("Test.csv"))
})

## cyto_channel_check ----------------------------------------------------------

test_that("cyto_channel_check flowFrame method", {
  expect_equal(cyto_channel_check(fs[[1]], 
                                  "CD4", 
                                  plot = TRUE), 
               "Alexa Fluor 700-A")
  expect_equal(cyto_channel_check(fs[[1]], 
                                  c("CD4", "Alexa Fluor 488-A"), 
                                  plot = TRUE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A"))
  expect_equal(cyto_channel_check(fs[[1]], 
                                  c("CD4", "CD8"), 
                                  plot = TRUE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A"))
  expect_error(cyto_channel_check(fs[[1]],
                                  c("CD4", "CD8", "CD69"), 
                                  plot = TRUE), 
               "Invalid number of supplied channels.", 
               fixed = TRUE)
  expect_equal(cyto_channel_check(fs[[1]], 
                                  c("CD4", "CD8", "CD69"), 
                                  plot = FALSE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "7-AAD-A"))
  expect_equal(cyto_channel_check(fs[[1]], 
                                  c("CD4", "Alexa Fluor 488-A", "CD69"), 
                                  plot = FALSE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "7-AAD-A"))
  expect_error(cyto_channel_check(fs[[1]], 
                                  c("Alexa Fluor 48-A", "CD69"), 
                                  plot = TRUE),
               "Alexa Fluor 48-A is not a valid channel for this flowFrame.")
  expect_equal(cyto_channel_check(fs[[1]], 
                                  channels = c(1, 2), 
                                  plot = TRUE), 
               colnames(fs[[1]])[1:2])
})

test_that("cyto_channel_check flowSet method", {
  expect_equal(cyto_channel_check(fs, 
                                  "CD4", 
                                  plot = TRUE), 
               "Alexa Fluor 700-A")
  expect_equal(cyto_channel_check(fs,
                                  c("CD4", "Alexa Fluor 488-A"), 
                                  plot = TRUE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A"))
  expect_equal(cyto_channel_check(fs, 
                                  c("CD4", "CD8"), 
                                  plot = TRUE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A"))
  expect_error(cyto_channel_check(fs, 
                                  c("CD4", "CD8", "CD69"), 
                                  plot = TRUE), 
               "Invalid number of supplied channels.", fixed = TRUE)
  expect_equal(cyto_channel_check(fs,
                                  c("CD4", "CD8", "CD69"), 
                                  plot = FALSE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "7-AAD-A"))
  expect_equal(cyto_channel_check(fs, 
                                  c("CD4", "Alexa Fluor 488-A", "CD69"), 
                                  plot = FALSE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "7-AAD-A"))
  expect_error(cyto_channel_check(fs,
                                  c("Alexa Fluor 48-A", "CD69"),
                                  plot = TRUE), 
               "Alexa Fluor 48-A is not a valid channel for this flowFrame.")
  expect_equal(cyto_channel_check(fs, 
                                  channels = c(1, 2), 
                                  plot = TRUE), 
               colnames(fs[[1]])[1:2])
})

test_that("cyto_channel_check GatingSet method", {
  expect_equal(cyto_channel_check(gs, 
                                  "CD4", 
                                  plot = TRUE), 
               "Alexa Fluor 700-A")
  expect_equal(cyto_channel_check(gs, 
                                  c("CD4", "Alexa Fluor 488-A"), 
                                  plot = TRUE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A"))
  expect_equal(cyto_channel_check(gs, 
                                  c("CD4", "CD8"), 
                                  plot = TRUE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A"))
  expect_error(cyto_channel_check(gs, 
                                  c("CD4", "CD8", "CD69"), 
                                  plot = TRUE), 
               "Invalid number of supplied channels.", fixed = TRUE)
  expect_equal(cyto_channel_check(gs,
                                  c("CD4", "CD8", "CD69"),
                                  plot = FALSE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "7-AAD-A"))
  expect_equal(cyto_channel_check(gs,
                                  c("CD4", "Alexa Fluor 488-A", "CD69"), 
                                  plot = FALSE), 
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "7-AAD-A"))
  expect_error(cyto_channel_check(gs,
                                  c("Alexa Fluor 48-A", "CD69"), 
                                  plot = TRUE), 
               "Alexa Fluor 48-A is not a valid channel for this flowFrame.")
  expect_equal(cyto_channel_check(gs,
                                  channels = c(1, 2), 
                                  plot = TRUE), 
               colnames(fs[[1]])[1:2])
})

## .cyto_gate_type_check -------------------------------------------------------

test_that(".cyto_gate_type_check", {
  expect_equal(.cyto_gate_type_check(type = "r", 
                                     alias = c("A", "B", "C")), 
               c("rectangle", "rectangle", "rectangle"))
  expect_equal(.cyto_gate_type_check(type = "p", 
                                     alias = c("A", "B", "C")), 
               c("polygon", "polygon", "polygon"))
  expect_equal(.cyto_gate_type_check(type = "q", 
                                     alias = c("A", "B", "C", "D")), 
               "quadrant")
  expect_error(.cyto_gate_type_check(type = "q", 
                                     alias = c("A", "B", "C")), 
               "Supply the names of 4 poulations to alias for quadrant gates.", 
               fixed = TRUE)
  expect_error(.cyto_gate_type_check(type = c("v", "j"), 
                                     alias = c("A", "B", "C")), 
               "v & j are not valid gate types for gate_draw!", 
               fixed = TRUE)
  expect_error(.cyto_gate_type_check(type = "z",
                                     alias = c("A", "B", "C")), 
               "z is not a valid gate type for gate_draw!", 
               fixed = TRUE)
})

## .cyto_alias_check -----------------------------------------------------------

test_that(".cyto_alias_check stops gating process if alias is incorrect", {
  expect_error(.cyto_alias_check(type = "web"), 
               "Supply names of populations to 'alias' for checking.", 
               fixed = TRUE)
  expect_error(.cyto_alias_check(type = "quadrant", 
                                 alias = "Cells"), 
          "Supply 4 population names to 'alias' to construct quadrant gates.", 
               fixed = TRUE)
  expect_error(.cyto_alias_check(type = c("rectangle", "polygon"),
                                 alias = "Cells"), 
           "Length of alias must be the same length as type for multi-gates.",
               fixed = TRUE)
})

## .cyto_gatingTemplate_check --------------------------------------------------

test_that(".cyto_gatingTemplate_check", {
  expect_error(
    .cyto_gatingTemplate_check(parent = "Cells", 
                              alias = "Single Cells", 
                          gatingTemplate = "Compensation-gatingTemplate.csv"), 
               "Supply another gatingTemplate or edit gate(s) using gate_edit.",
    fixed = TRUE)
})

## cyto_trans_check ------------------------------------------------------------

test_that("cyto_trans_check", {
  expect_error(cyto_trans_check(trans = "Test"), 
               "'trans' should be of class transformList or transformerList.")
})

## .cyto_overlay_check ---------------------------------------------------------

test_that(".cyto_overlay_check flowFrame method returns list of flowFrames", {
  expect_equal(.cyto_overlay_check(fs[[1]], 
                                   overlay = fs[[2]]), 
               list(fs[[2]]))
  expect_equal(.cyto_overlay_check(fs[[1]],
                                   overlay = list(fs[[2]])), 
               list(fs[[2]]))
  expect_equal(.cyto_overlay_check(fs[[1]], 
                                   overlay = fs), 
               list(fs[[1]], fs[[2]], fs[[3]], fs[[4]]))
  expect_equal(.cyto_overlay_check(fs[[1]], 
                                   overlay = list(fs)), 
               list(fs[[1]], fs[[2]], fs[[3]], fs[[4]]))

  expect_error(.cyto_overlay_check(fs[[1]], 
                                   overlay = list(fs, fs)), 
               paste("'overlay' must be a flowFrame, flowSet,",
                     "list of flowFrames or a list containing a flowSet.",
                     sep = " "))

})

test_that(".cyto_overlay_check flowSet method", {
  expect_equal(.cyto_overlay_check(fs, 
                                   overlay = fs[[2]]), 
              list(list(fs[[2]]), list(fs[[2]]), list(fs[[2]]), list(fs[[2]])))
  expect_equal(.cyto_overlay_check(fs, 
                                   overlay = list(fs[[2]])), 
              list(list(fs[[2]]), list(fs[[2]]), list(fs[[2]]), list(fs[[2]])))
  expect_equal(.cyto_overlay_check(fs, 
                                   overlay = fs), 
              list(list(fs[[1]]), list(fs[[2]]), list(fs[[3]]), list(fs[[4]])))
  expect_equal(.cyto_overlay_check(fs, 
                                   overlay = list(fs)), 
              list(list(fs[[1]]), list(fs[[2]]), list(fs[[3]]), list(fs[[4]])))
  expect_equal(.cyto_overlay_check(fs, 
                                   overlay = list(fs, fs)), 
               list(list(fs[[1]], fs[[1]]), 
                    list(fs[[2]], fs[[2]]), 
                    list(fs[[3]], fs[[3]]), 
                    list(fs[[4]], fs[[4]])))
  expect_equal(.cyto_overlay_check(fs, 
                                   overlay = list(fs[[1]], 
                                                  fs[[2]], 
                                                  fs[[3]], 
                                                  fs[[4]])), 
              list(list(fs[[1]]), list(fs[[2]]), list(fs[[3]]), list(fs[[4]])))

  expect_error(.cyto_overlay_check(fs, 
                                   overlay = exprs(fs[[1]])), 
               paste("'overlay' must be a flowFrame, flowSet,",
                     "list of flowFrames or a list of flowSets.",
                     sep = " "))
  
  expect_error(.cyto_overlay_check(fs, 
                                   overlay = list(fs[1:3])), 
               paste("Each flowSet in supplied list should be of the", 
                     "same length as the supplied flowSet.", sep = " "))
  expect_error(.cyto_overlay_check(fs, 
                                   overlay = list(fs[[1]], fs[[2]])),
               paste("Supplied list of flowFrames must be of the", 
                      "same length as the flowSet.", sep = " "))
  expect_error(.cyto_overlay_check(fs, 
                                   overlay = list(list(fs[[1]]), 
                                                  list(fs[[2]]))), 
            paste("'overlay' should be a list of flowFrames lists to overlay", 
                    "on each flowFrame in the flowSet.", sep = " "))
})

test_that(".cyto_overlay_check GatingHierarchy method", {
  expect_error(.cyto_overlay_check(gs[[1]], 
                                   overlay = "TEST"), 
               "'overlay' does not exist in the GatingHierarchy.", 
               fixed = TRUE)

  ov <- list(getData(gs, "T Cells")[[1]])
  names(ov) <- "T Cells"

  expect_equal(.cyto_overlay_check(gs[[1]], 
                                   overlay = "T Cells"), 
               ov)
  expect_equal(.cyto_overlay_check(gs[[1]], 
                                   overlay = getData(gs, "T Cells")[[4]]), 
               list(getData(gs, "T Cells")[[4]]))
})

test_that(".cyto_overlay_check GatingSet method", {
  expect_error(.cyto_overlay_check(gs, 
                                   overlay = "TEST"), 
               "overlay' does not exist in the GatingHierarchy.", 
               fixed = TRUE)

  TC <- getData(gs, "T Cells")
  ov <- list(list(TC[[1]]), list(TC[[2]]), list(TC[[3]]), list(TC[[4]]))

  expect_equivalent(.cyto_overlay_check(gs, overlay = "T Cells"), ov)
})

## .cyto_stat_check ------------------------------------------------------------

test_that(".cyto_stat_check", {
  expect_error(.cyto_stat_check("average"), 
               "Supplied statistic not supported.")
  expect_equal(.cyto_stat_check("Mean"), "mean")
  expect_equal(.cyto_stat_check("Median"), "median")
  expect_equal(.cyto_stat_check("Mode"), "mode")
  expect_equal(.cyto_stat_check("Count"), "count")
  expect_equal(.cyto_stat_check("Freq"), "freq")
  expect_equal(.cyto_stat_check("Geo mean"), "geo mean")
  expect_equal(.cyto_stat_check("cv"), "CV")
})

graphics.off()
