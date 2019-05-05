context("cyto_plot internal helpers")

# .cyto_plot_axes_limits -------------------------------------------------------

test_that(".cyto_plot_axes_limits", {
  
  # flowFrame ------------------------------------------------------------------
  # machine limits -
  lms <- list(c(0,262144), c(-0.4696335, 4.5))
  names(lms) <- c("FSC-A","PE-A")
  expect_equal(.cyto_plot_axes_limits(getData(gs, "root")[[32]],
                                      channels = c("FSC-A","PE-A"),
                                      limits = "machine"), 
               lms, tolerance = 0.0001)
  
  # data limits -
  lms <- list(c(0, 262143), c(-0.3862161, 3.8621611))
  names(lms) <- c("FSC-A","PE-A")
  expect_equal(.cyto_plot_axes_limits(getData(gs, "root")[[1]],
                                      channels = c("FSC-A","PE-A"),
                                      limits = "data"), 
  lms, tolerance = 0.0001)
  
  # overlay -
  lms <- list(c(-0.3044766, 3.0447665), c(-0.1917132, 1.9171320))
  names(lms) <- c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
  expect_equal(.cyto_plot_axes_limits(getData(gs,"CD4 T Cells")[[1]],
                                      channels = c("CD4","CD8"),
                                      overlay = getData(gs, "CD8 T Cells")[[1]],
                                      limits = "data"),
               lms, tolerance = 0.0001)
  
  # flowSet --------------------------------------------------------------------
  # machine limits -
  lms <- list(c(0,262144), c(-0.4653093, 4.5))
  names(lms) <- c("FSC-A","7-AAD-A")
  expect_equal(.cyto_plot_axes_limits(getData(gs, "root"),
                                      channels = c("FSC-A","CD69"),
                                      limits = "machine"), 
               lms, tolerance = 0.0001)
  
  # data limits -
  lms <- list(c(0,262143), c(-0.4651663, 4.49857))
  names(lms) <- c("FSC-A","7-AAD-A")
  expect_equal(.cyto_plot_axes_limits(getData(gs, "root"),
                                      channels = c("FSC-A","CD69"),
                                      limits = "data"), 
               lms, tolerance = 0.0001)
  
  # overlay -
  lms <- list(c(-0.3141541, 3.141541), c(-0.2255783, 2.2557833))
  names(lms) <- c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
  expect_equal(.cyto_plot_axes_limits(getData(gs,"CD4 T Cells"),
                                      channels = c("CD4","CD8"),
                                      overlay = getData(gs, "CD8 T Cells"),
                                      limits = "data"),
               lms, tolerance = 0.0001)
  
  # GatingSet ------------------------------------------------------------------
  lms <- list(c(-0.3141541, 3.141541), c(-0.2255783, 2.2557833))
  names(lms) <- c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
  expect_equal(.cyto_plot_axes_limits(gs,
                                      parent = "CD4 T Cells",
                                      channels = c("CD4","CD8"),
                                      overlay = "CD8 T Cells",
                                      limits = "data"),
               lms, tolerance = 0.0001)
  
})

# .cyto_plot_axes_text ---------------------------------------------------------

test_that(".cyto_plot_axes_text", {
  
  # flowFrame
  exp <- .cyto_plot_axes_text(fr_test,
                              channels = c("FSC-A","PE-A"),
                              trans = trans)
  ref <- list(label = expression(0, 10^2, 10^3, 10^4, 10^5),
              at = c(0.787854, 0.985926, 1.955290, 3.065260, 4.0803))
  ref <- list(NULL, ref)
  names(ref) <- c("FSC-A","PE-A")
  expect_equal(exp, ref, tolerance = 0.001)
  
  # GatingHierarchy
  exp <- .cyto_plot_axes_text(gs[[32]],
                              channels = c("FSC-A","PE-A"))
  ref <- list(label = expression(0, 10^2, 10^3, 10^4, 10^5),
              at = c(0.787854, 0.985926, 1.955290, 3.065260, 4.0803))
  ref <- list(NULL, ref)
  names(ref) <- c("FSC-A","PE-A")
  expect_equal(exp, ref, tolerance = 0.001)
  
})

# .cyto_plot_gate_count --------------------------------------------------------

test_that(".cyto_plot_gate_count", {
  
  # NULL
  expect_equal(.cyto_plot_gate_count(NULL), 0)
  
  # Gate object
  expect_equal(.cyto_plot_gate_count(rg), 1)
  expect_equal(.cyto_plot_gate_count(pg), 1)
  expect_equal(.cyto_plot_gate_count(eg), 1)
  
  # filters object
  expect_equal(.cyto_plot_gate_count(filters(getGate(gs,"CD4 T Cells"))), 33)
  
  # list of gate objects
  expect_equal(.cyto_plot_gate_count(list(rg,eg,pg)), 3)
  
  # list of gate object lists
  expect_equal(.cyto_plot_gate_count(list(getGate(gs,"CD4 T Cells"),
                                          getGate(gs, "CD8 T Cells"))), 66)
  
  # list gate objects and filters object
  expect_equal(.cyto_plot_gate_count(
    list(filters(getGate(gs,"CD4 T Cells")),
                 getGate(gs,"CD4 T Cells")[[1]])), 34)
  
})

# .cyto_plot_args --------------------------------------------------------------

test_that(".cyto_plot_args", {
  
  # Arguments for .cyto_plot_1d
  expect_equal(.cyto_plot_args("FSC-A"),
               c("xlab",
                 "ylab",
                 "title",
                 "title_text_font",
                 "title_text_size",
                 "title_text_col",
                 "axes_text_font",
                 "axes_text_size",
                 "axes_text_col",
                 "axes_label_text_font",
                 "axes_label_text_size",
                 "axes_label_text_col",
                 "border_line_type",
                 "border_line_width",
                 "border_line_col",
                 "legend",
                 "legend_text",
                 "legend_text_font",
                 "legend_text_size",            
                 "legend_text_col",
                 "label",            
                 "label_text",
                 "label_stat",
                 "label_text_font",
                 "label_text_size",
                 "label_text_col",
                 "label_box_x",
                 "label_box_y",
                 "label_box_alpha",
                 "gate_line_type", 
                 "gate_line_width", 
                 "gate_line_col",
                 "density_stack",
                 "density_fill",
                 "density_fill_alpha",
                 "density_line_type",
                 "density_line_width",
                 "density_line_col",            
                 "legend_line_col",
                 "legend_box_fill"))  
  
  # Arguments for .cyto_plot_2d
  expect_equal(.cyto_plot_args(c("FSC-A","PE-A")),
               c("xlab",
                 "ylab",
                 "title",
                 "title_text_font",
                 "title_text_size",
                 "title_text_col",
                 "axes_text_font",
                 "axes_text_size",
                 "axes_text_col",
                 "axes_label_text_font",
                 "axes_label_text_size",
                 "axes_label_text_col",
                 "border_line_type",
                 "border_line_width",
                 "border_line_col",
                 "legend",
                 "legend_text",
                 "legend_text_font",
                 "legend_text_size",            
                 "legend_text_col",
                 "label",            
                 "label_text",
                 "label_stat",
                 "label_text_font",
                 "label_text_size",
                 "label_text_col",
                 "label_box_x",
                 "label_box_y",
                 "label_box_alpha",
                 "gate_line_type", 
                 "gate_line_width", 
                 "gate_line_col",
                 "point_shape",
                 "point_size",
                 "point_col",
                 "point_alpha",
                 "contour_lines",
                 "contour_line_type",
                 "contour_line_width",
                 "contour_line_col",            
                 "legend_point_col"))
  
})

# .cyto_plot_args_split --------------------------------------------------------

test_that(".cyto_plot_args_split", {
  
  # List of args 1 in each 'category'
  args <- list("invalid" = 5,
               "xlab" = "CD44",
               "point_shape" = 16,
               "gate_line_type" = 2,
               "label_text_size" = 1.5)
  
  ref <- list("xlab" = list(`1` = "CD44", 
                            `2` = "CD44", 
                            `3` = "CD44"),
              "gate_line_type" = list(`1` = c(2,2), 
                                      `2` = c(2,2), 
                                      `3` = c(2,2)),
              "label_text_size" = list(`1` = rep(1.5,6), 
                                       `2` = rep(1.5,6), 
                                       `3` = rep(1.5,6)))

  # .cyto_plot_1d
  expect_equal(.cyto_plot_args_split(args,
                                     channels = "FSC-A",
                                     n = 9, # total layers
                                     plots = 3,
                                     layers = 3,
                                     gates = 2),
               ref)
  
  ref <- list("xlab" = list(`1` = "CD44",
                            `2` = "CD44",
                            `3` = "CD44"),
              "point_shape" = list(`1` = c(16,16,16),
                                   `2` = c(16,16,16),
                                   `3` = c(16,16,16)),
              "gate_line_type" = list(`1` = c(2,2),
                                      `2` = c(2,2),
                                      `3` = c(2,2)),
              "label_text_size" = list(`1` = c(1.5,1.5),
                                       `2` = c(1.5,1.5),
                                       `3` = c(1.5,1.5)))
  
  # .cyto_plot_2d
  expect_equal(.cyto_plot_args_split(args,
                                     channels = c("FSC-A","PE-A"),
                                     n = 9, # total layers
                                     plots = 3,
                                     layers = 3,
                                     gates = 2),
              ref)

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

# .cyto_plot_layout ------------------------------------------------------------

# .cyto_plot_margins -----------------------------------------------------------

# .cyto_plot_legend ------------------------------------------------------------

# .cyto_plot_title -------------------------------------------------------------

# .cyto_plot_axes_label --------------------------------------------------------

# .cyto_plot_density -----------------------------------------------------------
