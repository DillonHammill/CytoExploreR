context("cyto_plot")

# flowFrame Method -------------------------------------------------------------

test_that("cyto_plot flowFrame method", {
  
  expect_error(cyto_plot(fs[[1]],
                         parent = "root"
  ),
  "Supply channel/marker(s) to construct the plot.",
  fixed = TRUE)
  
  # 1D density distribution
  p <- function() cyto_plot(cyto_extract(gs,"T Cells")[[1]],
                            channels = "CD4",
                            overlay = list(cyto_extract(gs, "CD4 T Cells")[[1]],
                                           cyto_extract(gs, "CD8 T Cells")[[1]]))
  
  expect_doppelganger("cyto_plot-flowFrame-1D_001", p)
  
  # 2D scatter plot
  p <- function() cyto_plot(as(getData(gs, "T Cells"),"flowFrame"),
            channels = c("CD8","CD4"),
            gate = getGate(gs, "CD4 T Cells")[[1]],
            axes_trans = trans,
            legend = TRUE,
            contour_lines = 5,
            contour_line_col = "magenta",
            border_fill = "black",
            limits = "data")
  
  expect_doppelganger("cyto_plot-flowFrame-2D_001", p)
  
})

# GatingHierarchy Method -------------------------------------------------------

test_that("cyto_plot GatingHierarchy method", {
  
  # Errors
  expect_error(cyto_plot(gs[[1]],
                         channels = c("FSC-A", "SSC-A")
  ),
  "Supply the name of the 'parent' population to plot.",
  fixed = TRUE
  )
  expect_error(cyto_plot(gs[[1]],
                         parent = "root"
  ),
  "Supply channel/marker(s) to construct the plot.",
  fixed = TRUE
  )
  
  # 1D density distributions ---------------------------------------------------
  
  p <- function() cyto_plot(gs[[1]],
            parent = "T Cells",
            channels = "CD4",
            overlay = c("CD4 T Cells","CD8 T Cells"),
            gate = getGate(gs, "CD4 T Cells")[[1]],
            border_fill = "grey",
            legend = "line",
            density_line_type = c(1,2,3),
            density_line_width = 2,
            density_fill = "white",
            density_line_col = c("red","blue","green"),
            density_stack = 0.5,
            density_fill_alpha = c(0,0.5,1))
  
  expect_doppelganger("cyto_plot-GatingHierarchy-1D_001", p)
  
  # 2D scatter plots -----------------------------------------------------------
  
  # overlay by name
  p <- function() cyto_plot(gs[[1]],
            parent = "T Cells",
            alias = "",
            channels = c("CD44","CD69"),
            overlay = c("CD4 T Cells","CD8 T Cells"),
            point_size = 3,
            legend = TRUE,
            contour_lines = 15,
            border_fill = "black",
            point_col = c("red","orange","yellow"),
            gate_line_col = "magenta")
  
  expect_doppelganger("cyto_plot-GatingHierarchy-2D_001", p)
  
  # flowSet overlay
  p <- function() cyto_plot(gs[[1]],
            parent = "T Cells",
            alias = "",
            channels = c("CD4","CD8"),
            overlay = getData(gs, "CD4 T Cells")[1:2],
            point_size = 3,
            point_col_scale = c("green","yellow","magenta","purple"),
            legend = TRUE,
            axes_text = c(FALSE,FALSE),
            border_fill = "black")
  
  expect_doppelganger("cyto_plot-GatingHierarchy-2D_002", p)
  
  # flowFrame overlay
  p <- function() cyto_plot(gs[[1]],
            parent = "T Cells",
            alias = "",
            channels = c("CD4","CD8"),
            overlay = getData(gs, "CD4 T Cells")[[1]],
            point_size = 3)
  
  expect_doppelganger("cyto_plot-GatingHierarchy-2D_003", p)
  
})

# GatingSet Method -------------------------------------------------------------

test_that("cyto_plot GatingSet method", {

  # Errors
  expect_error(cyto_plot(gs,
    channels = c("FSC-A", "SSC-A")
  ),
  "Supply the name of the 'parent' population to plot.",
  fixed = TRUE
  )
  expect_error(cyto_plot(gs,
    parent = "root"
  ),
  "Supply channel/marker(s) to construct the plot.",
  fixed = TRUE
  )
  expect_error(
    cyto_plot(gs,
      parent = "root",
      channels = "CD4",
      density_layers = 8
    ),
    "Each plot must have the same number of layers!"
  )

  # 1D Density distributions ---------------------------------------------------

  # 1D stacked samples - labels without gates + stats
  p <- function() cyto_plot(gs[1:32],
      parent = "T Cells",
      channels = "CD4",
      density_smooth = 0.2,
      density_layers = 8,
      density_stack = 0.8,
      legend = TRUE,
      label = TRUE,
      label_text = paste0("layer", 1:8),
      label_stat = "count",
      label_box_x = 1.5
    )

  expect_doppelganger("cyto_plot-GatingSet-1D_001", p)

  # 1D grouped, stacked and gated
  p <- function() cyto_plot(gs[1:32],
      parent = "T Cells",
      alias = "",
      channels = "CD4",
      group_by = c("Treatment", "OVAConc"),
      density_smooth = 0.2,
      legend = TRUE,
      label_text_font = 1,
      label_text_size = 0.6,
      label_text_col = "red",
      label_box_alpha = 1,
      density_modal = FALSE,
      density_layers = 8,
      density_fill_alpha = 0.5,
      density_line_type = 2,
      density_line_width = 2.5,
      density_line_col = "blue",
      gate_line_col = "grey",
      gate_line_width = 3,
      gate_line_type = 2
    )

  expect_doppelganger("cyto_plot-GatingSet-1D_002", p)

  # 1D grouped, overlay and gated
  p <- function() cyto_plot(gs[1:32],
      parent = "T Cells",
      alias = "",
      channels = "CD4",
      group_by = c("Treatment", "OVAConc"),
      overlay = c("CD4 T Cells", "CD8 T Cells"),
      title_text_font = 1,
      title_text_size = 1,
      title_text_col = "deepskyblue",
      axes_text_font = 2,
      axes_text_size = 0.8,
      axes_text_col = "magenta",
      axes_label_text_font = 2,
      axes_label_text_size = 0.8,
      axes_label_text_col = "green",
      label_text_size = 0.4,
      border_line_type = 3,
      border_line_width = 2,
      border_line_col = "orange",
      border_fill = "black",
      border_fill_alpha = 0.2,
      density_line_col = "white"
    )

  expect_doppelganger("cyto_plot-GatingSet-1D_003", p)

  # Overlay with empty data

  p <- function() cyto_plot(gs[c(31, 33, 30, 32)],
      parent = "CD4 T Cells",
      alias = "",
      channels = "CD69",
      overlay = "CD8 T Cells"
    )

  expect_doppelganger("cyto_plot-GatingSet-1D_004", p)

  # 2D Scatter plots -----------------------------------------------------------

  # Multiple plots - check most arguments
  p <- function() cyto_plot(gs[1:4],
      parent = "T Cells",
      alias = "",
      channels = c("CD4", "CD8"),
      overlay = c("CD4 T Cells", "Dendritic Cells"),
      title = paste0("Test", 1:4),
      title_text_size = 0.8,
      title_text_font = 1,
      title_text_col = c("red", "brown", "orange", "green"),
      axes_text = c(FALSE, TRUE),
      axes_text_font = 3,
      axes_text_size = 0.8,
      axes_text_col = "deepskyblue",
      axes_label_text_font = 2,
      axes_label_text_size = 1.5,
      axes_label_text_col = "purple",
      gate_line_type = c(1, 2),
      gate_line_width = 3,
      gate_line_col = c("cyan", "green"),
      legend = TRUE,
      legend_text_col = c("green", "orange", "blue"),
      legend_text_size = 2,
      legend_text_font = 2,
      point_size = c(2, 3, 4),
      label_box_alpha = 0.2,
      label_text_size = 0.8,
      label_text_col = c("red", "blue"),
      label_text_font = 1,
      border_fill = "grey",
      border_line_col = "red",
      border_line_type = 2,
      border_line_width = 2
    )

  expect_doppelganger("cyto_plot-GatingSet-2D_001", p)

  # Grouped samples with contours
  p <- function() cyto_plot(gs[1:32],
      group_by = c("Treatment", "OVAConc"),
      parent = "T Cells",
      alias = "",
      channels = c("CD4", "CD8"),
      overlay = c("CD4 T Cells", "CD8 T Cells"),
      point_col = c("black", "red", "green4"),
      contour_lines = 15,
      contour_line_col = "blue"
    )

  expect_doppelganger("cyto_plot-GatingSet-2D_002", p)

  # Multiple with zero event sample
  p <- function() cyto_plot(gs[c(1, 33, 2:3, 4:8)],
      parent = "T Cells",
      channels = c("CD4", "CD8"),
      alias = "",
      overlay = c("CD4 T Cells", "CD8 T Cells"),
      point_size = 3,
      label_text_size = 0.6,
      contour_lines = 15,
      contour_line_col = "magenta",
      border_fill = "black",
      border_fill_alpha = 0.2
    )

  expect_doppelganger("cyto_plot-GatingSet-2D_003", p)
  
  # Supply gate objects through gate argument
  
  p <- function() cyto_plot(gs[1:4],
            parent = "T Cells",
            channels = c("CD4","CD8"),
            gate = list(getGate(gs, "CD4 T Cells")[[1]],
                        getGate(gs, "CD8 T Cells")[[1]]))

  expect_doppelganger("cyto_plot-GatingSet-2D_004", p)
  
})
