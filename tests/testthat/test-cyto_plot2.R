context("cyto_plot")

# GatingSet Method -------------------------------------------------------------

test_that("cyto_plot GatingSet method",{
  
  # Errors
  expect_error(cyto_plot(gs,
               channels = c("FSC-A","SSC-A")),
               "Supply the name of the 'parent' population to plot.", 
               fixed  = TRUE)
  expect_error(cyto_plot(gs,
               parent = "root"),
               "Supply channel/marker(s) to construct the plot.", 
               fixed = TRUE)
  expect_error(cyto_plot(gs,
                         parent = "root",
                         channels = "CD4",
                         density_layers = 8),
               "Each plot must have the same number of layers!")
  
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
            label_text = paste0("layer",1:8),
            label_stat = "count",
            label_box_x = 1.5)
  
  expect_doppelganger("cyto_plot-GatingSet-1D_001", p)
  
  # 1D grouped, stacked and gated
  p <- function() cyto_plot(gs[1:32],
            parent = "T Cells",
            alias = "",
            channels = "CD4",
            group_by = c("Treatment","OVAConc"),
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
            gate_line_type = 2)
  
  expect_doppelganger("cyto_plot-GatingSet-1D_002", p)
  
  # 1D grouped, overlay and gated
  p <- function() cyto_plot(gs[1:32],
            parent = "T Cells",
            alias = "",
            channels = "CD4",
            group_by = c("Treatment","OVAConc"),
            overlay = c("CD4 T Cells","CD8 T Cells"),
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
            density_line_col = "white")
  
  expect_doppelganger("cyto_plot-GatingSet-1D_003", p)
  
  # Overlay with empty data
  
  p <- function() cyto_plot(gs[c(31,33,30,32)],
                            parent = "CD4 T Cells",
                            alias = "",
                            channels = "CD69",
                            overlay = "CD8 T Cells")
  
  expect_doppelganger("cyto_plot-GatingSet-1D_004", p)
  
  # 2D Scatter plots -----------------------------------------------------------
  
  cyto_plot(gs[1:4],
            parent = "T Cells",
            alias = "",
            channels = c("CD4","CD8"),
            overlay = c("CD4 T Cells","Dendritic Cells"),
            label = FALSE,
            point_size = 3,
            legend = TRUE,
            border_fill = "black",
            contour_lines = 15,
            contour_line_col = "orange")
  
  cyto_plot(gs[1:4],
            parent = "T Cells",
            alias = "",
            channels = c("CD4","CD8"),
            overlay = rep(list(list(getData(gs, "Live Cells")[[1]],
                           getData(gs, "CD4 T Cells")[[1]])), 4),
            label = FALSE,
            point_size = 3,
            legend = TRUE)
  
  cyto_plot(gs[1], 
            parent = "T Cells", 
            channels = c("CD4","CD8"), 
            alias = "", 
            border_fill = "black", 
            point_size = 3, 
            point_col_scale = c("white","yellow","orange","red","brown"))
  
})
