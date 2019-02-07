context("cyto_plot")

# vdiffr::manage_cases() to validate and generate reference images

# cyto_plot flowFrame method ---------------------------------------------------

test_that("cyto_plot flowFrame method", {
  
  # 1-D density distribution errors
  expect_error(cyto_plot(Va2[[1]]), 
               "Supply channel/marker(s) to construct the plot.", 
               fixed = TRUE)
  
  # 1-D density distributions
  p <- function() cyto_plot(Va2[[1]], 
                             channels = "CD69", 
                             axes_trans = trans,
                             overlay = Va2[2:3],
                             gate = getGate(gs, "CD69+ CD4 T Cells")[[1]],
                             title = "Control",
                             xlab = "CD69",
                             ylab = "Density",
                             density_modal = FALSE,
                             density_smooth = 1.6,
                             density_stack = 0.5,
                             density_fill_alpha = c(0,0.5,1),
                             density_line_type = c(1,2,3),
                             density_line_width = 1.5,
                             density_line_col = "blue",
                             axes_text_font = 3,
                             axes_text_size = 1.5,
                             axes_text_col = "green4",
                             axes_label_text_font = 2,
                             axes_label_text_size = 1.5,
                             axes_label_text_col = "grey50",
                             title_text_font = 4,
                             title_text_size = 3,
                             title_text_col = "red",
                             border_line_type = 3,
                             border_line_width = 2.5,
                             border_line_col = "purple", 
                             legend = TRUE,
                             label_text = c("A","B","C"),
                             label_text_col = c("blue","black","green"))
  expect_doppelganger("cyto_plot-fr1", p)
  
})

# cyto_plot flowSet method -----------------------------------------------------

test_that("cyto_plot flowSet method", {
  
  # Channels error
  expect_error(cyto_plot(Va2), 
               "Supply channel/marker(s) to construct the plot.",
               fixed = TRUE)
  
  
})

# cyto_plot GatingHierarchy method ---------------------------------------------

test_that("cyto_plot GatingHierarchy method", {
  
  # Channels error
  expect_error(cyto_plot(gs[[1]], parent = "T Cells"), 
               "Supply channel/marker(s) to construct the plot.", 
               fixed = TRUE)
  
  # Parent error
  expect_error(cyto_plot(gs[[1]], 
                         channels = c("CD4","CD8")), 
               "Please supply the name of the parent population to plot.")
  
  # 2-D Scatter plots
  p <- function() cyto_plot(gs[[4]], 
                            parent = "T Cells",
                            alias = c("CD4 T Cells","CD8 T Cells"),
                            channels = c("CD4","CD8"), 
                            overlay = c("CD4 T Cells","CD8 T Cells"),
                            xlab = "CD4",
                            ylab = "CD8",
                            axes_text_font = 3,
                            axes_text_size = 1.5,
                            axes_text_col = "green4",
                            axes_label_text_font = 2,
                            axes_label_text_size = 1.5,
                            axes_label_text_col = c("red"),
                            title_text_font = 4,
                            title_text_size = 3,
                            title_text_col = "blue",
                            border_line_type = 3,
                            border_line_width = 2.5,
                            border_line_col = "purple", 
                            legend = TRUE,
                            label_text_col = c("blue","black"))
  expect_doppelganger("cyto_plot-gh1", p)
  
})

# cyto_plot GatingSet method ---------------------------------------------------

test_that("cyto_plot GatingSet method", {
  
  # Channels error
  expect_error(cyto_plot(gs, parent = "T Cells"), 
               "Supply channel/marker(s) to construct the plot.", 
               fixed = TRUE)
  
  # Parent error
  expect_error(cyto_plot(gs, 
                         channels = c("CD4","CD8")), 
               "Please supply the name of the parent population to plot.")
  
  # 1-D density distributions
  p <- function() cyto_plot(gs, 
                            parent = "T Cells",
                            channels = "CD8", 
                            axes_trans = trans,
                            overlay = c("CD4 T Cells","CD8 T Cells"),
                            gate = getGate(gs, "CD8 T Cells")[[1]],
                            xlab = "CD8 Expression",
                            ylab = "Density",
                            density_smooth = 1.6,
                            density_stack = 0.5,
                            density_fill = c("red",
                                             "blue",
                                             "green",
                                             "orange",
                                             "purple",
                                             "black"),
                            density_fill_alpha = c(0.2,0.5,1),
                            density_line_type = c(1,2,3),
                            density_line_width = 1.5,
                            density_line_col = "blue",
                            axes_text_font = 3,
                            axes_text_size = 1.5,
                            axes_text_col = "green4",
                            axes_label_text_font = 2,
                            axes_label_text_size = 1.5,
                            axes_label_text_col = c("red",
                                                    "yellow",
                                                    "darkorange",
                                                    "grey40"),
                            title_text_font = 4,
                            title_text_size = 3,
                            title_text_col = "red",
                            border_line_type = 3,
                            border_line_width = 2.5,
                            border_line_col = "purple", 
                            legend = TRUE,
                            label_text = c("A","B","C"),
                            label_text_col = c("blue",
                                               "black",
                                               "green"))
  expect_doppelganger("cyto_plot-gs1", p)
  
  # 1-D density distributions & group_by
  p <- function() cyto_plot(gs, 
                            parent = "T Cells",
                            alias = "CD8 T Cells",
                            group_by = "Treatment",
                            channels = "CD8", 
                            axes_trans = trans,
                            overlay = c("CD4 T Cells","CD8 T Cells"),
                            xlab = "CD8 Expression",
                            ylab = "Density",
                            density_smooth = 1.6,
                            density_stack = 0.5,
                            density_fill = c("red","blue","green"),
                            density_fill_alpha = c(0.2,0.5,1),
                            density_line_type = c(1,2,3),
                            density_line_width = 1.5,
                            density_line_col = "blue",
                            axes_text_font = 3,
                            axes_text_size = 1.5,
                            axes_text_col = "green4",
                            axes_label_text_font = 2,
                            axes_label_text_size = 1.5,
                            axes_label_text_col = c("red","darkorange"),
                            title_text_font = 4,
                            title_text_size = 3,
                            title_text_col = "red",
                            border_line_type = 3,
                            border_line_width = 2.5,
                            border_line_col = "purple", 
                            legend = TRUE)
  expect_doppelganger("cyto_plot-gs2", p)
  
  # 2-D Scatter plots
  p <- function() cyto_plot(gs, 
                            parent = "T Cells",
                            alias = c("CD4 T Cells","CD8 T Cells"),
                            channels = c("CD4","CD8"), 
                            overlay = c("CD4 T Cells","CD8 T Cells"),
                            xlab = "CD8 Expression",
                            ylab = "Density",
                            title = c("Title1", "Title2","Title3","Title4"),
                            axes_text_font = 3,
                            axes_text_size = 1.5,
                            axes_text_col = "green4",
                            axes_label_text_font = 2,
                            axes_label_text_size = 1.5,
                            axes_label_text_col = c("red",
                                                    "yellow",
                                                    "darkorange",
                                                    "grey40"),
                            title_text_font = 4,
                            title_text_size = 3,
                            title_text_col = "red",
                            border_line_type = 3,
                            border_line_width = 2.5,
                            border_line_col = "purple", 
                            legend = TRUE,
                            label_text = c("A","B","C"),
                            label_text_col = c("blue","black","green"))
  expect_doppelganger("cyto_plot-gs3", p)
  
})