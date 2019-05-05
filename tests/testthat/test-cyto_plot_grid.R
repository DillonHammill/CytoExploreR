context("cyto_plot_grid")

cyto_plot_save("GRID1.png")
cyto_plot_grid(getData(gs,"Live Cells"),
               channels = c("CD4","CD8"),
               axes_trans = trans,
               title = c("TESTX","TESTY"),
               panel_y_label = c(1,2),
               panel_x_label = c(3,4))

cyto_plot_save("GRID2.png")
cyto_plot_grid(getData(gs,"Live Cells"),
               channels = c("FSC-A","SSC-A"),
               title = c("TESTX","TESTY"),
               panel_y_label = c(1,2),
               panel_x_label = c(3,4))

cyto_plot_grid(getData(gs,"Live Cells"),
               channels = c("CD4","CD8"),
               group_by = c("Treatment","OVAConc"),
               axes_trans = trans,
               gate = list(getGate(gs,"CD4 T Cells")[[1]], 
                           getGate(gs,"CD8 T Cells")[[1]]),
               label_text = c("CD4 T Cells","CD8 T Cells"),
               point_size = 3,
               title_text_col = c("blue"),
               title = "OVA Concentration",
               overlay = getData(gs,"CD69+ CD4 T Cells"))

cyto_plot_grid(getData(gs,"Live Cells"),
               channels = c("CD4","CD8"),
               group_by = "OVAConc",
               axes_trans = trans,
               gate = list(getGate(gs,"CD4 T Cells")[[1]], 
                           getGate(gs,"CD8 T Cells")[[1]]),
               label_text = c("CD4 T Cells","CD8 T Cells"),
               point_size = 3,
               title_text_col = c("blue"),
               title = "OVA Concentration",
               overlay = getData(gs,"CD69+ CD4 T Cells"))