context("cyto_plot-internal update")

png("TestA.png",
    height = 7,
    width = 7,
    units = "in",
    res = 300)
.cyto_plot_1d_v2(Va2[[1]], 
                 channel = c("PE-A"),
                 overlay = Va2[2:32],
                 axes_trans = trans,
                 density_smooth = 0.5,
                 limits = "machine",
                 density_modal = FALSE,
                 legend = TRUE)
dev.off()

png("TestB.png",
    height = 7,
    width = 7,
    units = "in",
    res = 300)
.cyto_plot_1d_v2(Va2[[1]], 
                 channel = c("PE-A"),
                 overlay = Va2[2:32],
                 axes_trans = trans,
                 density_smooth = 0.5,
                 limits = "data",
                 density_modal = FALSE,
                 legend = TRUE)
dev.off()