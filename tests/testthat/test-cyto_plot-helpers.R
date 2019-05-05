context("cyto_plot-helpers")

# These tests can probably be removed later when higher level tests are added

# cyto_plot_empty --------------------------------------------------------------

test_that("cyto_plot_empty", {
  
  cyto_plot_empty(Va2[[1]],
                  channels = "CD4",
                  axes_trans = trans,
                  title = "Testing",
                  title_text_col = "red",
                  border_fill = "grey",
                  xlab = "FSC-A",
                  ylab = "Density")
  
})
