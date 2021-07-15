## CYTO_PLOT_PROFILE -----------------------------------------------------------

#' @export
cyto_plot_profile <- function(x,
                              channels = NULL,
                              layout = NULL,
                              header = NULL,
                              ...) {
  
  # CYTO_PLOT_EXIT -------------------------------------------------------------
  
  # SIGNAL CALL TO CYTO_PLOT
  if(is.null(cyto_option("cyto_plot_method"))) {
    cyto_option("cyto_plot_method", "profile")
    cyto_option("cyto_plot_data", TRUE)
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_par(reset = TRUE)
      cyto_option("cyto_plot_method", NULL)
      cyto_option("cyto_plot_data", FALSE)
    })
  }
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- names(cyto_markers(x))
  }
  
  # HEADER
  if(is.null(header)) {
    if(.all_na(overlay)) {
      header <- "Expression Profile"
    } else if(n > 1 & !.all_na(overlay)) {
      header <- paste(channels, "Expression")
    } else {
      header <- "Expression Profile"
    }
  }
  
  # REPEAT HEADER ARGUMENTS
  header <- rep(header, length.out = length(channels))
  
  
  
}