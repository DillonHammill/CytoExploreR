#' Set gloabl options and register plugins with openCyto on package loading
#'
#' Global options are set to NULL by default as these will not be provided to
#' mapply.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  
  # Turn off bell noises whilst gating
  options("locatorBell" = FALSE)
  
  # Interactive mode
  options("CytoExploreR_interactive" = TRUE)
  
  # Select active gatingTemplate
  options("CytoExploreR_gatingTemplate" = NULL)
  
  # Create custom theme for cyto_plot
  options("cyto_plot_theme" = NULL)
  
  # Signal cyto_plot_save method has been called
  options("cyto_plot_save" = FALSE)
  
  # Signal which cyto_plot method has been called
  options("cyto_plot_method" = NULL)
  
  # Signal if a custom plot is being constructed - require cyto_plot_complete
  options("cyto_plot_custom" = FALSE)
  
  # Save layout settings
  options("cyto_plot_layout" = NULL)
  
  # Register gating and preprocessing functions with openCyto
  suppressMessages({
    openCyto::register_plugins(fun = .cyto_gate_manual, 
                               methodName = "cyto_gate_manual")
    openCyto::register_plugins(fun = .cyto_gate_draw, 
                               methodName = "cyto_gate_draw")
    openCyto::register_plugins(fun = .pp_cyto_gate_draw, 
                               methodName = "pp_cyto_gate_draw", 
                               dep = NA, "preprocessing")
  })
  
}

# # Docker Message
# .onAttach <- function(libname, pkgname){
#   packageStartupMessage("CytoExploreR Docker v1.0.5")
# }