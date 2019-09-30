#' Set gloabl options and register plugins with openCyto on package loading
#'
#' Global options are set to NULL by default as these will not be provided to
#' mapply.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  
  # Turn off bell noises whilst gating
  options("locatorBell" = FALSE)
  
  # Select active gatingTemplate
  options("CytoExploreR_gatingTemplate" = NULL)
  
  # Bypass working directory checks for external files
  options("CytoExploreR_wd_check" = TRUE)
  
  # Signal when gate_draw has been called - turn off overlay sampling
  options("cyto_gate_draw" = FALSE)
  
  # Signals args called to cyto_plot - check if call is made twice
  options("cyto_plot_call" = NULL)
  
  # Create custom theme for cyto_plot
  options("cyto_plot_theme" = NULL)
  
  # Signal cyto_plot_save method has been called
  options("cyto_plot_save" = FALSE)
  
  # Signal which cyto_plot method has been called
  options("cyto_plot_method" = NULL)
  
  # Signal if a custom plot is being contructed - require cyto_plot_complete
  options("cyto_plot_custom" = FALSE)
  
  # Signal when cyto_plot_grid method is being called
  options("cyto_plot_grid" = FALSE)
  
  # Signal previous call to cyto_plot (same plot?)
  options("cyto_plot_call" = NULL)
  
  # Save label co-ordinates as list
  options("cyto_plot_label_coords" = NULL)
  
  # Register gating and preprocessing functions with openCyto
  openCyto::registerPlugins(fun = .gate_manual, 
                            methodName = "gate_manual")
  openCyto::registerPlugins(fun = .cyto_gate_draw, 
                            methodName = "cyto_gate_draw")
  openCyto::registerPlugins(fun = .pp_cyto_gate_draw, 
                            methodName = "pp_cyto_gate_draw", 
                            dep = NA, "preprocessing")
  
}