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
  options("CytoRSuite_gatingTemplate" = NULL)
  
  # Bypass working directory checks for external files
  options("CytoRSuite_wd_check" = TRUE)
  
  # Signal when gate_draw has been called - turn off overlay sampling
  options("CytoRSuite_gate_draw" = FALSE)
  
  # Create custom theme for cyto_plot
  options("CytoRSuite_cyto_plot_theme" = NULL)
  
  # Signal cyto_plot_save method has been called
  options("CytoRSuite_cyto_plot_save" = FALSE)
  
  # Signal which cyto_plot method has been called
  options("CytoRSuite_cyto_plot_method" = NULL)
  
  # Signal if a custom plot is being contructed - require cyto_plot_complete
  options("CytoRSuite_cyto_plot_custom" = FALSE)
  
  # Signal when cyto_plot_grid method is being called
  options("CytoRSuite_cyto_plot_grid" = FALSE)
  
  # Signal when axes text has already been calculated
  options("CytoRSuite_cyto_plot_axes_text" = NULL)
  
  # Register gating and preprocessing functions with openCyto
  openCyto::registerPlugins(fun = .gate_manual, 
                            methodName = "gate_manual")
  openCyto::registerPlugins(fun = .gate_draw, 
                            methodName = "gate_draw")
  openCyto::registerPlugins(fun = .pp_gate_draw, 
                            methodName = "pp_gate_draw", 
                            dep = NA, "preprocessing")
  
}