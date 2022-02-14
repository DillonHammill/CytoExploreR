#' Set gloabl options and register plugins with openCyto on package loading
#'
#' Global options are set to NULL by default as these will not be provided to
#' mapply.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  
  # GATING BELL NOISES
  options("locatorBell" = FALSE)
  
  # CYTOEXPLORER PARALLEL CONFIGURATION
  # options("CytoExploreR_parallel" = NULL)
  
  # INTERACTIVE MODE
  options("CytoExploreR_interactive" = TRUE)
  
  # ACTIVE GATINGTEMPLATE
  options("CytoExploreR_gatingTemplate" = NULL)
  
  # CUSTOM THEME FOR CYTO_PLOT()
  options("cyto_plot_theme" = NULL)
  
  # CYTO_PLOT_SAVE()
  options("cyto_plot_save" = FALSE)
  
  # CYTO_PLOT() METHOD
  options("cyto_plot_method" = NULL)
  
  # CYTO_PLOT() CUSTOM IN SHINY
  options("cyto_plot_shiny" = NULL)
  
  # CYTO_PLOT() SET GRAPHICAL PARAMETERS
  options("cyto_plot_par" = list())
  
  # CYTO_PLOT() RESET GRAPHICAL PARAMETERS
  options("cyto_plot_par_reset" = list())
  
  # CYTO_PLOT() DATA ALREADY PREPARED
  options("cyto_plot_data" = FALSE)
  
  # REGISTER OPENCYTO PREPROCESSING & GATING FUNCTIONS
  suppressMessages({
    # MANUAL GATING PLUGIN
    openCyto::register_plugins(
      fun = .cyto_gate_manual, 
      methodName = "cyto_gate_manual"
    )
    # CYTO_GATE_DRAW
    openCyto::register_plugins(
      fun = .cyto_gate_draw, 
      methodName = "cyto_gate_draw"
    )
    openCyto::register_plugins(
      fun = .pp_cyto_gate_draw, 
      methodName = "pp_cyto_gate_draw", 
      dep = NA,
      "preprocessing"
    )
    # CYTO_GATE_SAMPLE
    openCyto::register_plugins(
      fun = .pp_cyto_gate_sample, 
      methodName = "pp_cyto_gate_sample", 
      dep = NA,
      "preprocessing"
    )
    openCyto::register_plugins(
      fun = .cyto_gate_sample,
      methodName = "cyto_gate_sample"
    )
    # CYTO_GATE_CLEAN
    openCyto::register_plugins(
      fun = .cyto_gate_clean,
      methodName = "cyto_gate_clean"
    )
    # CYTO_GATE_CLUST
    openCyto::register_plugins(
      fun = .pp_cyto_gate_clust, 
      methodName = "pp_cyto_gate_clust", 
      dep = NA,
      "preprocessing"
    )
    openCyto::register_plugins(
      fun = .cyto_gate_clust,
      methodName = "cyto_gate_clust"
    )
  })
  
}

# # Docker Message
# .onAttach <- function(libname, pkgname){
#   packageStartupMessage("CytoExploreR Docker v1.0.5")
# }