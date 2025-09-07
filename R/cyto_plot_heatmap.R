## CYTO_PLOT_HEATMAP -----------------------------------------------------------

#' Plot heatmaps using HeatmapR
#'
#' A CytoExploreR wrapper for \code{HeatmapR::heat_map()} to allow visualisation
#' of matrices or data.frames as heatmaps.
#'
#' @param x object of class \code{"matrix"} or \code{"data.frame"}.
#' @param ... additional arguments passed to \code{HeatmapR::heat_map()}.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @return recorded heatmap.
#' 
#' @seealso \code{\link[HeatmapR]{heat_map}}
#' 
#' @examples
#' \dontrun{
#'   library(CytoExploreRData)
#'   gs <- GatingSet(Activation)
#'   spill <- cyto_spillover_extract(gs)[[1]]
#'   cyto_plot_heatmap(
#'     spill,
#'     cell_text = TRUE,
#'     title = "Spillover Matrix"
#'   )
#' }
#'
#' @export
cyto_plot_heatmap <- function(x,
                              ...) {
  
  # CYTOMETRY OBJECTS NOT SUPPORTED
  if(cyto_class(x, c("flowSet", "GatingSet"))) {
    stop(
      "cyto_plot_heatmap() expects objects of class matrix or data.frame!"
    )
  }
  
  # PLOT HEATMAP
  heat_map(
    x,
    ...
  )
  
}
