## CYTO_PLOT_GRID --------------------------------------------------------------

#' @name cyto_plot_grid
NULL

#' @noRd
#' @export
cyto_plot_grid <- function(x, 
                           ...) {
  UseMethod("cyto_plot_grid")
}

#' @rdname cyto_plot_grid
#' @export
cyto_plot_grid.GatingSet <- function(x) {
  
}

#' @rdname cyto_plot_grid
#' @export
cyto_plot_grid.flowSet <- function(x) {
  
}