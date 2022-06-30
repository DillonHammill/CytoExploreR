## CYTO_PLOT_SPECTRA -----------------------------------------------------------

#' Plot fluorescent spectrum for each fluorochrome across detectors
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied, set to the \code{"root"} node by
#'   default.
#' @param channels names of the channels or markers to be included in the
#'   spectral profiles, set to all channels except \code{"Time"}, \code{"FSC"}
#'   and \code{"SSC"} by default.
#' @param spectra_col_scale vector of ordered colours to use for the density
#'   colour gradient of spectra, matches \code{point_col_scale} by default.
#' @param spectra_col colour(s) to use for spectra, set to NA by default to use
#'   \code{spectra_col_scale}.
#' @param spectra_cols vector colours to draw from when selecting colours for
#'   spectra if none are supplied to \code{spectra_col}.
#' @param spectra_col_alpha numeric [0,1] to control point colour transparency
#'   of spectra, set to 1 by default to use solid colours.
#' @param ... additional arguments passed to \code{\link{cyto_plot}}.
#'
#' @return a list containing the recorded plots.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#' # gs - transformed GatingSet of unmixed spectral data
#' cyto_plot_sepctra(
#'   gs,
#'   parent = "root"
#' )
#' }
#'
#' @seealso \code{\link{cyto_plot}}
#' @seealso \code{\link{cyto_unmix_compute}}
#' @seealso \code{\link{cyto_unmix}}
#'
#' @export
cyto_plot_spectra <- function(x,
                              parent = "root",
                              channels = NULL,
                              spectra_col_scale = NA,
                              spectra_col = NA,
                              spectra_cols = NA,
                              spectra_col_alpha = 1,
                              ...) {
 
  # PULL DOWN ARGUMENTS
  args <- .args_list(...)
  
  # CHANNELS
  if(.empty(args$channels)) {
    args$channels <- cyto_channels(
      args$x,
      exclude = c(
        "FSC",
        "SSC",
        "Time"
      )
    )
  }
  
  # CALL CYTO_PLOT()
  cyto_func_call(
    "cyto_plot",
    args
  )
   
}
