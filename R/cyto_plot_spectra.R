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
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied data.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
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
#' @importFrom flowWorkspace cytoset
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
                              axes_trans = NA,
                              axes_limits = "machine",
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
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # AXES_TRANS
  if(.all_na(args$axes_trans)){
    args$axes_trans <- cyto_transformers_extract(
      args$x
    )
  }
  
  # EXTRACT PARENTAL POPULATIONS PER CONTROL
  if("parent" %in% colnames(pd)) {
    args$x <- cytoset(
      structure(
        lapply(
          seq_along(x),
          function(z) {
            # CYTOFRAME
            cyto_data_extract(
              x[z],
              parent = pd[z, "parent"],
              format = "cytoset",
              copy = FALSE
            )[[1]][[1]]
          }
        ),
        names = cyto_names(x)
      )
    )
  }
  
  # CALL CYTO_PLOT()
  cyto_func_call(
    "cyto_plot",
    args
  )
   
}
