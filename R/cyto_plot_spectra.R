## CYTO_PLOT_SPECTRA -----------------------------------------------------------

#' Plot fluorescent spectra for each sample or fluorochrome across detectors
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
#' @rdname cyto_plot_spectra
#'
#' @export
cyto_plot_spectra <- function(x, ...) {
  UseMethod("cyto_plot_spectra")
}

#' @rdname cyto_plot_spectra
#' @export
cyto_plot_spectra.default <- function(x,
                                      parent,
                                      channels,
                                      axes_trans = NA,
                                      axes_limits = "machine",
                                      spectra_col_scale = NA,
                                      spectra_col = NA,
                                      spectra_cols = NA,
                                      spectra_col_alpha = 1,
                                      ...) {
  
  # NOTE: PLOT SPECTRAL BANDS FOR CYTOSET/GATINGSET OBJECTS
  
  # TODO: SPACING BETWEEN BANDS & THRESHOLDING
  
  # PULL DOWN ARGUMENTS
  args <- .args_list(...)
  
  # CHANNELS
  if(.empty(args$channels)) {
    args$channels <- cyto_channels(
      args$x,
      exclude = c(
        "FSC",
        "SSC",
        "Time",
        "Event-ID",
        "Sample-ID"
      )
    )
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(args$x)
  
  # AXES_TRANS
  if(.all_na(args$axes_trans)){
    args$axes_trans <- cyto_transformers_extract(
      args$x
    )
  }
  
  # EXTRACT PARENTAL POPULATIONS PER CONTROL
  if(missing(parent)) {
    if("parent" %in% colnames(pd)) {
      parent <- pd[, "parent"]
    } else {
      parent <- "root"
    }
  }
  
  # PARENT PER CONTROL
  parent <- rep(
    parent,
    length.out = length(x)
  )
  
  # EXTRACT PARENT PER CONTROL
  args$x <- cytoset(
    structure(
      lapply(
        seq_along(x),
        function(z) {
          pop <- parent[z]
          # CHECK PARENT EXISTS
          if(cyto_class(x, "GatingSet")) {
            pop <- tryCatch(
              cyto_nodes_convert(
                x,
                nodes = pop,
                path = "auto"
              ),
              error = function(e) {
                return("root")
              }
            )
          }
          # CYTOFRAME
          cyto_data_extract(
            x[z],
            parent = pop,
            format = "cytoset",
            copy = FALSE
          )[[1]][[1]]
        }
      ),
      names = cyto_names(x)
    )
  )
  
  # TODO: EXTEND BOTTOM MARGIN TO ENCOMPASS CHANNEL NAMES
  
  # CALL CYTO_PLOT()
  cyto_func_call(
    "cyto_plot",
    args
  )
  
}

#' Plot fluorescent spectrum for each fluorochrome across detectors
#'
#' @param x object of class \code{"matrix"} or a list of matrices.
#' @param channels names of the channels or markers to be included in the
#'   spectral profiles, set to all channels contained in the column names of the
#'   supplied matrices by default.
#' @param select names of the dyes to extract from the rows of each matrix for
#'   plotting, set to NULL by default to include all rows.
#' @param ylim limits for the y axis in percentage units, default ylim is
#'   computed using the range across all matrices.
#' @param xlab label to use for the x axis, set to NA by default to remove x
#'   axis label making more space for x axis text.
#' @param ylab label to use for the y axis, set to Emission percentage by
#'   default.
#' @param legend_text text to use in the legend, defaults to the names of the
#'   matrices as supplied in \code{x}.
#' @param ... additional arguments passed to \code{\link{cyto_plot}} including
#'   all the \code{gate_} aesthetics arguments that can be used to customise the
#'   plotted spectra.
#'
#' @return a list containing the recorded plots.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace cytoset
#'
#' @seealso \code{\link{cyto_plot}}
#' @seealso \code{\link{cyto_unmix_compute}}
#' @seealso \code{\link{cyto_unmix}}
#'
#' @export
cyto_plot_spectra.list <- function(x,
                                   channels,
                                   select = NULL,
                                   ylim = c(NA, NA),
                                   xlab = NA,
                                   ylab = NA,
                                   legend_text = NA, 
                                   ...) {
  
  # NOTE: PLOT SPECTRAL PROFILES FOR SPILLOVER/UNMIXING MATRICES
  
  # TODO: ADD SUPPORT FOR LEGEND
  # TODO: ADD CUSTOMISATION
  
  # EXPECT MATRIX OR LIST OF MATRICES
  if(is.matrix(x)) {
    x <- list("spectra" = x)
  } else if(is.list(x)) {
    if(!all(LAPPLY(x, cyto_class, expect = "matrix", class = FALSE))) {
      stop(
        "'x' must be either a matrix or list of matrices!"
      )
    }
  } else {
    stop(
      "'x' must be either a matrix or list of matrices!"
    )
  }
  
  # PULL DOWN ARGUMENTS
  args <- .args_list(...)
  
  # CHANNELS
  if(.empty(args$channels)) {
    args$channels <- cyto_channels(
      args$x[[1]]
    )
  }
  
  # DYES MUST MATCH IN ALL MATRICES
  dyes <- lapply(
    args$x,
    rownames
  )
  dyes <- unique(dyes)
  if(length(dyes) != 1) {
    stop(
      "All matrices must contain the same set of dyes as rownames!"
    )
  }
  dyes <- dyes[[1]]
  
  # EXPERIMENT NAMES
  nms <- names(x)
  
  # YLIM -> CONVERTED TO PERCENTAGES
  rng <- do.call("range", args$x) * 100
  if(any(is.na(args$ylim))) {
    args$ylim[is.na(args$ylim)] <- rng[is.na(args$ylim)]
  }
  
  # CREATE (EMPTY) CYTOSET PER DYE
  cs <- cytoset(
    structure(
      lapply(
        dyes,
        function(dye) {
          cyto_empty(
            dye,
            args$channels
          )
        }
      ),
      names = dyes
    )
  )
  
  # FILTER DYES
  cs <- cyto_select(
    cs,
    select
  )
  
  # NOTE: DUPLICATE GATES WILL BE DROPPED SO WE NEED DISTINCT FILTERIDS
  # CREATE SPECTRAL GATES - PERCENTAGES
  args$gate <- lapply(
    cyto_names(cs),
    function(dye) {
      lapply(
        seq_along(args$x),
        function(id) {
          coords <- cbind(
            "x" = seq_along(args$channels),
            "y" = as.numeric(
                args$x[[id]][match(dye, rownames(args$x[[id]])), ]
            ) * 100
          )
          rownames(coords) <- args$channels
          spectralGate(
            coords = coords,
            filterId = if(is.null(nms[id])) {
              paste0(dye, "-", id)
            } else {
              paste0(dye, "-", nms[id])
            }
          )
        }
      )
    }
  )
  
  # LEGEND_TEXT
  if(.all_na(args[["legend_text"]])) {
    args$legend_text <- nms
  }
  
  # INSERT CYTOSET INTO PLOT ARGUMENTS
  args$x <- cs
  
  # SIGNAL SPECTRAL PLOT
  args$spectra <- TRUE
  
  # REMOVE GATE LABELS
  args$label <- FALSE
  
  # CALL CYTO_PLOT()
  cyto_func_call(
    "cyto_plot",
    args
  )
  
}

#' Internal spectralGate class to plot spectral profiles using cyto_plot
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setClass(
  "spectralGate",
  representation(
    coords = "matrix"
  ),
  contains = "parameterFilter",
  prototype = list(
    filterId = "defaultSpectralGate",
    coords = matrix(ncol = 2, nrow = 3)
  ),
  validity = function(object)
  {
    msg <- TRUE
    if(!is.matrix(object@coords) || nrow(object@coords) < 3 ||
       ncol(object@coords) != 2
    ) {
      stop("\nslot 'boundaries' must be a numeric matrix",
                   "of at least 3 rows and exactly 2 columns")
    }
  }
)

#' Constructor for internal spectralGate class
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
spectralGate <- function(coords,
                         filterId = "defaultSpectralGate") {
  # NOTE: COLUMN NAMES DON'T MATTER -> SHOULD BE XY
  if(ncol(coords) > 2) {
    stop(
      "Spectral gates must be constructed in using two parameters only!"
    )
  }
  if(is.null(rownames(coords))) {
    stop(
      "Channel names MUST be specified as rownames for 'coords'!"
    )
  }
  new(
    "spectralGate", 
    filterId = filterId, 
    parameters = rownames(coords),
    coords = coords
  )
}
