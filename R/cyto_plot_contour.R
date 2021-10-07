## CYTO_PLOT_CONTOUR -----------------------------------------------------------

#' Add contour lines to cyto_plot
#'
#' @param x either an object of class
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} or a list containing
#'   cytoframe or cytoset objects.
#' @param parent name of the parent population to extract from GatingHierarchy
#'   or GatingSt objects, set to the "root" node by default.
#' @param overlay list containing cytoframe or cytoset objects to overlay.
#' @param channels channels used to construct the existing cyto_plot.
#' @param events controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 15 by default.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 1 by default.
#' @param contour_line_col colour to use for contour lines, set to "black" by
#'   default.
#' @param contour_line_alpha numeric [0,1] to control transparency of contour
#'   lines, set to 1 by default to remove transparency.
#' @param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'   sampling is applied with each \code{\link{cyto_plot_contour}} call, set to
#'   an arbitrary numeric by default. This behaviour can be turned off by
#'   setting this argument to NULL.
#' @param ... not in use.
#'
#' @importFrom MASS kde2d
#' @importFrom graphics contour par
#' @importFrom grDevices adjustcolor
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(
#'   system.file(
#'     "extdata/Activation-GatingSet",
#'     package = "CytoExploreRData"
#'   )
#' )
#'
#' # Visualise T Cells
#' cyto_plot(gs[[1]],
#' parent = "T Cells",
#' channels = c("CD4", "CD8"))
#'
#' # Add contour lines
#' cyto_plot_contour(gs[[1]],
#' parent = "T Cells",
#' channels = c("CD4", "CD8"),
#' contour_lines = 15)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_plot_contour <- function(x,
                              parent = "root",
                              channels,
                              overlay = NA,
                              events = 1,
                              contour_lines = 15,
                              contour_line_type = 1,          
                              contour_line_width = 1,
                              contour_line_col = "black",
                              contour_line_alpha = 1,
                              seed = 42,
                              ...) {
  
  # TODO: BYPASS BKDE2D COMPUTATION
  
  # CHECKS -------------------------------------------------------------------
  
  # X - CYTOFRAME/CYTOSET/GATINGHIERARCHY/GATINGSET
  if(cyto_class(x, c("flowFrame", "flowSet", "GatingSet"))) {
    x <- cyto_data_extract(x,
                           parent = parent,
                           format = "cytoset", 
                           copy = FALSE)
    # OVERLAY
    if(!.all_na(overlay)) {
      if(!cyto_class(overlay, "list")) {
        overlay <- cyto_list(overlay)
      }
      x <- c(x, overlay)
    }
    # SAMPLING
    if(all(events != 1)) {
      x <- cyto_sample(x,
                       events = events,
                       seed = seed)
    }
  # CYTO_PLOT ARGUMENTS
  } else if(cyto_class(x, "cyto_plot")) { # not used - call point instead
    .args_update(x)
  # CHECK LISTS
  } else if(!all(LAPPLY(x, cyto_class, c("flowFrame", "flowSet")))) {
    stop("'x' must be a list of cytoframes or cytosets!")
  }
  
  
  # CONTOUR_LINES --------------------------------------------------------------
  
  invisible(mapply(
    function(cs,
             contour_lines,
             contour_line_type,
             contour_line_width,
             contour_line_col,
             contour_line_alpha) {
      # COMPUTE BKDE2D  
      if(contour_lines != 0) {
        # COMPUTE BKDE2D
        bkde2d <- cyto_apply(
          cs,
          "cyto_stat_bkde2d",
          input = "matrix",
          channels = channels,
          bins = c(250, 250),
          limits = list(.par("usr")[[1]][1:2],
                        .par("usr")[[1]][3:4]),
          copy = FALSE,
          simplify = FALSE
        )[[1]]
        # PLOT CONTOUR LINES
        graphics::contour(
          z = bkde2d$bkde,
          x = bkde2d$bins$x,
          y = bkde2d$bins$y,
          add = TRUE,
          drawlabels = FALSE,
          nlevels = contour_lines,
          col = adjustcolor(contour_line_col, contour_line_alpha),
          lwd = contour_line_width,
          lty = contour_line_type
        )
      }
    }, 
    x,
    contour_lines,
    contour_line_type,
    contour_line_width,
    contour_line_col,
    contour_line_alpha
  ))
  
}
