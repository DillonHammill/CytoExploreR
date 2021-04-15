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
#' @param display controls the number or percentage of events to display, set to
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
                              display = 1,
                              contour_lines = 15,
                              contour_line_type = 1,
                              contour_line_width = 1,
                              contour_line_col = "black",
                              contour_line_alpha = 1,
                              seed = 42,
                              ...) {
  
  # CHECKS -------------------------------------------------------------------
  
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
  # EXTRACT DATA
  if(cyto_class(x, "GatingSet")) {
    x <- cyto_data_extract(gs,
                           parent = parent,
                           channels = channels,
                           copy = FALSE)
  }
  
  # PREPARE DATA
  if(cyto_class(x,c("flowFrame", "flowSet"))) {
    # COMBINE LAYERS
    if(!.all_na(overlay)) {
      if(cyto_class(overlay, c("flowFrame", "flowSet"))) {
        overlay <- list(overlay)
      }
      x <- c(list(x), overlay)
    } else {
      x <- list(x)
    }
  # CYTO_PLOT ARGUMENTS
  } else if(cyto_class(x, "cyto_plot")) {
    .args_update(x)
  }
  
  # SAMPLING
  if(display != 1) {
    x <- cyto_sample(x,
                     display = display,
                     seed = seed)
  }
  
  # CONTOUR_LINES ------------------------------------------------------------
  
  invisible(mapply(
    function(cf,
             contour_lines,
             contour_line_type,
             contour_line_width,
             contour_line_col,
             contour_line_alpha) {
      # EXTRACT RAW DATA
      cf_exprs <- cyto_data_extract(cf, 
                                    format = "matrix", 
                                    channels = channels)[[1]]
      # MERGE CYTOSET MULTIPLE SAMPLES
      if(length(cf_exprs) > 1) {
        cf_exprs <- do.call("rbind", cf_exprs)
      } else {
        cf_exprs <- cf_exprs[[1]]
      }
      # ADD CONTOUR LINES
      if (contour_lines != 0) {
        # BYPASS INSUFFICIENT EVENTS
        if (nrow(cf_exprs) > 2) {
          # 2D KERNEL DENSITY
          z <- kde2d(
            x = cf_exprs[, 1],
            y = cf_exprs[, 2],
            n = 75,
            lims = par("usr")
          )
          # PLOT CONTOUR LINES
          graphics::contour(
            z = z$z,
            x = z$x,
            y = z$y,
            add = TRUE,
            drawlabels = FALSE,
            nlevels = contour_lines,
            col = adjustcolor(contour_line_col, contour_line_alpha),
            lwd = contour_line_width,
            lty = contour_line_type
          )
        }
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
