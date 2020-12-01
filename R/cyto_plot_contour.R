## CYTO_PLOT_CONTOUR -----------------------------------------------------------

#' Add contour lines to cyto_plot
#'
#' @param x either an object of class
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}} or a list of
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}} objects.
#' @param overlay optional argument if x is a cytoframe to overlay a list of
#'   cytoframe objects.
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
#' @param ... not in use.
#'
#' @importFrom MASS kde2d
#' @importFrom graphics contour par
#' @importFrom grDevices adjustcolor
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in Samples
#' fs <- Activation
#'
#' # Apply compensation
#' fs <- cyto_compensate(fs)
#'
#' # Transform fluorescent channels
#' trans <- cyto_transformer_logicle(fs)
#' fs <- cyto_transform(fs, trans)
#'
#' # Plot
#' cyto_plot(fs[[32]],
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   axes_trans = trans
#' )
#'
#' # Contour lines
#' cyto_plot_contour(fs[[32]],
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   contour_lines = 20
#' )
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_plot_contour <- function(x,
                              channels,
                              overlay = NA,
                              display = 1,
                              contour_lines = 15,
                              contour_line_type = 1,
                              contour_line_width = 1,
                              contour_line_col = "black",
                              contour_line_alpha = 1,
                              ...) {
  
  # CHECKS -------------------------------------------------------------------
  
  # X - FLOWFRAME/FLOWFRAME LIST/FLOWSET/CYTO_PLOT ARGS
  if(cyto_class(x, "flowFrame")) {
    overlay <- cyto_list(overlay)
    if(!.all_na(overlay)) {
      x <- c(structure(list(x), names = cyto_names(x)),
             overlay)
    } else {
      x <- cyto_list(x)
    }
  } else if(cyto_class(x, "flowSet")) {
    x <- cyto_list(x)
  } else if(cyto_class(x, "cyto_plot")) {
    .args_update(x)
  } else if(!all(LAPPLY(x, cyto_class, "flowFrame"))) {
    stop("'x' must be a list of cytoframe objects!")
  }
  
  # CHANNELS
  channels <- cyto_channels_extract(x[[1]], channels)
  
  # SAMPLING
  if(display != 1) {
    x <- cyto_sample(x,
                     display = display,
                     seed = 56)
  }
  
  # CONTOUR_LINES ------------------------------------------------------------
  
  invisible(mapply(
    function(fr,
             contour_lines,
             contour_line_type,
             contour_line_width,
             contour_line_col,
             contour_line_alpha) {
      fr_exprs <- cyto_data_extract(fr, 
                                    format = "matrix", 
                                    channels = channels)[[1]]
      if (contour_lines != 0) {
        # BYPASS INSUFFICIENT EVENTS
        if (nrow(fr_exprs) > 2) {
          # 2D KERNEL DENSITY
          z <- kde2d(
            x = fr_exprs[, 1],
            y = fr_exprs[, 2],
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
