## CYTO_PLOT_CONTOUR -----------------------------------------------------------

#' Add contour lines to cyto_plot
#'
#' @param x either an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}}, or
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or a list of
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}} objects.
#' @param channels channels used to construct the existing cyto_plot.
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
#' @importFrom flowCore exprs
#' @importFrom grDevices adjustcolor
#' @importFrom methods is
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
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_plot_contour
NULL

#' @rdname cyto_plot_contour
#' @export
cyto_plot_contour <- function(x, ...){
  UseMethod("cyto_plot_contour")
}

#' @rdname cyto_plot_contour
#' @export
cyto_plot_contour.flowFrame <- function(x,
                                        channels,
                                        contour_lines = 15,
                                        contour_line_type = 1,
                                        contour_line_width = 1,
                                        contour_line_col = "black",
                                        contour_line_alpha = 1,
                                        ...){
  
  # CHECKS -------------------------------------------------------------------
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
  # PREPARE DATA -------------------------------------------------------------
  fr <- x
  fr_exprs <- exprs(fr)[, channels]
  
  # CONTOUR_LINES ------------------------------------------------------------
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
  
}

#' @rdname cyto_plot_contour
#' @export
cyto_plot_contour.flowSet <- function(x,
                                      channels,
                                      contour_lines = 15,
                                      contour_line_type = 1,
                                      contour_line_width = 1,
                                      contour_line_col = "black",
                                      contour_line_alpha = 1,
                                      ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # LIST OF FLOWFRAMES
  fr_list <- lapply(seq_len(length(x)), function(z) {
    x[[z]]
  })
  
  # CONTOUR_LINES --------------------------------------------------------------
  invisible(mapply(
    function(fr,
             contour_lines,
             contour_line_type,
             contour_line_width,
             contour_line_col,
             contour_line_alpha) {
      cyto_plot_contour(fr,
                        channels = channels,
                        contour_lines = contour_lines,
                        contour_line_type = contour_line_type,
                        contour_line_width = contour_line_width,
                        contour_line_col = contour_line_col,
                        contour_line_alpha = contour_line_alpha
      )
    }, fr_list,
    contour_lines,
    contour_line_type,
    contour_line_width,
    contour_line_col,
    contour_line_alpha
  ))
}

#' @rdname cyto_plot_contour
#' @export
cyto_plot_contour.list <- function(x,
                                   channels,
                                   contour_lines = 15,
                                   contour_line_type = 1,
                                   contour_line_width = 1,
                                   contour_line_col = "black",
                                   contour_line_alpha = 1,
                                   ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # CLASS
  if (!all(lapply(x, is, "flowFrame"))) {
    stop("x should be a list of flowFrame objects.")
  }
  
  # CHANNELS 
  channels <- cyto_channels_extract(x[[1]], channels)
  
  # CONTOUR_LINES --------------------------------------------------------------
  invisible(mapply(
    function(fr,
             contour_lines,
             contour_line_type,
             contour_line_width,
             contour_line_col,
             contour_line_alpha) {
      cyto_plot_contour(fr,
                        channels = channels,
                        contour_lines = contour_lines,
                        contour_line_type = contour_line_type,
                        contour_line_width = contour_line_width,
                        contour_line_col = contour_line_col,
                        contour_line_alpha = contour_line_alpha
      )
    }, x,
    contour_lines,
    contour_line_type,
    contour_line_width,
    contour_line_col,
    contour_line_alpha
  ))
  
}
