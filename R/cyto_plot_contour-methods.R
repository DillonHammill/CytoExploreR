#' Contour Lines for cyto_plot
#'
#' @param x object to use for contour lines to overlay onto an existing
#'   cyto_plot.
#' @param ... additional method-specific arguments
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(
  name = "cyto_plot_contour",
  def = function(x, ...) {
    standardGeneric("cyto_plot_contour")
  }
)

#' Contour Lines for cyto_plot - flowFrame Method
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}} to
#'   use for contour lines to overlay onto an existing cyto_plot.
#' @param channels channels used to construct the existing cyto_plot.
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 15 by default.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 2 by default.
#' @param contour_line_col colour to use for contour lines, set to "black" by
#'   default.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom MASS kde2d
#' @importFrom graphics contour par
#' @importFrom flowCore exprs
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in Samples
#' fs <- Activation
#' 
#' # Apply compensation
#' fs <- compensate(fs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(fs[[4]], cyto_fluor_channels(fs))
#' fs <- transform(fs, trans)
#' 
#' # Plot
#' cyto_plot(fs[[4]],
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   axes_trans = trans
#' )
#' 
#' # Contour lines
#' cyto_plot_contour(fs[[4]],
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   contour_lines = 20
#' )
#' @export
setMethod(cyto_plot_contour,
  signature = "flowFrame",
  definition = function(x,
                          channels,
                          contour_lines = 15,
                          contour_line_type = 1,
                          contour_line_width = 1,
                          contour_line_col = "black") {

    # Assign x to fr
    fr <- x
    fr_exprs <- exprs(fr)[, channels]

    # Contours
    if (contour_lines != 0) {

      # Bypass contours if there are insufficient events
      if (nrow(fr_exprs) > 2) {

        # Calculate 2D kernel density using kde2d from MASS
        z <- MASS::kde2d(
          x = fr_exprs[, 1],
          y = fr_exprs[, 2],
          n = 75,
          lims = par("usr")
        )

        # Add contour lines to plot
        graphics::contour(
          z = z$z,
          x = z$x,
          y = z$y,
          add = TRUE,
          drawlabels = FALSE,
          nlevels = contour_lines,
          col = contour_line_col,
          lwd = contour_line_width,
          lty = contour_line_type
        )
      }
    }
  }
)

#' Contour Lines for cyto_plot - flowSet Method
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} to
#'   use for contour lines to overlay onto an existing cyto_plot.
#' @param channels channels used to construct the existing cyto_plot.
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 15 by default.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 2 by default.
#' @param contour_line_col colour to use for contour lines, set to "black" by
#'   default.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in Samples
#' fs <- Activation
#' 
#' # Apply compensation
#' fs <- compensate(fs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(fs[[4]], cyto_fluor_channels(fs))
#' fs <- transform(fs, trans)
#' 
#' # Plot
#' cyto_plot(fs,
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   group_by = "all", axes_trans = trans
#' )
#' 
#' # Contour lines for control sample
#' cyto_plot_contour(fs,
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   contour_lines = c(20, 0, 0, 0)
#' )
#' @export
setMethod(cyto_plot_contour,
  signature = "flowSet",
  definition = function(x,
                          channels,
                          contour_lines = 15,
                          contour_line_type = 1,
                          contour_line_width = 1,
                          contour_line_col = "black") {

    # Assign x to fs
    fs <- x

    # Convert to list of flowFrames
    fr.lst <- lapply(seq_len(length(fs)), function(x) {
      fs[[x]]
    })

    # Contours
    invisible(mapply(
      function(fr,
                     contour_lines,
                     contour_line_type,
                     contour_line_width,
                     contour_line_col) {
        cyto_plot_contour(fr,
          channels = channels,
          contour_lines = contour_lines,
          contour_line_type = contour_line_type,
          contour_line_width = contour_line_width,
          contour_line_col = "black"
        )
      }, fr.lst,
      contour_lines,
      contour_line_type,
      contour_line_width,
      contour_line_col
    ))
  }
)

#' Contour Lines for cyto_plot - list Method
#'
#' @param x a list of \code{flowFrame} objects to use for contour lines to
#'   overlay onto an existing cyto_plot.
#' @param channels channels used to construct the existing cyto_plot.
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 15 by default.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 2 by default.
#' @param contour_line_col colour to use for contour lines, set to "black" by
#'   default.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in Samples
#' fs <- Activation
#' 
#' # Apply compensation
#' fs <- compensate(fs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(fs[[4]], cyto_fluor_channels(fs))
#' fs <- transform(fs, trans)
#' 
#' # Plot
#' cyto_plot(fs,
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   group_by = "all", axes_trans = trans
#' )
#' 
#' # Contour lines for control sample
#' cyto_plot_contour(list(fs[[1]], fs[[2]], fs[[3]], fs[[4]]),
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   contour_lines = c(20, 0, 0, 0)
#' )
#' @export
setMethod(cyto_plot_contour,
  signature = "list",
  definition = function(x,
                          channels,
                          contour_lines = 15,
                          contour_line_type = 1,
                          contour_line_width = 1,
                          contour_line_col = "black") {
    # Check class of x
    if (!all(lapply(x, "class") == "flowFrame")) {
      stop("x should be a list of flowFrame objects.")
    }

    # Assign x to fr.lst
    fr.lst <- x

    # Contours
    invisible(mapply(
      function(fr,
                     contour_lines,
                     contour_line_type,
                     contour_line_width,
                     contour_line_col) {
        cyto_plot_contour(fr,
          channels = channels,
          contour_lines = contour_lines,
          contour_line_type = contour_line_type,
          contour_line_width = contour_line_width,
          contour_line_col = "black"
        )
      }, fr.lst,
      contour_lines,
      contour_line_type,
      contour_line_width,
      contour_line_col
    ))
  }
)
