#' Overlay Points on cyto_plot
#'
#' @param x object to overlay onto an existing cyto_plot.
#' @param ... additional method-specific arguments
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(
  name = "cyto_plot_overlay",
  def = function(x, ...) {
    standardGeneric("cyto_plot_overlay")
  }
)

#' Overlay Points on cyto_plot - flowFrame Method
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}} to
#'   overlay onto an existing cyto_plot.
#' @param channels channels used to construct the existing cyto_plot.
#' @param point_shape point character to use for points, set to "." by default
#'   to maximise plotting speed.
#' @param point_size numeric specifying the degree of character expansion for
#'   points, set to 2 by default.
#' @param point_col colours to use for points, set to \code{"black"} by default.
#' @param point_alpha numeric [0,1] used to control colour transparency, set to
#'   1 by default to remove transparency.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom graphics points
#' @importFrom flowCore exprs
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#' 
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#' 
#' # Gate using gate_draw
#' gating(Activation_gatingTemplate, gs)
#' 
#' # Plot
#' cyto_plot(gs[[4]],
#'   parent = "Live Cells",
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   axes_trans = trans
#' )
#' 
#' # Overlay CD4 T Cells
#' cyto_plot_overlay(getData(gs, "CD4 T Cells")[[4]],
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   point_col = "magenta"
#' )
#' @export
setMethod(cyto_plot_overlay,
          signature = "flowFrame",
          definition = function(x,
                                channels,
                                point_shape = ".",
                                point_size = 2,
                                point_col = "black",
                                point_alpha = 1) {
            
            # Assign x to fr
            fr <- x
            
            # Add points to existing plot
            graphics::points(
              x = exprs(fr)[, channels[1]],
              y = exprs(fr)[, channels[2]],
              pch = point_shape,
              col = adjustcolor(point_col, point_alpha),
              cex = point_size
            )
          }
)

#' Overlay Points on cyto_plot - flowSet Method
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} to
#'   overlay onto an existing cyto_plot.
#' @param channels channels used to construct the existing cyto_plot.
#' @param point_shape point character to use for points, set to "." by default
#'   to maximise plotting speed.
#' @param point_size numeric specifying the degree of character expansion for
#'   points, set to 2 by default.
#' @param point_col colours to use for points, set to \code{"black"} by default.
#' @param point_alpha numeric [0,1] used to control colour transparency, set to
#'   1 by default to remove transparency.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#' 
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#' 
#' # Gate using gate_draw
#' gating(Activation_gatingTemplate, gs)
#' 
#' # Plot
#' cyto_plot(gs[[4]],
#'   parent = "Live Cells",
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   axes_trans = trans
#' )
#' 
#' # Overlay CD4 T Cells
#' cyto_plot_overlay(getData(gs, "CD4 T Cells")[2:4],
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   point_col = c("magenta", "orange", "green")
#' )
#' @export
setMethod(cyto_plot_overlay,
          signature = "flowSet",
          definition = function(x,
                                channels,
                                point_shape = ".",
                                point_size = 2,
                                point_col = "black",
                                point_alpha = 1) {
            
            # Assign x to fs
            fs <- x
            
            # Convert to list of flowFrames
            fr.lst <- lapply(seq_len(length(fs)), function(x) {
              fs[[x]]
            })
            
            # Colours
            cols <- colorRampPalette(c(
              "darkorchid",
              "blueviolet",
              "magenta",
              "deeppink",
              "red4",
              "orange",
              "springgreen4"
            ))
            point_col <- c(point_col,cols(length(fr.lst)))[seq_len(length(fr.lst))]
            
            # Repeat other arguments
            point_shape <- rep(point_shape, length.out = length(fr.lst))
            point_size <- rep(point_size, length.out = length(fr.lst))
            point_alpha <- rep(point_alpha, length.out = length(fr.lst))
            
            # Add points to existing plot
            invisible(
              mapply(
                function(fr,
                         point_shape,
                         point_size,
                         point_col,
                         point_alpha) {
                  cyto_plot_overlay(
                    x = fr,
                    channels = channels,
                    point_shape = point_shape,
                    point_size = point_size,
                    point_col = point_col,
                    point_alpha = point_alpha
                  )
                }, fr.lst,
                point_shape,
                point_size,
                point_col,
                point_alpha
              )
            )
          }
)

#' Overlay Points on cyto_plot - list Method
#'
#' @param x list of \code{flowFrame} objects to overlay onto an existing
#'   cyto_plot.
#' @param channels channels used to construct the existing cyto_plot.
#' @param point_shape point character to use for points, set to "." by default
#'   to maximise plotting speed.
#' @param point_size numeric specifying the degree of character expansion for
#'   points, set to 2 by default.
#' @param point_col colours to use for points, set to \code{"black"} by default.
#' @param point_alpha numeric [0,1] used to control colour transparency, set to
#'   1 by default to remove transparency.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#' 
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#' 
#' # Gate using gate_draw
#' gating(Activation_gatingTemplate, gs)
#' 
#' # Plot
#' cyto_plot(gs[[4]],
#'   parent = "Live Cells",
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   axes_trans = trans
#' )
#' 
#' # Overlay CD4 T Cells & Dendritic Cells
#' cyto_plot_overlay(list(
#'   getData(gs, "CD4 T Cells")[[4]],
#'   getData(gs, "Dendritic Cells")[[4]]
#' ),
#' channels = c("APC-Cy7-A", "PE-A"),
#' point_col = c("magenta", "red")
#' )
#' @export
setMethod(cyto_plot_overlay,
          signature = "list",
          definition = function(x,
                                channels,
                                point_shape = ".",
                                point_size = 2,
                                point_col = "black",
                                point_alpha = 1) {
            
            # Check class of x
            if (!all(unlist(lapply(x, "class")) == "flowFrame")) {
              stop("x should be a list of flowFrame objects to overlay.")
            }
            
            # Assign x to fr.lst
            fr.lst <- x
            
            # Colours
            cols <- colorRampPalette(c(
              "darkorchid",
              "blueviolet",
              "magenta",
              "deeppink",
              "red4",
              "orange",
              "springgreen4"
            ))
            point_col <- c(point_col,cols(length(fr.lst)))[seq_len(length(fr.lst))]
            
            # Repeat other arguments
            point_shape <- rep(point_shape, length.out = length(fr.lst))
            point_size <- rep(point_size, length.out = length(fr.lst))
            point_alpha <- rep(point_alpha, length.out = length(fr.lst))
            
            # Add points to existing plot
            invisible(
              mapply(
                function(fr,
                         point_shape,
                         point_size,
                         point_col,
                         point_alpha) {
                  cyto_plot_overlay(
                    x = fr,
                    channels = channels,
                    point_shape = point_shape,
                    point_size = point_size,
                    point_col = point_col,
                    point_alpha = point_alpha
                  )
                }, fr.lst,
                point_shape,
                point_size,
                point_col,
                point_alpha
              )
            )
          }
)