#' Plot Gate Objects onto an Existing Plot
#'
#' @param x gate object of class
#'   \code{\link[flowCore:rectangleGate-class]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate-class]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}}, \code{list} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot_gate,rectangleGate-method}}
#' @seealso \code{\link{cyto_plot_gate,polygonGate-method}}
#' @seealso \code{\link{cyto_plot_gate,ellipsoidGate-method}}
#' @seealso \code{\link{cyto_plot_gate,list-method}}
#' @seealso \code{\link{cyto_plot_gate,filters-method}}
#'
#' @export
setGeneric(
  name = "cyto_plot_gate",
  def = function(x, ...) {
    standardGeneric("cyto_plot_gate")
  }
)

#' Plot rectangleGate Objects onto an Existing plot
#'
#' @param x an object of class
#'   \code{\link[flowCore:rectangleGate-class]{rectangleGate}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param gate_point logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param gate_point_shape integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param gate_point_size numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @return gate object with modified co-ordinates for plotting.
#'
#' @importFrom flowCore parameters
#' @importFrom graphics par rect lines abline points
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot_gate,polygonGate-method}}
#' @seealso \code{\link{cyto_plot_gate,ellipsoidGate-method}}
#' @seealso \code{\link{cyto_plot_gate,list-method}}
#' @seealso \code{\link{cyto_plot_gate,filters-method}}
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
#'   parent = "T Cells",
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   axes_trans = trans
#' )
#' 
#' # CD4 T Cells gate
#' cyto_plot_gate(getGate(gs, "CD4 T Cells")[[1]],
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   gate_line_col = "purple"
#' )
#' @export
setMethod(cyto_plot_gate, signature = "rectangleGate", definition = function(x,
                                                                             channels,
                                                                             gate_line_type = 1,
                                                                             gate_line_width = 2.5,
                                                                             gate_line_col = "red",
                                                                             gate_point = FALSE,
                                                                             gate_point_shape = 16,
                                                                             gate_point_size = 1) {

  # Assign x to gt
  gt <- x

  # Allow 1D gate plotted in 2D
  if (!missing(channels)) {
    if (length(channels) == 2 & length(parameters(gt)) == 1) {
      rg <- matrix(c(as.numeric(gt@min), as.numeric(gt@max), -Inf, Inf), ncol = 2, nrow = 2)
      colnames(rg) <- c(as.vector(parameters(gt)), channels[!channels == as.vector(parameters(gt))])
      rownames(rg) <- c("min", "max")
      gt <- rectangleGate(.gate = rg)
    }
  } else {
    channels <- as.vector(parameters(gt))
  }

  # 2D gatre in 1D
  if (length(channels) == 1 & length(parameters(gt)) == 2) {
    gt <- gt[channels]
  }

  if (!all(as.vector(parameters(gt)) %in% channels)) {
    stop("Channels used to construct the plot do not match those of the supplied gate.")
  }

  # 2D rectangleGate
  if (length(gt@min) == 2) {

    # Replace -Inf x values for plotting
    if (is.infinite(gt@min[channels[1]]) | gt@min[channels[1]] < (par("usr")[1] - 0.13 * par("usr")[1])) {
      gt@min[channels[1]] <- (par("usr")[1] - 0.13 * par("usr")[1])
    }

    # Replace Inf x values for plotting
    if (is.infinite(gt@max[channels[1]]) | gt@max[channels[1]] > 0.98 * par("usr")[2]) {
      gt@max[channels[1]] <- 0.98 * par("usr")[2]
    }

    # Replace -Inf y values for plotting
    if (is.infinite(gt@min[channels[2]]) | gt@min[channels[2]] < (par("usr")[3] - 0.13 * par("usr")[3])) {
      gt@min[channels[2]] <- (par("usr")[3] - 0.13 * par("usr")[3])
    }

    # Replace Inf y values for plotting
    if (is.infinite(gt@max[channels[2]]) | gt@max[channels[2]] > 0.98 * par("usr")[4]) {
      gt@max[channels[2]] <- 0.98 * par("usr")[4]
    }

    # Add points to gate
    if (gate_point == TRUE) {
      points(x = c(gt@min[channels[1]], gt@max[channels[1]]), y = c(gt@min[channels[2]], gt@max[channels[2]]), col = gate_line_col, pch = gate_point_shape, cex = gate_point_size)
    }

    rect(xleft = gt@min[channels[1]], ybottom = gt@min[channels[2]], xright = gt@max[channels[1]], ytop = gt@max[channels[2]], border = gate_line_col, lwd = gate_line_width, lty = gate_line_type)
  } else if (length(gt@min) == 1) {

    # Replace -Inf values for plotting
    if (is.infinite(gt@min[1]) | gt@min[1] < (par("usr")[1] - 0.13 * par("usr")[1])) {
      gt@min[1] <- (par("usr")[1] - 0.13 * par("usr")[1])
    }

    # Replace Inf values for plotting
    if (is.infinite(gt@max[1]) | gt@max[1] > 0.98 * par("usr")[2]) {
      gt@max[1] <- 0.98 * par("usr")[2]
    }

    # Add points (x1,hln) and (x2, hln)
    if (gate_point == TRUE) {
      points(x = c(gt@min[channels[1]], gt@min[channels[1]], gt@max[channels[1]], gt@max[channels[1]]), y = c(gt@min[channels[2]], gt@min[channels[2]], gt@max[channels[2]], gt@max[channels[2]]), col = gate_line_col, pch = gate_point_shape, cex = gate_point_size)
    }

    # Add rectangle
    rect(xleft = gt@min[channels[1]], ybottom = 0.6 * par("usr")[3], xright = gt@max[channels[1]], ytop = 0.985 * par("usr")[4], border = gate_line_col, lwd = gate_line_width, lty = gate_line_type)
  }

  invisible(gt)
})

#' Plot polygonGate Objects onto an Existing Plot
#'
#' @param x an object of class \code{\link[flowCore:polygonGate-class]{polygonGate}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param gate_point logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param gate_point_shape integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param gate_point_size numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @return gate object with modified co-ordinates for plotting.
#'
#' @importFrom flowCore parameters
#' @importFrom graphics par points polygon
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot_gate,rectangleGate-method}}
#' @seealso \code{\link{cyto_plot_gate,ellipsoidGate-method}}
#' @seealso \code{\link{cyto_plot_gate,list-method}}
#' @seealso \code{\link{cyto_plot_gate,filters-method}}
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
#'   parent = "root",
#'   channels = c("FSC-A", "SSC-A")
#' )
#' 
#' # Cells gate
#' cyto_plot_gate(getGate(gs, "Cells")[[1]],
#'   channels = c("FSC-A", "SSC-A"),
#'   gate_line_col = "purple"
#' )
#' @export
setMethod(cyto_plot_gate, signature = "polygonGate", definition = function(x,
                                                                           channels,
                                                                           gate_line_type = 1,
                                                                           gate_line_width = 2.5,
                                                                           gate_line_col = "red",
                                                                           gate_point = FALSE,
                                                                           gate_point_shape = 16,
                                                                           gate_point_size = 1) {

  # Assign x to gt
  gt <- x

  # Check Channels
  if (missing(channels)) {
    channels <- as.vector(parameters(gt))
  }

  if (!all(as.vector(parameters(gt)) %in% channels)) {
    stop("Channels used to construct the plot do not match those of the supplied gate.")
  }

  # Replace Inf values with plot limits
  if (!all(is.finite(gt@boundaries))) {
    cnt <- 0
    lapply(1:length(channels), function(x) {
      cnt <<- cnt + 1

      if (any(!is.finite(gt@boundaries[, channels[x]]) & any(gt@boundaries[, channels[x]] < 0))) {
        if (cnt == 1) {
          gt@boundaries[, channels[x]][which(gt@boundaries[, channels[x]] < 0)] <<- par("usr")[1]
        } else if (cnt == 2) {
          gt@boundaries[, channels[x]][which(gt@boundaries[, channels[x]] < 0)] <<- par("usr")[3]
        }
      }

      if (any(!is.finite(gt@boundaries[, channels[x]]) & any(!gt@boundaries[, channels[x]] < 0))) {
        if (cnt == 1) {
          gt@boundaries[, channels[x]][which(!is.finite(gt@boundaries[, channels[x]]) & !gt@boundaries[, channels[x]] < 0)] <<- par("usr")[2]
        } else if (cnt == 2) {
          gt@boundaries[, channels[x]][which(!is.finite(gt@boundaries[, channels[x]]) & !gt@boundaries[, channels[x]] < 0)] <<- par("usr")[4]
        }
      }
    })
  }

  # Plot Gate
  if (gate_point == TRUE) {
    points(x = c(gt@boundaries[, channels[1]]), y = c(gt@boundaries[, channels[2]]), pch = gate_point_shape, col = gate_line_col, cex = gate_point_size)
  }

  polygon(gt@boundaries[, channels[1]], gt@boundaries[, channels[2]], border = gate_line_col, lwd = gate_line_width, lty = gate_line_type)

  invisible(gt)
})

#' Plot ellipsoidGate Objects onto an Existing Plot
#'
#' @param x an object of class
#'   \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param gate_point logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param gate_point_shape integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param gate_point_size numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @return gate object with modified co-ordinates for plotting.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowCore parameters
#' @importFrom methods as
#' @importFrom graphics polygon
#'
#' @seealso \code{\link{cyto_plot_gate,rectangleGate-method}}
#' @seealso \code{\link{cyto_plot_gate,polygonGate-method}}
#' @seealso \code{\link{cyto_plot_gate,list-method}}
#' @seealso \code{\link{cyto_plot_gate,filters-method}}
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
#' # T Cells gate
#' cyto_plot_gate(getGate(gs, "T Cells")[[1]],
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   gate_line_col = "purple"
#' )
#' @export
setMethod(cyto_plot_gate, signature = "ellipsoidGate", definition = function(x,
                                                                             channels,
                                                                             gate_line_type = 1,
                                                                             gate_line_width = 2.5,
                                                                             gate_line_col = "red",
                                                                             gate_point = FALSE,
                                                                             gate_point_shape = 16,
                                                                             gate_point_size = 1) {

  # Assign x to gt
  gt <- x

  # Check Channels
  if (missing(channels)) {
    channels <- as.vector(parameters(gt))
  }

  if (!all(as.vector(parameters(gt)) %in% channels)) {
    stop("Channels used to construct the plot do not match those of the supplied gate.")
  }

  # Coerce to polygonGate
  gt <- as(gt, "polygonGate")

  # Plot gate
  polygon(gt@boundaries[, channels[1]], gt@boundaries[, channels[2]], border = gate_line_col, lwd = gate_line_width, lty = gate_line_type)

  invisible(gt)
})

#' Plot List of Gate Objects onto an Existing Plot
#'
#' @param x an object of class \code{list} containing objects of class
#'   \code{rectangleGate}, \code{polygonGate} or \code{ellipsoidGate}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param gate_point logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param gate_point_shape integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param gate_point_size numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @return gate object with modified co-ordinates for plotting.
#'
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot_gate,rectangleGate-method}}
#' @seealso \code{\link{cyto_plot_gate,polygonGate-method}}
#' @seealso \code{\link{cyto_plot_gate,ellipsoidGate-method}}
#' @seealso \code{\link{cyto_plot_gate,filters-method}}
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
#' # T Cells & Dendritic Cells gate
#' cyto_plot_gate(list(getGate(gs, "T Cells")[[1]], getGate(gs, "Dendritic Cells")[[1]]),
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   gate_line_col = c("purple", "magenta")
#' )
#' @export
setMethod(cyto_plot_gate, signature = "list", definition = function(x,
                                                                    channels,
                                                                    gate_line_type = 1,
                                                                    gate_line_width = 2.5,
                                                                    gate_line_col = "red",
                                                                    gate_point = FALSE,
                                                                    gate_point_shape = 16,
                                                                    gate_point_size = 1) {

  # Assign x to gts
  gts <- x

  # Check Channels
  if (missing(channels)) {
    channels <- unique(as.vector(sapply(gts, parameters)))
  }

  if (!any(unique(as.vector(sapply(gts, parameters))) %in% channels)) {
    stop("Channels used to construct the plot do not match those of the supplied gate.")
  }

  # Gate colours
  if (length(gate_line_col) != length(gts)) {
    if (length(gate_line_col) == 1) {
      gate_line_col <- rep(gate_line_col, length(gts))
    } else if (length(gate_line_col) > length(gts)) {
      gate_line_col <- gate_line_col[1:length(gts)]
    }
  }

  # Gate Line Type
  if (length(gate_line_type) != length(gts)) {
    if (length(gate_line_type) == 1) {
      gate_line_type <- rep(gate_line_type, length(gts))
    } else if (length(gate_line_type) > length(gts)) {
      gate_line_type <- gate_line_type[1:length(gts)]
    }
  }

  # Gate Line Width
  if (length(gate_line_width) != length(gts)) {
    if (length(gate_line_width) == 1) {
      gate_line_width <- rep(gate_line_width, length(gts))
    } else if (length(gate_line_width) > length(gts)) {
      gate_line_width <- gate_line_width[1:length(gts)]
    }
  }

  # Gate Point Character
  if (length(gate_point_shape) != length(gts)) {
    if (length(gate_point_shape) == 1) {
      gate_point_shape <- rep(gate_point_shape, length(gts))
    } else if (length(gate_point_shape) > length(gts)) {
      gate_point_shape <- gate_point_shape[1:length(gts)]
    }
  }

  # Gate Point Expansion
  if (length(gate_point_size) != length(gts)) {
    if (length(gate_point_size) == 1) {
      gate_point_size <- rep(gate_point_size, length(gts))
    } else if (length(gate_point_size) > length(gts)) {
      gate_point_size <- gate_point_size[1:length(gts)]
    }
  }

  # Plot Gates
  gts <- mapply(function(gt, gate_line_col, gate_line_type, gate_line_width, gate_point_shape, gate_point_size) {
    cyto_plot_gate(gt, channels = channels, gate_line_col = gate_line_col, gate_line_width = gate_line_width, gate_line_type = gate_line_type, gate_point = gate_point, gate_point_shape = gate_point_shape, gate_point_size = gate_point_size)
  }, gts, gate_line_col, gate_line_type, gate_line_width, gate_point_shape, gate_point_size)

  invisible(gts)
})

#' Plot filters List of Gate Objects onto an Existing Plot
#'
#' @param x an object of class \code{\link[flowCore:filters-class]{filters}}
#'   containing objects of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param gate_point logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param gate_point_shape integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param gate_point_size numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' 
#' @return gate object with modified co-ordinates for plotting.
#'
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot_gate,rectangleGate-method}}
#' @seealso \code{\link{cyto_plot_gate,polygonGate-method}}
#' @seealso \code{\link{cyto_plot_gate,ellipsoidGate-method}}
#' @seealso \code{\link{cyto_plot_gate,list-method}}
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
#' # T Cells & Dendritic Cells gate
#' cyto_plot_gate(filters(list(getGate(gs, "T Cells")[[1]], getGate(gs, "Dendritic Cells")[[1]])),
#'   channels = c("APC-Cy7-A", "PE-A"),
#'   gate_line_col = c("green4", "orange")
#' )
#' @export
setMethod(cyto_plot_gate, signature = "filters", definition = function(x,
                                                                       channels,
                                                                       gate_line_type = 1,
                                                                       gate_line_width = 2.5,
                                                                       gate_line_col = "red",
                                                                       gate_point = FALSE,
                                                                       gate_point_shape = 16,
                                                                       gate_point_size = 1) {

  # Assign x to gts
  gts <- x

  # Check Channels
  if (missing(channels)) {
    channels <- unique(as.vector(sapply(gts, parameters)))
  }

  if (!any(unique(as.vector(sapply(gts, parameters))) %in% channels)) {
    stop("Channels used to construct the plot do not match those of the supplied gate.")
  }

  # Gate colours
  if (length(gate_line_col) != length(gts)) {
    if (length(gate_line_col) == 1) {
      gate_line_col <- rep(gate_line_col, length(gts))
    } else if (length(gate_line_col) > length(gts)) {
      gate_line_col <- gate_line_col[1:length(gts)]
    }
  }

  # Gate Line Type
  if (length(gate_line_type) != length(gts)) {
    if (length(gate_line_type) == 1) {
      gate_line_type <- rep(gate_line_type, length(gts))
    } else if (length(gate_line_type) > length(gts)) {
      gate_line_type <- gate_line_type[1:length(gts)]
    }
  }

  # Gate Line Width
  if (length(gate_line_width) != length(gts)) {
    if (length(gate_line_width) == 1) {
      gate_line_width <- rep(gate_line_width, length(gts))
    } else if (length(gate_line_width) > length(gts)) {
      gate_line_width <- gate_line_width[1:length(gts)]
    }
  }

  # Gate Point Character
  if (length(gate_point_shape) != length(gts)) {
    if (length(gate_point_shape) == 1) {
      gate_point_shape <- rep(gate_point_shape, length(gts))
    } else if (length(gate_point_shape) > length(gts)) {
      gate_point_shape <- gate_point_shape[1:length(gts)]
    }
  }

  # Gate Point Expansion
  if (length(gate_point_size) != length(gts)) {
    if (length(gate_point_size) == 1) {
      gate_point_size <- rep(gate_point_size, length(gts))
    } else if (length(gate_point_size) > length(gts)) {
      gate_point_size <- gate_point_size[1:length(gts)]
    }
  }

  # Plot Gates
  gts <- mapply(function(gt, gate_line_col, gate_line_type, gate_line_width, gate_point_shape, gate_point_size) {
    cyto_plot_gate(gt, channels = channels, gate_line_col = gate_line_col, gate_line_width = gate_line_width, gate_line_type = gate_line_type, gate_point = gate_point, gate_point_shape = gate_point_shape, gate_point_size = gate_point_size)
  }, gts, gate_line_col, gate_line_type, gate_line_width, gate_point_shape, gate_point_size)

  invisible(gts)
})
