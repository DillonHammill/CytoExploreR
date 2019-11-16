## INTERNAL GATING FUNCTIONS ---------------------------------------------------

#' Draw Polygon Gate(s) Around Populations.
#'
#' \code{.cyto_gate_polygon_draw} constructs an interactive plotting window to
#' allow manual selection of the co-ordinates of a polygon gate(s) (through
#' mouse click) which are constructed into
#' \code{\link[flowCore:polygonGate-class]{polygonGate}} objects and stored in a
#' \code{\link[flowCore:filters-class]{filters}} list.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3","CD4")}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:polygonGate-class]{polygonGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, polygonGate, openCyto
#'
#' @importFrom flowCore polygonGate filters Subset
#' @importFrom graphics locator lines
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get polygonGate using .cyto_gate_polygon_draw
#' pg <- .cyto_gate_polygon_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("FSC-A", "SSC-A")
#' )
#'
#' # pg is a filters object - extract polygonGate using `[[`
#' pg[[1]]
#' }
#'
#' @noRd
.cyto_gate_polygon_draw <- function(fr,
                                    alias = NULL,
                                    channels,
                                    plot = TRUE,
                                    label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES ------------------------------------------------------------

  # GATES
  gates <- lapply(alias, function(alias) {
    message(
      paste(
        "Select at least 3 points to construct a polygon gate around the",
        alias, "population. \n"
      )
    )

    # GATE CO-ORDINATES
    options("show.error.messages" = FALSE)
    on.exit(options("show.error.messages" = TRUE))
    coords <- locator(
      type = "o",
      lwd = 2,
      pch = 16,
      col = "red"
    )

    # TOO FEW SELECTED POINTS
    if (length(coords$x) < 3) {
      stop("A minimum of 3 points is required to construct a polygon gate.")
    }

    # ADD GATE TO PLOT
    lines(
      x = coords$x[c(1, length(coords$x))],
      y = coords$y[c(1, length(coords$x))],
      lwd = 2.5,
      col = "red"
    )

    # TIDY COORDS
    coords <- as.data.frame(coords)
    coords <- as.matrix(coords)
    colnames(coords) <- channels

    # CONSTRUCT GATE
    gate <- flowCore::polygonGate(.gate = coords, filterId = alias)

    # LABEL GATED POPULATION - MANUAL FOR SPEED
    if (label == TRUE) {
      # GATE CENTER - LABEL POSITION
      gate_center <- .cyto_gate_center(gate,
        channels = channels
      )
      # GATE STAT
      gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
      gate_stat <- paste(.round(gate_stat), "%")
      # PLOT LABEL
      cyto_plot_labeller(
        label_text = paste(alias, gate_stat, sep = "\n"),
        label_text_size = 1,
        label_text_x = gate_center[, "x"],
        label_text_y = gate_center[, "y"]
      )
    }

    return(gate)
  })

  # RETURN CONSTRUCTED GATES
  gates <- filters(gates)
  return(gates)
}

#' Draw Rectangle Gate(s) Around Populations.
#'
#' \code{.cyto_gate_rectangle_draw} constructs an interactive plotting window to
#' allow manual selection of the co-ordinates of a rectangle gate(s) (through
#' mouse click) which are constructed into
#' \code{\link[flowCore:rectangleGate-class]{rectangleGate}} objects and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Simply select 2
#' diagonal co-ordinates to construct the rectangleGate(s).
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate-class]{rectangleGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, rectangleGate, openCyto
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom flowCore rectangleGate filters Subset
#' @importFrom flowCore exprs
#' @importFrom graphics locator rect
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get polygonGate using .cyto_gate_rectangle_draw - add contour lines
#' rg <- .cyto_gate_rectangle_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("FSC-A", "SSC-A"),
#'   contour_lines = 15
#' )
#'
#' # rg is a filters object - extract rectangleGate using `[[`
#' rg[[1]]
#' }
#'
#' @noRd
.cyto_gate_rectangle_draw <- function(fr,
                                      alias = NULL,
                                      channels,
                                      plot = TRUE,
                                      label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = TRUE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES ------------------------------------------------------------

  # GATES
  gates <- lapply(alias, function(alias) {
    message(
      paste(
        "Select 2 diagonal points to construct a rectangle gate around the",
        alias, "population. \n"
      )
    )

    # GATE COORDS
    options("show.error.messages" = FALSE)
    on.exit(options("show.error.messages" = TRUE))
    coords <- locator(
      n = 2,
      type = "p",
      lwd = 2,
      pch = 16,
      col = "red"
    )

    # TIDY GATE COORDS
    coords <- data.frame(coords)
    coords <- as.matrix(coords)
    colnames(coords) <- channels

    # CONSTRUCT GATE
    gate <- flowCore::rectangleGate(.gate = coords, filterId = alias)

    # PLOT GATE
    cyto_plot_gate(gate,
      channels = channels
    )

    # LABEL GATED POPULATION
    if (label == TRUE) {
      # GATE CENTER - LABEL POSITION
      gate_center <- .cyto_gate_center(gate,
        channels = channels
      )
      # GATE STAT
      gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
      gate_stat <- paste(.round(gate_stat), "%")
      # PLOT LABEL
      cyto_plot_labeller(
        label_text = paste(alias, gate_stat, sep = "\n"),
        label_text_size = 1,
        label_text_x = gate_center[, "x"],
        label_text_y = gate_center[, "y"]
      )
    }

    return(gate)
  })

  # RETURN CONSTRUCTED GATES
  gates <- filters(gates)
  return(gates)
}

#' Draw Interval Gate(s) Around Populations.
#'
#' \code{.cyto_gate_interval_draw} constructs an interactive plotting window for
#' user to select the lower and upper bounds of a population (through mouse
#' click) which is constructed into a
#' \code{\link[flowCore:rectangleGate-class]{rectangleGate}} object and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Both 1-D and 2-D
#' interval gates are supported, for 2-D interval gates an additional argument
#' \code{axis} must be supplied to indicate which axis should be gated.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate-class]{rectangleGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, rectangleGate, openCyto, interval
#'
#' @importFrom flowCore rectangleGate filters Subset
#' @importFrom graphics locator abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get 1-D interval gate using .cyto_gate_interval_draw - overlay control
#' ig <- .cyto_gate_interval_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = "PE-A",
#'   overlay = fs[[1]],
#'   density_stack = 0.5
#' )
#'
#' # ig is a filters object - extract rectangleGate using `[[`
#' ig[[1]]
#'
#' # Get 2-D interval gate on y axis using .cyto_gate_interval_draw
#' ig <- .cyto_gate_interval_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("PE-A", "Alexa Fluor 488-A"),
#'   axis = "y"
#' )
#'
#' # ig is a filters object - extract rectangleGate using `[[`
#' ig[[1]]
#' }
#'
#' @noRd
.cyto_gate_interval_draw <- function(fr,
                                     alias = NULL,
                                     channels,
                                     plot = TRUE,
                                     axis = "x",
                                     label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES
  gates <- lapply(alias, function(alias) {
    message(
      paste(
        "Select the lower and upper bounds of the",
        alias, "population to construct an interval gate. \n"
      )
    )

    # GATE COORDS
    options("show.error.messages" = FALSE)
    on.exit(options("show.error.messages" = TRUE))
    coords <- locator(
      n = 2,
      type = "o",
      lwd = 2.5,
      pch = 16,
      col = "red"
    )
    coords <- data.frame(coords)
    coords <- as.matrix(coords)

    if (length(channels) == 1) {
      colnames(coords) <- c(channels[1], "Density")
    } else {
      colnames(coords) <- channels
    }

    if (axis == "x") {
      abline(
        v = coords[, 1],
        lwd = 2.5,
        col = "red"
      )
    } else if (axis == "y") {
      abline(
        h = coords[, 2],
        lwd = 2.5,
        col = "red"
      )
    }

    if (axis == "x") {
      if (length(channels) == 1) {
        coords <- data.frame(x = coords[, 1])
        coords <- as.matrix(coords)
        colnames(coords) <- channels[1]
        rownames(coords) <- c("min", "max")
      } else if (length(channels) == 2) {
        coords <- data.frame(x = coords[, 1], y = c(-Inf, Inf))
        coords <- as.matrix(coords)
        colnames(coords) <- channels
        rownames(coords) <- c("min", "max")
      }
      gate <- rectangleGate(.gate = coords, filterId = alias)

    } else if (axis == "y") {
      if (length(channels) == 1) {
        stop("Cannot gate y axis if a single channel is supplied.")
      }
      coords <- data.frame(x = c(-Inf, Inf), y = coords[, 2])
      coords <- as.matrix(coords)
      colnames(coords) <- channels
      rownames(coords) <- c("min", "max")

      gate <- rectangleGate(.gate = coords, filterId = alias)
    }

    # LABEL GATED POPULATION
    if (label == TRUE) {
      # GATE CENTER - LABEL POSITION
      gate_center <- .cyto_gate_center(gate,
        channels = channels
      )
      # GATE STAT
      gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
      gate_stat <- paste(.round(gate_stat), "%")
      # PLOT LABEL
      cyto_plot_labeller(
        label_text = paste(alias, gate_stat, sep = "\n"),
        label_text_size = 1,
        label_text_x = gate_center[, "x"],
        label_text_y = gate_center[, "y"]
      )
    }

    return(gate)
  })

  # RETURN CONSTRUCTED GATES
  gates <- filters(gates)
  return(gates)
}

#' Draw Threshold Gate(s) Around Populations.
#'
#' \code{.cyto_gate_threshold_draw} constructs an interactive plotting window
#' for user to select the lower bound of a population which is constructed into
#' a \code{\link[flowCore:rectangleGate-class]{rectangleGate}} object and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Both 1-D and 2-D
#' threshold gates are supported, for 2-D threshold gates all events above the
#' select x and y coordinates are included in the gate. Multiple threshold gates
#' are not currently supported.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. Multiple
#'   \code{threshold} gates are not currently supported. \code{alias} is
#'   \code{NULL} by default which will halt the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate-class]{rectangleGate}}
#'   object.
#'
#' @keywords manual, gating, draw, rectangleGate, openCyto, threshold
#'
#' @importFrom flowCore rectangleGate filters Subset
#' @importFrom graphics locator rect abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get 1-D threshold gate using .cyto_gate_threshold_draw
#' tg <- .cyto_gate_threshold_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("PE-A")
#' )
#'
#' # tg is a filters object - extract rectangleGate using `[[`
#' tg[[1]]
#'
#' #' # Get 2-D threshold gate using .cyto_gate_threshold_draw - overlay control
#' tg <- .cyto_gate_threshold_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("Alexa Fluor 647-A", "7-AAD-A"),
#'   overlay = fs[[1]]
#' )
#'
#' # tg is a filters object - extract rectangleGate using `[[`
#' tg[[1]]
#' }
#'
#' @noRd
.cyto_gate_threshold_draw <- function(fr,
                                      alias = NULL,
                                      channels,
                                      plot = TRUE,
                                      label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES ------------------------------------------------------------

  # INSTRUCTIONS
  message(
    paste(
      "Select the lower bound of the",
      alias, "population to construct a threshold gate. \n"
    )
  )

  if (length(alias) > 1) {
    stop("Multiple threshold gates are not supported.")
  }

  # GATE COORDS
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  coords <- locator(
    n = 1,
    type = "p",
    lwd = 2.5,
    pch = 16,
    col = "red"
  )

  # TIDY GATE COORDS
  if (length(channels) == 1) {
    pts <- data.frame(x = c(coords$x, Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    rownames(pts) <- c("min", "max")
  } else if (length(channels) == 2) {
    pts <- data.frame(x = c(coords$x, Inf), y = c(coords$y, Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rownames(pts) <- c("min", "max")
  }

  # CONSTRUCT GATE
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  # PLOT GATE
  cyto_plot_gate(gate = gate,
                 channels = channels)

  # LABEL GATED POPULATION
  if (label == TRUE) {
    # GATE CENTER - LABEL POSITION
    gate_center <- .cyto_gate_center(gate,
      channels = channels
    )
    # GATE STAT
    gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
    gate_stat <- paste(.round(gate_stat), "%")
    # PLOT LABEL
    cyto_plot_labeller(
      label_text = paste(alias, gate_stat, sep = "\n"),
      label_text_size = 1,
      label_text_x = gate_center[, "x"],
      label_text_y = gate_center[, "y"]
    )
  }

  # REURN CONSTRUCTED GATES
  gates <- filters(list(gate))
  return(gates)
}

#' Draw Boundary Gate(s) Around Populations.
#'
#' \code{.cyto_gate_boundary_draw} constructs an interactive plotting window for
#' user to select the upper bound of a population which is constructed into a
#' \code{\link[flowCore:rectangleGate-class]{rectangleGate}} object and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Both 1-D and 2-D
#' boundary gates are supported, for 2-D boundary gates all events below the
#' select x and y coordinates are included in the gate. Multiple boundary gates
#' ares not currently supported.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. Multiple boundary
#'   gates ares not currently supported. \code{alias} is \code{NULL} by default
#'   which will halt the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate-class]{rectangleGate}}
#'   object.
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, boundary
#'
#' @importFrom flowCore rectangleGate filters Subset
#' @importFrom graphics locator rect abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get 1-D boundary gate using .cyto_gate_boundary_draw
#' bg <- .cyto_gate_boundary_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("PE-A")
#' )
#'
#' # bg is a filters object - extract rectangleGate using `[[`
#' bg[[1]]
#'
#' #' # Get 2-D boundary gate using .cyto_gate_boundary_draw
#' tg <- .cyto_gate_boundary_draw(fs[[2]],
#'   alias = "Cells",
#'   channels = c("PE-A", "Alexa Fluor 700-A")
#' )
#'
#' # bg is a filters object - extract rectangleGate using `[[`
#' bg[[1]]
#' }
#'
#' @noRd
.cyto_gate_boundary_draw <- function(fr,
                                     alias = NULL,
                                     channels,
                                     plot = TRUE,
                                     label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES ------------------------------------------------------------

  # INSTRUCTIONS
  message(
    paste(
      "Select the upper bound of the",
      alias, "population to construct a boundary gate. \n"
    )
  )

  if (length(alias) > 1) {
    stop("Multiple boundary gates are not supported.")
  }

  # GATE COORDS
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  coords <- locator(
    n = 1,
    type = "p",
    lwd = 2.5,
    pch = 16,
    col = "red"
  )

  # TIDY GATE COORDS
  if (length(channels) == 1) {
    pts <- data.frame(x = c(-Inf, coords$x))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    rownames(pts) <- c("min", "max")
  } else if (length(channels) == 2) {
    pts <- data.frame(x = c(-Inf, coords$x), y = c(-Inf, coords$y))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rownames(pts) <- c("min", "max")
  }

  # CONSTRUCT GATE
  gate <- rectangleGate(.gate = pts, filterId = alias)

  # PLOT GATE
  cyto_plot_gate(gate = gate,
                 channels = channels)
  
  # LABEL GATED POPULATION
  if (label == TRUE) {
    # GATE CENTER - LABEL POSITION
    gate_center <- .cyto_gate_center(gate,
      channels = channels
    )
    # GATE STAT
    gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
    gate_stat <- paste(.round(gate_stat), "%")
    # PLOT LABEL
    cyto_plot_labeller(
      label_text = paste(alias, gate_stat, sep = "\n"),
      label_text_size = 1,
      label_text_x = gate_center[, "x"],
      label_text_y = gate_center[, "y"]
    )
  }

  # RETURN CONSTRUCTED GATES
  gates <- filters(list(gate))
  return(gates)
}

#' Draw Ellipsoid Gate(s) Around Populations.
#'
#' \code{.cyto_gate_ellipse_draw} constructs an interactive plotting window for
#' user to select the limits of a population in 2 dimensions (4 points) which is
#' constructed into \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}}
#' object and stored in a \code{\link[flowCore:filters-class]{filters}} list.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3","CD4")}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, ellipsoidGate, openCyto, ellipse
#'
#' @importFrom flowCore ellipsoidGate filters Subset
#' @importFrom graphics locator polygon
#' @importFrom methods as
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get ellipsoidGate using .cyto_gate_ellipse_draw
#' eg <- .cyto_gate_ellipse_draw(fs[[4]],
#'   alias = "Cells",
#'   channels = c("PE-A", "Alexa Fluor 700-A"),
#'   overlay = fs[[1]]
#' )
#'
#' # eg is a filters object - extract ellipsoidGate using `[[`
#' eg[[1]]
#' }
#'
#' @noRd
.cyto_gate_ellipse_draw <- function(fr,
                                    alias = NULL,
                                    channels,
                                    plot = TRUE,
                                    label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES ------------------------------------------------------------

  # INSTRUCTIONS
  gates <- lapply(alias, function(alias) {
    message(
      paste(
        "Select 4 points to define the limits of the",
        alias, "population to construct an ellipsoid gate. \n"
      )
    )

    # GATE COORDS
    options("show.error.messages" = FALSE)
    on.exit(options("show.error.messages" = TRUE))
    coords <- locator(
      n = 4,
      type = "p",
      lwd = 2,
      pch = 16,
      col = "red"
    )

    coords <- data.frame(coords)

    # Find which points are on major axis
    dst <- as.matrix(stats::dist(coords))
    mj.pts <- coords[which(dst == max(dst), arr.ind = TRUE)[1, ], ]

    # Find which points are on minor axis
    mr.pts <- coords[!coords$x %in% mj.pts$x & !coords$y %in% mj.pts$y, ]

    # Find center of the major axis
    mj.center <- c(
      (sum(mj.pts$x) / nrow(mj.pts)),
      (sum(mj.pts$y) / nrow(mj.pts))
    )

    # Find center of all points
    center <- c(sum(c(mj.pts$x, mr.pts$x)) / 4, sum(c(mj.pts$y, mr.pts$y)) / 4)

    # Adjust mj.pts to fall on center
    adj <- c((mj.center[1] - center[1]), (mj.center[2] - center[2]))
    mj.pts$x <- mj.pts$x - adj[1]
    mj.pts$y <- mj.pts$y - adj[2]

    # Find major point which lies above center
    max.pt <- mj.pts[mj.pts$y > center[2], ]

    # Radius of the major axis
    a <- stats::dist(mj.pts) / 2

    # Radius of the minor axis
    b <- stats::dist(mr.pts) / 2

    # Angle between horizontal line through center and max.pt
    if (max.pt[1] > center[1]) { # angle < pi/2

      mj.pt.ct <- cbind(max.pt[1], center[2])
      colnames(mj.pt.ct) <- c("x", "y")
      adj <- stats::dist(rbind(center, mj.pt.ct))
      angle <- acos(adj / a)
    } else if (max.pt[1] <= center[1]) { # angle >= pi/2

      mj.pt.ct <- cbind(center[1], max.pt[2])
      colnames(mj.pt.ct) <- c("x", "y")
      opp <- stats::dist(as.matrix(rbind(max.pt, mj.pt.ct)))
      angle <- pi / 2 + asin(opp / a)
    }

    # Covariance matrix
    cinv <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
    cinv[1, 1] <- (((cos(angle) * cos(angle)) /
      (a^2)) + ((sin(angle) * sin(angle)) / (b^2)))
    cinv[2, 1] <- sin(angle) * cos(angle) * ((1 / (a^2)) - (1 / (b^2)))
    cinv[1, 2] <- cinv[2, 1]
    cinv[2, 2] <- (((sin(angle) * sin(angle)) / (a^2)) +
      ((cos(angle) * cos(angle)) / (b^2)))

    cvm <- solve(cinv)

    dimnames(cvm) <- list(channels, channels)

    gate <- ellipsoidGate(
      .gate = cvm,
      mean = center,
      filterId = alias
    )

    # PLOT GATE
    cyto_plot_gate(gate,
      channels = channels
    )

    # LABEL GATED POPULATION
    if (label == TRUE) {
      # GATE CENTER - LABEL POSITION
      gate_center <- .cyto_gate_center(gate,
        channels = channels
      )
      # GATE STAT
      gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
      gate_stat <- paste(.round(gate_stat), "%")
      # PLOT LABEL
      cyto_plot_labeller(
        label_text = paste(alias, gate_stat, sep = "\n"),
        label_text_size = 1,
        label_text_x = gate_center[, "x"],
        label_text_y = gate_center[, "y"]
      )
    }

    return(gate)
  })

  # RETURN FILTERS OBJECT
  gates <- filters(gates)
  return(gates)
}

#' Draw Quadrant Gates Around Populations.
#'
#' \code{.cyto_gate_quadrant_draw} constructs an interactive plotting window for
#' user to select the crosshair center of 4 populations which is used to
#' construct 4 \code{\link[flowCore:rectangleGate-class]{rectangleGate}} objects
#' which are stored in a\code{\link[flowCore:filters-class]{filters}}  list.
#' Populations are assigned in the following order: bottom left, bottom right,
#' top right and top left.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the 4 populations to be gated. \code{alias} is
#'   \code{NULL} by default which will halt the gating routine. \code{alias}
#'   must be supplied right to left and top to bottom (i.e. top right, top left,
#'   bottom right and bottom left).
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:quadGate]{quadGate}}.
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, quadrants
#'
#' @importFrom flowCore quadGate filters split
#' @importFrom graphics locator lines abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get quadrant gates using .cyto_gate_quadrant_draw
#' qg <- .cyto_gate_quadrant_draw(fs[[4]],
#'   alias = c("DN", "CD4", "DP", "CD8"),
#'   channels = c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
#' )
#'
#' # qg is a filters object - extract each rectangleGate using `[[`
#' qg[[4]]
#' }
#'
#' @noRd
.cyto_gate_quadrant_draw <- function(fr,
                                     alias = NULL,
                                     channels,
                                     plot = TRUE,
                                     label = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    alias <- c(
      paste0(channels[1], "+", channels[2], "+"),
      paste0(channels[1], "-", channels[2], "+"),
      paste0(channels[1], "+", channels[2], "-"),
      paste0(channels[1], "-", channels[2], "-")
    )
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  if (!length(alias) == 4) {
    stop("'alias' must contain 4 population names for quadrant gates.")
  }

  # INSTRUCTIONS
  message(
    paste("Select the center point to construct quadrant gates. \n")
  )

  # GATE COORDS
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  pts <- locator(
    n = 1,
    type = "o",
    lwd = 2,
    pch = 16,
    col = "red"
  )

  # CO-ORDINATES MATRIX
  pts <- matrix(unlist(pts), ncol = 2)
  colnames(pts) <- channels

  # QUADGATE CONSTRUCTION
  gate <- quadGate(.gate = pts, filterId = paste(alias, collapse = "|"))

  # PLOT GATE
  cyto_plot_gate(gate, channels = channels)

  # LABEL GATED POPULATION
  if (label == TRUE) {
    # GATE CENTER - LABEL POSITION
    gate_center <- .cyto_gate_center(gate,
      channels = channels
    )
    # GATE STAT
    gate_pops <- .cyto_label_pops(fr, gate)
    gate_stat <- LAPPLY(gate_pops, function(pop) {
      .cyto_count(pop) / .cyto_count(fr) * 100
    })
    gate_stat <- LAPPLY(gate_stat, function(z) {
      paste(.round(z), "%")
    })
    # PLOT LABEL
    cyto_plot_labeller(
      label_text = paste(alias, gate_stat, sep = "\n"),
      label_text_size = 1,
      label_text_x = gate_center[, "x"],
      label_text_y = gate_center[, "y"]
    )
  }

  # RETURN FILTERS OBJECT
  gate <- list(gate)
  return(gate)
}

#' Draw Web Gates Around Populations - EXPERIMENTAL
#'
#' \code{.cyto_gate_web_draw} is a variation of drawQuadrant which allows more
#' flexibility with gate co-ordinates (angled lines) and supports any number of
#' gates as indicated by the \code{alias} argument. To construct the gate simply
#' select the center point and surrounding divider points on plot edge.
#' \code{.cyto_gate_web_draw} will construct the
#' \code{\link[flowCore:polygonGate-class]{polygonGate}} objects and store them
#' in a \code{\link[flowCore:filters-class]{filters}} list.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3","CD4")}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine. Recommended for 3 or more populations.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:polygonGate-class]{polygonGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, polygonGate, openCyto, .cyto_gate_web_draw
#'
#' @importFrom flowCore polygonGate filters Subset
#' @importFrom flowCore exprs
#' @importFrom graphics locator lines par
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw}}
#'
#' @examples
#' \dontrun{
#'
#' # Copy and paste into console to interactively draw gates
#'
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Transform fluorescent channels
#' fs <- transform(fs, estimateLogicle(fs[[4]], cyto_fluor_channels(fs)))
#'
#' # Get web gates using .cyto_gate_web_draw
#' wg <- .cyto_gate_web_draw(fs[[4]],
#'   alias = c("DN", "CD4", "CD8"),
#'   channels = c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
#' )
#'
#' # wg is a filters object - extract each polygonGate using `[[`
#' wg[[4]]
#' }
#'
#' @noRd
.cyto_gate_web_draw <- function(fr,
                                alias = NULL,
                                channels,
                                plot = TRUE,
                                label = TRUE, ...) {

  # WARNING
  message("Web gates are an experimental feature - use at your own risk!")

  # CHECKS ---------------------------------------------------------------------

  # CHANNELS
  channels <- cyto_channels_extract(fr,
    channels = channels,
    plot = TRUE
  )

  # ALIAS
  if (is.null(alias)) {
    stop("Supply a name for the gated population(s) to the 'alias' argument.")
  }

  # CONSTRUCT PLOT -------------------------------------------------------------

  # PLOT
  if (plot == TRUE) {
    cyto_plot(fr,
      channels = channels,
      popup = TRUE,
      legend = FALSE,
      label = FALSE, ...
    )
  }

  # CONSTRUCT GATES ------------------------------------------------------------

  # Select center of the web gate
  message("Select the center of the web gate.")

  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  center <- locator(
    n = 1,
    type = "p",
    lwd = 2,
    pch = 16,
    col = "red"
  )

  # User Prompt
  message("Select surrounding co-ordinates on plot edges to draw a web gate.")

  # Minimum and maximum limits of plot
  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]

  # Get all gate co-ordinates - c(center, others)
  coords <- lapply(seq_len(length(alias)), function(x) {
    options("show.error.messages" = FALSE)
    on.exit(options("show.error.messages" = TRUE))
    pt <- locator(
      n = 1,
      type = "p",
      lwd = 2.5,
      pch = 16,
      col = "red"
    )

    lines(
      x = c(center$x, pt$x),
      y = c(center$y, pt$y),
      lwd = 2.5,
      col = "red"
    )

    return(c(pt$x, pt$y))
  })
  coords <- as.data.frame(do.call(rbind, coords))
  colnames(coords) <- c("x", "y")
  coords <- rbind(center, coords)

  # Determine which quadrants the points are in
  # bottom left anti-clockwise to top right (relative to center)
  quads <- c(0, rep(NA, length(alias)))
  for (i in seq_len(length(coords$x))[-1]) {

    # Bottom left Q1
    if (coords[i, ]$x < center$x & coords[i, ]$y <= center$y) {
      quads[i] <- 1

      # Bottom right Q2
    } else if (coords[i, ]$x >= center$x & coords[i, ]$y < center$y) {
      quads[i] <- 2

      # Top right Q3
    } else if (coords[i, ]$x > center$x & coords[i, ]$y >= center$y) {
      quads[i] <- 3

      # Top left Q4
    } else if (coords[i, ]$x <= center$x & coords[i, ]$y > center$y) {
      quads[i] <- 4
    }
  }
  coords[, "Q"] <- quads
  coords <- coords[with(coords, order(coords$Q)), ]

  # Push points to plot limits (intersection with plot limits)

  # Quadrant 1: find limit intercept and modify point co-ordinates
  if (1 %in% coords$Q) {
    q1 <- coords[coords$Q == 1, ]
    for (x in seq_len(length(q1$Q))) {

      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q1[x, "x"], q1[x, "y"]),
        c(xmin, center$y),
        c(xmin, ymin)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q1[x, "x"], q1[x, "y"]),
        c(center$x, ymin),
        c(xmin, ymin)
      )

      # Check which axis the point should be pushed onto
      if (vint[2] >= ymin) {
        q1[x, c("x", "y")] <- vint
      } else if (vint[2] < ymin) {
        q1[x, c("x", "y")] <- hint
      }
    }
    coords[coords$Q == 1, ] <- q1
  }

  # Quadrant 2: find limit intercept and modify point co-ordinates
  if (2 %in% coords$Q) {
    q2 <- coords[coords$Q == 2, ]
    for (x in seq_len(length(q2$Q))) {

      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q2[x, "x"], q2[x, "y"]),
        c(xmax, center$y),
        c(xmax, ymin)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q2[x, "x"], q2[x, "y"]),
        c(center$x, ymin),
        c(xmax, ymin)
      )

      # Check which axis the point should be pushed onto
      if (vint[2] >= ymin) {
        q2[x, c("x", "y")] <- vint
      } else if (vint[2] < ymin) {
        q2[x, c("x", "y")] <- hint
      }
    }
    coords[coords$Q == 2, ] <- q2
  }

  # Quadrant 3: find limit intercept and modify point co-ordinates
  if (3 %in% coords$Q) {
    q3 <- coords[coords$Q == 3, ]
    for (x in seq_len(length(q3$Q))) {

      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q3[x, "x"], q3[x, "y"]),
        c(xmax, ymax),
        c(xmax, center$y)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q3[x, "x"], q3[x, "y"]),
        c(center$x, ymax),
        c(xmax, ymax)
      )

      # Check which axis the point should be pushed onto
      if (vint[2] >= ymax) {
        q3[x, c("x", "y")] <- hint
      } else if (vint[2] < ymax) {
        q3[x, c("x", "y")] <- vint
      }
    }
    coords[coords$Q == 3, ] <- q3
  }

  # Quadrant 4: find limit intercept and modify point co-ordinates
  if (4 %in% coords$Q) {
    q4 <- coords[coords$Q == 4, ]
    for (x in seq_len(length(q4$Q))) {

      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q4[x, "x"], q4[x, "y"]),
        c(xmin, ymax),
        c(xmin, center$y)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q4[x, "x"], q4[x, "y"]),
        c(xmin, ymax),
        c(center$x, ymax)
      )

      # Check which axis the point should be pushed onto
      if (vint[2] >= ymax) {
        q4[x, c("x", "y")] <- hint
      } else if (vint[2] < ymax) {
        q4[x, c("x", "y")] <- vint
      }
    }
    coords[coords$Q == 4, ] <- q4
  }

  # If multiple points in same quadrant order anticlockwise Q1-Q4
  if (anyDuplicated(coords$Q) != 0) {

    # Quadrant 1
    if (1 %in% coords$Q[duplicated(coords$Q)]) {

      # Multiple points in Q1 - sort by -y then +x
      q1 <- coords[coords$Q == 1, ]
      q1 <- q1[with(q1, order(-q1$y, q1$x)), ]
      coords[coords$Q == 1, c("x", "y")] <- q1[, c("x", "y")]
    }

    # Quadrant 2
    if (2 %in% coords$Q[duplicated(coords$Q)]) {

      # Multiple points in Q2 - sort by +x then +y
      q2 <- coords[coords$Q == 2, ]
      q2 <- q2[with(q2, order(q2$x, q2$y)), ]
      coords[coords$Q == 2, c("x", "y")] <- q2[, c("x", "y")]
    }

    # Quadrant 3
    if (3 %in% coords$Q[duplicated(coords$Q)]) {

      # Multiple points in Q3 - sort by +y then -x
      q3 <- coords[coords$Q == 3, ]
      q3 <- q3[with(q3, order(q3$y, -q3$x)), ]
      coords[coords$Q == 3, c("x", "y")] <- q3[, c("x", "y")]
    }

    # Quadrant 4
    if (4 %in% coords$Q[duplicated(coords$Q)]) {

      # Multiple points in Q4 - sort by -x then -y
      q4 <- coords[coords$Q == 4, ]
      q4 <- q4[with(q4, order(-q4$x, -q4$y)), ]
      coords[coords$Q == 4, c("x", "y")] <- q4[, c("x", "y")]
    }
  }

  # Construct gates using input points
  # Duplicate first point after last point
  coords[(length(coords$Q) + 1), ] <- coords[2, ]
  coords[] <- lapply(coords, round, 4)

  # Gate coordinates using input points
  gates <- list()
  for (i in 2:(length(coords$Q) - 1)) {
    gates[[i - 1]] <- rbind(coords[1, ], coords[i, ], coords[i + 1, ])
  }

  # Check if a corner lies between the points - add as gate co-ordinate
  # Calculate corner points using min & max values
  Q1 <- c(xmin, ymin, 1)
  Q2 <- c(xmax, ymin, 2)
  Q3 <- c(xmax, ymax, 3)
  Q4 <- c(xmin, ymax, 4)
  Q <- matrix(c(Q1, Q2, Q3, Q4), byrow = TRUE, nrow = 4)
  colnames(Q) <- c("x", "y", "Q")
  Q <- data.frame(Q)

  indx <- 1:(length(alias) - 1)

  # Add corners to appropriate gates step-wise
  gates[indx] <- lapply(gates[indx], function(x) {

    # DUPLICATION - points in same quadrant
    if (any(duplicated(x$Q))) {

      # Quadrant 1
      if (1 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "x"] == xmin & x[3, "x"] != xmin) {

          # Include Q1 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, x[3, ])

          # Remove Q1 from Q
          if (1 %in% Q[, "Q"]) {
            Q <<- Q[-match(1, Q[, "Q"]), ]
          }
        }
      }

      # Quadrant 2
      if (2 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "y"] == ymin & x[3, "y"] != ymin) {

          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])

          # Remove Q2 from Q
          if (2 %in% Q[, "Q"]) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        }
      }

      # Quadrant 3
      if (3 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "x"] == xmax & x[3, "x"] != xmax) {

          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])

          # Remove Q3 from Q
          if (3 %in% Q[, "Q"]) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        }
      }

      # Quadrant 4
      if (4 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "y"] == ymax & x[3, "y"] != ymax) {

          # Include Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q4, x[3, ])

          # Remove Q4 from Q
          if (4 %in% Q[, "Q"]) {
            Q <<- Q[-match(4, Q[, "Q"]), ]
          }
        }
      }

      # ADJACENT - points in adjacent quadrants
    } else if (any(x[3, "Q"] - x[2, "Q"] == c(0, 1))) {

      # Q1-Q2
      if (x[2, "Q"] == 1 & x[3, "Q"] == 2) {
        if (x[2, "x"] == xmin & x[3, "x"] == xmax) {

          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, Q2, x[3, ])

          # Remove Q1 and Q2 from Q
          if (any(c(1, 2) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(1, 2), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] == xmin & x[3, "x"] != xmax) {

          # Include Q1 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, x[3, ])

          # Remove Q1 from Q
          if (any(1 %in% Q[, "Q"])) {
            Q <<- Q[-match(1, Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmin & x[3, "x"] == xmax) {

          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])

          # Remove Q2 from Q
          if (any(2 %in% Q[, "Q"])) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        }

        # Q2-Q3
      } else if (x[2, "Q"] == 2 & x[3, "Q"] == 3) {
        if (x[2, "y"] == ymin & x[3, "y"] == ymax) {

          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, x[3, ])

          # Remove Q2 and Q3 from Q
          if (any(c(2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] == ymin & x[3, "y"] != ymax) {

          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])

          # Remove Q2 from Q
          if (any(2 %in% Q[, "Q"])) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        } else if (x[2, "y"] != ymin & x[3, "y"] == ymax) {

          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])

          # Remove Q3 from Q
          if (any(3 %in% Q[, "Q"])) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        }

        # Q3-Q4
      } else if (x[2, "Q"] == 3 & x[3, "Q"] == 4) {
        if (x[2, "x"] == xmax & x[3, "x"] == xmin) {

          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, Q4, x[3, ])

          # Remove Q3 and Q4 from Q
          if (any(c(3, 4) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(3, 4), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] == xmax & x[3, "x"] != xmin) {

          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])

          # Remove Q3 from Q
          if (any(3 %in% Q[, "Q"])) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmax & x[3, "x"] == xmin) {

          # Include Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q4, x[3, ])

          # Remove Q4 from Q
          if (any(4 %in% Q[, "Q"])) {
            Q <<- Q[-match(4, Q[, "Q"]), ]
          }
        }
      }

      # SEPARATED - points separated by a quadrant
    } else if (x[3, "Q"] - x[2, "Q"] == 2) {

      # Q1-Q3
      if (x[2, "Q"] == 1 & x[3, "Q"] == 3) {
        if (x[2, "x"] == xmin & x[3, "y"] == ymax) {

          # Include Q1, Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, Q2, Q3, x[3, ])

          # Remove Q1, Q2 and Q3 from Q
          if (any(c(1, 2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(1, 2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] == xmin & x[3, "y"] != ymax) {

          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, Q2, x[3, ])

          # Remove Q1 and Q2 from Q
          if (any(c(1, 2) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(1, 2), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmin & x[3, "y"] == ymax) {

          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, x[3, ])

          # Remove Q2 and Q3 from Q
          if (any(c(2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmin & x[3, "y"] != ymax) {

          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])

          # Remove Q2 from Q
          if (any(2 %in% Q[, "Q"])) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        }

        # Q2-Q4
      } else if (x[2, "Q"] == 2 & x[3, "Q"] == 4) {
        if (x[2, "y"] == ymin & x[3, "x"] == xmin) {

          # Include Q2, Q3 & Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, Q4, x[3, ])

          # Remove Q2, Q3 and Q4 from Q
          if (any(c(2, 3, 4) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3, 4), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] == ymin & x[3, "x"] != xmin) {

          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, x[3, ])

          # Remove Q2 and Q3 from Q
          if (any(c(2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] != ymin & x[3, "x"] == xmin) {

          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, Q4, x[3, ])

          # Remove Q3 and Q4
          if (any(c(3, 4) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(3, 4), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] != ymin & x[3, "x"] != xmin) {

          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])

          # Remove Q3 from Q
          if (any(3 %in% Q[, "Q"])) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        }
      }
    }

    return(x)
  })

  # Last gate inherits remaining corners
  if (nrow(Q) != 0) {
    if (length(which(Q[, "Q"] >= gates[[length(alias)]][2, "Q"])) != 0) {
      g <- Q[which(Q[, "Q"] >= gates[[length(alias)]][2, "Q"]), ]
    }

    if (length(which(Q[, "Q"] < gates[[length(alias)]][2, "Q"])) != 0) {
      r <- Q[which(Q[, "Q"] < gates[[length(alias)]][2, "Q"]), ]
    }

    if (exists("g") & exists("r")) {
      Q <- rbind(g, r)
    } else if (exists("g") & !exists("r")) {
      Q <- g
    } else if (!exists("g") & exists("r")) {
      Q <- r
    }

    gates[[length(alias)]] <- rbind(
      gates[[length(alias)]][c(1, 2), ],
      Q,
      gates[[length(alias)]][3, ]
    )
  }


  # CONSTRUCT GATES
  gates <- lapply(seq(1, length(gates), 1), function(x) {
    coords <- as.matrix(gates[[x]])[, -3]
    colnames(coords) <- channels
    rownames(coords) <- NULL

    # CONSTRUCT GATE
    gate <- flowCore::polygonGate(.gate = coords, filterId = alias[x])

    # PLOT GATE
    cyto_plot_gate(gate, channels = channels)

    # LABEL GATED POPULATION
    if (label == TRUE) {
      # GATE CENTER - LABEL POSITION
      gate_center <- .cyto_gate_center(gate,
        channels = channels
      )
      # GATE STAT
      gate_stat <- .cyto_count(Subset(fr, gate)) / .cyto_count(fr) * 100
      gate_stat <- paste(.round(gate_stat), "%")
      # PLOT LABEL
      cyto_plot_labeller(
        label_text = paste(alias, gate_stat, sep = "\n"),
        label_text_size = 1,
        label_text_x = gate_center[, "x"],
        label_text_y = gate_center[, "y"]
      )
    }

    return(gate)
  })

  # RETURN CONSTRUCTED GATES
  gates <- filters(gates)
  return(gates)
}
