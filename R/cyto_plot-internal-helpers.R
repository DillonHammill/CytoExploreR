# AXES LIMITS ------------------------------------------------------------------

# This function can only be called on cyto_convert supported objects (i.e.
# flowFrames, flowSets, flowFrame lists and flowSet lists).

#' Get axes limits for cyto_plot
#'
#' Calculates the limits of a valid channel only. Need to calculate y axis
#' limits for density distributions separately in cyto_plot_empty.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:Gatinghierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels name of the channels or markers to be used to construct the
#'   plot.
#' @param parent name of the parental node to extract from GatingHierarchy or
#'   GatingSet.
#' @param overlay a \code{flowFrame}, \code{flowSet}, \code{list of flowFrames},
#'   \code{list of flowSets} or \code{list of flowFrame lists} containing
#'   populations to be overlayed onto the plot(s). Data for overlays will be
#'   merged with \code{x} prior to axis limit calculation to ensure that the
#'   axes limits are set based on all the data to be included in the plot.
#' @param limits indicates whether the limits of the "data" or limits of the
#'   "machine" should be returned. This argument will only influence the upper
#'   limit. The lower limit will always be set to 0, unless the data contains
#'   values below this limit. In such cases the lower limit of the data will be
#'   used instead. This argument is set to "machine" by default.
#' @importFrom flowCore exprs flowSet parameters
#' @importFrom flowWorkspace pData getData
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_axes_limits <- function(x, ...){
  UseMethod(".cyto_plot_axes_text")
}

#' @noRd
.cyto_plot_axes_limits.flowFrame <- function(x,
                                             channels,
                                             overlay = NA,
                                             limits = "machine"){
  
  # overlay should be a flowFrame, flowSet or flowFrame list
  
  # Incorrect limits argument
  if (!limits %in% c("data", "machine", "instrument")) {
    stop("Limits argument should be either 'data' or 'machine'.")
  }
  
  # Set limits to machine if set to instrument - flowCore compatibility
  if (limits == "instrument") {
    limits <- "machine"
  }
  
  # Check channels
  channels <- cyto_channels_extract(x, channels, plot = TRUE)
  
  # Combine x and overlay into flowFrame list
  if(!.all_na(overlay)){
    
    # Convert overlay to flowFrame list
    overlay <- cyto_convert(overlay, "flowFrame list")
    
    # Combine x and overlay into the same flowFrame list
    x <- c(list(x), overlay)
    
  }else{
    
    # Add x to list
    x <- list(x)
    
  }
  
  # Calculate ranges
  data_limits <- .cyto_range(x,
                             channels = channels,
                             limits = limits,
                             plot = TRUE
  )
  
  # Convert to list of x and y limits
  data_limits <- lapply(channels, function(x) {
    data_limits[, x]
  })
  names(data_limits) <- channels
  
  return(data_limits)
  
}

#' overlay either flowFrame, flowSet or flowFrame list - not list of flowFrame
#' lists
#' @noRd
.cyto_plot_axes_limits.flowSet <- function(x,
                                           channels,
                                           overlay = NA,
                                           limits = "machine"){
  
  # Get ranges for x (flowSet)
  data_limits <- .cyto_range(x,
                             channels = channels,
                             limits = limits,
                             plot = TRUE)
  
  # Convert overlay to flowFrame list
  if(!.all_na(overlay)){
    overlay <- cyto_convert(overlay, "flowSet")
    
    data_limits_overlay <- .cyto_range(x,
                                       channels = channels,
                                       limits = limits,
                                       plot = TRUE)
    
    # Replace values in data limits if min is lower or max is higher
    lapply(channels, function(y){
      
      # minimum for overly is less than x
      if(data_limits_overlay[,y][1] < data_limits[,y][1]){
        data_limits[,y][1] <<- data_limits_overlay[,y][1]
      }
      
      # maximum for overlay is greater than x
      if(data_limits_overlay[,y][2] > data_limits[,y][2]){
        data_limits[,y][2] <<- data_limits_overlay[,y][2]
      }
      
    })
    
  }
  
  # Convert to list of x and y limits
  data_limits <- lapply(channels, function(x) {
    data_limits[, x]
  })
  names(data_limits) <- channels
  
  return(data_limits)
  
}

#' overlay should be flowFrame list or flowSet
#' @noRd
.cyto_plot_axes_limits.GatingHierarchy <- function(x,
                                                   parent = "root",
                                                   channels,
                                                   overlay = NA,
                                                   limits = "machine"){
  
  # Convert x to flowFrame
  x <- cyto_convert(x, "flowFrame", parent)
  
  # Convert overlay to flowFrame list
  if(!.all_na(overlay)){
    overlay <- cyto_convert(overlay, "flowSet", parent)
  }
  
  # Make call to flowFrame method
  data_limits <- .cyto_plot_axes_limits(x,
                                        channels,
                                        overlay,
                                        limits)
  
  return(data_limits)
  
}

#' @noRd
.cyto_plot_axes_limits.GatingSet <- function(x,
                                             parent = "root",
                                             channels,
                                             overlay = NA,
                                             limits = "machine"){
  
  # Convert x to flowSet
  x <- cyto_convert(x, "flowSet", parent)
  
  # Convert overlay to flowSet
  if(!.all_na(overlay)){
    overlay <- cyto_convert(overlay, "flowSet", parent)
  }
  
  # Make call to flowSet method
  data_limits <- .cyto_plot_axes_limits(x,
                                        channels,
                                        overlay,
                                        limits)
  
  return(data_limits)
  
}

# AXES TEXT --------------------------------------------------------------------

#' Get Appropriate Axes Labels for Transformed Channels - flowWorkspace
#'
#' @param x object of class \code{flowFrame} or \code{GatingHierarchy}.
#' @param ... additional arguments.
#'
#' @return list containing axis labels and breaks.
#'
#' @noRd
.cyto_plot_axes_text <- function(x, ...) {
  UseMethod(".cyto_plot_axes_text")
}

#' Get Appropriate Axes Labels for Transformed Channels - flowFrame Method
#'
#' @param x an object of class \code{flowFrame}.
#' @param channels name(s) of the channel(s) used to construct the plot.
#' @param axes_trans object of class \code{"transformList"} or
#'   \code{"transformerList"} generated by estimateLogicle containing the
#'   transformations applied to the flowFrame.
#'
#' @return list containing axis labels and breaks.
#'
#' @importFrom flowCore transformList inverseLogicleTransform
#'
#' @noRd
.cyto_plot_axes_text.flowFrame <- function(x,
                                           channels,
                                           axes_trans = NA) {
  
  # Return NA if axes_trans is missing
  if (.all_na(axes_trans)) {
    return(NA)
  } else {

    # axes_trans of incorrect class
    if (!any(inherits(axes_trans, "transformList") |
      inherits(axes_trans, "transformerList"))) {
      stop("Supply a valid transformList/transformerList object to 'axes_trans'.")
    }
  }

  # Convert transformerList to transformList
  if (inherits(axes_trans, "transformerList")) {
    axes_trans <- cyto_transform_convert(axes_trans)
  }

  # Assign x to fr
  fr <- x

  # Get list of axis breaks and labels
  axs <- lapply(channels, function(channel) {

    # Channel not included in axes_trans
    if (!channel %in% names(axes_trans@transforms)) {
      return(NA)
    }

    # Range of values
    r <- as.vector(range(fr)[, channel])

    # Transformation Functions & Breaks
    trans.func <- axes_trans@transforms[[channel]]@f
    inv.func <- inverseLogicleTransform(axes_trans)@transforms[[channel]]@f
    raw <- inv.func(r)
    brks <- .cyto_plot_axes_breaks(raw, n = 5, equal.space = FALSE)


    pos <- signif(trans.func(brks))
    label <- .cyto_plot_axes_inverse(brks, drop.1 = TRUE)

    res <- list(label = label, at = pos)

    return(res)
  })
  names(axs) <- channels

  return(axs)
}

#' Get Appropriate Axes Labels for Transformed Channels - GatingHierarchy Method
#'
#' @param x \code{GatingHiearchy}.
#' @param channels \code{character} name(s) of the channel(S) used to construct
#'   the plot.
#'
#' @return when there is transformation function associated with the given
#'   channel, it returns a list of that contains positions and labels to draw on
#'   the axis otherwise returns NULL.
#'
#' @importFrom flowWorkspace getTransformations getData
#'
#' @noRd
.cyto_plot_axes_text.GatingHierarchy <- function(x,
                                                 channels) {

  # Assign x to gh
  gh <- x

  # Get list of axis breaks and labels
  axs <- lapply(channels, function(channel) {
    res <- gh@axis[[sampleNames(gh)]][[channel]]
    if (is.null(res)) {
      # try to grab trans and do inverse trans for axis label on the fly
      trans <- getTransformations(gh, channel, only.function = FALSE)
      if (is.null(trans)) {
        res <- NA
      } else {
        inv.func <- trans[["inverse"]]
        trans.func <- trans[["transform"]]
        brk.func <- trans[["breaks"]]

        fr <- getData(gh, use.exprs = FALSE)
        r <- as.vector(range(fr)[, channel]) # range
        raw <- inv.func(r)
        brks <- brk.func(raw)
        pos <- signif(trans.func(brks))
        # format it
        label <- trans[["format"]](brks)

        res <- list(label = label, at = pos)
      }
    } else {
      # use the stored axis label if exists
      res$label <- .cyto_plot_axes_inverse(as.numeric(res$label), drop.1 = TRUE)
    }

    return(res)
  })
  names(axs) <- channels

  return(axs)
}

#' Generate the breaks that makes sense for flow data visualization -
#' flowWorkspace
#'
#' @param x the raw data values
#' @param n desired number of breaks (the actual number will be different
#'   depending on the data range)
#' @param equal.space whether breaks at equal-spaced intervals
#' @param trans.fun the transform function (only needed when equal.space is
#'   TRUE)
#' @param inverse.fun the inverse function (only needed when equal.space is
#'   TRUE)
#'
#' @return either 10^n intervals or equal-spaced(after transformed) intervals in
#'   raw scale.
#'
#' @noRd
.cyto_plot_axes_breaks <- function(x,
                                   n = 6,
                                   equal.space = FALSE,
                                   trans.fun, inverse.fun) {
  rng.raw <- range(x, na.rm = TRUE)
  if (equal.space) {
    rng <- trans.fun(rng.raw)
    min <- floor(rng[1])
    max <- ceiling(rng[2])
    if (max == min) {
      return(inverse.fun(min))
    }
    by <- (max - min) / (n - 1)

    myBreaks <- inverse.fun(seq(min, max, by = by))
  } else {
    # log10 (e.g. 0, 10, 1000, ...)
    base10raw <- unlist(lapply(2:n, function(e) 10^e))
    base10raw <- c(0, base10raw)
    myBreaks <- base10raw[base10raw > rng.raw[1] & base10raw < rng.raw[2]]
  }

  myBreaks
}

# copy from sfsmisc/flowWorkspace package
# modified to handle NA values
.cyto_plot_axes_inverse <- function(x, drop.1 = FALSE, digits.fuzz = 7) {
  eT <- floor(log10(abs(x)) + 10^-digits.fuzz)
  mT <- signif(x / 10^eT, digits.fuzz)
  ss <- vector("list", length(x))

  for (i in seq(along = x)) ss[[i]] <- if (is.na(x[i])) {
      quote(NA)
    } else if (x[i] == 0) {
      quote(0)
    } else if (drop.1 && mT[i] == 1) {
      substitute(10^E, list(E = eT[i]))
    } else if (drop.1 && mT[i] == -1) {
      substitute(-10^E, list(E = eT[i]))
    } else {
      substitute(A %*% 10^E, list(A = mT[i], E = eT[i]))
    }

  do.call("expression", ss)
}

# GATES ------------------------------------------------------------------------

#' Get number of gates supplied to cyto_plot
#' @param x gate object, filters object or list object.
#' @return number of gate objects in x.
#' @noRd
.cyto_plot_gate_count <- function(x) {

  # Supported gate objects
  typs <- c("rectangleGate", "polygonGate", "ellipsoidGate", "filters")

  # No gates
  if (is.null(x)) {
    gate_count <- 0
  } else if (class(x) %in% typs) {
    gate_count <- length(x)
  } else if (class(x) == "list") {
    # list of gate objects
    if (all(unlist(lapply(x, "class")) %in% typs)) {
      gate_count <- length(unlist(x))
      # list of lists
    } else if (all(unlist(lapply(x, "class")) == "list")) {
      gate_count <- sum(unlist(lapply(x, function(y) {
        length(y)
      })))
    }
  } else {
    gate_count <- 0
  }
  return(gate_count)
}

# ARGUMENT HANDLERS ------------------------------------------------------------

#' Return a list of valid argument names for cyto_plot
#'
#' This returns a vector argument names currently accepted by
#' .cyto_plot_args_split.
#'
#' @param channels name(s) of channel(s) used to contruct the plot. This is used
#'   internally to filter arguments for .cyto_plot_1d or .cyto_plot_2d.
#'
#' @return vector of accepted .cyto_plot_args_split arguments
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_args <- function(channels) {

  # list of supported arguments
  args <- c(
    "xlab",
    "ylab",
    "title",
    "title_text_font",
    "title_text_size",
    "title_text_col",
    "axes_text_font",
    "axes_text_size",
    "axes_text_col",
    "axes_label_text_font",
    "axes_label_text_size",
    "axes_label_text_col",
    "border_line_type",
    "border_line_width",
    "border_line_col",
    "legend",
    "legend_text",
    "legend_text_font",
    "legend_text_size",
    "legend_text_col",
    "label",
    "label_text",
    "label_stat",
    "label_text_font",
    "label_text_size",
    "label_text_col",
    "label_box_x",
    "label_box_y",
    "label_box_alpha",
    "gate_line_type",
    "gate_line_width",
    "gate_line_col"
  )

  # Arguments only present in 1D plots
  if (length(channels) == 1) {
    args <- c(
      args,
      "density_stack",
      "density_fill",
      "density_fill_alpha",
      "density_line_type",
      "density_line_width",
      "density_line_col",
      "legend_line_col",
      "legend_box_fill"
    )
  }

  # Arguments only present in 2D plots
  if (length(channels) == 2) {
    args <- c(
      args,
      "point_shape",
      "point_size",
      "point_col",
      "point_alpha",
      "point_col_scale",
      "contour_lines",
      "contour_line_type",
      "contour_line_width",
      "contour_line_col",
      "legend_point_col"
    )
  }

  return(args)
}

#' Repeat and split arguments for use in cyto_plot
#'
#' @param x named list of arguments
#' @param channels vector of channels used to construct the plot
#' @param n total number of layers
#' @param plots number of plots
#' @param layers number of layers per plot
#' @param gates number of gates per plot
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_args_split <- function(x,
                                  channels,
                                  n,
                                  plots,
                                  layers,
                                  gates) {

  # Arguments per plot cyto_plot_1d/cyto_plot_2d
  args <- c(
    "xlab",
    "ylab",
    "title",
    "title_text_font",
    "title_text_size",
    "title_text_col",
    "density_stack",
    "axes_text_font",
    "axes_text_size",
    "axes_text_col",
    "axes_label_text_font",
    "axes_label_text_size",
    "axes_label_text_col",
    "legend",
    "legend_text_size",
    "label",
    "border_line_type",
    "border_line_width",
    "border_line_col",
    "contour_lines",
    "contour_line_type",
    "contour_line_width",
    "contour_line_col",
    "point_col_scale"
  )

  lapply(args, function(arg) {
    if (arg %in% names(x)) {
      res <- rep(x[[arg]], length.out = plots)

      if (plots == 1) {
        res <- list(res)
      } else {
        res <- split(res, rep(1:plots, length.out = plots))
      }

      x[[arg]] <<- res
    }
  })

  # Arguments of length 2 per plot
  args <- c("axes_text")

  lapply(args, function(arg) {
    if (arg %in% names(x)) {
      res <- rep(x[[arg]], length.out = plots * 2)

      if (plots == 1) {
        res <- list(res)
      } else {
        res <- split(res, rep(1:plots, length.out = plots * 2, each = 2))
      }

      x[[arg]] <<- res
    }
  })

  # Arguments per layer
  args <- c(
    "density_fill",
    "density_fill_alpha",
    "density_line_type",
    "density_line_width",
    "density_line_col",
    "legend_text",
    "legend_text_font",
    "legend_text_col",
    "legend_line_col",
    "legend_box_fill",
    "legend_point_col",
    "point_shape",
    "point_size",
    "point_col",
    "point_alpha"
  )

  lapply(args, function(arg) {
    if (arg %in% names(x)) {
      res <- rep(x[[arg]], length.out = n)

      if (plots == 1) {
        res <- list(res)
      } else {
        res <- split(res, rep(1:plots, length.out = n, each = layers))
      }

      x[[arg]] <<- res
    }
  })

  # Arguments per gate
  args <- c("gate_line_type", "gate_line_width", "gate_line_col")

  if (gates != 0) {
    lapply(args, function(arg) {
      if (arg %in% names(x)) {
        res <- rep(x[[arg]], length.out = gates * plots)

        if (plots == 1) {
          res <- list(res)
        } else {
          res <- split(res, rep(1:plots,
            length.out = gates * plots,
            each = gates
          ))
        }

        x[[arg]] <<- res
      }
    })
  }

  # cyto_plot_1d
  if (length(channels) == 1) {

    # Arguments per label
    args <- c(
      "label_text",
      "label_stat",
      "label_text_font",
      "label_text_size",
      "label_text_col",
      "label_box_x",
      "label_box_y",
      "label_box_alpha"
    )

    lapply(args, function(arg) {
      if (arg %in% names(x)) {
        if (gates != 0) {
          res <- rep(x[[arg]], length.out = gates * plots * layers)

          if (plots == 1) {
            res <- list(res)
          } else {
            res <- split(res, rep(1:plots,
              length_out = gates * plots * layers,
              each = gates * layers
            ))
          }

          x[[arg]] <<- res

          # labels without gates
        } else {
          res <- rep(x[[arg]], length.out = n)

          if (plots == 1) {
            res <- list(res)
          } else {
            res <- split(res, rep(1:plots,
              length_out = n,
              each = layers
            ))
          }

          x[[arg]] <<- res
        }
      }
    })

    # cyto_plot_2d
  } else if (length(channels) == 2) {

    # Arguments per gate
    args <- c(
      "label_text",
      "label_stat",
      "label_text_font",
      "label_text_size",
      "label_text_col",
      "label_box_x",
      "label_box_y",
      "label_box_alpha"
    )

    if (gates != 0) {
      lapply(args, function(arg) {
        if (arg %in% names(x)) {
          res <- rep(x[[arg]],
            length.out = gates * plots
          )

          if (plots == 1) {
            res <- list(res)
          } else {
            res <- split(res, rep(1:plots,
              length.out = gates * plots,
              each = gates
            ))
          }

          x[[arg]] <<- res
        }
      })

      # No gates expect 1 label per layer
    } else if (gates == 0) {
      lapply(args, function(arg) {
        if (arg %in% names(x)) {
          res <- rep(x[[arg]], length.out = n)

          if (plots == 1) {
            res <- list(res)
          } else {
            res <- split(res, rep(1:plots,
              length_out = n,
              each = layers
            ))
          }

          x[[arg]] <<- res
        }
      })
    }
  }

  return(x)
}

# POINT DENSITY COLOURS --------------------------------------------------------

#' Get density gradient colours for cyto_plot_2d
#'
#' @param point_cols vector of ordered colours to use for point density colour
#'   scale.
#'
#' @return a list of colorRampPalette functions to be used in densCols.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom grDevices colorRampPalette
#'
#' @noRd
.cyto_plot_point_cols <- function(point_cols = NA) {

  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Use default colour scale
  if (.all_na(args[["point_cols"]])) {
    args[["point_cols"]] <- c(
      "blue",
      "turquoise",
      "green",
      "yellow",
      "orange",
      "red",
      "darkred"
    )
  }

  return(args[["point_cols"]])
}

# GATES ------------------------------------------------------------------------

#' Check gate object supplied to cyto_plot
#'
#' @param gate gate object(s) to add to plot.
#' @param smp number of samples
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_gate_check <- function(gate, smp) {

  # A single gate supplied
  if (inherits(gate, "rectangleGate") |
    inherits(gate, "polygonGate") |
    inherits(gate, "ellipsoidGate")) {
    gate <- rep(list(list(gate)), smp)

    # A filters object supplied
  } else if (inherits(gate, "filters")) {
    gate <- rep(list(filters), smp)

    # A list of gate objects
  } else if (inherits(gate, "list")) {

    # List of individual gates
    if (all(unlist(lapply(gate, class)) %in%
      c("rectangleGate", "polygonGate", "ellipsoidGate"))) {

      # Assume 1 gate per sample
      if (length(gate) == smp) {
        gate <- lapply(gate, "list")

        # Assume list contains gates for each sample
      } else {
        gate <- rep(list(gate), smp)
      }

      # List of filters objects
    } else if (all(unlist(lapply(gate, class)) == "filters")) {

      # Assume 1 filters object per sample
      gate <- rep(gate, length.out = smp)
    }

    # filtersList
  } else if (inherits(gate, "filtersList")) {
    gate <- gate
  }

  return(gate)
}

# OVERLAY CONVERT --------------------------------------------------------------

#' Format overlay for use in cyto_plot
#'
#' \code{.cyto_plot_overlay_convert} will check whether the supplied overlay is
#' supported and convert it into an appropriate format for use in
#' \code{\link{cyto_plot}}.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param overlay object to overlay.
#' 
#' @return converted overlay
#' 
#' @importFrom flowCore fsApply
#' @importFrom flowWorkspace getData getNodes
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
cyto_plot_overlay_convert <- function(x, ...){
  UseMethod(".cyto_plot_overlay_convert")
}

#' @noRd
.cyto_plot_overlay_convert.flowFrame <- function(x, 
                                                 overlay = NA){
  
  # No overlay supplied
  if(.all_na(overlay)){
    return(overlay)
  }
  
  # Convert overlay to list of flowFrames
  
  # flowFrame -> flowFrame list
  if (inherits(x, "flowFrame")) {
      overlay <- cyto_convert(overlay,
                              "flowFrame list")
  # flowSet -> list of flowFrames
  } else if (inherits(overlay, "flowSet")) {
      overlay <- cyto_convert(overlay,
                              "list of flowFrames")
  # list of flowFrames
  } else if (all(unlist(lapply(overlay, class)) == "flowFrame")) {
    
  # flowSet list -> list of flowFrames
  } else if (all(unlist(lapply(overlay, function(x) {
    inherits(x, "flowSet")
  }))) &
  length(overlay) == 1) {
    
    overlay <- cyto_convert(overlay[[1]],
                            "list of flowFrames")
  } else {
    stop(paste("'overlay' must be a flowFrame, flowSet,",
               "list of flowFrames or a list containing a flowSet.",
               sep = " "
    ))
  }
  
  # return is a list of flowFrames to overlay
  return(overlay)
}

#' @noRd
.cyto_plot_overlay_convert.flowSet <- function(x,
                                               overlay = NA){
  
  # No overlay supplied
  if(.all_na(overlay)){
    return(overlay)
  }
  
  # Convert overlay to list of lists of flowFrames (1 list per element in x)
  
  # flowFrame
  if (class(overlay) == "flowFrame") {
    overlay <- lapply(rep(list(overlay), length(x)), "list")
  # flowSet
  } else if (inherits(overlay, "flowSet")) {
    overlay <- lapply(lapply(
      seq(1, length(overlay), 1),
      function(z) overlay[[z]]
    ), "list")
  # list of flowFrames
  } else if (all(unlist(lapply(overlay, class)) == "flowFrame")) {
    if (length(overlay) == 1) {
      overlay <- lapply(rep(list(overlay[[1]]), length(x)), "list")
    } else {
      if (length(overlay) != length(x)) {
        stop(paste("Supplied list of flowFrames must be of the",
                   "same length as the flowSet.",
                   sep = " "
        ))
      }
      overlay <- lapply(overlay, "list")
    }
  # list of flowSets
  } else if (all(unlist(lapply(overlay, function(x) {
    inherits(x, "flowSet")
  })))) {
    if (!all(unlist(lapply(overlay, length)) == length(x))) {
      stop(paste("Each flowSet in supplied list should be of the",
                 "same length as the supplied flowSet.",
                 sep = " "
      ))
    }
    
    # list of flowFrame lists
    overlay <- lapply(overlay, function(z) {
      lapply(seq(1, length(z), 1), function(y) z[[y]])
    })
    overlay <- lapply(seq_along(x), function(z) {
      lapply(overlay, `[[`, z)
    })
  } else if (all(do.call("rbind", lapply(overlay, function(z) {
    lapply(z, "class")
  })) == "flowFrame")) {
    if (length(overlay) != length(x)) {
      stop(
        paste("'overlay' should be a list of flowFrames lists",
              "to overlay on each flowFrame in the flowSet.",
              sep = " "
        )
      )
    }
  } else {
    stop(paste("'overlay' must be a flowFrame, flowSet,",
               "list of flowFrames or a list of flowSets.",
               sep = " "
    ))
  }
  
  # return is a list of flowFrame lists to overlay
  # 1 flowFrame list per flowFrame in x
  return(overlay)
}

#' @noRd
.cyto_plot_overlay_convert_GatingHierarchy <- function(x, 
                                                       overlay = NA){
  
  # Overlay could be character vector of population names
  if (inherits(overlay, "flowFrame") |
      inherits(overlay, "flowSet") |
      inherits(overlay, "list")) {
    
  } else if (inherits(overlay, "character")) {
    if (!all(overlay %in% basename(getNodes(x)))) {
      stop("'overlay' does not exist in the GatingHierarchy.")
    } else {
      nms <- overlay
      overlay <- lapply(overlay, function(z) {
        getData(x, z)
      })
      names(overlay) <- nms
    }
  }
  
  # .cyto_plot_overlay_convert to convert overlay to correct format
  .cyto_plot_overlay_convert(getData(x, "root"),
                             overlay = overlay)
  
}

#' @noRd
.cyto_plot_overlay_convert.GatingSet <- function(x, 
                                                 overlay = NA){
    # Overlay could be character vector of population names
    if (inherits(overlay, "flowFrame") |
      inherits(overlay, "flowSet") |
      inherits(overlay, "list")) {

    } else if (inherits(overlay, "character")) {
      if (!all(overlay %in% basename(getNodes(x)))) {
        stop("overlay' does not exist in the GatingHierarchy.")
      } else {
        nms <- overlay
        overlay <- lapply(overlay, function(z) {
          getData(x, z)
        })
        names(overlay) <- nms
      }
    }

    # .cyto_plot_overlay_convert to convert overlay to correct format
    .cyto_plot_overlay_convert(getData(x, "root"),
      overlay = overlay)
    
}

# LAYOUT -----------------------------------------------------------------------

#' Set plot layout
#'
#' @param x object to be plotted.
#' @param layout grid dimensions c(nr, nc), NULL or FALSE.
#' @param density_stack degree of offset.
#' @param denisity_layers number of layers per plot.
#'
#' @importFrom grDevices n2mfrow
#'
#' @noRd
.cyto_plot_layout <- function(x,
                              layout = NA,
                              density_stack = 0,
                              density_layers = 1) {

  # Number of samples
  smp <- length(x)

  # Stacking
  if (density_stack != 0) {
    if (density_layers == smp) {
      smp <- ceiling(smp / smp)
    } else {
      smp <- ceiling(smp / density_layers)
    }
  }

  # Plot layout
  if (.all_na(layout)) {
    if (smp > 1) {
      mfrw <- c(grDevices::n2mfrow(smp)[2], grDevices::n2mfrow(smp)[1])
    } else {
      mfrw <- c(1, 1)
    }
  } else if (!.all_na(layout)) {
    if (layout[1] == FALSE) {

      # Do nothing
    } else {
      mfrw <- layout
    }
  }

  return(mfrw)
}

# MARGINS ----------------------------------------------------------------------

#' Set plot margins
#'
#' @param x flowFrame or flowSet object to be plotted (post merging).
#' @param overlay object to overlay.
#' @param legend logical indicating whether a legend should be included in the
#'   plot.
#' @param legend_text text to be used in the legend, used to calculate required
#'   space.
#' @param title if NULL remove excess space above plot.
#' @param axes_text vector of logicals indicating whether the x and y axes
#'   should be included on the plot.
#'
#' @noRd
.cyto_plot_margins <- function(x,
                               overlay = NA,
                               legend = FALSE,
                               legend_text = NA,
                               title,
                               axes_text = c(TRUE,TRUE)) {

  # Bypass setting margins on cyto_plot_grid
  if (!getOption("CytoRSuite_cyto_plot_grid")) {
      
    # Pull down arguments to named lis
    args <- .args_list()
    
    # Default margins
    mar <- c(5.1,5.1,4.1,2.1)
      
    # Make space for legend text on right
    if(!.all_na(args[["overlay"]]) & 
       legend != FALSE & 
       !.all_na(legend_text)){
      mar[4] <- 7 + max(nchar(legend_text)) * 0.32
    }
      
    # Remove space above plot if no title
    if(.all_na(title)){
      mar[3] <- 2.1
    }
      
    # Remove space below plot if x axis is missing
    if(!axes_text[1]){
      mar[1] <- 4.1
    }
      
    # Remove space below plot if y axis is missing
    if(!axes_text[2]){
      mar[2] <- 4.1
    }
    
    # Set update graphics parameter
    par("mar" = mar)
  }
}

# LEGEND -----------------------------------------------------------------------

#' Create a legend for cyto_plot
#'
#' \code{.cyto_plot_margins} will handle setting the plot margins to make space
#' for the legend.
#'
#' @param channels name of the channels or markers to be used to construct the
#'   plot.
#' @param legend logical indicating whether a legend should be included for
#'   plots including overlays, set to FALSE by default.
#' @param legend_text vector of labels to use for the legend.
#' @param legend_text_font numeric indicating the font to use for legend text,
#'   set to 2 for bold font by default. See \code{\link[graphics:par]{?par}}
#'   font for details.
#' @param legend_text_size character expansion for legend text, set to 1 by
#'   default.
#' @param legend_text_col colour to use for legend text, set to "black by
#'   default.
#' @param legend_line_col vector of line colours to use for legend.
#' @param legend_box_fill vector of fill colours to use for legend.
#' @param legend_point_col vector of colours to use for points in legend.
#' @param density_fill colour(s) used to fill polygons.
#' @param density_fill_alpha numeric [0,1] used to control fill transparency,
#'   set to 1 by default to remove transparency.
#' @param density_line_type line type(s) to use for border(s), set to solid
#'   lines by default.
#' @param density_line_width line width for border.
#' @param density_line_col colour(s) for border line, set to "black" by default.
#' @param point_shape point character to use for points, set to "." by default
#'   to maximise plotting speed.
#' @param point_size numeric specifying the degree of character expansion for
#'   points, set to 2 by default.
#' @param point_col colours to use for points, set to NA by default to blue-red
#'   density colour scale.
#' @param point_alpha numeric [0,1] used to control colour transparency, set to
#'   1 by default to remove transparency.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom graphics legend
#'
#' @noRd
.cyto_plot_legend <- function(channels,
                              legend = "fill",
                              legend_text = NA,
                              legend_text_font = 1,
                              legend_text_size = 1,
                              legend_text_col = "black",
                              legend_line_col = NA,
                              legend_box_fill = NA,
                              legend_point_col = NA,
                              density_fill = NA,
                              density_fill_alpha = 1,
                              density_line_type = 1,
                              density_line_width = 1,
                              density_line_col = "black",
                              point_shape = ".",
                              point_size = 2,
                              point_col = NA,
                              point_col_alpha = 1) {
  
  # Legend for 1D density distributions
  if (length(channels) == 1) {

    # Set default legend type to fill
    if (legend == TRUE) {
      legend <- "fill"
    }
    
    # Reverse legend text order for legend
    legend_text <- rev(legend_text)

    # Line legend
    if (legend == "line") {

      # Revert to density_line_col if no colours supplied
      if (.all_na(legend_line_col)) {
        legend_line_col <- density_line_col
      }

      # Construct legend
      legend(
        x = "right",
        inset = c(-0.68,0),
        legend = legend_text,
        text_font = rev(legend_text_font),
        cex = legend_text_size,
        text.col = rev(legend_text_col),
        col = rev(legend_line_col),
        lty = rev(legend_line_type),
        lwd = rev(legend_line_width),
        xpd = TRUE,
        bty = "n",
        x.intersp = 0.5
      )
    } else if (legend == "fill") {

      # Revert to density_fill if no legend fill colours supplied
      if (.all_na(legend_box_fill)) {
        legend_box_fill <- density_fill
        
      # Alpha adjust colours if suppplied directly to legend_box_fill
      }else if (!.all_na(legend_box_fill) &
                !all(density_fill_alpha == 1)) {
        legend_box_fill <- mapply(
          function(legend_box_fill,
                             density_fill_alpha) {
            adjustcolor(legend_box_fill, density_fill_alpha)
          }, legend_box_fill, density_fill_alpha
        )
      }

      # Construct legend
      legend(
        x = "right", # right inside plot
        inset = c(-0.68,0), # move outside 0.47 graphics device widths
        legend = legend_text,
        fill = rev(legend_box_fill),
        xpd = TRUE,
        bty = "n",
        x.intersp = 0.5,
        cex = legend_text_size,
        text.col = rev(legend_text_col),
        text.font = rev(legend_text_font)
      )
    }

    # Legend for 2D scatter plot
  } else if (length(channels) == 2) {

    # 
    
    # Revert to point_col if no legen point cols supplied
    if(.all_na(legend_point_col)){
      legend_point_col <- point_col
    # Alpha adjust colours supplied directly to legend_point_col
    }else if (!.all_na(legend_point_col) &
              !all(point_col_alpha == 1)) {
      legend_point_col <- mapply(function(col, alpha) {
        adjustcolor(col, alpha)
      }, legend_point_col, point_col_alpha)
    }

    legend(
      x = "right",
      inset = c(-0.68,0),
      legend = rev(legend_text),
      col = rev(legend_point_col),
      pch = rev(point_shape),
      pt.cex = rev(2 * point_size),
      xpd = TRUE,
      bty = "n",
      x.intersp = 0.5
    )
  }
}

# THEME INHERIT ----------------------------------------------------------------

#' Inherit cyto_plot_theme arguments
#'
#' @param x list of named cyto_plot arguments.
#'
#' @return updated list of named arguments if cyto_plot_theme has been set.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_theme_inherit <- function(x) {
  
  # extract cyto_plot_theme arguments
  args <- getOption("CytoRSuite_cyto_plot_theme")
  
  if (!is.null(args)) {
    lapply(names(args), function(y){
      x[[y]] <<- args[[y]]
    })
  }

  return(x)
}

# TITLE ------------------------------------------------------------------------

#' Title for cyto_plot
#'
#' @param x flowFrame object.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_title <- function(x,
                             channels,
                             overlay = NA,
                             title) {
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # 1D density distributions
  if (length(args[["channels"]]) == 1) {

    # missing/empty replace with valid title
    if(.empty(args[["title"]])){
    
      # stacked/overlays lack a title
      if(.all_na(args[["overlay"]])){
        args[["title"]] <- identifier(args[["x"]])
        if(args[["title"]] == "anonymous"){
          args[["title"]] <- "Combined Events"
        }
      }else{
        args[["title"]] <- NA
      }
      
    # NA will remove title in cyto_plot_empty  
    }else if(.all_na(args[["title"]])){
      args[["title"]] <- NA
    }

    # 2D scatterplots
  } else if (length(args[["channels"]]) == 2) {

    # missing title replaced with sample name
    if(.empty(args[["title"]])){
      args[["title"]] <- identifier(args[["x"]])
      if(args[["title"]] == "anonymous"){
        args[["title"]] <- "Combined Events"
      }
      
    # NA will remove title in cyto_plot_empty  
    }else if(.all_na(args[["title"]])){
      args[["title"]] <- NA
    }
    
  }

  return(args[["title"]])
}

# AXES LABELS ------------------------------------------------------------------

#' Get axes titles for cyto_plot
#'
#' @param x flowFrame object.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_axes_label <- function(x,
                                  channels,
                                  xlab,
                                  ylab,
                                  density_modal = TRUE) {

  # Extract information about channels
  fr_data <- pData(parameters(x))
  fr_channels <- BiocGenerics::colnames(x)

  # 1D density distributions
  if (length(channels) == 1) {

    # x axis label
    if(missing(xlab) | .empty(xlab)){
      # Marker assigned to channel
      if (!is.na(fr_data$desc[which(fr_channels == channels)])) {
        xlab <- paste(fr_data$desc[which(fr_channels == channels)],
                      channels,
                      sep = " "
        )
        # No assigned marker to channel
      } else if (is.na(fr_data$desc[which(fr_channels == channels)])) {
        xlab <- paste(channels)
      }
    }else if(.all_na(xlab)){
      xlab <- NA
    }
    
    # y axis label
    if(missing(ylab) | .empty(ylab)){
      if (density_modal) {
        ylab <- "Density Normalised to Mode (%)"
      } else {
        ylab <- "Density"
      }
    }else if(.all_na(ylab)){
      ylab <- NA
    }
    
    # 2D scatterplots
  } else if (length(channels) == 2) {

    # x axis label
    if(missing(xlab) | .empty(xlab)){
      # Marker assigned to channel
      if (!is.na(fr_data$desc[which(fr_channels == channels[1])])) {
        xlab <- paste(fr_data$desc[which(fr_channels == channels[1])],
                      channels[1],
                      sep = " "
        )
        # No assigned marker to channel
      } else if (is.na(fr_data$desc[which(fr_channels == channels[1])])) {
        xlab <- paste(channels[1])
      }
    }else if(.all_na(xlab)){
      xlab <- NA
    }
    
    # y axis label
    if(missing(ylab) | .empty(ylab)){
      # Marker assigned to channel
      if (!is.na(fr_data$desc[which(fr_channels == channels[2])])) {
        ylab <- paste(fr_data$desc[which(fr_channels == channels[2])],
                      channels[2],
                      sep = " "
        )
        # No assigned marker to channel
      } else if (is.na(fr_data$desc[which(fr_channels == channels[2])])) {
        ylab <- paste(channels[2])
      }
    }else if(.all_na(ylab)){
      ylab <- NA
    }
    
  }

  return(list(xlab, ylab))
}

# DENSITY FILL -----------------------------------------------------------------

#' Get density fill colours for cyto_plot
#' 
#' @param x list of flowFrame or density objects.
#' @param density_fill vector of colours to use for each layer.
#' @param density_cols vector of colls to use to select density_fill colours.
#' 
#' @importFrom grDevices adjustcolor colorRampPalette
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
.cyto_plot_density_fill <- function(x,
                                    density_fill = NA,
                                    density_cols = NA,
                                    density_fill_alpha = 1) {
  
  # Expected number of colours
  n <- length(x)

  # Pull down arguments to named list
  args <- as.list(environment())
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # No density_cols supplied
  if(.all_na(args[["density_cols"]])){
    args[["density_cols"]] <- c("grey",
                                "bisque4",
                                "brown1",
                                "red",
                                "darkred",
                                "chocolate",
                                "orange",
                                "yellow",
                                "yellowgreen",
                                "green",
                                "aquamarine",
                                "cyan",
                                "cornflowerblue",
                                "blue",
                                "blueviolet",
                                "purple",
                                "magenta",
                                "deeppink")
  }

  
  # Make colorRampPalette
  if(class(args[["density_cols"]]) != "function"){
    cols <- colorRampPalette(args[["density_cols"]])
  }else{
    cols <- args[["density_cols"]]
  }

  # No colours supplied to density_fill either
  if(.all_na(args[["density_fill"]])){
      
    # Pull out a single colour per layer
    args[["density_fill"]] <- cols(n)
    
  # Colours supplied manually to density_fill 
  }else{
      
    # Too few colours supplied - pull others from cols
    if(length(args[["density_fill"]]) < n){
        
      args[["density_fill"]] <- c(args[["density_fill"]],
                                    cols(n - length(args[["density_fill"]])))
        
    # Too many colours supplied
    }else if(length(args[["density_fill"]]) > n){
        
      args[["density_fill"]] <- args[["density_fill"]][seq(1,n)]
        
    }
    
  }
  
  # Adjust colors by density_fill_alpha
  args[["density_fill"]] <- mapply(function(density_fill, density_fill_alpha){
    adjustcolor(density_fill, density_fill_alpha)
  }, args[["density_fill"]], args[["density_fill_alpha"]], USE.NAMES = FALSE)
  
  return(args[["density_fill"]])
  
}

# POINT COLOUR -----------------------------------------------------------------

#' Get point colours for cyto_plot
#'
#' @param x list of flowFrames.
#' @param channels used to construct the plot.
#' @param point_col_scale vector of colours to use for density gradient.
#' @param point_cols vector colours to select from when choosing a colour for
#'   each layer in x.
#' @param point_col vector of length x indicating colours to use for each layer.
#'   If NA set to default density gradient.
#' @param point_col_alpha transparency to use for point colours.
#'
#' @importFrom grDevices densCols colorRampPalette adjustcolor
#' @importFrom flowCore exprs
#' 
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_point_col <- function(x,
                                 channels,
                                 point_col_scale,
                                 point_cols,
                                 point_col,
                                 point_col_alpha){
  
  # Expected number of colours
  n <- length(x)
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # No colours supplied for density gradient
  if(.all_na(args[["point_col_scale"]])){
    args[["point_col_scale"]] <- c("blue",
                                   "turquoise",
                                   "green",
                                   "yellow",
                                   "orange",
                                   "red",
                                   "darkred")
  }
  
  # Make colorRampPalette
  if(class(args[["point_col_scale"]]) != "function"){
    col_scale <- colorRampPalette(args[["point_col_scale"]])
  }else{
    cols_scale <- args[["point_col_scale"]]
  }
  
  # No colours supplied for selection
  if(.all_na(args[["point_cols"]])){
    args[["point_cols"]] <- c("black",
                              "darkorchid",
                              "blueviolet",
                              "magenta",
                              "deeppink",
                              "red4",
                              "orange",
                              "springgreen4")
  }
  
  # Make colorRampPalette
  if(class(args[["point_cols"]]) != "function"){
    cols <- colorRampPalette(args[["point_cols"]])
  }else{
    cols <- args[["point_cols"]]
  }
  
  # Repeat point_col length of x and convert to list
  args[["point_col"]] <- rep(args[["point_col"]], 
                             length.out = length(x))
  args[["point_col"]] <- lapply(seq(1,length(x)), function(z){
                                args[["point_col"]][z]})
  
  # First layer contains density gradient in no other colour is designated
  if(.all_na(args[["point_col"]])){
    
    # Extract data
    fr_exprs <- exprs(x[[1]])[, channels]
    
    # Too few events for density computation
    if(nrow(fr_exprs) <= 2){
      # Get density colour for each point
      args[["point_col"]][[1]] <- densCols(fr_exprs, 
                                       colramp = col_scale)
    }else{
      args[["point_col"]][[1]] <- args[["point_col_scale"]][1]
    }

  }else if(any(is.na(args[["point_col"]]))){
    
    # Get position of NAs in point_col - replace with density scale
    ind <- which(is.na(args[["point_col"]]))
    
    # Run through each ind and get density gradient colours
    lapply(ind, function(z){
      
      # Extract data
      fr_exprs <- exprs(x[[z]])[,channels]
      
      # Too few events for density computation
      if(nrow(fr_exprs) <= 2){
        # Get density colour for each point
        args[["point_col"]][[z]] <<- densCols(fr_exprs,
                                              colramp = col_scale)
      }else{
        args[["point_col"]][[z]] <- args[["point_col_scale"]][1]
      }

    })
    
  }
  
  # Remaining colours are selected one per layer from point_cols
  if(any(is.na(args[["point_col"]]))){
    
    # Number of layers missing colours
    n <- length(args[["point_col"]][is.na(args[["point_col"]])])
    
    # Pull colours out of point_cols
    clrs <- cols(n) 
    
    # Replace NA values in point_col with selected colours
    args[["point_col"]][is.na(args[["point_col"]])] <- clrs
    
  }
  
  # Adjust colors by point_fill_alpha
  args[["point_col"]] <- mapply(function(point_col, point_col_alpha){
    adjustcolor(point_col, point_col_alpha)
  }, args[["point_col"]], args[["point_col_alpha"]], USE.NAMES = FALSE)
  
  return(args[["point_col"]])

}

# GATE STACKED DENSITY ---------------------------------------------------------

#' Gate 1D with overlays
#'
#' @param x flowFrame (base).
#' @param channels used in the plot.
#' @param overlay list of flowFrames to overlay.
#' @param gate gate object(s).
#' @param trans transform object used by cyto_plot_label to calculate
#'   statistics.
#' @param density_stack degree of stacking.
#' @param density_modal logical indicating whether y axis should be normalised
#'   to mode.
#' @param label_text text to use in label.
#' @param label_stat statistic to use in label.
#' @param gate_line_col gate(s) colour(s).
#' @param gate_line_width gate(s) line width(s).
#' @param gate_line_type gate(s) line type(s).
#' @param label logical indicating whether labels should be added.
#' @param label_text_font font(s) for labels.
#' @param label_text_size text size(s) for labels.
#' @param label_text_col text colour(s) for labels.
#' @param label_box_x x co-ordinate(s) for label(s).
#' @param label_box_y y co-ordinates for label(s).
#' @param label_box_alpha transparency for label(s).
#'
#' @importFrom flowCore parameters
#' @importFrom graphics par
#'
#' @noRd
.cyto_plot_overlay_gate <- function(x,
                                    channels,
                                    overlay = NA,
                                    gate = NA,
                                    trans = NA,
                                    density_stack = 0.6,
                                    density_modal = FALSE,
                                    label_text = NA,
                                    label_stat = "median",
                                    gate_line_col = "red",
                                    gate_line_width = 2.5,
                                    gate_line_type = 1,
                                    label = TRUE,
                                    label_text_font = 2,
                                    label_text_size = 0.8,
                                    label_text_col = "black",
                                    label_box_x = NA,
                                    label_box_y = NA,
                                    label_box_alpha = 0.6, ...) {
  
  # Add support for adding labels without gates
  
  # Check class of x
  if (!inherits(x, "flowFrame")) {
    stop("x should be a flowFrame object.")
  }
  
  # Samples
  smp <- length(overlay) + 1
  
  # Check channels
  channels <- cyto_channels_extract(
    x = x,
    channels = channels,
    plot = TRUE
  )
  
  # List of flowFrames for cyto_plot_label
  fr.lst <- c(list(x), overlay)
  
  # Gates supplied
  if(!.all_na(gate)){
  
    # list of gate
    if (inherits(gate, "filters")) {
      
      # Convert to list of gate
      gate <- lapply(seq_len(length(gate)), function(x) gate[[x]])
      
    } else if (inherits(gate, "list")) {
      
    } else if (inherits(gate, "rectangleGate") |
               inherits(gate, "polygonGate") |
               inherits(gate, "ellipsoidGate")) {
      gate <- list(gate)
    }
    
    # Plot gate
    gate <- cyto_plot_gate(gate,
                           channels = channels,
                           gate_line_col = gate_line_col,
                           gate_line_width = gate_line_width,
                           gate_line_type = gate_line_type
    )
    
    # Add labels to gated plot
    if(label == TRUE){
      
      # Repeat gate number of layers -
      gate <- do.call("rep", list(gate, smp))
    
      # Find center x co-ord for label position in each gate
      if (all(is.na(label_box_x))) {
        label_box_x <- unlist(lapply(unique(gate), function(x) {
          (unname(x@min) + unname(x@max)) / 2
        }))
      }
      label_box_x <- do.call("rep", list(label_box_x, smp))
    
      # Find y co-ord for each sample
      if (.all_na(label_box_y) & density_modal) {
        label_box_y <- unlist(lapply(rep(seq(1, smp), 
                                         length.out = length(gate), 
                                         each = length(unique(gate))), 
                                     function(x) {
                                     (0.5 * density_stack * 100) + ((x - 1) * density_stack * 100)
                                     }))
      
      # Too much computation required here - do this in cyto_plot_1d flowFrame
      }else if(.all_na(label_box_y) & !density_modal){
        stop("Need to supply y positions for labels")
      } 
    
      # Plot labels
      ind <- rep(seq_len(smp), 
                 each = length(unique(gate)), 
                 length.out = smp*length(unique(gate)))
      ind <- split(seq_len(smp*length(unique(gate))), ind)
    
      # Repeat arguments
      label_text <- rep(label_text, length.out = length(gate))
      label_stat <- rep(label_stat, length.out = length(gate))
      label_text_font <- rep(label_text_font, length.out = length(gate))
      label_text_size <- rep(label_text_size, length.out = length(gate))
      label_text_col <- rep(label_text_col, length.out = length(gate))
      label_box_alpha <- rep(label_box_alpha, length.out = length(gate))
    
      mapply(function(fr,x) {
        suppressMessages(cyto_plot_label(
          x = fr,
          channels = channels,
          gate = gate[x],
          trans = trans,
          text_x = label_box_x[x],
          text_y = label_box_y[x],
          text = label_text[x],
          stat = label_stat[x],
          text_font = label_text_font[x],
          text_col = label_text_col[x],
          text_size = label_text_size[x],
          box_alpha = label_box_alpha[x]
        ))
      
      }, fr.lst,
      ind)
    }
    
  # No gates - allow labels if label = TRUE
  }else{
    
    # Allow a single label per layer if label = TRUE
    if(label == TRUE){
      
      # Default x co-ord is 0.75 * x axis range
      
      # Find y co-ord for each sample
      if (.all_na(label_box_y) & density_modal) {
        label_box_y <- unlist(lapply(rep(seq(1, smp), 
                                         length.out = length(gate), 
                                         each = length(unique(gate))), 
                                     function(x) {
                                       (0.5 * density_stack * 100) + ((x - 1) * density_stack * 100)
                                     }))
        
        # Too much computation required here - do this in cyto_plot_1d flowFrame
      }else if(.all_na(label_box_y) & !density_modal){
        stop("Need to supply y positions for labels")
      } 
      
      # Add labels to plot
      mapply(function(fr,
                      label_box_x,
                      label_box_y,
                      label_text,
                      label_stat,
                      label_text_font,
                      label_text_size,
                      label_text_col,
                      label_box_alpha) {
        
        suppressMessages(cyto_plot_label(
          x = fr,
          channels = channels,
          gate = NULL,
          trans = trans,
          text_x = label_box_x,
          text_y = label_box_y,
          text = label_text,
          stat = label_stat,
          text_font = label_text_font,
          text_col = label_text_col,
          text_size = label_text_size,
          box_alpha = label_box_alpha
        ))
        
      }, fr.lst,
      label_box_x,
      label_box_y,
      label_text,
      label_stat,
      label_text_font,
      label_text_size,
      label_text_col,
      label_box_alpha)
      
    }
    
  }
  
}
