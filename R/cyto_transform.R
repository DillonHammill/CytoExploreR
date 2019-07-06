# CytoRSuite Implementation of Transformations

# Wrappers for flowCore transformation function which return transformerList
# objects and use cyto_plot to visualise the transformations.

# CYTO_TRANSFORM_ARCSINH ---------------------------------------------------------

#' Definition(s) of ArcSinh Transformation(s)
#'
#' CytoRSuite implementation of \code{flowWorkspace}
#' \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}} transformation which
#' always returns a \code{transformerList} object and displays the result of the
#' transformation(s) using \code{cyto_plot}. To combine different types of
#' transformations have a look at \code{cyto_transform_combine}.
#' \code{cyto_transform} should be used to apply the transformations to the
#' data.
#'
#' @param x an object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param channels name(s) of channel(s)/marker(s) for which transformation
#'   functions must be generated. Set to all the fluorescent channels by default
#'   if no channels/markers are supplied.
#' @param parent name of the parent population to use for generating the
#'   transformation(s).
#' @param select list of selection criteria passed to \code{cyto_select} to
#'   select a subset of samples for \code{cyto_plot}.
#' @param raw_max maximum value of the raw data.
#' @param width width of the transformed display in asymptotic decades.
#' @param trans_min minimum of transformed data in the display in asymptotic
#'   decades, set to 0 by default.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param equal_space logical indicating whether breaks should be equally
#'   spaced, set to FALSE by default.
#' @param breaks number of required breaks.
#' @param ... additional arguments passed to \code{cyto_plot}.
#'
#' @importFrom flowWorkspace asinh_Gml2 flow_trans transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh <- function(x, ...) {
  UseMethod("cyto_transform_arcsinh")
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.flowFrame <- function(x,
                                             channels,
                                             raw_max = 262144,
                                             width = 4.5,
                                             trans_min = 0,
                                             popup = FALSE,
                                             equal_space = FALSE,
                                             breaks = 6, ...) {

  # Prepare Channels
  if (missing(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }

  # Sort out transformations
  transform_list <- lapply(channels, function(z) {
    asinh_Gml2(
      T = raw_max,
      M = width,
      A = trans_min,
      inverse = FALSE
    )
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <- asinh_Gml2(
      T = raw_max,
      M = width,
      A = trans_min,
      inverse = TRUE
    )
    flow_trans(
      "arcsinh",
      z@.Data,
      inv@.Data,
      equal_space,
      breaks
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )

  print(transformer_list)
  
  # Apply transformations to data for visualisation
  transform_list <- cyto_transform_convert(transformer_list, inverse = FALSE)
  x <- transform(x, transform_list)

  # Set up plotting area
  cyto_plot_new(popup = popup)
  n <- length(channels)
  cyto_plot_layout(
    n2mfrow(n)[1],
    n2mfrow(n)[2]
  )

  # Generate plot for each channel
  lapply(channels, function(chan) {
    cyto_plot(x,
      channels = chan,
      axes_trans = transformer_list,
      ...
    )
  })

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.flowSet <- function(x,
                                           channels,
                                           select = NULL,
                                           raw_max = 262144,
                                           width = 4.5,
                                           trans_min = 0,
                                           popup = FALSE,
                                           equal_space = FALSE,
                                           breaks = 6, ...) {

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transform_arcsinh(x,
    channels = channels,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.GatingHierarchy <- function(x,
                                                   channels,
                                                   parent = "root",
                                                   raw_max = 262144,
                                                   width = 4.5,
                                                   trans_min = 0,
                                                   popup = FALSE,
                                                   equal_space = FALSE,
                                                   breaks = 6, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transform_arcsinh(x,
    channels = channels,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    quantile = quantile,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.GatingSet <- function(x,
                                             channels,
                                             parent = "root",
                                             select = NULL,
                                             raw_max = 262144,
                                             width = 4.5,
                                             trans_min = 0,
                                             popup = FALSE,
                                             equal_space = FALSE,
                                             breaks = 6, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transform_larcsinh(x,
    channels = channels,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

# CYTO_TRANSFORM_BIEX ----------------------------------------------------------

#' Definition(s) of Biexponential Transformation(s)
#'
#' CytoRSuite implementation of \code{flowWorkspace} flowJo biexponential
#' transformation which always returns a \code{transformerList}
#' object and displays in the result of the transformation(s) using
#' \code{cyto_plot}. To combine different types of transformations have a look
#' at \code{cyto_transform_combine}. \code{cyto_transform} should be used to
#' apply the transformations to the data.
#'
#' @param x an object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param channels name(s) of channel(s)/marker(s) for which transformation
#'   functions must be generated. Set to all the fluorescent channels by default
#'   if no channels/markers are supplied.
#' @param parent name of the parent population to use for generating the
#'   transformation(s).
#' @param select list of selection criteria passed to \code{cyto_select} to
#'   select a subset of samples for \code{cyto_plot}.
#' @param trans_max maximum value of the transformed data.
#' @param raw_max maximum value of the raw data.
#' @param width width of the transformed display in asymptotic decades.
#' @param trans_min minimum of transformed data in the display in asymptotic
#'   decades, set to 0 by default.
#' @param width_basis numeric to optimally position transformed data within
#'   display, set to -10 by default.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param equal_space logical indicating whether breaks should be equally
#'   spaced, set to FALSE by default.
#' @param breaks number of required breaks.
#' @param ... additional arguments passed to \code{cyto_plot}.
#'
#' @seealso \code{\link[flowWorkspace:flowJoTrans]{flowJoTrans}}
#' @seealso \code{\link[flowWorkspace:flowJo_biexp_trans]{flowJo_biexp_trans}}
#'
#' @return transformerList object.
#'
#' @importFrom flowWorkspace flowJoTrans flow_trans transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex <- function(x, ...) {
  UseMethod("cyto_transform_biex")
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.flowFrame <- function(x,
                                          channels,
                                          trans_max = 4096,
                                          raw_max = 262144,
                                          width = 4.5,
                                          trans_min = 0,
                                          width_basis = -10,
                                          popup = FALSE,
                                          equal_space = FALSE,
                                          breaks = 6, ...) {

  # Prepare Channels
  if (missing(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }

  # Sort out transformations
  transform_list <- lapply(channels, function(z) {
    flowJoTrans(
      channelRange = trans_max,
      maxValue = raw_max,
      pos = width,
      neg = trans_min,
      widthBasis = width_basis,
      inverse = FALSE
    )
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <- flowJoTrans(
      channelRange = trans_max,
      maxValue = raw_max,
      pos = width,
      neg = trans_min,
      widthBasis = width_basis,
      inverse = TRUE
    )
    flow_trans(
      "biexponential",
      z@.Data,
      inv@.Data,
      equal_space,
      breaks
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )

  # Apply transformations to data for visualisation
  transform_list <- cyto_transform_convert(transformer_list, inverse = FALSE)
  x <- transform(x, transform_list)

  # Set up plotting area
  cyto_plot_new(popup = popup)
  n <- length(channels)
  cyto_plot_layout(
    n2mfrow(n)[1],
    n2mfrow(n)[2]
  )

  # Generate plot for each channel
  lapply(channels, function(chan) {
    cyto_plot(x,
      channels = chan,
      axes_trans = transformer_list,
      ...
    )
  })

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.flowSet <- function(x,
                                        channels,
                                        select = NULL,
                                        trans_max = 4096,
                                        raw_max = 262144,
                                        width = 4.5,
                                        trans_min = 0,
                                        width_basis = -10,
                                        popup = FALSE,
                                        equal_space = FALSE,
                                        breaks = 6, ...) {

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transform_biex(x,
    channels = channels,
    trans_max = trans_max,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    width_basis = width_basis,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.GatingHierarchy <- function(x,
                                                channels,
                                                parent = "root",
                                                select = NULL,
                                                trans_max = 4096,
                                                raw_max = 262144,
                                                width = 4.5,
                                                trans_min = 0,
                                                width_basis = -10,
                                                popup = FALSE,
                                                equal_space = FALSE,
                                                breaks = 6, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transform_biex(x,
    channels = channels,
    trans_max = trans_max,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    width_basis = width_basis,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.GatingSet <- function(x,
                                          channels,
                                          parent = "root",
                                          select = NULL,
                                          trans_max = 4096,
                                          raw_max = 262144,
                                          width = 4.5,
                                          trans_min = 0,
                                          width_basis = -10,
                                          popup = FALSE,
                                          equal_space = FALSE,
                                          breaks = 6, ...) {
  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transform_biex(x,
    channels = channels,
    trans_max = trans_max,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    width_basis = width_basis,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

# CYTO_TRANSFORM_LOGICLE -------------------------------------------------------

#' Definition(s) of Logicle Transformation(s)
#'
#' \code{CytoRSuite} implementation of \code{flowCore} logicle transformation
#' which allows the use of multiple samples to estimate transformation
#' parameters, always returns a \code{transformerList} object and displays the
#' result of the transformation(s) using \code{cyto_plot}. To combine different
#' types of transformations have a look at \code{cyto_transform_combine}.
#' \code{cyto_transform} should be used to apply the transformations to the
#' data.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param channels name(s) of channel(s)/marker(s) for which transformation
#'   functions must be generated. Set to all the fluorescent channels by default
#'   if no channels/markers are supplied.
#' @param parent name of the parent population to use for generating the
#'   transformation(s).
#' @param select list of selection criteria passed to \code{cyto_select} to
#'   select a subset of samples for \code{cyto_plot}.
#' @param width width of transformed display in asymptotic decades.
#' @param raw_max maximum value of the raw data.
#' @param trans_min minimum of the transformed data in symptotic decades, set to
#'   0 by default.
#' @param quantile quantile of the negative data value (used to calculate
#'   linearisation width (w)).
#' @param type indicates whether the "instrument" or "data" limits should be
#'   used as the range for the data.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param equal_space logical indicating whether breaks should be equally
#'   spaced, set to FALSE by default.
#' @param breaks number of required breaks.
#' @param ... additional arguments passed to \code{cyto_plot}.
#'
#' @return transformerList object.
#'
#' @importFrom flowCore inverseLogicleTransform
#' @importFrom flowWorkspace flow_trans transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowCore:logicleTransform]{estimateLogicle}}
#' @seealso \code{\link[flowWorkspace:estimateLogicle]{estimateLogicle}}
#'
#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle <- function(x, ...) {
  UseMethod("cyto_transform_logicle")
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.flowFrame <- function(x,
                                             channels,
                                             width = 4.5,
                                             raw_max = 262144,
                                             trans_min = 0,
                                             quantile = 0.05,
                                             type = "instrument",
                                             popup = FALSE,
                                             equal_space = FALSE,
                                             breaks = 6, ...) {

  # Prepare Channels
  if (missing(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }

  # Sort out transformations
  transform_list <- flowCore:::.estimateLogicle(x,
    channels = channels,
    m = width,
    t = raw_max,
    a = trans_min,
    q = quantile,
    type = type
  )
  transformer_list <- lapply(transform_list, function(z) {
    inv <- inverseLogicleTransform(trans = z)
    flow_trans(
      "logicle",
      z@.Data,
      inv@.Data,
      equal_space,
      breaks
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )

  # Apply transformations to data for visualisation
  transform_list <- cyto_transform_convert(transformer_list, inverse = FALSE)
  x <- transform(x, transform_list)

  # Set up plotting area
  cyto_plot_new(popup = popup)
  n <- length(channels)
  cyto_plot_layout(
    n2mfrow(n)[1],
    n2mfrow(n)[2]
  )

  # Generate plot for each channel
  lapply(channels, function(chan) {
    cyto_plot(x,
      channels = chan,
      axes_trans = transformer_list,
      ...
    )
  })

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.flowSet <- function(x,
                                           channels,
                                           select = NULL,
                                           trans_max = 4.5,
                                           raw_max = 262144,
                                           trans_min = 0,
                                           quantile = 0.05,
                                           type = "instrument",
                                           popup = FALSE,
                                           equal_space = FALSE,
                                           breaks = 6, ...) {

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transform_logicle(x,
    channels = channels,
    trans_max = trans_max,
    raw_max = raw_max,
    trans_min = trans_min,
    quantile = quantile,
    type = type,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.GatingHierarchy <- function(x,
                                                   channels,
                                                   parent = "root",
                                                   trans_max = 4.5,
                                                   raw_max = 262144,
                                                   trans_min = 0,
                                                   quantile = 0.05,
                                                   type = "instrument",
                                                   popup = FALSE,
                                                   equal_space = FALSE,
                                                   breaks = 6, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transform_logicle(x,
    channels = channels,
    trans_max = trans_max,
    raw_max = raw_max,
    trans_min = trans_min,
    quantile = quantile,
    type = type,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.GatingSet <- function(x,
                                             channels,
                                             parent = "root",
                                             select = NULL,
                                             trans_max = 4.5,
                                             raw_max = 262144,
                                             trans_min = 0,
                                             quantile = 0.05,
                                             type = "instrument",
                                             popup = FALSE,
                                             equal_space = FALSE,
                                             breaks = 6, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transform_logicle(x,
    channels = channels,
    trans_max = trans_max,
    raw_max = raw_max,
    trans_min = trans_min,
    quantile = quantile,
    type = type,
    popup = popup,
    equal_space = equal_space,
    breaks = breaks, ...
  )

  # Return transformerList
  return(transformer_list)
}
