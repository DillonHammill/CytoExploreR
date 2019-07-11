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
#' @param equal_space logical indicating whether breaks should be equally
#'   spaced, set to FALSE by default.
#' @param breaks number of required breaks.
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... not in use.
#'
#' @importFrom flowWorkspace asinh_Gml2 flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}}
#' @seealso \code{\link{cyto_transform_biex}}
#' @seealso \code{\link{cyto_transform_logicle}}
#' @seealso \code{\link{cyto_transform_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh <- function(x, ...) {
  UseMethod("cyto_transform_arcsinh")
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.flowFrame <- function(x,
                                             channels = NULL,
                                             raw_max = 262144,
                                             width = 4.5,
                                             trans_min = 0,
                                             equal_space = FALSE,
                                             breaks = 5,
                                             plot = TRUE,
                                             popup = FALSE, ...) {

  # Prepare Channels
  if (is.null(channels)) {
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
  
  # Construct plots
  if(plot == TRUE){
      
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)

    # Old graphics parameters
    old_pars <- par("mfrow")
    
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
        axes_trans = transformer_list
      )
    })
    
    # Reset graphics parameters
    par(old_pars)
    
  }

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.flowSet <- function(x,
                                           channels = NULL,
                                           select = NULL,
                                           raw_max = 262144,
                                           width = 4.5,
                                           trans_min = 0,
                                           equal_space = FALSE,
                                           breaks = 5, 
                                           plot = TRUE,
                                           popup = FALSE,...) {

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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.GatingHierarchy <- function(x,
                                                   channels = NULL,
                                                   parent = "root",
                                                   raw_max = 262144,
                                                   width = 4.5,
                                                   trans_min = 0,
                                                   equal_space = FALSE,
                                                   breaks = 5,
                                                   plot = TRUE,
                                                   popup = FALSE) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transform_arcsinh(x,
    channels = channels,
    raw_max = raw_max,
    width = width,
    trans_min = trans_min,
    quantile = quantile,
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_arcsinh
#' @export
cyto_transform_arcsinh.GatingSet <- function(x,
                                             channels = NULL,
                                             parent = "root",
                                             select = NULL,
                                             raw_max = 262144,
                                             width = 4.5,
                                             trans_min = 0,
                                             equal_space = FALSE,
                                             breaks = 5,
                                             plot = TRUE,
                                             popup = FALSE) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
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
#' @param equal_space logical indicating whether breaks should be equally
#'   spaced, set to FALSE by default.
#' @param breaks number of required breaks.
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... not in use.
#'
#' @seealso \code{\link[flowWorkspace:flowJoTrans]{flowJoTrans}}
#' @seealso \code{\link[flowWorkspace:flowJo_biexp_trans]{flowJo_biexp_trans}}
#' @seealso \code{\link{cyto_transform_arcsinh}}
#' @seealso \code{\link{cyto_transform_logicle}}
#' @seealso \code{\link{cyto_transform_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @return transformerList object.
#'
#' @importFrom flowWorkspace flowJoTrans flow_trans transformerList
#' @importFrom graphics par
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
                                          channels = NULL,
                                          trans_max = 4096,
                                          raw_max = 262144,
                                          width = 4.5,
                                          trans_min = 0,
                                          width_basis = -10,
                                          equal_space = FALSE,
                                          breaks = 5, 
                                          plot = TRUE,
                                          popup = FALSE,...) {
  
  # Prepare Channels
  if (is.null(channels)) {
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

  # Construct plots
  if(plot == TRUE){
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)

    # Old graphics parameters
    old_pars <- par("mfrow")
     
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
        axes_trans = transformer_list
      )
    })
    
    # Reset graphics parameters
    par(old_pars)
  
  }

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.flowSet <- function(x,
                                        channels = NULL,
                                        select = NULL,
                                        trans_max = 4096,
                                        raw_max = 262144,
                                        width = 4.5,
                                        trans_min = 0,
                                        width_basis = -10,
                                        equal_space = FALSE,
                                        breaks = 5,
                                        plot = TRUE,
                                        popup = FALSE, ...) {
  
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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.GatingHierarchy <- function(x,
                                                channels = NULL,
                                                parent = "root",
                                                select = NULL,
                                                trans_max = 4096,
                                                raw_max = 262144,
                                                width = 4.5,
                                                trans_min = 0,
                                                width_basis = -10,
                                                equal_space = FALSE,
                                                breaks = 5,
                                                plot = TRUE,
                                                popup = FALSE) {
  
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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_biex
#' @export
cyto_transform_biex.GatingSet <- function(x,
                                          channels = NULL,
                                          parent = "root",
                                          select = NULL,
                                          trans_max = 4096,
                                          raw_max = 262144,
                                          width = 4.5,
                                          trans_min = 0,
                                          width_basis = -10,
                                          equal_space = FALSE,
                                          breaks = 5,
                                          plot = TRUE,
                                          popup = FALSE) {
  
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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
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
#' @param equal_space logical indicating whether breaks should be equally
#'   spaced, set to FALSE by default.
#' @param breaks number of required breaks.
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... not in use.
#'
#' @return transformerList object.
#'
#' @importFrom flowCore inverseLogicleTransform
#' @importFrom flowWorkspace flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowCore:logicleTransform]{estimateLogicle}}
#' @seealso \code{\link[flowWorkspace:estimateLogicle]{estimateLogicle}}
#' @seealso \code{\link{cyto_transform_arcsinh}}
#' @seealso \code{\link{cyto_transform_biex}}
#' @seealso \code{\link{cyto_transform_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle <- function(x, ...) {
  UseMethod("cyto_transform_logicle")
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.flowFrame <- function(x,
                                             channels = NULL,
                                             width = 4.5,
                                             raw_max = 262144,
                                             trans_min = 0,
                                             quantile = 0.05,
                                             type = "instrument",
                                             equal_space = FALSE,
                                             breaks = 5,
                                             plot = TRUE,
                                             popup = FALSE, ...) {

  # Prepare Channels
  if (is.null(channels)) {
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

  # Construct plots
  if(plot == TRUE){
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)

    # Old graphics parameters
    old_pars <- par("mfrow")
    
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
        axes_trans = transformer_list
      )
    })
    
    # Reset graphic parameters
    par(old_pars)
    
  }

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.flowSet <- function(x,
                                           channels = NULL,
                                           select = NULL,
                                           trans_max = 4.5,
                                           raw_max = 262144,
                                           trans_min = 0,
                                           quantile = 0.05,
                                           type = "instrument",
                                           equal_space = FALSE,
                                           breaks = 5,
                                           plot = FALSE,
                                           popup = FALSE, ...) {

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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.GatingHierarchy <- function(x,
                                                   channels = NULL,
                                                   parent = "root",
                                                   trans_max = 4.5,
                                                   raw_max = 262144,
                                                   trans_min = 0,
                                                   quantile = 0.05,
                                                   type = "instrument",
                                                   equal_space = FALSE,
                                                   breaks = 5,
                                                   plot = TRUE,
                                                   popup = FALSE) {

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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transform_logicle
#' @export
cyto_transform_logicle.GatingSet <- function(x,
                                             channels = NULL,
                                             parent = "root",
                                             select = NULL,
                                             trans_max = 4.5,
                                             raw_max = 262144,
                                             trans_min = 0,
                                             quantile = 0.05,
                                             type = "instrument",
                                             equal_space = FALSE,
                                             breaks = 5,
                                             plot = TRUE,
                                             popup = FALSE) {

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
    equal_space = equal_space,
    breaks = breaks,
    plot = plot,
    popup = popup
  )

  # Return transformerList
  return(transformer_list)
}

# CYTO_TRANSFORM_COMBINE -------------------------------------------------------

#' Combine Transformation Definitions
#'
#' \code{cyto_transform_combine} makes it easy to combine transformation
#' definitions obtained from \code{cyto_transform_arcsinh},
#' \code{cyto_transform_biex} and/or \code{cyto_transform_logicle} prior to
#' applying these transformations to the data using \code{cyto_transform}.
#'
#' @param ... objects of class \code{transformerList} to be combined into a
#'   single \code{transformerList} object.
#'
#' @importFrom flowWorkspace transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_transform_arcsinh}}
#' @seealso \code{\link{cyto_transform_biex}}
#' @seealso \code{\link{cyto_transform_logicle}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @export
cyto_transform_combine <- function(...) {
  
  # Combine transformerList objects
  transformer_list <- c(...)

  # Convert to transformerList
  transformer_list <- transformerList(names(transformer_list),
                                      transformer_list)
    
  # Return transformations in a single transformerList for cyto_transform
  return(transformer_list)
  
}
