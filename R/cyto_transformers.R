## CytoExploreR Implementation of Transformations ------------------------------

# Wrappers for flowCore/flowWorkspace transformation functions which return
# transformerList objects and use cyto_plot to visualise the transformations.

# All methods must have select as an argument to allowing passing of arguments
# from cyto_transform.

## CYTO_TRANSFORMER_LOG --------------------------------------------------------

#' Definition(s) of Log Transformation(s)
#'
#' CytoExploreR implementation of \code{flowWorkspace}
#' \code{\link[flowWorkspace:flowjo_flog]{flowjo_flog}} transformation which
#' always returns a \code{\link[flowWorkspace:transformerList]{transformerList}}
#' object and displays the result of the transformation(s) using
#' \code{\link{cyto_plot}}. To combine different types of transformations have a
#' look at \code{\link{cyto_transformer_combine}}. \code{\link{cyto_transform}}
#' should be used to apply the transformations to the data.
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
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:flowjo_flog]{flowjo_flog}}.
#'
#' @return a \code{transformerList} object.
#'
#' @importFrom flowWorkspace flowjo_flog flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowWorkspace:flowjo_flog]{flowjo_flog}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_logicle}}
#' @seealso \code{\link{cyto_transformer_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log <- function(x, ...){
  UseMethod("cyto_transformer_log")
}

#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log.GatingSet <- function(x,
                                           channels = NULL,
                                           parent = "root",
                                           select = NULL,
                                           plot = TRUE,
                                           popup = FALSE, ...){
  
  # Extract data
  x <- cyto_extract(x, parent = parent)
  
  # Select data
  x <- cyto_select(x, select)
  
  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")
  
  # Call to flowFrame method
  transformer_list <- cyto_transformer_log(x,
                                           channels = channels,
                                           plot = plot,
                                           popup = popup,
                                           ...
  )
  
  # Return transformerList
  return(transformer_list)
  
}

#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log.GatingHierarchy <- function(x,
                                                 channels = NULL,
                                                 parent = "root",
                                                 select = NULL,
                                                 plot = TRUE,
                                                 popup = FALSE, ...){
    
  # Extract data
  x <- cyto_extract(x, parent = parent)
  
  # Call to flowFrame method
  transformer_list <- cyto_transformer_log(x,
                                           channels = channels,
                                           plot = plot,
                                           popup = popup,
                                           ...
  )
  
  # Return transformerList
  return(transformer_list)
  
}

#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log.flowSet <- function(x,
                                         channels = NULL,
                                         select = NULL,
                                         plot = TRUE,
                                         popup = FALSE, ...){
  
  # Select data
  x <- cyto_select(x, select)
  
  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")
  
  # Call to flowFrame method
  transformer_list <- cyto_transformer_log(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )
  
  # Return transformerList
  return(transformer_list)
  
}

#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log.flowFrame <- function(x,
                                           channels = NULL,
                                           select = NULL,
                                           plot = TRUE,
                                           popup = FALSE, ...){
  
  # Prepare Channels
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }
  
  # Sort out transformations
  transform_list <- lapply(channels, function(z) {
    flowjo_flog(
      inverse = FALSE,
      ...
    )
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <- flowjo_flog(
      inverse = TRUE,
      ...
    )
    flow_trans(
      "log",
      z@.Data,
      inv@.Data
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )
  
  # Construct plots
  if(plot == TRUE){
    
    # Sample for speed
    x <- cyto_sample(x, display = 25000, seed = 56)
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)
    
    # Old graphics parameters
    old_pars <- .par("mfrow")
    
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
                title = NA
      )
    })
    
    # Reset graphics parameters
    par(old_pars)
    
  }
  
  # Return transformerList
  return(transformer_list)
  
}

## CYTO_TRANSFORMER_ARCSINH ------------------------------------------------------

#' Definition(s) of arcsinh transformation(s)
#'
#' CytoExploreR implementation of \code{flowWorkspace}
#' \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}} transformation which
#' always returns a \code{transformerList} object and displays the result of the
#' transformation(s) using \code{cyto_plot}. To combine different types of
#' transformations have a look at \code{cyto_transformer_combine}.
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
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}}.
#'
#' @return a \code{transformerList} object.
#'
#' @importFrom flowWorkspace asinh_Gml2 flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}}
#' @seealso \code{\link{cyto_transformer_log}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_logicle}}
#' @seealso \code{\link{cyto_transformer_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh <- function(x, ...) {
  UseMethod("cyto_transformer_arcsinh")
}

#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh.GatingSet <- function(x,
                                               channels = NULL,
                                               parent = "root",
                                               select = NULL,
                                               plot = TRUE,
                                               popup = FALSE, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_arcsinh(x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh.GatingHierarchy <- function(x,
                                                   channels = NULL,
                                                   parent = "root",
                                                   select = NULL,
                                                   plot = TRUE,
                                                   popup = FALSE, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transformer_arcsinh(x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh.flowSet <- function(x,
                                           channels = NULL,
                                           select = NULL,
                                           plot = TRUE,
                                           popup = FALSE,...) {

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_arcsinh(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh.flowFrame <- function(x,
                                               channels = NULL,
                                               select = NULL,
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
      inverse = FALSE,
      ...
    )
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <- asinh_Gml2(
      inverse = TRUE,
      ...
    )
    flow_trans(
      "arcsinh",
      z@.Data,
      inv@.Data
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )
  
  # Construct plots
  if(plot == TRUE){
      
    # Sample for speed
    x <- cyto_sample(x, display = 25000, seed = 56)
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)

    # Old graphics parameters
    old_pars <- .par("mfrow")
    
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
        title = NA
      )
    })
    
    # Reset graphics parameters
    par(old_pars)
    
  }

  # Return transformerList
  return(transformer_list)
}

## CYTO_TRANSFORMER_BIEX ---------------------------------------------------------

#' Definition(s) of biexponential transformation(s)
#'
#' CytoExploreR implementation of \code{flowWorkspace} flowJo biexponential
#' transformation which always returns a
#' \code{\link[flowWorkspace:transformerList]{transformerList}} object and
#' displays in the result of the transformation(s) using
#' \code{\link{cyto_plot}}. To combine different types of transformations have a
#' look at \code{\link{cyto_transformer_combine}}. \code{\link{cyto_transform}}
#' should be used to apply the transformations to the data.
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
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:flowjo_biexp]{flowjo_biexp}}.
#'
#' @return a \code{transformerList} object.
#'
#' @importFrom flowWorkspace flowjo_biexp flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link[flowWorkspace:flowjo_biexp]{flowjo_biexp}}
#' @seealso \code{\link[flowWorkspace:flowJo_biexp_trans]{flowJo_biexp_trans}}
#' @seealso \code{\link{cyto_transformer_log}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_logicle}}
#' @seealso \code{\link{cyto_transformer_combine}}
#' @seealso \code{\link{cyto_transform}}
#' 
#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex <- function(x, ...) {
  UseMethod("cyto_transformer_biex")
}

#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex.GatingSet <- function(x,
                                          channels = NULL,
                                          parent = "root",
                                          select = NULL,
                                          plot = TRUE,
                                          popup = FALSE, ...) {
  
  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_biex(x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex.GatingHierarchy <- function(x,
                                                channels = NULL,
                                                parent = "root",
                                                select = NULL,
                                                plot = TRUE,
                                                popup = FALSE) {
  
  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transformer_biex(x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex.flowSet <- function(x,
                                        channels = NULL,
                                        select = NULL,
                                        plot = TRUE,
                                        popup = FALSE, ...) {
  
  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_biex(x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex.flowFrame <- function(x,
                                          channels = NULL,
                                          select = NULL,
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
    flowjo_biexp(
      inverse = FALSE,
      ...
    )
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <- flowjo_biexp(
      inverse = TRUE,
      ...
    )
    flow_trans(
      "biexponential",
      z@.Data,
      inv@.Data
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )

  # Construct plots
  if(plot == TRUE){
    
    # Sample for speed
    x <- cyto_sample(x, display = 25000, seed = 56)
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)

    # Old graphics parameters
    old_pars <- .par("mfrow")
     
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
        title = NA
      )
    })
    
    # Reset graphics parameters
    par(old_pars)
  
  }

  # Return transformerList
  return(transformer_list)
}

## CYTO_TRANSFORMER_LOGICLE ----------------------------------------------------

#' Definition(s) of logicle transformation(s)
#'
#' \code{CytoExploreR} implementation of \code{flowCore} logicle transformation
#' which allows the use of multiple samples to estimate transformation
#' parameters, always returns a \code{transformerList} object and displays the
#' result of the transformation(s) using \code{cyto_plot}. To combine different
#' types of transformations have a look at \code{cyto_transformer_combine}.
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
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param popup logical indicating whether the plots should be constructed in a
#'   pop-up window.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:estimateLogicle.GatingHierarchy]{estimateLogicle}}.
#'
#' @return a \code{transformerList} object.
#'
#' @importFrom flowCore inverseLogicleTransform
#' @importFrom flowWorkspace flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowCore:logicleTransform]{estimateLogicle}}
#' @seealso \code{\link[flowWorkspace:estimateLogicle.GatingHierarchy]{estimateLogicle}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle <- function(x, ...) {
  UseMethod("cyto_transformer_logicle")
}

#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle.GatingSet <- function(x,
                                             channels = NULL,
                                             parent = "root",
                                             select = NULL,
                                             plot = TRUE,
                                             popup = FALSE, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_logicle(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle.GatingHierarchy <- function(x,
                                                   channels = NULL,
                                                   parent = "root",
                                                   select = NULL,
                                                   plot = TRUE,
                                                   popup = FALSE, ...) {

  # Extract data
  x <- cyto_extract(x, parent = parent)

  # Call to flowFrame method
  transformer_list <- cyto_transformer_logicle(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle.flowSet <- function(x,
                                           channels = NULL,
                                           select = NULL,
                                           plot = FALSE,
                                           popup = FALSE, ...) {

  # Select data
  x <- cyto_select(x, select)

  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_logicle(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle.flowFrame <- function(x,
                                             channels = NULL,
                                             select = NULL,
                                             plot = TRUE,
                                             popup = FALSE, ...) {

  # Prepare Channels
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }

  # Sort out transformations
  transform_list <- flowCore:::.estimateLogicle(
    x,
    channels = channels,
    ...
  )
  transformer_list <- lapply(transform_list, function(z) {
    inv <- inverseLogicleTransform(trans = z)
    flow_trans(
      "logicle",
      z@.Data,
      inv@.Data
    )
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )

  # Construct plots
  if(plot == TRUE){
    
    # Sample for speed
    x <- cyto_sample(x, display = 25000, seed = 56)
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)

    # Old graphics parameters
    old_pars <- .par("mfrow")
    
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
        title = NA
      )
    })
    
    # Reset graphic parameters
    par(old_pars)
    
  }

  # Return transformerList
  return(transformer_list)
}

## CYTO_TRANSFORMER_COMBINE ------------------------------------------------------

#' Combine cyto_transformer definitions
#'
#' \code{cyto_transformer_combine} makes it easy to combine transformation
#' definitions obtained from \code{\link{cyto_transformer_log}},
#' \code{\link{cyto_transformer_arcsinh}}, \code{\link{cyto_transformer_biex}}
#' and/or \code{\link{cyto_transformer_logicle}} prior to applying these
#' transformations to the data using \code{\link{cyto_transform}}.
#'
#' @param ... objects of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} to be
#'   combined into a single \code{transformerList} object.
#'
#' @importFrom flowWorkspace transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_transformer_log}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_logicle}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @export
cyto_transformer_combine <- function(...) {
  
  # Combine transformerList objects
  transformer_list <- c(...)

  # Convert to transformerList
  transformer_list <- transformerList(names(transformer_list),
                                      transformer_list)
    
  # Return transformations in a single transformerList for cyto_transform
  return(transformer_list)
  
}

## .CYTO_TRANSFORMER_COMPLETE --------------------------------------------------

#' Get complete transformerList for compensation functions
#' @noRd
.cyto_transformer_complete <- function(x, ...){
  UseMethod(".cyto_transformer_complete")
}

#' @noRd
.cyto_transformer_complete.flowFrame <- function(x,
                                                 axes_trans = NA) {
  
  # Fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # In case NULL axes_trans
  if(is.null(axes_trans)){
    axes_trans <- NA
  }
  
  # Transformations supplied
  if(!.all_na(axes_trans)){
    # Must be transformerList
    if(!inherits(axes_trans, "transformerList")){
      stop("'axes_trans' must be an object of class transformerList.")
    }
    # Some transformations are missing
    if(!all(channels %in% names(axes_trans))){
      # Missing channels
      trans_chans <- channels[channels %in% names(axes_trans)]
      linear_chans <- channels[!channels %in% names(axes_trans)]
      # Get remaining transformers
      trans <- cyto_transformer_biex(x,
                                     channels = linear_chans,
                                     plot = FALSE)
      # Combine transformers
      axes_trans <- cyto_transformer_combine(axes_trans[trans_chans], 
                                             trans)
    }
  # No transformations supplied
  }else{
    # Get biex transformers
    axes_trans <- cyto_transformer_biex(x, 
                                        channels = channels, 
                                        plot = FALSE)
    
  }
  
  # Return complete transformerList
  return(axes_trans)
  
}

#' @noRd
.cyto_transformer_complete.flowSet <- function(x, 
                                               axes_trans = NA){
  
  # Fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # In case NULL axes_trans
  if(is.null(axes_trans)){
    axes_trans <- NA
  }
  
  # Transformations supplied
  if(!.all_na(axes_trans)){
    # Must be transformerList
    if(!inherits(axes_trans, "transformerList")){
      stop("'axes_trans' must be an object of class transformerList.")
    }
    # Some stranformations are missing
    if(!all(channels %in% names(axes_trans))){
      # Missing channels
      trans_chans <- channels[channels %in% names(axes_trans)]
      linear_chans <- channels[!channels %in% names(axes_trans)]
      # Get remaining transformers
      trans <- cyto_transformer_biex(x,
                                     channels = linear_chans,
                                     plot = FALSE)
      # Combine transformers
      axes_trans <- cyto_transformer_combine(axes_trans[trans_chans], 
                                             trans)
    }
    # No transformations supplied
  }else{
    # Get biex transformers
    axes_trans <- cyto_transformer_biex(x, 
                                        channels = channels, 
                                        plot = FALSE)
  }
  
  # Return complete transformerList
  return(axes_trans)
  
}

#' @noRd
.cyto_transformer_complete.GatingHierarchy <- function(x,
                                                       axes_trans = NA){
  
  # Fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # Extract transformations from GatingHierarchy
  gh_trans <- x@transformation
  
  # In case NULL axes_trans
  if(is.null(axes_trans)){
    axes_trans <- NA
  }
  
  # No transformations found in GatingHierarchy
  if(length(gh_trans) == 0){
    # No transformations supplied
    if(.all_na(axes_trans)){
      # Get biex transformers for all channels
      axes_trans <- cyto_transformer_biex(x, 
                                          channels = channels, 
                                          plot = FALSE)
    # Transformations supplied
    }else if(!.all_na(axes_trans)){
      # All channels covered
      if(all(channels %in% names(axes_trans))){
        axes_trans <- axes_trans
      # Some channels covered
      }else if(any(channels %in% names(axes_trans))){
        # Channels not covered
        trans_chans <- channels[channels %in% names(axes_trans)]
        linear_chans <- channels[!channels %in% names(axes_trans)]
        # Get biex transformers for missing channels
        trans <- cyto_transformer_biex(x,
                                       channels = linear_chans,
                                       plot = FALSE)
        # Combine transformers into transformerList
        axes_trans <- cyto_transformer_combine(axes_trans[trans_chans],
                                               trans)
      # No channels covered
      }else if(!any(channels %in% names(axes_trans))){
        # Get biex transformers for all channels
        axes_trans <- cyto_transformer_biex(x,
                                            channels = channels,
                                            plot = FALSE)
      }
    }

  # Transformations found in GatingHierarchy
  }else{
    # All the transformations have been applied to GatingHierarchy
    if(all(channels %in% names(gh_trans))){
      # Regardless of axes_trans use complete gh_trans
      axes_trans <- cyto_transformer_combine(gh_trans[channels])
    # Some transformations have been applied to GatingHierachy
    }else if(any(channels %in% names(gh_trans))){
      # Channels not covered
      gh_trans_chans <- channels[channels %in% names(gh_trans)]
      linear_chans <- channels[!channels %in% names(gh_trans)]
      # No transformations supplied
      if(.all_na(axes_trans)){
        # Get biex transformers for missing channels
        trans <- cyto_transformer_biex(x,
                                       channels = linear_chans,
                                       plot = FALSE)
        # Combine transformers into transformerList
        axes_trans <- cyto_transformer_combine(gh_trans[gh_trans_chans], trans)
      # Some transformations supplied
      }else if(!.all_na(axes_trans)){
        # All missing transformations covered
        if(all(linear_chans %in% names(axes_trans))){
          axes_trans <- cyto_transformer_combine(gh_trans[gh_trans_chans],
                                                 axes_trans[linear_chans])
        # Some missing transformations covered
        }else if(any(linear_chans %in% names(axes_trans))){
          # Channels not covered
          trans_chans <- linear_chans[linear_chans%in% names(axes_trans)]
          extra_chans <- linear_chans[!linear_chans%in% names(axes_trans)]
          # Get biex transformers for missing channels
          trans <- cyto_transformer_biex(x,
                                         channels = extra_chans,
                                         plot = FALSE)
          # Combine transformers into transformerList
          axes_trans <- cyto_transformer_combine(gh_trans[gh_trans_chans],
                                                 axes_trans[trans_chans],
                                                 trans)
        # No missing transformations covered
        }else if(!any(linear_chans %in% names(axes_trans))){
          # Get biex transformers for all channels
          axes_trans <- cyto_transformer_biex(x,
                                              channels = linear_chans,
                                              plot = FALSE)
        }
      }
    # Applied transformations are not in fluorescent channels
    }else if(!any(channels %in% names(gh_trans))){
      # No transformations supplied
      if(.all_na(axes_trans)){
        # Get biex transformers for all channels
        axes_trans <- cyto_transformer_biex(x,
                                            channels = channels,
                                            plot = FALSE)
      # Some transformations supplied
      }else if(!.all_na(axes_trans)){
        # All channels covered by axes_trans
        if(all(channels %in% names(axes_trans))){
          axes_trans <- cyto_transformer_combine(axes_trans[channels])
        # Some channels are covered by axes_trans
        }else if(any(channels %in% names(axes_trans))){
          # Channels not covered
          trans_chans <- channels[channels %in% names(axes_trans)]
          linear_chans <- channels[!channels %in% names(axes_trans)]
          # Get biex transformers for missing channels
          trans <- cyto_transformer_biex(x,
                                         channels = linear_chans,
                                         plot = FALSE)
          # Combine transformers into transformerList
          axes_trans <- cyto_transformer_combine(axes_trans[trans_chans],
                                                 trans)
        # No channels covered by axes_trans
        } else if(!any(channels %in% names(axes_trans))){
          # Get biex transformers for all channels
          axes_trans <- cyto_transformer_biex(x,
                                              channels = channels,
                                              plot = FALSE)
        }
      }
    }
  }
  
  # Return complete transformerList
  return(axes_trans)
  
}

#' @noRd
.cyto_transformer_complete.GatingSet <- function(x,
                                                 axes_trans = NA){
  
  # Extract GatingHierachy
  x <- x[[1]]
  
  # Make call to GatingHierchy method
  axes_trans <- .cyto_transformer_complete(x,
                                           axes_trans = axes_trans)
  
  # Return complete transformerList
  return(axes_trans)
  
}
