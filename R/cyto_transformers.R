## CytoExploreR Implementation of Transformations ------------------------------

# Wrappers for flowCore/flowWorkspace transformation functions which return
# transformerList objects and use cyto_plot to visualise the transformations.

# All methods must have select as an argument to allow passing of arguments from
# cyto_transform.

## CYTO_TRANSFORMER_LOG --------------------------------------------------------

#' Definition(s) of Log Transformation(s)
#'
#' CytoExploreR implementation of \code{flowWorkspace}
#' \code{\link[flowWorkspace:flowjo_log_trans]{flowjo_log_trans}} transformation
#' which always returns a
#' \code{\link[flowWorkspace:transformerList]{transformerList}} object and
#' displays the result of the transformation(s) using \code{\link{cyto_plot}}.
#' To combine different types of transformations have a look at
#' \code{\link{cyto_transformer_combine}}. \code{\link{cyto_transform}} should
#' be used to apply the transformations to the data.
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
#' @param display number or frequency of events to display in plots, set to
#'   25000 events by default. See \code{\link{cyto_plot}} for details.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:flowjo_log_trans]{flowjo_log_trans}} and
#'   \code{\link{cyto_plot}}.
#'
#' @return a \code{transformerList} object.
#'
#' @importFrom flowWorkspace flowjo_log_trans flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowWorkspace:flowjo_log_trans]{flowjo_log_trans}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_logicle}}
#' @seealso \code{\link{cyto_transformer_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @name cyto_transformer_log
NULL

#' @noRd
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
                                           popup = FALSE,
                                           display = 25000,
                                           axes_limits = "machine", ...){
  
  # Extract data
  x <- cyto_extract(x, 
                    parent = parent, 
                    copy = TRUE)
  
  # Select data
  x <- cyto_select(x, select)
  
  # Coerce to flowFrame
  x <- cyto_convert(x, "flowFrame")
  
  # Call to flowFrame method
  transformer_list <- cyto_transformer_log(x,
                                           channels = channels,
                                           plot = plot,
                                           popup = popup,
                                           display = display,
                                           axes_limits = axes_limits,
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
                                                 popup = FALSE,
                                                 display = 25000,
                                                 axes_limits = "machine", ...){
    
  # Extract data
  x <- cyto_extract(x, 
                    parent = parent, 
                    copy = TRUE)
  
  # Call to flowFrame method
  transformer_list <- cyto_transformer_log(x,
                                           channels = channels,
                                           plot = plot,
                                           popup = popup,
                                           display = display,
                                           axes_limits = axes_limits,
                                           ...
  )
  
  # Return transformerList
  return(transformer_list)
  
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log.flowSet <- function(x,
                                         channels = NULL,
                                         parent = "root",
                                         select = NULL,
                                         plot = TRUE,
                                         popup = FALSE, 
                                         display = 25000,
                                         axes_limits = "machine", ...){
  
  # COPY
  x <- cyto_copy(x)
  
  # SELECT
  x <- cyto_select(x, select)
  
  # FLOWFRANE
  x <- cyto_convert(x, "flowFrame")
  
  # Call to flowFrame method
  transformer_list <- cyto_transformer_log(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )
  
  # Return transformerList
  return(transformer_list)
  
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_log
#' @export
cyto_transformer_log.flowFrame <- function(x,
                                           channels = NULL,
                                           parent = "root",
                                           select = NULL,
                                           plot = TRUE,
                                           popup = FALSE,
                                           display = 25000,
                                           axes_limits = "machine", ...){
  
  # PULL DOWN ALL ARGUMENTS
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  
  # COPY
  x <- cyto_copy(x)
  
  # PREPARE CHANNELS
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }
  
  # FLOWJO_LOG ARGUMENTS
  flowjo_log_args <- formalArgs("flowjo_log_trans")
  
  # TRANSFORMATIONS
  transformer_list <- lapply(channels, function(z) {
    do.call("flowjo_log_trans", 
            args[names(args) %in% flowjo_log_args])
  })
  names(transformer_list) <- channels
  transformer_list <- transformerList(
    from = channels,
    trans = transformer_list
  )
  
  # PLOTS
  if(plot == TRUE){
    
    # GRAPHICAL PARAMETERS
    old_pars <- .par("mfrow")
    on.exit(par(old_pars))
    
    # CYTO_PLOT ARGUMENTS
    cyto_plot_args <- formalArgs("cyto_plot.flowFrame")
    
    # REMOVE DISPLAY/AXES_TRANS/TITLE ARGUMENTS
    cyto_plot_args <- cyto_plot_args[-match(c("display",
                                              "axes_trans",
                                              "title"), 
                                            cyto_plot_args)]
    
    # SAMPLING
    x <- cyto_sample(x, display = display, seed = 56)
    
    # APPLY TRANSFORMATIONS
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)
    args[["x"]] <- x
    
    # PLOTTING AREA
    cyto_plot_new(popup = popup)
    n <- length(channels)
    cyto_plot_layout(
      n2mfrow(n)[1],
      n2mfrow(n)[2]
    )
    
    # PLOT PER CHANNEL
    lapply(channels, function(chan) {
      args[["channels"]] <- chan
      do.call("cyto_plot", 
              c(args[names(args) %in% cyto_plot_args],
                list("axes_trans" = transformer_list,
                     "title" = NA,
                     "display" = 1)))
    })
    
  }
  
  # Return transformerList
  return(transformer_list)
  
}

## CYTO_TRANSFORMER_ARCSINH ----------------------------------------------------

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
#' @param display number or frequency of events to display in plots, set to
#'   25000 events by default. See \code{\link{cyto_plot}} for details.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}} and
#'   \code{\link{cyto_plot}}.
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
#' @name cyto_transformer_arcsinh
NULL

#' @noRd
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
                                               popup = FALSE,
                                               display = 25000,
                                               axes_limits = "machine", ...) {

  # EXTRACT & COPY
  x <- cyto_extract(x, 
                    parent = parent, 
                    copy = TRUE)

  # SELECT
  x <- cyto_select(x, select)

  # FLOWFRAME
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_arcsinh(x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
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
                                                   popup = FALSE,
                                                   display = 25000,
                                                   axes_limits = "machine", 
                                                   ...) {

  # EXTRACT & COPY
  x <- cyto_extract(x, 
                    parent = parent,
                    copy = TRUE)

  # Call to flowFrame method
  transformer_list <- cyto_transformer_arcsinh(x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh.flowSet <- function(x,
                                           channels = NULL,
                                           parent = "root",
                                           select = NULL,
                                           plot = TRUE,
                                           popup = FALSE,
                                           display = 25000,
                                           axes_limits = "machine", ...) {

  # COPY
  x <- cyto_copy(x)
  
  # SELECT
  x <- cyto_select(x, select)

  # FLOWFRAME
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_arcsinh(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_arcsinh
#' @export
cyto_transformer_arcsinh.flowFrame <- function(x,
                                               channels = NULL,
                                               parent = "root",
                                               select = NULL,
                                               plot = TRUE,
                                               popup = FALSE,
                                               display = 25000,
                                               axes_limits = "machine", ...) {

  # PULL DOWN ALL ARGUMENTS
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  
  # COPY
  x <- cyto_copy(x)
  
  # CHANNELS
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }

  # ASINH_GML2 ARGUMENTS
  asinh_gml2_args <- formalArgs("asinh_Gml2")
  # REMOVE INVERSE ARGUMENT
  asinh_gml2_args <- asinh_gml2_args[-match("inverse", asinh_gml2_args)]
  
  # TRANSFORMATIONS
  transform_list <- lapply(channels, function(z) {
    do.call("asinh_Gml2",
            c(args[names(args) %in% asinh_gml2_args],
              list("inverse" = FALSE)))
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <- do.call("asinh_Gml2",
                   c(args[names(args) %in% asinh_gml2_args],
                     list("inverse" = TRUE)))
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
  
  # PLOTS
  if(plot == TRUE){
      
    # GRAPHICAL PARAMETERS
    old_pars <- .par("mfrow")
    on.exit(par(old_pars))
    
    # CYTO_PLOT ARGUMENTS
    cyto_plot_args <- formalArgs("cyto_plot.flowFrame")
    
    # REMOVE DISPLAY/AXES_TRANS/TITLE ARGUMENT
    cyto_plot_args <- cyto_plot_args[-match(c("display",
                                              "axes_trans",
                                              "title"), 
                                            cyto_plot_args)]
    
    # SAMPLING
    x <- cyto_sample(x, display = display, seed = 56)
    
    # Apply transformations to data for visualisation
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)
    args[["x"]] <- x

    # PLOTTING AREA
    cyto_plot_new(popup = popup)
    n <- length(channels)
    cyto_plot_layout(
      n2mfrow(n)[1],
      n2mfrow(n)[2]
    )
    
    # PLOT PER CHANNEL
    lapply(channels, function(chan) {
      args[["channels"]] <- chan
      do.call("cyto_plot", 
              c(args[names(args) %in% cyto_plot_args],
                list("axes_trans" = transformer_list,
                     "title" = NA,
                     "display" = 1)))
    })
    
  }

  # Return transformerList
  return(transformer_list)
}

## CYTO_TRANSFORMER_BIEX -------------------------------------------------------

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
#' @param display number or frequency of events to display in plots, set to
#'   25000 events by default. See \code{\link{cyto_plot}} for details.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:flowjo_biexp]{flowjo_biexp}} and
#'   \code{\link{cyto_plot}}.
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
#' @name cyto_transformer_biex
NULL

#' @noRd
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
                                          popup = FALSE,
                                          display = 25000,
                                          axes_limits = "machine", ...) {
  
  # EXTRACT & COPY
  x <- cyto_extract(x, 
                    parent = parent,
                    copy = TRUE)

  # SELECT
  x <- cyto_select(x, select)

  # FLOWFRAME
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_biex(x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
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
                                                popup = FALSE,
                                                display = 25000,
                                                axes_limits = "machine", ...) {
  
  # EXTRACT & COPY
  x <- cyto_extract(x, 
                    parent = parent,
                    copy = TRUE)

  # Call to flowFrame method
  transformer_list <- cyto_transformer_biex(x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex.flowSet <- function(x,
                                        channels = NULL,
                                        parent = "root",
                                        select = NULL,
                                        plot = TRUE,
                                        popup = FALSE,
                                        display = 25000,
                                        axes_limits = "machine", ...) {
  
  # COPY
  x <- cyto_copy(x)
  
  # SELECT
  x <- cyto_select(x, select)

  # FLOWFRAME
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_biex(x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_biex
#' @export
cyto_transformer_biex.flowFrame <- function(x,
                                          channels = NULL,
                                          parent = "root",
                                          select = NULL,
                                          plot = TRUE,
                                          popup = FALSE,
                                          display = 25000,
                                          axes_limits = "machine", ...) {
  
  # PULL DOWN ALL ARGUMENTS
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  
  # COPY
  x <- cyto_copy(x)
  
  # CHANNELS
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }

  # FLOWJO_BIEXP ARGUMENTS
  flowjo_biexp_args <- formalArgs("flowjo_biexp")
  # REMOVE INVERSE ARGUMENT
  flowjo_biexp_args <- flowjo_biexp_args[-match("inverse", flowjo_biexp_args)]
  # TRANSFORMATIONS
  transform_list <- lapply(channels, function(z) {
    do.call("flowjo_biexp",
            c(args[names(args) %in% flowjo_biexp_args],
              list("inverse" = FALSE)))
  })
  transformer_list <- lapply(transform_list, function(z) {
    inv <-  do.call("flowjo_biexp",
                    c(args[names(args) %in% flowjo_biexp_args],
                      list("inverse" = TRUE)))
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

  # PLOTS
  if(plot == TRUE){
    
    # GRAPHICAL PARAMETERS
    old_pars <- .par("mfrow")
    on.exit(par(old_pars))
    
    # CYTO_PLOT ARGUMENTS
    cyto_plot_args <- formalArgs("cyto_plot.flowFrame")
    
    # REMOVE DISPLAY ARGUMENT
    cyto_plot_args <- cyto_plot_args[-match(c("display",
                                              "axes_trans",
                                              "title"), 
                                            cyto_plot_args)]
    
    # SAMPLING
    x <- cyto_sample(x, display = display, seed = 56)
    
    # APPLY TRANSFORMATIONS
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)
    args[["x"]] <- x
     
    # PLOTTING AREA
    cyto_plot_new(popup = popup)
    n <- length(channels)
    cyto_plot_layout(
      n2mfrow(n)[1],
      n2mfrow(n)[2]
    )

    # PLOT PER CHANNEL
    lapply(channels, function(chan) {
      args[["channels"]] <- chan
      do.call("cyto_plot", 
              c(args[names(args) %in% cyto_plot_args],
                list("axes_trans" = transformer_list,
                     "title" = NA,
                     "display" = 1)))
    })
  
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
#' @param display number or frequency of events to display in plots, set to
#'   25000 events by default. See \code{\link{cyto_plot}} for details.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:estimateLogicle]{estimateLogicle}}
#'    and \code{\link{cyto_plot}}.
#'
#' @return a \code{transformerList} object.
#'
#' @importFrom flowCore inverseLogicleTransform CytoExploreR_.estimateLogicle
#' @importFrom flowWorkspace flow_trans transformerList
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowCore:logicleTransform]{estimateLogicle}}
#' @seealso
#' \code{\link[flowWorkspace:estimateLogicle]{estimateLogicle}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_combine}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @name cyto_transformer_logicle
NULL

#' @noRd
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
                                             popup = FALSE,
                                             display = 25000,
                                             axes_limits = "machine", ...) {

  # EXTRACT & COPY
  x <- cyto_extract(x, 
                    parent = parent,
                    copy = TRUE)

  # SELECT
  x <- cyto_select(x, select)

  # FLOWFRAME
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_logicle(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
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
                                                   popup = FALSE,
                                                   display = 25000,
                                                   axes_limits = "machine", 
                                                   ...) {

  # EXTRACT
  x <- cyto_extract(x, 
                    parent = parent,
                    copy = TRUE)

  # Call to flowFrame method
  transformer_list <- cyto_transformer_logicle(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle.flowSet <- function(x,
                                           channels = NULL,
                                           parent = "root",
                                           select = NULL,
                                           plot = TRUE,
                                           popup = FALSE,
                                           display = 25000,
                                           axes_limits = "machine", ...) {

  # COPY
  x <- cyto_copy(x)
  
  # SELECT
  x <- cyto_select(x, select)

  # FLOWFRAME
  x <- cyto_convert(x, "flowFrame")

  # Call to flowFrame method
  transformer_list <- cyto_transformer_logicle(
    x,
    channels = channels,
    plot = plot,
    popup = popup,
    display = display,
    axes_limits = axes_limits,
    ...
  )

  # Return transformerList
  return(transformer_list)
}

# Requires parent argument to allow passing from cyto_transform

#' @rdname cyto_transformer_logicle
#' @export
cyto_transformer_logicle.flowFrame <- function(x,
                                             channels = NULL,
                                             parent = "root",
                                             select = NULL,
                                             plot = TRUE,
                                             popup = FALSE,
                                             display = 25000,
                                             axes_limits = "machine", ...) {

  # COPY
  x <- cyto_copy(x)
  
  # PULL DOWN ALL ARGUMENTS
  args <- as.list(match.call(expand.dots = TRUE))[-1]  
  
  # CHANNELS
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)
  }
  args[["channels"]] <- channels

  # ESTIMATELOGICLE ARGUMENTS - CANNOT USE FORMALARGS (CYTOEXPLORER WRAPPER)
  estimateLogicle_args <- c("t","m","a","q")
  # TRANSFORMATIONS
  transform_list <- do.call("CytoExploreR_.estimateLogicle",
                            c(args[c("x","channels")],
                              args[names(args) %in% estimateLogicle_args]))
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

  # PLOTS
  if(plot == TRUE){
    
    # GRAPHICAL APAREMETERS
    old_pars <- .par("mfrow")
    on.exit(par(old_pars))
    
    # CYTO_PLOT ARGUMENTS
    cyto_plot_args <- formalArgs("cyto_plot.flowFrame")
    
    # REMOVE DISPLAY/AXES_TRANS/TITLE ARGUMENTS
    cyto_plot_args <- cyto_plot_args[-match(c("display",
                                              "axes_trans",
                                              "title"), 
                                            cyto_plot_args)]
    
    # SAMPLING
    x <- cyto_sample(x, display = display, seed = 56)
    
    # APPLY TRANSFORMATIONS
    transform_list <- cyto_transform_extract(transformer_list, inverse = FALSE)
    x <- transform(x, transform_list)
    args[["x"]] <- x

    # PLOTTING AREA
    cyto_plot_new(popup = popup)
    n <- length(channels)
    cyto_plot_layout(
      n2mfrow(n)[1],
      n2mfrow(n)[2]
    )

    # PLOT PER CHANNEL
    lapply(channels, function(chan) {
      args[["channels"]] <- chan
      do.call("cyto_plot", 
              c(args[names(args) %in% cyto_plot_args],
                list("axes_trans" = transformer_list,
                     "title" = NA,
                     "display" = 1)))
    })
    
  }

  # Return transformerList
  return(transformer_list)
}

## CYTO_TRANSFORMER_COMBINE ----------------------------------------------------

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
#' @importFrom methods is
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
    if(!is(axes_trans, "transformerList")){
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
    if(!is(axes_trans, "transformerList")){
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
  gh_trans <- gh_get_transformations(x,
                                     only.function = FALSE)
  if(length(gh_trans) != 0){
    gh_trans <- cyto_transformer_combine(gh_trans)
  }
  
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

## CYTO_TRANSFORMER_EXTRACT ----------------------------------------------------

#' Extract transformers from a GatingHierarchy or GatingSet
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}}
#'   or \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return transformerList or NA.
#'
#' @importFrom flowWorkspace gh_get_transformations
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' 
#' library(CytoExploreRData)
#' 
#' # Activation flowSet
#' fs <- Activation
#' 
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#' 
#' # Apply transformations
#' gs <- cyto_transform(gs, channels = c("CD4", "CD8"))
#'
#' # Extract transformations
#' trans <- cyto_transformer_extract(gs)
#'
#' @export
cyto_transformer_extract <- function(x){
  
  # GATINGSET
  if(is(x, "GatingSet")){
    x <- x[[1]]
  }
  
  # GATINGHIERARCHY
  if(is(x, "GatingHierarchy")){
    transformers <- gh_get_transformations(x, 
                                           only.function = FALSE)
    if(length(transformers) != 0){
      transformers <- cyto_transformer_combine(transformers)
    }else{
      transformers <- NA
    }
  # FLOWFRAME/FLOWSET
  }else{
    transformers <- NA
  }
  
  # RETURN TRANSFORMERLIST
  return(transformers)
}
