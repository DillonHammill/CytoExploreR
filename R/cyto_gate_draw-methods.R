# cyto_gate_draw ---------------------------------------------------------------

# NOTES ON SAMPLING:
# - flowFrame method - display applied to flowFrame and overlay

#' cyto_gate_draw
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' \code{cyto_gate_draw} is a convenient wrapper for the gating functions
#' shipped with \code{cytoRSuite} to facilitate analysis of flow cytometry by
#' gate drawing. Using \code{cyto_gate_draw} users can specify the type of
#' gate(s) to be constructed through the \code{type} argument and
#' \code{cyto_gate_draw} will automatically handle plotting the data and make
#' calls to the relevant gating function(s) to construct the gates around
#' populations of interest. \code{cyto_gate_draw} has methods for
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}} and
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} objects, refer to
#' their respective help pages for more information. The flowFrame and flowSet
#' methods simply return the constructed gates as a list of
#' \code{\link[flowCore:filters-class]{filters}}, whilst the GatingSet method
#' automatically applies the constructed gates to the GatingSet and saves the
#' constructed gates in an \code{openCyto}
#' \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}for future use.
#' See \code{\link{gate_edit}} and \code{\link{gate_remove}} to manipulate
#' constructed gates and modify their entries in the gatingTemplate.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_gate_draw,flowFrame-method}}
#' @seealso \code{\link{cyto_gate_draw,flowSet-method}}
#' @seealso \code{\link{cyto_gate_draw,GatingSet-method}}
#'
#' @export
cyto_gate_draw <- function(x, ...) {
  UseMethod("cyto_gate_draw")
}

#' cyto_gate_draw flowFrame Method.
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types include \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param overlay a \code{flowFrame}, \code{flowSet} or list of flowFrames to
#'   overlay. To improve plotting speed, overlays are subjected to the same
#'   \code{display} requirements as the underlying \code{flowFrame}.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a \code{\link[flowCore:filters-class]{filters}} list containing the
#'   drawn gate objects.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore filters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{cyto_gate_draw,flowSet-method}}
#' @seealso \code{\link{cyto_gate_draw,GatingSet-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # draw gates using cyto_gate_draw - add contour lines & overlay control
#' gt <- cyto_gate_draw(fs[[4]],
#'   channels = c("FSC-A", "SSC-A"),
#'   alias = "Cells",
#'   type = "ellipse",
#'   contour_lines = 15,
#'   overlay = fs[[1]]
#' )
#'
#' # gt is a filters object containing the contructed ellipsoidGate
#' gt[[1]]
#' }
#'
#' @export
cyto_gate_draw.flowFrame <- function(x,
                                     alias = NULL,
                                     channels = NULL,
                                     type = NULL,
                                     display = 1,
                                     overlay = NA,
                                     axis = "x",
                                     label = TRUE,
                                     plot = TRUE, ...) {

  # Default gate types
  if (is.null(type)) {
    if (length(channels) == 1) {
      if (is.null(type)) {
        type <- "interval"
      }
    } else if (length(channels) == 2) {
      if (is.null(type)) {
        type <- "polygon"
      }
    }
  }

  # Check type argument is valid
  type <- .cyto_gate_type_check(type = type, alias = alias)

  # Unsupported 1D gate types
  unspt <- c("polygon", "rectangle", "ellipse", "quadrant", "web")
  if (length(channels) == 1 &
    any(type %in% unspt)) {
    stop("Supported 1D gate types include interval, boundary and threshold.")
  }

  # Check alias is supplied correctly
  .cyto_alias_check(alias = alias, type = type)

  # Check supplied channel(s) are valid - for gating functions
  channels <- cyto_channels_extract(x,
                                    channels = channels,
                                    plot = TRUE
  )
  
  # Organise overlays to list of flowFrames
  if(!.all_na(overlay)){
    
    # flowFrame overlay added to list
    if(inherits(overlay, "flowFrame")){
      
      overlay <- list(overlay)
    
    # flowSet overlay converted to list of flowFrames  
    }else if(inherits(overlay, "flowSet")){
      
      overlay <- cyto_convert(overlay, "list of flowFrames")
    
    # flowFrame list as is - flowSet list use overlay[[1]]  
    }else if(inherits(overlay, "list")){
      
      # overlay should be list of flowFrames
      if (all(unlist(lapply(overlay, function(z) {
        inherits(z, "flowFrame")
      })))) {
        
        # overlay list of flowSets - use first fs convert to list of flowFrames
      } else if (all(unlist(lapply(overlay, function(z) {
        inherits(z, "flowSet")
      })))) {
        overlay <- overlay[[1]]
        overlay <- cyto_convert(overlay, "list of flowFrames")
        
        # overlay not supported
      } else {
        stop(paste(
          "'overlay' should be either the names of the populations to",
          "overlay, a flowFrame, a flowSet or a list of flowFrames."
        ))
      }
      
    }
    
    # Overlay sampling - list of flowFrames
    if(!is.null(display)){
      overlay <- lapply(overlay, function(z){
        cyto_sample(z, display)
      })
    }
    
  }
  
  # Sampling
  if(!is.null(display)){
    x <- cyto_sample(x, display)
  }
  
  # Generate plot - display already used - not passsed to cyto_plot
  if (plot == TRUE) {
    cyto_plot(x,
      channels = channels,
      overlay = overlay,
      popup = TRUE,
      legend = FALSE, ...
    )
  }

  # Construct gates save as filters object
  if (length(type) == 1 & type[1] == "quadrant") {
    gates <- .cyto_gate_quadrant_draw(
      x,
      channels = channels,
      alias = alias,
      plot = FALSE,
      label = label, ...
    )
  } else if (length(type) == 1 & type[1] == "web") {
    gates <- .cyto_gate_web_draw(
      x,
      channels = channels,
      alias = alias,
      plot = FALSE,
      label = label, ...
    )
  } else {
    gates <- mapply(function(type, alias) {
      if (type == "polygon") {
        .cyto_gate_polygon_draw(
          x,
          channels = channels,
          alias = alias,
          plot = FALSE,
          label = label, ...
        )
      } else if (type == "rectangle") {
        .cyto_gate_rectangle_draw(
          x,
          channels = channels,
          alias = alias,
          plot = FALSE,
          label = label, ...
        )
      } else if (type == "interval") {
        .cyto_gate_interval_draw(
          x,
          channels = channels,
          alias = alias,
          plot = FALSE,
          axis = axis,
          label = label, ...
        )
      } else if (type == "threshold") {
        .cyto_gate_threshold_draw(
          x,
          channels = channels,
          alias = alias,
          plot = FALSE,
          label = label, ...
        )
      } else if (type == "boundary") {
        .cyto_gate_boundary_draw(
          x,
          channels = channels,
          alias = alias,
          plot = FALSE,
          label = label, ...
        )
      } else if (type == "ellipse") {
        .cyto_gate_ellipse_draw(
          x,
          channels = channels,
          alias = alias,
          plot = FALSE,
          abel = label, ...
        )
      }
    }, type, alias)
  }

  gates <- filters(gates)
  return(gates)
}

#' cyto_gate_draw flowSet Method
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames} or \code{list of flowSets}
#'   containing populations to be overlaid onto the plot(s). Only overlaid
#'   flowSet objects are subjected to sampling by \code{display}.
#' @param group_by vector of pData column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to gating, set to "all" by default to construct a single gate for all
#'   samples. If group_by variables are supplied a different gate will be
#'   constructed for each group.
#' @param select designates which samples will be plotted and used for
#'   determining the best location to set the drawn gate(s). Filtering steps
#'   should be comma separated and wrapped in a list (e.g. list(Treatment ==
#'   "Stim-A", OVAConc %in% c(0,0.5))).
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowSet-method}}.
#'
#' @return a \code{\link[flowCore:filters-class]{filters}} object or list of
#'   \code{\link[flowCore:filters-class]{filters}} objects when grouping is
#'   applied.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore filters
#' @importFrom methods as
#' @importFrom purrr transpose
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot,flowSet-method}}
#' @seealso \code{\link{cyto_gate_draw,flowFrame-method}}
#' @seealso \code{\link{cyto_gate_draw,GatingSet-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # draw gates using cyto_gate_draw - add contour lines & overlay control
#' cyto_gate_draw(fs,
#'   channels = c("FSC-A", "SSC-A"),
#'   alias = "Cells",
#'   type = "polygon",
#'   contour_lines = 15,
#'   overlay = fs[[1]]
#' )
#'
#' # gt is a filters object containing the contructed polygonGate
#' gt[[1]]
#' }
#'
#' @export
cyto_gate_draw.flowSet <- function(x,
                                   alias = NULL,
                                   channels = NULL,
                                   type = NULL,
                                   display = NULL,
                                   overlay = NA,
                                   group_by = "all",
                                   select = NULL,
                                   axis = "x",
                                   label = TRUE,
                                   plot = TRUE, ...) {

  # Assign x to fs
  fs <- x  
  
  # Default gate types
  if (is.null(type)) {
    if (length(channels) == 1) {
      if (is.null(type)) {
        type <- "interval"
      }
    } else if (length(channels) == 2) {
      if (is.null(type)) {
        type <- "polygon"
      }
    }
  }

  # Check supplied channel(s) are valid - for gating functions
  channels <- cyto_channels_extract(fs,
                                    channels = channels,
                                    plot = TRUE
  )

  # Check type argument is valid
  type <- .cyto_gate_type_check(type = type, alias = alias)

  # Unsupported 1D gate types
  unspt <- c("polygon", "rectangle", "ellipse", "quadrant", "web")
  if (length(channels) == 1 &
    any(type %in% unspt)) {
    stop("Supported 1D gate types include interval, boundary and threshold.")
  }

  # Check alias is supplied correctly
  .cyto_alias_check(alias = alias, type = type)  
  
  # Group all samples together
  if (group_by[1] == "all") {

    # Select samples first based on selection criteria
    if (!is.null(select)) {
      fs <- cyto_select(fs, select)
      
    }
      
    # Restrict to 20 samples max for plotting - may need to refine over time
    if(length(fs) > 20){
      fs <- fs[sample(seq_len(length(fs)), 20)]
    }

    # Convert to merged flowFrame
    fr <- cyto_convert(fs, "flowFrame")
    
    # Sample for speed - only show mean number of events in a single sample
    if(is.null(display)){
      fr <- cyto_sample(fr, 1/length(fs))
    }else{
      fr <- cyto_sample(fr, display)
    }

    # Add merged sample to list
    fr_list <- list(fr)

    # Group samples by grouping variables
  } else {

    # Group samples by grouping variables
    fs_list <- cyto_group_by(fs, group_by)

    # Select samples per group if select supplied
    if (!is.null(select)) {
      fs_list <- lapply(fs_list, function(z) {
        # Select or return all samples if criteria not met
        tryCatch(cyto_select(z, select), error = function(e){z})
      })
    }
    
    # Maximum 20 samples per group for speed 
    fs_list <- lapply(fs_list, function(z){
      if(length(z) > 20){
        z[sample(seq_len(length(z)), 20)]
      }
      return(z)
    })

    # Merge each flowSet in fs_list
    fr_list <- lapply(fs_list, function(z) {
      # Number of samples per group
      n <- length(z)
      # Convert fs to flowFrame
      z <- cyto_convert(z, "flowFrame")
      # Sample for speed
      if(is.null(display)){
        return(cyto_sample(z, 1/n))
      }else{
        return(cyto_sample(z, display))
      }
    })
  }

  # Number of groups
  N <- length(fr_list)
  
  # Organise overlays - list of flowFrame lists of length(fr_list)
  if(!.all_na(overlay)){
    
    # flowFrame repeated in list - no sampling
    if(inherits(overlay, "flowFrame")){
      # Always show all events
      overlay <- rep(list(list(overlay)), N)
    # flowSet to lists of flowFrame lists
    }else if(inherits(overlay, "flowSet")){
      # Group by variables
      if(group_by[1] != "all"){
        # Grouping
        overlay <- cyto_group_by(overlay, group_by)
        # Coercion to list of flowFrames & sampling
        overlay <- lapply(overlay, function(z){
          # Select
          if(!is.null(select)){
            z <- tryCatch(cyto_select(z, select), error = function(e){z})
          }
          # Restrict
          if(length(z) > 20){
            z <- z[sample(seq_len(length(z)), 20)]
          }
          n <- length(z)
          y <- cyto_convert(z, "flowFrame")
          if(is.null(display)){
            return(cyto_sample(z, 1/n))
          }else{
            return(cyto_sample(y, display))
          }
        })
      # Group all samples togther
      }else{
        
        # Number of samples
        n <- length(overlay)
        
        # Select
        if(!is.null(select)){
          overlay <- cyto_select(overlay, select)
        }
        
        # Restrict
        if(length(overlay) > 20){
          overlay <- overlay[sample(seq_len(length(overlay)), 20)]
        }
        
        # Convert overlay to flowFrame
        overlay <- cyto_convert(overlay, "flowFrame")
        
        # Apply sampling
        if(is.null(display)){
          overlay <- cyto_sample(overlay, 1/n)
        }else{
          overlay <- cyto_sample(overlay, display)
        }
        overlay <- list(overlay)
      }
      # Convert list of flowFrames to list of flowFrame lists
      overlay <- lapply(overlay, function(z){
        list(z)
      })
    # Overlay is list of flowFrames or flowSets
    }else if(inherits(overlay, "list")){
      
      # List of flowFrames repeat fr_list times - no sampling or grouping
      if(all(lapply(overlay, function(z){
        inherits(z, "flowFrame")
      }))){
        
        overlay <- rep(list(overlay), N)
      
      # Allow list of flowFrame lists of length(fr_list)
      }else if (all(unlist(lapply(unlist(overlay), function(z) {
        inherits(z, "flowFrame")
      })))) {
        
        # Must be of same length as fr_list
        # No grouping, selecting or sampling - used as supplied
        if (length(overlay) != N) {
          stop(paste(
            "'overlay' must be a list of flowFrame lists -",
            "one flowFrame list per group."
          ))
        }
        
      # list of flowSets
      }else if(all(unlist(lapply(overlay, function(z){
        inherits(z, "flowSet")
      })))){
        
        # Group each flowSet, merge and sample
        overlay <- lapply(overlay, function(z){
          
          # Group by variables
          if(group_by[1] != "all"){
            # Grouping
            x <- cyto_group_by(z, group_by)
            # Coercion and sampling
            x <- lapply(x, function(y){
              # Selection
              if(!is.null(select)){
                y <- tryCatch(cyto_select(y, select), error = function(e){y})
              }
              # Restriction
              if(length(y) > 20){
                y <- y[sample(seq_len(length(y)), 20)]
              }
              n <- length(y)
              y <- cyto_convert(y, "flowFrame")
              if(is.null(display)){
                return(cyto_sample(y, 1/n))
              }else{
                return(cyto_sample(y, display))
              }
            })
          # Group all samples together
          }else{
            # Select
            if(!is.null(select)){
              z <- tryCatch(cyto_select(z, select), error = function(e){z})
            }
            # Restrict
            if(length(z) > 20){
              z <- z[sample(seq_len(length(z)), 20)]
            }
            # Coerce and sample
            n <- length(z)
            z <- cyto_convert(z, "flowFrame")
            if(is.null(display)){
              z <- cyto_sample(z, 1/n)
            }else{
              z <- cyto_sample(z, display)
            }
            z <- list(z)
            return(z)
          }
          
        })   
        
        # Overlay is a list of 
        overlay <- overlay %>% transpose()
      
      # Overlay is not supported   
      }else{
        stop(paste(
          "'overlay' should be either a flowFrame, a flowSet,",
          "list of flowFrames or a list of flowSets."
        ))
      }
    }
    
  }
  
  # Gate each group separately - list of filters
  filters_list <- lapply(seq_len(length(fr_list)), function(z){
    
    # Title 
    if(group_by[1] == "all"){
      title <- "Combined Events"
    }else{
      title <- names(fr_list)[z]
    }
    
    # Generate plot for gating - no popup or title control or overlay sampling
    if (plot == TRUE) {
      cyto_plot(fr_list[[z]],
                channels = channels,
                overlay = overlay[[z]],
                popup = TRUE,
                legend = FALSE,
                title = title, ...
      )
    }
    
    # Construct gates save as filters object
    if (length(type) == 1 & type[1] == "quadrant") {
      gates <- .cyto_gate_quadrant_draw(
        fr = fr_list[[z]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label, ...
      )
    } else if (length(type) == 1 & type[1] == "web") {
      gates <- .cyto_gate_web_draw(
        fr = fr_list[[z]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label, ...
      )
    } else {
      gates <- mapply(function(type, alias) {
        if (type == "polygon") {
          .cyto_gate_polygon_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "rectangle") {
          .cyto_gate_rectangle_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "interval") {
          .cyto_gate_interval_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            axis = axis,
            label = label, ...
          )
        } else if (type == "threshold") {
          .cyto_gate_threshold_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "boundary") {
          .cyto_gate_boundary_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "ellipse") {
          .cyto_gate_ellipse_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        }
      }, type, alias)
    }
    
  })
  names(filters_list) <- names(fr_list)

  # Combine gates in each filters_list element
  filters_list <- lapply(filters_list, function(z){
    filters(z)
  })
  
  # Group all samples returns filters object
  if(group_by[1] == "all"){
    filters_list <- filters_list[[1]]
  }
  
  return(filters_list)
}

#' cyto_gate_draw GatingSet Method
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param parent name of the \code{parent} population to extract for gating.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to be saved.
#' @param group_by vector of pData column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to gating, set to the length of x by default to construct a single
#'   gate for all samples. If group_by is supplied a different gate will be
#'   constructed for each group.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param select vector containing the indices of samples within gs to use for
#'   plotting. For large \code{flowSet} objects \code{select} is set to 20
#'   random \code{flowFrame} objects to improve processing speed.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot,GatingSet-method}}.
#'
#' @return drawn gates are applied to the
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and saved to a
#'   \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom openCyto add_pop
#' @importFrom methods as
#' @importFrom utils read.csv write.csv
#' @importFrom flowCore filters
#' @importFrom tools file_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot,GatingSet-method}}
#' @seealso \code{\link{cyto_gate_draw,flowFrame-method}}
#' @seealso \code{\link{cyto_gate_draw,flowSet-method}}
#'
#' @examples
#' \dontrun{
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
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(fs))
#' gs <- transform(gs, trans)
#'
#' # Gate using cyto_gate_draw
#' gating(Activation_gatingTemplate, gs)
#'
#' # draw gates using cyto_gate_draw - add contour lines & overlay control
#' cyto_gate_draw(gs,
#'   parent = "Dendritic Cells",
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   alias = c("CD8+ DC", "CD4+ DC"),
#'   gatingTemplate = "Example-gatingTemplate.csv",
#'   type = "rectangle",
#'   contour_lines = 15
#' )
#'
#' # Constructed gate applied directly to GatingSet
#' getNodes(gs)
#' }
#'
#' @export
cyto_gate_draw.GatingSet <- function(x,
                                     parent = "root",
                                     alias = NULL,
                                     channels = NULL,
                                     type = NULL,
                                     display = NULL,
                                     group_by = "all",
                                     select = NULL,
                                     gatingTemplate = NULL,
                                     axis = "x",
                                     label = TRUE,
                                     plot = TRUE, ...) {

  # Missing gatingtemplate
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # gatingTemplate still missing
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate csv file to save the gate(s),")
  }
  
  # File extension
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # Check gatingTemplate for existing entries
  .cyto_gatingTemplate_check(parent, alias, gatingTemplate)
  
  # Extract parent population
  fs <- cyto_extract(x, parent)

  # Default gate types
  if(length(channels) == 1){
    if(is.null(type)){
      type <- "interval"
    }
  }else if(length(channels) == 2){
   if(is.null(type)){
     type <- "polygon"
   } 
  }
  
  # Check supplied channels
  channels <- cyto_channels_extract(fs,
                                    channels = channels,
                                    plot = TRUE)
  
  # Check type argument is valid
  type <- .cyto_gate_type_check(type = type, alias = alias)
  
  # Unsupported 1D gate types
  unspt <- c("polygon","rectangle","ellipse","quadrant","web")
  if(length(channels) == 1 &
     any(type %in% unspt)){
    stop("Supported 1D gate types include interval, boundary and threshold.")
  }
  
  # Check alias is supplied correctly
  .cyto_alias_check(alias = alias, type = type)
  
  # Transformations
  axes_trans <- cyto_transform_convert(x[[1]]@transformation, inverse = FALSE)
  
  # Group all samples together
  if(group_by[1] == "all"){
    
    # Select samples first based on selection criteria
    if(!is.null(select)){
      fs <- cyto_select(fs, select)
    }
    
    # Restrict to 20 samples max for plotting - may need to refine over time
    if(length(fs) > 20){
      fs <- fs[sample(seq_len(length(fs)),20)]
    }
    
    # Convert to merged flowFrame
    fr <- cyto_convert(fs, "flowFrame")
    
    # Sample for speed - only show one sample's worth of events
    if(is.null(display)){
      fr <- cyto_sample(fr, 1/length(fs))
    }else{
      fr <- cyto_sample(fr, display)
    }
    
    # Add merged sample to list
    fr_list <- list(fr)
    
  # Group samples by grouping variables  
  }else{
    
    # Group samples by grouping variables
    fs_list <- cyto_group_by(fs, group_by)
    
    # Select samples per group if select supplied
    if(!is.null(select)){
      fs_list <- lapply(fs_list, function(z){
        # Select or return all samples if criteria not met
        tryCatch(cyto_select(z, select), error = function(e){z})
      })
    }
    
    # Restrict to maximum 20 samples per group for speed
    fs_list <- lapply(fs_list, function(z){
      if(length(z) > 20){
        z <- z[sample(seq_len(length(z)), 20)]
      }
      return(z)
    })
    
    # Merge each flowSet in fs_list
    fr_list <- lapply(fs_list, function(z){
      # Number of samples per group
      n <- length(z)
      # Convert fs to flowFrame
      z <- cyto_convert(z, "flowFrame")
      # Sample for speed
      if(is.null(display)){
        return(cyto_sample(z, 1/n))
      }else{
        return(cyto_sample(z, display))
      }
    })
    
  }
  
  # Gate each group separately - named list of filters
  filters_list <- lapply(seq_len(length(fr_list)), function(z){
    
    # Title
    if(group_by[1] == "all"){
      title <- paste("Combined Events" ,"\n", parent)
    }else{
      title <- paste(names(fr_list)[z], "\n", parent)
    }
    
    # Generate plot for gating - no popup or title control or overlay sampling
    if (plot == TRUE) {
      cyto_plot(fr_list[[z]],
                channels = channels,
                popup = TRUE,
                legend = FALSE,
                title = title,
                axes_trans = axes_trans, ...
      )
    }
    
    # Construct gates save as filters object
    if (length(type) == 1 & type[1] == "quadrant") {
      gates <- .cyto_gate_quadrant_draw(
        fr = fr_list[[z]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label, ...
      )
    } else if (length(type) == 1 & type[1] == "web") {
      gates <- .cyto_gate_web_draw(
        fr = fr_list[[z]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label, ...
      )
    } else {
      gates <- mapply(function(type, alias) {
        if (type == "polygon") {
          .cyto_gate_polygon_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "rectangle") {
          .cyto_gate_rectangle_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "interval") {
          .cyto_gate_interval_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            axis = axis,
            label = label, ...
          )
        } else if (type == "threshold") {
          .cyto_gate_threshold_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "boundary") {
          .cyto_gate_boundary_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "ellipse") {
          .cyto_gate_ellipse_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = alias,
            plot = FALSE,
            label = label, ...
          )
        }
      }, type, alias)
    }
    
  })
  names(filters_list) <- names(fr_list)
  
  # Combine gates in each list element
  filters_list <- lapply(filters_list, function(z){
    filters(z)
  })

  # format gates to be a list of alias lists
  # each of length group and appropriately named
  gates <- lapply(seq_len(length(alias)), function(y) {
    gates <- lapply(filters_list, function(x) {
      gts <- filters(list(x[[y]]))
    })
    names(gates) <- names(filters_list)
    return(gates)
  })

  # Prepare gatingTemplate entries - add negate support to cyto_plot
  pop <- "+"

  # Prepare group_by
  if (all(is.character(group_by))) {
    if(group_by[1] == "all"){
      group_by <- NA
    }else{
      group_by <- paste(group_by, collapse = ":")
    }
  } else if (all(is.na(group_by))) {
    group_by <- NA
  }

  # gatingTemplate not created yet
  if(!any(grepl(gatingTemplate, list.files()))){
    message(
      paste("Creating", gatingTemplate, "to save the constructed gate(s).")
      )
    cyto_gatingTemplate_create(gatingTemplate)
  }
  
  # Use add_pop to apply gates to GatingSet and construct gatingtemplate
  message(paste("Adding newly constructed gate(s) to", gatingTemplate, "."))
  
  gt <- read.csv(gatingTemplate, header = TRUE)
  
  pops <- list()
  for (i in seq_len(length(alias))) {
    pops[[i]] <- suppressWarnings(add_pop(
      gs = x,
      alias = alias[i],
      parent = parent,
      pop = pop,
      dims = paste(channels, collapse = ","),
      gating_method = "cyto_gate_draw",
      gating_args = list(gate = gates[[i]]),
      groupBy = group_by,
      collapseDataForGating = TRUE,
      preprocessing_method = "pp_cyto_gate_draw"
    ))
  }
  pops <- do.call("rbind", pops)
  gt <- rbind(gt, pops)
  
  write.csv(gt, gatingTemplate, row.names = FALSE)

  # Return new gatingTemplate entries
  invisible(pops)
}

#' Check Gate Type(s) Supplied to cyto_gate_draw.
#'
#' @param type vector indicating the types of gates to construct using
#'   \code{cyto_gate_draw}.
#' @param alias names of the populations to be gated.
#'
#' @return Stop gating process if type is incorrect or returns \code{type} as
#'   full lower case name(s). If a single type is supplied for multiple
#'   populations, the same type will be used for all populations.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .cyto_gate_type(type = "r", alias = c("A", "B", "C"))
#' @noRd
.cyto_gate_type_check <- function(type, alias) {
  if (all(type %in% c("q", "Q", "quadrant", "Quadrant")) & 
      length(alias) != 4) {
    stop("Supply the names of 4 poulations to alias for quadrant gates.")
  }
  
  gts <- c(
    "polygon",
    "Polygon",
    "p",
    "P",
    "rectangle",
    "Rectangle",
    "r",
    "R",
    "interval",
    "Interval",
    "i",
    "I",
    "threshold",
    "Threshold",
    "t",
    "T",
    "boundary",
    "Boundary",
    "b",
    "B",
    "ellipse",
    "Ellipse",
    "e",
    "E",
    "quadrant",
    "Quadrant",
    "q",
    "Q",
    "web",
    "Web",
    "w",
    "W"
  )
  
  if (!all(type %in% gts)) {
    if (length(type[type %in% gts == FALSE]) >= 2) {
      stop(
        paste(
          paste(type[type %in% gts == FALSE], collapse = " & "),
          "are not valid gate types for cyto_gate_draw!"
        )
      )
    } else {
      stop(paste(
        type[type %in% gts == FALSE],
        "is not a valid gate type for cyto_gate_draw!"
      ))
    }
  }
  
  type[type %in% c("polygon", "Polygon", "p", "P")] <- "polygon"
  
  type[type %in% c("rectangle", "Rectangle", "r", "R")] <- "rectangle"
  
  type[type %in% c("interval", "Interval", "i", "I")] <- "interval"
  
  type[type %in% c("threshold", "Threshold", "t", "T")] <- "threshold"
  
  type[type %in% c("boundary", "Boundary", "b", "B")] <- "boundary"
  
  type[type %in% c("ellipse", "Ellipse", "e", "E")] <- "ellipse"
  
  type[type %in% c("quadrant", "Quadrant", "q", "Q")] <- "quadrant"
  
  type[type %in% c("web", "Web", "w", "W")] <- "web"
  
  # Repeat type to equal length of alias
  if (length(type) != length(alias) &
      type[1] != "quadrant" &
      type[1] != "web") {
    type <- rep(type, length(alias))
  }
  
  return(type)
}

#' Check Alias Supplied to cyto_gate_draw
#'
#' @param alias vector indicating the names of the populations to be gated.
#' @param type vector indicating the type(s) of gate(s) to be constructed.
#'
#' @return Stops the gating process if alias is missing or of the incorrect
#'   length given the gate type.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .cyto_alias_check(alias = c("A", "B", "C", "D"), type = "q")
#' @noRd
.cyto_alias_check <- function(alias = NULL, type) {
  if (is.null(alias)) {
    stop("Supply names of populations to 'alias' for checking.")
  }
  
  if (type[1] == "quadrant" & length(alias) != 4) {
    stop("Supply 4 population names to 'alias' to construct quadrant gates.")
  }
  
  if (length(type) > 1) {
    if (length(alias) != length(type)) {
      stop("Length of alias must be the same length as type for multi-gates.")
    }
  }
}
