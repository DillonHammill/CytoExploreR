# CYTO_GATE_DRAW ---------------------------------------------------------------

#' cyto_gate_draw
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' \code{cyto_gate_draw} is a convenient wrapper for the gating functions
#' shipped with \code{CytoExploreR} to facilitate analysis of flow cytometry by
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
#' See \code{\link{cyto_gate_edit}} and \code{\link{cyto_gate_remove}} to
#' manipulate constructed gates and modify their entries in the gatingTemplate.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.flowSet}}
#' @seealso \code{\link{cyto_gate_draw.GatingSet}}
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
#' @seealso \code{\link{cyto_plot.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.flowSet}}
#' @seealso \code{\link{cyto_gate_draw.GatingSet}}
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
                                     overlay = NA,
                                     axis = "x",
                                     label = TRUE,
                                     plot = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # ALIAS
  alias <- .cyto_alias(alias, type)  
  
  # GATE TYPE
  type <- .cyto_gate_type(type, channels, alias)

  # CHANNELS
  channels <- cyto_channels_extract(x,
                                    channels = channels,
                                    plot = TRUE
  )
  
  # CONSTRUCT PLOT -------------------------------------------------------------
  
  # CYTO_PLOT
  if (plot == TRUE) {
    cyto_plot(x,
      channels = channels,
      overlay = overlay,
      popup = TRUE,
      legend = FALSE, ...
    )
  }

  # GATING ---------------------------------------------------------------------
  
  # CONSTRUCT GATE OBJECTS
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
          label = label, ...
        )
      }
    }, type, alias)
  }

  # RETURN GATE OBJECTS --------------------------------------------------------
  
  # GATES AS FILTERS OBJECTS
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
#' @param overlay  a \code{flowFrame}, \code{flowSet}, \code{list of flowFrames}
#'   or \code{list of flowSets} containing populations to be overlaid onto the
#'   plot(s). Only overlaid flowSet objects are subjected to sampling by
#'   \code{display}.
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
#' @seealso \code{\link{cyto_plot.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.GatingSet}}
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
                                   overlay = NA,
                                   group_by = "all",
                                   select = NULL,
                                   axis = "x",
                                   label = TRUE,
                                   plot = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # ASSIGN X TO FS
  fs <- x  
  
  # ALIAS
  alias <- .cyto_alias(alias, type) 
  
  # GATE TYPES
  type <- .cyto_gate_type(type, channels, alias)
  
  # Check supplied channel(s) are valid - for gating functions
  channels <- cyto_channels_extract(fs,
                                    channels = channels,
                                    plot = TRUE
  )

  # PREPARE SAMPLES ------------------------------------------------------------
  
  # GROUP ALL
  if (group_by[1] == "all") {
    # SELECT
    if (!is.null(select)) {
      fs <- cyto_select(fs, select)
    }
    # MERGED FLOWFRAME
    fr <- cyto_convert(fs, "flowFrame")
    # FLOWFRAME METHOD
    fr_list <- list(fr)
  # GROUP variables
  } else {
    # GROUPING
    fs_list <- cyto_group_by(fs, group_by)
    # SELECT PER GROUP
    if (!is.null(select)) {
      fs_list <- lapply(fs_list, function(z) {
        # Select or return all samples if criteria not met
        tryCatch(cyto_select(z, select), error = function(e){z})
      })
    }
    # MERGE EACH FLOWSET
    fr_list <- lapply(fs_list, function(z) {
      # Number of samples per group
      n <- length(z)
      # Convert fs to flowFrame
      z <- cyto_convert(z, "flowFrame")
      return(z)
    })
  }
  
  # GROUPS
  N <- length(fr_list)
  
  # PREPARE OVERLAY ------------------------------------------------------------
  
  # Organise overlays - list of flowFrame lists of length(fr_list)
  if(!.all_na(overlay)){
    # OVERLAY - FLOWFRAME
    if(inherits(overlay, "flowFrame")){
      # Always show all events
      overlay <- rep(list(list(overlay)), N)
    # flowSet to lists of flowFrame lists
    }else if(inherits(overlay, "flowSet")){
      # GROUP VARIABLES
      if(group_by[1] != "all"){
        # GROUPING
        overlay <- cyto_group_by(overlay, group_by)
        # LIST OF FLOWFRAMES
        overlay <- lapply(overlay, function(z){
          # SELECT
          if(!is.null(select)){
            z <- tryCatch(cyto_select(z, select), error = function(e){z})
          }
          # CONVERT
          z <- cyto_convert(z, "flowFrame")
          return(z)
        })
      # GROUP ALL
      }else{
        # SELECT
        if(!is.null(select)){
          overlay <- cyto_select(overlay, select)
        }
        # FLOWFRAME
        overlay <- cyto_convert(overlay, "flowFrame")
        # FLOWFRAME LIST
        overlay <- list(overlay)
      }
      # FLOWFRAME LIST TO LIST OF FLOWFRAME LISTS
      overlay <- lapply(overlay, function(z){
        list(z)
      })
    # OVERLAY - LIST OF FLOWFRAMES OR FLOWSETS
    }else if(inherits(overlay, "list")){
      # LIST OF FLOWFRAMES - REPEAT FR_LIST TIMES
      if(all(LAPPLY(overlay, function(z){
        inherits(z, "flowFrame")
      }))){
        overlay <- rep(list(overlay), N)
      # LIST FLOWFRAME LISTS OF LENGTH FR_LIST
      }else if (all(LAPPLY(unlist(overlay), function(z) {
        inherits(z, "flowFrame")
      }))) {
        # Must be of same length as fr_list
        # No grouping, selecting or sampling - used as supplied
        if (length(overlay) != N) {
          stop(paste(
            "'overlay' must be a list of flowFrame lists -",
            "one flowFrame list per group."
          ))
        }
      # LIST OF FLOWSETS
      }else if(all(LAPPLY(overlay, function(z){
        inherits(z, "flowSet")
      }))){
        # GROUP & MERGE EACH FLOWSET
        overlay <- lapply(overlay, function(z){
          # GROUP VARIABLES
          if(group_by[1] != "all"){
            # GROUPING
            x <- cyto_group_by(z, group_by)
            # Coercion and sampling
            x <- lapply(x, function(y){
              # SELECT
              if(!is.null(select)){
                y <- tryCatch(cyto_select(y, select), error = function(e){y})
              }
              # CONVERT
              y <- cyto_convert(y, "flowFrame")
              return(y)
            })
          # GROUP ALL
          }else{
            # SELECT
            if(!is.null(select)){
              z <- tryCatch(cyto_select(z, select), error = function(e){z})
            }
            # CONVERT
            z <- cyto_convert(z, "flowFrame")
            return(z)
          }
        })   
        # OVERLAY TRANSPOSE
        overlay <- overlay %>% transpose()
      # OVERLAY NOT SUPPORTED  
      }else{
        stop(paste(
          "'overlay' should be either a flowFrame, a flowSet,",
          "list of flowFrames or a list of flowSets."
        ))
      }
    }
  }
  
  # GATING ---------------------------------------------------------------------
  
  # GATE EACH GROUP - NAMED LIST OF FILTERS
  filters_list <- lapply(seq_len(length(fr_list)), function(z){
    # Title 
    if(group_by[1] == "all"){
      title <- "Combined Events"
    }else{
      title <- names(fr_list)[z]
    }
    # CONSTRUCT PLOT
    if (plot == TRUE) {
      if(!.all_na(overlay)){
        cyto_plot(fr_list[[z]],
                channels = channels,
                overlay = overlay[[z]],
                popup = TRUE,
                legend = FALSE,
                title = title, ...
        )
      }else{
        cyto_plot(fr_list[[z]],
                  channels = channels,
                  overlay = NA,
                  popup = TRUE,
                  legend = FALSE,
                  title = title, ...
        )
      }
    }
    
    # CONSTRUCT GATES
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

  # COMBINE GATES IN EACH LIST ELEMENT
  filters_list <- lapply(filters_list, function(z){
    filters(z)
  })
  
  # ALL GROUPED RETURN FILTERS OBJECT
  if(group_by[1] == "all"){
    filters_list <- filters_list[[1]]
  }
  
  # RETURN GATES
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
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames} or \code{list of flowSets}
#'   containing populations to be overlaid onto the plot(s). Only overlaid
#'   flowSet objects are subjected to sampling by \code{display}.
#' @param select vector containing the indices of samples within gs to use for
#'   plotting. For large \code{flowSet} objects \code{select} is set to 20
#'   random \code{flowFrame} objects to improve processing speed.
#' @param negate logical indicating whether a gatingTemplate entry should be
#'   made for the negated population (i.e. all events outside the constructed
#'   gates), set to FALSE by default. If negate is set to TRUE, a name for the
#'   negated population MUST be supplied at the end of the alias argument.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot.flowFrame}}.
#'
#' @return drawn gates are applied to the
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and saved to a
#'   \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom openCyto add_pop
#' @importFrom methods as
#' @importFrom utils read.csv write.csv
#' @importFrom flowCore filters split
#' @importFrom tools file_ext
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.flowFrame}}
#' @seealso \code{\link{cyto_gate_draw.flowSet}}
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
                                     gatingTemplate = NULL,
                                     overlay = NA,
                                     group_by = "all",
                                     select = NULL,
                                     negate = FALSE,
                                     axis = "x",
                                     label = TRUE,
                                     plot = TRUE, ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # ACTIVE GATINGTEMPLATE
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # MISSING GATINGTEMPLATE
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate csv file to save the gate(s),")
  }
  
  # MISSING GATINGTEMPLATE FILE EXTENSION
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # EXISTING ENTRIES IN GATINGTEMPLATE
  .cyto_gatingTemplate_check(parent, alias, gatingTemplate)
  
  # ALIAS
  alias <- .cyto_alias(alias, type, negate)
  
  # GATE TYPES
  type <- .cyto_gate_type_check(type, channels, alias, negate)
  
  # NEGATED ALIAS
  if(negate == TRUE){
    ALIAS <- alias[-length(alias)]
  }else{
    ALIAS <- alias
  }
  
  # CHANNELS
  channels <- cyto_channels_extract(fs,
                                    channels = channels,
                                    plot = TRUE)
  
  # TRANSFORMATIONS
  axes_trans <- x[[1]]@transformation
  
  # PREPARE SAMPLES ------------------------------------------------------------
  
  # EXTRACT PARENT POPULATION
  fs <- cyto_extract(x, parent)

  # GROUP ALL
  if (group_by[1] == "all") {
    # SELECT
    if (!is.null(select)) {
      fs <- cyto_select(fs, select)
    }
    # MERGED FLOWFRAME
    fr <- cyto_convert(fs, "flowFrame")
    # FLOWFRAME METHOD
    fr_list <- list(fr)
    # GROUP variables
  } else {
    # GROUPING
    fs_list <- cyto_group_by(fs, group_by)
    # SELECT PER GROUP
    if (!is.null(select)) {
      fs_list <- lapply(fs_list, function(z) {
        # Select or return all samples if criteria not met
        tryCatch(cyto_select(z, select), error = function(e){z})
      })
    }
    # MERGE EACH FLOWSET
    fr_list <- lapply(fs_list, function(z) {
      # Number of samples per group
      n <- length(z)
      # Convert fs to flowFrame
      z <- cyto_convert(z, "flowFrame")
      return(z)
    })
  }
  
  # GROUPS
  N <- length(fr_list)
  
  # PREPARE OVERLAY ------------------------------------------------------------
  
  # Organise overlays - list of flowFrame lists of length(fr_list)
  if(!.all_na(overlay)){
    # OVERLAY - POPUALTION NAMES
    if (is.character(overlay)) {
      # VALID OVERLAY
      if (all(overlay %in% basename(cyto_nodes(x)))) {
        # EXTRACT POPULATIONS
        nms <- overlay
        overlay <- lapply(overlay, function(z) {
          cyto_extract(x, z)
        })
        names(overlay) <- nms
      }
    }
    # OVERLAY - FLOWFRAME
    if(inherits(overlay, "flowFrame")){
      # Always show all events
      overlay <- rep(list(list(overlay)), N)
      # flowSet to lists of flowFrame lists
    }else if(inherits(overlay, "flowSet")){
      # GROUP VARIABLES
      if(group_by[1] != "all"){
        # GROUPING
        overlay <- cyto_group_by(overlay, group_by)
        # LIST OF FLOWFRAMES
        overlay <- lapply(overlay, function(z){
          # SELECT
          if(!is.null(select)){
            z <- tryCatch(cyto_select(z, select), error = function(e){z})
          }
          # CONVERT
          z <- cyto_convert(z, "flowFrame")
          return(z)
        })
        # GROUP ALL
      }else{
        # SELECT
        if(!is.null(select)){
          overlay <- cyto_select(overlay, select)
        }
        # FLOWFRAME
        overlay <- cyto_convert(overlay, "flowFrame")
        # FLOWFRAME LIST
        overlay <- list(overlay)
      }
      # FLOWFRAME LIST TO LIST OF FLOWFRAME LISTS
      overlay <- lapply(overlay, function(z){
        list(z)
      })
      # OVERLAY - LIST OF FLOWFRAMES OR FLOWSETS
    }else if(inherits(overlay, "list")){
      # LIST OF FLOWFRAMES - REPEAT FR_LIST TIMES
      if(all(LAPPLY(overlay, function(z){
        inherits(z, "flowFrame")
      }))){
        overlay <- rep(list(overlay), N)
        # LIST FLOWFRAME LISTS OF LENGTH FR_LIST
      }else if (all(LAPPLY(unlist(overlay), function(z) {
        inherits(z, "flowFrame")
      }))) {
        # Must be of same length as fr_list
        # No grouping, selecting or sampling - used as supplied
        if (length(overlay) != N) {
          stop(paste(
            "'overlay' must be a list of flowFrame lists -",
            "one flowFrame list per group."
          ))
        }
        # LIST OF FLOWSETS
      }else if(all(LAPPLY(overlay, function(z){
        inherits(z, "flowSet")
      }))){
        # GROUP & MERGE EACH FLOWSET
        overlay <- lapply(overlay, function(z){
          # GROUP VARIABLES
          if(group_by[1] != "all"){
            # GROUPING
            x <- cyto_group_by(z, group_by)
            # Coercion and sampling
            x <- lapply(x, function(y){
              # SELECT
              if(!is.null(select)){
                y <- tryCatch(cyto_select(y, select), error = function(e){y})
              }
              # CONVERT
              y <- cyto_convert(y, "flowFrame")
              return(y)
            })
            # GROUP ALL
          }else{
            # SELECT
            if(!is.null(select)){
              z <- tryCatch(cyto_select(z, select), error = function(e){z})
            }
            # CONVERT
            z <- cyto_convert(z, "flowFrame")
            return(z)
          }
        })   
        # OVERLAY TRANSPOSE
        overlay <- overlay %>% transpose()
        # OVERLAY NOT SUPPORTED  
      }else{
        stop(paste(
          "'overlay' should be either a flowFrame, a flowSet,",
          "list of flowFrames or a list of flowSets."
        ))
      }
    }
  }
  
  # GATING ---------------------------------------------------------------------
  
  # GATE EACH GROUP - NAMED LIST OF FILTERS
  filters_list <- lapply(seq_len(length(fr_list)), function(z){
    # TITLE
    if(group_by[1] == "all"){
      title <- paste("Combined Events" ,"\n", parent)
    }else{
      title <- paste(names(fr_list)[z], "\n", parent)
    }
    # CONSTRUCT PLOT
    if (plot == TRUE) {
      if(!.all_na(overlay)){
        cyto_plot(fr_list[[z]],
                channels = channels,
                overlay = overlay[[z]],
                popup = TRUE,
                legend = FALSE,
                title = title,
                axes_trans = axes_trans, ...
        )
      }else{
        cyto_plot(fr_list[[z]],
                  channels = channels,
                  overlay = NA,
                  popup = TRUE,
                  legend = FALSE,
                  title = title,
                  axes_trans = axes_trans, ...
        )
      }
    }
    # CONSTRUCT GATES
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
      gates <- mapply(function(type, ALIAS) {
        if (type == "polygon") {
          .cyto_gate_polygon_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = ALIAS,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "rectangle") {
          .cyto_gate_rectangle_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = ALIAS,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "interval") {
          .cyto_gate_interval_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = ALIAS,
            plot = FALSE,
            axis = axis,
            label = label, ...
          )
        } else if (type == "threshold") {
          .cyto_gate_threshold_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = ALIAS,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "boundary") {
          .cyto_gate_boundary_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = ALIAS,
            plot = FALSE,
            label = label, ...
          )
        } else if (type == "ellipse") {
          .cyto_gate_ellipse_draw(
            fr = fr_list[[z]],
            channels = channels,
            alias = ALIAS,
            plot = FALSE,
            label = label, ...
          )
        }
      }, type, ALIAS)
      # NEGATED POPULATION
      negate_filter <- do.call("|", unlist(gates))
      NP <- split(fr_list[[z]], negate_filter)[[2]]
      # NEGATED LABEL XCOORD
      negate_text_x <- suppressMessages(
        .cyto_mode(NP, channels = channels[1])
      )
      # NEGATED LABEL YCOORD
      if(length(channels) == 1){
        negate_text_x <- mean(par("usr")[c(3,4)])
      }else if(length(channels) == 2){
        negate_text_y <- suppressMessages(
          .cyto_mode(NP, channels = channels[2])
        )
      }
      # NEGATED STATISTIC
      negate_stat <- .cyto_count(NP)/.cyto_count(fr_list[[z]]) * 100
      negate_stat <- paste(.round(negate_stat, 2), "%")
      # NEGATED LABEL
      cyto_plot_labeller(label_text = paste(alias[length(alias)],
                                            negate_stat,
                                            sep = "\n"),
                         label_text_x = negate_text_x,
                         label_text_y = negate_text_y,
                         label_text_size = 1)
    }
    
  })
  names(filters_list) <- names(fr_list)
  
  # COMBINE GATES IN EACH LIST ELEMENT
  filters_list <- lapply(filters_list, function(z){
    filters(z)
  })

  # FORMAT GATES - LIST OF ALIAS LISTS - EACH LENGTH GROUP & NAMED
  gates <- lapply(seq_len(length(alias)), function(y) {
    gates <- lapply(filters_list, function(x) {
      gts <- filters(list(x[[y]]))
    })
    names(gates) <- names(filters_list)
    return(gates)
  })

  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # POP
  pop <- "+"

  # GROUP_BY
  if (all(is.character(group_by))) {
    if(group_by[1] == "all"){
      group_by <- NA
    }else{
      group_by <- paste(group_by, collapse = ":")
    }
  } else if (all(is.na(group_by))) {
    group_by <- NA
  }

  # GATINGTEMPLATE NOT CREATED YET
  if(!any(grepl(gatingTemplate, list.files()))){
    message(
      paste("Creating", gatingTemplate, "to save the constructed gate(s).")
      )
    cyto_gatingTemplate_create(gatingTemplate)
  }
  
  # ADD_POP - GATINGTEMPLATE ENTRY & APPLY TO GATINGSET
  message(paste("Adding newly constructed gate(s) to", gatingTemplate, "."))
  
  # READ IN GATINGTEMPLATE
  gt <- read.csv(gatingTemplate, header = TRUE)
  
  # NO NEGATED POPULATIONS
  if(negate == FALSE){
    # GATED POPULATIONS
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
  # NEGATED POPULATIONS
  }else if(negate == TRUE){
    # GATED POPULATIONS
    pops <- list()
    for (i in seq_len(length(ALIAS))) {
      pops[[i]] <- suppressWarnings(add_pop(
        gs = x,
        alias = ALIAS[i],
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
    # NEGATED POPULATION
    pops[[length(pops) + 1]] <- suppressWarnings(
      add_pop(
        gs =x,
        alias = alias[length(alias)],
        parent = parent,
        pop = pop,
        dims = paste(channels, collapse = ","),
        gating_method = "boolGate",
        gating_args = paste(paste0("!", ALIAS), collapse = "&"),
        groupBy = group_by,
        collapseDataForGating = TRUE,
        preprocessing_method = NA
      )
    )
    pops <- do.call("rbind", pops)
  }
  # COMBINE ROWS
  gt <- rbind(gt, pops)
  
  # SAVE UPDATED GATINGTEMPLATE
  write.csv(gt, gatingTemplate, row.names = FALSE)

  # RETURN GATINGTEMPLATE ENTRIES
  invisible(pops)
}
