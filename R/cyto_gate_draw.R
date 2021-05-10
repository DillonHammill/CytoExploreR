## CYTO_GATE_DRAW --------------------------------------------------------------

#' Interactively draw and construct gate objects
#'
#' \code{cyto_gate_draw} allows to users to interactive draw gates around
#' populations which are returned as \code{flowCore} gate objects. The flowFrame
#' and flowSet methods simply return the constructed gates as a list of
#' flowCore-compatible gate objects, whilst the GatingSet method automatically
#' applies the constructed gates to the GatingSet and saves the constructed
#' gates in an \code{openCyto}
#' \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}} for future use.
#' See \code{\link{cyto_gate_edit}}, \code{\link{cyto_gate_remove}} and
#' \code{\link{cyto_gate_rename}} to manipulate constructed gates and modify
#' their entries in the gatingTemplate.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param parent name of the \code{parent} population to extract for gating when
#'   a \code{GatingSet} object is supplied.
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
#' @param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'   \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'   saved.
#' @param group_by vector of \code{\link{cyto_details}} column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to gating, set to the length of x by default to construct a single
#'   gate for all samples. If group_by is supplied a different gate will be
#'   constructed for each group. This argument only applies to \code{flowSet} or
#'   \code{GatingSet} objects.
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames} or \code{list of flowSets}
#'   containing populations to be overlaid onto the plot(s).
#' @param select designates which samples will be plotted and used for
#'   determining the best location to set the drawn gate(s). Filtering steps
#'   should be comma separated and wrapped in a list. Refer to
#'   \code{\link{cyto_select}}.
#' @param negate logical indicating whether a gatingTemplate entry should be
#'   made for the negated population (i.e. all events outside the constructed
#'   gates), set to FALSE by default. If negate is set to TRUE, a name for the
#'   negated population MUST be supplied at the end of the alias argument.
#' @param display fraction or number of events to display in the plot during the
#'   gating process, set to 25 000 events by default.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param popup logical indicating whether the plot should be constructed in a
#'   pop-up window, set to TRUE by default.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param gate_point_shape shape to use for selected gate points, set to
#'   \code{16} by default to use filled circles. See
#'   \code{\link[graphics:par]{pch}} for alternatives.
#' @param gate_point_size numeric to control the size of the selected gate
#'   points, set to 1 by default.
#' @param gate_point_col colour to use for the selected gate points, set to
#'   "red" by default.
#' @param gate_point_col_alpha numeric [0,1] to control the transparency of the
#'   selected gate points, set to 1 by default to use solid colours.
#' @param gate_line_type integer [0,6] to control the line type of gates, set to
#'   \code{1} to draw solid lines by default. See
#'   \code{\link[graphics:par]{lty}} for alternatives.
#' @param gate_line_width numeric to control the line width(s) of gates, set to
#'   \code{2.5} by default.
#' @param gate_line_col colour to use for gates, set to \code{"red"} by default.
#' @param gate_line_col_alpha numeric [0,1] to control the transparency of the
#'   selected gate lines, set to 1 by default to use solid colours.
#' @param ... additional arguments for \code{\link{cyto_plot}}.
#'
#' @return \code{cytoframe} and \code{cytoset} methods return a list of flowCore
#'   gate objects per group in group_by, whilst the \code{GatingSet} applies the
#'   constructed gates directly to the \code{GatingSet} and adds appropriate
#'   entries into the specified \code{gatingTemplate}. The \code{GatingSet}
#'   method does not return the constructed gates but instead invisibly returns
#'   the \code{gatingTemplate} entries.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom openCyto gs_add_gating_method
#' @importFrom methods as
#' @importFrom flowCore filters split
#' @importFrom flowWorkspace gs_pop_get_children gh_pop_get_descendants
#' @importFrom graphics par
#' @importFrom purrr transpose
#' @importFrom magrittr %>%
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot}}
#' @seealso \code{\link{cyto_gate_edit}}
#' @seealso \code{\link{cyto_gate_remove}}
#' @seealso \code{\link{cyto_gate_rename}}
#'
#' @examples
#' \dontrun{
#' # Gate drawing requires an interactive R session
#' library(CytoExploreRData)
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
#' gt_gating(Activation_gatingTemplate, gs)
#'
#' # draw gates using cyto_gate_draw
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
#' cyto_nodes(gs)
#' }
#'
#' @name cyto_gate_draw
NULL

#' @noRd
#' @export
cyto_gate_draw <- function(x, ...) {
  UseMethod("cyto_gate_draw")
}

#' @rdname cyto_gate_draw
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
                                     display = 25000,
                                     axis = "x",
                                     label = TRUE,
                                     plot = TRUE,
                                     popup = TRUE,
                                     axes_limits = "machine",
                                     gate_point_shape = 16,
                                     gate_point_size = 1,
                                     gate_point_col = "red",
                                     gate_point_col_alpha = 1,
                                     gate_line_type = 1,
                                     gate_line_width = 2.5,
                                     gate_line_col = "red",
                                     gate_line_col_alpha = 1, ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # ALIAS MISSING
  if (is.null(alias)) {
    stop("Supply the name(s) for the gated population(s) to 'alias'.")
  }
  
  # ACTIVE GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate csv file to save the gate(s),")
  }
  
  # EXISTING ENTRIES IN GATINGTEMPLATE
  gt <- .cyto_gatingTemplate_check(parent, alias, gatingTemplate)
  
  # GATE TYPES (REPEAT BASED ON INPUT ALIAS)
  type <- .cyto_gate_type(type, channels, alias, negate)
  
  # ALIAS (LENGTH AS EXPECTED BASED ON GATE TYPES)
  alias <- .cyto_alias(alias, type)
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels = channels, plot = TRUE)
  
  # TRANSFORMATIONS
  axes_trans <- cyto_transformers_extract(x)
  
  # NODES
  nds <- cyto_nodes(x, path = "auto")
  
  # PREPARE SAMPLES ------------------------------------------------------------
  
  # EXTRACT PARENT POPULATION
  fs <- cyto_extract(x, parent)
  
  # GROUPING (MERGE_BY)
  fr_list <- cyto_merge_by(fs,
                           merge_by = group_by,
                           select = select
  )
  
  # GROUPS
  N <- length(fr_list)
  GRPS <- names(fr_list)
  
  # PREPARE OVERLAY ------------------------------------------------------------
  
  # Organise overlays - list of flowFrame lists of length(fr_list)
  if (!.all_na(overlay)) {
    # OVERLAY - POPUALTION NAMES
    if (is.character(overlay)) {
      # OVERLAY DESCENDANTS
      if (any(grepl("descendants", overlay))) {
        overlay <- tryCatch(gh_pop_get_descendants(x[[1]],
                                                   parent,
                                                   path = "auto"
        ),
        error = function(e) {
          NA
        }
        )
        # OVERLAY CHILDREN
      } else if (any(grepl("children", overlay))) {
        overlay <- tryCatch(gs_pop_get_children(x,
                                                parent,
                                                path = "auto"
        ),
        error = function(e) {
          NA
        }
        )
      }
      # CHECK OVERLAY - MAY BE NA ABOVE
      if (!.all_na(overlay)) {
        # VALID OVERLAY
        if (all(overlay %in% nds)) {
          # EXTRACT POPULATIONS - LIST OF FLOWSETS
          nms <- overlay
          overlay <- lapply(overlay, function(z) {
            cyto_extract(x, z)
          })
          names(overlay) <- nms
        }
      }
    }
    # OVERLAY - FLOWFRAME
    if (is(overlay, "flowFrame")) {
      # Always show all events
      overlay <- rep(list(list(overlay)), N)
      # flowSet to lists of flowFrame lists
    } else if (is(overlay, "flowSet")) {
      # GROUPING (MERGE_BY) - LIST OF FLOWFRAMES
      overlay <- cyto_merge_by(overlay,
                               merge_by = group_by,
                               select = select
      )
      # FLOWFRAME LIST TO LIST OF FLOWFRAME LISTS
      overlay <- lapply(overlay, function(z) {
        list(z)
      })
      # OVERLAY - LIST OF FLOWFRAMES OR FLOWSETS
    } else if (is(overlay, "list")) {
      # LIST OF FLOWFRAMES - REPEAT FR_LIST TIMES
      if (all(LAPPLY(overlay, function(z) {
        is(z, "flowFrame")
      }))) {
        overlay <- rep(list(overlay), N)
        # LIST FLOWFRAME LISTS OF LENGTH FR_LIST
      } else if (all(LAPPLY(unlist(overlay), function(z) {
        is(z, "flowFrame")
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
      } else if (all(LAPPLY(overlay, function(z) {
        is(z, "flowSet")
      }))) {
        # GROUP & MERGE EACH FLOWSET
        overlay <- lapply(overlay, function(z) {
          # GROUPING (MERGE_BY)
          cyto_merge_by(z,
                        merge_by = group_by,
                        select = select
          )
        })
        # OVERLAY TRANSPOSE
        overlay <- overlay %>% transpose()
        # OVERLAY NOT SUPPORTED
      } else {
        stop(paste(
          "'overlay' should be either a flowFrame, a flowSet,",
          "list of flowFrames or a list of flowSets."
        ))
      }
    }
  } else {
    overlay <- rep(list(list(NA)), N)
  }
  
  # COMBINE FR_LIST & OVERLAY
  fr_list <- lapply(seq_along(fr_list), function(z) {
    return(c(fr_list[z], overlay[[z]]))
  })
  names(fr_list) <- GRPS
  
  # SAMPLING - APPLY SAME SAMPLING PER LAYER
  fr_list <- lapply(seq_along(fr_list), function(z) {
    lapply(seq_along(fr_list[[z]]), function(y){
      # WATCH OUT NA OVERLAY
      if(!is.logical(fr_list[[z]][[y]])){
        cyto_sample(fr_list[[z]][[y]],
                    display = display,
                    seed = 56,
                    plot = FALSE)
      }else{
        return(fr_list[[z]][[y]])
      }
    })
  })
  names(fr_list) <- GRPS
  
  # GATING ---------------------------------------------------------------------
  
  # CONSTRUCT GATES PER GROUP
  gate_list <- lapply(seq_along(fr_list), function(z) {
    
    # CYTO_PLOT - PROPER TITLES
    if (plot == TRUE) {
      
      # PREPARE TITLE - KEEP VALID PARENT FOR GATINGTEMPLATE
      if (parent == "root") {
        prnt <- "All Events"
      } else {
        prnt <- parent
      }
      
      # TITLE
      title <- paste(GRPS[z], "\n", prnt)
      
      # PLOT
      cyto_plot(fr_list[[z]][[1]],
                channels = channels,
                overlay = fr_list[[z]][seq_along(fr_list[[z]])[-1]],
                display = 1,
                popup = popup,
                legend = FALSE,
                title = title,
                axes_trans = axes_trans,
                axes_limits = axes_limits, ...
      )
    }
    
    # CYTO_GATE_DRAW FLOWFRAME METHOD
    cyto_gate_draw(fr_list[[z]][[1]],
                   alias = unlist(alias),
                   channels = channels,
                   type = type,
                   negate = negate,
                   display = 1,
                   axis = axis,
                   label = label,
                   plot = FALSE,
                   gate_point_shape = gate_point_shape,
                   gate_point_size = gate_point_size,
                   gate_point_col = gate_point_col,
                   gate_point_col_alpha = gate_point_col_alpha,
                   gate_line_type = gate_line_type,
                   gate_line_width = gate_line_width,
                   gate_line_col = gate_line_col,
                   gate_line_col_alpha = gate_line_col_alpha,
    )
  })
  names(gate_list) <- names(fr_list)
  
  # TRANSPOSE GATE_LIST - LIST LENGTH ALIAS - EACH LENGTH GROUP
  gate_list <- transpose(gate_list)
  
  # LIST OF FILTERS OF LENGTH ALIAS - NAMES IMPORTANT!
  gate_list <- lapply(seq_along(gate_list), function(z) {
    gates <- lapply(seq_along(gate_list[[z]]), function(y) {
      if (!is(gate_list[[z]][[y]], "quadGate")) {
        filters(gate_list[[z]][y])
      }else{
        gate_list[[z]][y]
      }
    })
    names(gates) <- GRPS
    return(gates)
  })
  
  # GATE NAMES - IMPORTANT
  if(negate == TRUE){
    names(gate_list) <- unlist(alias)[-length(alias)]
  }else{
    if(all(type %in% "quadrant")){
      names(gate_list) <- paste(unlist(alias), collapse = "|")
    }else{
      names(gate_list) <- unlist(alias)
    }
  }
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # GROUP_BY
  if (all(is.character(group_by))) {
    if (group_by[1] == "all") {
      group_by <- "NA"
    } else {
      group_by <- paste(group_by, collapse = ":")
    }
  } else if (all(is.na(group_by))) {
    group_by <- "NA"
  }
  
  # GATINGTEMPLATE NOT CREATED YET
  if (is.null(gt)) {
    message(
      paste("Creating", gatingTemplate, "to save the constructed gate(s).")
    )
    cyto_gatingTemplate_create(gatingTemplate)
  }
  
  # gs_add_gating_method - GATINGTEMPLATE ENTRY & APPLY TO GATINGSET
  message(paste("Adding newly constructed gate(s) to", gatingTemplate, "."))
  
  # READ IN GATINGTEMPLATE
  if(is.null(gt)){
    gt <- read_from_csv(gatingTemplate)
  }
  
  # GATED POPULATIONS
  pops <- lapply(seq_len(length(alias[!is.na(type)])), function(z) {
    # PREPARE POP
    if (length(alias[[z]]) == 1 | type[1] == "web") {
      pop <- "+"
    } else {
      pop <- "*"
    }
    # QUADRANT GATES - NOT FILTERS
    if (type[z] == "quadrant") {
      gate_list[[z]] <- unlist(gate_list[[z]])
    }
    # GATED POPULATIONS
    suppressWarnings(
      gs_add_gating_method(
        gs = x,
        alias = paste(alias[[z]], collapse = ","),
        parent = parent,
        pop = pop,
        dims = paste(channels, collapse = ","),
        gating_method = "cyto_gate_draw",
        gating_args = list(gate = gate_list[[z]],
                           openCyto.minEvents = -1),
        groupBy = group_by,
        collapseDataForGating = TRUE,
        preprocessing_method = "pp_cyto_gate_draw"
      )
    )
  })
  
  # NEGATED POPULATIONS
  if (negate == TRUE) {
    pops[[length(pops) + 1]] <- suppressMessages(
      gs_add_gating_method(
        gs = x,
        alias = alias[[which(is.na(type))]],
        parent = parent,
        pop = "+",
        dims = paste(channels, collapse = ","),
        gating_method = "boolGate",
        gating_args = paste(paste0("!", unlist(alias[which(!is.na(type))])),
                            collapse = "&!"
        ),
        groupBy = group_by,
        collapseDataForGating = TRUE,
        preprocessing_method = NA
      )
    )
  }
  
  # RBIND GATINGTEMPLATE ENTRIES
  pops <- do.call("rbind", pops)
  
  # COMBINE NEW ENTRIES WITH EXISTING ONES
  gt <- rbind(gt, pops)
  
  # SAVE UPDATED GATINGTEMPLATE
  write_to_csv(gt, gatingTemplate)
  
  # RETURN GATINGTEMPLATE ENTRIES
  invisible(pops)
}

#' @rdname cyto_gate_draw
#' @export
cyto_gate_draw.flowSet <- function(x,
                                   alias = NULL,
                                   channels = NULL,
                                   type = NULL,
                                   overlay = NA,
                                   group_by = "all",
                                   select = NULL,
                                   negate = FALSE,
                                   display = 25000,
                                   axis = "x",
                                   label = TRUE,
                                   plot = TRUE,
                                   popup = TRUE,
                                   axes_limits = "machine",
                                   gate_point_shape = 16,
                                   gate_point_size = 1,
                                   gate_point_col = "red",
                                   gate_point_col_alpha = 1,
                                   gate_line_type = 1,
                                   gate_line_width = 2.5,
                                   gate_line_col = "red",
                                   gate_line_col_alpha = 1, ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # ALIAS MISSING
  if (is.null(alias)) {
    stop("Supply the name(s) for the gated population(s) to 'alias'.")
  }
  
  # ASSIGN X TO FS
  fs <- x
  
  # PREPARE SAMPLES ------------------------------------------------------------
  
  # GROUPING (MERGE BY)
  fr_list <- cyto_merge_by(fs,
                           merge_by = group_by,
                           select = select)
  
  # GROUPS
  N <- length(fr_list)
  GRPS <- names(fr_list)
  
  # PREPARE OVERLAY ------------------------------------------------------------
  
  # Organise overlays - list of flowFrame lists of length(fr_list)
  if (!.all_na(overlay)) {
    # OVERLAY - FLOWFRAME
    if (is(overlay, "flowFrame")) {
      # Always show all events
      overlay <- rep(list(list(overlay)), N)
      # flowSet to lists of flowFrame lists
    } else if (is(overlay, "flowSet")) {
      # GROUPING (MERGE BY) - LIST OF FLOWFRAMES
      overlay <- cyto_merge_by(overlay,
                               merge_by = group_by,
                               select = select
      )
      # FLOWFRAME LIST TO LIST OF FLOWFRAME LISTS
      overlay <- lapply(overlay, function(z) {
        list(z)
      })
      # OVERLAY - LIST OF FLOWFRAMES OR FLOWSETS
    } else if (is(overlay, "list")) {
      # LIST OF FLOWFRAMES - REPEAT FR_LIST TIMES
      if (all(LAPPLY(overlay, function(z) {
        is(z, "flowFrame")
      }))) {
        overlay <- rep(list(overlay), N)
        # LIST FLOWFRAME LISTS OF LENGTH FR_LIST
      } else if (all(LAPPLY(unlist(overlay), function(z) {
        is(z, "flowFrame")
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
      } else if (all(LAPPLY(overlay, function(z) {
        is(z, "flowSet")
      }))) {
        # GROUP & MERGE EACH FLOWSET
        overlay <- lapply(overlay, function(z) {
          # GROUPING (MERGE BY)
          cyto_merge_by(z,
                        merge_by = group_by,
                        select = select
          )
        })
        # OVERLAY TRANSPOSE
        overlay <- overlay %>% transpose()
        # OVERLAY NOT SUPPORTED
      } else {
        stop(paste(
          "'overlay' should be either a flowFrame, a flowSet,",
          "list of flowFrames or a list of flowSets."
        ))
      }
    }
  } else {
    overlay <- rep(list(list(NA)), N)
  }
  
  # COMBINE FR_LIST & OVERLAY
  fr_list <- lapply(seq_along(fr_list), function(z) {
    return(c(fr_list[z], overlay[[z]]))
  })
  names(fr_list) <- GRPS
  
  # GATING ---------------------------------------------------------------------
  
  # CONSTRUCT GATES PER GROUP
  gate_list <- lapply(seq_along(fr_list), function(z) {
    
    # CYTO_GATE_DRAW FLOWFRAME METHOD
    cyto_gate_draw(fr_list[[z]][[1]],
                   alias = alias,
                   channels = channels,
                   type = type,
                   overlay = fr_list[[z]][seq_along(fr_list[[z]])[-1]],
                   negate = negate,
                   display = display,
                   axis = axis,
                   label = label,
                   plot = plot,
                   popup = popup,
                   axes_limits = axes_limits,
                   gate_point_shape = gate_point_shape,
                   gate_point_size = gate_point_size,
                   gate_point_col = gate_point_col,
                   gate_point_col_alpha = gate_point_col_alpha,
                   gate_line_type = gate_line_type,
                   gate_line_width = gate_line_width,
                   gate_line_col = gate_line_col,
                   gate_line_col_alpha = gate_line_col_alpha,
                   ...
    )
  })
  names(gate_list) <- names(fr_list)
  
  # RETURN GATES
  return(gate_list)
}

#' @rdname cyto_gate_draw
#' @export
cyto_gate_draw.flowFrame <- function(x,
                                     alias = NULL,
                                     channels = NULL,
                                     type = NULL,
                                     overlay = NA,
                                     negate = FALSE,
                                     display = 25000,
                                     axis = "x",
                                     label = TRUE,
                                     plot = TRUE,
                                     popup = TRUE,
                                     axes_limits = "machine",
                                     gate_point_shape = 16,
                                     gate_point_size = 1,
                                     gate_point_col = "red",
                                     gate_point_col_alpha = 1,
                                     gate_line_type = 1,
                                     gate_line_width = 2.5,
                                     gate_line_col = "red",
                                     gate_line_col_alpha = 1, ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # ALIAS MISSING
  if (is.null(alias)) {
    stop("Supply the name(s) for the gated population(s) to 'alias'.")
  }
  
  # GATE TYPE
  type <- .cyto_gate_type(type, channels, alias, negate)
  
  # ALIAS
  alias <- .cyto_alias(alias, type)
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels = channels, plot = TRUE)
  
  # PREPARE SAMPLES & OVERLAY --------------------------------------------------
  
  # X FLOWFRAME LIST
  fr_list <- list(x)
  
  # OVERLAY FLOWFRAME LIST
  if (!.all_na(overlay)) {
    # OVERLAY FLOWFRAME
    if (is(overlay, "flowFrame")) {
      overlay_list <- list(overlay)
      # OVERLAY FLOWSET
    } else if (is(overlay, "flowSet")) {
      overlay_list <- cyto_convert(overlay, "list of flowFrames")
      # OVERLAY FLOWFRAME LIST
    } else if (is(overlay, "list")) {
      # overlay should be list of flowFrames
      if (all(LAPPLY(overlay, function(z) {
        is(z, "flowFrame")
      }))) {
        overlay_list <- overlay
        # OVERLAY LIST OF FLOWSETS
      } else if (all(LAPPLY(overlay, function(z) {
        is(z, "flowSet")
      }))) {
        overlay <- overlay[[1]]
        overlay_list <- cyto_convert(overlay, "list of flowFrames")
        # OVERLAY NOT SUPPORTED
      } else {
        stop(paste(
          "'overlay' should be either the names of the populations to",
          "overlay, a flowFrame, a flowSet or a list of flowFrames."
        ))
      }
    }
  } else {
    overlay_list <- rep(list(NA), length(fr_list))
  }
  
  # COMBINE FR_LIST & OVERLAY
  fr_list <- c(fr_list, overlay_list)
  
  # SAMPLING - SAME SAMPLING PER LAYER
  fr_list <- lapply(seq_along(fr_list), function(z){
    # WATCH OUT NA OVERLAY
    if(!is.logical(fr_list[[z]])){
      cyto_sample(fr_list[[z]],
                  display = display,
                  seed = 56,
                  plot = FALSE)
    }else{
      return(fr_list[[z]])
    }
  })
  
  # CONSTRUCT PLOT -------------------------------------------------------------
  
  # CYTO_PLOT
  if (plot == TRUE) {
    cyto_plot(fr_list[[1]],
              channels = channels,
              overlay = fr_list[seq(2, length(fr_list))],
              display = 1,
              popup = popup,
              legend = FALSE,
              axes_limits = axes_limits,
              ...
    )
  }
  
  # GATING ---------------------------------------------------------------------
  
  # CONSTRUCT GATE OBJECTS
  gates <- mapply(function(type,
                           alias) {
    if (type == "polygon") {
      .cyto_gate_polygon_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "rectangle") {
      .cyto_gate_rectangle_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "interval") {
      .cyto_gate_interval_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        axis = axis,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "threshold") {
      .cyto_gate_threshold_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "boundary") {
      .cyto_gate_boundary_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "ellipse") {
      .cyto_gate_ellipse_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "quadrant") {
      .cyto_gate_quadrant_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    } else if (type == "web") {
      .cyto_gate_web_draw(
        fr_list[[1]],
        channels = channels,
        alias = alias,
        plot = FALSE,
        label = label,
        gate_point_shape = gate_point_shape,
        gate_point_size = gate_point_size,
        gate_point_col = gate_point_col,
        gate_point_col_alpha = gate_point_col_alpha,
        gate_line_type = gate_line_type,
        gate_line_width = gate_line_width,
        gate_line_col = gate_line_col,
        gate_line_col_alpha = gate_line_col_alpha
      )
    }
  }, type[!is.na(type)], alias[!is.na(type)], USE.NAMES = FALSE)
  
  # LABEL NEGATED POPULATION
  if (negate == TRUE) {
    # GATE
    if (length(unlist(gates)) == 1) {
      NP <- split(fr_list[[1]], unlist(gates)[[1]])[[2]]
    } else {
      negate_filter <- do.call("|", unname(unlist(gates)))
      NP <- split(fr_list[[1]], negate_filter)[[2]]
    }
    # NEGATED LABEL XCOORD
    negate_text_x <- suppressMessages(
      # WATCH FOR EMPTY FLOWFRAMES
      if (cyto_stat_count(NP) < 2) {
        mean(par("usr")[c(1, 2)])
      } else {
        cyto_apply(NP, 
                   "cyto_stat_mode", 
                   channels = channels[1],
                   input = "matrix",
                   inverse = FALSE,
                   copy = FALSE)[, 1]
      }
    )
    # NEGATED LABEL YCOORD
    if (length(channels) == 1) {
      negate_text_y <- mean(par("usr")[c(3, 4)])
    } else if (length(channels) == 2) {
      negate_text_y <- suppressMessages(
        # WATCH OUT FOR EMPTY FLOWFRAMES
        if (cyto_stat_count(NP) < 2) {
          mean(par("usr")[c(3, 4)])
        } else {
          cyto_apply(NP, 
                     "cyto_stat_mode", 
                     channels = channels[2],
                     input = "matrix",
                     inverse = FALSE,
                     copy = FALSE)[, 1]
        }
      )
    }
    # NEGATED STATISTIC
    negate_stat <- cyto_stat_count(NP) / cyto_stat_count(fr_list[[1]]) * 100
    negate_stat <- paste(.round(negate_stat, 2), "%")
    # NEGATED LABEL
    cyto_plot_labeller(
      label_text = paste(alias[length(alias)],
                         negate_stat,
                         sep = "\n"
      ),
      label_text_x = negate_text_x,
      label_text_y = negate_text_y,
      label_text_size = 1
    )
  }
  
  # RETURN GATE OBJECTS --------------------------------------------------------
  
  # GATE OBJECTS LIST - EXCLUDE NEGATED ENTRIES
  gates <- unlist(gates)
  return(gates)
}
