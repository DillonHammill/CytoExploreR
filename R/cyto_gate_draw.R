## CYTO_GATE_DRAW --------------------------------------------------------------

#'Interactively draw and construct gate objects
#'
#'\code{cyto_gate_draw} allows to users to interactive draw gates around
#'populations which are returned as \code{flowCore} gate objects. The flowSet
#'methods simply return the constructed gates as a list of flowCore-compatible
#'gate objects, whilst the GatingSet method automatically applies the
#'constructed gates to the GatingSet and saves the constructed gates in an
#'\code{openCyto} \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}
#'for future use. See \code{\link{cyto_gate_edit}},
#'\code{\link{cyto_gate_copy}}, \code{\link{cyto_gate_remove}} and
#'\code{\link{cyto_gate_rename}} to manipulate constructed gates and modify
#'their entries in the gatingTemplate.
#'
#'@param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'  \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'  \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'@param channels vector of channel names to use for plotting, can be of length
#'  1 for 1-D density histogram or length 2 for 2-D scatter plot.
#'@param parent name of the \code{parent} population to extract for gating when
#'  a \code{GatingSet} object is supplied.
#'@param alias the name(s) of the populations to be gated. If multiple
#'  population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates will
#'  be returned. \code{alias} is \code{NULL} by default which will halt the
#'  gating routine.
#'@param type vector of gate type names used to construct the gates. Multiple
#'  gate types are supported but should be accompanied with an \code{alias}
#'  argument of the same length (i.e. one \code{type} per \code{alias}).
#'  Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'  boundary, interval, quadrant and web} which can be abbreviated as upper or
#'  lower case first letters as well. Default \code{type} is \code{"interval"}
#'  for 1D gates and \code{"polygon"} for 2D gates.
#'@param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'  \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'  saved.
#'@param merge_by vector of \code{\link{cyto_details}} column names (e.g.
#'  c("Treatment","Concentration") indicating how the samples should be grouped
#'  prior to gating, set to "all" by default to apply the same gate(s) to all
#'  samples. If group_by is supplied a different gate will be constructed for
#'  each merged group.
#'@param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'  \code{flowSet}, \code{list of flowFrames} or \code{list of flowSets}
#'  containing populations to be overlaid onto the plot(s).
#'@param select designates which samples will be plotted and used for
#'  determining the best location to set the drawn gate(s). Filtering steps
#'  should be comma separated and wrapped in a list. Refer to
#'  \code{\link{cyto_select}}.
#'@param negate logical indicating whether a gatingTemplate entry should be made
#'  for the negated population (i.e. all events outside the constructed gates),
#'  set to FALSE by default. If negate is set to TRUE, a name for the negated
#'  population MUST be supplied at the end of the alias argument.
#'@param display fraction or number of events to display in the plot during the
#'  gating process, set to 25 000 events by default.
#'@param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'  gated for 2-D interval gates.
#'@param label logical indicating whether to include
#'  \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'  default.
#'@param plot logical indicating whether a plot should be drawn, set to
#'  \code{TRUE} by default.
#'@param title title to use for the plot, set to the name of the sample and
#'  parent population by default. Title can be removed by setting this argument
#'  to \code{NA}.
#'@param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied data. \code{cyto_plot} does
#'   not support in-line transformations and as such the transformations should
#'   be applied to the data prior to plotting. The transformerList is used
#'   internally to ensure that the axes on the constructed plots are
#'   appropriately labelled.
#'@param axes_limits options include \code{"auto"}, \code{"data"} or
#'  \code{"machine"} to use optimised, data or machine limits respectively. Set
#'  to \code{"machine"} by default to use entire axes ranges. Fine control over
#'  axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'  arguments.
#'@param gate_point_shape shape to use for selected gate points, set to
#'  \code{16} by default to use filled circles. See
#'  \code{\link[graphics:par]{pch}} for alternatives.
#'@param gate_point_size numeric to control the size of the selected gate
#'  points, set to 1 by default.
#'@param gate_point_col colour to use for the selected gate points, set to "red"
#'  by default.
#'@param gate_point_col_alpha numeric [0,1] to control the transparency of the
#'  selected gate points, set to 1 by default to use solid colours.
#'@param gate_line_type integer [0,6] to control the line type of gates, set to
#'  \code{1} to draw solid lines by default. See \code{\link[graphics:par]{lty}}
#'  for alternatives.
#'@param gate_line_width numeric to control the line width(s) of gates, set to
#'  \code{2.5} by default.
#'@param gate_line_col colour to use for gates, set to \code{"red"} by default.
#'@param gate_line_col_alpha numeric [0,1] to control the transparency of the
#'  selected gate lines, set to 1 by default to use solid colours.
#'@param label_text_font numeric to control the font of text in plot labels, set
#'  to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'  alternatives.
#'@param label_text_size numeric to control the size of text in the plot labels,
#'  set to 1 by default.
#'@param label_text_col colour(s) to use for text in plot labels, set to
#'  \code{"black"} by default.
#'@param label_text_col_alpha numeric [0, 1] to control the transparency of the
#'  text colour, set to 1 by default to remove transparency.
#'@param label_fill fill colour(s) to use for labels, set to "white" by default.
#'@param label_fill_alpha numeric to control background fill transparency of
#'  label, set to 0.6 by default to introduce some transparency.
#'@param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'  sampling is applied with each \code{\link{cyto_plot}} call, set to an
#'  arbitrary numeric by default. This behaviour can be turned off by setting
#'  this argument to NULL.
#'@param ... additional arguments for \code{\link{cyto_plot}}.
#'
#'@return \code{cytoframe} and \code{cytoset} methods return a list of flowCore
#'  gate objects per group in group_by, whilst the \code{GatingSet} applies the
#'  constructed gates directly to the \code{GatingSet} and adds appropriate
#'  entries into the specified \code{gatingTemplate}. The \code{GatingSet}
#'  method does not return the constructed gates but instead invisibly returns
#'  the \code{gatingTemplate} entries.
#'
#'@importFrom BiocGenerics colnames
#'@importFrom openCyto gs_add_gating_method
#'@importFrom methods as
#'@importFrom flowCore filters split
#'@importFrom flowWorkspace gs_pop_get_children gh_pop_get_descendants
#'@importFrom graphics par
#'@importFrom purrr transpose
#'@importFrom methods is
#'
#'@author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#'@seealso \code{\link{cyto_plot}}
#'@seealso \code{\link{cyto_gate_edit}}
#'@seealso \code{\link{cyto_gate_remove}}
#'@seealso \code{\link{cyto_gate_rename}}
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
#' gs <- compensate(gs)
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
#'@export
cyto_gate_draw <- function(x,
                           parent = "root",
                           alias = NULL,
                           channels = NULL,
                           type = NULL,
                           gatingTemplate = NULL,
                           overlay = NA,
                           merge_by = "all",
                           select = NULL,
                           negate = FALSE,
                           display = 50000,
                           axis = "x",
                           label = TRUE,
                           plot = TRUE,
                           title = "",
                           axes_trans = NA,
                           axes_limits = "machine",
                           gate_point_shape = 16,
                           gate_point_size = 1,
                           gate_point_col = "red",
                           gate_point_col_alpha = 1,
                           gate_line_type = 1,
                           gate_line_width = 2.5,
                           gate_line_col = "red",
                           gate_line_col_alpha = 1,
                           label_text_size = 1,
                           label_text_font = 2,
                           label_text_col = "black",
                           label_text_col_alpha = 1,
                           label_fill = "white",
                           label_fill_alpha = 0.75,
                           seed = 42,
                           ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # FLAG - SKIP DATA PREPARTION IN CYTO_PLOT
  cyto_option("cyto_plot_data", TRUE)
  on.exit({
    cyto_option("cyto_plot_data", FALSE)
  })
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
  # GATINGTEMPLATE CHECKS
  if(cyto_class(x, "GatingSet")) {
    # ACTIVE GATINGTEMPLATE
    if (is.null(gatingTemplate)) {
      gatingTemplate <- cyto_gatingTemplate_active()
    }
    # MISSING GATINGTEMPLATE
    if (is.null(gatingTemplate)) {
      stop(
        "Supply the name of the gatingTemplate csv file to save the gate(s)."
        )
    }
    # CHECK EXISTING ENTRIES IN GATINGTEMPLATE
    gt <- .cyto_gatingTemplate_check(parent, 
                                     alias, 
                                     gatingTemplate)
    # CREATE GATINGTEMPLATE
    if (is.null(gt)) {
      message(
        paste("Creating", gatingTemplate, "to save the constructed gate(s).")
      )
      cyto_gatingTemplate_create(gatingTemplate) # ACTIVE?
      gt <- cyto_gatingTemplate_read(gatingTemplate, data.table = TRUE)
    }
  }
  
  # TYPE -> LIST
  type <- .cyto_gate_type(type = type,
                          channels = channels,
                          alias = alias,
                          negate = negate)
  
  # ALIAS MISSING
  if(is.null(alias)) {
    stop("Supply the name(s) for the gated population(s) to 'alias'.")
  # ALIAS -> LIST
  } else {
    alias <- .cyto_gate_alias(alias, type)
  }
  
  # TRANSFORMATIONS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }

  # PREPARE DATA ---------------------------------------------------------------
  
  # LIST OF CYTOSET LISTS
  cs_lists <- .cyto_plot_data(x,
                              parent = parent,
                              overlay = overlay,
                              merge_by = merge_by,
                              select = select,
                              display = display,
                              seed = seed)

  # GATING ---------------------------------------------------------------------
  
  # TITLE
  title <- rep(c(title, rep("", length(cs_lists))), 
               length.out = length(cs_lists))
  
  # LOOP THROUGH GROUPS
  gate_list <- structure(
    lapply(seq_along(cs_lists), function(z){
      # CYTOSET LIST
      cs_list <- cs_lists[[z]]
      # TITLE
      if(.empty(title[z])) {
        # GROUP
        title[z] <<- names(cs_lists)[z]
        # PARENT POPULATION
        if(!is.null(names(cs_lists[[z]])[1])) {
          title[z] <<- paste(title[z], 
                             names(cs_lists[[z]])[1], 
                             sep = "\n")
        }
      }
      # PLOT
      if(plot == TRUE) {
        # MERGED CYTOSET
        cyto_plot(cs_list[[1]],
                  channels = channels,
                  overlay = ifelse(
                    length(cs_list) > 1,
                    cs_list[seq_along(cs_list)[-1]],
                    NA
                  ),
                  legend = FALSE,
                  axes_limits = "machine",
                  title = title[z],
                  ...)
      }
      # GATE
      gt <- list()
      structure(
        mapply(function(alias,
                        type, 
                        axis,
                        gate_point_shape,
                        gate_point_size,
                        gate_point_col,
                        gate_point_col_alpha,
                        gate_line_type,
                        gate_line_width,
                        gate_line_col,
                        gate_line_col_alpha,
                        label,
                        label_text_size,
                        label_text_font,
                        label_text_col,
                        label_text_col_alpha,
                        label_fill,
                        label_fill_alpha){
          # DATA TO GATE
          x <- cs_list[[1]]
          # ARGUMENTS
          args <- .args_list()
          args$channels <- channels
          # GATE
          if(!.all_na(type)) {
            gate <- .cyto_gate_draw_dispatch(args) # ... not used
            # NEGATE
          } else {
            if(length(gt) == 1) {
              gate <- !gt[[1]]
            } else {
              gate <- !do.call("|", unname(unlist(gt)))
            }
          }
          # ADD GATE
          args$gate <- gate
          gt <<- c(gt, structure(
            list(gate),
            names = paste(alias, collapse = "|") # quad gates
          ))
          # LABELS - GATE BASE LAYER
          do.call(".cyto_gate_draw_label", args) # ... not used
        },
        alias,
        type, 
        axis,
        gate_point_shape,
        gate_point_size,
        gate_point_col,
        gate_point_col_alpha,
        gate_line_type,
        gate_line_width,
        gate_line_col,
        gate_line_col_alpha,
        label,
        label_text_size,
        label_text_font,
        label_text_col,
        label_text_col_alpha,
        label_fill,
        label_fill_alpha),
        names = alias
      )
      return(gt)
    }),
    names = names(cs_lists)
  )
  
  # CYTOSET METHOD RETURNS GATES
  if(!cyto_class(x, "GatingSet")) {
    return(gate_list)
  }
  
  # PREPARE GATES --------------------------------------------------------------
  
  # TRANSPOSE GATES - LIST OF LENGTH ALIAS - EACH OF LENGTH GROUP
  gate_list <- transpose(gate_list)
  
  # WRAP GATES AS FILTERS LIST
  gate_list <- structure(
    lapply(seq_along(gate_list), function(z){
      gates <- structure(
        lapply(seq_along(gate_list[[z]]), function(y){
          if(!cyto_class(gate_list[[z]][[y]], "quadGate")) {
            filters(gate_list[[z]][y])
          } else {
            gate_list[[z]][y]
          }
        }),
        names = names(gate_list[[z]])
      )
    }), 
    names = names(gate_list)
  )
  
  # REMOVE NEGATED GATE FROM LIST
  if(negate == TRUE) {
    gate_list <- gate_list[!names(gate_list) %in% alias[[length(alias)]]]
  }
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # PREPARE GROUP_BY
  if(all(is.character(merge_by))) {
    if(all(merge_by == "all")) {
      group_by <- "NA"
    } else {
      group_by <- paste(merge_by, collapse = ":")
    }
  } else {
    group_by <- "NA"
  }
  
  # gs_add_gating_method - GATINGTEMPLATE ENTRY & APPLY TO GATINGSET
  message(paste("Adding newly constructed gate(s) to", gatingTemplate, "."))
  
  # ADD GATED POPULATIONS TO GATINGTEMPLATE
  pops <- lapply(seq_along(alias[!is.na(type)]), function(z) {
    # PREPARE POP
    if (length(alias[[z]]) == 1 | type[[z]] == "web") {
      pop <- "+"
    } else {
      pop <- "*"
    }
    # QUADRANT GATES - NOT FILTERS
    if (type[[z]] == "quadrant") {
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

  # ADD NEGATED POPULATIONS TO GATINGTEMPLATE
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
  cyto_gatingTemplate_write(gt, gatingTemplate)
  
  # RETURN GATINGTEMPLATE ENTRIES
  return(x)
  
}