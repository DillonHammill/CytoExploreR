#' @export
cyto_gate_edit <- function(x,
                           parent = NULL,
                           alias = NULL,
                           channels = NULL,
                           type = NULL,
                           gatingTemplate = NULL,
                           overlay = NA,
                           merge_by = "all",
                           select = NULL,
                           negate = FALSE,
                           logic = NULL,
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
                           seed = 42,
                           ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # FLAG - SKIP DATA PREPARTION IN CYTO_PLOT
  cyto_option("cyto_plot_data", TRUE)
  on.exit({
    cyto_option("cyto_plot_data", FALSE)
  })
  
  # PARENT
  if(is.null(parent)) {
    stop("Supply the name of the parental population to 'parent'.")
  } else {
    parent <- cyto_nodes_convert(x, 
                                 nodes = parent,
                                 path = "auto")
  }
  
  # ALIAS
  if(is.null(alias)) {
    stop("Supply the name(s) of the gates to edit to 'alias'.")
  } else {
    alias <- cyto_nodes_convert(x, 
                                 nodes = alias,
                                 path = "auto")
  }
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to edit gate(s).")
  }
  
  # AXES_TRANS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # PREPARE GATINGTEMPLATE -----------------------------------------------------
  
  # READ GATINGTEMPLATE
  gtf <- cyto_gatingTemplate_read(gatingTemplate,
                                  data.table = FALSE)
  
  # REPLACE PARENT WITH AUTO PATHS (IN CASE)
  gtf$parent <- cyto_nodes_convert(x, 
                                   nodes = gtf$parent,
                                   path = "auto")
  
  # RESTRICT ENTRIES TO PARENT
  gtf_parent <- gtf[gtf$parent == parent, ]
  
  # RESTRICT BY PARENT & ALIAS
  gtf_alias <- lapply(gtf_parent$alias, function(z) {
    unlist(strsplit(as.character(z), ","))
  })
  ind <- LAPPLY(gtf_alias, function(z) {
    if (any(alias %in% z)) {
      TRUE
    } else {
      FALSE
    }
  })
  # UPDATE gtf_alias
  gtf_alias <- gtf_alias[ind]
  
  # PULL OUT RELEVENT GATINTEMPLATE ENTRIES
  gtf_chunk <- gtf_parent[ind, ]

  # EXTRACT CHANNELS
  if (is.null(channels)) {
    channels <- unique(
      unlist(
        strsplit(
          as.character(gtf_chunk[gtf_chunk$alias %in% alias, "dims"]), ",",
          fixed = TRUE
        )
      )
    )
  }
  
  # EXTRACT CHANNELS - REQUIRED FOR QUADGATES
  channels <- cyto_channels_extract(gs, channels)
  
  # EXTRACT GROUPBY
  gtf_groupBy <- unique(gtf_chunk[, "groupBy"])
  
  # EXTRACT GATING METHOD
  gtf_gating_method <- unique(as.character(gtf_chunk[, "gating_method"]))
  
  # UNSUPPORTED GATING METHOD
  if (any(gtf_gating_method %in% c("refGate", "boolGate"))) {
    stop(
      "Use cyto_gatingTemplate_edit() to modify refGate and boolGate entries."
      )
  }
  
  # REMOVE BOOLEAN GATES
  bool_alias <- gtf_chunk$alias[gtf_chunk$gating_method == "boolGate"]
  if(length(bool_alias) > 0) {
    alias <- alias[!alias %in% bool_alias]
  }
  
  # REFERENCE GATES
  if(any(gtf_chunk$gating_method == "refGate")) {
    message(
      "Reference gate(s) will be replaced with their own gate(s)."
    )
  }
  
  # EXTRACT GATES FROM GATINGSET DIRECTLY --------------------------------------
  
  # GATES PER GROUP
  gates_gs <- cyto_gate_extract(x,
                                parent = parent,
                                alias = alias,
                                select = select,
                                merge_by = merge_by)
  
  # PREPARE ALIAS & TYPE -------------------------------------------------------
  
  # GET GATE TYPE FROM EXISTING GATES
  if (is.null(type)) {
    type <- cyto_gate_type(gates_gs[[1]])
  }
  
  # TYPE
  type <- .cyto_gate_type(type, channels, alias, negate = FALSE)
  
  # ALIAS - LIST PER GATE TYPE
  alias <- .cyto_gate_alias(alias, type)
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # LIST OF CYTOSET LISTS
  cs_lists <- .cyto_plot_data(x,
                              parent = parent,
                              overlay = overlay,
                              merge_by = merge_by,
                              select = select,
                              display = display,
                              seed = seed)
  
  # SELECT GROUP(S) TO EDIT ----------------------------------------------------
  
  # MENU - SELECT GROUPS TO EDIT
  if(merge_by[1] == "all") {
    grps <- "all"
  # SKIP MENU IF NEW GROUPING ADDED
  } else if(.all_na(gtf_groupBy) & merge_by[1] != "all") {
    grps <- names(cs_lists)
  # SKIP MENU IF GROUP HAS CHANGED
  } else if(!.all_na(gtf_groupBy) &
            !setequal(group_by, unlist(strsplit(gtf_groupBy, ":")))) {
    grps <- names(cs_lists)
  # GROUPING REMAINS UNCHANGED
  } else {
    message("Select the group(s) to edit:")
    # INTERACTIVE GROUP SELECTION
    if(interactive()) {
      grps <- data.frame("group" = names(cs_lists),
                         "select" = NA,
                         stringsAsFactors = FALSE)
      grps <- data_edit(grps,
                        title = "Group Selection",
                        logo = CytoExploreR_logo(),
                        col_edit = FALSE,
                        row_edit = FALSE,
                        col_options = list("select" = c(TRUE, FALSE)),
                        colnames = c("select", "group"),
                        col_readonly = "group",
                        viewer = "pane")
      grps <- grps[, "group"][which(grps[, "select"] == 1)]
      # NA TO CHARACTERS
      if(any(is.na(grps))){
        grps[is.na(grps)] <- "NA"
      }
    # OTHERWISE ALL GROUPS
    } else {
      grps <- names(cs_lists)
    }
  }
  
  # GROUPS
  cyto_groups <- cyto_groups()
  
  # # EXPERIMENT DETAILS
  # pd <- cyto_details(x)
  # 
  # # GROUPING
  # if (group_by[1] == "all") {
  #   pd$groupby <- rep("all", nrow(pd))
  # } else {
  #   pd$groupby <- do.call(paste, pd[, group_by, drop = FALSE])
  # }
  
  # GATE EDITING ---------------------------------------------------------------
  
  # TITLE
  if(plot){
    title <- rep(c(title, rep("", length(grps))), 
                 length.out = length(grps))
  }

  # LOOP THROUGH GROUPS
  w <- 0
  lapply(match(grps, names(cs_lists)), function(y) {
    w <<- w + 1
    # CYTOSET LIST
    cs_list <- cs_lists[[y]]
    # EXISTING GATE LIST
    gate <- gates_gs[[y]]
    print(gate)
    # CYTO_PLOT
    if(plot == TRUE) {
      # TITLE
      if(.empty(title[w])) {
        # GROUP
        title[w] <<- names(cs_lists)[y]
        # PARENT POPULATION
        if(!is.null(names(cs_lists[[y]])[1])) {
          title[w] <<- paste(title[w], 
                             names(cs_lists[[y]])[1], 
                             sep = "\n")
        }
      }
      # CYTO_PLOT
      cyto_plot(cs_list[[1]],
                overlay = ifelse(length(cs_list) > 1,
                                 cs_list[seq_along(cs_list)[-1]],
                                 NA),
                channels = channels,
                axes_trans = axes_trans,
                axes_limits = axes_limits,
                legend = FALSE,
                gate = gate,
                gate_line_width = 2.5,
                gate_line_col = "magenta",
                gate_line_alpha = 0.8,
                label = FALSE,
                title = title,
                ...)
    }
    
    # 2D INTERVAL - AXIS ARGUMENT
    if("interval" %in% type) {
      intvl <- rbind(
        gate[[match("interval", type)[1]]]@min,
        gate[[match("interval", type)[1]]]@max
      )
      if (all(is.finite(intvl[, 1]))) {
        axis <- "x"
      } else if (all(is.finite(intvl[, 2]))) {
        axis <- "y"
      }
    }
    
    # DRAW NEW GATES
    gate_new <- cyto_gate_draw(cs_list[[1]],
                               alias = unlist(alias),
                               channels = channels,
                               type = type,
                               negate = FALSE,
                               axis = axis,
                               plot = FALSE,
                               gate_point_shape = gate_point_shape,
                               gate_point_size = gate_point_size,
                               gate_point_col = gate_point_col,
                               gate_point_col_alpha = gate_point_col_alpha,
                               gate_line_type = gate_line_type,
                               gate_line_width = gate_line_width,
                               gate_line_col = gate_line_col,
                               gate_line_col_alpha = gate_line_col_alpha)
    
    # NAME GATES
    names(gate_new) <- LAPPLY(alias, function(z) {
      paste(z, collapse = ",")
    })
    
    gates_gs[[y]] <<- gate_new
    gates_gT[[y]] <<- gate_new
  })
  
  # PREPARE GATES FOR GATINGTEMPLATE
  gates_gT_transposed <- transpose(gates_gT)
  
  # WRAP GATES IN FILTER LIST
  gates_gT_transposed <- structure(
    lapply(seq_along(gates_gT_transposed), function(z){
      gates <- structure(
        lapply(seq_along(gates_gT_transposed[[z]]), function(y){
          if(!cyto_class(gates_gT_transposed[[z]][[y]], "quadGate")) {
            filters(gates_gT_transposed[[z]][y])
          } else {
            gates_gT_transposed[[z]][y]
          }
        }),
        names = names(gates_gT_transposed[[z]])
      )
    }), 
    names = names(gates_gT_transposed)
  )
  
  # DATA.TABLE FRIENDLY NAMES
  prnt <- parent
  als <- alias
  gtmd <- "cyto_gate_draw"
  ppmd <- "pp_cyto_gate_draw"
  
  # REMOVE NEGATED POPULATIONS FROM ALS (NEGATED ENTRY CANNOT BE MODIFIED)
  als <- lapply(als, function(z) {
    paste(z, collapse = ",")
  })
  
  # REMOVE NEGATED REFERENCE
  if (negate == TRUE) {
    als <- als[-length(als)]
  }
  
  # FIND & EDIT GATINGTEMPLATE ENTRIES
  gt <- cyto_gatingTemplate_read(gatingTemplate,
                                 data.table = TRUE)
  
  # DATA.TABLE R CMD CHECK NOTE
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  preprocessing_args <- NULL
  .SD <- NULL
  
  # GROUP_BY
  if (group_by[1] == "all") {
    group_by <- "NA"
  } else {
    group_by <- paste(group_by, collapse = ":")
  }
  
  # MODIFY GATINGTEMPLATE - GATED POPULATIONS ONLY (IGNORE NEGATED ENTRIES)
  for (i in seq_len(length(als))) {
    gt[parent == prnt & alias == als[i], gating_method := gtmd]
    gt[
      parent == prnt & alias == als[i],
      gating_args := CytoExploreR_.argDeparser(list(
        gate = gates_gT_transposed[[i]],
        openCyto.minEvents = -1
      ))
    ]
    gt[parent == prnt & alias == als[i], collapseDataForGating := TRUE]
    # groupBy must be character class
    gt[, groupBy := lapply(.SD, as.character), .SDcols = "groupBy"]
    gt[parent == prnt & alias == als[i], groupBy := group_by]
    gt[parent == prnt & alias == als[i], preprocessing_method := ppmd]
    gt[parent == prnt & alias == als[i], preprocessing_args := as.logical(NA)]
  }
  
  # UPDATE BOOLEAN LOGIC
  if(length(bool_alias) > 0 & !is.null(logic)) {
    
  }
  
  
  # SAVE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(gt,
                            save_as = gatingTemplate)
  
  # CONVERT ALIAS BACK TO VECTOR
  alias <- as.character(unlist(alias))
  
  # CONVERT QUAD TO RECTANGLE FOR GATINGSET SAVING (EASIEST WAY)
  gates_gs <- lapply(gates_gs, function(z){
    if(cyto_class(z[[1]], "quadGate")){
      .cyto_gate_quad_convert(z[[1]], channels)
    }else{
      z
    }
  })
  
  # APPLY NEW GATES TO GATINGSET (INCLUDES NEGATED ALIAS)
  lapply(grps, function(z) {
    # MATCHING GROUPS
    ind <- which(pd$groupby == z)
    lapply(seq_along(alias), function(y) {
      # SET GATES - GATED POPULATIONS ONLY
      if (y < length(alias) | (negate == FALSE & y == length(alias))) {
        # SET UPDATED GATES- IGNORE NEGATED GATE
        # REPEAT GATES FOR EACH SAMPLE IN GROUP
        gates_update <- rep(
          list(gates_gs[[z]][[alias[y]]]),
          length(gs[ind])
        )
        names(gates_update) <- cyto_names(gs[ind])
        suppressMessages(gs_pop_set_gate(
          gs[ind],
          paste(parent, alias[y], sep = "/"),
          gates_update
        ))
      }
      # RECOMPUTE STATISTICS (INCLUDE NEGATED POPULATION)
      suppressMessages(recompute(gs[ind], paste(parent, alias[y], sep = "/")))
    })
  })
  
  # UPDATE GATINGSET GLOBALLY
  assign(deparse(substitute(x)), gs, envir = globalenv())
  
}