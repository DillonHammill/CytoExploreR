## CYTO_GATE_CLUST -------------------------------------------------------------

#' Use clustering algorithm to gate populations
#' 
#' @param x object of class \code{GatingSet}.
#' 
#' @export
cyto_gate_clust <- function(x,
                            parent = NULL,
                            alias = NULL,
                            channels = NULL,
                            type = "FlowSOM", 
                            merge_by = "all",
                            gatingTemplate = NULL,
                            group_by = NULL,
                            slot = NULL,
                            input = "flowFrame",
                            inverse = FALSE,
                            ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET REQUIRED
  if(!cyto_class(x, "GatingSet")) {
    stop(
      paste0(
        "'cyto_gate_clust()' only supports GatingHierarchy and GatingSet ",
        "objects!"
      )
    )
  }
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    channels <- channels[!grepl("-H$|-W$", channels, ignore.case = TRUE)]
  } else {
    channels <- cyto_channels_extract(x, channels = channels)
  }
  
  # CURRENT NODES
  nodes <- cyto_nodes(x, path = "auto")
  
  # GATINGTEMPLATE CHECKS
  if(cyto_class(x, "GatingSet")) {
    # ACTIVE GATINGTEMPLATE
    if (is.null(gatingTemplate)) {
      gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
    }
    # CHECK EXISTING ENTRIES IN GATINGTEMPLATE
    gt <- .cyto_gatingTemplate_check(
      parent, 
      alias, 
      gatingTemplate
    )
    # CREATE GATINGTEMPLATE
    if (is.null(gt)) {
      message(
        paste("Creating", gatingTemplate, "to save the constructed gate(s).")
      )
      cyto_gatingTemplate_create(gatingTemplate, active = TRUE)
      gt <- cyto_gatingTemplate_read(gatingTemplate, data.table = TRUE)
    }
  }
  
  # ALGORITHM REFERENCES -------------------------------------------------------
  
  # GATING MESSAGE
  message(
    paste0(
      "Applying ", cyto_func_name(type), "() ",
      "clustering algorithm to the ", parent, " population..."
    )
  )
  
  # DEFAULT ALGORITHM TYPES
  if(is.character(type)) {
    # NIMBUS
    if(grepl("^nimbus", type, ignore.case = TRUE)) {
      # LOAD NIMBUS
      cyto_require(
        "nimbus",
        source = "GitHub",
        repo = "DillonHammill/nimbus"
      )
    # FLOWSOM
    } else if(grepl("^FlowSOM$", type, ignore.case = TRUE)) {
      # LOAD FLOWSOM
      cyto_require(
        "FlowSOM",
        source = "BioC",
        repo = "SofieVG/FlowSOM",
        ref = paste0(
          "Van Gassen S, et al. (2015). FlowSOM: Using self-organising maps ",
          "for visualisation and interpretation of cytometry data. Cytometry ",
          "A 87(7) 10.1002/cyto.a.22625"
        )
      )
    # IMMUNOCLUST
    } else if(grepl("^ImmunoClust$", type, ignore.case = TRUE)) {
      # LOAD IMMUNOCLUST
      cyto_require(
        "immunoClust",
        source = "BioC",
        ref = paste0(
          "Sorensen T, et al. (2015). immunoClust--An automated analysis ",
          "pipeline for the identification of immunophenotypic signatures ",
          "in high-dimensional cytometric datasets. Cytometry A 87(7) ",
          "10.1002/cyto.a.22626"
        )
      )
    # KMEANS
    } else if(grepl("^kmeans$", type, ignore.case = TRUE)) {
      # USES STATS::KMEANS()
    }
  }
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # GROUP_BY -> MERGE_BY
  if(length(group_by) > 0) {
    message(
      "Support for 'group_by' is ending. Please use 'merge_by' instead."
    )
    merge_by <- group_by
  }
  
  # PREPARE GROUP_BY
  if(all(is.character(merge_by))) {
    if(all(merge_by == "all")) {
      group_by <- NA
    } else {
      group_by <- paste(merge_by, collapse = ":")
    }
  } else {
    group_by <- NA
  }

  # ADD POPULATION(S) TO GATINGSET
  pop <- suppressWarnings(
    suppressMessages(
      gs_add_gating_method(
        gs = x,
        alias = "*",
        parent = parent,
        pop = "*",
        dims = "", # EMPTY DIMS
        gating_method = "cyto_gate_clust",
        gating_args = list(
          "type" = type,
          "input" = input,
          "params" = channels,
          "slot" = slot,
          "inverse" = inverse,
          "alias" = alias,
          "openCyto.minEvents" = -1,
          ...
        ),
        groupBy = group_by,
        collapseDataForGating = TRUE, # REQUIRED
        preprocessing_method = "pp_cyto_gate_clust",
        preprocessing_args = list(
          "parent" = parent
        )
      )
    )
  )
  
  # NEW NODES
  new_nodes <- cyto_nodes(x, path = "auto")
  new_nodes <- new_nodes[!new_nodes %in% nodes]
  
  # EXTRACT GATING ARGUMENTS
  gating_args <- eval(
    parse(
      text = paste0("list(", pop$gating_args, ")")
    )
  )
  
  # UPDATE ALIAS IN GATING_ARGS
  gating_args$alias <- paste0(new_nodes, collapse = ",")
  pop$gating_args <- CytoExploreR_.argDeparser(gating_args)
  
  # ADD POPULATIONS TO GATINGTEMPLATE
  gt <- rbind(gt, pop)
  
  # WRITING NEW GATINGTEMPLATE ENTRIES
  message(paste("Re-writing", gatingTemplate, "with new gating entries..."))
  
  # SAVE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(gt, gatingTemplate)
  
  # RETURN GATINGSET
  return(x)
  
}

#' Preprocessing method to pass sample counts to clustering plugin
#' @importFrom flowWorkspace gs_pop_get_children
#' @noRd
.pp_cyto_gate_clust <- function(fs,
                                gs,
                                gm,
                                channels,
                                groupBy = NA,
                                isCollapse = NA,
                                parent  = "root",
                                ...) {
  
  # NOTE: OUTPUT NAMES = GROUPS -> COLLAPSEDATAFORGATING = TRUE
  
  # SAME GROUP
  if(groupBy %in% "") {
    groupBy <- NA 
  }
  
  # GROUPING VARIABLES SEPARATED BY :
  if(is.character(groupBy)) {
    groupBy <- strsplit(groupBy, ":")[[1]]
  }
  
  # GROUP NAME
  grp <- cyto_groups(
    fs,
    group_by = groupBy,
    details = FALSE,
    sep = ":"
  )
  
  # TODO: OPENCYTO COMPATIBILITY - COMBINED EVENTS GROUP NAME
  if(grp %in% "Combined Events") {
    grp <- "all"
  }
  
  # PP_RES - TRANSFORMERS | COUNTS | INDEX
  return(
    list(
      "trans" = cyto_transformers_extract(gs),
      "counts" =   cyto_apply(
        fs,
        FUN = "cyto_stat_count",
        input = "matrix",
        channels = cyto_channels(fs)[1],
        copy = FALSE,
        simplify = FALSE
      ),
      "nodes" = gs_pop_get_children(
        gs, 
        parent,
        path = "auto",
        showHidden = TRUE
      )
    )
  )
  
}

#' openCyto plugin
#' @noRd
.cyto_gate_clust <- function(fr,
                             pp_res,
                             channels,
                             type = "FlowSOM",
                             input = "flowFrame",
                             params = NULL,
                             slot = NULL,
                             inverse = FALSE,
                             alias = "cluster",
                             ...){
  
  # TODO: ADD FRAMEID TO FILTERRESULTS?
  # TODO: APPLY GATE TO NEW DATA - RESET ALIAS?
  
  # BYPASS EMPTY CYTOFRAME
  if(nrow(fr) == 0) {
    return(
      structure(
        lapply(
          names(pp_res$counts),
          function(z) {
            res <- as(factor(), "filterResult")
            res@filterId <- "cyto_gate_clust"
            return(res)
          }
        ),
        names = names(pp_res$counts)
      )
    )
  }
  
  # CHANNELS FOR GATING PASSED THROUGH PARAMS ARGUMENT
  if(is.null(params)) {
    params <- cyto_channels(
      fr,
      exclude = c("Event", "Sample", "Time")
    )
  }
  
  # DEFAULT CLUSTERING METHODS
  if(is.character(type)) {
    # NIMBUS
    if(grepl("^nimbus$", type, ignore.case = TRUE)) {
      
    # FLOWSOM
    } else if(grepl("^FlowSOM$", type, ignore.case = TRUE)) {
      # TODO: CHECK DIMENSIONS OF SOM REQUIRED?
      # FLOWSOM DEFAULT ARGUMENTS
      args_default <- list(
        compensate = FALSE,
        transform = FALSE,
        scale = FALSE,
        silent = FALSE,
        seed = 56
      )
      # COMBINE INCOMING ARGUMENTS
      args <- list(...)
      args <- c(
        args,
        args_default[!names(args_default) %in% names(args)]
      )
      # PREPARE DATA - INVERSE DATA OPTIONAL
      fr <- cyto_data_extract(
        fr,
        format = "flowFrame",
        channels = params,
        trans = pp_res$trans,
        inverse = ifelse(args$transform, TRUE, FALSE),
        copy = FALSE
      )[[1]][[1]]
      # ADD DATA TO ARGUMENTS
      args <- c(
        list(fr),
        args
      )
      # RUN FLOWSOM
      res <- cyto_func_call(
        "FlowSOM::FlowSOM",
        args
      )
      # CLUSTERS
      gate <- cyto_func_call(
        "FlowSOM::GetMetaclusters",
        list(res)
      )
    # IMMUNOCLUST
    } else if(grepl("^ImmunoClust$", type, ignore.case = TRUE)) {
      # PREPARE DATA - INVERSE DATA REQUIRED
      fr <- cyto_data_extract(
        fr,
        format = "flowFrame",
        channels = params,
        trans = pp_res$trans,
        inverse = TRUE,
        copy = FALSE
      )[[1]][[1]]
      # IMMUNOCLUST ARGUMENTS
      args <- list(
        fr,
        ...
      )
      # RUN IMMUNOCLUST
      gate <- cyto_func_call(
        "immunoClust::cell.process",
        args
      )@label
    # KMEANS  
    } else if(grepl("^kmeans$", type, ignore.case = TRUE)) {
      # PREPARE DATA - MATRIX REQUIRED
      fr <- cyto_data_extract(
        fr,
        format = "matrix",
        channels = params,
        trans = pp_res$trans,
        inverse = FALSE,
        copy = FALSE
      )[[1]][[1]]
      # ARGUMENTS
      args <- list(
        fr,
        ...
      )
      # RUN KMEANS
      gate <- cyto_func_call(
        "stats::kmeans",
        args
      )$cluster
    # CONVERT TYPE TO FUNCTION
    } else {
      type <- cyto_func_match(type)
    }
  }
  
  # NON-STANDARD CLUSTERING ALGORITHM
  if(is.function(type)) {
    # CONVERT DATA TO REQUIRED FORMAT - DATA RESTRICTED TO CHANNELS
    if(!cyto_class(fr, input, TRUE)) {
      fr <- cyto_data_extract(
        fr,
        format = input,
        channels = params,
        trans = pp_res$trans,
        inverse = inverse,
        copy = FALSE,
      )[[1]][[1]]
    }
    # APPLY CLUSTERING ALGORITHM
    gate <- cyto_slot(
      cyto_func_call(
        type,
        list(fr, ...)
      ),
      slot = slot
    )
  }
  
  # CONVERT GATE TO FACTOR IF REQUIRED
  if(!is.factor(gate)) {
    gate <- as.factor(gate)
  }
  
  # CONVERT NA -> 0
  if(any(is.na(gate))) {
    `levels<-` (addNA(gate), c(levels(gate), 0))
  }
  
  # ALIAS
  if(length(alias) == 0) {
    alias <- cyto_func_name(type)
  } else if(length(alias) == 1) {
    alias <- unlist(strsplit(alias, ","))
  } 
  
  # ALIAS LENGTH CHANGED - APPLY GATINGTEMPLATE TO NEW DATA
  if(length(alias) > 1 & length(alias) != length(levels(gate))) {
    warning(
      paste0(
        "'alias' is of the incorrect length - naming the clustered ",
        "populations with the name of the clustering function instead."
      )
    )
    alias <- cyto_func_name(type)
  }
  
  # PREPARE ALIAS
  if(length(alias) != length(levels(gate))) {
    # CHECK EXISTING INDICES
    index <- 0
    # PULL OUT INDICES FROM EXISTING NODES
    if(length(pp_res$nodes) > 0) {
      nodes <- gsub(
        paste0("^", alias, "-([0-9]+)$"),
        "\\1",
        pp_res$nodes
      )
      nodes <- nodes[!nodes %in% pp_res$nodes]
      # INDEX
      if(length(nodes) > 0) {
        index <- max(as.numeric(nodes))
      }
    }
    # ALIAS - INDEX
    alias <- paste0(alias, "-", seq(index + 1, index + length(levels(gate))))
  }
  
  # LABEL CLUSTERS IN SIZE ORDER -> LARGEST TO SMALLEST
  levels(gate)[order(table(gate))] <- alias
  
  # SPLIT GATE PER SAMPLE
  gate <- split(
    gate,
    rep(names(pp_res$counts), unlist(pp_res$counts))
  )
  
  # EMPTY SAMPLES
  nms <- names(pp_res$counts)[!names(pp_res$counts) %in% names(gate)]
  if(length(nms) > 0) {
    gate <- c(
      gate,
      structure(
        lapply(
          nms,
          function(z) {
            g <- factor()
            levels(g) <- levels(gate[[1]])
            return(g)
          }
        )
      )
    )
  }
  
  # ORDER GATES
  gate <- gate[names(pp_res$counts)]
  
  # CONVERT TO FILTERRESULT - RETAINS SET POPULATION NAMES
  gate <- structure(
    lapply(
      gate,
      function(z) {
        res <- as(z, "filterResult")
        res@filterId <- "cyto_gate_clust"
        return(res)
      }
    ),
    names = names(gate)
  )
  
  # RETURN CLUSTERING FILTER
  return(gate)
  
}
