## CYTO_GATE_AUTO --------------------------------------------------------------

#' CytoExploreR wrapper around openCyto automated gating methods
#'
#' Perform automated using the algorithms supported by \code{openCyto} without
#' the need to write any plugins. The algorithms used within
#' \code{cyto_gate_auto()} behave in much the same way as they do within the
#' \code{openCyto} package except that any range related arguments (e.g. min,
#' max and target) expect values on the linear scale. \code{cyto_gate_auto()}
#' will automatically map these values to the transformed scale prior to passing
#' them to the specified \code{openCyto} gating function. Eventually, automated
#' gates will be edited through the \code{cyto_gate_edit()} but in the meantime,
#' users should use \code{cyto_gatingTemplate_edit()} instead.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to gate using the specified
#'   automated gating algorithm supplied to \code{type}.
#' @param alias names(s) to use for the newly gated populations.
#' @param channels names of the channels or markers to be be passed to the
#'   automated gating algorithm.
#' @param type the name of the automated gating algorithm to use, options
#'   include \code{"mindensity"}, \code{"quantileGate"}, \code{"tailgate"},
#'   \code{"boundary"}, \code{"singletGate"}, \code{"flowClust"},
#'   \code{"quadGate.tmix"} or \code{"flowDensity"}, set to \code{"flowClust"}
#'   by default.
#' @param merge_by a vector of experiment variables to merge the data into
#'   groups prior to applying the automated gating algorithm supplied to
#'   \code{type}, set to \code{"all"} by default to merge all samples prior to
#'   gating. If \code{merge_by} is used, the specified automated gating
#'   algorithm will be applied separately to each group of merged samples.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'   \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'   saved, set to \code{cyto_gatingTemplate_active()} by default.
#' @param group_by included for compatibility with \code{openCyto} but is
#'   treated in the same way as \code{merge_by}.
#' @param slot a character string or function to apply to the output of a custom
#'   automated gating algorithm to extract the gate object after gating. The
#'   purpose of \code{slot} is to ensure that a valid gate type (i.e. a
#'   rectangleGate) is passed to \code{openCyto} during the gating step.
#' @param events indicates the number of events to pass to the automated gating
#'   algorithm.
#' @param input indicates the format in which the data should be passed to the
#'   automated gating algorithm passed to \code{type}, options include
#'   \code{"flowFrame"}, \code{"cytoframe"}, \code{"flowSet"}, \code{"cytoset"}
#'   or \code{"matrix"}.
#' @param inverse logical to indicate whether inverse data transformations
#'   should be applied to the data prior to passing it to the automated gating
#'   algorithm specified by \code{type}, set to FALSE by default.
#' @param scale optional argument to scale each channel prior to passing data to
#'   the specified automated gating algorithm, options include \code{"range"},
#'   \code{"mean"}, \code{"median"} or \code{"zscore"}. Set to \code{FALSE} by
#'   default to bypass the data scaling step.
#' @param seed numeric to set seed when \code{events != 1}.
#' @param plot logical indicating whether the gating result should be displayed
#'   using \code{cyto_plot()} set to TRUE by default.
#' @param ... additional arguments passed to the automated gating function or
#'   \code{openCyto::gs_add_gating_method()}.
#'
#' @importFrom openCyto gs_add_gating_method
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @return a GatingHierarchy or GatingSet with automated gates added and updated
#'   in the gatingTemplate.
#'
#' @seealso \code{\link{cyto_gate_draw}}
#' @seealso \code{\link{cyto_gate_clust}}
#' @seealso \code{\link{cyto_gatingTemplate_edit}}
#' @seealso \code{\link{cyto_gate_remove}}
#' @seealso \code{\link{cyto_gate_rename}}
#'
#' @export
cyto_gate_auto <- function(x,
                           parent = NULL,
                           alias = NULL,
                           channels = NULL,
                           type = "flowClust",
                           merge_by = "all",
                           gatingTemplate = NULL,
                           group_by = NULL,
                           slot = NULL,
                           events = 1,
                           input = "flowFrame",
                           inverse = FALSE,
                           scale = FALSE,
                           seed = 2022,
                           plot = TRUE,
                           ...) {
  
  # TODO: ADD SUPPORT FOR MULTIPLE GATES
  # TODO: ADD STANDARD MANUAL GATING TYPES -= ELLIPSE, RECTANGLE ETC
  
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
    stop(
      "Supply the 'channels' to gate!"
    )
  }
  channels <- cyto_channels_extract(
    x,
    channels = channels
  )
  
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

  # APPLY AUTOMATED GATING ALGORITHM -------------------------------------------
  
  # GATING MESSAGE
  message(
    paste0(
      "Applying ", cyto_func_name(type), "() ",
      "automated gating algorithm to the ", parent, " population..."
    )
  )
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # GROUP_BY -> MERGE_BY
  if(length(group_by) > 0) {
    message(
      "Support for 'group_by' is ending. PLease use 'merge_by' instead."
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
  
  # DEFAULT GATING METHOD ARGUMENTS
  gm <- list(
    gs = x,
    alias = paste0(
      alias,
      collapse = ","
    ),
    parent = paste0(
      parent,
      collapse = ","
    ),
    dims = paste0(
      channels,
      collapse = ","
    ),
    pop = "+",
    gating_method = "cyto_gate_auto",
    gating_args = list(
      "parent" = parent,
      "alias" = alias,
      "type" = type,
      "input" = input,
      "inverse" = inverse,
      "events" = events,
      "scale" = scale,
      "seed" = seed,
      "slot" = slot,
      "plot" = plot
    ),
    collapseDataForGating = TRUE,
    groupBy = group_by,
    preprocessing_method = "pp_cyto_gate_auto",
    preprocessing_args = NA,
    strip_extra_quotes = FALSE
  )
  
  # PREPARE ARGUMENTS
  args <- list(...)
  
  # REMOVE ILLEGAL ARGUMENTS
  args <- args[
    !names(args) %in% c(
      "gs",
      "dims",
      "groupBy",
      "gating_method", 
      "preprocessing_method",
      "preprocessing_args"
    )
  ]
  
  # PARALLEL GATING ARGUMENTS
  if(any(c("mc.cores", "parallel_type", "cl") %in% names(args))) {
    # APPEN GATING METHOD ARGUMENTS
    gm <- c(
      gm,
      args[
        names(args) %in% c(
          "mc.cores",
          "parallel_type",
          "cl"
        )
      ]
    )
    # EXCLUDE FROM GATING ARGUMENTS
    args <- args[
      !names(args) %in% c(
        "mc.cores", 
        "parallel_type",
        "cl"
      )
    ]
  }
  
  # GATING METHOD ARGUMENTS
  if(any(names(args) %in% names(gm))) {
    nms <- names(gm)[names(gm) %in% names(args)]
    gm[nms] <- args[nms]
    args <- args[!names(args) %in% nms]
  }
  
  # DEPARSE GATING ARGUMENTS
  if(length(args) > 0) {
    gm[["gating_args"]] <- c(
      gm[["gating_args"]],
      args
    )
  }
  
  print(gm)
  
  # ADD POPULATIONS TO GATINGSET
  pop <- cyto_func_call(
    "gs_add_gating_method",
    gm
  )
  
  # ADD GATES TO GATINGTEMPLATE
  gt <- rbind(gt, pop)
  
  # WRITING NEW GATINGTEMPLATE ENTRIES
  message(
    paste(
      "Re-writing",
      gatingTemplate,
      "with new gating entries..."
    )
  )
  
  # SAVE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(gt, gatingTemplate)
  
  # RETURN GATINGSET
  return(x)
  
}

#' Preprocessing method to pass transformers to gating function
#' @noRd
.pp_cyto_gate_auto <- function(fs, 
                               gs,
                               gm,
                               channels,
                               groupBy = NA,
                               isCollapse = NA,
                               ...) {
 
  # GROUPBY
  if(is.character(groupBy) & nchar(groupBy) == 0) {
    groupBy <- NA
  }
  
  # GROUPING VARIABLES SEPARATED BY COLON
  if(is.character(groupBy)) {
    groupBy <- unlist(strsplit(groupBy, ":"))
  }
  
  # GROUP NAME
  grp <- cyto_groups(
    fs,
    group_by = groupBy,
    details = FALSE
  )
  
  # RETURN - GROUP & TRANSFORMERS
  return(
    list(
      "group" = grp,
      "trans" = cyto_transformers_extract(gs)
    )
  )
   
}

#' openCyto plugin
#' @noRd
.cyto_gate_auto <- function(fr,
                            pp_res,
                            channels,
                            parent = "root",
                            alias = NULL,
                            type = "flowClust",
                            input = "flowFrame",
                            inverse = FALSE,
                            events = 1,
                            scale = FALSE,
                            seed = 2022,
                            slot = NULL,
                            plot = TRUE,
                            ...) {
  
  # TODO: INVERSE CAUSE ISSUES WITH GATE CO-ORDINATES!
  
  # PRECAUTIONARY - DON'T MODIFY DATA IN PLACE
  if(cyto_class(fr, "cytoframe")) {
    fr <- cyto_copy(fr)
  }
  
  # TRANSFORMERS
  trans <- pp_res$trans
  
  # GROUP NAME
  grp <- pp_res$group
  
  # INVERSE TRANSFORM
  if(inverse & !.all_na(trans)) {
    fr <- cyto_transform(
      fr,
      trans = trans,
      inverse = inverse,
      plot = FALSE,
      quiet = TRUE
    )
  }
  
  # SCALE DATA PRIOR TO GATING
  if(!scale %in% FALSE) {
    cyto_exprs(fr) <- cyto_stat_scale(
      cyto_exprs(
        fr,
        channels = channels,
        drop = FALSE
      ),
      type = scale
    )
  }
  
  # EXTRACT DATA
  fr <- cyto_data_extract(
    fr,
    format = input,
    channels = channels,
    trans = trans,
    inverse = inverse,
    events = events,
    seed = seed,
    copy = FALSE
  )[[1]][[1]]
  
  # PULL DOWN ARGUMENTS
  args <- list(fr, ...)
  
  # DEFAULT AUTOMATED GATING METHODS
  if(is.character(type)) {
    # FLOWCLUST
    if(grepl("^flowClust", type, ignore.case = TRUE)) {
      type <-  if(length(channels) == 1) {
        "openCyto:::.gate_flowclust_1d"
      } else {
        "openCyto:::.gate_flowclust_2d"
      }
    # MINDENSITY
    } else if(grepl("^mindensity", type, ignore.case = TRUE)) {
      type <- "openCyto:::.gate_mindensity"
    # TAILGTE
    } else if(grepl("^tailgate", type, ignore.case = TRUE)) {
      type <- "openCyto:::.gate_tail"
    # SINGLETGATE
    } else if(grepl("^singletGate", type, ignore.case = TRUE)) {
      type <- "openCyto:::.singletGate"
    # QUANTILEGATE
    } else if(grepl("^quantileGate", type, ignore.case = TRUE)) {
      type <- "openCyto:::.singletGate"
    # RANGEGATE
    } else if(grepl("^rangeGate", type, ignore.case = TRUE)) {
      type <- "openCyto:::.rangeGate"
    # BOUNDARY
    } else if(grepl("^boundary", type, ignore.case = TRUE)) {
      type <- "openCyto:::.boundary"
    # SEQUENTIAL QUADRANT GATE
    } else if(grepl("^q.*seq", type, ignore.case = TRUE)) {
      type <- "openCyto:::.gate_quad_sequential"
    # TMIX QUADRANT GATE
    } else if(grepl("^q.*tmix", type, ignore.case = TRUE)) {
      type <- "openCyto:::.gate_quad_tmix"
    # FLOWDENSITY
    } else if(grepl("^flowdensity", type, ignore.case = TRUE)) {
      cyto_require(
        "flowDensity",
        source = "BioC"
      )
      message(
        paste0(
          "Malek M, Taghiyar MJ, Chong L, Finak G, Gottardo R, Brinkman RR. ",
          "flowDensity: reproducing manual gating of flow cytometry data by ",
          "automated density-based cell population identification. ",
          "Bioinformatics. 2015 Feb 15;31(4):606-7. ",
          "doi: 10.1093/bioinformatics/btu677."
        )
      )
      if(length(channels)==2) {
        type <- "flowDensity:::.flowDensity.2d"
      } else {
        type <- "flowDensity:::.flowDensity.1d"
      }
    }
    # OPENCYTO EXPECTS ARGUMENTS ON TRANSFORMED SCALE
    if(!.all_na(trans) & grepl("^openCyto", type, ignore.case = TRUE)) {
      # ARGUMENTS TO TRANSFORM LINEAR -> TRANSFORMED SCALE
      ind <- grep(
        "^min$|^max$|^target$|_min$|_max$",
        names(args)
      )
      # ARGUMENTS SHOULD BE OF LENGTH CHANNELS
      lapply(
        ind,
        function(z) {
          for(i in seq_along(channels)) {
            if(channels[i] %in% names(trans)) {
              args[[z]] <<- .cyto_transform(
                args[[z]][i],
                channel = channels[i],
                trans = trans,
                inverse = FALSE
              )
            }
          }
        }
      )
    }
    # CONVERT TYPE TO FUNCION
    type <- cyto_func_match(type)
  }
  
  # CUSTOM AUTOMATED GATING ALGORITHM
  if(is.function(type)) {
    # APPLY AUTOMATED GATING ALGORITHM
    gate <- cyto_slot(
      cyto_func_call(
        type,
        args = if("channels" %in% cyto_func_args(type)) {
          args <- c(args, list("channels" = channels))
        } else {
          args
        }
      )
    )
    # PLOT
    if(plot) {
      # TODO: CHECK FOR SUPPORTED GATE TYPES
      cyto_plot(
        fr,
        channels = channels,
        gate = gate,
        trans = trans,
        title = paste0(
          grp, "\n", parent
        )
      )
    }
  # UNSUPPORTED AUTOMATED GATING ALGORITHM
  } else {
    stop(
      "Invalid 'type' argument refere to docs for supported algorithms"
    )
  }
  
  # TODO: HANDLING OF MULTIGATE RETURNS?
  
  # RETURN GATE(S)
  return(gate)
  
}
