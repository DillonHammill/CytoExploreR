## CYTO_PLOT INTERNAL FUNCTIONS ------------------------------------------------

# These functions are used within cyto_plot() to prepare arguments for use in
# cyto_plot().

## BUILD -----------------------------------------------------------------------

#' Build cyto_plot from arguments
#' 
#' @param args list of arguments prepared by cyto_plot.
#' 
#' @noRd
.cyto_plot_build <- function(args){
  
  # ARGUMENTS - CYTO_PLOT CLASS
  class(args) <- "cyto_plot"
  
  # CYTO_PLOT_EMPTY
  cyto_plot_empty(args)
  
  # HISTOGRAMS
  if(length(args$channels) == 1) {
    cyto_plot_hist(args)
  # POINTS & CONTOURS
  } else {
    cyto_plot_point(args)
  }
  
  # AUTOMATIC LABEL POSITION - BYPASS ON SAVING - MUST BE AFTER CYTO_PLOT_EMPTY
  if(!getOption("cyto_plot_save") &
     args$label == TRUE &
     args$label_position == "auto") {
    label_text_xy <- do.call(".cyto_plot_label_coords", args)
    args$label_text_x <- label_text_xy[, "x"]
    args$label_text_y <- label_text_xy[, "y"]
  }
  
  # GATES & LABELS
  if (!.all_na(args$gate)) {
    # PLOT GATE & ASSOCIATED LABELS
    label_text_xy <- .cyto_plot_gate_label(args)
    args$label_text_x <- label_text_xy[, "x"]
    args$label_text_y <- label_text_xy[, "y"]
    # LABELS
  } else {
    if (args$label == TRUE) {
      label_text_xy <- do.call("cyto_plot_labeller", args)
      args$label_text_x <- label_text_xy[, "x"]
      args$label_text_y <- label_text_xy[, "y"]
    }
  }
  
  return(args)
  
}

## DATA ------------------------------------------------------------------------

#' Prepare data for cyto_plot
#' @return list of singular cytosets for each plot
#' @importFrom flowWorkspace cytoset cf_get_uri
#' @importFrom purrr transpose
#' @noRd
.cyto_plot_data <- function(x,
                            parent = "root",
                            overlay = NA,
                            merge_by = "name",
                            select = NULL,
                            channels = NULL,
                            events = 50000,
                            hist_layers = NA,
                            hist_stack = 0,
                            barcode = FALSE,
                            seed = 42,
                            ...) {

  # TODO: COMPUTE NROW ON CYTOFRAME DIRECTLY NOT EXPRS
  
  # DATA PREPARED ALREADY ------------------------------------------------------
  
  # BYPASS DATA PREPARTION
  if(cyto_class(x, "list", TRUE)) {
    # LIST OF CYTOSET LISTS
    if(all(LAPPLY(x, "cyto_class", "flowSet"))) {
      return(
        list(x)
      )
    }
    return(x)
  # .CYTO_PLOT_DATA() CALLED
  } else if(!is.null(cyto_option("cyto_plot_method")) & 
            cyto_option("cyto_plot_data")) {
    if(.all_na(overlay)) {
      return(
        list(
          structure(
            list(x),
            names = cyto_names(x)
          )
        )
      )
    } else {
      return(
        list(
          c(
            structure(
              list(x),
              names = cyto_names(x)
            ),
            overlay
          )
        )
      )
    }
  }
  
  # CHECKS ---------------------------------------------------------------------
  
  # INITIALISE GS/GH
  gs <- NULL
  gh <- NULL
  
  # GATINGSET
  if(cyto_class(x, "GatingSet", TRUE)) {
    gs <- x
    gh <- x[[1]]
  # GATINGHIERARCHY
  } else if(cyto_class(x, "GatingHierarchy", TRUE)) {
    gh <- x
  # CYTOSET
  } else if(cyto_class(x, "cytoset")) {
    # DO NOTHING
  # COERCE MATRIX/DATA.FRAME/CYTOFRAME TO CYTOSET
  } else {
    x <- tryCatch(
      as(x, "cytoset"),
      error = function(e) {
        stop(
          paste0(
            "cyto_plot() only supports cytoset, GatingHierarchy or GatingSet ",
            "objects!"
          )
        )
      }
    )
  }
  
  # EXTRACT DATA
  x <- cyto_data_extract(
    x,
    parent = parent,
    select = select,
    format = "cytoset",
    copy = FALSE
  )
  parent <- names(x) # path
  
  # GROUP
  x <- cyto_group_by(x[[1]], merge_by)
  
  # PREPARE X
  x <- structure(
    lapply(x, function(z){
      if(is.null(gh)) {
        return(
          list(z)
        )
      } else {
        return(
          structure(
            list(z), 
            names = parent
          )
        )
      }
    }),
    names = names(x)
  )
  
  # OVERLAY - POPULATION NAMES
  if(!.all_na(overlay)) {
    # OVERLAY - POPULATION NAMES
    if(cyto_class(overlay, "character")) {
      # GATINGHIERRACHY REQUIRED
      if(!is.null(gh)) {
        # DESCENDANTS
        if(any(grepl("^descendants", overlay, ignore.case = TRUE))) {
          overlay <- overlay[!grepl("^descendants", 
                                    overlay, 
                                    ignore.case = TRUE)]
          overlay <- c(overlay, tryCatch(
            gh_pop_get_descendants(
              gh,
              parent,
              path = "auto"
            ),
            error = function(e){
              return(NA)
            }
          ))
        } 
        # CHILDREN
        if(any(grepl("^children", overlay, ignore.case = TRUE))) {
          overlay <- overlay[!grepl("^children",
                                    overlay,
                                    ignore.case = TRUE)]
          overlay <- c(overlay, tryCatch(
            gh_pop_get_children(
              gh,
              parent,
              path = "auto"
            ),
            error = function(e){
              return(NA)
            }
          ))
        }
        # REMOVE DUPLICATES
        overlay <- unique(overlay)
        # EXTRACT DATA
        overlay <- structure(
          lapply(overlay, function(z) {
            if(.all_na(z)) {
              return(NA)
            } else {
              cyto_data_extract(
                if(is.null(gs)){
                  gh
                }else{
                  gs
                }, # GatingHierarchy/GatingSet
                parent = cyto_nodes_convert(
                  gh,
                  nodes = z,
                  anchor = parent
                ),
                copy = FALSE,
                format = "cytoset"
              )[[1]]
            }
          }),
          names = overlay
        )
        # INVALID OVERLAY
      } else {
        overlay <- NA
      }
    }
  }
  # OVERLAY - CYTOSET/LIST OF CYTOSETS
  if(!.all_na(overlay)) {
    # OVERLAY -> LIST OF CYTOSETS
    if(!cyto_class(overlay, "list", TRUE)) {
      overlay <- structure(
        list(overlay),
        names = ifelse(
          length(overlay) == 1, 
          cyto_names(overlay), 
          NULL
        )
      )
    }
    # OVERLAY - SELECT & GROUP
    overlay  <- structure(
      lapply(seq_along(overlay), function(z){
        cs <- overlay[[z]]
        # CHECK
        if(!cyto_class(cs, "cytoset")) {
          cs <- tryCatch(
            as(cs, "cytoset"),
            error = function(e) {
              stop("cyto_plot only supports cytosets for overlays!")
            }
          )
        }
        # SELECT - IF POSSIBLE
        cs <- tryCatch(
          cyto_select(cs, select),
          error = function(e){return(cs)}
        )
        # GROUP - ONLY IF SAME GROUPS
        if(setequal(cyto_groups(cs, merge_by), names(x))) {
          cs <- cyto_group_by(cs, merge_by)
        } else {
          cs <- structure(
            rep(list(cs), length(x)),
            names  = names(x)
          )
        }
        # LAYER
        return(cs)
      }),
      names = names(overlay)
    )
    # TRANSPOSE OVERLAY - GROUPS
    overlay <- transpose(overlay)
    # OVERLAY EACH PLOT - OVERLAY MUST BE SAME LENGTH AS X FOR TRANSPOSE
    if(length(overlay) != length(x)) {
      overlay <- rep(overlay, length.out = length(x))
    }
    # COMBINE X & OVERLAY
    x <- structure(
      lapply(
        seq_along(x), 
        function(z) {
          c(x[[z]], overlay[[z]])
        }
      ),
      names = names(x)
    )
  }
  
  # SAMPLE & MERGE
  x <- structure(
    # LOOP THROUGH EACH PLOT
    lapply(seq_along(x), function(z){
      # PLOT LAYERS - LIST OF CYTOSETS
      cs_list <- x[[z]]
      # GROUP NAME
      grp <- names(x)[z]
      # CHECK EACH LAYER FOR EVENTS IN PREVIOUS LAYERS - CANNOT STORE GROUP NAME
      i <- lapply(seq_along(cs_list) , function(w){
        # CYTOSET
        cs <- cs_list[[w]]
        # IDENTIFIERS
        ids <- cyto_names(cs)
        # CYTOFRAME
        structure(
          lapply(seq_along(cs), function(q) {
            # IDENTIFIER
            id <- ids[q]
            l <- list(
              "total" = nrow(cyto_exprs(cs[[q]], drop = FALSE)),
              "layer" = NULL,
              "match" = 0,
              "sample" = NULL
            )
            for(v in seq_len(w-1)) {
              # FIRST LAYER
              if(v == 0){
                break()
                # SUBSEQUENT LAYERS
              } else {
                # SAMPLE ID MATCH PREVIOUS LAYER
                if(id %in% cyto_names(cs_list[[v]])) {
                  # MATCHING EVENT IDS?
                  m <- sum(
                    cyto_exprs(cs[[q]], "Event-ID") %in%
                      cyto_exprs(cs_list[[v]][[match(id, ids)]], "Event-ID")
                  ) 
                  # ANY MATCHING EVENTS
                  if(m > 0) {
                    # UPDATE LIST
                    l$layer <- v
                    l$match <- m
                    break()
                  }
                }
              }
            }
            return(l)
          }),
          names = ids
        )
      })
      
      # SAMPLE - SAVE TO CS_SUB_LIST - NEED CS_LIST CANNOT UPDATE IN PLACE
      cs_sub_list <- list()
      lapply(seq_along(cs_list), function(w){
        # CYTOSET
        cs <- cs_list[[w]]
        # IDENTIFIERS
        ids <- cyto_names(cs)
        # ALL CYTOFRAMES CONTAIN NEW EVENTS
        if(is.null(LAPPLY(i[[w]], `[[`, "layer"))) {
          # COMPUTE SAMPLE SIZES
          n <- cyto_sample_n(
            cs,
            events = events
          )
          for(id in ids) {
            i[[w]][[id]]$sample <<- n[id]
          }
          cs_sub_list[[w]] <<- cyto_sample(
            cs,
            events = n,
            seed = seed
          )
        # SOME/ALL CYTOFRAMES CONTAIN EVENTS ON PLOT
        } else {
          # FIND WHICH CYTOFRAMES HAVE EVENTS ON PLOT
          m <- LAPPLY(i[[w]], function(s){
            !is.null(s$layer)
          })
          # STORE CYTOFRAMES IN LIST
          cf_list <- list()
          # CYTOFRAMES WITH EVENTS ON PLOT
          for(id in names(m[m])) {
            # MATCHING LAYER
            l <- i[[w]][[id]]$layer
            # ORIGINAL SIZE OF REFERNCE LAYER
            t <- i[[l]][[id]]$total
            # SAMPLE SIZE OF REFERENCE LAYER
            s <- i[[l]][[id]]$sample
            # CYTOFRAME SUBSET OF PREVIOUS LAYER
            if(i[[w]][[id]]$match == nrow(cyto_exprs(cs[[id]], drop = FALSE))) {
              # USE EVENT IDs FROM REFERENCE LAYER
              cf_list[[id]] <- cs[[id]][
                cyto_exprs(cs[[id]], "Event-ID") %in%
                  cyto_exprs(cs_sub_list[[l]][[id]], "Event-ID")
                ,]
              i[[w]][[id]]$sample <<- nrow(
                cyto_exprs(cf_list[[id]], drop = FALSE)
              )
            # CYTOFRAME OVERLAP WITH PREVIOUS LAYER
            } else {
              # EVENT IDs FROM LAYER SAMPLE
              e <- cyto_exprs(cs[[id]], "Event-ID")[
                cyto_exprs(cs[[id]], "Event-ID") %in%
                  cyto_exprs(cs_sub_list[[l]][[id]], "Event-ID")
              ]
              # NON-OVERLAPPING PORTION OF CYTOFRAME
              cf_new <- cs[[id]][
                !cyto_exprs(cs[[id]], "Event-ID") %in%
                  cyto_exprs(cs_list[[l]][[id]], "Event-ID"),
              ]
              n <- nrow(cyto_exprs(cf_new, drop = FALSE))
              # SAMPLE NON-OVERLAPPING PORTION OF CYTOFRAME
              e <- c(e,
                     cyto_exprs(cf_new)[
                       sample(
                         nrow(cf_new),
                         round(
                           (n/i[[l]][[id]]$total)*(i[[l]][[id]]$sample)
                           )
                         )
                       , "Event-ID"])
              # COMPLETE SAMPLE CYTOFRAME
              cf_list[[id]] <- cs[[id]][
                cyto_exprs(cs[[id]], "Event-ID") %in% e
              , ]
              # STORE SAMPLE SIZE
              i[[w]][[id]]$sample <<- length(e)
            }
          }
          # CYTOFRAMES WITHOUT EVENTS ON PLOT
          # NEW LAYER - ALL NEW CYTOFRAMES
          if(length(m[!m]) == length(cs)) {
            # COMPUTE SAMPLE SIZES
            n <- cyto_sample_n(
              cs,
              events = events
            )
            # SAMPLE
            for(id in names(m[!m])) {
              i[[w]][[id]]$sample <<- n[id]
              cf_list[[id]] <- cyto_sample(
                cs[[id]],
                events = n[id],
                seed = seed
              )
            }
          # SOME NEW CYTOFRAME(S)
          } else if(length(m[!m]) != 0 & length(m[!m]) < length(cs)) {
            # RATIO SAMPLING TO COMPUTE N FOR NEW CYTOFRAMES
            n <- round(
              sum(LAPPLY(i[[w]][names(m[m])], `[[`, "sample")) *
              (sum(LAPPLY(i[[w]][names(m[!m])], `[[`, "total")) /
                 sum(LAPPLY(i[[w]][names(m[m])], `[[`, "total")))
            )
            # COMPUTE SAMPLE SIZES FOR NEW CYTOFRAMES
            n <- cyto_sample_n(cs[names(m[!m])],
                               events = n)
            # SAMPLE
            for(id in names(m[!m])) {
              i[[w]][[id]]$sample <<- n[id]
              cf_list[[id]] <- cyto_sample(cs[[id]],
                                           events = n[id])
            }
          }
          # ADD CYTOFRAMES TO NEW CYTOSET
          cs_sub_list[[w]] <<- cytoset(cf_list)
        }
      })
      # COERCE
      structure(
        lapply(seq_along(cs_sub_list), function(r){
          if(length(cs_sub_list[[r]]) > 1) {
            cytoset(
              structure(
                list(
                  as(cs_sub_list[[r]], "cytoframe")
                ),
                names = grp
              )
            )
          } else {
            return(cs_sub_list[[r]])
          }
        }),
        names = names(cs_list)
      )
    }),
    names = names(x)
  )
  
  # FORMAT DATA FOR HISTOGRAMS
  if(length(channels) == 1 & all(LAPPLY(x, "length") == 1)) {
    # HIST_LAYERS
    if(.all_na(hist_layers)) {
      if(hist_stack == 0) {
        hist_layers <- rep(1, length(x))
      } else {
        hist_layers <- length(x)
      }
    }
    # CHECK HIST_LAYERS
    if(sum(hist_layers) != length(x)) {
      # SAME NUMBER OF LAYERS PER PLOT
      if(length(x) %% hist_layers[1] != 0) {
        stop(
          "Each plot must have the same number of layers!"
        )
      }
      # REPEAT HIST_LAYERS
      hist_layers <- rep(
        hist_layers, 
        length.out = length(x)/hist_layers[1]
      )
    }
    # UNPACK X
    x <- structure(
      unlist(x),
      names = names(x)
    )
    # FORMAT X
    L <- split(seq_along(x), 
               rep(1:(length(x)/hist_layers[1]),
               each = hist_layers[1]))
    x <- structure(
      lapply(
        L,
        function(z) {
          x[z]
        }
      ),
      names = NULL
    )
  }
  
  # DATA READY FOR CYTO_PLOT
  return(x)
  
}

## GATES -----------------------------------------------------------------------

#' Prepare gates for cyto_plot
#' @param x cytoset, GatingHierarchy or GatingSet.
#' @param parent name of the parent population
#' @param channels 
#' @param alias names of gated populations
#' @param gate gate objects
#' @param channels channels used to construct plot (required for gate
#'   conversion)
#' @param merge_by how the data will be merged within cyto_plot (need gates per
#'   group)
#' @return (named) list of gate lists
#' @noRd
.cyto_plot_gates <- function(x,
                             parent = "root",
                             channels,
                             alias = NA, 
                             gate = NA,
                             merge_by = "name",
                             select = NULL,
                             negate = FALSE) {
  
  # TODO: BELOW STEP IN CYTO_PLOT()? GROUPS OUT OF SYNC IN CYTO_PLOT_DATA()?
  # PREPARE X - NON-STANDARD DATA STRUCTURES -> CYTOSET
  if(cyto_class(x, "list", TRUE)) {
    x <- structure(
      lapply(
        x,
        function(z) {
          if(!cyto_class(z, c("flowSet", "GatingSet"))){
            z <- as(z, "cytoset")
          }
          return(z)
        }
      )
    )
  } else if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    x <- as(x, "cytoset")
  }
  
  # X - PREPARED LIST of CYTOSETS PER GROUP (DATA PREPARED ALREADY)
  if(cyto_class(x, "list")) {
    grps <- names(x)
    names(grps) <- grps
  # DATA NOT PREPARED YET
  } else {
    # EXPERIMENTAL GROUPS
    grps <- cyto_groups(
      x,
      select = select,
      group_by = merge_by,
      details = TRUE
    )
  }
  
  # ALIAS - GATINGHIERARCHY/GATINGSET
  if(!.all_na(alias)) {
    # CYTOSET
    if(!cyto_class(x, "GatingSet")) {
      message(
        "'alias' is only supported for GatingHierarchy and GatingSet objects"
      )
      alias <- NA
    # GATINGHIERARCHY/GATINGSET
    } else {
      # GATINGHIERARCHY
      gh <- x[[1]]
      # PARENT - AUTO PATH
      parent <- cyto_nodes_convert(
        gh,
        nodes = parent,
        path = "auto"
      )
      # GATINGTEMPLATE - AUTO PATHS
      gt <- gh_generate_template(gh)
      # NO GATES EXIST
      if(nrow(gt) == 0) {
        stop(
          paste0(
            "This ", 
            cyto_class(x),
            "doesn't contain any gated populations!"
          )
        )
      }
      # ALIAS - MUST BE BEFORE PARENT
      gt$alias <- cyto_nodes_convert(
        gh,
        nodes = paste0(gt$parent, "/", gt$alias),
        path = "auto"
      )
      # PARENT
      gt$parent <- cyto_nodes_convert(
        gh,
        nodes = gt$parent,
        path = "auto"
      )
      # EMPTY ALIAS - BOOLEAN FILTERS NOT SUPPORTED (LACK CHANNELS)
      if(any(LAPPLY(alias, ".empty"))) {
        # REMOVE EMPTY ALIAS
        alias <- alias[!LAPPLY(alias, ".empty")]
        # 2D - MATCH BOTH CHANNELS
        if(length(channels) == 2) {
          pops <- gt$alias[
            gt$dims == paste(channels, collapse = ",") |
              gt$dims == paste(rev(channels), collapse = ",")
          ]
          # 1D - SINGLE CHANNEL MATCH
        } else if(length(channels) == 1) {
          pops <- gt$alias[
            LAPPLY(gt$dims, function(z){
              grepl(channels, z)
            })
          ]
        }
        # SEARCH POSSIBLE BOOLEAN GATES
        pops <- c(pops,
                  add = gt$alias[LAPPLY(gt$dims, ".empty")])
        # UPDATE ALIAS
        alias <- c(alias, pops)
      }
      # REMOVE BOOLEAN GATES FROM ALIAS
      bool <- LAPPLY(alias, function(z){
        if(.empty(gt[gt$alias == z, "dims"])) {
          return(z)
        }else {
          return(NULL)
        }
      })
      if(length(bool) > 0){
        alias <- alias[-match(bool, alias)]
      }
      # CHECK BOOLEAN GATES IN ALIAS
      if(length(bool) > 0) {
        bool <- LAPPLY(bool, function(z){
          # CHECK ADDED GATES - MUST ANCHOR TO PARENT (BYPASS)
          if(grepl("add", names(bool)[z])) {
            z <- tryCatch(
              cyto_nodes_convert(
                gh,
                nodes = z,
                anchor = parent
              ),
              error = function(e){
                return(NULL)
              })
            # INVALID BOOLEAN GATE
            if(is.null(z)) {
              return(z)
            }
          }
          # EXTRACT GATE
          gate <- gh_pop_get_gate(
            gh,
            cyto_nodes_convert(
              gh,
              nodes = z,
              anchor = parent
            )
          )
          # BOOLEAN GATE
          if(cyto_class(gate, "booleanFilter")) {
            # BOOLEAN LOGIC
            logic <- gate@deparse
            # ONLY AND/NOT SUPPORTED
            if(!grepl("^!", logic)) {
              if(!grepl("add", names(bool)[z])) {
                message("Only NOT boolean gates are supported!")
              }
              return(NULL)
            } else if(grepl("|", logic, fixed = TRUE)) {
              if(!grepl("add", names(bool)[z])){
                message("Only AND NOT boolean gates are supported!")
              }
              return(NULL)
            } else {
              # STRIP &! - CHECK ALIAS
              bool_alias <- gsub("^!", "", logic)
              bool_alias <- unlist(strsplit(bool_alias, "&!"))
              # BOOLEAN GATE MUST REFERENCE EVERY POPULATION IN ALIAS
              if(!all(bool_alias %in% alias)) {
                return(NULL)
              } else {
                return(z)
              }
              # # BOOL ALIAS MUST BE IN ALIAS
              # if(!all(bool_alias %in% alias)) {
              #   alias <<- unique(c(alias, bool_alias))
              # }
            }
            # UNSUPPORTED GATE
          } else {
            return(NULL)
          }
        })
        # COMBINE ALIAS & BOOLEAN GATES
        alias <- unique(c(alias, bool))
      }
      # GATES PER GROUP - USE FIRST GH - SELECT HANDLED ABOVE
      alias <- cyto_gate_extract(
        x,
        parent = parent,
        alias = alias,
        merge_by = merge_by
      )
    }
  }
  
  # GATE SUPPLIED MANUALLY -> LIST OF GATE OBJECT LISTS
  if(!.all_na(gate)) {
    # GATE OBJECTS SUPPLIED (NOT A LIST)
    if(!cyto_class(gate, "list", TRUE)) {
      # FILTERS -> GATE OBJECT LIST
      if(cyto_class(gate, "filters")) {
        gate <- unlist(gate)
      # GATE -> GATE OBJECT LIST
      } else {
        gate <- list(gate)
      }
      # STORE FILTER NAMES IN LIST NAMES
      ids <- LAPPLY(gate, function(z){
        tryCatch(
          z@filterId,
          error = function(e){
            return(NA)
          })
      })
      names(gate)[!is.na(ids)] <- ids[!is.na(ids)]
      # REPEAT GATES PER GROUP
      gate <- structure(
        lapply(
          names(grps), 
          function(z){
            return(gate)
          }
        ),
        names = names(grps)
      )
    # GATE OBJECTS SUPPLIED IN LIST
    } else {
      # LIST OF GATE OBJECT LISTS
      if(cyto_class(gate[[1]], "list", TRUE)) {
        # LIST OF GATE OBJECTS PER GROUP
        if(length(gate)!= length(grps)) {
          gate <- rep(gate, length.out = length(grps))
        }
        # NEGATE & FILTERID
        gate <- structure(
          lapply(
            gate, 
            function(w){
              # EXTRACT FILTERS
              gate_list <- unlist(gate[[w]])
              ids <- LAPPLY(
                gate_list, 
                function(z){
                  tryCatch(
                    z@filterId,
                    error = function(e){
                      return(NA)
                    }
                  )
                }
              )
              names(gate_list)[!is.na(ids)] <- ids[!is.na(ids)]
              return(gate_list)
            }
          ),
          names = names(grps)
        )
      # LIST OF GATE OBJECTS  
      } else {
        # EXTRACT FILTERS
        gate <- unlist(gate)
        # STORE FILTER NAMES IN LIST NAMES
        ids <- LAPPLY(gate, function(z){
          tryCatch(z@filterId,
                   error = function(e){
                     return(NA)
                   })
        })
        names(gate)[!is.na(ids)] <- ids[!is.na(ids)]
        # REPEAT GATES PER GROUP
        gate <- structure(
          lapply(
            names(grps), 
            function(z){
              return(gate)
            }
          ),
          names = names(grps)
        )
      }
    }
    # COMBINE GATE & ALIAS
    gate <- structure(
      lapply(
        seq_along(gate), 
        function(z){
          # COMBINE ALIAS
          if(!.all_na(alias)) {
            alias <<- NA
            gates <- c(alias[[z]], gate[[z]])
          } else {
            gates <- gate[[z]]
          }
          # NEGATE
          if(negate) {
            if(length(gates) == 1) {
              gates <- c(gates,
                         list("negate" = !gates[[1]]))
            } else {
              gates <- c(gates,
                         list("negate" = !do.call("|", unname(unlist(gates)))))
            }
          }
          return(gates)
        }
      ),
      names = ifelse(.all_na(alias), names(gate), names(alias))
    )
  # GATE - NA
  } else {
    # ALIAS INSTEAD OF GATE
    if(!.all_na(alias)) {
      gate <- structure(
        lapply(
          alias, 
          function(z){
            # NEGATE
            if(negate) {
              if(length(z) == 1) {
                z <- c(z,
                       list("negate" = !z[[1]]))
              } else {
                z <- c(z,
                       list("negate" = !do.call("|", unname(unlist(z)))))
              }
            }
            return(z)
          }
        ),
        names = names(alias)
      )
    }
  }
  
  # PREPARED GATE OBJECTS
  return(gate)
  
}

## AXES ------------------------------------------------------------------------

## AXES LIMITS ----

#' Compute axes limits for cyto_plot
#'
#' The lower limit is always set to zero unless there is data below this limit.
#' The upper limit is determined by the axes_limits argument.
#'
#' @param x cytometry object(s) which require range calculation.
#' @param parent name of parent population to extract for range calculation.
#' @param channels name of a channel.
#' @param gate list of gates to plot, gate co-ordinates used to compute
#'   axes_limits.
#' @param axes_limits either "auto", "trim", "data" or "machine". "auto" use
#'   data limits but always includes zero.
#' @param plot logical indicating whether a check should be performed for
#'   channel length.
#' @param buffer fraction indicating the amount of buffering to be added on top
#'   of the upper and lower limit, set to 0.03 by default.
#' @param anchor logical indicating if the lower limit should be anchored to
#'   zero if the data range is above zero, set to TRUE by default.
#'
#' @return vector containing minimum and maximum values.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_axes_limits <- function(x,
                                   parent = NULL,
                                   channels = NA,
                                   gate = NA,
                                   axes_limits = "auto",
                                   plot = FALSE,
                                   buffer = 0.03,
                                   anchor = TRUE){
  
  # COMPATIBILITY --------------------------------------------------------------
  
  # FLOWCORE USES INSTRUMENT
  if(axes_limits == "instrument"){
    axes_limits <- "machine"
  }
  
  # DATA -----------------------------------------------------------------------
  
  # LIST CYTOFRAMES/CYTOSETS
  if(cyto_class(x, "list")){
    x <- unlist(x) 
    x <- lapply(x, function(z){
      # GATINGHIERARCHY/GATINGSET
      if(cyto_class(z, "GatingSet")){
        z <- cyto_data_extract(z, 
                               parent = parent,
                               copy = FALSE)[[1]]
      }
      return(z)
    })
  # GATINGHIERARCHY/GATINGSET
  }else if(cyto_class(x, "GatingSet")){
    x <- list(cyto_data_extract(x, 
                                parent = parent,
                                copy = FALSE)[[1]])
  # CYTOFRAME/CYTOSET
  } else {
    x <- list(x)
  }
  
  # CHANNELS -------------------------------------------------------------------
  
  # MARKERS TO CHANNELS
  if(!.all_na(channels)){
    channels <- cyto_channels_extract(x[[1]], channels, plot)
  }else{
    channels <- cyto_channels(x[[1]])
  }
  
  # TIME CHANNEL ALWAYS USE DATA LIMITS
  if(grepl("^Time", channels, ignore.case = TRUE)){
    axes_limits <- "data"
  } 
  
  # RANGE ----------------------------------------------------------------------
  
  # DATA RANGE
  data_range <- lapply(x, function(z){
    if(cyto_class(z, "flowFrame")){
      if(nrow(z) == 0){
        type <- "instrument"
      }else{
        type <- "data"
      }
      # QUANTILE TRIM 1%
      if(type == "trim") {
        rng <- cyto_stat_quantile(
          cyto_exprs(z,
                     channels = channels,
                     drop = FALSE),
          probs = c(0.01, 1)
        )
        rownames(rng) <- c("min", "max")
      } else {
        rng <- suppressWarnings(
          range(z[, channels],
                type = type)
        )
      }
    }else if(cyto_class(z, "flowSet")){
      suppressWarnings(
        cyto_apply(z, function(y){
          if(nrow(y) == 0){
            type <- "instrument"
          }else{
            type <- "data"
          }
          if(axes_limits == "trim") {
            rng <- cyto_stat_quantile(
              cyto_exprs(y,
                         channels = channels,
                         drop = FALSE),
              probs = c(0.01, 1)
            )
            rownames(rng) <- c("min", "max")
          } else {
            rng <- suppressWarnings(
              range(y, type = type)
            )
          }
          return(rng)
        }, 
        input = "cytoframe", 
        channels = channels,
        inverse = FALSE)
      )
    }
  })
  data_range <- do.call("rbind", unname(data_range))
  
  # MIN/MAX DATA RANGE
  data_range <- lapply(seq_len(ncol(data_range)), function(z){
    c(min(data_range[, z]), max(data_range[, z]))
  })
  data_range <- do.call("cbind", data_range)
  colnames(data_range) <- channels
  rownames(data_range) <- c("min", "max")
  
  # MACHINE RANGE
  if(axes_limits == "machine"){
    machine_range <- lapply(x, function(z){
      if(cyto_class(z, "flowFrame")){
        rng <- suppressWarnings(
          range(z,
                type = "instrument")[, channels, drop = FALSE]
        )
      }else if(cyto_class(z, "flowSet")){
        rng <- suppressWarnings(
          cyto_apply(z,
                     "range",
                     type = "instrument",
                     input = "cytoframe",
                     channels = channels,
                     inverse = FALSE)
        )
      }
      return(rng)
    })
    machine_range <- do.call("rbind", machine_range)
    # MIN/MAX DATA RANGE
    machine_range <- lapply(seq_len(ncol(machine_range)), function(z){
      c(min(machine_range[, z]), max(machine_range[, z]))
    })
    machine_range <- do.call("cbind", machine_range)
    colnames(machine_range) <- channels
    rownames(machine_range) <- c("min", "max")
    # REPLACE MAX DATA RANGE 
    data_range["max", ] <- machine_range["max", ]
  }
  
  # REPLACE LOWER LIMIT IF > 0 - AUTO
  if(axes_limits != "data" & anchor == TRUE){
    if(any(data_range[1,] > 0)){
      data_range[1, data_range[1,] > 0] <- 0
    }
  }
  
  # GATES ----------------------------------------------------------------------
  
  # MAKE SURE ALL GATES ARE VISIBLE
  if(!.all_na(gate)) {
    # GATE CO-ORDINATES
    coords <- .cyto_gate_coords(gate, channels)
    # DROP NON-FINITE CO-ORDINATES
    coords <- coords[is.finite(coords[, 1]), , drop = FALSE]
    # BYPASS INFINITE CO-ORDINATES
    if(nrow(coords) > 0) {
      # MINIMUM
      if(min(coords) < min(data_range)) {
        data_range[1, 1] <- min(coords)
      }
      # MAXIMUM
      if(max(coords) > max(data_range)) {
        data_range[2, 1] <- max(coords)
      }
    }
  }
  
  # BUFFER ---------------------------------------------------------------------
  
  # ADD 2% BUFFER EITHER SIDE - 4% TOTAL IN PLOT
  if(buffer != 0){
    # BUFFER
    for(z in channels) {
      # RANGE
      RNG <- data_range[, z][2] - data_range[, z][1]
      # ADD BUFFER
      data_range[, z][1] <- data_range[, z][1] - buffer * RNG
      data_range[, z][2] <- data_range[, z][2] + buffer * RNG
    }
  }
  
  return(data_range)
  
}

## AXES TEXT ----

#' Axes ticks and text
#'
#' @param x list of cytoframes/cytosets.
#' @param channels name(s) of the channel(s) used to construct the plot.
#' @param axes_trans transformerList.
#' @param axes_range named list of axes limits for each each axis (i.e.
#'   list(xlim,ylim)).
#' @param axes_limits either "auto", "data" or "machine".
#'
#' @return list containing axis labels and breaks.
#'
#' @importFrom grDevices axisTicks
#'
#' @noRd
.cyto_plot_axes_text <- function(x,
                                 channels,
                                 axes_trans = NA,
                                 axes_range = list(NA, NA),
                                 axes_limits = "data",
                                 axes_limits_buffer = 0.03,
                                 rescale = c(NA, NA),
                                 format = FALSE) {
  
  # AXES_RANGE IS PADDED BUT RESCALE IS NOT

  # AXES RANGES ----------------------------------------------------------------
  
  # PAD AXES RANGES - 4% BASE GRAPHICS
  axes_range <- structure(
    lapply(axes_range, function(z){
      if(!.all_na(z)) {
        # RANGE
        rng <- diff(z)
        pad <- (rng - rng / 1.04) / 2
        return(
          c(z[1] - pad,
            z[2] + pad)
        )
      } else {
        return(z)
      }
    }
    ), 
    names = names(axes_range)
  )
  
  # AXES TICKS & LABELS --------------------------------------------------------
  
  # LOOP THROUGH CHANNELS
  axes_text <- lapply(seq_along(channels), function(z) {
    # CHANNEL
    chan <- channels[z]
    # AXES RANGE
    rng <- axes_range[[z]]
    # NO CHANNEL OR AXES RANGE SUPPLIED - BYPASS
    if(.all_na(chan) & .all_na(rng)) {
      return(NA)
    }
    # LINEAR SCALE
    if(!chan %in% names(axes_trans)) {
      # COMPUTE MAJOR TICKS
      axis_major_ticks <- axisTicks(
        rng,
        log = FALSE
      )
      names(axis_major_ticks) <- trimws(
        format(
          unname(axis_major_ticks), 
          scientific = FALSE
        )
      )
      # COMPUTE MINOR TICKS - WITHIN MAJOR TICK RANGE
      axis_minor_ticks <- LAPPLY(
        seq_len(length(axis_major_ticks) - 1),
        function(w){
          pretty(
            axis_major_ticks[c(w, w + 1)],
            n = 5
          )
        }
      )
      # COMPUTE MINOR TICKS - OUTSIDE MAJOR TICK RANGE
      axis_minor_ticks <- c(
        seq(
          min(axis_major_ticks),
          min(rng),
          diff(axis_minor_ticks[2:1])
        ),
        axis_minor_ticks,
        seq(
          max(axis_major_ticks),
          max(rng),
          diff(axis_minor_ticks[1:2])
        )
      )
      axis_minor_ticks <- axis_minor_ticks[
        !axis_minor_ticks %in% axis_major_ticks
      ]
      names(axis_minor_ticks) <- rep("", length(axis_minor_ticks))
      # COMBINED AXIS TICKS
      axis_ticks <- sort(
        c(
          axis_minor_ticks,
          axis_major_ticks
        )
      )
      # RESCALE TICK LOCATIONS - KEY MAPPING TO Y AXIS
      if(!.all_na(rescale)) {
        # RESCALE TICK LOCATIONS
        axis_ticks <- structure(
          LAPPLY(
            axis_ticks,
            function(v) {
              min(rescale) + ((v - min(rng))/diff(rng)) * diff(rescale)
            }
          ),
          names = names(axis_ticks)
        )
      }
      # FORMAT LABELS - NOT MAPPED TO AXIS - ABBREVIATE LABELS
      if(format) {
        names(axis_ticks) <- LAPPLY(
          names(axis_ticks),
          function(v) {
            if(!.empty(v)) {
              v <- as.numeric(v)
              if(v / 1000 >= 1) {
                if(v / 1000000 >= 1) {
                  return(
                    paste0(
                      round(v / 1000000, 1),
                      "M"
                    )
                  )
                } else {
                  return(
                    paste0(
                      round(v / 1000, 1),
                      "K"
                    )
                  )
                }
              } else {
                return(v)
              }
            } else {
              return(v)
            }
          }
        )
      }
      # RETURN PREPARED AXIS LOCATIONS & LABELS
      return(
        list(
          "label" = names(axis_ticks),
          "at" = axis_ticks
        )
      )
    # TRANSFORMED SCALE
    } else {
      # AXIS TICKS
      axis_ticks <- c(
        sort(
          LAPPLY(
            1 * 10^seq_len(9),
            function(v) {
              -seq(90, 10, -10) * v
            }
          )
        ),
        seq(-9, 9, 1),
        LAPPLY(
          1 * 10^seq_len(9),
          function(v) {
            seq(10, 90, 10) * v
          }
        )
      )
      # AXIS LABELS
      axis_labels <- lapply(axis_ticks, function(v) {
        if (v != 0) {
          pwr <- log10(abs(v))
        }
        if (v == 0) {
          quote(0)
        } else if (pwr == 0) {
          quote("")
        } else if (abs(pwr) %% 1 == 0) {
          if(v < 0) {
            substitute(-10^pwr)
          } else {
            substitute(10^pwr)
          }
        } else {
          quote("")
        }
      })
      # AXES TEXT
      trans_func <- axes_trans[[chan]]$transform
      inv_func <- axes_trans[[chan]]$inverse
      # AXES RANGE ON LINEAR SCALE
      if (!.all_na(rng)) {
        rng <- inv_func(rng)
      } else {
        rng <- inv_func(
          .cyto_plot_axes_limits(
            x,
            channels = chan,
            axes_limits = axes_limits,
            buffer = axes_limits_buffer
          )[, chan])
      }
      # RESTRICT AXES TICKS & LABELS BY RANGE
      ind <- which(axis_ticks > min(rng) & axis_ticks < max(rng))
      axis_ticks <- axis_ticks[ind]
      axis_labels <- axis_labels[ind]
      # BREAKS - TRANSFORMED SCALE
      axis_ticks <- signif(trans_func(axis_ticks))
      # RANGE - TRANSFORMED SCALE
      rng <- trans_func(rng)
      # REMOVE OVERLAPPING LABELS AROUND ZERO 
      if("0" %in% axis_labels) {
        # BUFFERING AROUND ZERO
        zero_ind <- match_ind("0", axis_labels)
        zero_break <- axis_ticks[zero_ind]
        zero_rng <- c(
          zero_break - 0.025 * diff(rng),
          zero_break + 0.025 * diff(rng)
        )
        # TICKS WITHIN RANGE - EMPTY LABELS
        ind <- which(
          axis_ticks > min(zero_rng) &
            axis_ticks < max(zero_rng)
        )
        lapply(ind, function(q){
          axis_labels[[q]] <<- ""
        })
        axis_labels[[zero_ind]] <- "0"
      }
      axis_labels <- do.call("expression", axis_labels)
      # RESCALE TICK LOCATIONS - KEY MAPPING TO Y AXIS
      if(!.all_na(rescale)) {
        # RESCALE TICK LOCATIONS
        axis_ticks <- LAPPLY(
          axis_ticks,
          function(v) {
            min(rescale) + 
              ((v - min(rng))/diff(rng)) * diff(rescale)
          }
        )
      }
      # BREAKS & LABELS
      return(
        list(
          "label" = axis_labels, 
          "at" = axis_ticks
        )
      )
    }
  })
  
}

## AXES LABELS ----

#' Get axes titles for cyto_plot
#'
#' @param x list of cytoframes/cytosets.
#' @param channels used to construct the plot.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param hist_stat "percent", "count" or "density".
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_axes_label <- function(x,
                                  channels,
                                  xlab,
                                  ylab,
                                  hist_stat = "count") {
  
  # CHANNELS
  channels <- c(channels, NA)[1:2]
  
  # MARKERS
  markers <- cyto_markers(x[[1]])
  
  # AXES_LABELS
  cnt <- 0
  axes_labels <- lapply(channels, function(z) {
    cnt <<- cnt + 1
    # LABELS
    if (cnt == 1) {
      lab <- xlab
    } else {
      lab <- ylab
    }
    # XLAB OR SCATTER YLAB
    if (!.all_na(z)) {
      # XLAB/YLAB
      if (missing(lab) | .empty(lab)) {
        # MARKER ASSIGNED - DIFFERENT FROM CHANNEL
        ind <- match_ind(z, names(markers))
        if (length(ind) > 0 & !z %in% markers) {
          return(paste(markers[ind], z))
        } else {
          return(z)
        }
      } else if (.all_na(lab)) {
        return(NA)
      }
      # HISTOGRAM YLAB
    } else {
      if (missing(lab) | .empty(lab)) {
        # COUNT
        if (hist_stat == "count") {
          return("Count")
          # PERCENT
        } else if (hist_stat == "percent") {
          return("% of Mode")
          # DENSITY
        } else if (hist_stat == "density") {
          return("Density")
        }
      } else if (.all_na(lab)) {
        return(NA)
      }
    }
  })
  names(axes_labels) <- c("xlab", "ylab")
  
  # RETURN AXES LABELS
  return(axes_labels)
}

## LAYOUT/MARGINS --------------------------------------------------------------

## LAYOUT ----

#' Set plot layout
#'
#' @param x list of  cytoframe/cytoset lists.
#' @param layout grid dimensions c(nr, nc), NA or FALSE.
#'
#' @importFrom grDevices n2mfrow
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_layout <- function(x,
                              layout = NULL) {
  
  # LAYOUT
  if (is.null(layout) | .empty(layout)) {
    layout <- rev(n2mfrow(length(x)))
  }
  return(layout)
  
}

## MARGINS ----

#' Set plot margins
#'
#' @param x list of flowFrames or density objects to plot.
#' @param legend logical indicating whether a legend should be included in the
#'   plot.
#' @param title if NULL remove excess space above plot.
#' @param axes_text vector of logicals indicating whether the x and y axes
#'   should be included on the plot.
#' @param margins a vector of length 4 to control the margins around the bottom,
#'   left, top and right of the plot, set to c(NA, NA, NA, NA) by default to let
#'   `cyto_plot` compute optimal margins.
#' @param point_col required to make space for point_col_scale legend.
#'
#' @importFrom methods is
#'
#' @noRd
.cyto_plot_margins <- function(x,
                               channels = NULL,
                               legend = FALSE,
                               legend_text = NA,
                               legend_text_size = 1,
                               title,
                               axes_text = list(TRUE, TRUE),
                               margins = c(NA, NA, NA, NA),
                               point_col = NA,
                               key = "both",
                               key_scale = "fixed") {
  
  # ARGUMENTS
  args <- .args_list()

  # MARGINS
  if (is.null(margins) | .empty(margins) | .all_na(margins)) {
    
    # DEFAULT
    mar <- c(5.1, 5.1, 4.1, 2.1)
    
    # KEY 
    if(length(channels) == 2 & 
       (any(is.na(point_col)) | 
        any(point_col %in% c(cyto_channels(x), cyto_markers(x)))) &
       !key %in% "none") {
      # KEY SCALE & TEXT
      if(key %in% "both") {
        mar[4] <- mar[4] + 2.9
      } else {
        mar[4] <- mar[4] + 1
      }
    }
    
    # LEGEND TEXT
    if (legend != FALSE & !.all_na(legend_text)) {
      # KEY REQUIRED
      if(length(channels) == 2 & 
         (any(is.na(point_col)) | 
          any(point_col %in% c(cyto_channels(x), cyto_markers(x)))) &
         !key %in% "none") {
        # KEY SCALE & TEXT
        if(key %in% "both" & length(channels) == 2) {
          mar[4] <- 8.8 + max(nchar(legend_text))*0.32*mean(legend_text_size)
        # KEY SCALE ONLY
        } else {
          mar[4] <- 7.5 + max(nchar(legend_text))*0.32*mean(legend_text_size)
        }
      # NO KEY
      } else {
        mar[4] <- 7 + max(nchar(legend_text))*0.32*mean(legend_text_size)
      }
    }
    
    # TITLE
    if (.all_na(title)) {
      mar[3] <- 2.2
    }
    
    # X AXIS
    if (!all(is(axes_text[[1]], "list"))) {
      if (.all_na(axes_text[[1]])) {
        # NA == FALSE returns NA not T/F
      } else if (all(axes_text[[1]] == FALSE)) {
        mar[1] <- 4.1
      }
    }
    
    # Y AXIS
    if (!all(is(axes_text[[2]], "list"))) {
      if (.all_na(axes_text[[2]])) {
        # NA == FALSE return NA not T/F
      } else if (all(axes_text[[2]] == FALSE)) {
        mar[2] <- 4.1
      }
    }
  } else {
    if (length(margins) != 4) {
      stop("'margins' must be a vector with 4 elements.")
    }
    mar <- margins
  }
  
  # SET MAR PARAMETER
  cyto_plot_par(mar = mar)

}

## LABELS ----------------------------------------------------------------------

## LABEL_TEXT & LABEL_STAT ARGUMENTS ----

#' Prepare label_text and label_stat arguments
#' @noRd
.cyto_plot_label_args <- function(x,
                                  channels = NULL,
                                  label = "",
                                  label_text = "",
                                  label_stat = "",
                                  gate = NA,
                                  hist_stack = 0,
                                  ...) {
  
  
  if(any(is.function(label_stat))) {
    stop("label_stat must be a character string!")
  }
  
  # POPULATIONS PER LAYER
  NP <- .cyto_gate_count(gate)
  
  # TOTAL POPULATIONS
  SMP <- length(x)
  TNP <- NP * SMP
  TNP_split <- split(seq_len(TNP), rep(seq_len(SMP), each = NP))
  
  # LABEL_STAT FILLED WITH EMPTY
  if (!all(LAPPLY(label_stat, ".empty")) &
      any(LAPPLY(label_stat, ".empty"))) {
    label_stat[LAPPLY(label_stat, ".empty")] <- NA
  }
  
  # LABEL_STAT
  # 1D PLOT NO STACK
  if (length(channels) == 1 & hist_stack == 0) {
    # LABEL_STAT MISSING
    if (all(LAPPLY(label_stat, ".empty"))) {
      # GATE - FREQ STAT
      if (!.all_na(gate)) {
        # LABEL_STAT - BASE LAYER ONLY
        label_stat <- c(
          rep("freq", NP),
          rep(NA, TNP - NP)
        )
        # NO GATE - NO STAT
      } else {
        # LABEL_STAT REMOVED
        label_stat <- rep(NA, TNP)
      }
      # LABEL_STAT SUPPLIED - FILL WITH NA
    } else {
      # GATE - BASE LAYER ONLY
      if (!.all_na(gate)) {
        if (length(label_stat) == 1) {
          label_stat <- rep(label_stat, length.out = NP)
        }
        label_stat <- rep(c(
          label_stat,
          rep(NA, length.out = TNP)
        ),
        length.out = TNP
        )
        # NO GATE
      } else {
        # LABEL EACH LAYER
        label_stat <- rep(label_stat, length.out = TNP)
      }
    }
    # 1D PLOT STACK
  } else if (length(channels) == 1 & hist_stack != 0) {
    # LABEL_STAT MISSING
    if (all(LAPPLY(label_stat, ".empty"))) {
      # GATE - FREQ STAT
      if (!.all_na(gate)) {
        # LABEL_STAT - ALL LAYERS
        label_stat <- rep("freq", length.out = TNP)
        # NO GATE
      } else {
        # LABEL_STAT REMOVED
        label_stat <- rep(NA, length.out = TNP)
      }
      # LABEL_STAT SUPPLIED - FILL WITH NA
    } else {
      label_stat <- rep(label_stat, length.out = TNP)
    }
    # 2D PLOT
  } else if (length(channels) == 2) {
    # LABEL_STAT MISSING
    if (all(LAPPLY(label_stat, ".empty"))) {
      # GATE - FREQ STAT
      if (!.all_na(gate)) {
        # LABEL_STAT - BASE LAYER ONLY
        label_stat <- c(
          rep("freq", NP),
          rep(NA, TNP - NP)
        )
        # NO GATE - NO STAT
      } else {
        # LABEL_STAT REMOVED
        label_stat <- rep(NA, length.out = TNP)
      }
      # LABEL_STAT SUPPLIED- FILL WITH NA
    } else {
      # GATE - BASE LAYER ONLY
      if (!.all_na(gate)) {
        if (length(label_stat) == 1) {
          label_stat <- rep(label_stat, length.out = NP)
        }
        label_stat <- rep(c(
          label_stat,
          rep(NA, length.out = TNP)
        ),
        length.out = TNP
        )
        # NO GATE
      } else {
        # LABEL EACH LAYER
        label_stat <- rep(label_stat, length.out = TNP)
      }
    }
  }
  
  # LABEL_TEXT
  if(any(LAPPLY(label_text, ".empty"))) {
    # POPULATION NAMES FROM GATES - WATCH OUT QUADGATES
    if(!.all_na(gate)) {
      pops <- rep(
        LAPPLY(
          seq_along(gate),
          function(z) {
            # USE NAMES OF GATE LIST
            if(!is.null(names(gate)[z])) {
              unlist(
                strsplit(
                  names(gate)[z],
                  "\\|"
                )
              )
            # USE FILTERIDs
            } else {
              unlist(
                strsplit(
                  gate[[z]]@filterId,
                  "\\|"
                )
              )
            }
          }
        ),
        length(x)
      )
      label_text <- LAPPLY(
        seq_along(label_text),
        function(z) {
          if(!.all_na(label_stat[z])) {
            if(.empty(label_text[z])) {
              return(pops[z])
            } else {
              return(label_text[z])
            }
          } else {
            return(NA)
          }
        }
      )
    # REMOVE LABEL_TEXT
    } else {
      label_text[LAPPLY(label_text, ".empty")] <- NA
    }
  }
  
  # LABEL
  if (all(LAPPLY(label, ".empty"))) {
    # TURN LABELS ON
    if (!.all_na(c(label_text, label_stat))) {
      label <- TRUE
      # TURN LABELS OFF
    } else {
      label <- FALSE
    }
  }
  
  # RETURN LABEL ARGUMENTS
  return(list("label" = label,
              "label_text" = label_text,
              "label_stat" = label_stat))
  
}

## LABEL POPULATIONS ----

#' Get a list of gated populations to label
#'
#' @param x list of cytoframes/cytosets per layer to be gated
#' @param gate list of gate objects to apply to each element of x, gates only
#'   required for base layer.
#' @param label_stat required to identify which layers require gating to get
#'   statistics
#' @param ... not in use.
#'
#' @return list of flowFrames lists
#'
#' @importFrom flowWorkspace flowSet_to_cytoset
#' @importFrom flowCore Subset split quadGate
#' @importFrom methods is
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_pops <- function(x,
                                  gate = NA,
                                  label = TRUE,
                                  label_stat = NA,
                                  ...) {
  
  # FLOWFRAME LIST
  if(cyto_class(x, c("flowFrame", "flowSet"))){
    x <- list(x)
  }
  
  # NO GATES -------------------------------------------------------------------
  
  # RETURN X
  if (.all_na(gate) | label == FALSE) {
    pops <- lapply(seq_len(length(x)), function(z){
      if(!.all_na(gate)) {
        rep(x[z], length(gate))
      } else {
        x[z]
      }
    })
    return(pops)
  }
  
  # PREPARE GATES --------------------------------------------------------------
  
  # LIST OF GATE OBJECT LISTS
  if(cyto_class(gate, "list") &
     all(LAPPLY(gate, "cyto_class", "list", TRUE))){
    # USE BASE LAYER GATES
    gate <- gate[[1]]
  }
  
  # LABEL_STAT -----------------------------------------------------------------
  
  # LAYERS
  layers <- length(x)
  
  # LABELS PER LAYER (SAME LENGTH AS GATE WHICH INCLUDES NEGATE)
  labels_per_layer <- split(label_stat,
                            rep(seq_len(layers), 
                                each = length(label_stat)/layers))
  
  # POPULATIONS ----------------------------------------------------------------
  
  # CAREFUL CANNOT NEGATE INDIVIDUAL QUADRANTS EITHER
  
  # ARGUMENTS
  args <- .args_list()
  
  # GATING PER LAYER (list of pops per layer)
  pops <- lapply(seq_along(x), function(z){
    cs <- x[[z]]
    # NO LABEL - NO GATING REQUIRED (CHECK WHOLE LAYER)
    if(.all_na(labels_per_layer[[z]])) {
      return(rep(list(cs), length(label_stat)/layers))
    # LABEL - GATING REQUIRED
    } else {
      # LIST OF GATED POPULATIONS PER GATE -> LIST OF POPS (CYTOSETS)
      unlist(
        cyto_gate_apply(cs,
                        gate = gate)
      )
    }
  })
  
  # RETURN LIST OF GATED POPULATIONS
  return(pops)
}

## LABEL STATISTICS ----

#' Compute and prepare statistics for labels
#'
#' @param x list of parental flowFrame objects.
#' @param pops list of population lists to label.
#' @param channels vector of channels used to construct the plot.
#' @param label_stat names of statistics to include in labels, supplied per
#'   layer.
#'
#' @return computed statistics to include in labels.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_stat <- function(x,
                                  pops,
                                  channels,
                                  axes_trans = NA,
                                  label_stat,
                                  gate = NA,
                                  hist_smooth = 1,
                                  hist_bins = 256,
                                  xlim = c(NA,NA),
                                  ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # NO LABEL_STAT
  if (.all_na(label_stat)) {
    return(label_stat)
  }
  
  # FLOWFRAME LIST
  if(cyto_class(x, c("flowFrame", "flowSet"))){
    x <- list(x)
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # LABEL_STAT
  label_stat <- split(label_stat,
                      rep(seq_len(length(x)),
                          each  = length(pops[[1]])))
  
  # COMPUTE STAT AGAINST BASE LAYER (NO GATES)
  if(.all_na(gate)){
    x <- rep(x[1], length(x))
  }
  
  # COMPUTE STATISTICS ---------------------------------------------------------
  
  # STATISTICS
  LABEL_STAT <- LAPPLY(seq_len(length(x)), function(z) {
    # LABEL_STAT
    ST <- lapply(seq_len(length(pops[[z]])), function(y) {
      # STATISTIC SUPPLIED
      if (!.all_na(label_stat[[z]][y])) {
        # DISPATCH
        label_stat_fun <- .cyto_stat_dispatch(label_stat[[z]][y])
        # FREQUENCY
        if(grepl("freq$", label_stat_fun, ignore.case = TRUE) |
           grepl("percent$", label_stat_fun, ignore.case = TRUE)) {
          res <- cyto_stats_compute(pops[[z]][[y]],
                                    channels = channels,
                                    parent = x[z],
                                    stat = label_stat_fun,
                                    round = 2,
                                    input = "matrix",
                                    markers = FALSE,
                                    details = FALSE)
        # HISTOGRAM STATISTICS - ADDITIONAL ARGUMENTS
        } else if(any(LAPPLY(c("mode",
                               "auc"), function(w){
                                 grepl(paste0(w, "$"), 
                                       label_stat_fun, 
                                       ignore.case = TRUE)
                               }))){
          res <- cyto_stats_compute(pops[[z]][[y]],
                                    channels = channels,
                                    stat = label_stat_fun,
                                    inverse = TRUE,
                                    trans = axes_trans,
                                    round = 2,
                                    input = "matrix",
                                    smooth = hist_smooth,
                                    bins = hist_bins,
                                    limits = xlim,
                                    markers = FALSE,
                                    details = FALSE)
        # OTHER STATISTICS - GATED POPS
        } else {
          res <- cyto_stats_compute(pops[[z]][[y]],
                                    channels = channels,
                                    stat = label_stat_fun,
                                    inverse = TRUE,
                                    trans = axes_trans,
                                    round = 2,
                                    input = "matrix",
                                    markers = FALSE,
                                    details = FALSE)
        }
        # DROP NAME & ALIAS COLUMNS
        res <- res[, -which(c("name", "alias") %in% colnames(res)), 
                   drop = FALSE]
        # EXTRACT STATISTICS - EITHER SINGLE OR DOUBLE
        res <- res[1, , drop = TRUE]
        # ROUND 2 DECIMAL PLACES
        if(!grepl("count$", label_stat_fun)) {
          res <- .round(res)
        }
        # STATSTICS REQUIRE %
        if(any(LAPPLY(c("freq",
                        "percent",
                        "cv",
                        "rcv"), function(w){
                          grepl(paste0(w, "$"), 
                                label_stat_fun, 
                                ignore.case = TRUE)
                        }))) {
          res <- paste(res, "%")
        }
        # STATISTIC PER CHANNEL - 2D PLOT
        if(length(res) > 1) {
          res <- paste("x =", res[1], "\n", "y =", res[2])
        }
      } else {
        res <- NA
      }
      return(res)
    })
    return(ST)
  })
  
  # RETURN COMPUTED STATISTICS -------------------------------------------------
  return(LABEL_STAT)
}

## PREPARE LABEL_TEXT ----

#' Prepare text for labels to include statistics
#'
#' @param label_text text to include in the labels.
#' @param label_stat vector of computed statistics to include in the labels.
#'
#' @return vector of finalised labels incorporating both label_text and
#'   label_stat.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_text <- function(label_text,
                                  label_stat,
                                  ...) {
  
  # MERGE LABEL_TEXT & LABEL_STAT
  for(z in seq_len(length(label_text))) {
    # NO LABEL_TEXT & NO LABEL_STAT
    if (.all_na(label_text[z]) & .all_na(label_stat[z])) {
      label_text[z] <- NA
      # NO LABEL_TEXT & LABEL_STAT
    } else if (.all_na(label_text[z]) & !.all_na(label_stat[z])) {
      label_text[z] <- label_stat[z]
      # LABEL_TEXT & NO LABEL_STAT
    } else if (!.all_na(label_text[z]) & .all_na(label_stat[z])) {
      label_text[z] <- label_text[z]
      # LABEL_TEXT & LABEL_STAT
    } else if (!.all_na(label_text[z]) & !.all_na(label_stat[z])) {
      label_text[z] <- paste(label_text[z], label_stat[z], sep = "\n")
    }
  }
  
  # RETURN LABEL_TEXT
  return(label_text)
}

## LABEL CO-ORDINATES ----

#' Compute offset label co-ordinates
#'
#' Used internally within cyto_plot to compute offset label co-ordinates. Only
#' called if label_position and label are set to TRUE.
#'
#' @param args list of named cyto_plot arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @noRd
.cyto_plot_label_coords <- function(x,
                                    channels = NULL,
                                    pops = NULL,
                                    axes_limits = c(NA, NA),
                                    gate = NA,
                                    label_text = NA,
                                    label_text_x = NA,
                                    label_text_y = NA,
                                    label_text_size = 1.8,
                                    hist_stack = 0,
                                    hist_layers = NA,
                                    hist_smooth = 1,
                                    hist_bins = 256,
                                    d = NULL,
                                    ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # X - CYTO_PLOT ARGUMENTS
  if(class(x) == "cyto_plot") {
    .args_update(x)
  }
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04) / 2 # 2% BUFFER
  xmin <- xmin + xpad
  xmax <- xmax - xpad
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04) / 2 # 2% BUFFER
  ymin <- ymin + ypad
  ymax <- ymax - ypad
  yrng <- ymax - ymin
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # LAYERS & POPULATIONS
  L <- length(pops)
  GNP <- length(pops[[1]])
  GC <- .cyto_gate_count(gate, 
                         drop = TRUE,
                         total = TRUE)
  
  # PREPARE LABEL ARGUMENTS PER LAYER
  label_text <- split(label_text, rep(seq_len(L), each = GNP))
  label_text_x <- split(label_text_x, rep(seq_len(L), each = GNP))
  label_text_y <- split(label_text_y, rep(seq_len(L), each = GNP))
  label_text_size <- split(label_text_size, rep(seq_len(L), each = GNP))
  
  # hist_stack --------------------------------------------------------------
  
  # 1D PLOT
  if (length(channels) == 1) {
    # STACKING
    if (hist_stack != 0 &
        ifelse(.all_na(hist_layers), TRUE, hist_layers != 1)) {
      stk <- LAPPLY(d, function(z){
        z$range
      })
      y_coords <- stk + 0.5 * (
        d[[2]]$range[1] - d[[1]]$range[1]
      )
      # REPEAT (Y_COORDS/LAYER)
      y_coords <- rep(y_coords, each = GNP)
      # NO STACKING
    } else {
      # DEFAULT Y COORD - 50% Y RANGE
      y_coords <- rep(0.5 * yrng, L * GNP)
    }
    y_coords <- split(y_coords, rep(seq_len(L), each = GNP))
  }
  
  # COMPUTE LABEL CO-ORDINATES -------------------------------------------------
  
  # GATE CENTERS - (IGNORE SUPPLIED LABEL COORDS)
  if (!.all_na(gate)) {
    gate_centers <- do.call(
      "rbind",
      .cyto_gate_center(gate,
                        channels = channels,
                        text_x = rep(NA, GC),
                        text_y = rep(NA, GC))
    )
  # NO GATES - NO CENTERS
  } else {
    gate_centers <- matrix(NA,
                           ncol = 2,
                           nrow = sum(GNP),
                           dimnames = list(NULL, c("x", "y")))
  }
  
  # LABEL_TEXT_XY - MATRIX
  label_text_xy <- lapply(seq_len(L), function(z) {
    # TEMPORARY STORAGE VECTORS
    text_x <- c()
    text_y <- c()
    # COMPUTE LABEL CO-ORDINATES
    for(y in seq_len(GNP)) {
      # LABEL
      if (!.all_na(label_text[[z]][y])) {
        # ID PLOT - CENTER OF RANGE
        if (length(channels) == 1) {
          # GATE CENTER
          gate_center <- gate_centers[y, ]
          # X COORD REQUIRED
          if(.all_na(label_text_x[[z]][y])) {
            # GATE CENTER AVAILABLE
            if(!.all_na(gate_center[1])) {
              text_x[y] <- gate_center[1]
              # NO GATE CENTER
            } else {
              # NO EVENTS
              if (cyto_apply(pops[[z]][[y]],
                             "nrow",
                             input = "matrix",
                             copy = FALSE)[, 1] == 0) {
                # RANGE CENTER - PLOT LIMITS
                text_x[y] <- mean(c(xmin, xmax))
              } else {
                # RANGE
                rng <- suppressMessages(
                  .cyto_plot_axes_limits(pops[[z]][[y]],
                                         channels = channels[1],
                                         axes_limits = axes_limits,
                                         buffer = 0,
                                         anchor = FALSE
                  )[, channels]
                )
                # UNIMODAL - 50% RANGE
                if (abs(diff(rng)) <= 0.6 * xrng) {
                  text_x[y] <- quantile(rng, 0.5)
                  # UMULTIMODAL - 56% RANGE
                } else {
                  text_x[y] <- quantile(rng, 0.56)
                }
              }
            }
            # X COORD - MANUALLY SUPPLIED
          } else {
            text_x[y] <- label_text_x[[z]][y]
          }
          # Y COORD - STACKS/LIMITS
          if (.all_na(label_text_y[[z]][y])) {
            text_y[y] <- y_coords[[z]][y]
            # Y COORD MANUALLY SUPPLIED
          } else if (!.all_na(label_text_y[[z]][y])) {
            text_y[y] <- label_text_y[[z]][y]
          }
          # 2D PLOT - MODE
        } else if (length(channels) == 2) {
          # GATE CENTER
          gate_center <- gate_centers[y, ]
          # X COORD REQUIRED
          if(.all_na(label_text_x[[z]][y])) {
            # GATE CENTER AVAILABLE
            if(!.all_na(gate_center[1])) {
              text_x[y] <- gate_center[1]
              # NO GATE CENTER
            } else {
              # NO EVENTS
              if (cyto_apply(pops[[z]][[y]],
                             "nrow",
                             input = "matrix",
                             copy = FALSE)[, 1] < 2) {
                # RANGE CENTER
                text_x[y] <- mean(c(xmin, xmax))
              } else {
                # MODE
                text_x[y] <- cyto_apply(pops[[z]][[y]],
                                        "cyto_stat_mode",
                                        input = "matrix",
                                        channels = channels[1],
                                        smooth = hist_smooth,
                                        bins = hist_bins,
                                        inverse = FALSE)[, 1]
              }
            }
          # X COORD - MANUALLY SUPPLIED
          } else {
            text_x[y] <- label_text_x[[z]][y]
          }
          # Y COORD REQUIRED
          if(.all_na(label_text_y[[z]][y])) {
            # GATE CENTER AVAILABLE
            if(!.all_na(gate_center[2])) {
              text_y[y] <- gate_center[2]
            # NO GATE CENTER
            } else {
              # NO EVENTS
              if (cyto_apply(pops[[z]][[y]],
                             "nrow",
                             input = "matrix",
                             copy = FALSE)[, 1] == 0) {
                # RANGE CENTER
                text_y[y] <- mean(c(ymin, ymax))
              } else {
                # MODE
                text_y[y] <- cyto_apply(pops[[z]][[y]],
                                        "cyto_stat_mode",
                                        input = "matrix",
                                        channels = channels[2],
                                        smooth = hist_smooth,
                                        bins = hist_bins,
                                        inverse = FALSE)[, 1]
              }
            }
          # Y COORD - MANUALLY SUPPLIED
          } else {
            text_y[y] <- label_text_y[[z]][y]
          }
        }
        # NO LABEL
      } else if (.all_na(label_text[[z]][y])) {
        text_x[y] <- NA
        text_y[y] <- NA
      }
    }
    
    # MATRIX
    text_xy <- matrix(c(text_x, text_y),
                      ncol = 2,
                      nrow = GNP,
                      byrow = FALSE
    )
    colnames(text_xy) <- c("x", "y")
    return(text_xy)
  })
  
  # UPDATE LABEL_COORDS
  label_text_x <- lapply(label_text_xy, function(z){
    unlist(z[, "x", drop = TRUE])
  })
  label_text_y <- lapply(label_text_xy, function(z){
    unlist(z[, "y", drop = TRUE])
  })
  
  # OFFSET LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL DIMENSIONS
  label_dims <- lapply(seq_len(L), function(z){
    lapply(seq_len(GNP), function(y){
      # COMPUTE LABEL DIMENSIONS
      if (!.all_na(label_text[[z]][y])) {
        .cyto_plot_label_dims(
          label_text = label_text[[z]][y],
          label_text_x = label_text_x[[z]][y],
          label_text_y = label_text_y[[z]][y],
          label_text_size = label_text_size[[z]][y]
        )
      } else {
        matrix(rep(NA, 4),
               ncol = 2,
               dimnames = list(NULL, c("x", "y"))
        )
      }
    })
  })
  
  # OFFSET BY LAYER
  if (length(channels) == 1 & hist_stack != 0) {
    # OFFSET PER LAYER
    for(z in seq_len(L)) {
      # LABEL OVERLAP
      if (.cyto_plot_label_overlap(label_dims[[z]])) {
        # LABEL HEIGHT - OFFSETTING
        label_height <- max(LAPPLY(label_dims[[z]], function(y) {
          max(y[, "y"]) - min(y[, "y"])
        }), na.rm = TRUE)
        # LABEL HEIGHT BUFFERING
        label_height <- 1.18 * label_height
        # Y COORDS TO OFFSET
        text_y <- label_text_y[[z]]
        # OFFSET Y CO-ORDINATES (EXCLUDE NA)
        label_text_y[[z]][!is.na(text_y)] <- tryCatch(
          {
            .suppress_all_messages(
              .spread.labels(text_y[!is.na(text_y)],
                             mindiff = label_height,
                             min = ymin,
                             max = ymax
              )
            )
          },
          error = function(e) {
            return(text_y[!is.na(text_y)])
          }
        )
      }
    }
    # OFFSET ALL LABELS
  } else {
    # COLLAPSE LABEL_DIMS
    label_dims <- unlist(label_dims, recursive = FALSE)
    # LABEL OVERLAP
    if (.cyto_plot_label_overlap(label_dims)) {
      # LABEL HEIGHT - OFFSETTING
      label_height <- max(LAPPLY(label_dims, function(y) {
        max(y[, "y"]) - min(y[, "y"])
      }), na.rm = TRUE)
      # LABEL HEIGHT BUFFERING
      label_height <- 1.18 * label_height
      # OFFSET Y CO-ORDINATES (EXCLUDE NA)
      text_y <- unlist(label_text_y)
      text_y[!is.na(text_y)] <-
        .suppress_all_messages(
          .spread.labels(text_y[!is.na(text_y)],
                         mindiff = label_height,
                         min = ymin,
                         max = ymax
          )
        )
      # UPDATE LABEL_TEXT_Y
      label_text_y <- split(text_y,
                            rep(seq_len(L), each = GNP))
    }
  }
  
  # RETURN COMPUTED LABEL CO-ORDINATES -----------------------------------------
  
  # LABEL CO-ORDINATE MATRIX
  label_text_xy <- matrix(c(unlist(label_text_x), 
                            unlist(label_text_y)),
                          ncol = 2,
                          byrow = FALSE
  )
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
}

## LABEL DIMENSIONS ----

#' Compute label dimensions
#'
#' Labels should already contain statistic as well. Co-ordiante for each label
#' must already be computed.
#'
#' @importFrom graphics strwidth strheight
#'
#' @return upper left and bottom right x and y co-ordinates of labels
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_dims <- function(label_text,
                                  label_text_x,
                                  label_text_y = NA,
                                  label_text_size = 1,
                                  ...) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # RESET CEX & XPD
  old_pars <- par(c("cex", "xpd"))
  on.exit(par(old_pars))
  
  # SET CEX & XPD
  par(cex = label_text_size)
  par(xpd = TRUE)
  
  # LABEL DIMENSIONS -----------------------------------------------------------
  
  # ARGUMENTS
  xpad <- 1.2
  ypad <- 1.2
  adj <- 0.5
  
  # BOX ADJUSTMENT
  box_adj <- adj + (xpad - 1) * label_text_size * (0.5 - adj)
  
  # BOX DIMENSIONS
  lwidths <- strwidth(label_text)
  rwidths <- lwidths * (1 - box_adj)
  lwidths <- lwidths * box_adj
  bheights <- theights <- strheight(label_text) * 0.5
  
  # BOX X COORDS
  xr <- label_text_x - lwidths * xpad
  xl <- label_text_x + lwidths * xpad
  
  # BOX Y COORDS
  yb <- label_text_y - bheights * ypad
  yt <- label_text_y + theights * ypad
  
  # LABEL DIMENSIONS MATRIX ----------------------------------------------------
  
  # MATRIX - TOP LEFT THEN BOTTOM RIGHT
  coords <- matrix(c(
    min(c(xl, xr)),
    max(c(yb, yt)),
    max(c(xl, xr)),
    min(c(yb, yt))
  ),
  ncol = 2,
  byrow = TRUE
  )
  colnames(coords) <- c("x", "y")
  
  # RETURN LABEL DIMENSIONS ----------------------------------------------------
  return(coords)
}

## LABEL OVERLAP ----

#' Check if any labels are overlapping.
#'
#' @param x list of label dimensions.
#'
#' @return TRUE or FALSE based on whether any overlapping labels are found.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_overlap <- function(x) {
  
  # For each rectangle in x
  overlaps <- LAPPLY(seq_len(length(x)), function(y) {
    
    # Check if other rectangles overlap
    LAPPLY(seq_len(length(x))[-y], function(z) {
      
      # Co-ordinates of reference label
      x1 <- x[[y]][, "x"]
      y1 <- x[[y]][, "y"]
      
      # Co-ordinates of comparison label
      x2 <- x[[z]][, "x"]
      y2 <- x[[z]][, "y"]
      
      # MISSING LABELS - NO OVERLAP
      if (any(is.na(c(x1, x2, y1, y2)))) {
        return(FALSE)
      }
      
      # X co-ordinates are overlapping
      if (min(x2) >= min(x1) & min(x2) <= max(x1) |
          max(x2) >= min(x1) & max(x2) <= max(x1)) {
        # Y co-ordinates are also overlapping
        if (min(y2) >= min(y1) & min(y2) <= max(y1) |
            max(y2) >= min(y1) & max(y2) <= max(y1)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
      
      # Non-overlapping x and y co-ordinates
      return(FALSE)
    })
  })
  
  # RETURN TRUE OR FALSE
  return(any(overlaps))
}

## GATE LABELS -----------------------------------------------------------------

#' Add labels to gates as they are plotted
#' 
#' @param args list of cyto_plot arguments
#'
#' @noRd
.cyto_plot_gate_label <- function(args) {
  
  # GENERAL --------------------------------------------------------------------
  
  # SAMPLES
  L <- length(args$x)
  
  # POPULATIONS PER LAYER
  NP <- length(args$pops[[1]])
  
  # TOTAL POPULATIONS
  TNP <- L * NP
  
  # REMOVE ANY FILTERS
  args$gate <- structure(
    lapply(args$gate, function(z){
      # KEEP GATE OBJECTS
      if(grepl("gate", cyto_class(z), ignore.case = TRUE)) {
        return(z)
      } else {
        return(NULL)
      }
    }),
    names = names(args$gate)
  )
  args$gate[sapply(args$gate, "is.null")] <- NULL
  
  # GATES
  NG <- length(args$gate)
  
  # POPULATIONS PER GATE
  P <- .cyto_gate_count(args$gate,
                        total = FALSE)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # CYTO_PLOT_GATE & CYTO_PLOT_LABELLER ARGUMENTS
  gate_args <- formalArgs("cyto_plot_gate.list")
  label_args <- formalArgs("cyto_plot_labeller")
  
  # RESTRICT ARGUMENTS
  args <- args[names(args) %in% c(gate_args, "label", label_args)]
  
  # SPLIT GATE_FILL ARGUMENTS BY POPULATIONS
  for(z in gate_args[which(grepl("gate_fill", gate_args))]){
    if(z %in% names(args)) {
      args[[z]] <- split(args[[z]], rep(seq_len(NG), times = P))
    }
  }
  
  # SPLIT GATE_LINE ARGUMENTS BY GATE
  for(z in gate_args[which(grepl("gate_line", gate_args))]){
    if(z %in% names(args)) {
      args[[z]] <- split(args[[z]], seq_len(NG))
    }
  }
  
  # SPLIT LABEL ARGUMENTS BY POPULATIONS PER LAYER
  for(z in label_args) {
    if(z %in% names(args)) {
      args[[z]] <- rep(args[[z]], length.out = TNP)
      args[[z]] <- split(args[[z]], rep(seq_len(L), each = NP))
      # RE-ARRANGE LABEL COORDS PER GATE
      args[[z]] <- lapply(seq_len(NP), function(y){
        LAPPLY(args[[z]], `[[`, y)
      })
    }
  }
  
  # GATE ARGUMENTS
  gate_args <- args[names(args) %in% gate_args]
  
  # LABEL ARGUMENTS
  label_args <- args[names(args) %in% label_args]
  
  # GATE & ASSOCIATED LABELS ---------------------------------------------------
  
  # GATES & LABELS
  label_text_xy <- lapply(seq_len(NP), function(z) {
    # GATED POPULATION(S) - GATE & LABEL(S)
    if (z <= NG) {
      # PLOT GATE
      do.call(
        "cyto_plot_gate",
        c(
          list("channels" = gate_args$channels),
          lapply(
            gate_args[!grepl("channels", names(gate_args))],
            `[[`, z
          )
        )
      )
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (args$label == TRUE & !.all_na(args$label_text[[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(args$label_text[z])),
                          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
      # NEGATED POPULATION(S) - LABEL(s) ONLY
    } else {
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (args$label == TRUE & !.all_na(args$label_text[[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(args$label_text[[z]])),
                          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
    }
    return(text_xy)
  })
  label_text_xy <- do.call("rbind", label_text_xy)
  
  # RE-ARRANGE LABEL ARGUMENTS -------------------------------------------------
  
  # UPDATE LABEL_TEXT_X & LABEL_TEXT_Y
  args$label_text_x <- label_text_xy[, "x"]
  args$label_text_y <- label_text_xy[, "y"]
  
  # REVERT LABEL_TEXT_X & LABEL_TEXT_Y TO ORIGINAL FORMAT
  if (L > 1) {
    args$label_text_x <- LAPPLY(seq_len(L), function(z) {
      args$label_text_x[names(args$label_text_x) == z]
    })
    args$label_text_y <- LAPPLY(seq_len(L), function(z) {
      args$label_text_y[names(args$label_text_y) == z]
    })
  }
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL_COORDS MATRIX
  label_text_xy <- matrix(c(args$label_text_x, args$label_text_y),
                          ncol = 2,
                          byrow = FALSE)
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
  
}

## LEGEND ----------------------------------------------------------------------

## LEGEND ----

#' Create a legend for cyto_plot
#'
#' \code{.cyto_plot_margins} will handle setting the plot margins to make space
#' for the legend.
#'
#' @param x list of flowFrame objects to include in the plot.
#' @param channels name of the channels or markers to be used to construct the
#'   plot.
#' @param legend logical indicating whether a legend should be included for
#'   plots including overlays, set to FALSE by default.
#' @param legend_text vector of labels to use for the legend.
#' @param legend_text_font numeric indicating the font to use for legend text,
#'   set to 2 for bold font by default. See \code{\link[graphics:par]{?par}}
#'   font for details.
#' @param legend_text_size character expansion for legend text, set to 1 by
#'   default.
#' @param legend_text_col colour to use for legend text, set to "black by
#'   default.
#' @param legend_line_col vector of line colours to use for legend.
#' @param legend_box_fill vector of fill colours to use for legend.
#' @param legend_point_col vector of colours to use for points in legend.
#' @param hist_cols vector colours to draw from when selecting density fill
#'   colours if none are supplied to hist_fill.
#' @param hist_fill colour(s) used to fill polygons.
#' @param hist_fill_alpha numeric [0,1] used to control fill transparency,
#'   set to 1 by default to remove transparency.
#' @param hist_line_type line type(s) to use for border(s), set to solid
#'   lines by default.
#' @param hist_line_width line width for border.
#' @param hist_line_col colour(s) for border line, set to "black" by default.
#' @param point_shape point character to use for points, set to "." by default
#'   to maximise plotting speed.
#' @param point_size numeric specifying the degree of character expansion for
#'   points, set to 2 by default.
#' @param point_col_scale vector of colours to use for density gradient.
#' @param point_cols vector colours to draw from when selecting colours for
#'   points if none are supplied to point_col.
#' @param point_col colours to use for points, set to NA by default to blue-red
#'   density colour scale.
#' @param point_col_alpha numeric [0,1] used to control colour transparency, set to
#'   1 by default to remove transparency.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom graphics legend strheight
#' @importFrom grDevices adjustcolor colorRamp rgb
#'
#' @noRd
.cyto_plot_legend <- function(x,
                              channels,
                              legend = "fill",
                              legend_text = NA,
                              legend_text_font = 1,
                              legend_text_size = 1,
                              legend_text_col = "black",
                              legend_line_type = NA,
                              legend_line_width = NA,
                              legend_line_col = NA,
                              legend_box_fill = NA,
                              legend_point_col = NA,
                              hist_cols = NA,
                              hist_fill = NA,
                              hist_fill_alpha = 1,
                              hist_line_type = 1,
                              hist_line_width = 1,
                              hist_line_col = "black",
                              point_shape = ".",
                              point_size = 2,
                              point_col_scale = NA,
                              point_cols = NA,
                              point_col = NA,
                              point_col_alpha = 1,
                              key = "both") {
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # LEGEND_TEXT ----------------------------------------------------------------
  
  # ESTIMATE LEGEND HEIGHT
  lgnd <- paste(legend_text, collapse = " \n ")
  lgnd_height <- strheight(
    lgnd,
    cex = legend_text_size,
    font = legend_text_font
  )
  
  # LEGEND POSITION ------------------------------------------------------------
  
  # PLOT CENTER
  cnt <- .par("usr")[[1]][3] + (.par("usr")[[1]][4] - .par("usr")[[1]][3]) / 2
  
  # LEGEND - HISTOGRAMS - NO KEY
  if (length(channels) == 1) {
    # DEFAULT - FILL
    if (legend == TRUE) {
      legend <- "fill"
    }
    # REVERSE TEXT ORDER FOR LEGEND
    legend_text <- rev(legend_text)
    # LEGEND - LINE
    if (legend == "line") {
      # USE HIST_LINE_COL
      if (.all_na(legend_line_col)) {
        legend_line_col <- hist_line_col
      }
      # USE HIST_LINE_TYPE
      if (.all_na(legend_line_type)) {
        legend_line_type <- hist_line_type
      }
      # USE HIST_LINE_WIDTH
      if (.all_na(legend_line_width)) {
        legend_line_width <- hist_line_width
      }
      # LEGEND
      legend(
        x = 1.07 * .par("usr")[[1]][2],
        y = cnt + 0.52 * lgnd_height,
        legend = legend_text,
        text.font = rev(legend_text_font),
        cex = legend_text_size,
        text.col = rev(legend_text_col),
        col = rev(legend_line_col),
        lty = rev(legend_line_type),
        lwd = rev(legend_line_width),
        xpd = TRUE,
        bty = "n",
        x.intersp = 0.5
      )
    # LEGEND - FILL
    } else if (legend == "fill") {
      # COLOURS
      hist_fill <- .cyto_plot_hist_fill(
        x,
        hist_fill = hist_fill,
        hist_cols = hist_cols,
        hist_fill_alpha = 1
      )
      # USE HIST_FILL
      if (.all_na(legend_box_fill)) {
        legend_box_fill <- hist_fill
      }
      # ALPHA ADJUST COLOURS SUPPLIED TO LEGEND_BOX_FILL
      if (!.all_na(legend_box_fill) &
          !all(hist_fill_alpha == 1)) {
        legend_box_fill <- mapply(
          function(legend_box_fill,
                   hist_fill_alpha) {
            adjustcolor(legend_box_fill, hist_fill_alpha)
          }, legend_box_fill, hist_fill_alpha
        )
      }
      # LEGEND - BOXES - POINT KEY TOO SMALL
      legend(
        x = 1.08 * .par("usr")[[1]][2],
        y = cnt + 0.6 * lgnd_height,
        legend = paste0("   ", legend_text), # HACKY WAY TO ALIGN TEXT
        col = rev(legend_box_fill),
        pch = 15, # BOX - MATCH POINT LEGEND
        pt.cex = 2.4,
        xpd = TRUE,
        bty = "n",
        x.intersp = 0.5, # MOVE DIAGONALLY
        y.intersp = 1.2,
        adj = c(0, 0.49), # XY ADJUSTMENT
        cex = legend_text_size,
        text.col = rev(legend_text_col),
        text.font = rev(legend_text_font)
      )
    }
  # LEGEND - SCATTERPLOTS - KEY
  } else if (length(channels) == 2) {
    # CYTO_PLOT_POINT_COL_SCALE
    point_col_scale <- .cyto_plot_point_col_scale(point_col_scale)
    # DENSITY/MARKER EXPRESSION - MINIMUM COLOUR - KEY
    if(.all_na(point_col[1]) | point_col[1] %in% c(cyto_channels(x[[1]]),
                                                   cyto_markers(x[[1]]))) {
      point_col_scale <- colorRamp(point_col_scale)
      point_col[1] <- rgb(
        point_col_scale(0), 
        maxColorValue = 255
      )
      # ADJUST LEGEND POSITION TO KEY 
      if(key == "none") {
        legend_x <- 1.08 * .par("usr")[[1]][2]
      # KEY - SCALE & TEXT
      } else if(key == "both") {
        legend_x <- 1.15 * .par("usr")[[1]][2]
      # KEY - SCALE ONLY
      } else {
        legend_x <- 1.1 * .par("usr")[[1]][2]
      }
    # NO KEY
    } else {
      legend_x <- 1.08 * .par("usr")[[1]][2]
    }
    # GET REMAINING COLOURS
    point_col <- .cyto_plot_point_col(
      x,
      channels = channels,
      point_col_scale = point_col_scale,
      point_cols = point_cols,
      point_col = point_col,
      point_col_alpha = 1
    )
    # USE POINT_COL
    if (.all_na(legend_point_col)) {
      legend_point_col <- unlist(point_col)
    }
    # ALPHA ADJUST COLOURS SUPPLIED DIRECTLY TO LEGEND_POINT_COL
    if (!.all_na(legend_point_col) &
        !all(point_col_alpha == 1)) {
      legend_point_col <- mapply(function(col, alpha) {
        adjustcolor(col, alpha)
      }, legend_point_col, point_col_alpha)
    }
    # USE SQUARE PCH FOR "."
    point_shape[point_shape == "."] <- 15
    point_shape <- as.numeric(point_shape)
    # LEGEND - BOXES - POINT KEY TOO SMALL
    legend(
      x = legend_x,
      y = cnt + 0.6 * lgnd_height,
      legend = paste0("   ", rev(legend_text)), # HACKY WAY TO ALIGN TEXT
      col = rev(legend_point_col),
      pch = rev(point_shape),
      pt.cex = 2.4,
      xpd = TRUE,
      bty = "n",
      x.intersp = 0.5, # MOVE DIAGONALLY
      y.intersp = 1.2,
      adj = c(0, 0.49), # XY ADJUSTMENT
      cex = legend_text_size,
      text.col = rev(legend_text_col),
      text.font = rev(legend_text_font)
    )
  }
}

## KEY -------------------------------------------------------------------------

#' Add a key for point colour scale to cyto_plot
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @importFrom graphics text rect lines
#' @importFrom grDevices adjustcolor colorRamp
#' 
#' @noRd
.cyto_plot_key <- function(x,
                           channels,
                           xlim = c(NA, NA),
                           ylim = c(NA, NA),
                           point_col_scale = NA,
                           point_col_alpha = 1,
                           point_col = NA,
                           key = "both",
                           key_scale = "fixed",
                           key_text_font = 1,
                           key_text_size = 0.9,
                           key_text_col = "black",
                           key_text_col_alpha = 1,
                           key_title = "",
                           key_title_text_font = 1,
                           key_title_text_size = 1,
                           key_title_text_col = "black",
                           key_title_text_col_alpha = 1,
                           axes_trans = NA) {
  
  # CONSTRUCT KEY
  if(!.all_na(key) | !key == "none") {
    # COMPUTE KEY LOCATION
    usr <- .par("usr")[[1]]
    key_x <- c(
      usr[2] + 0.005 * diff(usr[1:2]),
      usr[2] + 0.035 * diff(usr[1:2])
    )
    key_y <- usr[3:4]
    # KEY AXIS & TITLE
    if(key == "both") {
      # KEY_SCALE REQUIRED
      if(!cyto_class(key_scale, "cyto_plot_key")) {
        key_scale <- .cyto_plot_key_scale(
          list(x),
          channels = channels,
          xlim = xlim,
          ylim = ylim,
          point_col = point_col,
          key_scale = key_scale,
          axes_trans = axes_trans
        )[[1]]
      }
      # KEY AXIS TEXT - BYPASS EMPTY SAMPLES FREE SCALE
      if(!.all_na(key_scale$at)) {
        # KEY AXIS
        lapply(
          seq_along(key_scale$at),
          function(z) {
            # MAJOR TICK
            if(as.character(key_scale$label[z]) != "") {
              # MAJOR TICK
              lines(
                x = c(key_x[2], key_x[2] + 0.6 * diff(key_x)),
                y = c(key_scale$at[z], key_scale$at[z]),
                lwd = 1,
                lty = 1,
                col = "black",
                xpd = TRUE
              )
              # TEXT
              text(
                x = key_x[2] + 0.8 * diff(key_x),
                y = key_scale$at[z],
                pos = 4,
                offset = 0,
                labels = key_scale$label[z],
                font = key_text_font,
                cex = key_text_size,
                col = adjustcolor(key_text_col, key_text_col_alpha),
                xpd = TRUE
              )
              # MINOR TICK
            } else {
              # MINOR TICK
              lines(
                x = c(key_x[2], key_x[2] + 0.3 * diff(key_x)),
                y = c(key_scale$at[z], key_scale$at[z]),
                lwd = 1,
                lty = 1,
                col = "black",
                xpd = TRUE
              )
            }
          }
        )
        # KEY TITLE
        if(!.all_na(key_title)) {
          # DEFAULT TITLE
          if(.empty(key_title)) {
            # COUNTS
            if(.all_na(point_col[1])) {
              key_title <- "count"
            } else {
              # TRY MARKER FIRST (SHORTER)
              key_title <- tryCatch(
                cyto_markers_extract(
                  x,
                  channels = point_col[1]
                ),
                error = function(e){
                  # RESORT TO CHANNEL
                  cyto_channels_extract(
                    x,
                    channels = point_col[1]
                  )
                }
              )
            }
          }
          # TITLE TEXT
          text(
            x = key_x[2] + 0.99 * diff(key_x),
            y = key_y[2] + 0.02 * diff(key_y),
            pos = 3,
            offset = 0,
            labels = key_title,
            font = key_title_text_font,
            cex = key_title_text_size,
            col = adjustcolor(key_title_text_col, key_title_text_col_alpha),
            xpd = TRUE
          )
        }
      }
    }
    # KEY SCALE
    key_boxes <- 60
    # POINT_COL_ALPHA
    key_col_scale <- .cyto_plot_point_col_scale(point_col_scale)
    key_col_ramp <- colorRamp(key_col_scale)
    key_cols <- seq(0, 1, 1/key_boxes)
    key_cols <- key_col_ramp(key_cols)
    key_cols <- rgb(
      key_cols[, 1],
      key_cols[, 2],
      key_cols[, 3],
      maxColorValue = 255
    )
    key_cols <- adjustcolor(key_cols, point_col_alpha[1])
    # KEY BORDER
    rect(
      key_x[1],
      key_y[1],
      key_x[2],
      key_y[2],
      xpd = TRUE,
      lwd = 1
    )
    # KEY BOXES
    key_box_x <- key_x
    key_box_y <- seq(
      usr[3],
      usr[4],
      diff(usr[3:4]) / key_boxes
    )
    lapply(
      seq_len(key_boxes),
      function(z){
        rect(
          key_box_x[1],
          key_box_y[z],
          key_box_x[2],
          key_box_y[z + 1],
          xpd = TRUE,
          col = key_cols[z],
          border = NA
        )
      }
    )
  }

  invisible(NULL)
  
}

#' Prepare scales for cyto_plot keys
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @importFrom flowCore exprs
#' 
#' @noRd
.cyto_plot_key_scale <- function(x,
                                 channels,
                                 xlim = c(NA, NA),
                                 ylim = c(NA, NA),
                                 point_col = NA,
                                 key_scale = "fixed",
                                 axes_trans = NA) {

  # X AXES LIMITS REQUIRED - CROP GRID
  if(any(is.na(xlim))) {
    xlim <- .par("usr")[[1]][1:2]
  } else {
    # BUFFER XLIM - 4% MATCH PLOT LIMITS
    xlim <- c(
      min(xlim) - (diff(xlim) - diff(xlim) / 1.04) / 2,
      max(xlim) + (diff(xlim) - diff(xlim) / 1.04) / 2
    )
  }
  
  # Y AXES LIMITS REQUIRED - CROP GRID
  if(any(is.na(ylim))) {
    ylim <- .par("usr")[[1]][3:4]
  } else {
    # BUFFER YLIM - 4% MATCH PLOT LIMITS
    ylim <- c(
      min(ylim) - (diff(ylim) - diff(ylim) / 1.04) / 2,
      max(ylim) + (diff(ylim) - diff(ylim) / 1.04) / 2
    )
  }
  
  # STORE RANGES FOR SCALE
  key_range <- NULL
  
  # KEY RANGES MANUALLY SUPPLIED ON LINEAR SCALE - NO NA SAFEGUARD
  if(cyto_class(key_scale, "list", TRUE) | is.numeric(key_scale)) {
    # FIXED SCALE - VECTOR
    if(is.numeric(key_scale)) {
      key_range <- structure(
        rep(
          list(key_scale),
          length(x)
        ),
        names = names(x)
      )
    # SCALE - LIST
    } else {
      key_range <- structure(
        rep(
          key_scale,
          length.out = length(x)
        ),
        names = names(x)
      )
    }
  }
  
  # COLOUR SCALE - COUNTS
  if(is.na(point_col[1])) {
    # SCALE MANUALLY SUPPLIED
    if(is.null(key_range)) {
      # COMPUTE BINNED COUNTS
      key_range <- structure(
        lapply(
          seq_along(x),
          function(z) {
            # COMPUTE BINNED COUNTS
            range(
              cyto_apply(
                x[[z]][[1]],
                "cyto_stat_bkde2d",
                input = "matrix",
                channels = channels,
                limits = list(xlim, ylim),
                smooth = FALSE,
                copy = FALSE,
                simplify = FALSE
              )[[1]]$counts,
              na.rm = TRUE
            )
          }
        ),
        names = names(x)
      )
      # FIXED SCALE
      if(key_scale == "fixed") {
        key_range  <- structure(
          rep(
            list(range(unlist(key_range), na.rm = TRUE)),
            length(x)
          ),
          names = names(key_range)
        )
      }
    }
    # PREPARE KEYS
    key <- structure(
      lapply(
        key_range,
        function(z){
          # NA BKDE - TOO FEW EVENTS
          if(.all_na(z)) {
             k <-list(
              "range" = NA,
              "label" = NA,
              "at" = NA
            )
          # PREPARE AXIS TICKS & LABELS
          } else {
            k <- .cyto_plot_axes_text(
              channels = NA,
              axes_range = list(z),
              rescale = ylim,
              format = TRUE
            )[[1]]
            k$range <- z
          }
          class(k) <- "cyto_plot_key"
          return(k)
        }
      )
    )
  # COLOUR SCALE - MARKER EXPRESSION
  } else if(point_col[1] %in% c(cyto_markers(x[[1]][[1]]),
                                cyto_channels(x[[1]][[1]]))) {
    # CONVERT POINT_COL TO CHANNEL
    point_col[1] <- cyto_channels_extract(
      x[[1]][[1]], 
      channels = point_col[1]
    )
    
    # KEY RANGES REQUIRED
    if(is.null(key_range)) {
      # CALIBRATION SETTINGS
      cyto_cal <- .cyto_plot_calibrate_recall()
      # FIXED KEY SCALE
      if(key_scale == "fixed") {
        # TRY CALIBARTION SETTINGS
        if(point_col[1] %in% colnames(cyto_cal)) {
          key_range <- structure(
            rep(
              list(cyto_cal[, point_col[1]]),
              length(x)
            ),
            names = names(x)
          )
        }
      }
      # FREE KEY SCALE - DATA RANGES
      if(is.null(key_range)) {
        # COMPUTE DATA RANGES
        key_range <- structure(
          lapply(
            seq_along(x),
            function(z) {
              # COMPUTE DATA RANGE
              cyto_apply(
                x[[z]][[1]],
                "cyto_stat_range",
                input = "matrix",
                channels = point_col[1],
                copy = FALSE
              )[, 1]
            }
          ),
          names = names(x)
        )
        # EXCLUDE INFINITE LIMITS
        key_range <- lapply(key_range, function(v){
          if(!any(is.finite(v))) {
            v[!is.finite(v)] <- NA
          }
          return(v)
        })
        # FIXED KEY SCALE - DATA RANGE
        if(key_scale == "fixed") {
          key_range  <- structure(
            rep(
              list(range(unlist(key_range), na.rm = TRUE)),
              length(x)
            ),
            names = names(key_range)
          )
        }
      }
    # TRANSFORM MANUALLY SUPPLIED SCALE
    } else {
        # TRANSFORM SCALES
        if(point_col[1] %in% names(axes_trans)) {
          key_range <- structure(
            lapply(
              seq_along(key_range),
              function(z) {
                .cyto_transform(key_range[[z]],
                                trans = axes_trans,
                                channel = point_col[1],
                                inverse = FALSE)
              }
            ),
            names = names(key_range)
          )
        }
    }
    # PREPARE KEYS
    key <- structure(
      lapply(
        key_range,
        function(z){
          # INFINITE RANGE
          if(.all_na(z) | any(!is.finite(z))) {
            k <- list(
              "range" = NA,
              "label" = NA,
              "at" = NA
            )
            # PREPARE AXIS TICKS & LABELS
          } else {
            k <- .cyto_plot_axes_text(
              channels = point_col[1],
              axes_range = list(z),
              axes_trans = axes_trans,
              rescale = ylim,
              format = TRUE
            )[[1]]
            k$range <- z
          }
          class(k) <- "cyto_plot_key"
          return(k)
        }
      )
    )
  # NO COLOUR SCALE REQUIRED
  } else {
    k <- list(
      "range" =  NA,
      "label" = NA,
      "at" = NA
    )
    class(k) <- "cyto_plot_key"
    key <- rep(
      list(k),
      length(x)
    )
    names(key) <- names(x)
  }
  
  # RETURN KEYS
  return(key)
  
}

## TITLE -----------------------------------------------------------------------

#' Title for cyto_plot
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_title <- function(x,
                             channels,
                             overlay = NA,
                             title = "") {
  
  # x can be a list
  if (class(x) == "list") {
    if (length(x) > 1) {
      overlay <- x[2:length(x)]
      x <- x[[1]]
    }
  }
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Update arguments
  .args_update(args)
  
  # 1D density distributions
  if (length(channels) == 1) {
    # missing/empty replace with valid title
    if (.empty(title)) {
      # stacked/overlays lack a title
      if (.all_na(overlay)) {
        title <- cyto_names(x)
      } else {
        title <- NA
      }
    # NA will remove title in cyto_plot_empty
    } else if (.all_na(title)) {
      title <- NA
    }
  # 2D scatterplots
  } else if (length(channels) == 2) {
    # missing title replaced with sample name
    if (.empty(title)) {
      title <- cyto_names(x)
      # NA will remove title in cyto_plot_empty
    } else if (.all_na(title)) {
      title <- NA
    }
  }
  
  # COMBINED/ALL EVENTS
  title <- LAPPLY(title, function(z){
    if(!.all_na(z)) {
      z <- gsub("^all", "Combined Events", z)
      z <- gsub("root", "All Events", z)
    }
    return(z)
  })
  
  return(title)
}

## HEADER ----------------------------------------------------------------------

#' Add header to plot
#' 
#' @importFrom graphics mtext
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_header <- function(header,
                              header_text_font = 2,
                              header_text_size = 1,
                              header_text_col = "black") {
  
  # HEADER - OMA MUST BE PRESET - c(0,0,3,0)
  mtext(header, 
        outer = TRUE, 
        font = header_text_font,
        cex = header_text_size,
        col = header_text_col)
  
}

## POINTS ----------------------------------------------------------------------

## POINT COLOUR SCALE ----

#' Get density gradient colours for cyto_plot
#'
#' @param point_col_scale vector of ordered colours to use for point density
#'   colour scale.
#'
#' @return a list of colorRampPalette functions to be used in densCols.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_point_col_scale <- function(point_col_scale = NA) {
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Use default colour scale
  if (.all_na(args[["point_col_scale"]])) {
    args[["point_col_scale"]] <- .cyto_plot_colour_palette("point_col_scale")
  }
  
  return(args[["point_col_scale"]])
}

## POINT COLOUR ----

#' Get point colours for cyto_plot
#'
#' @param x list of flowFrames.
#' @param channels used to construct the plot.
#' @param point_col_scale vector of colours to use for density gradient.
#' @param point_cols vector colours to select from when choosing a colour for
#'   each layer in x.
#' @param point_col vector of length x indicating colours to use for each layer.
#'   If NA set to default density gradient.
#' @param point_col_alpha transparency to use for point colours.
#' @param point_col_smooth logical indicating whether the 2D binned counts
#'   should be smoothed using kernel density estimates prior to selecting
#'   colours from \code{point_col_scale}, set to TRUE by default. Setting
#'   \code{point_col_smooth} to FALSE will significantly improve plotting speed
#'   on less powerful machines but produce more granular plots.
#' @param ... additional arguments inherited from \code{cyto_plot_point()} or
#'   \code{cyto_plot()} including \code{bkde2d} and \code{key_scale}
#'
#' @importFrom grDevices colorRampPalette adjustcolor colorRamp rgb
#' @importFrom flowCore exprs
#' @importFrom graphics par
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_point_col <- function(x,
                                 channels,
                                 point_col_scale,
                                 point_cols,
                                 point_col,
                                 point_col_alpha = 1,
                                 point_col_smooth = TRUE,
                                 ...) {
  
  # PULL DOWN ARGUMENTS --------------------------------------------------------
  
  # ARGUMENTS - CARRY ADDITIONAL ARGUMENTS FROM CYTO_PLOT()
  args <- .args_list(...)
  
  # INHERIT CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # POINT_COL_SCALE
  if(.all_na(args$point_col_scale)) {
    args$point_col_scale <- .cyto_plot_colour_palette(
      type = "point_col_scale"
    )
  }
  
  # CREATE COLORRAMP USING POINT_COL_SCALE
  if(!cyto_class(args$point_col_scale, "function")) {
    args$point_col_scale <- colorRamp(
      args$point_col_scale
    )
  }
  
  # USE DEFUALT COLOURS FOR POINTS
  if(.all_na(args$point_cols)) {
    args$point_cols <- .cyto_plot_colour_palette(
      type = "point_cols"
    )
  }
  
  # CREATE COLORRAMPPALETTE USING POINT_COLS
  if (!cyto_class(args$point_cols, "function")) {
    args$point_cols <- colorRampPalette(
      args$point_cols
    )
  }
  
  # REPEAT POINT COLOUR ARGUMENTS
  lapply(
    c("point_col", "point_col_alpha"),
    function(z) {
      args[[z]] <<- rep(args[[z]], length.out = length(args$x))
    }
  )
  
  # POINT_COL -> LIST
  if(!cyto_class(args$point_col, "list")) {
    args$point_col <- lapply(
      seq_along(args$point_col),
      function(z) {
        args$point_col[z]
      }
    )
  }
  
  # BASE LAYER - NA/MARKER/COLOUR
  if(length(args$point_col[[1]]) == 1) {
    # DENSITY GRADIENT
    if(.all_na(args$point_col[[1]])) {
      # NUMBER OF EVENTS
      N <- nrow(cyto_exprs(x[[1]][[1]], drop = FALSE))
      # INSUFFICIENT EVENTS
      if(N < 2) {
        # USE MINIMUM COLOUR FOR POINTS
        args$point_col[[1]] <- rgb(
          args$point_col_scale(rep(0, N)),
          maxColorValue = 255
        )
        # COMPUTE BKDE - MAP TO ROWS
      } else {
        # BKDE COUNTS REQUIRED
        if(!"bkde2d" %in% names(args)) {
          args$bkde2d <- cyto_apply(
            args$x[[1]],
            "cyto_stat_bkde2d",
            input = "matrix",
            channels = args$channels,
            limits = list(.par("usr")[[1]][1:2],
                          .par("usr")[[1]][3:4]),
            smooth = args$point_col_smooth,
            copy = FALSE,
            simplify = FALSE
          )[[1]]
        }
        # SCALE RANGE
        if(!"key_scale" %in% names(args)) {
          # USE RANGE OF COUNTS FOR SCALING
          args$key_scale <- list(
            range = range(args$bkde2d$counts),
            label = NA,
            at = NA
          )
        }
        # RE-SCALE BINNED COUNTS - RANGE [0,1]
        args$bkde2d$counts <- cyto_stat_rescale(
          args$bkde2d$counts,
          scale = matrix(
            args$key_scale$range,
            nrow = 2,
            ncol = ncol(args$bkde2d$counts)
          )
        )
        # RE-SCALE BKDE
        if(args$point_col_smooth) {
          # STORE RE-SCALED BKDE IN COUNTS SLOT
          args$bkde2d$counts <- min(args$bkde2d$counts) +
            ((args$bkde2d$bkde - min(args$bkde2d$bkde))/
               diff(range(args$bkde2d$bkde))) * diff(range(args$bkde2d$counts))
        }
        # MAP BKDE TO ROWS & ASSIGN COLOURS
        args$point_col[[1]] <- cyto_apply(
          args$x[[1]],
          channels = args$channels,
          input = "matrix", 
          copy = FALSE,
          simplify = FALSE,
          bkde = args$bkde,
          point_col_scale = args$point_col_scale,
          FUN = function(
            z,
            bkde,
            point_col_scale) {
            # COMPUTE BREAKS FOR CUT - SEE DENSCOLS()
            mkBreaks <-  function(u){
              u - diff(range(u))/(length(u)-1)/2
            }
            # COMPUTE BINS
            b <- do.call(
              "cbind",
              lapply(
                seq_len(2),
                function(w){
                  cut(
                    z[, w],
                    breaks = mkBreaks(
                      bkde$bins[[w]]
                    ),
                    labels = FALSE
                  )
                }
              )
            )
            # MAP DATA TO (SMOOTHED) COUNTS
            b <- bkde$counts[b]
            b[is.na(b)] <- min(bkde$counts, na.rm = TRUE) # EVENTS OUTSIDE GRID
            # CONVERT TO COLOURS
            return(
              rgb(
                point_col_scale(b),
                maxColorValue = 255
              )
            )
          }
        )[[1]]
      }
      # BASE LAYER - MARKER EXPRESSION - DATA SORTED IN CYTO_PLOT_POINT()
    } else if(all(args$point_col[[1]] %in% c(cyto_channels(args$x[[1]]),
                                             cyto_markers(args$x[[1]])))) {
      # CONVERT POINT_COL TO CHANNEL
      args$point_col[[1]] <- cyto_channels_extract(
        args$x[[1]],
        channels = args$point_col[[1]]
      )
      # INSUFFICIENT EVENTS
      N <- nrow(cyto_exprs(x[[1]][[1]]))
      if(N < 2) {
        # USE MINIMUM COLOUR FOR POINTS
        args$point_col[[1]] <- rgb(
          args$point_col_scale(rep(0, N)),
          maxColorValue = 255
        )
      } else {
        # KEY_SCALE
        if(!"key_scale" %in% names(args)) {
          # CALIBRATION SETTINGS
          cal <- .cyto_plot_calibrate_recall()
          # USED STORED SETTINGS
          if(args$point_col[[1]] %in% colnames(cal)) {
            args$key_scale <- list(
              range = cal[, args$point_col[[1]]],
              label = NA,
              at = NA
            )
            # USE DATA RANGE
          } else {
            args$key_scale <- list(
              range = cyto_apply(
                args$x[[1]],
                "cyto_stat_range",
                channels = args$point_col[[1]],
                copy = FALSE,
                input = "matrix"
              )[, 1],
              label = NA,
              at = NA
            )
          }
        }
        # RE-SCALE DATA TO RANGE & ASSIGN COLOURS
        args$point_col[[1]] <- cyto_apply(
          args$x[[1]],
          channels = args$point_col[[1]],
          input = "column",
          copy = FALSE,
          scale = args$key_scale$range,
          FUN = function(z,
                         scale) {
            rgb(
              args$point_col_scale(
                cyto_stat_rescale(
                  z, 
                  scale = scale
                )
              ),
              maxColorValue = 255
            )
          }
        )
      }
    }
  }

  # REMAINING LAYERS - SELECT FROM POINT_COLS
  if (any(LAPPLY(args$point_col, ".all_na"))) {
    # NUMBER OF REQUIRED COLOURS
    n <- length(
      args$point_col[
        LAPPLY(args$point_col, ".all_na")
      ]
    )
    # GET MISSING COLOURS FROM POINT_COLS
    clrs <- args$point_cols(n)
    # UPDATE COLOURS IN POINT_COL
    args$point_col[
      LAPPLY(args$point_col, ".all_na")
    ] <- clrs
  }
  
  # POINT COLOUR ALPHA ADJUSTMENT
  lapply(seq_along(args$point_col), function(z) {
    args$point_col[[z]] <<- adjustcolor(
      args$point_col[[z]], 
      args$point_col_alpha[z]
    )
  })
  
  # RETURN LIST OF COLOURS
  return(args$point_col)
}

## HISTOGRAMS ------------------------------------------------------------------

## DENSITY ----

#' Prepare histograms for cyto_plot
#' @param x list of cytosets
#' @importFrom graphics par
#' @noRd
.cyto_plot_hist <- function(x,
                            channels = NULL,
                            hist_stat = "count",
                            hist_smooth = 1,
                            hist_bins = 256,
                            hist_stack = 0,
                            xlim = c(NA, NA),
                            ...) {
  
  # COMPUTE BANDWIDTH BY BINNING XLIM (256 BINS - FLOWJO)
  # INSTRUMENT RANGE IS INCONSISTENT - FLOWWORKSPACE ISSUE #348
  
  # NOTE: NON-FINITE FROM ERROR - DUE TO AXES LIMITS
  
  # KERNEL DENSITY - LIST OF DENSITY LISTS PER SAMPLE
  d <- structure(
    lapply(x, function(cs){
      # HISTOGRAMS - INSUFFICIENT EVENTS WARNING
      suppressWarnings(
        cyto_apply(
          cs, 
          "cyto_stat_density",
          input = "matrix",
          channels = channels,
          copy = FALSE,
          stat = hist_stat,
          bins = hist_bins,
          smooth = hist_smooth,
          limits = matrix(
            xlim,
            ncol = 1,
            dimnames = list(
              c("min", "max"),
              channels
            )
          ),
          simplify = FALSE
        )[[1]][[1]]
      )
    }), names = names(x)
  )
  
  # STACKING - RANGE STORED IN DENSITY OBJECT -  CANNOT ROUND DENSITY
  hist_heights <- LAPPLY(
    d, 
    function(D){
      if(.all_na(D)){
        return(NA)
      }else{
        D$range[2]
      }
    }
  )
  
  # NO EVENTS TO PLOT IN ANY LAYERS
  if(.all_na(hist_heights)) {
    # DEFAULT COUNT SCALE
    if(hist_stat == "count") {
      hist_heights <- structure(
        rep(1000, length(hist_heights)),
        names = names(hist_heights)
      )
    # DEFAULT PERCENT SCALE
    } else if(hist_stat == "percent") {
      hist_heights <- structure(
        rep(100, length(hist_heights)),
        names = names(hist_heights)
      )
    # DEFAULT DENSITY SCALE
    } else {
      hist_heights <- structure(
        rep(1, length(hist_heights)),
        names = names(hist_heights)
      )
    }
  # SOME LAYERS CONTAIN NO EVENTS
  } else if(any(is.na(hist_heights))) {
    # USE MEAN HEIGHT FOR EMPTY LAYERS
    hist_heights[is.na(hist_heights)] <- mean(hist_heights, na.rm = TRUE)
  }
  hist_stack <- max(hist_heights, na.rm = TRUE) * hist_stack
  hist_levels <- c(0, seq_along(x)) * hist_stack
  d <- structure(
    lapply(
      seq_along(d), 
      function(z){
        D <- d[[z]]
        if(z > 1 & !.all_na(D)){
          D$y <- D$y + hist_levels[z]
        }
        return(D)
      }
    ),
    names = names(d)
  )
  
  # STORE NEW RANGES IN RANGE SLOT
  d <- structure(
    lapply(
      seq_along(d), 
      function(z){
        # STORE RNAGE IN LIST
        if(.all_na(d[[z]])) {
          return(
            list(
              NA,
              range = c(hist_levels[z],
                        hist_levels[z] + hist_heights[z])
            )
          )
        # APPEND RANGE TO DENSITY OBJECTS
        } else {
          d[[z]]$range <- c(hist_levels[z],
                            hist_levels[z] + hist_heights[z])
          return(d[[z]])
        }
      }
    ), 
    names = names(d)
  )
  
  # STACKED HISTOGRAMS
  return(d)
  
}

## HISTOGRAM FILL ----

#' Get density fill colours for cyto_plot
#'
#' @param x list of flowFrame or density objects.
#' @param density_fill vector of colours to use for each layer.
#' @param density_cols vector of colls to use to select density_fill colours.
#'
#' @importFrom grDevices adjustcolor colorRampPalette
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_hist_fill <- function(x,
                                 hist_fill = NA,
                                 hist_cols = NA,
                                 hist_fill_alpha = 1) {
  
  # INHERIT CYTO_PLOT_THEME ----------------------------------------------------
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Update arguments
  .args_update(args)
  
  # GENERAL --------------------------------------------------------------------
  
  # Expected number of colours
  SMP <- length(x)
  
  # DENSITY_FILL ---------------------------------------------------------------
  
  # No hist_cols supplied
  if (.all_na(hist_cols)) {
    hist_cols <- .cyto_plot_colour_palette(type = "hist_cols")
  }
  
  # Make colorRampPalette
  if (class(hist_cols) != "function") {
    cols <- colorRampPalette(hist_cols)
  } else {
    cols <- hist_cols
  }
  
  # No colours supplied to hist_fill either
  if (.all_na(hist_fill)) {
    
    # Pull out a single colour per layer
    hist_fill <- cols(SMP)
    
    # Colours supplied manually to hist_fill
  } else {
    
    # Too few colours supplied - pull others from cols
    if (length(hist_fill) < SMP) {
      hist_fill <- c(
        hist_fill,
        cols(SMP - length(hist_fill))
      )
      
      # Too many colours supplied
    } else if (length(hist_fill) > SMP) {
      hist_fill <- hist_fill[seq_len(SMP)]
    }
  }
  
  # Adjust colors by hist_fill_alpha
  hist_fill <- mapply(function(hist_fill, hist_fill_alpha) {
    if (hist_fill_alpha != 1) {
      adjustcolor(hist_fill, hist_fill_alpha)
    } else {
      hist_fill
    }
  }, hist_fill, hist_fill_alpha, USE.NAMES = FALSE)
  
  return(hist_fill)
}

## THEME -----------------------------------------------------------------------

## INHERIT THEME ----

#' Inherit cyto_plot_theme arguments
#'
#' @param x list of named cyto_plot arguments.
#'
#' @return updated list of named arguments if cyto_plot_theme has been set.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_theme_inherit <- function(x) {
  
  # extract cyto_plot_theme arguments
  args <- getOption("cyto_plot_theme")
  
  if (!is.null(args)) {
    for(y in names(args)) {
      x[[y]] <- args[[y]]
    }
  }
  
  return(x)
}

## PALETTES --------------------------------------------------------------------

#' cyto_plot colour palette
#'
#' @param type indicates whether to return the "point_cols", "point_col_scale"
#'   or "hist_cols" colour palette.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_colour_palette <- function(type = "point_cols") {
  # POINT COLOUR PALETTE
  if (type == "point_cols") {
    pal <- c(
      "grey25",
      "bisque4",
      "brown1",
      "red",
      "darkred",
      "chocolate",
      "orange",
      "yellow",
      "yellowgreen",
      "green",
      "limegreen",
      "turquoise",
      "aquamarine",
      "cyan",
      "cornflowerblue",
      "blue",
      "blueviolet",
      "purple4",
      "purple",
      "magenta",
      "deeppink"
    )
    
    # POINT COLOUR SCALE
  } else if (type == "point_col_scale") {
    
    # # NEW VIRIDIS PALETTE
    # viridis_pal <- viridis(12, option = "D")
    # plasma_pal <- viridis(12, option = "C")
    # viridis_pal <- viridis_pal[-12]
    # plasma_pal <- plasma_pal[-c(1:4, 12)]
    # custom_pal <- c(viridis_pal, "#E1DD37FF", rev(plasma_pal))
    # scales::show_col(custom_pal)
    
    pal <- c(
      "#440154FF",
      "#482173FF",
      "#433E85FF",
      "#38598CFF",
      "#2D708EFF",
      "#25858EFF",
      "#1E9B8AFF",
      "#2BB07FFF",
      "#51C56AFF",
      "#85D54AFF",
      "#C2DF23FF",
      "#E1DD37FF",
      "#FCD225FF",
      "#FDAD32FF",
      "#F58C46FF",
      "#E76F5AFF",
      "#D5546EFF",
      "#C03A83FF",
      "#A62098FF"
    )
    
    # DENSITY COLOUR PALETTE
  } else if (type == "hist_cols") {
    pal <- c(
      "grey50",
      "bisque4",
      "brown1",
      "red",
      "darkred",
      "chocolate",
      "orange",
      "yellow",
      "yellowgreen",
      "green",
      "limegreen",
      "turquoise",
      "aquamarine",
      "cyan",
      "cornflowerblue",
      "blue",
      "blueviolet",
      "purple4",
      "purple",
      "magenta",
      "deeppink"
    )
  }
  
  return(pal)
}
