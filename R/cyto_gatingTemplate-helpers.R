## CYTO_GATINGTEMPLATE_EXTRACT -------------------------------------------------

#' Extract a gatingTemplate from a GatingHierrachy or GatingSet
#'
#' Simply a modified wrapper around
#' \code{\link[openCyto:gh_generate_template]{gh_generate_template()}} to
#' extract a gatingTemplate from either a GatingHierarchy or GatingSet.
#'
#' @param x an object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select passed to \code{cyto_select()} to select the GatingHierarchy
#'   from which the gatingTemplate should be extracted, set to 1 by default to
#'   select the first GatingHierarchy.
#' @param sort logical indicating whether the entries in the gatingTemplate
#'   should be sorted to match the order of the nodes in the GatingHierarchy or
#'   GatingSet, set to TRUE by default.
#' @param bool logical indicating whether attempts should be made to supply
#'   \code{dims} for boolean filters, set to FALSE by default. Basically, we
#'   check whether the boolean gate was defined based on gates constructed on
#'   the same \code{parent} in the same channels.
#' @param drop logical indicating whether entries with empty \code{dims} should
#'   be excluded from the extracted gatingTemplate, set to FALSE by default.
#' @param data.table logical indicating whether the extracted gatingTemplate
#'   should be returned as a \code{data.table}, set to FALSE by default.
#' @param ... not in use.
#'
#' @return gatingTemplate as either a \code{data.frame} or \code{data.table}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom openCyto gh_generate_template
#' @importFrom data.table data.table
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#' gs <- cyto_compensate(gs)
#' gs <- cyto_transform(gs)
#' gs <- cyto_gatingTemplate_apply(gs)
#'
#' # Extract gatingTemplate
#' cyto_gatingTemplate_extract(gs)
#'
#' @export
cyto_gatingTemplate_extract <- function(x,
                                        select = 1,
                                        sort = TRUE,
                                        bool = FALSE,
                                        drop = FALSE,
                                        data.table = FALSE,
                                        ...) {
  
  # SELECT GATINGHIERARCHY
  gh <- cyto_select(
    x, 
    select
  )
  
  # EXTRACT GATINGTEMPLATE
  gt <- gh_generate_template(gh)
  
  # CONVERT EMPTY DIMS -> NA
  gt$dims[!nzchar(gt$dims)] <- NA
  
  # TRY SETTING BOOLEAN FILTER DIMS - CYTO_PLOT_GATING_SCHEME()
  if(any(is.na(gt$dims)) & bool) {
    # ENTRIES WITH EMPTY DIMS
    ind <- which(is.na(gt$dims))
    # CHECK ENTRIES
    lapply(
      ind,
      function(z) {
        # EXTRACT GATE
        gate <- cyto_gate_extract(
          gh,
          alias = cyto_nodes_convert(
            gh,
            nodes = gt$alias[z],
            anchor = gt$parent[z],
            path = "auto"
          ),
          bool = FALSE
        )[[1]][[1]]
        # PARSE DIMS FROM VALID BOOLEANFILTERS
        if(cyto_class(gate, "booleanFilter")) {
          # EXTRACT POPULATIONS
          pops <- unlist(
            strsplit(
              gate@deparse,
              "[^[:alnum:][:space:]]"
            )
          )
          pops <- pops[nzchar(pops)]
          pops <- tryCatch(
            cyto_nodes_convert(
              gh,
              nodes = pops,
              anchor = gt$parent[z],
              path = "auto"
            ),
            error = function(e) {
              return(NULL)
            }
          )
          # BOOLENAFILTER REFERENCING VALID POPULATIONS
          if(!is.null(pops)) {
            # SET GATING METHOD
            gt$gating_method[z] <<- "boolGate"
            # NEGATED BOOLEANFILTERS ONLY (FOR NOW)
            ops <- unlist(
              strsplit(
                gate@deparse,
                "[[:alnum:][:space:]]"
              )
            )
            ops <- ops[nzchar(ops)]
            # POPULATION DIMS
            pop_dims <-  unique(
              gt$dims[gt$alias %in% pops & gt$parent %in% gt$parent[z]]
            )
            # NEGATED BOOLEANFILTER BASED ON POPS IN SAME DIMS
            if(length(pop_dims) == 1 & 
               ops[1] %in% "!") {
              if(length(ops) > 1) {
                if(all(ops[-1] %in% c("!&", "&!"))) {
                  gt$dims[z] <<- pop_dims
                }
              } else {
                gt$dims[z] <<- pop_dims
              }
            }
          }
        }
      }
    )
  }
  
  # DROP EMPTY DIMS ENTRIES
  if(drop) {
    gt <- gt[which(!is.na(gt$dims)), , drop = FALSE]
  }
  
  # SORT ENTRIES
  if(sort) {
    # AUTO PATHS FOR ENTRIES
    pops <- LAPPLY(
      seq_len(nrow(gt)),
      function(z) {
        cyto_nodes_convert(
          gh,
          nodes = gt$alias[z],
          anchor = gt$parent[z],
          path = "auto"
        )
      }
    )
    # SORT GATINGTEMPLATE
    gt <- gt[
      order(
        match(
          pops, 
          cyto_nodes(gh, path = "auto")[-1]
        )
      ), , drop = FALSE
    ]
  }
  
  # DATA.TABLE
  if(data.table) {
    data.table(gt)
  }
  
  # GATINGTEMPLATE
  return(gt)
  
}

## GATINGTEMPLATE GENERATE -----------------------------------------------------

#' Generate a CytoExploreR gatingTemplate for a GatingHierarchy or GatingSet
#'
#' @param x an object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select sample selection criteria passed to \code{cyto_select} to
#'   identify the samples to included in the generated gatingTemplate, set to
#'   NULL by default to use all samples in \code{x}.
#' @param merge_by vector of experiment variable names to indicate how samples
#'   should be merged prior to creating the gatingTemplate entries, set to
#'   \code{"name"} by default to store a unique gate for every sample in
#'   \code{x}. If \code{merge_by} is specified, a consensus gate will be created
#'   for each sample group and stored in the gatingTemplate.
#' @param nodes vector of population names to include in the gatingTemplate, set
#'   to NULL by default to include all nodes.
#' @param gatingTemplate logical indicating whether the generated gatingTemplate
#'   should be returned as a \code{gatingTemplate} object instead of a
#'   \code{data.frame}, set to TRUE by default.
#' @param save_as name of a CSV file to which the generated gatingTemplate
#'   should be written.
#' @param active logical to indicate whether the generated gatingTemplate should
#'   be set to the active gatingTemplate, set to FALSE by default.
#' @param inverse placeholder for future support for storing gates in the
#'   gatingTemplate on the linear scale.
#' @param ... additional arguments passed to \code{gatingTemplate}.
#'
#' @return write generated gatingTemplate to designated CSV file and return the
#'   gatingTemplate object.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom data.table as.data.table
#' @importFrom grDevices chull
#' @importFrom flowCore rectangleGate polygonGate quadGate filters
#' @importFrom flowWorkspace gh_pop_is_negated gh_pop_get_gate gh_pop_get_data
#' @importFrom openCyto CytoExploreR_.argDeparser gatingTemplate
#'
#' @export
cyto_gatingTemplate_generate <- function(x,
                                         select = NULL,
                                         merge_by = "name",
                                         nodes = NULL,
                                         gatingTemplate = TRUE,
                                         save_as = NULL,
                                         active = FALSE,
                                         inverse = FALSE,
                                         ...){
  
  # TODO: INVERSE IS PLACEHOLDER FOR WHEN GATE TRANSFORM SUPPORTED IN FLOWCORE
  if(inverse) {
    warning(
      "'inverse' is not currently supported."
    )
  }
  
  # GATINGSET | GATINGHIERARCHY ONLY
  if(!cyto_class(x, "GatingSet")) {
    stop(
      "cyto_gatingTemplate_generate() requires either a GatingHierarchy or ",
      "GatingSet!"
    )
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # TRANSFORMERS
  trans <- cyto_transformers_extract(x)
  
  # FULL NODE PATHS
  nodes_full <- cyto_nodes(
    x,
    path = "full"
  )
  
  # ALL NODES
  if(is.null(nodes)) {
    nodes <- nodes_full
  # SUBSET NODES
  } else {
    nodes <- cyto_nodes_convert(
      x, 
      nodes = nodes,
      path = "full"
    )
  }
  
  # GROUPBY - SPLIT SAMPLES INTO GROUPS
  pd_groups <- cyto_groups(
    x, 
    group_by = merge_by,
    details = TRUE,
    select = select
  )
  
  # COMPUTE INSTRUMENT LIMITS - REPLACE INFINITE COORDS
  rng <- apply(
    do.call(
      "rbind",
      lapply(
        seq_along(x),
        function(id) {
          cf <- gh_pop_get_data(
            x[[id]]
          )
          # USE BOTH DATA & INSTRUMENT LIMITS
          rbind(
            range(cf, type = "instrument"),
            range(cf, type = "data")
          )
        }
      )
    ),
    2,
    "range"
  )
  
  # INITIALIZE GATINGTEMPLATE VARIABLES - R CMD CHECK NOTE
  parent <- NULL
  alias <- NULL
  pop <- NULL
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  
  # GENERATE BARE BONES GATINGTEMPLATE
  gt <- cyto_gatingTemplate_extract(x)
  
  gt$dims[is.na(gt$dims)] <- "NA"
  
  # SPLIT GATINGTEMPLATE BY PARENT & DIMS
  gt_chunks <- split(
    gt,
    gt[, c("parent", "dims")],
    sep = "|",
    drop = TRUE
  )
  
  # ORDER GT_CHUNKS BY NODES - CORRECT ORDER IN GATINGTEMPLATE
  idx <- sapply(
    gt_chunks,
    function(gt_chunk) {
      parent <- gt_chunk$parent
      parent[parent == "root"] <- ""
      alias <- unlist(strsplit(gt_chunk$alias, ","))
      min(
        match(
          paste0(
            parent,
            "/",
            alias
          ),
          nodes
        )
      )
    }
  )
  
  # DROP UNWANTED NODES
  gt_chunks <- gt_chunks[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  
  # REORDER GATINGTEMPLATE CHUNKS
  gt_chunks <- gt_chunks[order(idx)]
  
  # CONSTRUCT GATINGTEMPLATE
  gt_new <- do.call(
    "rbind",
    structure(
      lapply(
        gt_chunks,
        function(gt_chunk) {
          # PARENT
          parent <- unique(gt_chunk$parent)
          # CHANNELS
          channels <- strsplit(unique(gt_chunk$dims), ",")[[1]]
          # ALIAS
          alias <- gt_chunk$alias
          # REQUIRE: LIST OF FILTERS PER POPULATION IN EACH SAMPLE
          group_gate_list <- structure(
            lapply(
              pd_groups,
              function(pd_group) {
                # SAMPLE INDICES
                idx <- match(pd_group$name, pd$name)
                # EXTRACT GATES PER SAMPLE IN GROUP
                sample_gates <- lapply(
                  idx,
                  function(id) {
                    # EXTRACT GATE OBJECTS
                    gate_list <- structure(
                      lapply(
                        alias,
                        function(pop) {
                          gate <- gh_pop_get_gate(
                            x[[id]],
                            paste0(
                              if(parent == "root") {
                                ""
                              } else {
                                parent
                              },
                              "/",
                              pop
                            )
                          )
                          # INHERIT NEGATED FLAG
                          attributes(gate)$negated <- gh_pop_is_negated(
                            x[[id]],
                            paste0(
                              if(parent == "root") {
                                ""
                              } else {
                                parent
                              },
                              "/",
                              pop
                            )
                          )
                          return(gate)
                        }
                      ),
                      names = alias
                    )
                    # FORMAT QUADRANT GATES
                    quad_idx <- which(
                      sapply(
                        gate_list,
                        function(gate){
                          if(inherits(gate, "rectangleGate")) {
                            any(grepl("quad", names(attributes(gate))))
                          } else {
                            FALSE
                          }
                        }
                      )
                    )
                    # QUADRANTS LOCATED
                    if(length(quad_idx) > 0) {
                      # EXTRACT RECTANGLEGATES
                      quads <- gate_list[quad_idx]
                      # REMOVE QUADRANTS FROM GATE_LIST
                      gate_list <- gate_list[-quad_idx]
                      # QUADRANT CO-ORDINATES
                      coords <- do.call(
                        "rbind",
                        lapply(
                          quads,
                          function(quad) {
                            structure(
                              c(quad@min, quad@max),
                              names = parameters(quad)
                            )
                          }
                        )
                      )
                      # QUADRANT CENTER
                      coords <- unlist(
                        sapply(
                          channels,
                          function(z) {
                            unique(coords[, z][is.finite(coords[, z])])
                          }
                        )
                      )
                      names(coords) <- channels
                      # CONSTRUCT QUADRANTGATE
                      quads <- list(
                        quadGate(
                          filterId = paste0(
                            names(attributes(quads[[1]])[["quadrants"]]),
                            collapse = "|"
                          ),
                          .gate = coords
                        )
                      )
                      attributes(quads[[1]])$negated <- FALSE
                      gate_list <- c(quads, gate_list)
                    }
                    names(gate_list) <- sapply(
                      gate_list,
                      function(gate) {
                        gate@filterId
                      }
                    )
                    # # APPLY INVERSE TRANSFORMERS
                    # if(inverse & length(trans) > 0) {
                    #   # EXTRACT INVERSE TRANSFORMLIST
                    #   inv_trans <- transformList(
                    #     names(trans),
                    #     lapply(
                    #       trans,
                    #       `[[`,
                    #       "inverse"
                    #     )
                    #   )
                    #   # apply inverse transformers to gates
                    #   gate_list <- structure(
                    #     lapply(
                    #       gate_list,
                    #       function(gate) {
                    #         if(any(parameters(gate) %in% names(trans))) {
                    #           # excess transformers may be supplied
                    #           gate <- transform(
                    #             gate,
                    #             inv_trans
                    #           )
                    #         }
                    #         return(gate)
                    #       }
                    #     ),
                    #     names = names(gate_list)
                    #   )
                    # }
                    return(gate_list)
                  }
                )
                # CONSTRUCT GATINGTEMPLATE ENTRIES
                structure(
                  lapply(
                    names(sample_gates[[1]]),
                    function(gate) {
                      # EXTRACT UNIQUE GATES ACROSS SAMPLES IN GROUP
                      gates <- unique(
                        lapply(
                          sample_gates,
                          `[[`,
                          gate
                        )
                      )
                      # SAMPLES IN GROUP HAVE DIFFERENT GATES
                      if(length(gates) > 1) {
                        # MULTIRANGEGATE - USE FIRST GATE
                        if(cyto_class(gates[[1]],c("booleanFilter", "multiRangeGate"))) {
                          return(gates[[1]])
                        # QUADRANTGATE
                        } else if(cyto_class(gates[[1]], "quadGate")) {
                          coords <- colMeans(
                            do.call(
                              "rbind",
                              lapply(
                                gates,
                                function(pop) {
                                  pop@boundary
                                }
                              )
                            )
                          )
                          return(
                            quadGate(
                              .gate = coords,
                              filterId = gate
                            )
                          )
                        # 1D RECTANGLEGATE
                        } else if(inherits(gates[[1]], "rectangleGate") &
                                  length(parameters(gates[[1]])) == 1) {
                          # COMBINE GATE COORDINATES
                          coords <- do.call(
                            "rbind",
                            lapply(
                              gates,
                              function(pop) {
                                rbind(
                                  pop@min,
                                  pop@max
                                )
                              }
                            )
                          )
                        # COERCE COORDS & CHULL
                        } else {
                          coords <- do.call(
                            "rbind",
                            lapply(
                              gates,
                              function(pop) {
                                as(pop, "polygonGate")@boundaries
                              }
                            )
                          )
                        }
                        # REPLACE INFINITE COORDS WITH INSTRUMENT LIMITS
                        lapply(
                          colnames(coords),
                          function(chan) {
                            coords[, chan][
                              is.infinite(coords[, chan]) & coords[, chan] < 0
                            ] <<- min(rng[, chan])
                            coords[, chan][
                              is.infinite(coords[, chan]) & coords[, chan] > 0
                            ] <<- max(rng[, chan])
                          }
                        )
                        # 1D RANGE GATE
                        if(ncol(coords) == 1) {
                          # RANGE ACROSS GATES
                          chull <- rectangleGate(
                            .gate = apply(
                              coords,
                              2,
                              "range"
                            ),
                            filterId = gate
                          )
                        # 2D POLYGON
                        } else {
                          # CHULL POLYGONGATE
                          chull <- polygonGate(
                            .gate = coords[chull(coords), ],
                            filterId = gate
                          )
                        }
                        # INHERIT NEGATED FLAG
                        attributes(chull)$negated <- 
                          attributes(gates[[1]])$negated
                        # COSENSUS GATE
                        return(chull)
                      } else {
                        return(gates[[1]])
                      }
                    }
                  ),
                  names = names(sample_gates[[1]])
                )
              }
            )
          )
          # CONASTRUCT GATINGTEMPLATE ENTRIES
          do.call(
            "rbind",
            lapply(
              names(group_gate_list[[1]]),
              function(alias) {
                # GATE FLAGS
                bool <- FALSE
                quad <- FALSE
                negated <- FALSE
                # GATES PER GROUP
                gates <- structure(
                  lapply(
                    group_gate_list,
                    function(gate_list) {
                      gate <- gate_list[[alias]]
                      # NEGATED FLAG
                      negated <<- attributes(gate)$negated
                      # BOOLEAN FLAG
                      if(inherits(gate, "booleanFilter")) {
                        bool <<- TRUE
                      }
                      # DON'T WRAP QUADGATES IN FILTERS
                      if(inherits(gate, "quadGate")) {
                        quad <<- TRUE
                        gate_list[[alias]]
                      } else {
                        filters(gate_list[alias])
                      }
                    }
                  ),
                  names = names(group_gate_list)
                )
                # GATINGTEMPLATE ENTRY
                gt_entry <- as.data.table(gt_chunk[1, , drop = FALSE])
                # ALIAS
                gate <- gsub("\\|", ",", alias)
                # POP
                if(negated) {
                  gt_entry[, pop := "-"]
                } else {
                  gt_entry[, pop := "+"]
                }
                # ALIAS
                gt_entry[, alias := gate]
                # PREPROCESSING_METHOD
                gt_entry[, preprocessing_method := "pp_cyto_gate_draw"]
                # GATING_METHOD
                gt_entry[, gating_method := "cyto_gate_draw"]
                # COLLAPSEDATAFORGATING REQUIRED - OR NEED GATE FOR ALL SAMPLES
                gt_entry[, collapseDataForGating := TRUE]
                # GROUPBY
                if(all(is.na(merge_by))) {
                  group_by <- "NA"
                } else {
                  group_by <- paste0(
                    merge_by,
                    collapse = ":"
                  )
                }
                gt_entry[, groupBy := group_by]
                # BOOLEANFILTER
                if(bool) {
                  gt_entry[, groupBy := NA]
                  gt_entry[, preprocessing_method := NA]
                  gt_entry[, gating_method := "boolGate"]
                  gt_entry[, gating_args := gates[[1]][[1]]@deparse]
                  gt_entry[, collapseDataForGating := NA]
                # OTHER GATES PREPARED ABOVE
                } else {
                  if(quad) {
                    gt_entry[, pop := "*"]
                  }
                  gt_entry[, gating_args := CytoExploreR_.argDeparser(
                    list(
                      gate = gates,
                      openCyto.minEvents  = -1
                    )
                  )]
                }
                return(gt_entry)
              }
            )
          )
        }
      ),
      names = names(gt_chunks)
    )
  )
  
  # BACKWARDS COMPATIBLE GATINGTEMPLATE
  if(is.character(gatingTemplate)) {
    save_as <- gatingTemplate
    gatingTemplate <- TRUE
  }
  
  # WRITE GATINGTEMPLATE TO CSV FILE
  if(!is.null(save_as)) {
    write_to_csv(
      gt_new,
      save_as
    )
    if(active) {
      cyto_gatingTemplate_active(
        save_as
      )
    }
  }
  
  # CONSTRUCT GATINGTEMPLATE OBJECT
  if(gatingTemplate) {
    gt_new <- tryCatch(
      gatingTemplate(
        gt_new,
        ...
      ),
      error = function(e) {
        return(gt_new)
      }
    )
  }
  
  # return constructed gatingTemplate
  return(gt_new)
  
}

## GATINGTEMPLATE ACTIVE -------------------------------------------------------

#' Retrieve or set active gatingTemplate
#'
#' @param x name of the gatingTemplate to set as the active gatingTemplate. If
#'   not supplied, the name of the current active gatingTemplate will be
#'   returned.
#' @param reset logical indicating whether the active gatingTemplate setting
#'   should be reset so that no active gatingTemplate is set, set to FALSE by
#'   default.
#' @param ask logical indicating whether the user should be asked to supply a
#'   name for the active gatingTemplate if one has not been set, set to FALSE
#'   by default.
#'
#' @return retrieve or set current active gatingTemplate.
#'
#' @seealso \code{\link{cyto_gatingTemplate_create}}
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_gatingTemplate_active <- function(x = NULL, 
                                       reset = FALSE,
                                       ask = FALSE) {
  if(reset) {
    message(
      paste0(
        "Resetting the active gatingTemplate..."
      )
    )
    cyto_option("CytoExploreR_gatingTemplate", NULL)
  } else {
    if(is.null(x)) {
      gt <- cyto_option("CytoExploreR_gatingTemplate")
      # ACTIVE GATINGTEMPLATE NOT SET
      if(is.null(gt) & ask == TRUE) {
        if(interactive()){
          x <- cyto_enquire(
            "Please supply a name for the active gatingTemplate:"
          )
        } else {
          stop(
            "Please supply the name of a gatingTemplate CSV file."
          )
        }
      # ACTIVE GATINGTEMPLATE SET
      } else {
        return(gt)
      }
    }
    # SET NEW ACTIVE GATINGTEMPLATE
    x <- file_ext_append(x, ".csv")
    message(
      paste0(
        "Setting ", x, " as the active gatingTemplate..."
      )
    )
    cyto_option("CytoExploreR_gatingTemplate", x)
    # ACTIVE GATINGTEMPLATE
    invisible(x)
  }
}

#' @export
cyto_gatingTemplate_select <- function(...){
  .Defunct("cyto_gatingTemplate_active")
}

## GATINGTEMPLATE CREATE -------------------------------------------------------

#' Create an empty gatingTemplate csv file
#'
#' @param gatingTemplate name of the gatingTemplate csv file to create.
#' @param active logical indicating whether the created gatingTemplate should be
#'   set as the active gatingTemplate, set to FALSE by default.
#' @param data.table logical passed to \code{cyto_gatingTemplate_read()} to
#'   control whether the empty gatingTemplate is returned as a \code{data.table}
#'   or \code{data.frame}, set to FALSE by default.
#'
#' @return an empty gatingTemplate in the form of a \code{data.table} or
#'   \code{data.frame}.
#'
#' @importFrom utils write.csv
#' @importFrom stats setNames
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_gatingTemplate_create <- function(gatingTemplate = NULL,
                                       active = FALSE,
                                       data.table = FALSE) {
  
  # ACTIVE GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  
  # CREATE GATINGTEMPLATE
  cols <- c(
    "alias",
    "pop",
    "parent",
    "dims",
    "gating_method",
    "gating_args",
    "collapseDataForGating",
    "groupBy",
    "preprocessing_method",
    "preprocessing_args"
  )
  # EMPTY GATINGTEMPLATE
  gt <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
  # FILE EXTENSION
  gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
  # GATINGTEMPLATE EXISTS
  if(file_exists(gatingTemplate)) {
    if(!cyto_enquire(
      paste0(
        gatingTemplate, 
        " already exists. Do you want to replace it? (Y/N)"
      ),
      options = c("T", "Y")
    )) {
      invisible(NULL)
    }
  }
  # MESSAGE
  message(paste("Creating", gatingTemplate, "..."))
  # WRITE GATINGTEMPLATE
  write_to_csv(gt, gatingTemplate)
  # SET ACTIVE GATINGTEMPLATE
  if(active){
    cyto_gatingTemplate_active(gatingTemplate)
  }
  # RETURN EMPTY DATA.TABLE
  invisible(
    cyto_gatingTemplate_read(
      gatingTemplate,
      data.table = data.table,
      parse = FALSE
    )
  )
}

## GATINGTEMPLATE_EDIT ---------------------------------------------------------

#' Interactively Edit a gatingTemplate
#'
#' \code{cyto_gatingTemplate_edit} provides an interactive interface to editing
#' the openCyto gatingTemplate. This function is intended to aid in adding
#' boolean and reference gates to the gatingTemplate. Users should NOT modify
#' existing entries in the gatingTemplate.
#'
#' @param x object of class \code{GatingSet} which has been gated with the
#'   supplied gatingTemplate.
#' @param gatingTemplate name of the gatingTemplate csv file to be edited.
#' @param ... not in use.
#'
#' @return update gatingTemplate csv file and update the GatingSet accordingly.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom flowWorkspace recompute
#' @importFrom DataEditR data_edit
#'
#' @examples
#' \dontrun{
#' library(CytoExploreR)
#'
#' # gs is a GatingSet object
#' cyto_gatingTemplate_edit(gs, "gatingTemplate.csv")
#' }
#'
#' @export
cyto_gatingTemplate_edit <- function(x, 
                                     gatingTemplate = NULL,
                                     ...) {
  
  # x of wrong class
  if (!cyto_class(x, "GatingSet")) {
    stop("'x' should be a GatingSet object.")
  }
  
  # Assign x to gs
  gs <- x
  
  # Missing gatingTemplate - check global ooption
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # gatingTemplate still missing
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate csv file to edit.")
  }
  
  # Read in gatingTemplate
  gt <- read_from_csv(gatingTemplate)
  
  # Edit gatingTemplate
  message("Do not modify existing gatingTemplate entries!")
  message("Add new rows to add boolean or reference gates to the GatingSet.")
  
  # Edit gatingTemplate
  gt <- data_edit(
    gt, 
    title = "gatingTemplate Editor",
    logo = CytoExploreR_logo(),
    col_edit = FALSE,
    col_names = colnames(gt),
    quiet = TRUE,
    hide = TRUE,
    viewer = "pane",
    ...
  )
  
  # BE NICE - REMOVE WHITESPACE FROM DIMS (COMMON MISTAKE)
  gt$dims <- gsub(" *, *", ",", gt$dims)
  
  # WRITE TO FILE
  write_to_csv(gt, gatingTemplate)
  
  # Read in updated template to gatingTemplate object
  gt <- gatingTemplate(gatingTemplate)
  
  # Re-apply template to GatingSet
  suppressMessages(gt_gating(gt, gs))
  
  # Recompute statistics
  suppressMessages(recompute(gs))
  
  # Update GatingSet globally
  assign(deparse(substitute(x)), gs, envir = globalenv())
  
  # Return gatingTemplate object
  invisible(return(gt))
}

## GATINGTEMPLATE_APPLY ---------------------------------------------------------

#' Apply gates saved in gatingTemplate to a GatingHierarchy or GatingSet
#'
#' A convenient wrapper for
#' \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}} and
#' \code{\link[openCyto:gt_gating]{gt_gating}} to apply a gatingTemplate csv
#' file to a GatingSet.
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param gatingTemplate name of the gatingTemplate csv file which contains the
#'   gates to be applied to \code{x}. This can alos be an object of class
#'   \code{gatingTemplate} if a \code{GatingSet} object is supplied.
#' @param active logical indicating whether the applied gatingTemplate should be
#'   set as the current active gatingTemplate set to TRUE by default.
#' @param overwrite logical indicating whether existing gates should be removed
#'   prior to applying the gatingTemplate, only required when running
#'   CytoExploreR non-interactively.
#' @param ... additional arguments passed to
#'   \code{\link[openCyto:gt_gating]{gt_gating()}} including \code{start} and
#'   \code{stop.at} to only apply particular gates from the gatingTemplate.
#'
#' @return NULL and update the GatingHierarchy/ GatingSet in the global
#'   environment with the gates in the gatingTemplate.
#'
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
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
#' # Apply gatingTemplate in working directory to GatingSet
#' cyto_gatingTemplate_apply(gs, "gatingTemplate.csv")
#' }
#'
#' @export
cyto_gatingTemplate_apply <- function(x,
                                      gatingTemplate = NULL,
                                      active = TRUE,
                                      overwrite = NULL,
                                      ...) {
  
  # GATINGHIERARCHY/GATINGSET REQUIRED
  if (missing(x)) {
    stop("Supply a GatingHierarchy or GatingSet to 'x'.")
  } else {
    if (!cyto_class(x, "GatingSet")) {
      stop("'x' must be an object of class GatingHierarchy or GatingSet.")
    }
  }
  
  # GATINGTEMPLATE SUPPLIED
  if (!is.null(gatingTemplate)) {
    # GATINGTEMPLATE
    if (cyto_class(gatingTemplate, "gatingTemplate")) {
      gt <- gatingTemplate
    # GATINGTEMPLATE FILE
    } else {
      # GATINGTEMPLATE FILE
      if(is.character(gatingTemplate)) {
        # READ GATINGTEMPLATE FROM FILE
        gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
        file_exists(gatingTemplate, error = TRUE)
      # COERCE TO DATA.TABLE
      } else {
        gatingTemplate <- as(gatingTemplate, "data.table")
      }
      # GATINGTEMPLATE FILE -> GATINGTEMPLATE
      gt <- suppressMessages(gatingTemplate(gatingTemplate))
    }
  # NO GATINGTEMPLATE
  } else if (is.null(gatingTemplate)) {
    # ACTIVE GATINGTEMPLATE
    if (is.null(gatingTemplate)) {
      gatingTemplate <- cyto_gatingTemplate_active()
    }
    # GATINGTEMPLATE MISSING
    if (is.null(gatingTemplate)) {
      stop("Supply the name of the gatingTemplate csv file to apply.")
    }
    # READ GATINGTEMPLATE
    gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
    file_exists(gatingTemplate, error = TRUE)
    # GATINGTEMPLATE FILE -> GATINGTEMPLATE
    gt <- suppressMessages(gatingTemplate(gatingTemplate))
  }
  
  # GATINGSET CONTAINS GATES
  if(length(cyto_nodes(x)) > 1) {
    # EMPTY OVERWRITE & INTERACTIVE
    if(!is.logical(overwrite)) {
      # INTERACTIVE MODE
      if(cyto_option("CytoExploreR_interactive")) {
        overwrite <- cyto_enquire(
          "Overwrite existing gates in the GatingSet? (Y/N): ",
          options = c("T", "Y")
        )
      # NON-INTERACTIVE MODE
      } else {
        overwrite <- FALSE
      }
    }
    # OVERWRITE GATES
    if(overwrite) {
      # MESSAGE
      message(
        "Removing existing gates from the GatingSet..."
      )
      # REMOVE GATES
      lapply(
        cyto_nodes_kin(
          x,
          nodes = "root",
          type = "children",
          path = "full"
        ),
        function(z) {
          gs_pop_remove(
            x,
            z
          )
        }
      )
    # WARNING IN CASE
    } else {
      warning(
        paste0(
          "Applying a gatingTemplate to a GatingSet that already contains ",
          "gates is not recommended - see overwrite argument."
        )
      )
    }
  }
  
  # MESSAGE - APPLY GATINGTEMPLATE FILE
  message(
    paste(
      "Applying", 
      if(is.character(gatingTemplate)) {
        gatingTemplate
      } else {
        "gatingTemplate"
      },
      "to the GatingSet..."
    )
  )
  
  # APPLY GATINGTEMPLATE
  suppressWarnings(
    gt_gating(
      gt,
      x,
      ...
    )
  )
  
  # ACTIVE GATINGTEMPLATE
  if(active) {
    if(is.character(gatingTemplate)) {
      cyto_gatingTemplate_active(gatingTemplate)
    }
  }
  
  # GATINGHIERARCHY/GATINGSET
  return(x)
  
}

## GATINGTEMPLATE UPDATE ------------------------------------------------------

#' Convert CytoRSuite gatingTemplate to be compatible with CytoExoloreR
#'
#' @param gatingTemplate name of the gatingTemplate csv file to convert.
#' @param save_as name for the updated gatingTemplate csv file, set to "Updated
#'   gatingTemplate.csv" by default.
#'
#' @return converted gatingTemplate saved to .csv file.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_gatingTemplate_update <- function(gatingTemplate = NULL,
                                       save_as = "Updated-gatingTemplate.csv") {
  
  # READ IN GATINGTEMPLATE
  gt <- read_from_csv(gatingTemplate)
  
  # CONVERT OLD GATING METHOD
  gt[gt$gating_method == "gate_draw", "gating_method"] <- "cyto_gate_draw" 
  
  # CONVERT OLD PREPROCESSING  METHOD
  gt[gt$preprocessing_method == "pp_gate_draw",
     "preprocessing_method"] <- "pp_cyto_gate_draw" 
  
  # SAVE UPDATED GATINGTEMPLATE
  write_to_csv(gt, save_as)
  
}

## GATINGTEMPLATE CHECK --------------------------------------------------------

#' Check gatingTemplate for Existing Entry
#'
#' @param parent name of the parent population.
#' @param alias name of the population of interest.
#' @param gatingTemplate csv file name of the gatingTemplate.
#' @param data.table logical indicating whether the checked gatingTemplate
#'   should be returned as a \code{data.table}, set to TRUE by default.
#'
#' @return stops the gating process if an entry already exists in the
#'   gatingTemplate for the supplied alias.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_gatingTemplate_check <- function(parent,
                                       alias,
                                       gatingTemplate = NULL,
                                       data.table = TRUE) {
  
  # BYPASS CHECKS - GATINGTEMPLATE DOES NOT EXIST
  gt <- tryCatch(
    suppressWarnings(
      cyto_gatingTemplate_read(
        gatingTemplate,
        parse = TRUE,
        data.table = TRUE
      )
    ),
    error = function(e) {
      NULL
    }
  )
  
  # MATCHING PARENT & ALIAS
  if (!is.null(gt)) {
    if(nrow(gt) > 0 ){
      gt_chunk <- gt[
        LAPPLY(
          seq_len(nrow(gt)),
          function(z) {
            all(parent %in% unlist(strsplit(gt$parent[z], ",")))
          }
        ),
        , drop = FALSE
      ]
      # CHECK ALIAS
      if(nrow(gt_chunk) > 0) {
        lapply(
          seq_len(nrow(gt_chunk)),
          function(z) {
            pops <- unlist(strsplit(gt_chunk$alias[z], ","))
            if(!is.null(alias)) {
              if(any(alias %in% pops)) {
                # POPULATION(S) EXIST
                message(
                  paste0(
                    "The following population(s) already exist in this ",
                    ifelse(
                      is.character(gatingTemplate),
                      gatingTemplate,
                      "gatingTemplate"
                    ), ": \n",
                    paste0(
                      alias[alias %in% pops],
                      collapse = "\n"
                    )
                  )
                )
                # ERROR
                stop(
                  paste0(
                    "Supply another gatingTemplate or edit existing gate(s) ",
                    "using cyto_gate_edit()."
                  )
                )
              }
            }
          }
        )
      }
    }
    # RETURN UNPARSED GATINGTEMPLATE
    gt <- cyto_gatingTemplate_read(
      gatingTemplate,
      parse = FALSE,
      data.table = data.table
    )
  }
  
  # GATINGTEMPLATE
  return(gt)
  
}

## CYTO_GATINGTEMPLATE_READ ----------------------------------------------------

#' Read in a gatingTemplate object or gatingTemplate csv file
#'
#' @param gatingTemplate gatingTemplate object or the name of a gatingTemplate
#'   csv file to read in.
#' @param data.table logical indicating whether the gatingTemplate should be
#'   read into a \code{data.table}, set to FALSE by default to use
#'   \code{data.frame}.
#' @param parse logical indicating whether attempts should be made to parse
#'   \code{alias = "*"} into a valid set of comma separated aliases, set to
#'   FALSE by default. Required to appropriately parse gating entries created by
#'   \code{cyto_gate_clust()}.
#' @param ... additional arguments passed to
#'   \code{\link[data.table:fread]{fread}}.
#'
#' @importFrom data.table as.data.table 
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # gatingTemplate
#' gt <- cyto_gatingTemplate_read(Activation_gatingTemplate)
#'
#' # write gatingTemplate
#' cyto_gatingTemplate_write(gt, "gatingTemplate.csv")
#'
#' # gatingTemplate csv file
#' gt <- cyto_gatingTemplate_read("gatingTemplate.csv")
#' }
#'
#' @export
cyto_gatingTemplate_read <- function(gatingTemplate = NULL,
                                     data.table = FALSE, 
                                     parse = FALSE,
                                     ...){
  
  # NO GATINGTEMPLATE
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
    if(is.null(gatingTemplate)){
      stop("Supply the name of the gatingTemplate csv file to read in.")
    }
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if(cyto_class(gatingTemplate, "character")){
    gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
    file_exists(gatingTemplate, error = TRUE)
  }
  
  # GATINGTEMPLATE
  if(cyto_class(gatingTemplate, "gatingTemplate")){
    gt <- as.data.table(gatingTemplate)
    # DROP EMPTY POP
    gt <- gt[!gt[["pop"]] %in% c("", NA), , drop = FALSE]
    # PARSE OUT POP = * WITH SINGLE ALIAS
    ind <- which(
      gt[["pop"]] %in% "*" & !grepl("/|\\*", gt[["alias"]])
    )
    # ROWS TO PARSE
    if(length(ind) > 0) {
      gt_parse <- gt[ind, ]
      gt_parse <- cbind(gt[ind, ], "index" = ind)
      gt_parse <- split(
        gt_parse,
        as.list(
          gt_parse
        )[
          colnames(gt_parse)[
            -match(
              c("alias", "pop", "index"),
              colnames(gt_parse)
            )
          ]
        ],
        drop = TRUE
      )
      lapply(
        gt_parse,
        function(z) {
          # COLLAPSE ALIAS
          gt_alias <- paste0(
            z[["alias"]],
            collapse = ","
          )
          # REPLACE ALIAS
          gt[min(z[["index"]]), "alias"] <<- gt_alias
          # KEEP FIRST ROW
          ind <- z[["index"]][!z[["index"]] == min(z[["index"]])]
          # REMOVE EXCESS ROWS
          gt <<- gt[-ind, ]
        }
      )
    }
  # DATA.TABLE OR DATA.FRAME
  } else if( cyto_class(gatingTemplate, c("data.table", "data.frame"))) {
    gt <- gatingTemplate
    if(!cyto_class(gt, "data.table")){
      gt <- as.data.table(gt)
    } 
  # FILE NAME
  } else if(cyto_class(gatingTemplate, "character")) {
    gt <- read_from_csv(
      gatingTemplate,
      data.table = TRUE,
      ...
    )
  # UNSUPPORTED GATINGTEMPLATE
  } else {
    stop(
      paste0(
        "'gatingTemplate' must be either a name of a CSV file, a ",
        "gatingTremplate or a data.table."
      )
    )
  }
  
  # PARSE GATINGTEMPLATE
  if(parse) {
    gt <- cyto_gatingTemplate_parse(gt)
  }
  
  # FORMAT GATINGTEMPLATE - DATA.FRAME REQUIRED
  if(!data.table) {
    gt <- as.data.frame(gt)
  }
  
  # RETURN GATINGTEMPLATE
  return(gt)
  
}

## CYTO_GATINGTEMPALTE_WRITE ---------------------------------------------------

#' Write gatingTemplate to csv file
#'
#' @param gatingTemplate data.frame, data.table, gatingTemplate object or the
#'   name of a gatingTemplate csv file.
#' @param save_as name of the csv file to write the gatingtemplate to.
#' @param ... additional arguments passed to
#'   \code{\link[data.table:fwrite]{fwrite}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # gatingTemplate
#' gt <- cyto_gatingTemplate_write(Activation_gatingTemplate)
#'
#' # write gatingTemplate
#' cyto_gatingTemplate_write(gt, "gatingTemplate.csv")
#' }
#'
#' @export
cyto_gatingTemplate_write <- function(gatingTemplate = NULL,
                                      save_as = NULL,
                                      ...){
  
  # SAVE_AS
  if(is.null(save_as)){
    stop(
      "Supply the name of the csv file to 'save_as'."
    )
  }else{
    save_as <- file_ext_append(
      save_as,
      ".csv"
    )
  }
  
  # READ GATINGTEMPLATE
  if(cyto_class(gatingTemplate, c("character", "gatingTemplate"))){
    gatingTemplate <- cyto_gatingTemplate_read(
      gatingTemplate,
      data.table = TRUE
    )
  }
  
  # HANDLE NAs
  gatingTemplate$groupBy[
    LAPPLY(
      gatingTemplate$groupBy,
      ".empty"
    )
  ] <- "NA"
  gatingTemplate$preprocessing_args[
    LAPPLY(
      gatingTemplate$preprocessing_args,
      ".empty"
    )
  ] <- "NA"
  
  # WRITE GATINGTEMPLATE
  write_to_csv(
    gatingTemplate,
    file = save_as, 
    na = "NA",
    ...
  )
  
}

## CYTO_GATINGTEMPLATE_PARSE ---------------------------------------------------

#' Convert parsed gatingTemplate into a format compatible with CytoExploreR
#'
#' @param gatingTemplate gatingTemplate object or the name of a gatingTemplate
#'   csv file to read in and parse.
#' @param data.table logical indicating whether the parsed gatingTemplate should
#'   be returned as a \code{data.table}, set to TRUE by default.
#' @param ... not in use
#'
#' @return parsed gatingTemplate as either a \code{data.table} or
#'   \code{data.frame}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#' library(CytoExploreRData)
#' # Parse Activation gatingTemplate
#' cyto_gatingTemplate_parse(Activation_gatingTemplate)
#' }
#'
#' @export
cyto_gatingTemplate_parse <- function(gatingTemplate,
                                      data.table = TRUE,
                                      ...) {
  
  # TODO: UPDATE .PREPROCESS_ROW IMPORT WHEN AVAILABLE IN OPENCYTO
  
  # READ IN GATINGTEMPLATE
  gt <- cyto_gatingTemplate_read(
    gatingTemplate,
    data.table = TRUE
  )
  
  # PARSE GATINGTEMPLATE ENTRIES
  lapply(
    seq_len(nrow(gt)),
    function(z) {
      # OPENCYTO PARSING|EXPANSION
      gt_chunk <- openCyto:::.preprocess_row(
        gt[z, ]
      )
      # CYTO_GATE_CLUST - PULL ALIAS FROM GATING_ARGS
      if(all(gt_chunk$alias %in% "*" & 
             gt_chunk$gating_method %in% "cyto_gate_clust")) {
        # EXTRACT ALIAS FROM GATING ARGUMENTS
        gating_args <- eval(
          parse(
            text = paste0(
              "list(",
              gt_chunk$gating_args,
              ")"
            )
          )
        )
        # REPLACE * ALIAS ENTRY
        gt_chunk$alias <- rep(
          paste0(
            gating_args$alias,
            collapse = ","
          ), 
          nrow(gt_chunk)
        )
      }
      # PREPARE ALIAS - UPDATE GATINGTEMPLATE
      gt$alias[z] <<- paste0(
        unique(
          unlist(
            strsplit(
              gt_chunk$alias,
              ","
            )
          )
        ), 
        collapse = ","
      )
    }
  )
  
  # RETURN PARSED GATINGTEMPLATE
  if(!data.table) {
    gt <- as.data.frame(gt)
  }
  return(gt)
  
}
