## CYTO_GATE_APPLY -------------------------------------------------------------

#' Apply a list of gates to a cytoframe or cytoset
#'
#' \code{cyto_gate_apply()} is a helper function used within CytoExploreR to
#' handle all gating operations. In particular it ensures appropriate handling
#' of \code{quadGate} and \code{rectangleGates} associated with
#' \code{quadGates}.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}} or
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}.
#' @param gate a filters, rectangleGate, polygonGate, ellipsoidGate, quadGate or
#'   list of these objects.
#' @param negate logical flag indicating whether events outside of the supplied
#'   gates should be included as a separate population, set to FALSE by default.
#' @param skip logical indicating whether gating should be skipped for gates
#'   when\code{negate = TRUE}, set to FALSE by default. Setting this option to
#'   TRUE will return a list containing only the negated population.
#'
#' @return a list of cytoframes or cytosets per gate.
#'
#' @importFrom flowCore Subset split quadGate
#' @importFrom flowWorkspace flowSet_to_cytoset
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_gate_apply <- function(x,
                            gate = NA,
                            negate = FALSE,
                            skip = FALSE) {
  
  # GATE - LIST
  if(!cyto_class(gate, "list", TRUE)) {
    gate <- unlist(list(gate))
  # EXTRACT FILTERS
  } else if(cyto_class(gate, "filters", TRUE)) {
    gate <- unlist(gate)
  # EXTRACT ANY FILTERS IN GATE LIST
  } else {
    gate <- unlist(gate)
  }
  
  # RECTANGLEGATES TO QUADGATE
  quad_order <- NULL
  if (length(gate) == 4 & all(LAPPLY(gate, function(z) {
    cyto_class(z, "rectangleGate") & any(grepl("quad", names(attributes(z))))
  }))) {
    # CHANNELS
    chans <- as.character(parameters(gate[[1]]))
    quad_order <- LAPPLY(gate, function(z) {
      z@filterId
    })
    gate <- list(.cyto_gate_quad_convert(gate, channels = chans))
  }
  
  # NEGATE - APPEND FILTER INCLUDING ALL GATES
  if(negate == TRUE) {
    if(length(gate) == 1) {
      negate_filter <- !gate[[1]]
    } else {
      negate_filter <- !do.call("|", unname(unlist(gate)))
    }
    if(!skip) {
      gate <- c(gate, negate_filter)
    } else {
      gate <- list("negate" = negate_filter)
    }
  }
  
  # GATED POPULATIONS
  alias <- c()
  pops <- lapply(seq_along(gate), function(z){
    # NO GATE
    if(.all_na(gate[[z]])) {
      return(list(x))
    }
    # QUADGATES - MULTIPLE POPULATIONS
    if(cyto_class(gate[[z]], "quadGate")) {
      # ORDER
      quads <- unlist(strsplit(gate[[z]]@filterId, "\\|"))
      # ALIAS
      alias <<- c(alias, quads)
      # FIX ORDER FROM ABOVE
      if (!is.null(quad_order)) {
        pop <- tryCatch(
          # FIX ORDER
          split(x, gate[[z]])[c(2, 1, 3, 4)][match(quad_order, quads)],
          error = function(e){
            return(
              structure(
                list(x, x, x, x),
                names = quads
              )
            )
          }
        )
      } else {
        pop <- tryCatch(
          split(x, gate[[z]])[c(2, 1, 3, 4)], # FIX ORDER
          error = function(e) {
            return(
              structure(
                list(x, x, x, x),
                names = quads
              )
            )
          }
        )
      }
      # *** CYTOSET CONVERSION ***
      pop <- structure(
        lapply(pop, function(q){
          cyto_convert(q)
        }),
        names = names(pop)
      )
      # SINGLE POPULATIONS
    } else {
      # RECTANGLE BELONGS TO QUADGATE
      if (cyto_class(gate[[z]], "rectangleGate") &
          any(grepl("quad", names(attributes(gate[[z]]))))) {
        q <- names(attributes(gate[[z]])[["quadrants"]])
        coords <- .cyto_gate_coords(gate[z],
                                    channels = as.character(
                                      parameters(gate[[z]]))
        )
        chans <- colnames(coords)
        coords <- lapply(colnames(coords), function(y) {
          unique(coords[, y][is.finite(coords[, y])])
        })
        names(coords) <- chans
        qg <- quadGate(
          filterId = paste(q, collapse = "|"),
          .gate = coords
        )
        p <- tryCatch(
          split(x, qg)[c(2, 1, 3, 4)], # FIX ORDER
          error = function(e) {
            return(
              structure(
                list(x, x, x, x),
                names = q
              )
            )
          }
        )
        names(p) <- q
        # ALIAS
        alias <<- c(alias, q)
        # *** CYTOSET CONVERSION ***
        p <- structure(
          lapply(p, function(b){
            cyto_convert(b)
          }),
          names = names(p)
        )
        pop <- p[[match(gate[[z]]@filterId, names(p))]]
      } else {
        # alias
        alias <<- c(alias, gate[[z]]@filterId)
        pop <- tryCatch(
          list(Subset(x, gate[[z]])),
          error = function(e) {
            return(
              list(x)
            )
          }
        )
        names(pop) <- gate[[z]]@filterId
      }
    }
    return(pop)
  })
  
  # POPULATION NAMES
  if(length(unique(alias)) == length(pops)) {
    names(pops) <- alias
  }
  
  # LIST OF GATED POPULATIONS PER GATE
  return(pops)
  
}

## CYTO_GATE_REMOVE ------------------------------------------------------------

#' Remove Gate(s) and Edit gatingTemplate csv File
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param parent name of the parent population from which to remove the gate.
#'   This argument is not necessary but is included to allow conversion of
#'   \code{cyto_gate_draw} code to \code{cyto_gate_remove} code by simply
#'   changing \code{"draw"} to \code{"remove"}.
#' @param alias name(s) of the population(s) to remove (e.g. "Single Cells"). By
#'   default all descendant populations will be removed as well.
#' @param channels names of the channel(s) used to gate the population. This
#'   argument is not necessary but is included to allow conversion of
#'   \code{cyto_gate_draw} code to \code{cyto_gate_remove} code by simply
#'   changing \code{"draw"} to \code{"remove"}.
#' @param type gate type(s) used to for the gates to be removed. This argument
#'   is not necessary but is included to allow conversion of
#'   \code{cyto_gate_draw} code to \code{cyto_gate_remove} code by simply
#'   changing \code{"draw"} to \code{"remove"}.
#' @param gatingTemplate name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv").
#' @param ... additional \code{cyto_gate_draw} arguments that will be ignored.
#'
#' @return an object of class \code{GatingSet} with gate and children removed
#'   and updated gatingTemplate to reflect these changes.
#'
#' @importFrom flowWorkspace gh_pop_get_descendants gs_pop_remove
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
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
#' gs <- cyto_compensate(gs)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # Remove T Cells population - replace gatingTemplate name
#' cyto_gate_remove(gs, "T Cells", gatingTemplate = "gatingTemplate.csv")
#' }
#' @export
cyto_gate_remove <- function(x,
                             parent = NULL,
                             alias = NULL,
                             channels = NULL,
                             type = NULL,
                             gatingTemplate = NULL, 
                             ...) {
  
  # ALIAS MISSING
  if (is.null(alias)) {
    stop("Supply the name of the population to be removed to 'alias'.")
  }
  
  # GATINGTEMPLATE MISSING
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  
  # PARSE GATINGTEMPLATE
  gt <- cyto_gatingTemplate_parse(
    gatingTemplate,
    data.table = FALSE
  )
  
  # PREPARE GATINGTEMPLATE ALIAS
  gt_alias <- lapply(
    seq_len(nrow(gt)), 
    function(z) {
      # PARSE PARENT
      prnt <- unlist(strsplit(as.character(gt$parent[z]), ","))
      # PARSE ALIAS
      pops <- unlist(strsplit(as.character(gt$alias[z]), ","))
      # ANCHOR ALIAS
      LAPPLY(
        prnt,
        function(v) {
          cyto_nodes_convert(
            x,
            nodes = pops,
            anchor = v,
            path = "auto"
          )
        }
      )
    }
  )
  
  # ALIAS
  alias <- cyto_nodes_convert(
    x,
    nodes = alias,
    path = "auto",
    anchor = parent
  )
  
  # MULTIPLE ALIAS - REMOVE ALL ASSOCIATED NODES
  alias <- unique(
    LAPPLY(
      seq_len(
        length(alias)
      ), 
      function(z) {
        LAPPLY(
          gt_alias, 
          function(y) {
            if (alias[z] %in% y) {
              y
            } else {
              NA
            }
          }
        )
      }
    )
  )
  alias <- alias[!is.na(alias)]
  
  # CHILDREN
  chldrn <- LAPPLY(
    alias,
    function(z) {
      gh_pop_get_descendants(x[[1]], z, path = "auto")
    }
  )
  chldrn <- unlist(chldrn, use.names = FALSE)
  chldrn <- unique(c(alias, chldrn))
  
  # REMOVE ROWS ALIAS == CHILDREN
  ind <- LAPPLY(
    gt_alias,
    function(z) {
      any(chldrn %in% z)
    }
  )
  
  # UNPARSED GATINGTEMPLATE
  gt <- cyto_gatingTemplate_read(
    gatingTemplate,
    # data.table = FALSE,
    parse = FALSE
  )
  gt <- gt[!ind, , drop = FALSE]
  
  # MESSAGE
  message(
    paste0(
      "Removing the following gates from this ",
      cyto_class(x), " and ", 
      ifelse(is.character(gatingTemplate) , gatingTemplate, "gatingTemplate"), 
      ": \n",
      paste0(
        unlist(gt_alias[ind]),
        collapse = "\n"
      )
    )
  )
  
  # REMOVE NODES FROM GATINGSET
  for (i in seq_len(length(alias))) {
    if (alias[i] %in% cyto_nodes(x, path = "auto")) {
      suppressMessages(gs_pop_remove(x, alias[i]))
    }
  }
  
  # UPDATE GATINGTEMPLATE
  cyto_gatingTemplate_write(
    gt,
    save_as = gatingTemplate
  )
  
  # RETURN GATINGSET
  return(x)
  
}

## CYTO_GATE_RENAME ------------------------------------------------------------

#' Rename Gates
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param alias current names of the gates to be changed.
#' @param names new names to use for \code{alias}.
#' @param gatingTemplate name of the gatingTemplate csv file to be edited.
#' @param ... not in use.
#'
#' @return update gate names in \code{GatingSet} or \code{GatingHierarchy} and
#'   update the gatingTemplate accordingly.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace gh_pop_set_name gs_pop_set_name
#'
#' @export
cyto_gate_rename <- function(x,
                             alias = NULL,
                             names = NULL,
                             gatingTemplate = NULL,
                             ...) {
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  
  # ALIAS INVALID
  alias <- cyto_nodes_convert(
    x,
    nodes = alias,
    path = "auto"
  )
  
  # RENAME GATES IN GATINGHIERARCHY/GATINGSET
  mapply(
    function(alias, name) {
      if (cyto_class(x, "GatingHierarchy")) {
        gh_pop_set_name(x, alias, name)
      } else {
        gs_pop_set_name(x, alias, name)
      }
    },
    alias, 
    names
  )
  
  # READ IN GATINGTEMPLATE
  gt <- cyto_gatingTemplate_read(gatingTemplate)
  
  # UPDATE PARENTAL NAMES
  lapply(
    seq_along(alias), 
    function(z) {
      # CHANGES DESCENDANT NODES AS WELL (eg CD4 T Cells & CD69+ CD4 T Cells)
      if (any(grepl(alias[z], gt$parent, fixed = TRUE))) {
        gt[grepl(alias[z], gt$parent, fixed = TRUE), "parent"] <<- names[z]
      }
    }
  )
  
  # UPDATE ALIAS NAMES
  gt_alias <- lapply(gt$alias, function(z) {
    unlist(strsplit(as.character(z), ","))
  })
  gt_alias <- LAPPLY(
    gt_alias, 
    function(z) {
      lapply(
        seq_len(length(alias)), 
        function(y) {
          # CHANGES DESCENDANT NODES
          if (any(grepl(alias[y], z, fixed = TRUE))) {
            z[grepl(alias[y], z, fixed = TRUE)] <<- names[y]
          }
        }
      )
      # RE-COLLAPSE ALIAS
      z <- paste(z, collapse = ",")
      return(z)
    }
  )
  gt[, "alias"] <- gt_alias
  
  # UPDATE GATING ARGUMENTS - CAPTURES BOOLEAN | CLUSTERS
  lapply(
    seq_len(nrow(gt)), 
    function(z) {
      gating_args <- gt[z, "gating_args"]
      lapply(
        seq_along(alias), 
        function(y) {
          # DESCENDANT NODES CHANGE
          if (any(grepl(alias[y], gating_args, fixed = TRUE))) {
            gating_args <<- gsub(alias[y], names[y], gating_args)
          }
        }
      )
      gt[z, "gating_args"] <<- gating_args
    }
  )
  
  # SAVE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(
    gt,
    save_as = gatingTemplate
  )
  
  # RETURN GATINGSET
  return(x)
  
}

## CYTO_GATE_COPY --------------------------------------------------------------

#' Copy gates to other nodes in GatingSet and update gatingTemplate
#'
#' @param x object of class \code{GatingSet}.
#' @param parent name of the parental node to which the gates should be copied.
#' @param alias vector of names for the copied gates, must be the same length as
#'   copy.
#' @param copy vector of names to copy to new parent.
#' @param gatingTemplate name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the new entries should be saved.
#' @param ... not in use.
#'
#' @return object of class \code{GatingSet} with copied gates applied and update
#'   gatingTemplate csv file with appropriate entries.
#'
#' @importFrom openCyto gs_add_gating_method
#' @importFrom flowCore parameters
#' @importFrom flowWorkspace gs_pop_get_gate
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_gate_copy <- function(x,
                           parent = NULL,
                           alias = NULL,
                           copy = NULL,
                           gatingTemplate = NULL,
                           ...){
  
  # TODO: CHECK IF GATE ALREADY EXISTS - MESSAGE TO USE CYTO_GATE_EDIT()
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET
  if(!cyto_class(x, "GatingSet", TRUE)){
    stop("'x' must be an object of class GatingSet.")
  }
  
  # PARENT
  if (is.null(parent)) {
    stop("Supply the name of the parent population.")
  } else {
    parent <- cyto_nodes_convert(x,
                                 nodes = parent,
                                 path = "auto")
  }
  
  # ALIAS
  if (is.null(copy)) {
    stop("Supply the name of the reference node to 'copy'.")
  } else {
    copy <- cyto_nodes_convert(x,
                               nodes = copy,
                               path = "auto")
  }
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  
  # UPDATE GATINGTEMPLATE ------------------------------------------------------
  
  # NEW ENTRIES
  gt_entries <- lapply(seq_along(alias), function(z){
    suppressMessages(
      gs_add_gating_method(
        gs = x,
        parent = parent,
        alias = alias[z],
        dims = paste(
          parameters(
            gs_pop_get_gate(x, copy[z])[[1]]
          ),
          collapse = ","),
        gating_method = "refGate",
        gating_args = copy[z]
      )
    )
  })
  
  # COMBINE GATINGTEMPLATE ENTRIES
  gt_entries <- do.call("rbind", 
                        gt_entries)
  
  # LOAD GATINGTEMPLATE
  gt <- cyto_gatingTemplate_read(gatingTemplate)
  
  # UPDATE GATINGTEMPLATE
  gt <- rbind(gt, 
              gt_entries)
  
  # WRITE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(gt,
                            save_as = gatingTemplate)
  
  # RETURN UPDATED GATINGSET 
  return(x)
  
}

## CYTO_GATE_BOOL --------------------------------------------------------------

#' Add boolean gate to GatingSet and gatingTemplate
#'
#' @param x object of class \code{GatingSet}.
#' @param parent name of the parent population to which the boolean gates should
#'   be added, set to the most recent common ancestor by default.
#' @param alias vector of names for the boolean populations.
#' @param logic vector of logic to define each of the boolean populations.
#' @param gatingTemplate name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the new entries should be saved.
#' @param ... not in use.
#'
#' @return object of class \code{GatingSet} with new boolean gates and updated
#'   gatingTemplate csv file with appropriate entries.
#'
#' @importFrom flowWorkspace booleanFilter gs_pop_add
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transformations
#' gs <- cyto_transform(gs)
#'
#' # Gating
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#'
#' # Add boolean gate
#' gs <- cyto_gate_bool(gs,
#' alias = "CD4+CD8",
#' logic = "CD4 T Cells|CD8 T Cells")
#' }
#'
#' @export
cyto_gate_bool <- function(x,
                           parent = NULL,
                           alias = NULL,
                           logic = NULL,
                           gatingTemplate = NULL,
                           ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET
  if(!cyto_class(x, "GatingSet", TRUE)){
    stop("'x' must be an object of class GatingSet.")
  }
  
  # ALIAS
  if(is.null(alias)){
    stop("Supply name(s) for the boolean populations to 'alias'.")
  }
  
  # LOGIC
  if(is.null(logic)){
    stop("Supply boolean logic to 'logic' to gate boolean populations.")
  }
  
  # ALIAS & LOGIC
  if(length(alias) != length(logic)){
    stop("'alias' and 'logic' must be the same length.")
  }
  
  # LOGIC - STRIP WHITE SPACE
  logic <- lapply(seq_along(logic), function(z){
    logic_to_clean <- logic[z]
    logic_vector <- LAPPLY(seq_len(nchar(logic_to_clean)), function(z){
      substr(logic_to_clean, z, z)
    })
    logic_match <- which(grepl("\\||\\&|\\!", logic_vector))
    # NO LEADING BOOLEAN OPERATOR
    if(logic_match[1] != 1){
      logic_match <- c(0, logic_match)
    }
    logic_split <- LAPPLY(seq_along(logic_match), function(z){
      # FIRST CHUNK NO OPERATOR
      if(logic_match[z] == 0){
        return(
          trimws(
            substr(logic_to_clean,
                   1,
                   logic_match[z + 1] - 1)
          )
        )
      }
      # CHUNK BOUND BY OPERATORS
      if(length(logic_match) > z){
        # ADJACENT OPERATORS
        if(logic_match[z + 1] - logic_match[z] == 1){
          return(
            substr(logic_to_clean,
                   logic_match[z],
                   logic_match[z])
          )
          # OPERATOR - CHUNK - OPERATOR  
        }else{
          return(
            paste0(
              c(substr(logic_to_clean,
                       logic_match[z],
                       logic_match[z]),
                trimws(
                  substr(logic_to_clean,
                         logic_match[z] + 1,
                         logic_match[z + 1] - 1)
                )), collapse = ""
            )
          )
        }
      }
      # LAST OPERATOR
      if(length(logic_match) == z){
        return(
          paste0(
            c(substr(logic_to_clean,
                     logic_match[z],
                     logic_match[z]),
              trimws(
                substr(logic_to_clean,
                       logic_match[z] + 1,
                       nchar(logic_to_clean))
              )), collapse = ""
          )
        )
      }
    })
    return(logic_split)
  })
  logic <- LAPPLY(logic, function(z){
    paste0(z, collapse = "")
  })
  
  # LOGIC - POPULATIONS
  logic_pop_list <- lapply(logic, function(z){
    p <- unlist(strsplit(z, "\\||\\&|\\!"))
    p <- p[!LAPPLY(p, ".empty")]
    LAPPLY(p, function(y){
      cyto_nodes_convert(x,
                         nodes = y,
                         path = "auto")
    })
  })
  
  # PARENT
  if(is.null(parent)){
    # COMMON ANCESTOR
    parent <- LAPPLY(logic_pop_list, function(pops){
      cyto_nodes_ancestor(x,
                          nodes = pops)
    })
  }else{
    parent <- rep(parent, length.out = length(alias))
  }  
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  
  # READ GATINGTEMPLATE
  gt <- cyto_gatingTemplate_read(gatingTemplate)
  gt_alias <- unlist(strsplit(gt$alias, ","))
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # EMPTY GATINGTEMPLATE ENTRY
  gt_entry <- gt[1, , drop = FALSE]
  gt_entry[1, ] <- rep(NA, ncol(gt_entry))
  
  # ADD BOOLEAN ENTRIES
  gt_entries <- structure(
    lapply(seq_along(alias), function(z){
      # BOOLEAN GATE ALREADY EXISTS - UPDATE LOGIC?
      if(alias[z] %in% gt_alias) {
        # DON'T UPDATE EXISTING BOOLEAN LOGIC
        if(!cyto_enquire(
          paste0(
            alias[z], " already exists in the gatingTemplate. Do you want ",
            "to update the boolean logic for this population? (Y/N)"
          ),
          options = c("Y", "T")
        )) {
          return(NULL)
        # UPDATE EXISTING BOOLEAN LOGIC
        } else {
          # REMOVE EXISTING ENTRY FROM GATINGTEMPLATE
          gt <<- gt[gt$alias != alias[z], ]
        }
      }
      # NEW GATINGTEMPLATE ENTRIES
      gt_pop <- gt_entry
      gt_pop[, "alias"] <- alias[z]
      gt_pop[, "pop"] <- "+"
      gt_pop[, "parent"] <- parent[z]
      gt_pop[, "gating_method"] <- "boolGate"
      gt_pop[, "gating_args"] <- logic[z]
      # CREATE BOOLEANFILTER
      gate <- eval(
        substitute(
          booleanFilter(v, filterId = alias[z]),
          list(v = as.symbol(logic[z]))
        )
      )
      # ADD BOOLEANFILTER TO GATINGSET
      gs_pop_add(x,
                 gate = gate,
                 parent = parent[z])
      # RETURN NEW GATINGTEMPLATE ENTRY
      return(gt_pop)
    }), 
    names = alias)
  gt_entries <- do.call("rbind", gt_entries)
  
  # UPDATE GATINGTEMPLATE
  gt <- rbind(gt, gt_entries)
  cyto_gatingTemplate_write(gt,
                            save_as = gatingTemplate)
  
  # RETURN GATINGSET
  return(x)
  
}

## CYTO_GATE_EXTRACT -----------------------------------------------------------

#' Extract Saved Gate(s) from GatingHierarchy, GatingSet or gatingTemplate.
#'
#' @param x object of class \code{GatingHierarchy}, \code{GatingSet},
#'   \code{gatingTemplate} or the name of a gatingTemplate CSV file.
#' @param parent name of the parental population.
#' @param alias name of the population for which the gate must be extracted.
#' @param select passed to \code{cyto_select}.
#' @param merge_by passed to \code{cyto_groups} to split \code{GatingSet} into
#'   groups, gates will be extracted from the first GatingHierarchy within each
#'   group.
#' @param path passed to \code{cyto_nodes()}.
#' @param bool logical indicating whether booleanFilters should be converted
#'   into complementFilters, set to TRUE by default. Only negated booleanFilters
#'   are supported.
#' @param ... not in use.
#'
#' @importFrom openCyto gt_get_gate gatingTemplate
#' @importFrom flowCore filters parameters<-
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # gatingTemplate
#' gtfile <- system.file("extdata",
#'   "Activation-gatingTemplate.csv",
#'   package = "CytoExploreRData"
#' )
#'
#' # Extract T Cells gate
#' cyto_gate_extract("Live Cells", "T Cells", gatingTemplate = gtfile)
#'
#' @export
cyto_gate_extract <- function(x,
                              parent = NULL,
                              alias,
                              select = NULL,
                              merge_by = "name",
                              path = "auto",
                              bool = TRUE,
                              ...) {
  
  # ALIAS MISSING
  if(missing(alias)) {
    stop("Supply the name(s) of the gate(s) to extract to 'alias'.")
  }
  
  # GATINGHIERARCHY/GATINGSET
  if(cyto_class(x, "GatingSet")) {
    # ANCHOR ALIAS TO PARENT
    alias <- cyto_nodes_convert(
      x,
      nodes = alias,
      anchor = parent,
      path = path
    )
    # SPLIT INTO GROUPS
    x <- cyto_group_by(
      x, 
      group_by = merge_by,
      select = select
    )
    # EXTRACT GATES PER GROUP
    gates <- structure(
      lapply(x, function(z){
        # GATINGHIERARCHY - FIRST GATINGHIERARCHY
        gh <- z[[1]]
        # EXTRACT GATES
        gates <- structure(
          lapply(alias, function(w){
            gate <- gh_pop_get_gate(gh, w)
            # CONVERT BOOLEANFILTER TO VALID FILTER
            if(cyto_class(gate, "booleanFilter") & bool) {
              logic <- gate@deparse
              pops <- unlist(strsplit(logic, "[^[:alnum:][:space:]]"))
              pops <- pops[!LAPPLY(pops, ".empty")]
              # BOOLEAN GATE REFERS TO OTHER POPULATIONS
              if(is.null(
                tryCatch(
                  cyto_nodes_convert(
                    gh,
                    nodes = pops,
                    path = "auto"
                  ),
                  error = function(e) {
                    NULL
                  }
                )
              )) {
                return(gate)
              }
              # EXTRACT OPERATORS
              ops <- unlist(strsplit(logic, "[[:alnum:][:space:]]"))
              ops <- ops[!LAPPLY(ops, ".empty")]
              # EXTRACT GATES
              gts <- structure(
                lapply(
                  pops, 
                  function(q) {
                    gh_pop_get_gate(gh, q)
                  }
                ),
                names = pops
              )
              # CREATE COMBINED COMPLEX FILTER
              gate <- eval(
                parse(
                  text = paste0(ops, "gts[['", pops, "']]", collapse = "")
                )
              )
            }
            return(gate)
          }),
          names = alias
        )
        # RECTANGLEGATES BELONG TO QUADGATE?
        quad_ind <- which(LAPPLY(gates, function(v){
          cyto_class(v, "rectangleGate") &
            any(grepl("quad", names(attributes(v))))
        }))
        # ALL RECTANGLES MUST BE PRESENT
        if(length(quad_ind) == 4) {
          # CREATE QUADGATE
          quad <- .cyto_gate_quad_convert(
            gates[quad_ind]
          )
          # REMOVE RECTANGLEGATES
          gates[quad_ind] <- NULL
          gates[[quad@filterId]] <- quad
        }
        return(gates)
      }),
      names = names(x)
    )
  # GATINGTEMPLATE  
  } else {
    # MISSING PARENT
    if (is.null(parent)) {
      stop("Supply the name of the parent population.")
    }
    # GATINGTEMPLATE FILENAME
    if(is.character(gatingTemplate)) {
      gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
      gatingTemplate <- suppressMessages(gatingTemplate(gatingTemplate))
    }
    # EXTRACT NODES - GATINGTEMPLATE
    nds <- cyto_nodes(gatingTemplate, path = "auto")
    # PARENTAL NODE
    parent <- names(nds[match(parent, nds)])
    # EXTRACT GATES GIVEN PARENTAL & CHILD NODES
    gates <- structure(
      lapply(alias, function(x) {
        # ALIAS NODE
        ind <- LAPPLY(
          seq_len(length(nds)), 
          function(z) {
            if (x %in% nds[[z]]) {
              z
            } else {
              NA
            }
          }
        )
        ind <- ind[!is.na(ind)][1]
        alias <- names(nds[ind])
        gm <- gt_get_gate(gatingTemplate, parent, alias)
        # GATE OBJECT
        if("gate" %in% names(parameters(gm))) {
          gate <- unlist(eval(parameters(gm)$gate))[[1]]
        # BOOLEANFILTER
        } else if(bool) {
          logic <- as.character(parameters(gm)[[1]])
          pops <- unlist(strsplit(logic, "[^[:alnum:][:space:]]"))
          pops <- pops[!LAPPLY(pops, ".empty")]
          # BOOLEAN GATE REFERS TO OTHER POPULATIONS
          if(is.null(
            tryCatch(
              cyto_nodes_convert(
                x,
                nodes = pops,
                path = "auto"
              ),
              error = function(e){
                NULL
              }
            )
          )) {
            return(gate)
          }
          # OPERATORS
          ops <- unlist(strsplit(logic, "[[:alnum:][:space:]]"))
          ops <- ops[!LAPPLY(ops, ".empty")]
          # EXTRACT GATES
          gts <- structure(
            lapply(
              pops, 
              function(q) {
                gt_gates <- gt_get_gate(
                  gatingTemplate, 
                  parent, 
                  names(nds)[match(q, nds)]
                )
                unlist(eval(parameters(gt_gates)$gate))[[1]]
              }
            ),
            names = pops
          )
          # CREATE COMBINED FILTER
          gate <- eval(
            parse(
              text = paste0(ops, "gts[['", pops, "']]", collapse = "")
            )
          )
        }
        return(gate)
      }), 
      names = alias
    )
  }
  
  # EXTRACTED GATES
  return(gates)
}

## CYTO_GATE_EDIT --------------------------------------------------------------

#' Edit gates that exist in a GatingSet and gatingTemplate
#'
#' @param x an object of class \code{GatingSet}.
#' @param parent name of the parental population.
#' @param alias name(s) of the gate to edit (e.g. "Single Cells").
#' @param channels name(s) of the channel(s) used to construct the gate(s). This
#'   argument is not necessary but is included to allow conversion of
#'   \code{cyto_gate_draw} code to \code{cyto_gate_remove} code by simply
#'   changing \code{"draw"} to \code{"remove"}.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   \code{types} are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported \code{gate_types} are \code{polygon, rectangle, ellipse,
#'   threshold, boundary, interval, quadrant and web} which can be abbreviated
#'   as upper or lower case first letters as well. Default \code{type} is
#'   \code{"polygon"}.
#' @param gatingTemplate name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the gate is saved.
#' @param merge_by vector of pData column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to gating, set to the length of x by default to construct a single
#'   gate for all samples. If merge_by is supplied a different gate will be
#'   constructed for each group.
#' @param group_by same as \code{merge_by} included for backwards compatibility
#'   with older versions of CytoExploreR, \code{group_by} was renamed to
#'   \code{merge_by} in CytoExploreR v2.0.0 to maintain consistency with
#'   \code{cyto_plot}. Users should therefore use \code{merge_by} as support for
#'   \code{group_by} will be ended in the future.
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames} or \code{list of flowSets}
#'   containing populations to be overlaid onto the plot(s). 
#' @param select vector containing the indices of samples within gs to use for
#'   plotting.
#' @param events fraction or number of events to display in the plot during the
#'   gating process, set to 50000 events by default.
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
#' @param title title to use for the plot, set to the name of the group by
#'   default. Title can be removed by setting this argument to \code{NA}.
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied data. \code{cyto_plot} does not
#'   support in-line transformations and as such the transformations should be
#'   applied to the data prior to plotting. The transformerList is used
#'   internally to ensure that the axes on the constructed plots are
#'   appropriately labelled.
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
#' @param label_text_font numeric to control the font of text in plot labels,
#'   set to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param label_text_size numeric to control the size of text in the plot
#'   labels, set to 1 by default.
#' @param label_text_col colour(s) to use for text in plot labels, set to
#'   \code{"black"} by default.
#' @param label_text_col_alpha numeric [0, 1] to control the transparency of the
#'   text colour, set to 1 by default to remove transparency.
#' @param label_fill fill colour(s) to use for labels, set to "white" by
#'   default.
#' @param label_fill_alpha numeric to control background fill transparency of
#'   label, set to 0.6 by default to introduce some transparency.
#' @param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'   sampling is applied with each \code{\link{cyto_plot}} call, set to an
#'   arbitrary numeric by default. This behaviour can be turned off by setting
#'   this argument to NULL.
#' @param ... additional arguments for \code{\link{cyto_plot}}.
#'
#' @return an object of class \code{GatingSet} with edited gate applied, as well
#'   as gatingTemplate file with edited gate saved.
#'
#' @importFrom flowWorkspace gh_pop_get_gate gs_pop_set_gate recompute
#'   gh_pop_get_descendants
#' @importFrom flowCore parameters filterList
#' @importFrom openCyto gatingTemplate CytoExploreR_.argDeparser
#' @importFrom data.table as.data.table fread :=
#' @importFrom methods as
#' @importFrom purrr transpose
#' @importFrom methods is
#' @importFrom DataEditR data_edit
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
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
#' gs <- cyto_compensate(gs)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # Edit CD4 T Cells Gate - replace gatingTemplate name
#' cyto_gate_edit(gs,
#'   parent = "T Cells",
#'   alias = "CD4 T Cells",
#'   gatingTemplate = "gatingTemplate.csv"
#' )
#' }
#'
#' @export
cyto_gate_edit <- function(x,
                           parent = NULL,
                           alias = NULL,
                           channels = NULL,
                           type = NULL,
                           gatingTemplate = NULL,
                           overlay = NA,
                           merge_by = "all",
                           group_by = NULL,
                           select = NULL,
                           negate = FALSE,
                           events = 50000,
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
                                 path = "full")
  }
  
  # ALIAS
  if(is.null(alias)) {
    stop("Supply the name(s) of the gates to edit to 'alias'.")
  } else {
    # ALIAS MAY INCLUDE NEW NAME FOR BOOLEAN GATE (NEGATE NOW TRUE)
    alias_new <- NULL
    alias <- LAPPLY(seq_along(alias), function(z){
      pop <- tryCatch(
        cyto_nodes_convert(x,
                           nodes = alias[z],
                           path = "auto"),
        error = function(e){
          return(NULL)
        }
      )
      if(is.null(pop)) {
        alias_new <<- c(alias_new, alias[z])
      }
      return(pop)
    })
  }
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  
  # AXES_TRANS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # GROUP_BY - BACKWARDS COMPATIBILITY
  if(!is.null(group_by)) {
    message(
      "Support for 'group_by' is ending. Please use 'merge_by' instead."
    )
    merge_by <- group_by
  }
  
  # PREPARE GATINGTEMPLATE -----------------------------------------------------
  
  # READ GATINGTEMPLATE
  gt <- cyto_gatingTemplate_read(gatingTemplate,
                                 data.table = FALSE)
  
  # REPLACE PARENT WITH FULL PATHS (IN CASE)
  gt$parent <- cyto_nodes_convert(x, 
                                  nodes = gt$parent,
                                  path = "full")
  
  # RESTRICT GATINGTEMPLATE
  gt_chunk <- 
    gt[gt$parent == parent & 
         LAPPLY(gt$alias, function(z){
           any(alias %in% unlist(strsplit(z, ",")))
         }), ]
  
  # INVALID ALIAS/PARENT COMBINATION
  if(nrow(gt_chunk) == 0) {
    stop(
      "Populations in 'alias' were not gated on the specified 'parent'!"
    )
  }
  
  # EXTRACT CHANNELS
  if (is.null(channels)) {
    channels <- unique(
      unlist(
        strsplit(gt_chunk$dims, ",")
      )
    )
  }
  
  # EXTRACT CHANNELS - REQUIRED FOR QUADGATES
  channels <- cyto_channels_extract(x, channels)
  
  # EXTRACT GROUPBY
  gt_groupBy <- unique(gt_chunk[, "groupBy"])
  
  # EXTRACT GATING METHOD
  gt_gating_method <- unique(gt_chunk[, "gating_method"])
  
  # BOOLEAN GATE(S)
  gt_bool_chunk <- gt_chunk[gt_chunk$gating_method == "boolGate", ]
  if(nrow(gt_bool_chunk) > 0) {
    # REMOVE FROM ALIAS
    bool_alias <- gt_bool_chunk$alias
    alias <- alias[!alias %in% bool_alias]
  } else {
    bool_alias <- NULL
  }
  
  # NAME FOR NEW BOOLEAN GATE
  if(is.null(bool_alias) & negate == TRUE) {
    if(is.null(alias_new) | length(alias_new) > 1) {
      alias_new <- cyto_enquire(
        "Provide a name for the new boolean gate:"
      )
    }
  }
  
  # REFERENCE GATES - MESSAGE
  if(any(gt_chunk$gating_method == "refGate")) {
    message(
      "Reference gate(s) will be replaced with their own gate(s)."
    )
  }
  
  # EXTRACT GATES FROM GATINGSET DIRECTLY --------------------------------------
  
  # GATES PER GROUP - BOOLEAN GATES EXCLUDED
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
                              events = events,
                              seed = seed)
  
  # SELECT GROUP(S) TO EDIT ----------------------------------------------------
  
  # GROUPS
  pd_groups <- cyto_groups(x,
                           select = select,
                           group_by = merge_by,
                           details = TRUE)
  
  # PREPARE GROUP NAMES
  grps <- names(cs_lists)
  grps[is.na(grps)] <- "NA"
  
  # MENU - SELECT GROUPS TO EDIT
  if(merge_by[1] == "all") {
    merge_by <- "Combined Events" # match cyto_groups
    grps <- "Combined Events"
    # SKIP MENU IF NEW GROUPING ADDED
  } else if(.all_na(gt_groupBy) & merge_by[1] != "all") {
    grps <- grps
    # SKIP MENU IF GROUP HAS CHANGED
  } else if(!.all_na(gt_groupBy) &
            !setequal(merge_by, unlist(strsplit(gt_groupBy, ":")))) {
    grps <- grps
    # GROUPING REMAINS UNCHANGED
  } else {
    message("Select the group(s) to edit:")
    # INTERACTIVE GROUP SELECTION
    if(interactive()) {
      grps <- data.frame("group" = grps,
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
                        viewer = "pane",
                        hide = TRUE,
                        quiet = TRUE)
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
    # EXISTING GATE LIST - EXLUDES BOOLEAN GATES
    gate <- gates_gs[[y]]
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
                overlay = if(length(cs_list) > 1){
                  cs_list[seq_along(cs_list)[-1]]
                } else {
                  NA
                },
                channels = channels,
                axes_trans = axes_trans,
                axes_limits = axes_limits,
                legend = FALSE,
                gate = gate,
                gate_line_width = 2.5,
                gate_line_col = "magenta",
                gate_line_alpha = 0.8,
                label = FALSE,
                title = title[w],
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
                               alias = if(negate){
                                 c(unlist(alias), bool_alias, alias_new)
                               } else {
                                 unlist(alias)
                               },
                               channels = channels,
                               type = type,
                               negate = negate,
                               axis = axis,
                               plot = FALSE,
                               gate_point_shape = gate_point_shape,
                               gate_point_size = gate_point_size,
                               gate_point_col = gate_point_col,
                               gate_point_col_alpha = gate_point_col_alpha,
                               gate_line_type = gate_line_type,
                               gate_line_width = gate_line_width,
                               gate_line_col = gate_line_col,
                               gate_line_col_alpha = gate_line_col_alpha,
                               label_text_size = label_text_size,
                               label_text_font = label_text_font,
                               label_text_col = label_text_col,
                               label_text_col_alpha = label_text_col_alpha,
                               label_fill = label_fill,
                               label_fill_alpha = label_fill_alpha)[[1]]
    
    # REMOVE BOOLEAN GATE
    if(negate == TRUE) {
      if(!is.null(bool_alias) | !is.null(alias_new)) {
        gate_new <- gate_new[-length(gate_new)]
      }
    }
    
    # UPDATE GATES
    gates_gs[[y]] <<- gate_new
  })
  
  # PREPARE GATES FOR GATINGTEMPLATE
  gates_gs_transposed <- transpose(gates_gs)
  
  # WRAP GATES IN FILTER LIST
  gates_gs_transposed <- structure(
    lapply(seq_along(gates_gs_transposed), function(z){
      gates <- structure(
        lapply(seq_along(gates_gs_transposed[[z]]), function(y){
          if(!cyto_class(gates_gs_transposed[[z]][[y]], "quadGate")) {
            filters(gates_gs_transposed[[z]][y])
          } else {
            gates_gs_transposed[[z]][[y]] # here
          }
        }),
        names = names(gates_gs_transposed[[z]])
      )
    }), 
    names = names(gates_gs_transposed)
  )
  
  # DATA.TABLE FRIENDLY NAMES
  prnt <- parent
  als <- alias
  gtmd <- "cyto_gate_draw"
  ppmd <- "pp_cyto_gate_draw"
  
  # FIND & EDIT GATINGTEMPLATE ENTRIES
  gt <- cyto_gatingTemplate_read(gatingTemplate, data.table = TRUE)
  
  # DATA.TABLE R CMD CHECK NOTE
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  preprocessing_args <- NULL
  .SD <- NULL
  
  # GROUP_BY
  if (merge_by[1] == "Combined Events") {
    group_by <- "NA"
  } else {
    group_by <- paste(merge_by, collapse = ":")
  }
  
  # HANDLE MULTIPLE POPULATIONS
  als <- lapply(als, "paste", collapse = ",")

  # INDEX PARENT IN GATINTEMPLATE (FULL PATH)
  parent_ind <- match(
    prnt,
    cyto_nodes_convert(x, 
                       nodes = gt$parent,
                       path = "full")
  )
  # MODIFY GATINGTEMPLATE - GATED POPULATIONS ONLY
  for (i in seq_len(length(als))) {
    alias_ind <- match(
      als[i],
      gt$alias
    )
    ind <- intersect(parent_ind, alias_ind)
    gt[ind, gating_method := gtmd]
    gt[
      ind,
      gating_args := CytoExploreR_.argDeparser(list(
        gate = gates_gs_transposed[[i]],
        openCyto.minEvents = -1
      ))
    ]
    gt[ind, collapseDataForGating := TRUE]
    gt[ind, preprocessing_method := ppmd]
    gt[ind, preprocessing_args := as.logical(NA)]
    # groupBy must be character class
    gt[, groupBy := lapply(.SD, as.character), .SDcols = "groupBy"]
    gt[ind, groupBy := group_by]
  }
  
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
  
  # APPLY NEW GATES TO GATINGSET
  lapply(grps, function(z) {
    # MATCHING GROUPS
    ind <- rownames(pd_groups[[z]])
    # SET GATES - GATED POPULATIONS ONLY
    lapply(seq_along(alias), function(y) {
      # REPEAT GATES FOR EACH SAMPLE IN GROUP
      gates_update <- rep(
        list(gates_gs[[z]][[alias[y]]]),
        length(x[ind])
      )
      names(gates_update) <- cyto_names(x[ind])
      # ANCHOR
      alias_anchor <- cyto_nodes_convert(x,
                                         nodes = alias[y],
                                         anchor = parent)
      # UPDATE GATE
      suppressMessages(gs_pop_set_gate(
        x[ind],
        alias_anchor,
        gates_update
      ))
      # RECOMPUTE STATISTICS
      suppressMessages(recompute(x[ind], alias_anchor))
    })
  })
  
  # BOOLEAN GATE REQUIRED
  if(negate == TRUE) {
    # BOOLEAN GATE ADDED
    if(is.null(bool_alias)) {
      if(is.null(alias_new)) {
        alias_new <- cyto_enquire(
          "Please provide a name for new new boolean gate:"
        )
      }
      # GATINGSET
      pop <- suppressMessages(
        gs_add_gating_method(
          gs = x,
          alias = alias_new,
          parent = parent,
          pop = "+",
          dims = paste(channels, collapse = ","),
          gating_method = "boolGate",
          gating_args = paste0("!", paste(unlist(alias), collapse = "&!")),
          groupBy = group_by,
          collapseDataForGating = TRUE,
          preprocessing_method = NA
        )
      )
      # GATINGTEMPLATE
      gt <- rbind(gt, pop)
      # UPDATE EXISTING BOOLEAN GATE  
    } else{
      suppressMessages(recompute(x, cyto_nodes_convert(x,
                                                       nodes = alias_new,
                                                       anchor = parent)))
    }
    # REMOVE BOOLEAN GATE
  } else {
    # BOOLEAN ENTRY EXISTS
    if(nrow(gt_bool_chunk) > 0) {
      # REMOVE BOOLEAN GATE(S)?
      if(cyto_enquire(
        paste0(
          "Do you want to remove the boolean gate(s) and their descendants ",
          "from the GatingSet and gatingTemplate? (Y/N)"
        ),
        options = c("Y", "T")
      )) {
        for(i in 1:nrow(gt_bool_chunk)) {
          # DESCENDANTS
          bool_desc <- gh_pop_get_descendants(x[[1]], 
                                              gt_bool_chunk[i, "alias"],
                                              path = "auto")
          # REMOVE BOOLEAN GATE & DESCENDANTS
          gt <- gt[!alias %in% c(gt_bool_chunk[i, "alias"], bool_desc), ]
          # REMOVE GATE(S) FROM GATINGSET
          suppressMessages(gs_pop_remove(x, gt_bool_chunk[i, "alias"]))
        }
      }
    }
  }
  
  # # BOOLEAN GATE(S) - REMOVE OR MODIFY LOGIC
  # if(length(bool_ind) > 0) {
  #   # UPDATE BOOLEAN GATE LOGIC
  #   if(!is.null(logic)) {
  #     for(i in 1:nrow(gt_bool)) {
  #       # GATINGTEMPLATE
  #       gt[parent == gt_bool_chunk[i, "parent"] & 
  #            alias == gt_bool_chunk[i, "alias"], 
  #          gating_args := gt_bool_logic[i]]
  #       # GATINGSET
  #       gs_pop_set_gate(x,
  #                       cyto_nodes_convert(x,
  #                                         nodes = gt_bool_chunk[i, "alias],
  #                                         anchor = gt_bool_chunk[i, "parent]),
  #                       eval(
  #                         substitute(
  #                           booleanFilter(v, 
  #                                         filterId = gt_bool_chunk[i, "alias"]),
  #                           list(v = as.symbol(logic[i]))
  #                         )
  #                       )
  #       )
  #     }
  #   }
  # }
  
  # SAVE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(gt,
                            save_as = gatingTemplate)
  
  # UPDATE GATINGSET GLOBALLY
  return(x)
  
}

## CYTO_GATE_TYPE --------------------------------------------------------------

#' Get Gate Type from Saved Gate.
#'
#' @param gates an object of class \code{filters} containing the gates from
#'   which the \code{type(s)} must be obtained.
#'
#' @return vector of gate type names for supplied gates.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # Get gate type used for T Cells gate
#' cyto_gate_type(gs_pop_get_gate(gs, "Cells")[[1]])
#' cyto_gate_type(gs_pop_get_gate(gs, "T Cells")[[1]])
#' cyto_gate_type(gs_pop_get_gate(gs, "CD69+ CD4 T Cells")[[1]])
#' @export
cyto_gate_type <- function(gates) {
  
  # One gate supplied
  if (length(gates) == 1) {
    if (cyto_class(gates, c("filters", "list"), TRUE)) {
      gates <- gates[[1]]
    }
    # ELLIPSE
    if (cyto_class(gates, "ellipsoidGate")) {
      types <- "ellipse"
      # RECTANGLE/BOUNDARY/INTERVAL/THRESHOLD
    } else if (cyto_class(gates, "rectangleGate")) {
      # Includes rectangle, interval, threshold and boundary gate_types
      if (length(parameters(gates)) == 1) {
        
        # Gate in One Dimension
        if (is.infinite(gates@min)) {
          types <- "boundary"
        } else if (is.infinite(gates@max)) {
          types <- "threshold"
        } else if (is.finite(gates@min) & is.finite(gates@max)) {
          types <- "interval"
        }
      } else if (length(parameters(gates)) == 2) {
        # Gate in 2 Dimensions
        if (is.infinite(gates@min[1]) & is.infinite(gates@min[2])) {
          types <- "boundary"
        } else if (is.infinite(gates@max[1]) & is.infinite(gates@max[2])) {
          types <- "threshold"
        } else if (all(is.infinite(c(gates@min[1], gates@max[1]))) |
                   all(is.infinite(c(gates@min[2], gates@max[2])))) {
          types <- "interval"
        } else {
          types <- "rectangle"
        }
      }
      # POLYGON
    } else if (cyto_class(gates, "polygonGate")) {
      types <- "polygon"
      # QUADRANT
    } else if (cyto_class(gates, "quadGate")) {
      types <- "quadrant"
    }
    # Multiple gates supplied
  } else if (length(gates) > 1) {
    # Get classes of gates
    classes <- LAPPLY(gates, function(x) {
      cyto_class(x, class = TRUE)
    })
    # All gates are of the same class
    if (all(classes[1] == classes)) {
      # Gates are all ellipses
      if (classes[1] == "ellipsoidGate") {
        types <- rep("ellipse", length(gates))
        # Gates are all rectangles
      } else if (classes[1] == "rectangleGate") {
        # Combine gate co-ordinates
        pts <- lapply(gates, function(x) {
          rbind(x@min, x@max)
        })
        pts <- do.call(rbind, pts)
        # if 4 gates are supplied - type may be "quadrant"
        if (length(gates) == 4) {
          # Quadrant gates should have finite and infinite values in all gates
          # and all finite co-ordinates should be the same
          if (sum(is.finite(pts[, 1])) == 4 &
              sum(is.infinite(pts[, 1])) == 4 &
              sum(duplicated(pts[, 1][is.finite(pts[, 1])])) == 3) {
            types <- "quadrant"
            # Each gate could be either rectangle, interval, threshold, boundary
          } else {
            types <- LAPPLY(gates, function(x) {
              # Includes rectangle, interval, threshold and boundary gate_types
              if (length(parameters(x)) == 1) {
                # Gate in One Dimension
                if (is.infinite(x@min)) {
                  types <- "boundary"
                } else if (is.infinite(x@max)) {
                  types <- "threshold"
                } else if (is.finite(x@min) & is.finite(x@max)) {
                  types <- "interval"
                }
              } else if (length(parameters(x)) == 2) {
                # Gate in 2 Dimensions
                if (is.infinite(x@min[1]) & is.infinite(x@min[2])) {
                  types <- "boundary"
                } else if (is.infinite(x@max[1]) & is.infinite(x@max[2])) {
                  types <- "threshold"
                } else if (all(is.infinite(c(x@min[1], x@max[1]))) |
                           all(is.infinite(c(x@min[2], x@max[2])))) {
                  types <- "interval"
                } else {
                  types <- "rectangle"
                }
              }
            })
          }
        } else {
          types <- LAPPLY(gates, function(x) {
            # Includes rectangle, interval, threshold and boundary gate_types
            if (length(parameters(x)) == 1) {
              # Gate in One Dimension
              if (is.infinite(x@min)) {
                types <- "boundary"
              } else if (is.infinite(x@max)) {
                types <- "threshold"
              } else if (is.finite(x@min) & is.finite(x@max)) {
                types <- "interval"
              }
            } else if (length(parameters(x)) == 2) {
              # Gate in 2 Dimensions
              if (is.infinite(x@min[1]) & is.infinite(x@min[2])) {
                types <- "boundary"
              } else if (is.infinite(x@max[1]) & is.infinite(x@max[2])) {
                types <- "threshold"
              } else if (all(is.infinite(c(x@min[1], x@max[1]))) |
                         all(is.infinite(c(x@min[2], x@max[2])))) {
                types <- "interval"
              } else {
                types <- "rectangle"
              }
            }
          })
        }
        # Gates are all polygons
      } else if (classes[1] == "polygonGate") {
        # Combine gate co-ordinates
        pts <- lapply(gates, function(x) {
          rbind(x@boundaries)
        })
        pts <- do.call(rbind, pts)
        dupl <- pts[pts[, 1] == pts[which(duplicated(pts))[1], ][1], ]
        # May be type == "web" need to see if any points are conserved
        if (length(dupl[, 1]) == length(gates)) {
          types <- "web"
        } else {
          types <- rep("polygon", length(gates))
        }
      }
      # Not all supplied gates are of the same class - treat separately
    } else {
      types <- LAPPLY(gates, function(x) {
        # ELLIPSE
        if (class(x) == "ellipsoidGate") {
          types <- "ellipse"
          # RECTANGLE/BOUNDARY/THRESHOLD/INTERVAL
        } else if (class(x) == "rectangleGate") {
          # Includes rectangle, interval, threshold and boundary gate_types
          if (length(parameters(x)) == 1) {
            # Gate in One Dimension
            if (is.infinite(x@min)) {
              types <- "boundary"
            } else if (is.infinite(x@max)) {
              types <- "threshold"
            } else if (is.finite(x@min) & is.finite(x@max)) {
              types <- "interval"
            }
          } else if (length(parameters(x)) == 2) {
            # Gate in 2 Dimensions
            if (is.infinite(x@min[1]) & is.infinite(x@min[2])) {
              types <- "boundary"
            } else if (is.infinite(x@max[1]) & is.infinite(x@max[2])) {
              types <- "threshold"
            } else if (all(is.infinite(c(x@min[1], x@max[1]))) |
                       all(is.infinite(c(x@min[2], x@max[2])))) {
              types <- "interval"
            } else {
              types <- "rectangle"
            }
          }
          # POLYGON
        } else if (class(x) == "polygonGate") {
          types <- "polygon"
          # QUADRANT
        } else if (class(x) == "quadGate") {
          types <- "quadrant"
        }
      })
    }
  }
  
  return(types)
}

## CYTO_GATE_CONVERT -----------------------------------------------------------

# CYTO_GATE_CONVERT is designed to convert constructed gate objects to the
# number of dimensions specified by channels. Primary use is to easily convert
# between 1D and 2D gate objects.

#' Convert Between 1D and 2D Gate Objects
#'
#' Useful function to convert between 1D and 2D gate objects.
#'
#' @param x gate object(s) to be converted.
#' @param channels indicates the required dimensions of the gate for plotting.
#' @param ... not in use.
#'
#' @return modified gate object with appropriate dimensions.
#'
#' @importFrom flowCore parameters rectangleGate
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_gate_convert
NULL

#' @noRd
#' @export
cyto_gate_convert <- function(x, ...) {
  UseMethod("cyto_gate_convert")
}

#' @noRd
#' @export
cyto_gate_convert.default <- function(x,
                                      channels = NULL,
                                      ...) {
  
  # # Invalid gate object
  # if (!cyto_class(x, c("list",
  #                      "filters",
  #                      "rectangleGate",
  #                      "polygonGate",
  #                      "ellipsoidGate",
  #                      "quadGate"), TRUE)) {
  #   stop(paste(
  #     "'x' should contain either filters, rectangleGate, polygonGate,",
  #     "ellipsoidGate or quadGate objects."
  #   ))
  # }
  
  return(x)
  
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.rectangleGate <- function(x,
                                            channels = NULL,
                                            ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (is.null(channels)) {
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # 1D PLOT
  if (length(channels) == 1) {
    # 1D GATE IN 1D PLOT - CORRECT CHANNEL
    if (length(chans) == 1 & all(chans %in% channels)) {
      return(x)
      # 1D GATE IN 1D PLOT - INCORRECT CHANNEL
    } else if (length(chans) == 1 & all(!chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channel.")
      # 2D GATE IN 1D PLOT - CORRECT CHANNEL
    } else if (length(chans) == 2 & any(chans == channels)) {
      return(x[channels])
      # 2D GATE IN 1D PLOT - INCORRECT CHANNEL
    } else if (length(chans) == 2 & !any(chans == channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channel.")
    }
    # 2D PLOT
  } else if (length(channels) == 2) {
    # 2D GATE IN 2D PLOT - CORRECT CHANNELS
    if (length(chans) == 2 & all(chans %in% channels)) {
      return(x)
      # 2D GATE in 2D PLOT - ONE CHANNEL MATCH
    } else if (length(chans) == 2 & any(chans %in% channels)) {
      # FITERID
      ID <- x@filterId
      # CHANNEL INDEX
      ind <- match_ind(chans, channels)
      # COORDS - MISSING CHANNEL
      coords <- matrix(c(x@min[channels[ind]], x@max[channels[ind]], -Inf, Inf),
                       ncol = 2,
                       nrow = 2,
                       byrow = FALSE
      )
      colnames(coords) <- c(channels[ind], channels[-ind])
      x <- rectangleGate(.gate = coords)
      return(x)
      # 2D GATE IN 2D PLOT - NO CHANNELS MATCH
    } else if (length(chans) == 2 & !any(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
      # 1D GATE IN 2D PLOT - CORRECT CHANNEL
    } else if (length(chans) == 1 & all(chans %in% channels)) {
      # FITERID
      ID <- x@filterId
      # CHANNEL INDEX
      ind <- match_ind(chans, channels)
      # COORDS - MISSING CHANNEL
      coords <- matrix(c(x@min, x@max, -Inf, Inf),
                       ncol = 2,
                       nrow = 2,
                       byrow = FALSE
      )
      colnames(coords) <- c(channels[ind], channels[-ind])
      x <- rectangleGate(.gate = coords)
      return(x)
      # 1D GATE IN 2D PLOT - INCORRECT CHANNEL
    } else if (length(chans) == 1 & !all(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.polygonGate <- function(x,
                                          channels = NULL,
                                          ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (is.null(channels)) {
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # 1D PLOT - CONVERT TO RECTANGLEGATE
  if (length(channels) == 1) {
    # 2D GATE IN 1D PLOT - CORERCT CHANNEL
    if (any(chans %in% channels)) {
      # FILTERID
      ID <- x@filterId
      # COORDS IN CHANNEL
      coords <- x@boundaries[, channels[channels %in% chans]]
      coords <- c(min(coords), max(coords))
      # RECTANGLEGATE
      coords <- matrix(coords,
                       ncol = 1,
                       nrow = 2
      )
      colnames(coords) <- channels[channels %in% chans]
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 1D RECTANGLEGATE
      return(x)
      # 2D GATE IN 1D PLOT - INCORRECT CHANNELS
    } else if (!any(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
    # 2D PLOT
  } else if (length(channels) == 2) {
    # 2D GATE IN 2D PLOT - CORRECT CHANNELS
    if (all(chans %in% channels)) {
      return(x)
      # 2D GATE IN 2D PLOT - SINGLE CHANNEL
    } else if (any(chans %in% channels)) {
      # FILTERID
      ID <- x@filterId
      # COORDS IN CHANNEL
      coords <- x@boundaries[, channels[channels %in% chans]]
      coords <- c(min(coords), max(coords))
      # COORDS in MISSING CHANNEL
      coords <- c(coords, -Inf, Inf)
      # RECTANGLEGATE
      coords <- matrix(coords,
                       ncol = 2,
                       nrow = 2,
                       byrow = FALSE
      )
      colnames(coords) <- c(
        channels[channels %in% chans],
        !channels[channels %in% chans]
      )
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 2D RECTANGLEGATE
      return(x)
      # 2D GATE IN 2D PLOT - INCORRECT CHANNELS
    } else if (!any(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.ellipsoidGate <- function(x,
                                            channels = NULL,
                                            ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (is.null(channels)) {
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # 1D PLOT
  if (length(channels) == 1) {
    # 2D GATE IN 1D PLOT - CHANNEL MATCH
    if (any(chans %in% channels)) {
      # FILTERID
      ID <- x@filterId
      # POLYGONGATE
      x <- as(x, "polygonGate")
      # COORDS IN CHANNEL
      coords <- x@boundaries[, channels[channels %in% chans]]
      coords <- c(min(coords), max(coords))
      # RECTANGLEGATE
      coords <- matrix(coords,
                       ncol = 1,
                       nrow = 2
      )
      colnames(coords) <- channels[channels %in% chans]
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 1D RECTANGLEGATE
      return(x)
      # 2D GATE IN 1D PLOT - NO CHANNEL MATCH
    } else if (!any(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
    # 2D PLOT
  } else if (length(channels) == 2) {
    # 2D GATE IN 2D PLOT - CORRECT CHANNELS
    if (all(chans %in% channels)) {
      return(x)
      # 2D GATE IN 2D PLOT - SINGLE CHANNEL MATCH
    } else if (any(chans %in% channels)) {
      # FILTERID
      ID <- x@filterId
      # POLYGONGATE
      x <- as(x, "polygonGate")
      # COORDS IN CHANNEL
      coords <- x@boundaries[, channels[channels %in% chans]]
      coords <- c(min(coords), max(coords))
      # COORDS IN MISSING CHANNEL
      coords <- c(coords, -Inf, Inf)
      # RECTANGLEGATE
      coords <- matrix(coords,
                       ncol = 2,
                       nrow = 2,
                       byrow = FALSE
      )
      colnames(coords) <- c(
        channels[channels %in% chans],
        channels[!channels %in% chans]
      )
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 2D RECTANGLEGATE
      return(x)
      # 2D GATE IN 2D PLOT - NO CHANNEL MATCH
    } else if (!any(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.quadGate <- function(x,
                                       channels = NULL,
                                       ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (is.null(channels)) {
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)
  
  # PREPARE GATES --------------------------------------------------------------
  
  # 1D PLOT
  if (length(channels) == 1) {
    stop("Quadrant gates cannot be converted into 1D gate objects.")
    # 2D PLOT
  } else if (length(channels) == 2) {
    # 2D GATES IN 2D PLOT - CORRECT CHANNELS
    if (all(chans %in% channels)) {
      return(x)
      # 2D GATES IN 2D PLOT - INCORRECT CHANNELS
    } else if (!any(chans %in% channels)) {
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.filters <- function(x,
                                      channels = NULL,
                                      ...) {
  
  # GATE OBJECT LIST -----------------------------------------------------------
  x <- unlist(x)
  
  # LIST METHOD CALL -----------------------------------------------------------
  x <- cyto_gate_convert(
    x,
    channels
  )
  
  # RETURN CONVERTED GATES -----------------------------------------------------
  return(x)
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.list <- function(x,
                                   channels = NULL,
                                   ...) {
  
  # GATE OBJECT LIST -----------------------------------------------------------
  x <- unlist(x)
  
  # CONVERT GATES --------------------------------------------------------------
  x <- structure(
    lapply(x, function(z) {
      cyto_gate_convert(z, channels)
    }),
    names = names(x)
  )
  
  # RETURN CONVERTED GATES -----------------------------------------------------
  return(x)
}

## CYTO_GATE_PREPARE -----------------------------------------------------------

#' Prepare a list of gate objects
#'
#' Returns a list of modified gate objects with appropriate dimensions.
#'
#' @param x gate object(s) to be converted.
#' @param channels indicates the required dimensions of the gate for plotting.
#'
#' @return a list of unique modified gate objects.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_gate_convert
cyto_gate_prepare <- function(x,
                              channels = NULL) {
  
  # LIST OF GATES --------------------------------------------------------------
  
  # PREPARE GATE LIST
  if (cyto_class(x, "list", TRUE)) {
    if (!all(LAPPLY(x, "cyto_class") %in% c(
      "rectangleGate",
      "polygonGate",
      "ellipsoidGate",
      "quadGate",
      "filters"
    ))) {
      x <- unlist(x)
    }
  } else if (cyto_class(x, "filters", TRUE)) {
    x <- unlist(x)
  } else if (cyto_class(x, c("rectangleGate",
                             "polygonGate",
                             "ellipsoidGate",
                             "quadGate",
                             "filters"), TRUE)) {
    x <- list(x)
  }
  
  # UNIQUE GATE LIST - UNIQUE() DROP NAMES
  x <- x[!duplicated(x)]
  
  # DIMENSIONS
  x <- cyto_gate_convert(x, channels = channels)
  
  # RETURN PREPARED GATE LIST
  return(x)
}

## CYTO_GATE_TRANSFORM ---------------------------------------------------------

#' Transform gates objects
#'
#' @param x object of class rectangleGate, polygonGate, ellipsoidGate, quadGate,
#'   filters or a list of these gate objects.
#' @param trans list of transformers to apply to the gate co-ordinates.
#' @param inverse logical indicating whether the inverse transformations should
#'   be applied, set to FALSE by default.
#' @param ... not in use.
#'
#' @return transformed gate object.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowCore parameters rectangleGate polygonGate ellipsoidGate
#'   quadGate
#' @importFrom stats dist
#'
#' @name cyto_gate_transform
NULL

#' @noRd
#' @export
cyto_gate_transform <- function(x, ...){
  UseMethod("cyto_gate_transform")
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.rectangleGate <- function(x,
                                              trans = NULL,
                                              inverse = FALSE,
                                              ...){
  
  # TRANSFORMERS
  if(is.null(trans) | !cyto_class(trans, "transformerList")){
    stop("Supply a list of transformers to transform the gate co-ordinates.")
  }
  
  # GATE FILTERID
  gate_name <- x@filterId
  
  # GATE CHANNELS
  gate_channels <- parameters(x) 
  
  # GATE COORDS
  gate_coords <- .cyto_gate_coords(list(x), channels = gate_channels)
  
  # CHANNELS TO TRANSFORM
  gate_trans_channels <- gate_channels[which(gate_channels %in% names(trans))]
  
  # TRANSFORMATIONS
  lapply(gate_trans_channels, function(z){
    # TRANSFORM COORDS - WATCH OUT INF COORDS
    if(inverse == FALSE){
      gate_coords[, z] <<- LAPPLY(gate_coords[,z], function(y){
        if(is.finite(y)){
          trans[[z]]$transform(y)
        }else{
          y
        }
      })
    }else{
      gate_coords[, z] <<- LAPPLY(gate_coords[,z], function(y){
        if(is.finite(y)){
          trans[[z]]$inverse(y)
        }else{
          y
        }
      })
    }
  })
  
  # UPDATE GATE
  x <- rectangleGate(filterId = gate_name,
                     .gate = gate_coords)
  return(x)
  
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.polygonGate <- function(x,
                                            trans = NULL,
                                            inverse = FALSE,
                                            ...){
  
  # TRANSFORMERS
  if(is.null(trans) | !cyto_class(trans, "transformerList")){
    stop("Supply a list of transformers to transform the gate co-ordinates.")
  }
  
  # GATE FILTERID
  gate_name <- x@filterId 
  
  # GATE CHANNELS
  gate_channels <- parameters(x) 
  
  # GATE COORDS
  gate_coords <- .cyto_gate_coords(list(x), channels = gate_channels)
  
  # CHANNELS TO TRANSFORM
  gate_trans_channels <- gate_channels[which(gate_channels %in% names(trans))]
  
  # TRANSFORMATIONS
  lapply(gate_trans_channels, function(z){
    # TRANSFORM COORDS - WATCH OUT INF COORDS
    if(inverse == FALSE){
      gate_coords[, z] <<- LAPPLY(gate_coords[,z], function(y){
        if(is.finite(y)){
          trans[[z]]$transform(y)
        }else{
          y
        }
      })
    }else{
      gate_coords[, z] <<- LAPPLY(gate_coords[,z], function(y){
        if(is.finite(y)){
          trans[[z]]$inverse(y)
        }else{
          y
        }
      })
    }
  })
  
  # UPDATE GATE
  x <- polygonGate(filterId = gate_name,
                   .gate = gate_coords)
  return(x)
  
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.ellipsoidGate <- function(x,
                                              trans = NULL,
                                              inverse = FALSE,
                                              ...) {
  
  # TODO: APPROACH DOES NOT WORK SINCE SYMMETRY IS LOST DURING TRANSFORMATION
  
  # TRANSFORMERS
  if(is.null(trans) | !cyto_class(trans, "transformerList")){
    stop("Supply a list of transformers to transform the gate co-ordinates.")
  }
  
  # GATE FILTERID
  name <- x@filterId  
  
  # GATE CHANNELS
  channels <- parameters(x)  
  
  # GATE COVARIANCE MATRIX
  cov <- x@cov
  
  # SORT CHANNELS BY COVARIANCE
  chans <- names(sort(diag(cov), decreasing = FALSE))
  
  # GATE CENTER
  center <- x@mean
  
  # DISTANCE
  dist <- x@distance
  
  # EIGEN DECOMPOSITION
  eig <- eigen(cov)
  
  # EIGEN VALUES - MINOR = LARGER EIGENVALUE
  evals <- eig$values
  names(evals) <- chans
  
  # EIGEN VECTORS
  evecs <- eig$vectors
  colnames(evecs) <- chans
  
  # MAJOR AXIS LENGTH
  a2 <- evals[2]*(dist^2)
  a <- sqrt(a2)
  a_vec <- a*evecs[, 2]
  
  # MINOR AXIS LENGTH
  b2 <- evals[1]*(dist^2)
  b <- sqrt(b2)
  b_vec <- b*evecs[, 1]
  
  # ANTIPODAL CO-ORDINATES
  coords <- matrix(
    c(
      c(center[1] + a_vec[1], center[2] + a_vec[2]),
      c(center[1] - a_vec[1], center[2] - a_vec[2]),
      c(center[1] + b_vec[1], center[2] + b_vec[2]),
      c(center[1] - b_vec[1], center[2] - b_vec[2])
    ),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, chans)
  )
  
  # APPEND CENTER
  coords <- rbind(coords, center[chans])

  # TRANSFORM CO-ORDINATES
  cnt <- 1
  coords <- apply(
    coords,
    2,
    function(z) {
      if(colnames(coords)[cnt] %in% names(trans)) {
        if(!inverse) {
          z <- trans[[colnames(coords)[cnt]]]$transform(z)
        } else {
          z <- trans[[colnames(coords)[cnt]]]$inverse(z)
        }
      }
      cnt <<- cnt + 1
      return(z)
    }
  )

  
  # NEW TRANSFORMED CENTER
  center <- coords[nrow(coords), ]
  coords <- coords[-nrow(coords), ]
  
  # RECOMPUTE COVARIANCE MATRIX
  cov <- .cyto_ellipse_cov(coords)

  # UPDATE GATE
  x <- ellipsoidGate(
    filterId = name,
    .gate = cov,
    mean = center
  )
  return(x)
  
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.quadGate <- function(x, 
                                         trans = NULL,
                                         inverse = FALSE,
                                         ...){
  
  # TRANSFORMERS
  if(is.null(trans) | !cyto_class(trans, "transformerList")){
    stop("Supply a list of transformers to transform the gate co-ordinates.")
  }
  
  # GATE FILTERID
  gate_name <- x@filterId
  
  # GATE CHANNELS
  gate_channels <- parameters(x) 
  
  # GATE COORDS
  gate_coords <- .cyto_gate_coords(list(x), channels = gate_channels)
  
  # CHANNELS TO TRANSFORM
  gate_trans_channels <- gate_channels[which(gate_channels %in% names(trans))]
  
  # TRANSFORMATIONS
  lapply(gate_trans_channels, function(z){
    # TRANSFORM COORDS - WATCH OUT INF COORDS
    if(inverse == FALSE){
      gate_coords[, z] <<- LAPPLY(gate_coords[,z], function(y){
        if(is.finite(y)){
          trans[[z]]$transform(y)
        }else{
          y
        }
      })
    }else{
      gate_coords[, z] <<- LAPPLY(gate_coords[,z], function(y){
        if(is.finite(y)){
          trans[[z]]$inverse(y)
        }else{
          y
        }
      })
    }
  })
  
  # UPDATE GATE
  x <- quadGate(filterId = gate_name,
                .gate = structure(as.list(gate_coords),
                                  names = colnames(gate_coords)))
  return(x)
  
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.filters <- function(x,
                                        trans = NULL,
                                        inverse = FALSE,
                                        ...){
  
  cyto_gate_transform(unlist(x),
                      trans = trans,
                      inverse = inverse)
  
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.list <- function(x,
                                     trans = NULL,
                                     inverse = FALSE,
                                     ...){
  
  # TRANSFORM GATES
  gates <- structure(
    lapply(x, function(z){
      cyto_gate_transform(
        z,
        trans = trans,
        inverse = inverse
      )
    }), names = names(x)
  )

  
  # RETURN UPDATED GATES
  return(gates)
  
}
