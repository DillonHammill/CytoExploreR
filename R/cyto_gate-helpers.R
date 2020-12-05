## CYTO_GATE_REMOVE ------------------------------------------------------------

#' Remove Gate(s) and Edit gatingTemplate csv File
#'
#' @param gs an object of class \code{GatingSet}.
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
#' @importFrom utils read.csv write.csv
#' @importFrom tools file_ext
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
cyto_gate_remove <- function(gs,
                             parent = NULL,
                             alias = NULL,
                             channels = NULL,
                             type = NULL,
                             gatingTemplate = NULL, ...) {

  # ALIAS MISSING
  if (is.null(alias)) {
    stop("Supply the name of the population to be removed to 'alias'.")
  }

  # ALIAS
  if (!all(alias %in% cyto_nodes(gs, path = "auto"))) {
    stop("Supplied alias does not exist in the GatingSet.")
  }

  # GATINGTEMPLATE MISSING
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }

  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to remove gate(s).")
  }

  # GATINGTEMPLATE FILE EXTENSION
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # READ IN GATINGTEMPLATE
  gt <- read.csv(gatingTemplate,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # PREPARE GATINGTEMPLATE ALIAS
  gt_alias <- lapply(gt$alias, function(z) {
    unlist(strsplit(as.character(z), ","))
  })

  # MULTIPLE ALIAS - REMOVE ALL ASSOCIATED NODES
  alias <- unique(LAPPLY(seq_len(length(alias)), function(z) {
    LAPPLY(gt_alias, function(y) {
      if (alias[z] == y) {
        y
      } else {
        NA
      }
    })
  }))
  alias <- alias[!is.na(alias)]

  # CHILDREN
  chldrn <- LAPPLY(
    alias,
    function(x) {
      if (!is.null(parent)) {
        pop <- paste0(parent, "/", x)
      } else {
        pop <- x
      }
      basename(gh_pop_get_descendants(gs[[1]], pop))
    }
  )
  chldrn <- unlist(chldrn, use.names = FALSE)
  chldrn <- unique(c(alias, chldrn))

  # REMOVE ROWS ALIAS == CHILDREN
  ind <- LAPPLY(gt_alias, function(z) {
    any(chldrn %in% z)
  })
  gt <- gt[!ind, ]

  # MESSAGE
  message(paste0(
    "Removing gate(s) from the GatingSet and ",
    gatingTemplate, "."
  ))

  # REMOVE NODES FROM GATINGSET
  for (i in seq_len(length(alias))) {
    if (alias[i] %in% cyto_nodes(gs, path = "auto")) {
      if (!is.null(parent)) {
        pop <- paste0(parent, "/", alias[i])
      } else {
        pop <- alias[i]
      }
      suppressMessages(gs_pop_remove(gs, pop))
    }
  }

  # UPDATE GATINGTEMPLATE
  write.csv(gt, gatingTemplate, row.names = FALSE)
}

## CYTO_GATE_RENAME ------------------------------------------------------------

#' Rename Gates
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param alias current names of the gates to be changed.
#' @param names new names to use for \code{alias}.
#' @param gatingTemplate name of the gatingTemplate csv file to be edited.
#'
#' @return update gate names in \code{GatingSet} or \code{GatingHierarchy} and
#'   update the gatingTemplate accordingly.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace gh_pop_set_name gs_pop_set_name
#' @importFrom utils read.csv write.csv
#' @importFrom tools file_ext
#'
#' @export
cyto_gate_rename <- function(x,
                             alias = NULL,
                             names = NULL,
                             gatingTemplate = NULL) {

  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }

  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to rename gate(s).")
  }

  # GATINGTEMPLATE FILE EXTENSION
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # ALIAS INVALID
  if (!all(alias %in% cyto_nodes(x, path = "auto"))) {
    stop(paste0("Supplied gate(s) do not exist in this ", class(x), "."))
  }

  # RENAME GATES IN GATINGHIERARCHY/GATINGSET
  mapply(function(alias, name) {
    if (is(x, "GatingHierarchy")) {
      gh_pop_set_name(x, alias, name)
    } else {
      gs_pop_set_name(x, alias, name)
    }
  }, alias, names)

  # READ IN GATINGTEMPLATE
  gt <- read.csv(gatingTemplate,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # UPDATE PARENTAL NAMES
  lapply(seq_along(alias), function(z) {
    # CHANGES DESCENDANT NODES AS WELL (eg CD4 T Cells & CD69+ CD4 T Cells)
    if (any(grepl(alias[z], gt$parent, fixed = TRUE))) {
      gt[grepl(alias[z], gt$parent, fixed = TRUE), "parent"] <<- names[z]
    }
  })

  # UPDATE ALIAS NAMES
  gt_alias <- lapply(gt$alias, function(z) {
    unlist(strsplit(as.character(z), ","))
  })
  gt_alias <- LAPPLY(gt_alias, function(z) {
    lapply(seq_len(length(alias)), function(y) {
      # CHANGES DESCENDANT NODES AS WELL (eg CD4 T Cells & CD69+ CD4 T Cells)
      if (any(grepl(alias[y], z, fixed = TRUE))) {
        z[grepl(alias[y], z, fixed = TRUE)] <<- names[y]
      }
    })
    # RE-COLLAPSE ALIAS
    z <- paste(z, collapse = ",")
    return(z)
  })
  gt[, "alias"] <- gt_alias

  # UPDATE GATING ARGUMENTS
  lapply(seq_len(nrow(gt)), function(z) {
    gating_args <- gt[z, "gating_args"]
    lapply(seq_along(alias), function(y) {
      # DESCENDANT NODES WILL ALSO CHANGE (eg CD4 T Cells & CD69+ CD4 T Cells)
      if (any(grepl(alias[y], gating_args, fixed = TRUE))) {
        gating_args <<- gsub(alias[y], names[y], gating_args)
      }
    })
    gt[z, "gating_args"] <<- gating_args
  })

  # SAVE UPDATED GATINGTEMPLATE
  write.csv(as.data.frame(gt), gatingTemplate, row.names = FALSE)
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
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_gate_copy <- function(x,
                           parent = NULL,
                           alias = NULL,
                           copy = NULL,
                           gatingTemplate = NULL){
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET
  if(!is(x, "GatingSet")){
    stop("'x' must be an object of class GatingSet.")
  }
  
  # NODES
  nodes_full <- cyto_nodes(x, path = "full")
  nodes_auto <- cyto_nodes(x, path = "auto")
  
  # PARENT
  if (is.null(parent)) {
    stop("Supply the name of the parent population.")
  } else if (!parent %in% c(nodes_full, nodes_auto)) {
    stop("Supplied parent does not exist in the GatingSet.")
  }
  
  # ALIAS
  if (is.null(copy)) {
    stop("Supply the name of the reference node to 'copy'.")
  } else if (!all(copy %in% c(nodes_full, nodes_auto))) {
    stop("Reference node does not exist in the GatingSet.")
  }
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to edit gate(s).")
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # UPDATE GATINGTEMPLATE ------------------------------------------------------
  
  # NEW ENTRIES
  gt_entries <- lapply(seq_along(alias), function(z){
    suppressMessages(
      gs_add_gating_method(
        gs = x,
        parent = parent,
        alias = alias[z],
        gating_method = "refGate",
        gating_args = copy[z]
      )
    )
  })
  
  # COMBINE GATINGTEMPLATE ENTRIES
  gt_entries <- do.call("rbind", 
                        gt_entries)
  
  # LOAD GATINGTEMPLATE
  gt <- read.csv(gatingTemplate, 
                 header = TRUE)
  
  # UPDATE GATINGTEMPLATE
  gt <- rbind(gt, 
              gt_entries)
  
  # WRITE UPDATED GATINGTEMPLATE
  write.csv(gt, 
            gatingTemplate, 
            row.names = FALSE)
  
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
                           gatingTemplate = NULL){
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET
  if(!is(x, "GatingSet")){
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
  
  # NODES
  nodes_full <- cyto_nodes(x, path = "full")
  nodes_auto <- cyto_nodes(x, path = "auto")
  
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
      if(!y %in% c(nodes_full, nodes_auto)){
        stop(paste(y, "does not exist or is not unique in the GatingSet."))
      }
    })
    return(p)
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
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to edit gate(s).")
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # CHECK GATINGTEMPLATE
  tryCatch(.cyto_gatingTemplate_check(parent,
                                      alias,
                                      gatingTemplate),
           error = function(e){
             stop(
               paste("Boolean populations already exist in gatingTemplate.",
                     "Remove them with cyto_gate_remove first to update them.")
             )
           })
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # GATINGTEMPLATE
  gt <- read.csv(gatingTemplate,
                 header = TRUE)
  
  # ENTRY
  gt_entry <- gt[1, , drop = FALSE]
  gt_entry[1, ] <- rep(NA, ncol(gt_entry))
  
  # ADD BOOLEAN ENTRIES
  gt_entries <- lapply(seq_along(alias), function(z){
    gt_pop <- gt_entry
    gt_pop[, "alias"] <- alias[z]
    gt_pop[, "pop"] <- "+"
    gt_pop[, "parent"] <- parent[z]
    gt_pop[, "gating_method"] <- "boolGate"
    gt_pop[, "gating_args"] <- logic[z]
    return(gt_pop)
  })
  gt_entries <- do.call("rbind", gt_entries)
  
  # UPDATE GATINGTEMPLATE
  gt <- rbind(gt, gt_entries)
  write.csv(gt,
            gatingTemplate, 
            row.names = FALSE)
  
  # CREATE BOOLEAN GATES
  gates <- lapply(seq_along(alias), function(z){
    call <- substitute(booleanFilter(v, filterId = alias[z]), 
                       list(v = as.symbol(logic[z])))
    eval(call)
  })
  names(gates) <- alias
  
  # APPLY BOOLEAN GATES
  lapply(seq_along(gates), function(z){
    gs_pop_add(x, 
               gate = gates[[z]])
  })
  
  # RETURN GATINGSET
  return(x)
  
}

## CYTO_GATE_EXTRACT -----------------------------------------------------------

#' Extract Saved Gate(s) from gatingTemplate.
#'
#' @param parent name of the parental population.
#' @param alias name of the population for which the gate must be extracted.
#' @param gatingTemplate name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the gate(s) are saved.
#' @param ... not used.
#'
#' @importFrom openCyto gt_gating gt_get_nodes gt_get_gate
#' @importFrom flowCore filters parameters<-
#' @importFrom tools file_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Bypass working directory check for external files
#' options("CytoExploreR_wd_check" = FALSE)
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
#' # Reset CytoExploreR_wd_check to default
#' options("CytoExploreR_wd_check" = TRUE)
#' @export
cyto_gate_extract <- function(parent,
                              alias,
                              gatingTemplate = NULL, ...) {

  # MISSING PARENT
  if (missing(parent)) {
    stop("Supply the name of the parent population.")
  }

  # MISSING ALIAS
  if (missing(alias)) {
    stop("Supply the name(s) of the alias to extract.")
  }

  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }

  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to extract gate(s).")
  }

  # GATINGTEMPLATE FILE EXTENSION
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # READ IN GATINGTEMPLATE
  gt <- suppressMessages(gatingTemplate(gatingTemplate))

  # EXTRACT NODES - GATINGTEMPLATE
  nds <- gt_get_nodes(gt, only.names = TRUE)

  # PARENTAL NODE
  parent <- names(nds[match(parent, nds)])

  # EXTRACT GATES GIVEN PARENTAL & CHILD NODES
  gates <- lapply(alias, function(x) {
    # ALIAS NODE
    ind <- LAPPLY(seq_len(length(nds)), function(z) {
      if (x %in% nds[[z]]) {
        z
      } else {
        NA
      }
    })
    ind <- ind[!is.na(ind)][1]
    alias <- names(nds[ind])
    gm <- gt_get_gate(gt, parent, alias)
    gate <- eval(parameters(gm)$gate)
    return(gate)
  })

  return(gates)
}

## CYTO_GATE_EDIT --------------------------------------------------------------

#' Edit Existing Gate(s).
#'
#' @param x an object of class \code{GatingSet}.
#' @param parent name of the parental population.
#' @param alias name(s) of the gate to edit (e.g. "Single Cells").
#' @param channels name(s) of the channel(s) used to construct the gate(s). This
#'   argument is not necessary but is included to allow conversion of
#'   \code{cyto_gate_draw} code to \code{cyto_gate_remove} code by simply changing
#'   \code{"draw"} to \code{"remove"}.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   \code{types} are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported \code{gate_types} are \code{polygon, rectangle, ellipse,
#'   threshold, boundary, interval, quadrant and web} which can be abbreviated
#'   as upper or lower case first letters as well. Default \code{type} is
#'   \code{"polygon"}.
#' @param gatingTemplate name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the gate is saved.
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
#'   plotting.
#' @param display fraction or number of events to display in the plot during the
#'   gating process, set to 25 000 events by default.
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
#' @param ... additional arguments for \code{\link{cyto_plot.flowFrame}}.
#'
#' @return an object of class \code{GatingSet} with edited gate applied, as well
#'   as gatingTemplate file with edited gate saved.
#'
#' @importFrom flowWorkspace gh_pop_get_gate gs_pop_set_gate recompute
#' @importFrom flowCore parameters filterList
#' @importFrom openCyto gatingTemplate CytoExploreR_.argDeparser
#' @importFrom data.table as.data.table fread :=
#' @importFrom methods as
#' @importFrom utils select.list write.csv
#' @importFrom tools file_ext
#' @importFrom purrr transpose
#' @importFrom methods is
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

  # PARENT
  if (is.null(parent)) {
    stop("Supply the name of the parent population.")
  } else if (!parent %in% cyto_nodes(x, path = "auto")) {
    stop("Supplied parent does not exist in the GatingSet.")
  }

  # ALIAS
  if (is.null(alias)) {
    stop("Supply the name(s) of the gates to edit to 'alias'.")
  } else if (!all(alias %in% cyto_nodes(x, path = "auto"))) {
    stop("Supplied alias does not exist in the GatingSet.")
  }

  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }

  # GATINGTEMPLATE STILL MISSING
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate to edit gate(s).")
  }

  # GATINGTEMPLATE FILE EXTENSION
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # ASSIGN X TO GS
  gs <- x

  # TRANSFORMATIONS
  axes_trans <- cyto_transformer_extract(gs)

  # NODES
  nds <- cyto_nodes(gs, path = "auto")

  # EXTRACT INFORMATION FROM GATINGTEMPLATE FILE -------------------------------

  # READ IN GATINGTEMPLATE
  gtf <- read.csv(gatingTemplate, header = TRUE, stringsAsFactors = FALSE)

  # RESTRICT BY PARENT
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
    stop("Use cyto_gatingTemplate_edit to modify refGate and boolGate entries.")
  }

  # EXTRACT EXISTING GATES FROM GATINGSET DIRECTLY -----------------------------

  # GATINGSET LIST - NEW GROUPING - EXTRACT EXISTING GATES FROM EACH GROUP
  gs_list <- cyto_group_by(gs, 
                           group_by = group_by)

  # GROUPS
  N <- length(gs_list)
  GRPS <- names(gs_list)

  # GATES FORMATTED FOR GATINGSET (NO QUADGATES) - EXISTING GATES
  gates_gs <- lapply(gs_list, function(gs) {
    gates <- lapply(alias, function(z) {
      gh_pop_get_gate(gs[[1]], paste(parent, z, sep = "/"))
    })
    names(gates) <- alias
    return(gates)
  })
  names(gates_gs) <- names(gs_list)
  
  # GATES FORMATTED FOR GATINGTEMPLATE (QUADGATES) - EXISTING GATES
  gates_gT <- lapply(gates_gs, function(z){
    if(all(LAPPLY(z, function(y){
      is(y, "rectangleGate") & any(grepl("quad", names(attributes(y))))
    }))){
      return(.cyto_gate_quad_convert(z, channels))
    }else{
      return(z)
    }
  })
  
  # PREPARE TYPE & ALIAS -------------------------------------------------------

  # GET GATE TYPE FROM EXISTING GATES
  if (is.null(type)) {
    type <- cyto_gate_type(gates_gs[[1]])
  }

  # TYPE
  type <- .cyto_gate_type(type, channels, alias, negate)

  # ALIAS - LIST PER GATE TYPE
  alias <- .cyto_alias(alias, type)

  # PREPARE SAMPLES & OVERLAYS -------------------------------------------------

  # EXTRACT PARENT POPULATION
  fs <- cyto_extract(x, parent)

  # GROUPING & MERGING
  fr_list <- cyto_merge_by(fs,
    merge_by = group_by,
    select = select
  )
  names(fr_list) <- GRPS

  # PREPARE OVERLAY ------------------------------------------------------------

  # Organise overlays - list of flowFrame lists of length(fr_list)
  if (!.all_na(overlay)) {
    # OVERLAY - POPULATION NAMES
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
          # EXTRACT POPULATIONS
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
  }else{
    overlay <- rep(list(list(NA)), N)
  }

  # COMBINE FR_LIST & OVERLAY
  fr_list <- lapply(seq_along(fr_list), function(z) {
    return(c(fr_list[z], overlay[[z]]))
  })
  names(fr_list) <- GRPS

  # SAMPLING - SAME SAMPLING PER LAYER
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
  
  # SELECT GROUP(S) TO EDIT ----------------------------------------------------

  # MENU - SELECT GROUPS TO EDIT
  if (group_by[1] == "all") {
    grps <- "all"
  # SKIP MENU IF GROUPING ADDED HERE
  }else if(is.na(gtf_groupBy) & group_by[1] != "all"){
    grps <- names(fr_list)
  # SKIP MENU IF GROUPING HAS CHANGED  
  }else if(!is.na(gtf_groupBy) & 
           !identical(group_by, unlist(strsplit(gtf_groupBy, ":")))){  
    grps <- names(fr_list)
  # GROUPING REMAINS THE SAME
  } else {
    message("Select the group(s) to edit:")
    grps <- data_editor(data.frame("group" = names(fr_list),
                                   "select" = NA,
                                   stringsAsFactors = FALSE),
                        title = "Group Selection",
                        type = "menu")
    grps <- grps[, "group"][which(grps[, "select"] == 1)]
    # NA TO CHARCTERS
    if(any(is.na(grps))){
      grps[is.na(grps)] <- "NA"
    }
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)

  # GROUPING
  if (group_by[1] == "all") {
    pd$groupby <- rep("all", nrow(pd))
  } else {
    pd$groupby <- do.call(paste, pd[, group_by, drop = FALSE])
  }

  # GATE EDITING ---------------------------------------------------------------

  # PLOT/EDIT/SAVE EACH GROUP - ORGANISE FOR GS & GT
  lapply(match(grps, names(fr_list)), function(y) {

    # CYTO_PLOT - PROPER TITLES
    if(plot == TRUE){
      
      # EXISTING GATES TO PLOT
      gate <- gates_gT[[y]]
      
      # PARENT TITLE
      if (parent == "root") {
        prnt <- "All Events"
      } else {
        prnt <- parent
      }
      
      # TITLE
      if(!is.na(GRPS[y]) & GRPS[y] == "all"){
        grp <- "Combined Events"
      }else{
        grp <- GRPS[y]
      }
      title <- paste(grp, "\n", prnt)

      # CYTO_PLOT
      cyto_plot(fr_list[[y]][[1]],
        channels = channels,
        overlay = fr_list[[y]][seq_along(fr_list[[y]])[-1]],
        display = 1,
        legend = FALSE,
        gate = gate,
        gate_line_col = "magenta",
        gate_line_alpha = 0.8,
        axes_trans = axes_trans,
        label = FALSE,
        title = title,
        gate_line_width = 2.5,
        popup = popup,
        axes_limits = axes_limits, ...
      )
      
    }
    
    # 2D Interval gates require axis argument
    if ("interval" %in% type) {
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
    
    # DRAW NEW GATES - FILTERS OBJECT OF LENGTH ALIAS (NEGATE LABEL)
    gate_new <- cyto_gate_draw(fr_list[[y]][[1]],
      alias = unlist(alias),
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
    )
    
    # NAME GATES
    names(gate_new) <- LAPPLY(alias, function(z) {
      paste(z, collapse = ",")
    })

    # MODIFY EXISTING GATES - GATINGSET
    gates_gs[[y]] <<- gate_new

    # MODIFY EXISTING GATES - GATINGTEMPLATE
    gates_gT[[y]] <<- gate_new
  })
  
  # PREPARE GATES FOR GATINGTEMPLATE - LIST OF FILTERS OF LENGTH ALIAS
  gates_gT_transposed <- transpose(gates_gT)
  
  # ADD GATES TO FILTERS LISTS
  gates_gT_transposed <- lapply(seq_along(gates_gT_transposed), 
                               function(z){
    gates <- lapply(seq_along(gates_gT_transposed[[z]]), function(y){
      if(!is(gates_gT_transposed[[z]][[y]], "quadGate")){
        filters(gates_gT_transposed[[z]][y])
      }else{
        gates_gT_transposed[[z]][[y]]
      }
    })
    names(gates) <- GRPS
    return(gates)
  })
  names(gates_gT_transposed) <- alias
  
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
  gt <- data.table::fread(gatingTemplate)

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

  # SAVE UPDATED GATINGTEMPLATE
  write.csv(as.data.frame(gt), gatingTemplate, row.names = FALSE)

  # CONVERT ALIAS BACK TO VECTOR
  alias <- as.character(unlist(alias))

  # CONVERT QUAD TO RECTANGLE FOR GATINGSET SAVING (EASIEST WAY)
  gates_gs <- lapply(gates_gs, function(z){
    if(is(z[[1]], "quadGate")){
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
    if (is(gates)[1] %in% c("filters", "list")) {
      gates <- gates[[1]]
    }
    # ELLIPSE
    if (is(gates, "ellipsoidGate")) {
      types <- "ellipse"
      # RECTANGLE/BOUNDARY/INTERVAL/THRESHOLD
    } else if (is(gates, "rectangleGate")) {
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
    } else if (is(gates, "polygonGate")) {
      types <- "polygon"
      # QUADRANT
    } else if (is(gates, "quadGate")) {
      types <- "quadrant"
    }
    # Multiple gates supplied
  } else if (length(gates) > 1) {

    # Get classes of gates
    classes <- LAPPLY(gates, function(x) {
      class(x)
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
  
  # Invalid gate object
  if (!is(x)[1] %in% c(
    "list",
    "filters",
    "rectangleGate",
    "polygonGate",
    "ellipsoidGate",
    "quadGate"
  )) {
    stop(paste(
      "'x' should contain either filters, rectangleGate, polygonGate,",
      "ellipsoidGate or quadGate objects."
    ))
  }
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
  x <- lapply(x, function(z) {
    cyto_gate_convert(z, channels)
  })
  
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
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_gate_convert
cyto_gate_prepare <- function(x,
                              channels = NULL) {

  # LIST OF GATES --------------------------------------------------------------

  # PREPARE GATE LIST
  if (is(x)[1] == "list") {
    if (all(LAPPLY(x, "is") %in% c(
      "rectangleGate",
      "polygonGate",
      "ellipsoidGate",
      "quadGate",
      "filters"
    ))) {
      x <- unlist(x)
    }
  } else if (is(x)[1] == "filters") {
    x <- unlist(x)
  } else if (is(x)[1] %in% c(
    "rectangleGate",
    "polygonGate",
    "ellipsoidGate",
    "quadGate",
    "filters"
  )) {
    x <- list(x)
  }

  # UNIQUE GATE LIST
  x <- unique(x)

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
  if(is.null(trans) | !is(trans, "transformerList")){
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
  if(is.null(trans) | !is(trans, "transformerList")){
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
  
  # TRANSFORMERS
  if(is.null(trans) | !is(trans, "transformerList")){
    stop("Supply a list of transformers to transform the gate co-ordinates.")
  }
  
  # GATE FILTERID
  gate_name <- x@filterId  
  
  # GATE CHANNELS
  gate_channels <- parameters(x)  
  
  # GATE COORDS (POLYGON)
  gate_coords <- .cyto_gate_coords(list(x), channels = gate_channels)
  
  # GATE CENTER (ELLIPSE)
  gate_center <- x@mean  
  
  # ADD GATE CENTER TO GATE COORDS
  gate_coords <- rbind(gate_coords, gate_center)
  
  # GATE COORDS (ELLIPSE)
  gate_dist <- as.matrix(dist(round(gate_coords, 2)))
  gate_dist <- gate_dist[, ncol(gate_dist)]
  gate_dist <- round(gate_dist[gate_dist != 0], 2)
  names(gate_dist) <- seq_along(gate_dist) # ORIGINAL ORDER
  gate_dist_sort <- sort(gate_dist)
  
  # MINOR AXIS - WATCH OUT FOR DUPLICATES
  ind <- as.numeric(names(gate_dist_sort[c(1,3)]))
  gate_min <- gate_coords[ind, ] # closest 2 points
  
  # MAJOR AXIS - WATCH OUT FOR DUPLICATES
  ind <- as.numeric(names(gate_dist_sort[c(length(gate_dist_sort),
                                           length(gate_dist_sort) - 3)]))
  gate_max <- gate_coords[ind, ] # closest 2 points
  gate_coords <- rbind(gate_min, gate_max, gate_center)
  
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
  
  # RADIUS MINOR AXIS
  b <- dist(gate_coords[c(1,2), ])/2  
  
  # RADIUS MAJOR AXIS
  a <- dist(gate_coords[c(3,4), ])/2
  
  # ANGLE - HORIZONTAL THROUGH CENTER TO MAX POINT
  max_point <- gate_coords[which.max(gate_coords[, 2]), ]
  if (max_point[1] > gate_center[1]) { # angle < pi/2
    mj.pt.ct <- cbind(max_point[1], gate_center[2])
    colnames(mj.pt.ct) <- c("x", "y")
    adj <- stats::dist(rbind(gate_center, mj.pt.ct))
    angle <- acos(adj / a)
  } else if (max_point[1] <= gate_center[1]) { # angle >= pi/2
    mj.pt.ct <- cbind(gate_center[1], max_point[2])
    colnames(mj.pt.ct) <- c("x", "y")
    opp <- stats::dist(as.matrix(rbind(max_point, mj.pt.ct)))
    angle <- pi / 2 + asin(opp / a)
  }
  
  # COVARIANCE MATRIX
  cinv <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
  cinv[1, 1] <- (((cos(angle) * cos(angle)) /
                    (a^2)) + ((sin(angle) * sin(angle)) / (b^2)))
  cinv[2, 1] <- sin(angle) * cos(angle) * ((1 / (a^2)) - (1 / (b^2)))
  cinv[1, 2] <- cinv[2, 1]
  cinv[2, 2] <- (((sin(angle) * sin(angle)) / (a^2)) +
                   ((cos(angle) * cos(angle)) / (b^2)))
  cvm <- solve(cinv)
  dimnames(cvm) <- list(gate_channels, gate_channels)
  
  # UPDATE GATE
  x <- ellipsoidGate(filterId = gate_name,
                     .gate = cvm,
                     mean = gate_center)
  return(x)
  
}

#' @rdname cyto_gate_transform
#' @export
cyto_gate_transform.quadGate <- function(x, 
                                         trans = NULL,
                                         inverse = FALSE,
                                         ...){
  
  # TRANSFORMERS
  if(is.null(trans) | !is(trans, "transformerList")){
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
  nms <- names(x)
  gates <- lapply(x, function(z){
    cyto_gate_transform(z,
                        trans = trans,
                        inverse = inverse)
  })
  if(!is.null(nms)){
    names(gates) <- nms
  }
  
  # RETURN UPDATED GATES
  return(gates)
  
}
