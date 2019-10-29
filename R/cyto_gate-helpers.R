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
#' @importFrom flowWorkspace getDescendants Rm getNodes
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
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gating(gt, gs)
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
    stop("Please supply the name of the population to be removed.")
  }
  
  # ALIAS
  if (!all(alias %in% basename(getNodes(gs)))) {
    stop("Supplied alias does not exist in the GatingSet.")
  }
  
  # GATINGTEMPLATE MISSING
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to remove gate(s).")
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # READ IN GATINGTEMPLATE
  gt <- read.csv(gatingTemplate, header = TRUE) 
  
  # PREPARE GATINGTEMPLATE ALIAS
  gt_alias <- lapply(gt$alias, function(z){
    unlist(strsplit(as.character(z), ","))
  })
  
  # QUADGATES - REMOVE ALL NODES
  alias <- unique(LAPPLY(seq_len(length(alias)), function(z){
    LAPPLY(gt_alias, function(y){
      if(any(grepl(alias[z], y, fixed = TRUE))){
        y
      }
    })
  }))
  
  # CHILDREN
  chldrn <- LAPPLY(
    alias,
    function(x) basename(getDescendants(gs[[1]], paste0(parent,"/",x)))
  )
  chldrn <- unlist(chldrn, use.names = FALSE)
  chldrn <- unique(c(alias, chldrn))
  
  # REMOVE ROWS ALIAS == CHILDREN
  ind <- LAPPLY(gt_alias, function(z){
    LAPPLY(chldrn, function(y){
      any(grepl(y, z), fixed = TRUE)
    })
  })
  gt <- gt[!ind, ]
  
  # MESSAGE
  message(paste0("Removing gate(s) from the GatingSet and ",
                gatingTemplate,"."))
  
  # REMOVE NODES FROM GATINGSET
  for (i in seq_len(length(alias))) {
    if (alias[i] %in% basename(getNodes(gs))) {
      suppressMessages(Rm(paste0(parent,"/", alias[i]), gs))
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
#' @importFrom flowWorkspace getNodes setNode
#' @importFrom utils read.csv write.csv
#' @importFrom tools file_ext
#'
#' @export
cyto_gate_rename <- function(x,
                             alias = NULL,
                             names = NULL,
                             gatingTemplate = NULL){
  
  # MMISSING GATINGTEMPLATE
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to rename gate(s).")
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # ALIAS INVALID
  if(!all(alias %in% basename(getNodes(x)))){
    stop(paste0("Supplied gate(s) do not exist in this ", class(x), "."))
  }
  
  # RENAME GATES IN GATINGHIERARCHY/GATINGSET
  mapply(function(alias, name){
    setNode(x, alias, name)
  }, alias, names)
  
  # READ IN GATINGTEMPLATE
  gt <- read.csv(gatingTemplate, 
                 header = TRUE, 
                 stringsAsFactors = FALSE)
  
  # UPDATE PARENTAL NAMES
  if(any(alias %in% gt$parent)){
    gt[gt$parent %in% alias,"parent"] <- names
  }
  
  # UPDATE ALIAS NAMES
  gt_alias <- lapply(gt$alias, function(z){
    unlist(strsplit(as.character(z), ","))
  })
  gt_alias <- LAPPLY(gt_alias, function(z){
    lapply(seq_len(length(alias)), function(y){
      if(alias[y] %in% z){
        z[z == alias[y]] <<- names[y]
      }
    })
    # RE-COLLAPSE ALIAS
    z <- paste(z, collapse = ",")
    return(z)
  })
  gt[, "alias"] <- gt_alias
  
  # SAVE UPDATED GATINGTEMPLATE
  write.csv(as.data.frame(gt), gatingTemplate, row.names = FALSE)
  
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
#' @importFrom flowWorkspace getGate getNodes
#' @importFrom openCyto gating
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
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#' 
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gating(gt, gs)
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
    stop("Please supply the name of the parent population.")
  }
  
  # MISSING ALIAS
  if (missing(alias)) {
    stop("Please supply the name(s) of the alias to extract.")
  }
  
  # MISSING GATINGTEMPLATE
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to extract gate(s).")
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # READ IN GATINGTEMPLATE
  gt <- suppressMessages(gatingTemplate(gatingTemplate))
  
  # EXTRACT NODES - GATINGTEMPLATE
  nds <- getNodes(gt, only.names = TRUE)
  
  # PARENTAL NODE
  parent <- names(nds[match(parent, nds)])
  
  # EXTRACT GATES GIVEN PARENTAL & CHILD NODES
  gates <- lapply(alias, function(x) {
    # ALIAS NODE
    ind <- LAPPLY(seq_len(length(nds)), function(z){
      if(x %in% nds[[z]]){
        z
      }else{
        NA
      }
    })
    ind <- ind[!is.na(ind)][1]
    alias <- names(nds[ind])
    gm <- getGate(gt, parent, alias)
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
#' @param select vector containing the indicies of samples within gs to use for
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
#' @param ... additional arguments for \code{\link{cyto_plot.flowFrame}}.
#'
#' @return an object of class \code{GatingSet} with edited gate applied, as well
#'   as gatingTemplate file with edited gate saved.
#'
#' @importFrom flowWorkspace getData getTransformations GatingSet getGate
#'   setGate recompute pData
#' @importFrom flowCore parameters filterList
#' @importFrom openCyto gatingTemplate CytoExploreR_.argDeparser
#' @importFrom data.table as.data.table fread :=
#' @importFrom methods as
#' @importFrom utils select.list write.csv
#' @importFrom tools file_ext
#' @importFrom magrittr %>%
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
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gating(gt, gs)
#'
#' # Edit CD4 T Cells Gate - replace gatingTemplate name
#' gate_edit(gs,
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
                           label = TRUE, ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # PARENT
  if (is.null(parent)) {
    stop("Please supply the name of the parent population.")
  } else if (!parent %in% cyto_nodes(x, path = "auto")) {
    stop("Supplied parent does not exist in the GatingSet.")
  }
  
  # ALIAS
  if (is.null(alias)) {
    stop("Please supply the name(s) of the gates to edit to 'alias'.")
  } else if (!all(alias %in% cyto_nodes(x, path = "auto"))) {
    stop("Supplied alias does not exist in the GatingSet.")
  }
  
  # MISSING GATINGTEMPLATE
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # GATINGTEMPLATE STILL MISSING
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to edit gate(s).")
  }
  
  # GATINGTEMPLATE FILE EXTENSION
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # ASSIGN X TO GS
  gs <- x
  
  # TRANSFORMATIONS
  axes_trans <- gs@transformation
  if(length(axes_trans) != 0){
    axes_trans <- axes_trans[[1]]
  }else{
    axes_trans <- NA
  }  
  
  # EXTRACT INFORMATION FROM GATINGTEMPLATE FILE -------------------------------
  
  # READ IN GATINGTEMPLATE
  gtf <- read.csv(gatingTemplate, header = TRUE, stringsAsFactors = FALSE)
  
  # RESTRICT BY PARENT
  gtf_parent <- gtf[gtf$parent == parent, ]
  
  # RESTRICT BY PARENT & ALIAS
  gtf_alias <- lapply(gtf_parent$alias, function(z){
    unlist(strsplit(as.character(z), ","))
  })
  ind <- LAPPLY(gtf_alias, function(z){
    if(any(alias %in% z)){
      TRUE
    }else{
      FALSE
    }
  })
  # UPDATE gtf_alias
  gtf_alias <- gtf_alias[ind]

  # PULL OUT RELEVENT GATINTEMPLATE ENTRIES
  gtf_chunk <- gtf_parent[ind, ]

  # EXTRACT CHANNELS
  if(is.null(channels)){
    channels <- unique(
      unlist(
        strsplit(
          as.character(gtf_chunk[gtf_chunk$alias %in% alias, "dims"]), ",",
          fixed = TRUE
        )
      )
    )  
  }
  
  # EXTRACT GROUPBY
  gtf_groupBy <- unique(gtf_chunk[, "groupBy"])
  
  # EXTRACT GATING METHOD
  gtf_gating_method <- unique(as.character(gtf_chunk[, "gating_method"]))
  
  # UNSUPPORTED GATING METHOD
  if(any(gtf_gating_method %in% c("refGate", "boolGate"))){
    stop("Use cyto_gatingTemplate_edit to modify refGate and boolGate entries.")
  }
  
  # EXTRACT EXISTING GATES FROM GATINGSET DIRECTLY -----------------------------
  
  # GATINGSET LIST - NEW GROUPING - EXTRACT EXISTING GATES FROM EACH GROUP
  gs_list <- cyto_group_by(gs, group_by = group_by)
  
  # GROUPS
  N <- length(gs_list)  
  
  # GATES FORMATTED FOR GATINGSET (NO QUADGATES) - EXISTING GATES
  gates_gs <- lapply(gs_list, function(gs){
    gates <- lapply(alias, function(z){
      getGate(gs[[1]], paste(parent, z, sep = "/"))
    })
    names(gates) <- alias
    return(gates)
  })
  names(gates_gs) <- names(gs_list)
  
  # GATES FORMATTED FOR GATINGTEMPLATE (QUADGATES) - EXISTING GATES
  gates_gT <- lapply(gates_gs, function(z){
    if(cyto_gate_type(z) == "quadrant"){
      z <- list(.cyto_gate_quad_convert(z, channels))
    }
    return(z)
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
        tryCatch(cyto_select(z, select), error = function(e) {
          z
        })
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
  names(fr_list) <- names(gs_list)

  # PREPARE OVERLAY ------------------------------------------------------------
  
  # Organise overlays - list of flowFrame lists of length(fr_list)
  if (!.all_na(overlay)) {
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
    if (inherits(overlay, "flowFrame")) {
      # Always show all events
      overlay <- rep(list(list(overlay)), N)
      # flowSet to lists of flowFrame lists
    } else if (inherits(overlay, "flowSet")) {
      # GROUP VARIABLES
      if (group_by[1] != "all") {
        # GROUPING
        overlay <- cyto_group_by(overlay, group_by)
        # LIST OF FLOWFRAMES
        overlay <- lapply(overlay, function(z) {
          # SELECT
          if (!is.null(select)) {
            z <- tryCatch(cyto_select(z, select), error = function(e) {
              z
            })
          }
          # CONVERT
          z <- cyto_convert(z, "flowFrame")
          return(z)
        })
        # GROUP ALL
      } else {
        # SELECT
        if (!is.null(select)) {
          overlay <- cyto_select(overlay, select)
        }
        # FLOWFRAME
        overlay <- cyto_convert(overlay, "flowFrame")
        # FLOWFRAME LIST
        overlay <- list(overlay)
      }
      # FLOWFRAME LIST TO LIST OF FLOWFRAME LISTS
      overlay <- lapply(overlay, function(z) {
        list(z)
      })
      # OVERLAY - LIST OF FLOWFRAMES OR FLOWSETS
    } else if (inherits(overlay, "list")) {
      # LIST OF FLOWFRAMES - REPEAT FR_LIST TIMES
      if (all(LAPPLY(overlay, function(z) {
        inherits(z, "flowFrame")
      }))) {
        overlay <- rep(list(overlay), N)
        # LIST FLOWFRAME LISTS OF LENGTH FR_LIST
      } else if (all(LAPPLY(unlist(overlay), function(z) {
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
      } else if (all(LAPPLY(overlay, function(z) {
        inherits(z, "flowSet")
      }))) {
        # GROUP & MERGE EACH FLOWSET
        overlay <- lapply(overlay, function(z) {
          # GROUP VARIABLES
          if (group_by[1] != "all") {
            # GROUPING
            x <- cyto_group_by(z, group_by)
            # Coercion and sampling
            x <- lapply(x, function(y) {
              # SELECT
              if (!is.null(select)) {
                y <- tryCatch(cyto_select(y, select), error = function(e) {
                  y
                })
              }
              # CONVERT
              y <- cyto_convert(y, "flowFrame")
              return(y)
            })
            # GROUP ALL
          } else {
            # SELECT
            if (!is.null(select)) {
              z <- tryCatch(cyto_select(z, select), error = function(e) {
                z
              })
            }
            # CONVERT
            z <- cyto_convert(z, "flowFrame")
            return(z)
          }
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
  }
  
  # COMBINE FR_LIST & OVERLAY
  if(!.all_na(overlay)){
    fr_list <- lapply(N, function(z){
      c(fr_list[z], overlay[[z]])
    })
    names(fr_list) <- names(gs_list)
  }
  
  # SELECT GROUP(S) TO EDIT ----------------------------------------------------
  
  # MENU - SELECT GROUPS TO EDIT
  if(group_by[1] == "all"){
    grps <- "all"
  }else{
    # GROUPED SAMPLES
    grps <- select.list(names(fr_list),
                        multiple = TRUE,
                        graphics = TRUE,
                        title = "Select the group(s) to edit:")
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # GROUPING
  if(group_by[1] == "all"){
    pd$groupby <- rep("all", nrow(pd))
  }else{
    pd$groupby <- do.call(paste, pd[, group_by, drop = FALSE])
  }
  
  # GATE EDITING ---------------------------------------------------------------
  
  # PLOT/EDIT/SAVE EACH GROUP - COORGANISE FOR GS & GT
  lapply(match(grps, names(fr_list)), function(y) {
    
    # SAMPLING
    FR_LIST <- cyto_sample(fr_list[y], display = display, seed = 56)
    
    # EXISTING GATES TO PLOT
    gate <- gates_gT[[y]]
    
    # Title
    if(group_by[1] == "all"){
      title <- paste("Combined Events" ,"\n", parent)
    }else{
      title <- paste(names(fr_list)[y], "\n", parent)
    }
    
    # Call to cyto_plot - careful about overlay
    if(!length(FR_LIST) == 1){
      cyto_plot(FR_LIST[[1]],
                channels = channels,
                overlay = FR_LIST[seq(2, length(FR_LIST))],
                legend = FALSE,
                gate = gate,
                gate_line_col = "magenta",
                gate_line_alpha = 0.8,
                axes_trans = axes_trans,
                label = FALSE,
                title = title,
                gate_line_width = 2.5,
                popup = TRUE, ...
      )
    }else{
      cyto_plot(FR_LIST[[1]],
                channels = channels,
                overlay = NA,
                legend = FALSE,
                gate = gate,
                gate_line_col = "magenta",
                gate_line_alpha = 0.8,
                axes_trans = axes_trans,
                label = FALSE,
                title = title,
                gate_line_width = 2.5,
                popup = TRUE, ...
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
    gate_new <- cyto_gate_draw(FR_LIST[[1]],
                                alias = unlist(alias),
                                channels = channels,
                                type = type,
                                negate = negate,
                                axis = axis,
                                plot = FALSE
    )
    
    # REMOVE FILTERS LAYER
    gate_new <- unlist(gate_new)
    names(gate_new) <- LAPPLY(alias, function(z){paste(z, collapse = ",")})
    
    # PREPARE GATES FOR GATINGSET (QUADGATE -> REACTANGLES)
    gate_prep <- lapply(gate_new, function(z){
      if(is(z, "quadGate")){
        z <- .cyto_gate_quad_convert(z, channels)
      }
      return(z)
    })
    gate_prep <- unlist(gate_prep)
    names(gate_prep) <- unlist(alias)
    
    # MODIFY EXISTING GATES - GATINGSET
    gates_gs[[y]] <<- gate_prep
    
    # MODIFY EXISTING GATES - GATINGTEMPLATE
    gates_gT[[y]] <<- gate_new
    
  })
  
  # PREPARE GATES FOR GATINGTEMPLATE - LIST OF FILTERS OF LENGTH ALIAS
  gates_gT_transposed <- gates_gT %>% transpose()
  
  # ADD GATES TO FILTERS LISTS
  gates_gT_transposed <- lapply(gates_gT_transposed, function(z){
    if(!any(LAPPLY(z, function(y){is(y, "quadGate")}))){
      z <- filters(z)
    }
    return(z)
  })
  
  # DATA.TABLE FRIENDLY NAMES
  prnt <- parent
  als <- alias
  gtmd <- "cyto_gate_draw"
  ppmd <- "pp_cyto_gate_draw"
  
  # REMOVE NEGATED POPULATIONS FROM ALS (NEGATED ENTRY CANNOT BE MODIFIED)
  als <- lapply(als, function(z){
      paste(z, collapse = ",")
  })
  
  # REMOVE NEGATED REFERENCE
  if(negate == TRUE){
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
  
  # GROUP_BY
  if(group_by[1] == "all"){
    group_by <- NA
  }else{
    group_by <- paste(group_by, collapse = ":")
  }
  
  # MODIFY GATINGTEMPLATE - GATED POPULATIONS ONLY (IGNORE NEGATED ENTRIES)
  for (i in seq_len(length(als))) {
    gt[parent == prnt & alias == als[i], gating_method := gtmd]
    gt[parent == prnt & alias == als[i], 
       gating_args := CytoExploreR_.argDeparser(list(
      gate = gates_gT_transposed[[i]]
    ))]
    gt[parent == prnt & alias == als[i], collapseDataForGating := TRUE]
    gt[parent == prnt & alias == als[i], groupBy := group_by]
    gt[parent == prnt & alias == als[i], preprocessing_method := ppmd]
    gt[parent == prnt & alias == als[i], preprocessing_args := as.logical(NA)]
  }

  # SAVE UPDATED GATINGTEMPLATE
  write.csv(as.data.frame(gt), gatingTemplate, row.names = FALSE)
  
  # CONVERT ALIAS BACK TO VECTOR
  alias <- as.character(unlist(alias))

  # APPLY NEW GATES TO GATINGSET (INCLUDES NEGATED ALIAS)
  lapply(grps, function(z){
    # MATCHING GROUPS
    ind <- which(pd$groupby == z)
    lapply(seq_along(alias), function(y){
      # SET GATES - GATED POPULATIONS ONLY
      if(y < length(alias) | (negate == FALSE & y == length(alias))){
        # SET UPDATED GATES- IGNORE NEGATED GATE
        # REPEAT GATES FOR EACH SAMPLE IN GROUP
        gates_update <- rep(list(gates_gs[[z]][[alias[y]]]), 
                            length(gs[ind]))
        names(gates_update) <- cyto_names(gs[ind])
        suppressMessages(setGate(gs[ind], 
                                 paste(parent, alias[y], sep = "/"),
                                 gates_update))
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
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#' 
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#' 
#' # Get gate type used for T Cells gate
#' cyto_gate_type(cyto_extract(gs, "Cells")[[1]])
#' cyto_gate_type(cyto_extract(gs, "T Cells")[[1]])
#' cyto_gate_type(cyto_extract(gs, "CD69+ CD4 T Cells")[[1]])
#' @export
cyto_gate_type <- function(gates) {
  
  # One gate supplied
  if (length(gates) == 1) {
    if (class(gates) %in% c("filters", "list")) {
      gates <- gates[[1]]
    }
    # ELLIPSE
    if (class(gates) == "ellipsoidGate") {
      types <- "ellipse"
    # RECTANGLE/BOUNDARY/INTERVAL/THRESHOLD
    } else if (class(gates) == "rectangleGate") {
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
    } else if (class(gates) == "polygonGate") {
      types <- "polygon"
    # QUADRANT
    }else if(class(gates) == "quadGate"){
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
        } else if(class(x) == "quadGate"){
          types <- "quadrant"
        }
      })
    }
  }
  
  return(types)
}

## CYTO_GATE_CONVERT -----------------------------------------------------------

# CYTO_GATE_CONVERT is designed to convert constrcted gate objects to the number
# of dimensions specified by channels. Primary use is to easily convert between
# 1D and 2D gate objects.

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
cyto_gate_convert <- function(x, ...){
  UseMethod("cyto_gate_convert")
}

#' @noRd
#' @export
cyto_gate_convert.default <- function(x, 
                                      channels = NULL,
                                      ...){
  
  # Invalid gate object
  if(!is(x) %in% c("list", 
                      "filters",
                      "rectangleGate",
                      "polygonGate",
                      "ellipsoidGate",
                      "quadGate")){
    stop(paste(
      "'x' should contain either filters, rectangleGate, polygonGate,",
      "ellipsoidGate or quadGate objects."))
  }
  
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.rectangleGate <- function(x, 
                                            channels = NULL,
                                            ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)){
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)

  # PREPARE GATE ---------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # 1D GATE IN 1D PLOT - CORRECT CHANNEL
    if(length(chans) == 1 & all(chans %in% channels)){
      return(x)
    # 1D GATE IN 1D PLOT - INCORRECT CHANNEL
    }else if(length(chans) == 1 & all(!chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channel.")
    # 2D GATE IN 1D PLOT - CORRECT CHANNEL
    }else if(length(chans) == 2 & any(chans == channels)){
      return(x[channels])
    # 2D GATE IN 1D PLOT - INCORRECT CHANNEL
    }else if(length(chans) == 2 & !any(chans == channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channel.")
    }
  # 2D PLOT
  }else if(length(channels) == 2){
    # 2D GATE IN 2D PLOT - CORRECT CHANNELS
    if(length(chans) == 2 & all(chans %in% channels)){
      return(x)
    # 2D GATE in 2D PLOT - ONE CHANNEL MATCH
    }else if(length(chans) == 2 & any(chans %in% channels)){
      # FITERID
      ID <- x@filterId
      # CHANNEL INDEX
      ind <- which(chans %in% channels)
      # COORDS - MISSING CHANNEL
      coords <- matrix(c(x@min[channels[ind]], x@max[channels[ind]], -Inf, Inf),
                       ncol = 2,
                       nrow = 2,
                       byrow = FALSE)
      colnames(coords) <- c(channels[ind], channels[-ind])
      x <- rectangleGate(.gate = coords)
      return(x)
    # 2D GATE IN 2D PLOT - NO CHANNELS MATCH
    }else if(length(chans) == 2 & ! any(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    # 1D GATE IN 2D PLOT - CORRECT CHANNEL
    }else if(length(chans) == 1 & all(chans %in% channels)){
      # FITERID
      ID <- x@filterId
      # CHANNEL INDEX
      ind <- which(chans %in% channels)
      # COORDS - MISSING CHANNEL
      coords <- matrix(c(x@min, x@max, -Inf, Inf),
                       ncol = 2,
                       nrow = 2,
                       byrow = FALSE)
      colnames(coords) <- c(channels[ind], channels[-ind])
      x <- rectangleGate(.gate = coords)
      return(x)
    # 1D GATE IN 2D PLOT - INCORRECT CHANNEL  
    }else if(length(chans) == 1 & !all(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.polygonGate <- function(x, 
                                          channels = NULL,
                                          ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)){
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)

  # PREPARE GATE ---------------------------------------------------------------
  
  # 1D PLOT - CONVERT TO RECTANGLEGATE
  if(length(channels) == 1){
    # 2D GATE IN 1D PLOT - CORERCT CHANNEL
    if(any(chans %in% channels)){
      # FILTERID
      ID <- x@filterId
      # COORDS IN CHANNEL
      coords <- x@boundaries[, channels[channels %in% chans]]
      coords <- c(min(coords), max(coords))
      # RECTANGLEGATE
      coords <- matrix(coords,
                       ncol = 1,
                       nrow = 2)
      colnames(coords) <- channels[channels %in% chans]
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 1D RECTANGLEGATE
      return(x)
    # 2D GATE IN 1D PLOT - INCORRECT CHANNELS
    }else if(!any(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  # 2D PLOT
  }else if(length(channels) == 2){
    # 2D GATE IN 2D PLOT - CORRECT CHANNELS
    if(all(chans %in% channels)){
      return(x)
    # 2D GATE IN 2D PLOT - SINGLE CHANNEL 
    }else if(any(chans %in% channels)){
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
                       byrow = FALSE)
      colnames(coords) <- c(channels[channels %in% chans],
                            !channels[channels %in% chans])
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 2D RECTANGLEGATE
      return(x)
    # 2D GATE IN 2D PLOT - INCORRECT CHANNELS
    }else if(!any(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.ellipsoidGate <- function(x, 
                                            channels = NULL,
                                            ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)){
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)

  # PREPARE GATE ---------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # 2D GATE IN 1D PLOT - CHANNEL MATCH
    if(any(chans %in% channels)){
      # FILTERID
      ID <- x@filterId
      # POLYGONGATE
      x <- as(x, "polygonGate")
      # COORDS IN CHANNEL
      coords <- x@boundaries[, channels[channels %in% chans]]
      coords <- c(min(coords), max(coords))
      print(coords)
      # RECTANGLEGATE
      coords <- matrix(coords,
                       ncol = 1,
                       nrow = 2)
      colnames(coords) <- channels[channels %in% chans]
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 1D RECTANGLEGATE
      return(x)
    # 2D GATE IN 1D PLOT - NO CHANNEL MATCH
    }else if(!any(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  # 2D PLOT
  }else if(length(channels) == 2){
    # 2D GATE IN 2D PLOT - CORRECT CHANNELS
    if(all(chans %in% channels)){
      return(x)
    # 2D GATE IN 2D PLOT - SINGLE CHANNEL MATCH
    }else if(any(chans %in% channels)){
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
                       byrow = FALSE)
      colnames(coords) <- c(channels[channels %in% chans],
                            channels[!channels %in% chans])
      x <- rectangleGate(.gate = coords, filterId = ID)
      # RETURN 2D RECTANGLEGATE
      return(x)
    # 2D GATE IN 2D PLOT - NO CHANNEL MATCH 
    }else if(!any(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.quadGate <- function(x, 
                                       channels = NULL,
                                       ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)){
    stop("Supply the channels in which the gate will be plotted.")
  }
  
  # GATE CHANNELS
  chans <- parameters(x)
  
  # PREPARE GATES --------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    stop("Quadrant gates cannot be converted into 1D gate objects.")
  # 2D PLOT
  }else if(length(channels) == 2){
    # 2D GATES IN 2D PLOT - CORRECT CHANNELS
    if(all(chans %in% channels)){
      return(x)
    # 2D GATES IN 2D PLOT - INCORRECT CHANNELS
    }else if(!any(chans %in% channels)){
      stop("Supplied gate must contain co-ordinates in the supplied channels.")
    }
    
  }
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.filters <- function(x, 
                                      channels = NULL,
                                      ...){
  
  # GATE OBJECT LIST -----------------------------------------------------------
  x <- unlist(x)
  
  # LIST METHOD CALL -----------------------------------------------------------
  x <- cyto_gate_convert(x,
                         channels)
  
  # RETURN CONVERTED GATES -----------------------------------------------------
  return(x)
  
}

#' @rdname cyto_gate_convert
#' @export
cyto_gate_convert.list <- function(x, 
                                   channels = NULL,
                                   ...){
  
  # GATE OBJECT LIST -----------------------------------------------------------
  x <- unlist(x)
  
  # CONVERT GATES --------------------------------------------------------------
  x <- lapply(x, function(z){
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
                              channels = NULL){
  
  # LIST OF GATES --------------------------------------------------------------
  
  # PREPARE GATE LIST
  if(is(x)[1] == "list"){
    if(all(LAPPLY(x, "is") %in% c("rectangleGate",
                                  "polygonGate",
                                  "ellipsoidGate",
                                  "quadGate",
                                  "filters"))){
      x <- unlist(x)
    }
  }else if(is(x)[1] == "filters"){
    x <- unlist(x)
  }else if(is(x)[1] %in% c("rectangleGate",
                          "polygonGate",
                          "ellipsoidGate",
                          "quadGate",
                          "filters")){
    x <- list(x)
  }
  
  # UNIQUE GATE LIST
  x <- unique(x)
  
  # DIMENSIONS
  x <- cyto_gate_convert(x, channels = channels)
  
  # RETURN PREPARED GATE LIST
  return(x)
  
}
