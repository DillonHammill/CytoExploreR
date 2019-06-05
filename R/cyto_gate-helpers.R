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
#' library(CytoRSuiteData)
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
  
  # Supply alias
  if (is.null(alias)) {
    stop("Please supply the name of the population to be removed.")
  }
  
  # Check Alias
  if (!all(alias %in% basename(getNodes(gs)))) {
    stop("Supplied alias does not exist in the GatingSet.")
  }
  
  # Missing gatingTemplate
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # gatingTemplate still missing
  if(is.null(gatingTemplate)){
    stop("Suppply the name of the gatingTemplate to remove gate(s).")
  }
  
  # File extension
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # Get children from GatingSet
  chldrn <- unlist(lapply(
    alias,
    function(x) basename(getDescendants(gs[[1]], paste0(parent,"/",x)))
  ))
  chldrn <- unlist(chldrn, use.names = FALSE)
  chldrn <- c(alias, unique(chldrn))
  
  # Read in gatingTemplate
  gt <- read.csv(gatingTemplate, header = TRUE)
  
  # Remove all rows with alias = chldrn
  gt <- gt[!gt$alias %in% chldrn, ]
  
  # Message
  message(paste0("Removing gate(s) from the GatingSet and ",
                gatingTemplate,"."))
  
  # Remove nodes from GatingSet
  for (i in seq_len(length(alias))) {
    if (alias[i] %in% basename(getNodes(gs))) {
      suppressMessages(Rm(paste0(parent,"/", alias[i]), gs))
    }
  }
  
  write.csv(gt, gatingTemplate, row.names = FALSE)
}

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
#' library(CytoRSuiteData)
#' 
#' # Bypass working directory check for external files
#' options("CytoRSuite_wd_check" = FALSE)
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
#'   package = "CytoRSuiteData"
#' )
#' 
#' # Extract T Cells gate
#' gate_extract("Live Cells", "T Cells", gatingTemplate = gtfile)
#' 
#' # Reset CytoRSuite_wd_check to default
#' options("CytoRSuite_wd_check" = TRUE)
#' @export
cyto_gate_extract <- function(parent,
                         alias,
                         gatingTemplate = NULL, ...) {
  
  # Parent
  if (missing(parent)) {
    stop("Please supply the name of the parent population.")
  }
  
  # Alias
  if (missing(alias)) {
    stop("Please supply the name(s) of the alias to extract.")
  }
  
  # Missing gatingTemplate
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # gatingTemplate still missing
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to extract gate(s).")
  }
  
  # File extension
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # Load gatingTemplate into gatingTemplate
  gt <- suppressMessages(gatingTemplate(gatingTemplate))
  
  # Extract population nodes from gt
  nds <- getNodes(gt, only.names = TRUE)
  
  # Parent Node
  parent <- names(nds[match(parent, nds)])
  
  # Extract gates given parent and child node(s)
  gates <- lapply(alias, function(x) {
    
    # Alias node
    alias <- names(nds[match(x, nds)])
    
    gm <- getGate(gt, parent, alias)
    gate <- eval(parameters(gm)$gate)
    names(parameters(gate[[1]][[1]])) <- parameters(gate[[1]][[1]])
    return(gate)
  })
  
  return(gates)
}

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
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param select vector containing the indicies of samples within gs to use for
#'   plotting.
#' @param overlay name(s) of the population(s) to overlay onto the plot.
#' @param ... additional arguments passed to cyto_plot, see ?cyto_plot for
#'   details.
#'
#' @return an object of class \code{GatingSet} with edited gate applied, as well
#'   as gatingTemplate file with edited gate saved.
#'
#' @importFrom flowWorkspace getData getTransformations GatingSet getGate
#'   setGate recompute pData
#' @importFrom flowCore parameters filterList
#' @importFrom openCyto gatingTemplate
#' @importFrom data.table as.data.table fread :=
#' @importFrom methods as
#' @importFrom utils select.list write.csv
#' @importFrom tools file_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
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
                      display = NULL,
                      select = NULL,
                      overlay = NULL, ...) {
  
  # Parent
  if (is.null(parent)) {
    stop("Please supply the name of the parent population.")
  } else if (!parent %in% basename(getNodes(x))) {
    stop("Supplied parent does not exist in the GatingSet.")
  }
  
  # Alias
  if (is.null(alias)) {
    stop("Please supply the name(s) of the gates to edit to 'alias'.")
  } else if (!all(alias %in% basename(getNodes(x)))) {
    stop("Supplied alias does not exist in the GatingSet.")
  }
  
  # Missing gatingTemplate
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # gatingTemplate still missing
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to edit gate(s).")
  }
  
  # File extension
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # Assign x to gs
  gs <- x
  
  # Extract pData information
  pd <- cyto_details(gs)
  
  # Extract transList from gs
  if (length(getTransformations(gs[[1]])) != 0) {
    transList <- transformList(
      names(getTransformations(gs[[1]])),
      getTransformations(gs[[1]])
    )
  } else {
    transList <- NULL
  }
  
  # Extract gates from gt
  gT <- suppressMessages(gatingTemplate(gatingTemplate))
  
  # Extract population nodes from gt
  nds <- getNodes(gT, only.names = TRUE)
  
  # Parent Node
  prnt <- names(nds)[match(parent, nds)]
  
  # Read in gatingTemplate
  gt <- read.csv(gatingTemplate, header = TRUE)
  
  # Get groupBy from gatingTemplate
  grpby <- as.character(gt[gt$parent == parent &
                             gt$alias == alias[1], "groupBy"])
  
  # groupBy is NA
  if (all(is.na(grpby))) {
    grpby <- length(gs)
  } else {
    grpby <- unlist(strsplit(grpby, ":"))
  }
  
  # groupBy is numeric
  if (!all(grepl("^[A-Za-z]+$", grpby))) {
    grpby <- as.numeric(grpby)
    if (grpby > length(gs)) {
      grpby <- length(gs)
    } else if (grpby < length(gs)) {
      stop("Numeric groupBy must be equal to the length of the GatingSet.")
    }
  }
  
  # Get channels from gatingTemplate
  channels <- unique(
    unlist(
      strsplit(
        as.character(gt[gt$parent == parent &
                          gt$alias %in% alias, "dims"]), ",",
        fixed = TRUE
      )
    )
  )
  
  # Menu to select which groups require editing
  if (is.numeric(grpby)) {
    
    # All samples in same group
    if (grpby >= length(gs)) {
      
      # no selection required - grps indicates which filters object to edit
      grps <- 1
      pd$groupby <- rep(1, length(gs))
    } else {
      stop("Numeric groupBy is not currently supported.")
    }
  } else if (is.character(grpby)) {
    vrs <- unlist(strsplit(grpby, ",", fixed = TRUE))
    opts <- unique(do.call(paste, pd[, grpby, drop = FALSE]))
    
    pd$groupby <- do.call(paste, pd[, grpby, drop = FALSE])
    
    grps <- select.list(opts,
                        multiple = TRUE,
                        graphics = TRUE,
                        title = "Select the group(s) to edit:"
    )
  }
  
  # Extract gates given parent and child node(s)
  gt_gates <- lapply(seq_along(alias), function(x) {
    
    # Alias node
    als <- names(nds[match(alias[x], nds)])
    gm <- getGate(gT, prnt, als)
    gate <- eval(parameters(gm)$gate)
    names(gate) <- unique(pd$groupby)
    
    return(gate)
  })
  names(gt_gates) <- alias
  # gt_gates is list (length(alias)) of list of filters (1x filters per group)
  
  # Split GatingSet into list of GatingSet groups
  if (is.numeric(grpby)) {
    gs.lst <- list(gs)
    names(gs.lst) <- 1
  } else if (is.character(grpby)) {
    gs.lst <- lapply(grps, function(x) gs[which(pd$groupby == x)])
    names(gs.lst) <- grps
  }
  
  # plot, edit and save gate for each group
  new_gates <- lapply(gs.lst, function(grp) {
    
    # Extract parent population for plotting
    fs <- suppressMessages(getData(grp, parent))
    fr <- as(fs, "flowFrame")
    
    # Remove "Original" column introduced by coercion
    if (is.na(match("Original", BiocGenerics::colnames(fr))) == FALSE) {
      fr <- suppressWarnings(
        fr[, -match("Original", BiocGenerics::colnames(fr))]
      )
    }
    
    # Display
    if (is.null(display)) {
      display <- 1 / length(fs)
    }
    
    # Overlay
    if (!is.null(overlay)) {
      
      # Extract populations to overlay - list of flowSets
      if (class(overlay) == "character") {
        overlay <- lapply(overlay, function(overlay) {
          getData(grp, overlay)
        })
        
        overlay <- lapply(overlay, function(x) {
          fr <- as(x, "flowFrame")
          
          if (is.na(match("Original", BiocGenerics::colnames(fr))) == FALSE) {
            fr <- suppressWarnings(
              fr[, -match("Original", BiocGenerics::colnames(fr))]
            )
          }
          
          return(fr)
        })
      }
    }
    
    # Turn off overlay sampling
    options("CytoRSuite_overlay_display" = FALSE)
    
    # Extract gate(s) for plotting
    gates <- filters(lapply(alias, function(x) {
      getGate(grp[[1]], x)
    }))
    
    # Plot data and existing gates
    if (is.numeric(grpby)) {
      if (parent == "root") {
        pnt <- "All Events"
      } else {
        pnt <- parent
      }
      
      if (grpby >= length(gs)) {
        main <- paste("Combined Events", "\n", pnt)
      }
    } else if (is.character(grpby)) {
      if (parent == "root") {
        pnt <- "All Events"
      } else {
        pnt <- parent
      }
      
      main <- paste(
        unique(pd[pd$name %in% sampleNames(grp), "groupby"]),
        "\n",
        pnt
      )
    }
    
    if (getOption("CytoRSuite_interact") == TRUE) {
      cyto_plot(fr,
                channels = channels,
                overlay = overlay,
                display = display,
                popup = TRUE,
                legend = FALSE,
                gate = gates,
                gate_line_col = "magenta",
                axes_trans = transList,
                label = FALSE,
                title = main,
                gate_line_width = 2.5,
                legend_text = NA, ...
      )
    } else {
      cyto_plot(fr,
                channels = channels,
                overlay = overlay,
                display = display,
                legend = FALSE,
                gate = gates,
                gate_line_col = "magenta",
                axes_trans = transList,
                label = FALSE,
                title = main,
                gate_line_width = 2.5,
                legend_text = NA, ...
      )
    }
    
    # Reset overlay sampling
    options("CytoRSuite_overlay_display" = TRUE)
    
    # If no type supplied determine using gate_type
    if (is.null(type)) {
      type <- gate_type(gates)
    }
    
    # Check type argument is valid
    type <- .cyto_gate_type_check(type = type, alias = alias)
    
    # Check alias is supplied correctly
    .cyto_alias_check(alias = alias, type = type)
    
    # Extract channels from gates
    channels <- unique(unlist(lapply(gates, parameters)))
    
    # 2D Interval gates require axis argument
    if ("interval" %in% type) {
      intvl <- rbind(
        gates[[match("interval", type)[1]]]@min,
        gates[[match("interval", type)[1]]]@max
      )
      
      if (all(is.finite(intvl[, 1]))) {
        axis <- "x"
      } else if (all(is.finite(intvl[, 2]))) {
        axis <- "y"
      }
    }
    
    # Draw new gates - set plot to FALSE (filters object of length alias)
    new_gates <- cyto_gate_draw(fr,
                           alias = alias,
                           channels = channels,
                           type = type,
                           axis = axis,
                           plot = FALSE
    )
    names(new_gates) <- alias
    
    # Modify existing gate(s)
    lapply(seq_along(alias), function(pop) {
      gt_gates[[alias[pop]]][[match(
        pd[pd$name %in% sampleNames(grp), "groupby"][1], names(gt_gates[[pop]])
      )]] <<- filters(list(new_gates[[alias[pop]]]))
    })
    
    return(new_gates)
  })
  
  
  # Re-name parent and pop to be data.table friendly
  prnt <- parent
  als <- alias
  gtmd <- "cyto_gate_draw"
  ppmd <- "pp_cyto_gate_draw"
  
  # Find and Edit gatingTemplate entries - each alias and gate separate
  gt <- data.table::fread(gatingTemplate)
  
  # data.table R CMD Check NOTE
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  preprocessing_args <- NULL
  
  # Prepare grpby
  if (length(grpby) == 1) {
    if (!all(grepl("^[A-Za-z]+$", grpby)) & grpby[1] == length(gs)) {
      grpby <- as.logical(NA)
    }
  } else {
    grpby <- paste(grpby, collapse = ":")
  }
  
  # Modify template
  for (i in seq_len(length(alias))) {
    gt[parent == prnt & alias == als[i], gating_method := gtmd]
    gt[parent == prnt & alias == als[i], gating_args := .argDeparser(list(
      gate = gt_gates[[i]]
    ))]
    gt[parent == prnt & alias == als[i], collapseDataForGating := TRUE]
    gt[parent == prnt & alias == als[i], groupBy := grpby]
    gt[parent == prnt & alias == als[i], preprocessing_method := ppmd]
    gt[parent == prnt & alias == als[i], preprocessing_args := as.logical(NA)]
  }
  
  # Save updated gatingTemplate
  write.csv(as.data.frame(gt), gatingTemplate, row.names = FALSE)
  
  # Apply new gates to GatingSet by group
  lapply(grps, function(grp) {
    lapply(seq_along(alias), function(pop) {
      fltrs <- rep(
        list(gt_gates[[alias[pop]]][[grp]][[1]]),
        length(gs[which(pd$groupby == grp)])
      )
      names(fltrs) <- as.character(sampleNames(gs[which(pd$groupby == grp)]))
      
      suppressMessages(setGate(gs[which(pd$groupby == grp)], alias[pop], fltrs))
      suppressMessages(recompute(gs[which(pd$groupby == grp)], alias[pop]))
    })
  })
  
  assign(deparse(substitute(x)), gs, envir = globalenv())
}

#' Rename Gates
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param gates current names of the gates to be changed.
#' @param names new names to use for \code{gates}.
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
                        gates = NULL,
                        names = NULL,
                        gatingTemplate = NULL){
  
  # Missing gatingTemplate
  if(is.null(gatingTemplate)){
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # gatingTemplate still missing
  if(is.null(gatingTemplate)){
    stop("Supply the name of the gatingTemplate to rename gate(s).")
  }
  
  # File extension
  if(.empty(file_ext(gatingTemplate))){
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }
  
  # Gates do not exist
  if(!all(gates %in% basename(getNodes(x)))){
    stop(paste0("Supplied gate(s) do not exist in this ", class(x), "."))
  }
  
  # Rename gates in GatingHierarchy/GatingSet
  mapply(function(gate, name){
    
    setNode(x, gate, name)
    
  }, gates, names)
  
  # Modify names in gatingTemplate
  gt <- read.csv(gatingTemplate, 
                 header = TRUE, 
                 stringsAsFactors = FALSE)
  
  # Update parental names
  if(any(gates %in% gt$parent)){
    
    gt[gt$parent %in% gates,"parent"] <- names
  
  }
  
  # Update alias names
  gt[gt$alias %in% gates,"alias"] <- names
  
  # Re-write updated gatingTemplate
  write.csv(as.data.frame(gt), gatingTemplate, row.names = FALSE)
  
}

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
#' library(CytoRSuiteData)
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
#' # Get gate type used for T Cells gate
#' gate_type(getGate(gs, "Cells")[[1]])
#' gate_type(getGate(gs, "T Cells")[[1]])
#' gate_type(getGate(gs, "CD69+ CD4 T Cells")[[1]])
#' @export
cyto_gate_type <- function(gates) {
  
  # One gate supplied
  if (length(gates) == 1) {
    if (class(gates) %in% c("filters", "list")) {
      gates <- gates[[1]]
    }
    
    if (class(gates) == "ellipsoidGate") {
      
      # type == "ellipse"
      types <- "ellipse"
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
    } else if (class(gates) == "polygonGate") {
      
      # type == "polygon"
      types <- "polygon"
    }
    
    # Multiple gates supplied
  } else if (length(gates) > 1) {
    
    # Get classes of gates
    classes <- unlist(lapply(gates, function(x) {
      class(x)
    }))
    
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
            types <- unlist(lapply(gates, function(x) {
              
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
            }))
          }
        } else {
          types <- unlist(lapply(gates, function(x) {
            
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
          }))
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
      types <- unlist(lapply(gates, function(x) {
        if (class(x) == "ellipsoidGate") {
          
          # type == "ellipse"
          types <- "ellipse"
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
        } else if (class(x) == "polygonGate") {
          
          # type == "polygon"
          types <- "polygon"
        }
      }))
    }
  }
  
  return(types)
}
