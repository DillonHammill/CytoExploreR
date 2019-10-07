# GATINGTEMPLATE SELECT ---------------------------------------------------------

#' Select a gatingTemplate for downstream analyses
#'
#' \code{cyto_gatingTemplate_select} will designate which gatingTemplate should
#' be used for storing downstream gating results.
#' \code{cyto_gatingTemplate_select} should therefore be called prior to gating
#' or \code{gate_draw} will resort to using \code{"gatingTemplate.csv"} to save
#' the constructed gates. \code{cyto_gatingTemplate_select} also provides an
#' easy way to switch between gatingTemplates when multiple templates are
#' required.
#'
#' @param x name of the gatingTemplate csv file to assign as the active
#'   gatingTemplate for downstream analyses.
#'
#' @return select a gatingTemplate as the active gatingTemplate for downstream
#'   analyses.
#'
#' @importFrom tools file_ext
#'
#' @examples
#' cyto_gatingTemplate_select("Activation_gatingTemplate")
#' @seealso \code{\link{cyto_gatingTemplate_edit}}
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_gatingTemplate_select <- function(x) {

  # Name does not include a file extension
  if (.empty(file_ext(x))) {
    x <- paste0(x, ".csv")
  }

  # Name does not include a csv file extension
  if (file_ext(x) != "csv") {
    x <- paste0(x, ".csv")
  }

  # Set new gatingTemplate as active
  options("CytoExploreR_gatingTemplate" = x)
}

# GATINGTEMPLATE ACTIVE --------------------------------------------------------

#' Active gatingTemplate
#'
#' @return name of the last and therefore active gatingTemplate assigned by
#'   \code{cyto_gatingTemplate_select}.
#'
#' @seealso \code{\link{cyto_gatingTemplate_select}}
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_gatingTemplate_active <- function() {
  getOption("CytoExploreR_gatingTemplate")
}

# GATINGTEMPLATE CREATE --------------------------------------------------------

#' Create an empty gatingTemplate csv file
#'
#' @param gatingTemplate name of the gatingTemplate csv file to create.
#'
#' @importFrom utils write.csv
#' @importFrom stats setNames
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_gatingTemplate_create <- function(gatingTemplate = NULL) {

  # Inherit global option
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }

  # Create gatingTemplate for saving
  if (!is.null(gatingTemplate)) {
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

    # Make empty gatingTemplate
    gt <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)

    # Write to csv file
    write.csv(gt, gatingTemplate, row.names = FALSE)
  }
}

# GATINGTEMPLATE EDIT ----------------------------------------------------------

#' Interactively Edit a gatingTemplate
#'
#' \code{cyto_gatingTemplate_edit} provides an interactive interface to editing the
#' openCyto gatingTemplate. This function is intended to aid in adding boolean
#' and reference gates to the gatingTemplate. Users should NOT modify existing
#' entries in the gatingTemplate.
#'
#' @param x object of class \code{GatingSet} which has been gated with the
#' supplied gatingTemplate.
#' @param gatingTemplate name of the gatingTemplate csv file to be edited.
#'
#' @return update gatingTemplate csv file and update the GatingSet accordingly.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom openCyto gatingTemplate gating
#' @importFrom flowWorkspace recompute
#' @importFrom utils edit read.csv write.csv
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#' library(CytoRSuite)
#'
#' # gs is a GatingSet object
#' cyto_gatingTemplate_edit(gs, "gatingTemplate.csv")
#' }
#'
#' @export
cyto_gatingTemplate_edit <- function(x, gatingTemplate = NULL) {

  # x of wrong class
  if (!inherits(x, "GatingSet")) {
    stop("'x' should be either a GatingSet object.")
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

  # File extension
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # Working directory check
  if (getOption("CytoExploreR_wd_check") == TRUE) {
    if (file_wd_check(gatingTemplate) == FALSE) {
      stop(paste(gatingTemplate, "is not in this working directory."))
    }
  }

  # Read in gatingTemplate
  gt <- read.csv(gatingTemplate, header = TRUE)

  # Edit gatingTemplate
  message("Do not modify existing gatingTemplate entries!")
  message("Add new rows to add boolean or reference gates to the GatingSet.")
  gt <- suppressMessages(edit(gt))

  # Write updated template to csv file
  write.csv(gt, gatingTemplate, row.names = FALSE)

  # Read in updated template to gatingTemplate object
  gt <- gatingTemplate(gatingTemplate)

  # Re-apply template to GatingSet
  suppressMessages(gating(gt, gs))

  # Recompute statistics
  suppressMessages(recompute(gs))

  # Update GatingSet globally
  assign(deparse(substitute(x)), gs, envir = globalenv())

  # Return gatingTemplate object
  invisible(return(gt))
}

# GATINGTEMPLATE APPLY ---------------------------------------------------------

#' Apply gates saved in gatingTemplate to a GatingHierarchy or GatingSet
#'
#' A convenient wrapper for
#' \code{\link[openCyto:gatingTemplate]{gatingTemplate}} and
#' \code{\link[openCyto:gating]{gating}} to apply a gatingTemplate csv file to a
#' GatingSet.
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param gatingTemplate name of the gatingTemplate csv file which contains the
#'   gates to be applied to \code{x}. This can alos be an object of class
#'   \code{gatingTemplate} if a \code{GatingSet} object is supplied.
#' @param ... additional arguments passed to gating.
#'
#' @return NULL and update the GatingHierarchy/ GatingSet in the global
#'   environment with the gates in the gatingTemplate.
#'
#' @importFrom openCyto gatingTemplate gating
#' @importFrom tools file_ext
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
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(fs))
#' gs <- transform(gs, trans)
#'
#' # Apply gatingTemplate in working directory to GatingSet
#' cyto_gatingTemplate_apply(gs, "gatingTemplate.csv")
#' }
#'
#' @export
cyto_gatingTemplate_apply <- function(x,
                                      gatingTemplate = NULL, ...) {

  # No GatingSet supplied
  if (missing(x)) {
    stop("Supply a GatingHierarchy or GatingSet to 'x'.")
  } else {
    if (!any(c(
      inherits(x, "GatingSet"),
      inherits(x, "gatinghierarchy")
    ))) {
      stop("'x' must be an object of class GtaingHierarchy or GatingSet.")
    }
  }

  # gatingTemplate supplied
  if(!is.null(gatingTemplate)){
    
    # gatingTemplate object
    if(inherits(gatingTemplate, "gatingTemplate")){
    
      # Apply gatingTemplates to GatingSets only
      if(!inherits(x, "GatingSet")){
        stop("gatingTemplates can only be applied to GatingSet objects.")
      }
      
    # name of gatingTemplate csv file  
    }else{
      
      # File extension if missing
      if (.empty(file_ext(gatingTemplate))) {
        gatingTemplate <- paste0(gatingTemplate, ".csv")
      }
      
      # Working directory check
      if (getOption("CytoExploreR_wd_check") == TRUE) {
        if (file_wd_check(gatingTemplate) == FALSE) {
          stop(paste(gatingTemplate, "is not in this working directory."))
        }
      }
      
      # Message to indicate which gatingTemplate is being applied
      message(paste("Applying", gatingTemplate, "to the GatingSet..."))
      
      # Read in gatingTemplate csv file to gatingTemplate object
      gt <- suppressMessages(gatingTemplate(gatingTemplate))
      
    }
  # no gatingTemplate supplied
  }else if(is.null(gatingTemplate)){
    
    # Missing gatingTemplate - check global option
    if (is.null(gatingTemplate)) {
      gatingTemplate <- cyto_gatingTemplate_active()
    }
    
    # gatingTemplate still missing
    if (is.null(gatingTemplate)) {
      stop("Supply the name of the gatingtemplate csv file to apply.")
    }
    
    # File extension if missing
    if (.empty(file_ext(gatingTemplate))) {
      gatingTemplate <- paste0(gatingTemplate, ".csv")
    }
    
    # Working directory check
    if (getOption("CytoExploreR_wd_check") == TRUE) {
      if (file_wd_check(gatingTemplate) == FALSE) {
        stop(paste(gatingTemplate, "is not in this working directory."))
      }
    }
    
    # Message to indicate which gatingTemplate is being applied
    message(paste("Applying", gatingTemplate, "to the GatingSet..."))
    
    # Read in gatingTemplate csv file to gatingTemplate object
    gt <- suppressMessages(gatingTemplate(gatingTemplate))
    
  }
  
  # Apply gatingTemplate to GatingHierarchy/GatingSet
  suppressWarnings(gating(gt, x, ...))
  
}

# GATINGTEMPLATE CONVERT -------------------------------------------------------

#' Convert Old gatingTemplates to New Format
#'
#' @param gs object of class GatingSet.
#' @param gatingTemplate name of the gatingTemplate csv file to convert.
#'
#' @return converted gatingTemplate saved to .csv file.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom data.table fread fwrite
#' @importFrom openCyto gatingTemplate
#' @importFrom flowWorkspace getGate getNodes
#' @importFrom flowCore parameters filters
#' @importFrom utils read.csv write.csv
#' @importFrom tools file_ext
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
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(fs))
#' gs <- transform(gs, trans)
#'
#' # Convert gatingTemplate CytoRSuite <= 0.9.5 to new format
#' cyto_gatingTemplate_convert(gs, "gatingTemplate.csv")
#'
#' # Updated gatingTemplate will work as expected
#' gt <- gatingTemplate("gatingTemplate.csv")
#' gating(gt, gs)
#' }
#'
#' @export
cyto_gatingTemplate_convert <- function(gs, gatingTemplate = NULL) {

  # data.table R CMD Check NOTE
  alias <- NULL
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  preprocessing_args <- NULL

  # Number of samples
  smp <- length(gs)

  # Missing gatingTemplate - check global option
  if (is.null(gatingtemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }

  # gatingtemplate sill missing
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate csv file to convert.")
  }

  # File extension
  if (.empty(file_ext(gatingTemplate))) {
    gatingTemplate <- paste0(gatingTemplate, ".csv")
  }

  # Working directory check
  if (getOption("CytoExploreR_wd_check") == TRUE) {
    if (file_wd_check(gatingTemplate) == FALSE) {
      stop(paste(gatingTemplate, "is not in this working directory."))
    }
  }

  # Read in gatingTemplate
  gt <- data.table::fread(gatingTemplate)

  # Modify template
  for (i in seq_len(nrow(gt))) {

    # Load gatingTemplate into gatingTemplate
    gT <- suppressMessages(gatingTemplate(gatingTemplate))
    pops <- basename(gT@nodes)[-1]

    # Extract population nodes from gT
    nds <- getNodes(gT, only.names = TRUE)

    als <- names(nds[match(pops[i], nds)])

    # Parent Node
    parent <- gt[alias == pops[i], "parent"]
    parent <- as.character(parent)

    prnt <- names(nds[match(parent, nds)])
    gm <- getGate(gT, prnt, als)
    gate <- eval(parameters(gm)$gate)
    gts <- list(filters(list(gate)))

    gtmd <- "gate_draw"
    ppmd <- "pp_gate_draw"
    als <- pops[i]

    gt[alias == als, gating_method := gtmd]
    gt[alias == als, gating_args := .argDeparser(list(gate = gts))]
    gt[alias == als, collapseDataForGating := TRUE]
  }

  # Save updated gatingTemplate
  data.table::fwrite(gt, gatingTemplate)

  gt <- read.csv(gatingTemplate, header = TRUE)
  gt$collapseDataForGating <- rep(TRUE, nrow(gt))
  gt$groupBy <- rep(smp, nrow(gt))
  gt$preprocessing_method <- rep("pp_gate_draw", nrow(gt))
  gt$preprocessing_args <- rep(NA, nrow(gt))

  write.csv(gt, gatingTemplate, row.names = FALSE)
}

# GATINGTEMPLATE CHECK ---------------------------------------------------------

#' Check gatingTemplate for Existing Entry
#'
#' @param parent name of the parent population.
#' @param alias name of the population of interest.
#' @param gatingTemplate csv file name of the gatingTemplate.
#'
#' @return stops the gating process if an entry already exists in the
#'   gatingTemplate for the supplied alias.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom utils read.csv
#'
#' @noRd
.cyto_gatingTemplate_check <- function(parent, 
                                       alias, 
                                       gatingTemplate = NULL) {
  
  # GATINGTEMPLATE CSV FILE
  if (inherits(gatingTemplate, "gatingTemplate")) {
    stop("'gatingTemplate' should be the name of the gatingTemplate csv file.")
  } else {
    # BYPASS CHECKS - GATINGTEMPLATE DOES NOT EXIST
    gt <- tryCatch(suppressWarnings(read.csv(gatingTemplate, header = TRUE)),
      error = function(e) {
        NULL
      }
    )
    # MATCHING PARENT & ALIAS
    if (!is.null(gt)) {
      if(parent %in% gt$parent){
        gt <- gt[gt$parent == parent, ]
        gt_alias <- lapply(gt$alias, function(z){
          unlist(strsplit(as.character(z), ","))
        })
        lapply(gt_alias, function(z){
          if(any(alias %in% z)){
            message(paste(alias[alias %in% z], collapse = " & "),
                    " already exists in ",
                    gatingTemplate, ".")
            stop(paste(
              "Supply another gatingTemplate",
              "or edit gate(s) using cyto_gate_edit."
            ))
          }
        })
      }
    }
  }
}