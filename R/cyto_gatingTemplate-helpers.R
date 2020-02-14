## GATINGTEMPLATE SELECT -------------------------------------------------------

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

## GATINGTEMPLATE ACTIVE -------------------------------------------------------

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

## GATINGTEMPLATE CREATE -------------------------------------------------------

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

## GATINGTEMPLATE EDIT ---------------------------------------------------------

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
#'
#' @return update gatingTemplate csv file and update the GatingSet accordingly.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom flowWorkspace recompute
#' @importFrom utils edit read.csv write.csv
#' @importFrom tools file_ext
#' @importFrom methods is
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
cyto_gatingTemplate_edit <- function(x, gatingTemplate = NULL) {

  # x of wrong class
  if (!is(x, "GatingSet")) {
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
#' @param ... additional arguments passed to gating.
#'
#' @return NULL and update the GatingHierarchy/ GatingSet in the global
#'   environment with the gates in the gatingTemplate.
#'
#' @importFrom openCyto gatingTemplate gt_gating
#' @importFrom tools file_ext
#' @importFrom methods is
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
                                      gatingTemplate = NULL, ...) {

  # No GatingSet supplied
  if (missing(x)) {
    stop("Supply a GatingHierarchy or GatingSet to 'x'.")
  } else {
    if (!any(c(
      is(x, "GatingSet"),
      is(x, "GatingHierarchy")
    ))) {
      stop("'x' must be an object of class GatingHierarchy or GatingSet.")
    }
  }

  # gatingTemplate supplied
  if (!is.null(gatingTemplate)) {

    # gatingTemplate object
    if (is(gatingTemplate, "gatingTemplate")) {

      # Apply gatingTemplates to GatingSets only
      if (!is(x, "GatingSet")) {
        stop("gatingTemplates can only be applied to GatingSet objects.")
      }

      # Assign gatingTemplate to gt
      gt <- gatingTemplate
      
      # name of gatingTemplate csv file
    } else {

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
  } else if (is.null(gatingTemplate)) {

    # Missing gatingTemplate - check global option
    if (is.null(gatingTemplate)) {
      gatingTemplate <- cyto_gatingTemplate_active()
    }

    # gatingTemplate still missing
    if (is.null(gatingTemplate)) {
      stop("Supply the name of the gatingTemplate csv file to apply.")
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
  suppressWarnings(gt_gating(gt, x, ...))
  
  # Return GatingSet
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
#' @importFrom utils read.csv write.csv
#'
#' @export
cyto_gatingTemplate_update <- function(gatingTemplate = NULL,
                                       save_as = "Updated-gatingTemplate.csv") {

  # READ IN GATINGTEMPLATE
  gt <- read.csv(gatingTemplate, 
                 header = TRUE,
                 stringsAsFactors = FALSE)
  
  # CONVERT OLD GATING METHOD
  gt[gt$gating_method == "gate_draw", "gating_method"] <- "cyto_gate_draw" 
  
  # CONVERT OLD PREPROCESSING  METHOD
  gt[gt$preprocessing_method == "pp_gate_draw",
     "preprocessing_method"] <- "pp_cyto_gate_draw" 

  # SAVE UPDATED GATINGTEMPLATE
  write.csv(gt, save_as, row.names = FALSE)
}

## GATINGTEMPLATE CHECK --------------------------------------------------------

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

  # BYPASS CHECKS - GATINGTEMPLATE DOES NOT EXIST
  gt <- tryCatch(suppressWarnings(read.csv(gatingTemplate, header = TRUE)),
    error = function(e) {
      NULL
    }
  )
  # MATCHING PARENT & ALIAS
  if (!is.null(gt)) {
    if (parent %in% gt$parent) {
      gt <- gt[gt$parent == parent, ]
      gt_alias <- lapply(gt$alias, function(z) {
        unlist(strsplit(as.character(z), ","))
      })
      lapply(gt_alias, function(z) {
        if (any(alias %in% z)) {
          message(
            paste(alias[alias %in% z], collapse = " & "),
            " already exists in ",
            gatingTemplate, "."
          )
          stop(paste(
            "Supply another gatingTemplate",
            "or edit gate(s) using cyto_gate_edit."
          ))
        }
      })
    }
  }
}
