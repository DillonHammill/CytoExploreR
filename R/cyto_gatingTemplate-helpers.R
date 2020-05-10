## GATINGTEMPLATE GENERATE -----------------------------------------------------

#' Generate a CytoExploreR gatingTemplate for a GatingHierarchy or GatingSet
#'
#' @param x object of class GatingHierarchy or GatingSet.
#' @param gatingTemplate name of a CSV file to which the generated
#'   gatingTemplate should be saved.
#' @param ... additional arguments passed to \code{gatingTemplate}.
#'
#' @return write generated gatingTemplate to designated CSV file and return the
#'   gatingTemplate object.
#'
#' @importFrom tools file_ext
#' @importFrom utils write.csv
#' @importFrom data.table as.data.table
#' @importFrom openCyto gh_generate_template gatingTemplate
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @name cyto_gatingTemplate_generate
NULL

#' @noRd
#' @export
cyto_gatingTemplate_generate <- function(x, 
                                         gatingTemplate = NULL,
                                         ...) {
  UseMethod("cyto_gatingTemplate_generate")
}

#' @rdname cyto_gatingTemplate_generate
#' @export
cyto_gatingTemplate_generate.GatingHierarchy <- function(x,
                                                         gatingTemplate = NULL,
                                                         ...){
  
  # GATINGTEMPLATE MISSING
  if(is.null(gatingTemplate)){
    stop("Supply the name of a csv file to 'gatingTemplate'.")
  }else{
    # FILE EXTENSION MISSING
    if(file_ext(gatingTemplate) != "csv"){
      gatingTemplate <- paste0(gatingTemplate, ".csv")
    }
  }
  
  # DATA.TABLE R CMD CHECK NOTE
  parent <- NULL
  alias <- NULL
  pop <- NULL
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  
  # GENERATE GATINGTEMPLATE - DATA FRAME
  gt <- gh_generate_template(x)

  # BOOLEAN ENTRIES EXIST
  if(any(LAPPLY(gt$dims, ".empty"))){
    # PREPARE BOOLEAN ENTRIES
    gt_bool_chunk <- gt[LAPPLY(gt$dims, ".empty"), ]
    gt_bool_chunk$parent <- cyto_nodes_convert(x,
                                               nodes = gt_bool_chunk$parent,
                                               path = "auto")
    gt_bool_chunk$gating_args <- LAPPLY(
      seq_len(nrow(gt_bool_chunk)), function(z){
         gate <- gh_pop_get_gate(x,
                                 gt_bool_chunk$alias[z])
        return(gate@deparse)
      }
    )
    gt_bool_chunk$gating_method <- rep("boolGate", nrow(gt_bool_chunk))
    gt_bool_chunk <- as.data.table(gt_bool_chunk)
    # RESTRICT GATINGTEMPLATE
    gt <- gt[!LAPPLY(gt$dims, ".empty"), ]
  # NO BOOLEAN ENTRIES EXIST
  }else{
    gt_bool_chunk <- NULL
  }
  
  # SORT ENTRIES BASED ON PARENT
  gt <- gt[order(nchar(gt$parent)), ]
  
  # GATINGTEMPLATE INFO
  gt_parents <- unique(gt$parent)
  
  # GATINGTEMPLATE ENTRIES PER PARENT
  gt_entries <- lapply(gt_parents, function(gt_parent){
    # GATINGTEMPLATE CHUNK
    gt_parent_chunk <- gt[gt$parent == gt_parent, ]
    gt_parent_chunk_channels <- unique(gt_parent_chunk$dims)
    # PREPARE ENTRIES (PER PLOT)
    gt_entries <- lapply(gt_parent_chunk_channels, function(y){
      gt_chunk <- gt_parent_chunk[gt_parent_chunk$dims == y, ]
      gt_alias <- gt_chunk$alias
      gt_channels <- unlist(strsplit(y, ","))
      # GATES
      gates <- lapply(gt_alias, function(w){
        gh_pop_get_gate(x,
                        cyto_nodes_convert(x,
                                           nodes = w,
                                           anchor = gt_parent,
                                           path = "auto"))
      })
      names(gates) <- gt_alias
      # GATES BELONG TO QUADGATE
      if(length(gates) == 4 &
         all(LAPPLY(gates, "class") == "rectangleGate")){
        # RECTANGLES BELONG TO QUADGATE
        if(all(LAPPLY(gates, function(v){
          any(grepl("quad", names(attributes(v))))
        }))){
          gates <- list(.cyto_gate_quad_convert(gates,
                                                channels = gt_channels))
          names(gates) <- paste0(gt_alias, collapse = ",")
        }
      }
      # GATINGTEMPLATE CHUNK DATA TABLE 
      gt_chunk_dt <- as.data.table(gt_chunk)
      gt_entries <- lapply(seq_along(gates), function(q){
        # GATE
        gate <- gates[[q]]
        # QUAD GATE
        if(is(gate, "quadGate")){
          # USE FIRST ROW ONLY
          gt_entry <- gt_chunk_dt[1, ]
          # ALIAS
          gt_entry[, alias := cyto_nodes_convert(x,
                                                 nodes = unlist(
                                                   strsplit(names(gates), ",")
                                                 ),
                                                 anchor = gt_parent,
                                                 path = "auto")]
          # POP
          gt_entry[, pop := "*"]
          # PARENT
          gt_entry[, parent := cyto_nodes_convert(x,
                                                  nodes = gt_parent,
                                                  path = "auto")]
          # GATING METHOD
          gt_entry[, gating_method := "cyto_gate_draw"]
          # GATING ARGS
          gt_entry[, gating_args := CytoExploreR_.argDeparser(
            list(gate = gates,
                 openCyto.minEvents  = -1)
          )]
          # COLLAPSE DATA FOR GATING
          gt_entry[, collapseDataForGating := TRUE]
          # GROUP BY
          gt_entry[, groupBy := "NA"]
          # PREPROCESSING METHOD
          gt_entry[, preprocessing_method := "pp_cyto_gate_draw"]
        # OTHER GATE
        }else{
          # GATINGTEMPLATE ENTRY
          gt_entry <- gt_chunk_dt[gt_chunk_dt$alias == names(gates)[q], ]
          # ALIAS
          gt_entry[, alias := cyto_nodes_convert(x,
                                                 nodes = names(gates)[q],
                                                 anchor = gt_parent,
                                                 path = "auto")]
          # PARENT
          gt_entry[, parent := cyto_nodes_convert(x,
                                                  nodes = gt_parent,
                                                  path = "auto")]
          # GATING METHOD
          gt_entry[, gating_method := "cyto_gate_draw"]
          # GATING ARGS
          gt_entry[, gating_args := CytoExploreR_.argDeparser(
            list(gate = list(filters(list(gate))),
                 openCyto.minEvents  = -1)
          )]
          # COLLAPSE DATA FOR GATING
          gt_entry[, collapseDataForGating := TRUE]
          # GROUP BY
          gt_entry[, groupBy := "NA"]
          # PREPROCESSING METHOD
          gt_entry[, preprocessing_method := "pp_cyto_gate_draw"]
        }
        return(gt_entry)
      })
      gt_entries <- do.call("rbind", gt_entries)
      return(gt_entries)
    })
    gt_entries <- do.call("rbind", gt_entries)
    return(gt_entries)
  })
  gt_entries <- do.call("rbind", gt_entries)

  # ADD BOOLEAN ENTRIES
  if(!is.null(gt_bool_chunk)){
    lapply(seq_len(nrow(gt_bool_chunk)), function(z){
      gate_pops <- gt_bool_chunk$gating_args[z]
      gate_pops <- unlist(strsplit(gate_pops, "\\||\\&|\\!"))
      gate_pops <- gate_pops[!LAPPLY(gate_pops, ".empty")]
      ind <- LAPPLY(gate_pops, function(pop){match(pop, gt_entries$alias)})
      ind <- max(ind)
      gt_entries <<- rbind(gt_entries[1:ind, ],
                           gt_bool_chunk[z, ],
                           gt_entries[ind+1:nrow(gt_entries), , drop = TRUE])
    })
  }
  
  # REMOVE EMPTY ROWS
  gt_entries <- gt_entries[LAPPLY(seq_len(nrow(gt_entries)), function(z){
    !.all_na(gt_entries[z, ])
  }), ]
  
  # SAVE GATINGTEMPLATE
  write.csv(gt_entries, gatingTemplate, row.names = FALSE)
  
  # GATINGTEMPLATE OBJECT
  gt <- gatingTemplate(gatingTemplate, ...)
  
  # RETURN GATINGTEMPLATE
  return(gt)
  
}

#' @rdname cyto_gatingTemplate_generate
#' @export
cyto_gatingTemplate_generate.GatingSet <- function(x,
                                                   gatingTemplate = NULL,
                                                   ...){
  
  # CALL GATINGHIERARCHY METHOD
  gt <- cyto_gatingTemplate_generate(x[[1]],
                                     gatingTemplate = gatingTemplate,
                                     ...)
  
  # RETURN GATINGTEMPLATE OBJECT
  return(gt)
  
}

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
  gt <- suppressMessages(
    data_editor(gt, 
                title = "gatingTemplate Editor",
                save_as = gatingTemplate)
    )

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
