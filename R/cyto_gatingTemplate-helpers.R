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
        if(cyto_class(gate, "quadGate")){
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
  
  # SORT GATINGTEMPLATE ENTRIES
  gt_entries <- gt_entries[
    order(match(gt_entries$parent, cyto_nodes(x, path = "auto"))), ]
  
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
  write_to_csv(
    gt_entries,
    gatingTemplate
  )
  
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
  gt <- cyto_gatingTemplate_generate(
    x[[1]],
    gatingTemplate = gatingTemplate,
    ...
  )
  
  # RETURN GATINGTEMPLATE OBJECT
  return(gt)
  
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
#'
#' @return stops the gating process if an entry already exists in the
#'   gatingTemplate for the supplied alias.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_gatingTemplate_check <- function(parent,
                                       alias,
                                       gatingTemplate = NULL) {
  
  # BYPASS CHECKS - GATINGTEMPLATE DOES NOT EXIST
  gt <- tryCatch(
    suppressWarnings(
      cyto_gatingTemplate_read(gatingTemplate)
      ),
    error = function(e) {
      NULL
    }
  )
  # MATCHING PARENT & ALIAS
  if (!is.null(gt)) {
    if (parent %in% gt$parent) {
      gt_chunk <- gt[gt$parent == parent, ]
      gt_alias <- lapply(gt_chunk$alias, function(z) {
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
            "or edit gate(s) using cyto_gate_edit()."
          ))
        }
      })
    }
  }
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
  } else {
    gt <- read_from_csv(
      gatingTemplate,
      data.table = TRUE,
      ...
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
    stop("Supply the name of the csv file to 'save_as'.")
  }else{
    save_as <- file_ext_append(save_as, ".csv")
  }
  
  # READ GATINGTEMPLATE
  if(cyto_class(gatingTemplate, c("character", "gatingTemplate"))){
    gatingTemplate <- cyto_gatingTemplate_read(gatingTemplate,
                                               data.table = TRUE)
  }
  
  # HANDLE NAs
  gatingTemplate$groupBy[LAPPLY(gatingTemplate$groupBy, ".empty")] <- "NA"
  gatingTemplate$preprocessing_args[
    LAPPLY(gatingTemplate$preprocessing_args, ".empty")] <- "NA"
  
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
#' @param gatingtemplate gatingTemplate object or the name of a gatingTemplate
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
#' @importFrom openCyto CytoExploreR_.preprocess_row
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
                                      data.table = TRUE) {
  
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
      gt_chunk <- CytoExploreR_.preprocess_row(
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
          paste0(gating_args$alias,  collapse = ","), 
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
