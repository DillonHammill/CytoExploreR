## CYTO_UNMIX_COMPUTE ----------------------------------------------------------

#' Compute spectral unmixing matrix
#'
#' \code{cyto_unmix_compute()} will guide the user through the steps to required
#' to compute the spectral unmixing matrix using single colour control samples.
#'
#' \code{cyto_unmix_compute()} does not perform any gating internally so samples
#' should be gated beforehand to isolate homogeneous populations of cells and/or
#' beads. \code{cyto_unmix_compute()} also requires annotation of some
#' additional information in \code{cyto_details(x)} prior to computing the
#' unmixing matrix. This includes \itemize{\item{group}{indicates the type of
#' each sample (i.e. cells or beads)}\item{parent}{name of the population to
#' extract from each sample (e.g. single cells or single beads) when a
#' GatingHierarchy or Gatingset is supplied}\item{label}{a name for the new
#' unmixed parameter for each control (e.g. CD4 FITC) - labels for unstained
#' controls should be left empty}}. Each \code{group} must contain a reference
#' unstained population in order to compute the unmixing matrix.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing single
#'   colour controls.
#' @param parent name of the population to extract from each control to compute
#'   the spectral unmixing coefficients, set to \code{"root"} by default.
#' @param channels names of the channels/markers for which unmixing coefficients
#'   should be computed, set to all fluorescent channels by default.
#' @param select passed to \code{cyto_select()} to extract a subset of samples
#'   for computation of unmixing matrix.
#' @param auto indicates which group (see \code{group_by}) should be used to get
#'   the levels of background autofluorescence.
#' @param save_as name of a CSV to which the computed unmixing matrix should be
#'   written, set to \code{Spectral-Unmixing-Matrix.csv} prefixed with the date
#'   by default. Set this argument to NA if you don't want to write the unmixing
#'   matrix to file.
#' @param ... not in use.
#'
#' @return computed unmixing coefficients in the form of a matrix.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' # GatingSet of controls
#' unmix <- cyto_unmix_compute(
#'   gs,
#'   group_by = "type",
#'   auto = "cells",
#'   save_as = "Spectral-Unmixing-Matrix.csv"
#' )
#' }
#'
#' @export
cyto_unmix_compute <- function(x,
                               parent = "root",
                               channels = NULL,
                               select = NULL,
                               auto = NULL,
                               save_as = NULL,
                               ...) {
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # SELECT
  x <- cyto_select(x, select)
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  }
  
  # UNMIXING DETAILS
  pd <- cyto_details(x)
  
  # REQUIRED VARIABLES
  vars <- c(
    "name",
    "group",
    "parent",
    "label"
  )
  
  # APPEND REQUIRED COLUMNS
  vars_req <- vars[!vars %in% colnames(pd)]
  if(length(vars_req) > 0) {
    pd <- cbind(
      pd,
      matrix(
        NA,
        ncol = length(vars_req),
        nrow = nrow(pd),
        dimnames = list(
          rownames(pd),
          vars_req
        )
      )
    )
  }
  
  # IMPORT DATA FROM FILE
  pd_new <- cyto_file_search(
    "-Unmixing-Details.*\\.csv$",
    colnames = c("name", "group", "parent", "label"),
    rownames = rownames(cyto_details(x)),
    ignore.case = TRUE,
    data.table = FALSE,
    type = "spectral unmixing details",
    files = NULL
  )
  
  # EXPERIMENT DETAILS INDICES
  pd_ind <- match(rownames(pd), rownames(pd_new[[1]]))
  
  # INHERIT DETAILS FROM FILE - FILE MAY CONTAIN EXTRA ROWS
  if(!is.null(pd_new)) {
    # APPEND EXTRA DETAILS NOT IN FILE
    vars_req <- colnames(pd)[!colnames(pd) %in% colnames(pd_new[[1]])]
    if(length(vars_req) > 0) {
      lapply(
        vars_req,
        function(z) {
          pd_new[[1]] <<- cbind(
            pd_new[[1]],
            structure(
              list(
                rep(NA, nrow(pd_new[[1]])),
                names = z
              )
            )
          )
          pd_new[[1]][pd_ind, z] <<- pd[, z]
        }
      )
    }
    pd <- pd_new[[1]]
    file <- names(pd_new)
  } else {
    file <- NULL
  }
  
  # DEFAULT PARENT - EDUCATED GUESS
  if(cyto_class(x, "GatingSet")) {
    # MISSING PARENTS - AVAILABLE SAMPLES
    ind <- which(
      rownames(pd) %in% rownames(cyto_details(x)) & 
        is.na(pd$parent)
    )
    # DEFAULT PARENT
    if(length(ind) > 0) {
      # TERMINAL NODES
      pops <- cyto_nodes(
        x[ind],
        terminal = TRUE,
        path = "auto"
      )
      # COMPUTE COUNTS FOR EACH TERMINAL NODE
      pop_stats <- cyto_apply(
        x[ind],
        parent = pops,
        channels = channels[1],
        input = "matrix",
        FUN = "cyto_stat_count",
        copy = FALSE
      )
      if(cyto_class(pop_stats, "list", TRUE)) {
        pop_stats <- do.call("cbind", pop_stats)
        dimnames(pop_stats) <- list(rownames(pop_stats), pops)
      }
      # PARENT - TERMINAL NODE MOST EVENTS
      pd$parent[ind] <- pops[
        apply(
          pop_stats,
          1,
          "which.max"
        )
      ]
    }
    # GROUPS MISSING - DEFAULT TO PARENTS
    ind <- which(is.na(pd$group))
    if(length(ind) > 0) {
      pd$group[ind] <- pd$parent[ind]
    }
  }
  
  # EDIT DETAILS - ROWNAMES CANNOT BE EDITED
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    # REMOVE ROW NAMES
    cyto_names <- rownames(pd)
    rownames(pd) <- NULL
    # PARENT - DROPDOWN COLUMN
    if(cyto_class(x, "GatingSet")) {
      col_options <- list(
        parent = cyto_nodes(x)
      )
    } else {
      col_options = NULL
    }
    # EDIT
    pd <- data_edit(
      pd,
      logo = CytoExploreR_logo(),
      title = "Spectral Unmixing Details Editor",
      row_edit = FALSE,
      # col_readonly = "name",
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane",
      col_options = col_options,
      ...
    )
    # REPLACE ROW NAMES - REQUIRED TO UPDATE CYTO_DETAILS
    rownames(pd) <- cyto_names
  }
  
  # DETAILS REQUIRE ROWNAMES FOR UPDATING - (ROWNAMES MISSING IN FILE)
  if(is.null(rownames(pd))) {
    rownames(pd) <- pd[, "name"]
  }
  
  # UPDATE DETAILS - PD MAY CONTAIN EXTRA INFORMATION
  cyto_details(x) <- pd[
    match_ind(
      rownames(
        cyto_details(x)
      ), 
      rownames(pd)
    ), 
    , 
    drop = FALSE
  ]
  
  # SAVE UNMIXING DETAILS TO FILE
  if(is.null(file)) {
    file <- cyto_file_name(
      paste0(
        format(
          Sys.Date(), 
          "%d%m%y"
        ), 
        "-Spectral-Unmixing-Details.csv"
      )
    )
  }
  
  # EXPORT UNMIXINING DETAILS
  write_to_csv(
    pd,
    file = file,
    row.names = TRUE
  )

  # RESTRICTED EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # SPLIT DATA INTO GROUPS
  cs_list <- cyto_group_by(
    x,
    group_by = "group"
  )
  
  # AUTOFLUORESCENCE
  if(is.null(auto)) {
    # SINGLE GROUP
    if(length(cs_list) == 1) {
      auto <- names(cs_list)
    # MULTIPLE GROUPS - SELECTION REQUIRED
    } else {
      # ENQUIRE
      if(interactive() & cyto_option("CytoExploreR_interactive")) {
        message(
          paste0(
            "Which of the following groups do you want to use for the ",
            "autofluorescence parameter? \n"
          )
        )
        message(
          paste0(
            1:length(cs_list), 
            ": ", 
            names(cs_list),
            sep = "\n"
          )
        )
        auto <- cyto_enquire(NULL)
        if(!is.na(suppressWarnings(as.numeric(auto)))) {
          auto <- names(cs_list)[as.numeric(auto)]
        }
        # AUTO REQUIRED
      } else {
        stop(
          paste0(
            "Supply the name of the group to use for autofluorescence ",
            "parameter to 'auto'."
          )
        )
      }
    }
  }
  
  # EXTRACT DATA - EACH GROUP SAME PARENT
  cs_list <- structure(
    lapply(
      seq_along(cs_list),
      function(z) {
        cyto_data_extract(
          cs_list[[z]],
          parent = cyto_details(cs_list[[z]])$parent[1],
          channels = channels
        )[[1]]
      }
    ),
    names = names(cs_list)
  )
  
  # COMPUTE SPECTRAL UNMIXING MATRIX -------------------------------------------
  
  # COMPUTE MEDFI - ALL CHANNELS
  medFI_list <- structure(
    lapply(
      cs_list,
      function(z) {
        do.call(
          "rbind",
          cyto_apply(
            z,
            FUN = "cyto_stat_median",
            input = "matrix",
            simplify = FALSE
          )
        )
      }
    ),
    names = names(cs_list)
  )
  
  # SWEEP OUT UNSTAINED SIGNAL PER GROUP
  unmix <- do.call(
    "rbind",
    structure(
      lapply(
        seq_along(medFI_list),
        function(z) {
          # UNSTAINED INDEX
          ind <- which.min(rowSums(medFI_list[[z]]))
          # SWEEP OUT UNSTAINED SIGNAL
          m <- sweep(
            medFI_list[[z]],
            2,
            medFI_list[[z]][ind, , drop = FALSE],
            "-"
          )[-ind, , drop = FALSE]
          # AUTOFLUORESCENCE
          if(auto %in% names(medFI_list)[z]) {
            m <- rbind(
              m, 
              "*auto*" = medFI_list[[z]][ind, ]
            )
          }
          return(m)
        }
      ),
      names = names(medFI_list)
    )
  )
  
  # PROPORTION OF SIGNAL IN EACH DETECTOR
  unmix <- t(
    apply(
      unmix,
      1,
      function(z){
        z / max(z)
      }
    )
  )
  
  # REMOVE NEGATIVE VALUES
  unmix[unmix < 0] <- 0
  
  # LABELS
  rownames(unmix) <- LAPPLY(
    rownames(unmix),
    function(z) {
      ind <- tryCatch(
        cyto_match(
          x,
          z,
          exact = TRUE
        ),
        error = function(e) {
          return(NULL)
        }
      )
      if(is.null(ind)) {
        return(
          "Autofluorescence"
        )
      } else {
        return(
          pd$label[ind]
        )
      }
    }
  )
  
  # DEFAULT FILENAME
  if(is.null(save_as)) {
    save_as <- cyto_file_name(
      paste0(
        format(
          Sys.Date(), 
          "%d%m%y"
        ), 
        "-Spectral-Unmixing-Matrix.csv"
      )
    )
  }
  
  # SAVE MATRIX
  if(!.all_na(save_as)) {
    write_to_csv(
      unmix, 
      save_as,
      row.names = TRUE
    )
  }

  # SPECTRAL UNMIXING MATRIX
  return(unmix)
  
}

## CYTO_UNMIX ------------------------------------------------------------------

#' Apply spectral unmixing matrix to spectral data
#'
#' @param x an object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing unmixed
#'   spectral data.
#' @param unmix an unmixing matrix or the name of a CSV file containing the
#'   unmixing matrix.
#' @param copy logical indicating whether unmixing should be applied to a copy
#'   of the original data, set to FALSE by default.
#' @param ... not in use.
#'
#' @return unmixed \code{cytoframe}, \code{cytoset}, \code{GatingHierarchy} or
#'   \code{GatingSet}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom stats lsfit
#'
#' @examples
#' \dontrun{
#' # Compute unmixing matrix
#' unmix <- cyto_unmix_compute(
#'   gs,
#'   group_by = "type",
#'   auto = "cells",
#'   save_as = "Spectral-Unmixing-Matrix.csv"
#' )
#' # Apply unmixing matrix
#' gs <- cyto_unmix(
#'   gs,
#'   unmix = unmix
#' )
#' }
#'
#' @name cyto_unmix
NULL

#' @noRd
#' @export
cyto_unmix <- function(x,
                       ...) {
  UseMethod("cyto_unmix")
}

#' @export
cyto_unmix.default <- function(x,
                               unmix = NULL,
                               copy = FALSE,
                               ...) {
  
  # APPLY SPECTRAL UNMIXING
  cs <- cyto_apply(
    x,
    parent = "root",
    FUN = "cyto_unmix",
    input = "cytoframe",
    unmix = unmix,
    copy = copy
  )
  
  # UPDATE DATA IN GATINGHIERARCHY/GATINGSET
  if(cyto_class(x, "GatingSet")) {
    gs_cyto_data(x) <- cs
    return(x)
  # UNMIXED CYTOSET
  } else {
    return(cs)
  }
  
}

#' @export
cyto_unmix.flowFrame <- function(x, 
                                 unmix = NULL,
                                 copy = FALSE,
                                 ...){
  
  # MISSING UNMIX
  if(is.null(unmix)) {
    stop(
      "Supply an unmixing matrix to 'unmix'."
    )
  # PREPARE UNMIX
  } else {
    # UNMIX FILE
    if(is.character(unmix)) {
      unmix <- read_from_csv(
        unmix
      )
    }
  }
  
  # CHANNELS
  channels <- cyto_channels(x)
  
  # CHANNELS UNMIX
  exprs <- cyto_exprs(
    x,
    channels = colnames(unmix),
    drop = FALSE
  )
  
  # LEAST SQUARES FIT
  ls <- lsfit(
    x = t(unmix),
    y = t(exprs),
    intercept = FALSE
  )
  
  # UNMIX
  unmix_coef <- t(ls$coefficients)
  
  # REMOVE UNMIXED PARAMETERS (EXTRA PARAMETERS REQUIRED OR CBIND NOT WORK)
  rm <- cyto_channels(x) %in% colnames(unmix)
  
  # ALL CHANNELS (CAUSES ISSUES WITH CYTO_CBIND())
  if(length(rm) == length(channels)) {
    # KEEP A CHANNEL FOR CBINDING
    x <- x[, rm[1]]
    # ADD UNMIXED PARAMETERS
    x <- cyto_cbind(
      x,
      unmix
    )
    # COPY & REMOVE EXCESS PARAMETER
    return(
      cyto_copy(
        x[, -match(rm[1], cyto_channels(x))]
      )
    )
  # SOME EXTRA CHANNELS
  } else {
    # REMOVE RAW PARAMETERS
    x <- x[, !cyto_channels(x) %in% colnames(unmix)]
    # ADD UNMIXED PARAMETERS
    x <- cyto_cbind(
      x,
      unmix_coef
    )
    # RETURN UNMIXED CYTOFRAME
    return(x)
  }

}
