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
  
  # REQUIRED COLUMNS
  vars_req <- c()
  lapply(
    vars,
    function(z) {
      if(any(grepl(z, colnames(pd), ignore.case = TRUE))) {
        colnames(pd)[grep(z, colnames(pd), ignore.case = TRUE)[1]] <<- z
      } else {
        vars_req <<- c(vars_req, z)
      }
    }
  )
  
  # FILE MUST REFLECT ANY CHANGES MADE TO DATA
  # SEARCH FOR UNMIXING DETAILS
  pd_new <- cyto_file_search(
    "-Unmixing-Details.csv$",
    colnames = c("name", "group", "parent", "label"),
    rownames = rownames(cyto_details(x)),
  )
  # NO UNMIXING DETAILS FOUND
  if(length(pd_new) == 0) {
    # UNMIXING DETAILS FILE
    file <- NULL
    # PREPARE MISSING VARIABLES
    pd_new <- matrix(
      NA,
      nrow = nrow(pd),
      ncol = length(vars_req),
      dimnames = list(
        rownames(pd),
        vars_req
      )
    )
    # SET DEFAULT FOR PARENT
    if("parent" %in% colnames(pd_new)) {
      pd_new[, "parent"] <- rep(parent, length.out = nrow(pd_new))
    }
    # APPEND MISSING VARIABLES
    pd <- cbind(
      pd, 
      pd_new
    )
  # UNMIXING DETAILS FOUND
  } else {
    # MULTIPLE COMPENSATION DETAILS FILES FOUND
    if(length(pd_new) > 1) {
      # ENQUIRE
      if(interactive() & cyto_option("CytoExploreR_interactive")) {
        message(
          "Multiple files found with spectral unmixing details for this ",
          cyto_class(x),
          ". Which file would you like to import unmixing details from?"
        )
        message(
          paste0(
            paste0(
              1:length(pd_new),
              ": ",
              names(pd_new)
            ),
            sep = "\n"
          )
        )
        opt <- cyto_enquire(NULL)
        opt <- tryCatch(
          as.numeric(opt),
          warning = function(w) {
            return(
              match(opt, names(pd_new))
            )
          }
        )
        # COMPENSATION DETAILS FILE
        pd_new <- pd_new[opt]
      } else {
        pd_new <- pd_new[1]
      }
    }
    # COMPENSATION DETAILS FILE
    file <- names(pd_new)
    message(
      paste0(
        "Importing saved unmixing details from ",
        file,
        "..."
      )
    )
    # TODO: DO WE WANT TO APPEND ALL DETAILS HERE?
    # APPEND MISSING VARIABLES
    pd <- cbind(
      pd,
      pd_new[[file]][
        match_ind(
          rownames(pd),
          rownames(pd_new[[file]])
        ), vars_req]
    )
  }
  
  # UPDATE EXPERIMENT DETAILS
  cyto_details(x) <- pd
  
  # INTERACTIVELY EDIT UNMIXING DETAILS
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    # EDIT UNMIXING DETAILS - PARENT DROPDOWN
    x <- cyto_details_edit(
      x,
      file = NA,
      col_options = if(cyto_class(x, "GatingSet")) {
        list(
          "parent" = cyto_nodes(x, path = "auto")
        )
      } else {
        NULL
      }
    )
    # UPDATE UNMIXING DETAILS
    pd <- cyto_details(x)
  }
  
  # MISSING VARIABLES
  if(!all(vars %in% colnames(pd))) {
    stop(
      paste0(
        "cyto_details(x) is missing required variables: ",
        paste0(
          vars[!vars %in% colnames(pd)],
          collapse = " & "
        )
      )
    )
  }
  
  # FILE - SAVE UNMIXING DETAILS
  if(is.null(file)) {
    file <- paste0(
      format(
        Sys.Date(), 
        "%d%m%y"
      ), 
      "-Spectral-Unmixing-Details.csv"
    )
    pd_new <- pd
  # UPDATE DETAILS IN EXISTING FILE
  } else {
    # ORIGINAL DETAILS STORED IN PD_NEW[[FILE]] - NEW DETAILS STORED IN PD
    # ADD NEW COLUMNS TO ORIGINAL DETAILS
    ind <- which(!colnames(pd) %in% colnames(pd_new[[file]]))
    if(length(ind) > 0) {
      pd_new[[file]] <- do.call(
        "cbind",
        c(
          list(
            pd_new[[file]]
          ),
          structure(
            lapply(
              ind,
              function(z) {
                rep(NA, nrow(pd_new[[file]]))
              }
            ),
            names = colnames(pd)[ind]
          )
        )
      )
    }
    # UPDATE ROWS
    ind <- match_ind(rownames(pd), rownames(pd_new[[file]]))
    pd_new[[file]][ind, colnames(pd)] <- pd
    pd_new <- pd_new[[file]]
  }
  
  # EXPORT UPDATED UNMIXING DETAILS
  write_to_csv(
    pd_new,
    file,
    row.names = TRUE
  )
  
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
          auto <- names(cs_list)[auto]
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
    save_as <- paste0(
      format(
        Sys.Date(), 
        "%d%m%y"
      ), 
      "-Spectral-Unmixing-Matrix.csv"
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
