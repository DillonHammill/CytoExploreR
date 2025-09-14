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
#' controls should be set to \code{"Unstained"}}}. Each \code{group} must
#' contain a reference unstained population in order to compute the unmixing
#' matrix.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing single
#'   colour controls.
#' @param parent  name of the population to use for the unmixing calculation
#'   when a GatingSet object is supplied, set to the last node of the GatingSet
#'   by default (e.g. "Single Cells"). For greater flexibility, users can
#'   specify a parent population for each control, which will be extracted for
#'   the unmixing calculation (e.g. Lymphocytes for CD4 APC or Myeloid Cells for
#'   CD11b FITC). The parent populations for each control can also be specified
#'   in a \code{parent} column in the channel match CSV file or in
#'   \code{cyto_details}.
#' @param select passed to \code{\link{cyto_select}} to select the samples
#'   required to compute the spillover matrix.
#' @param channels names of the channels/markers for which unmixing coefficients
#'   should be computed, set to all fluorescent channels by default.
#' @param gate indicates whether to \code{draw} or use \code{auto} gating
#'   methods to gate the negative and positive populations for each control.
#'   Gating is optional for \code{"Autospill"} and \code{"CytoDecode"} methods
#'   but is recommended to achieve better results and ensure compatibility with
#'   \code{cyto_panel_design()}. Set to \code{"draw"} by default to allow for
#'   more flexibility over gating. For \code{gate = "draw"} users can also
#'   specify whether to gate the 1D histograms or 2D scatter plots 1D or 2D
#'   suffix to \code{gate}. For example, to gate in 2D set \code{gate =
#'   "draw-2D"} which is the new default manual gating method.
#' @param type options include \code{"Bagwell" or "autocomp"}, \code{"Roca" or
#'   "autospill"} or \code{"CytoDecode"} to indicate which method to use when
#'   computing the spillover matrix, set to \code{"CytoDecode"} by default. 
#' @param save_as name of a CSV to which the computed unmixing matrix should be
#'   written, set to \code{Spectral-Unmixing-Matrix.csv} prefixed with the date
#'   by default. Set this argument to NA if you don't want to write the unmixing
#'   matrix to file.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'   \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'   saved, set to \code{cyto_gatingTemplate_active()} by default.
#' @param axes_trans object of class \code{transformerList} containing the
#'   transformer definitions for transformations applied to the supplied data,
#'   only required when \code{cytoset} objects are supplied.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param heatmap logical indicating whether the computed spectral unmixing
#'   matrix should be displayed in a heatmap, set to TRUE by default.
#' @param events number of events to extract from each control prior to fitting
#'   RLM models, set to 500 events by default.
#' @param model indicates the type of model to use when running autospill,
#'   options include \code{"rlm"} for robust linear models or \code{"vwqr"} for
#'   variance weighted quantile regression, set to \code{"rlm"} by default.
#' @param max_iter indicates the maximum number of allowable iterations for
#'   refining the spillover coefficients, set to 20 by default.
#' @param trim proportion of events to exclude from the top and bottom of each
#'   scale when fitting models in \code{AutoSpill} and \code{CytoDecode}
#'   methods, set to 0.001 by default.
#' @param tol tolerance level for convergence of vwqr model in \code{CytoDecode}
#'   method, set to \code{1e-6} by default.
#' @param grid_size size of the grid to use for wvqr model in \code{CytoDecode}
#'   method, set to 100 by default. Smaller grid sizes will provide a speed
#'   boost at the cost of accuracy.
#' @param resid indicates whether to use \code{"lower"}, \code{"upper"} or
#'   \code{"both"} residuals when fitting variance models in vwqr used in
#'   \code{CytoDecode}, set to \code{"both"} by default.
#' @param details name of a CSV file to which the details of the compensation
#'   controls should be saved, set to NULL by default to use
#'   \code{date-Control-Details.csv}. Setting this argument to \code{NA}
#'   will prevent details from being written to a CSV file.
#' @param ... additional arguments passed to \code{cyto_plot()} to allow
#'   customisation of plots used for gating.
#'
#' @return unmixing matrix and write unmixing matrix to csv file named in
#'   accordance with \code{save_as}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' # GatingSet of controls
#' unmix <- cyto_unmix_compute(
#'   gs,
#'   auto = "cells",
#'   save_as = "Spectral-Unmixing-Matrix.csv"
#' )
#' }
#'
#' @export
cyto_unmix_compute <- function(x,
                               parent = "root",
                               select = NULL,
                               channels = NULL,
                               gate = "draw-2d",
                               type = "CytoDecode",
                               save_as = NULL,
                               gatingTemplate = NULL,
                               axes_trans = NA,
                               axes_limits = "machine",
                               heatmap = TRUE,
                               events = 500,
                               model = "rlm",
                               max_iter = 20,
                               trim = 0.001,
                               tol = 1e-6,
                               grid_size = 100,
                               resid = "both",
                               details = NULL,
                               ...) {
 
  # PULL DOWN ARGUMENTS
  args <- CytoExploreR:::.args_list(...)
  args[["unmix"]] <- TRUE
  
  # COMPUTE UNMIXING MATRIX
  unmix <- cyto_func_call(
    "cyto_spillover_compute",
    args
  )
  
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
#' @param parent name of the parental population in the supplied
#'   \code{GatingHierarchy} or \code{GatingSet} to unmix, set to \code{"root"}
#'   by default.
#' @param unmix an unmixing matrix or the name of a CSV file containing the
#'   unmixing matrix.
#' @param auto either \code{FALSE} to subtract autofluorescence prior to
#'   unmixing or \code{TRUE} to include autofluorescence spectrum in the
#'   unmixing matrix, set to \code{FALSE} by default. Autofluorescence is
#'   automatically computed for each parent population using the unstained or
#'   autofluorescence control in each sample group.
#' @param type method to use for unmxing, options include \code{"OLS"}.
#' @param select sample selection criteria to select a subset of samples for
#'   unmixing, see \code{\link{cyto_select}} for details.
#' @param drop logical indicating whether the original channels used compute the
#'   unmixing matrix should be dropped from the unmixed data, set to TRUE by
#'   default.
#' @param trans object of class \code{transformerList} containing the
#'   definitions of transformers already applied to the input \code{cytoframe}
#'   or \code{cytoset}, set to NA by default.
#' @param inverse logical to indicate whether inverse data transformations
#'   should be applied to the data prior to unmixing, set to TRUE by default.
#' @param save_as name of the output directory to use when writing new unmixed
#'   FCS files, must be manually supplied by the user.
#' @param ... additional arguments passed to \code{cyto_save} such as
#'   \code{name} to specify how the unmixed data should be re-written to new FCS
#'   files.
#'
#' @return unmixed \code{cytoframe}, \code{cytoset}, \code{GatingHierarchy} or
#'   \code{GatingSet} with fluorescent values on the linear scale.
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
                               parent = "root",
                               unmix = NULL,
                               auto = FALSE,
                               type = "OLS",
                               select = NULL,
                               drop = TRUE,
                               trans = NA,
                               inverse = TRUE,
                               save_as = NULL,
                               ...) {
  
  # TODO: ENFORCE SAVE_AS - AVOID TEMPFILES SO METADATA MATCHES
  
  # TODO: DEFAULT PARENT IS LEAF NODES OF GATINGSET
  
  # TODO: AUTOFLUORESCENCE IS NOT STORED IN UNMIXING MATRIX
  
  # NOTE: WE RETURN LINEAR UNMIXED DATA + LINEAR RAW DATA (CYTOSET | GATINGSET)
  # NOTE: WE CAN'T KEEP TRANSFORMERS AS UNMIXED PARAMETERS WILL ALSO BE TRANSFORMED
  # DOWNSTREAM - SO WE DON'T TRANSFER GATES
  
  # MISSING UNMIX
  if(is.null(unmix)) {
    stop(
      "Supply an unmixing matrix to 'unmix'."
    )
    # PREPARE UNMIX
  } else {
    # UNMIX FILE
    if(is.character(unmix)) {
      # ROWNAMES REQUIRED
      unmix <- read_from_csv(
        unmix,
        data.table = FALSE
      )
    }
  }
  
  # SAVE_AS REQUIRED
  if(!is.null(save_as)) {
    if(!dir.exists(save_as)) {
      dir.create(save_as)
    }
  }
  
  # EXTRACT TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # CHECK CHANNELS
  channels <- colnames(unmix)
  rm_idx <- which(
    !channels %in% cyto_channels(x)
  )
  if(length(rm_idx) > 0) {
    message(
      "'unmix' references channels absent in 'x': \n",
      paste(
        colnames(unmix)[rm_idx],
        collapse = "\n"
      )
    )
  }
  
  # CHECK PARENT
  if(cyto_class(x, "GatingSet")) {
    parent <- cyto_nodes_convert(
      x,
      nodes = parent,
      path = "auto"
    )
  }

  # MULTI-PARENT ONLY SUPPORTED FOR GATINGSETS
  if(cyto_class(x, "flowSet")) {
    parent <- "root"
  }
  names(parent) <- parent
  
  # SAMPLE SELECTION
  if(!is.null(select)) {
    x <- cyto_select(
      x,
      select
    )
  }
  
  # METADATA
  pd <- cyto_details(x)
  
  # UNMIXING REQUIRES GROUP VARIABLE
  if(any(!c("group", "stain") %in% colnames(pd))) {
    if(!"group" %in% colnames(pd)) {
      pData(x) <- cbind(
        pData(x),
        "group" =  rep(NA, length(x))
      )
    }
    # UNMIXING REQUIRES STAIN VARIABLE
    if(!"stain" %in% colnames(pd)) {
      pData(x) <- cbind(
        pData(x),
        "stain" =  rep(NA, length(x))
      )
    }
    # UPDATE METADATA
    message(
      "Update 'group' and 'stain' variables to unmix data..."
    )
    cyto_details_edit(
      x
    )
  }
  
  # SPLIT DATA INTO GROUPS
  cgs_list <- cyto_group_by(
    x,
    group_by = "group"
  )
  
  # UNMIXED PARAMETERS
  chans <- gsub(
    "^<(.*)>(.*)", 
    "\\2",
    rownames(unmix)
  )
  chans <- trimws(chans, "both")
  
  # UNMIXING PROGRESS BAR
  pb <- cyto_progress(
    label = "cyto_unmix",
    total = length(x),
    clear = FALSE
  )
  
  # UNMIX EACH PARENT IN EACH GROUP
  cf_list <- lapply(
    cgs_list,
    function(cgs) {
      # GROUP METADATA
      pd <- cyto_details(cgs)
      # CHECK FOR PARENT VARIABLE
      if("parent" %in% colnames(pd)) {
        parent <- unique(
          pd$parent
        )
        # MULTIPLE PARENTS - COMMA SEPARATED - NOT DOCUMENTED
        parent <- unlist(
          strsplit(
            parent, 
            ","
          )
        )
        names(parent) <- parent
      }
      # LOCATE GROUP UNSTAINED
      unst_idx <- grep(
        "unst|auto|nil",
        pd$stain,
        ignore.case = TRUE
      )
      # NOTE: MULTIPLE UNSTAINED WILL BE COERCED
      # ESTIMATE AUTOFLUORESCENCE FOR EACH PARENT
      if(length(unst_idx) == 0) {
        AF <- NULL
      } else {
        AF <- do.call(
          "rbind",
          lapply(
            parent,
            function(pop) {
              # COMPUTE MEDFI ACROSS DETECTORS - LINEAR SCALE
              res <- cyto_apply(
                cgs[unst_idx],
                parent = parent,
                channels = channels,
                FUN = "cyto_stat_median",
                input = "matrix",
                coerce = if(length(unst_idx) > 1) {
                  TRUE
                } else {
                  FALSE
                },
                barcode = FALSE,
                trans = trans,
                inverse = TRUE,
                copy = TRUE
              )
              if(length(parent) > 1) {
                res <- do.call(
                  "rbind",
                  res
                )
              }
              rownames(res) <- parent
              return(res)
            }
          )
        )
      }
      # UNMIX EACH PARENT for EACH SAMPLE
      idx <- structure(
        seq_along(cgs), 
        names = gsub(
          "\\.fcs",
          "_unmixed.fcs",
          cyto_names(cgs)
        )
      )
      # RETURN A LIST OF UNMIXED CYTOFRAMES
      res <- lapply(
        idx,
        function(id) {
          # UNMIX EACH PARENT CYTOFRAME THEN MERGE
          # NOTE: COERCED DATA READ FROM TEMPFILE SO NAME VARIABLE WON'T MATCH
          cf <- as(
            cytoset(
              lapply(
                parent,
                function(pop) {
                  # UPDATE UNMIXING MATRIX TO INCLUDE AUTOFLUORESCENCE
                  if(!is.null(AF)) {
                    unmix <- rbind(
                      unmix,
                      "Autofluorescence" = as.numeric(
                        AF[
                          match(pop, rownames(AF)),
                          colnames(unmix), 
                          drop = TRUE
                        ]
                      )
                    )
                  }
                  # SPECTRAL UNMIXING
                  cf <- cyto_unmix(
                    cyto_data_extract(
                      cgs[id],
                      parent = pop,
                      channels = NULL,
                      format = "cytoframe",
                      trans = trans,
                      inverse = TRUE,
                      copy = TRUE
                    )[[1]][[1]],
                    unmix = unmix,
                    auto = auto,
                    type = type,
                    drop = drop,
                    copy = FALSE,
                    trans = NA,
                    save_as = NULL
                  )
                  # SET NEW FILENAME
                  cyto_keyword(
                    cf,
                    "GUID",
                    gsub(
                      "\\.fcs",
                      "_unmixed.fcs",
                      cyto_names(cgs)[id]
                    )
                  )
                  cyto_keyword(
                    cf,
                    "$FIL",
                    gsub(
                      "\\.fcs",
                      "_unmixed.fcs",
                      cyto_names(cgs)[id]
                    )
                  )
                  # UPDATE PNTYPE
                  params <- pData(parameters(cf))
                  params <- rownames(params)[
                    match(
                      chans,
                      params$name
                    )
                  ]
                  for(param in params) {
                    cyto_keyword(
                      cf,
                      paste0(param, "TYPE"),
                      "Unmixed_Fluorescence"
                    )
                  }
                  # INCREMENT PROGRESS BAR
                  cyto_progress(pb)
                  # RETURN UNMIXED CYTOFRAME
                  return(cf)
                }
              )
            ),
            "cytoframe"
            # if(!is.null(save_as)) {
            #   paste0(save_as, "/", names(id))
            # } else{
            #   save_as
            # },
            # ...
          )
          # RETURN MERGED CYTOFRAME
          return(cf)
        }
      )
      # RETURN UNMIXED DATA
      return(res)
    }
  )
  
  # UNMIXED CYTOSET
  cs <- cytoset(
    do.call(
      "c", 
      unname(
        cf_list
      )
    )
  )
  
  # SAVE_AS
  if(!is.null(save_as)) {
    cyto_save(
      cs,
      save_as = save_as,
      names = cyto_names(cs)
    )
  }
  
  # INHERIT METDATA
  pd <- cyto_details(x)[
    match(
      gsub(
        "_unmixed\\.fcs$",
        ".fcs",
        cyto_names(cs)
      ), 
      cyto_names(x)
    ),
    ,
    drop = FALSE
  ]
  pd$name <- cyto_details(cs)$name
  rownames(pd) <- cyto_names(cs)
  cyto_details(cs) <- pd
  
  # PREPARE LINEAR GATINGSET - ALL DATA SAME SCALE
  if(cyto_class(x, "GatingSet")) {
    cs <- GatingSet(cs)
    if(cyto_class(x, "GatingHierarchy")) {
      cs <- cs[[1]]
    }
  }
  
  # RETURN UNMIXED DATA
  return(cs)
  
}

#' @export
cyto_unmix.flowFrame <- function(x, 
                                 unmix = NULL,
                                 auto = FALSE,
                                 type = "OLS",
                                 drop = TRUE,
                                 trans = NA,
                                 inverse = TRUE,
                                 save_as = NULL,
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
      # ROWNAMES REQUIRED
      unmix <- read_from_csv(
        unmix,
        data.table = FALSE
      )
    }
  }
  
  # SAVE_AS REQUIRED
  if(!is.null(save_as)) {
    if(!dir.exists(save_as)) {
      dir.create(save_as)
    }
  }
  
  # CHANNELS
  channels <- cyto_channels(x)
  
  # INVERSE TRANSFORMATIONS - ALL CHANNELS ON LINEAR SCALE
  if(!.all_na(trans) & inverse) {
    trans <- cyto_transformers_combine(
      trans[names(trans) %in% colnames(unmix)]
    )
    x <- cyto_transform(
      x,
      trans = trans,
      inverse = TRUE,
      copy = TRUE
    )
  }
  
  # CHANNELS UNMIX
  exprs <- cyto_exprs(
    x,
    channels = colnames(unmix),
    drop = FALSE
  )
  
  # SUBTRACT AUTOFLUORESCENCE
  auto_idx <- grep(
    "Autofluorescence",
    rownames(unmix),
    ignore.case = TRUE
  )
  if(length(auto_idx) > 0 & !auto) {
    # NOTE: ONLY ONE AF SPECTRUM SUPPORTED
    exprs <- sweep(
      exprs,
      2,
      as.numeric(
        unmix[
          auto_idx[1], 
          colnames(exprs), 
          drop = TRUE
        ]
      ),
      "-"
    )
    # DROP AUTOFLUORESCENCE SPECTRA
    unmix <- unmix[
      -auto_idx, , 
      drop = FALSE
    ]
  }
  
  # ROW SUM-TO-ONE
  unmix <- sweep(
    unmix,
    1,
    rowSums(unmix),
    "/"
  )
  
  # NON-NEGATIVE LEAST SQUARES UNMIXING
  if(grepl("^n", type, ignore.case = TRUE)) {
    # NNLS PACKAGE REQUIRED
    cyto_require(
      "nnls"
    )
    # UNMIX IN PARALLEL
    unmix_coef <- t(
      future_apply(
        exprs,
        1, 
        function(z) {
          cyto_func_call(
            "nnls::nnls",
            args = list(
              rbind(
                t(unmix),
                rep(1, nrow(unmix)) # weight = 1 for sum-to-one
              ),
              c(z, 1) # weight = 1 for sum-to-one
            )
          )$x
        }
      )
    )
  # FULLY CONSTRAINED LEAST SQUARES
  } else if(grepl("^f", type, ignore.case = TRUE)) {
    # QUADPROG PACKAGE REQUIRED
    cyto_require(
      "quadprog"
    )
    # PREAPRE MATRICES
    Dmat <- 2 * tcrossprod(unmix)
    Amat <- t(
      rbind(
        rep(1, nrow(unmix)),
        diag(nrow(unmix))
      )
    )
    bvec <- c(1, rep(0, nrow(unmix)))
    unmix_coef <- t(
      future_apply(
        exprs,
        1, 
        function(z) {
          dvec <- 2 * tcrossprod(z, unmix)
          cyto_func_call(
            "quadprog::solve.QP",
            args = list(
              Dmat,
              dvec,
              Amat,
              bvec,
              meq = 1 # first contraint = equality
            )
          )$solution
        }
      )
    )
  # SOLVE - FAST
  } else if(grepl("^s", type, ignore.case = TRUE)) {
    unmix_coef <- t(
      solve(
        tcrossprod(unmix),
        t(
          tcrossprod(
            exprs,
            unmix
          )
        )
      )
    )
  # ORDINARY LEAST SQUARES
  } else {
    unmix_coef <- t(
      lsfit(
        x = t(unmix),
        y = t(exprs),
        intercept = FALSE
      )$coefficients
    )
  }
  
  # SET ENDMEMBER NAMES
  colnames(unmix_coef) <- rownames(unmix)
  
  # PREPARE NEW PARAMETER NAMES
  mrks <- gsub("^<(.*)>(.*)", "\\1", rownames(unmix))
  chans <- gsub("^<(.*)>(.*)", "\\2", rownames(unmix))
  chans <- trimws(chans, "both")
  mrks <- trimws(mrks, "both")
  mrks[mrks == chans] <- NA
  
  # SET NEW CHANNEL NAMES
  colnames(unmix_coef) <- chans
  rownames(unmix_coef) <- NULL
  
  # APPEND UNMIXED PARAMETERS
  if(!drop) {
    # ADD UNMIXED PARAMETERS
    x <- cyto_cbind(
      x,
      unmix_coef
    )
  # DROP ORIGINAL PARAMETERS
  } else {
    rm <- colnames(unmix)
    # REMOVE ALL CHANNELS
    if(length(rm) == length(channels)) {
      # KEEP A CHANNEL FOR CBINDING
      x <- x[, rm[1]]
      # ADD UNMIXED PARAMETERS
      x <- cyto_cbind(
        x,
        unmix_coef
      )
      # DROP EXCESS PARAMETERS
      x <- cyto_copy(
        x[, -match(rm[1], cyto_channels(x))]
      )
    # REMOVE SOME CHANNELS
    } else {
      # ADD UNMIXED PARAMETERS
      x <- cyto_cbind(
        x[, -match(rm, cyto_channels(x))],
        unmix_coef
      )
    }
  }
  
  # SET MARKER NAMES
  if(!.all_na(mrks)) {
    idx <- which(!is.na(mrks))
    mrks <- mrks[idx]
    names(mrks) <- chans[idx]
    cyto_markers(x) <- mrks
  }
  
  # SET PNTYPE FOR UNMIXED PARAMETERS
  pd <- pData(parameters(x))
  lapply(
    chans,
    function(channel) {
      # PARAMETER IDENTIFIER
      id <- rownames(
        pd[pd$name %in% channel, , drop = FALSE]
      )
      # SET PNR KEYWORD
      cyto_keyword(
        x,
        paste0(id, "TYPE"),
        "Unmixed_Fluorescence"
      )
    }
  )
  
  # UPDATE FILENAME KEYWORD
  cyto_keyword(
    x,
    "$FIL",
    gsub(
      "\\.fcs",
      "_unmixed.fcs",
      cyto_keyword(x, "$FIL")
    )
  )

  # UPDATE GUID
  id <- cyto_keyword(
    x,
    c("GUID", "$GUID")
  )
  id <- id[!sapply(id, "is.null")]
  cyto_keyword(
    x,
    names(id),
    gsub(
      "\\.fcs",
      "_unmixed.fcs",
      id[[1]]
    )
  )
  
  # WRITE UNMIXED CYTOFRAME
  if(!is.null(save_as)) {
    cyto_save(
      x,
      save_as = save_as,
      ...
    )
  }

  # RETURN UNMIXED CYTOFRAME
  return(x)

}
