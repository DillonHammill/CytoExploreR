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
#'   the unmixing calculation (e.g. Lymphocytes for CD4 APC or Myeloid Cells
#'   for CD11b FITC). The parent populations for each control can also be
#'   specified in a \code{parent} column in the channel match CSV file or in
#'   \code{cyto_details}.
#' @param select passed to \code{\link{cyto_select}} to select the samples
#'   required to compute the spillover matrix.
#' @param channels names of the channels/markers for which unmixing coefficients
#'   should be computed, set to all fluorescent channels by default.
#' @param type options include \code{"Bagwell"}, \code{"Roca"} or
#'   \code{"hybrid"} to indicate which method to use when computing the
#'   spillover matrix, set to \code{"Roca"} by default. The \code{"hybrid"}
#'   method computes the spillover coefficients using the \code{Bagwell}
#'   approach (no RLM) and refines the coefficients using the \code{Autospill}
#'   approach. Refer to \code{references} section for more details about each
#'   method.
#' @param auto logical indicating whether autofluorescence subtraction should be
#'   performed, set to TRUE by default. Users will be asked which group to use
#'   for autofluorescence subtraction or the name of the group can be supplied
#'   manually here.
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
#' @param iter indicates the maximum number of allowable iterations for refining
#'   the spillover coefficients, set to 100 by default.
#' @param trim proportion of events to exclude from the top and bottom of each
#'   scale when fitting RLM models in autospill method, set to 0.001 by default.
#' @param details name of a CSV file to which the details of the compensation
#'   controls should be saved, set to NULL by default to use
#'   \code{date-Compensation-Details.csv}. Setting this argument to \code{NA}
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
                               type = "roca",
                               auto = FALSE,
                               save_as = NULL,
                               gatingTemplate = NULL,
                               axes_trans = NA,
                               axes_limits = "machine",
                               heatmap = TRUE,
                               events = 500,
                               iter = 100,
                               trim = 0.001,
                               details = NULL,
                               ...) {
 
  # TODO: LOOSEN MATCHING ON UNSTAINED INCLUDE NA | NIL | AUTO
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(
      x, 
      select
    )
  }
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x,
      exclude = c(
        "FSC",
        "SSC",
        "Time",
        "Event-ID",
        "Sample-ID",
        "\\-H$",
        "\\-W$"
      )
    )
  } else {
    channels <- cyto_channels_extract(
      x,
      channels = channels
    )
  }
  
  # MATCH CHANNELS - LABEL VARIABLE REQUIRED
  pd <- cyto_channel_match(
    x,
    channels = channels,
    label = TRUE,
    save_as = details
  )
  pd <- pd[
    match(
      rownames(cyto_details(x)),
      rownames(pd)
    ), 
    ,
    drop = FALSE
  ]
  
  # CHECK COMPENSATION DETAILS
  lapply(
    seq_len(nrow(pd)),
    function(z) {
      # CHANNEL
      if(!pd$channel[z] %in% c("Unstained",
                               "unstained",
                               cyto_channels(x))) {
        stop(
          paste0(
            pd$channel[z], 
            " is not a valid channel for this ", 
            cyto_class(x), 
            "!"
          )
        )
        # PARENT
        if(cyto_class(x, "GatingSet")) {
          cyto_nodes_convert(x, pd$parent[z]) # ERRORS IF MISSING
        }
      }
    }
  )
  
  # MISSING VARIABLES
  vars <- c("name", "group", "parent", "channel", "label")
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
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # TRANSFORMED DATA REQUIRED FOR GATING
  if(any(!channels %in% names(axes_trans))) {
    # DEFINE NEW TRANSFORMERS
    trans_new <- cyto_transformers_define(
      x, 
      channels = channels[!channels %in% names(axes_trans)],
      type = "biex",
      widthBasis = -10,
      plot = FALSE,
      progress = FALSE
    )
    # APPLY NEW TRANSFORMERS
    x <- suppressWarnings(
      cyto_transform(
        x,
        trans = trans_new,
        copy = TRUE,
        plot = FALSE,
        quiet = TRUE
      )
    )
    # COMBINE TRANSFROMERS
    if(.all_na(axes_trans)) {
      axes_trans <- trans_new
    } else {
      axes_trans <- cyto_transformers_combine(axes_trans, trans_new)
    }
  }
  
  # RESTRICT TRANSFORMERS TO CHANNELS
  if(!.all_na(axes_trans)) {
    if(any(channels %in% names(axes_trans))) {
      axes_trans <- cyto_transformers_combine(
        axes_trans[names(axes_trans) %in% channels]
      )
    } else {
      axes_trans <- NA
    }
  }

  # SPLIT DATA INTO GROUPS
  groups <- unique(pd$group)
  
  # CHECK AUTOFLUORESCENCE
  if(!all(auto %in% c(TRUE, FALSE))) {
    auto <- auto[
      auto %in% groups
    ]
    if(length(auto) == 0) {
      warning(
        "'auto' must be the name(s) of a sample group as in ",
        "cyto_details(x)$group!"
      )
    }
  } else if(all(auto %in% TRUE)) {
    auto <- groups
  } else {
    auto <- NULL
  }
  
  # COMPUTE SPECTRAL UNMIXING MATRIX -------------------------------------------
  
  # STAINED INDICES IN X
  idx <- which(
    !grepl(
      "Unstained",
      pd$channel,
      ignore.case = TRUE
    )
  )
  
  # LOOP THROUGH STAINED CONTROLS
  unmix <- do.call(
    "rbind",
    lapply(
      idx, 
      function(id) {
        # LOCATE UNSTAINED CONTROL(S)
        unst_idx <- which(
          grepl(
            "Unstained",
            pd$channel,
            ignore.case = TRUE
          ) &
            pd$group %in% pd$group[id]
        )
        # COMPUTE SPILLOVER COEFFICIENTS
        spill <- cyto_spillover_compute(
          x[c(unst_idx, id)],
          channels = channels,
          type = type,
          save_as = NA,
          gatingTemplate = gatingTemplate,
          axes_trans = axes_trans,
          axes_limits = axes_limits,
          heatmap = FALSE,
          events = events,
          iter = iter,
          trim = trim,
          details = NA, # DON'T EXPORT DETAILS - DONE ALREADY
          unmix = TRUE,
          ...
        )
        # EXTRACT NON-EMPTY ROW
        spill <- spill[
          apply(
            spill, 
            1,
            function(z) {
              !all(z %in% c(0, 1))
            }
          ),
          , 
          drop = FALSE
        ]
        # SET LABEL AS ROWNAMES
        rownames(spill) <- pd$label[id]
        
        # RETURN COMPUTED COEFFICIENTS
        return(spill)
      }
    )
  )
  
  # AUTOFLUORESCENCE
  if(length(auto) > 0) {
    # INDICES OF AUTOFLUORESCENCE CONTROLS
    auto_idx <- which(
      pd$group %in% auto &
        grepl("Unstained", pd$channel, ignore.case = TRUE)
    )
    # USE FIRST UNSTAINED CONTROL PER CGROUP
    auto_idx <- auto_idx[!is.duplicated(pd$group[auto_idx])]
    # COMPUTE AUTOFLUORESCENCE SPECTRA
    auto <- do.call(
      "rbind",
      lapply(
        auto_idx,
        function(id) {
          # EXTRACT CYTOSET
          auto_cs <- cyto_data_extract(
            x[id],
            parent =  pd$parent[id],
            channels = channels,
            format = "cytoset",
            inverse = TRUE,
            copy = TRUE,
            trans = axes_trans
          )[[1]]
          # COMPUTE MEDIAN IN ALL LINEAR CHANNELS 
          auto_spec <- cyto_apply(
            auto_cs,
            FUN = "cyto_stat_median",
            channels = channels,
            input = "matrix",
            simplify = TRUE
          )
          # SWEEP OUT MAX
          auto_spec <- auto_spec/max(auto_spec)
          # SET ROWNAMES
          rownames(auto_spec) <- paste0(
            "<",
            basename(pd$parent[id]),
            "> ",
            "Autofluorescence"
          )
          # RETURN AUTOFLUORESCENCE SPECTRA
          return(auto_spec)
        }
      )
    )
    # APPEND AUTOFLUORESCENCE TO UNMXING MATRIX
    unmix <- rbind(unmix, auto)
  }
  
  # REMOVE NEGATIVE VALUES
  unmix[unmix < 0] <- 0
  
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

  # HEATMAP
  if(heatmap) {
    um <- unmix
    # CONSTRUCT HEATMAP - USE CYTO_PLOT_HEATMAP?
    heat_map(
      um,
      tree_y = TRUE,
      cell_col_scale = .cyto_plot_point_col_scale(),
      title = "Spectral Unmixing Matrix"
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
#' @param drop logical indicating whether the original channels used compute the
#'   unmixing matrix should be dropped from the unmixed data, set to TRUE by
#'   default.
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
                               drop = TRUE,
                               ...) {
  
  # APPLY SPECTRAL UNMIXING
  cs <- cyto_apply(
    x,
    parent = "root",
    FUN = "cyto_unmix",
    input = "cytoframe",
    unmix = unmix,
    copy = copy,
    drop = drop
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
                                 drop = TRUE,
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
  colnames(unmix_coef) <- rownames(unmix)
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
      # DROP EXCESS PARAMETER
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
  
  # RETURN UNMIXED CYTOFRAME
  return(x)

}
