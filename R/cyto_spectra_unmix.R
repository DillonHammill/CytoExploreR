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
#' @param parent name of the population to extract from each control to compute
#'   the spectral unmixing coefficients, set to \code{"root"} by default.
#' @param channels names of the channels/markers for which unmixing coefficients
#'   should be computed, set to all fluorescent channels by default.
#' @param select passed to \code{cyto_select()} to extract a subset of samples
#'   for computation of unmixing matrix.
#' @param auto logical indicating whether autofluorescence subtraction should be
#'   performed, set to TRUE by default. Users will be asked which group to use
#'   for autofluorescence subtraction or the name of the group can be supplied
#'   manually here.
#' @param save_as name of a CSV to which the computed unmixing matrix should be
#'   written, set to \code{Spectral-Unmixing-Matrix.csv} prefixed with the date
#'   by default. Set this argument to NA if you don't want to write the unmixing
#'   matrix to file.
#' @param axes_trans object of class \code{transformerList} containing the
#'   transformer definitions for transformations applied to the supplied data,
#'   only required when \code{cytoset} objects are supplied.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param ... additional arguments passed to \code{cyto_plot()} to allow
#'   customisation of plots used for gating.
#'
#' @return computed unmixing coefficients in the form of a matrix.
#'
#' @importFrom flowCore rectangleGate
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
                               channels = NULL,
                               select = NULL,
                               auto = TRUE,
                               save_as = NULL,
                               axes_trans = NA,
                               axes_limits = "machine",
                               ...) {

  # PREPARE DATA ---------------------------------------------------------------
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(x, select)
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
        "Sample-ID"
      )
    )
  } else {
    channels <- cyto_channels_extract(
      x,
      channels = channels
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
  
  # SPLIT DATA INTO GROUPS
  cgs_list <- cyto_group_by(
    x,
    group_by = "group"
  )
  
  # CHECK AUTOFLUORESCENCE
  if(!auto %in% c(TRUE, FALSE)) {
    if(!auto %in% names(cgs_list)) {
      auto <- NULL
    }
  } else if(auto %in% TRUE) {
    auto <- names(cgs_list)
  } else {
    auto <- NULL
  }
  
  # SELECT AUTOFLUORESCENCE GROUP
  if(length(auto) > 1) {
    # INTERCTAIVE ENQUIRE
    if(interactive() & cyto_option("CytoExploreR_interactive")) {
      message(
        paste0(
          "Which of the following groups do you want to use for the ",
          "autofluorescence parameter? \n"
        )
      )
      message(
        paste0(
          1:length(auto), 
          ": ", 
          auto,
          sep = "\n"
        )
      )
      opt <- cyto_enquire(NULL)
      if(!is.na(suppressWarnings(as.numeric(opt)))) {
        auto <- auto[as.numeric(opt)]
      } else {
        auto <- opt
      }
    # NON-INTERACTIVE 
    } else {
      auto <- auto[1]
    }
  }
  
  # COMPUTE SPECTRAL UNMIXING MATRIX -------------------------------------------
  
  # PREPARE DATA
  unmix <- do.call(
    "rbind",
    lapply(
      seq_along(cgs_list),
      function(z) {
        # CYTOSET | GATINGSET
        cgs <- cgs_list[[z]]
        # GROUP DETAILS
        pd_grp <- cyto_details(cgs)
        # LOCATE UNSTAINED CONTROLS
        unst <- NULL
        unst_ind <- which(
          grepl(
            "Unstained|NIL|Auto",
            pd_grp[, "label"],
            ignore.case = TRUE
          )
        )
        # DROP EXCESS UNSTAINED CONTRTOLS
        if(length(unst_ind) > 0) {
          if(length(unst_ind) == 1) {
            # EXTRACT UNSTAINED SAMPLE
            unst <- cgs[unst_ind]
            # DROP UNSTAINED SAMPLE
            cgs <- cgs[-unst_ind]
            # UPDATE DETAILS
            pd_grp <- pd_grp[-unst_ind, , drop = FALSE]
          } else {
            unst_rm <- unst_ind[
              -which.min(
                apply(
                  do.call(
                    "rbind",
                    lapply(
                      unst_ind,
                      function(w) {
                        # CYTOSET COPY
                        cs <- cyto_data_extract(
                          cgs[w],
                          parent = pd_grp[, "parent"][w],
                          channels = channels,
                          format = "cytoset",
                          trans = axes_trans,
                          inverse = TRUE,
                          copy = TRUE
                        )[[1]]
                        # COMPUTE CHANNEL CV
                        cyto_apply(
                          cs,
                          channels = channels,
                          input = "matrix",
                          FUN = "cyto_stat_cv",
                          trans = axes_trans,
                          inverse = FALSE,
                          copy = FALSE
                        )
                      }
                    )
                  ),
                  1,
                  "mean"
                )
              )
            ]
            # EXTRACT UNSTAINED SAMPLE
            unst <- cgs[unst_ind[!unst_ind %in% unst_rm]]
            # DROP UNSTAINED SAMPLE
            cgs <- cgs[-unst_ind]
            # UPDATE DETAILS
            pd_grp <- pd_grp[-unst_ind, , drop = FALSE]
          }
        }
        # COMPUTE PEAK EMISSION ALL DETECTORS
        peak_em <- structure(
          lapply(
            seq_along(cgs),
            function(w) {
              # UNSTAINED PEAK EM
              if(!is.null(unst)) {
                # UNSTAINED MATCH PARENT
                cs <- cyto_data_extract(
                  unst,
                  parent = pd_grp[, "parent"][w],
                  channels = channels,
                  format = "cytoset",
                  trans = axes_trans,
                  inverse = FALSE,
                  copy = FALSE
                )[[1]]
                # COMPUTE MAXIMA
                neg_em <- cyto_apply(
                  cs,
                  channels = channels,
                  input = "matrix",
                  FUN = "cyto_stat_quantile",
                  probs = c(0.99),
                  trans = axes_trans,
                  inverse = TRUE,
                  copy = TRUE
                )
              } else {
                neg_em <- NULL
              }
              # POSITIVE PEAK EM - LIENAR COPY
              cs <- cyto_data_extract(
                cgs[w],
                parent = pd_grp[, "parent"][w],
                channels = channels,
                format = "cytoset",
                trans = axes_trans,
                inverse = FALSE,
                copy = FALSE
              )[[1]]
              # COMPUTE MAXIMA
              pos_em <- cyto_apply(
                cs,
                channels = channels,
                input = "matrix",
                FUN = "cyto_stat_quantile",
                probs = 0.99,
                trans = axes_trans,
                inverse = TRUE,
                copy = TRUE
              ) 
              # RETURN NEGATIVE & POSITIVE EM
              list(
                "-" = neg_em,
                "+" = pos_em
              )
            }
          ),
          names = cyto_names(cgs)
        )
        # COMBINE PEAK EM
        peak_em <- structure(
          lapply(
            c("-", "+"),
            function(w) {
              do.call(
                "rbind",
                lapply(
                  peak_em,
                  `[[`,
                  w
                )
              )
            }
          ),
          names = c("-", "+")
        )
        # DROP DUPLICATE CONTROLS BY LABEL
        ind_rm <- LAPPLY(
          unique(pd_grp[, "label"]),
          function(w) {
            ind <- which(pd_grp[, "label"] %in% w)
            pd_chunk <- pd_grp[ind, ]
            if(nrow(pd_chunk) > 1) {
              return(
                ind[
                  -which.max(
                    rowSums(
                      peak_em[["+"]][ind, , drop = FALSE]
                    )
                  )
                ]
              )
            } else {
              return(NULL)
            }
          }
        )
        # DROP DUPLICATES
        if(length(ind_rm) > 0) {
          # REMOVE DUPLICATES FROM PEAK EM
          peak_em[["+"]] <- peak_em[["+"]][-ind_rm, ]
          peak_em[["-"]] <- peak_em[["-"]][-ind_rm, ]
          # DROP DUPLICATE SAMPLES
          cgs <- cgs[-ind_rm]
          # UPDATE DETAILS
          pd_grp <- pd_grp[-ind_rm, , drop = FALSE]
        }
        # SUBTRACT UNSTAINED PROFILE
        if(!is.null(peak_em[["-"]])) {
          # SWEEP OUT UNSTAINED
          peak_em <- peak_em[["+"]] - peak_em[["-"]]
        } else {
          peak_em <- peak_em[["+"]]
        }
        # LOCATE PEAK EMISSION CHANNEL
        pd_grp[, "peak"] <- channels[
          apply(
            peak_em,
            1,
            "which.max"
          )
        ]
        # PERFORM GATING IN PEAK EMISSION CHANNEL - POSITIVE + NEGTAIVE EVENTS
        cs_list <- structure(
          lapply(
            seq_along(cgs),
            function(w) {
              # UNSTAINED CYTOSET - TRANSFORMED COPY
              if(!is.null(unst)) {
                neg_events <- cyto_data_extract(
                  unst,
                  parent = pd_grp[w, "parent"],
                  channels = channels,
                  format = "cytoset",
                  trans = axes_trans,
                  inverse = FALSE,
                  copy = FALSE
                )[[1]]
              } else {
                neg_events <- NULL
              }
              # STAINED CYTOSET - TRANSFORMED COPY
              pos_events <- cyto_data_extract(
                cgs[w],
                parent = pd_grp[w, "parent"],
                channels = channels,
                format = "cytoset",
                trans = axes_trans,
                inverse = FALSE,
                copy = FALSE
              )[[1]]
              # GATING
              cyto_plot(
                if(length(neg_events) == 0){
                  pos_events
                } else {
                  neg_events
                },
                channels = pd_grp[w, "peak"],
                overlay = if(!is.null(neg_events)){
                  pos_events
                } else {
                  NA
                },
                hist_stack = 0,
                hist_fill = if(is.null(neg_events)){
                  "dodgerblue"
                } else {
                  c("red", "dodgerblue")
                },
                hist_fill_alpha = 0.6,
                title = cyto_names(pos_events),
                axes_limits = axes_limits, 
                axes_trans = axes_trans,
                legend = FALSE,
                ...
              )
              # INTERACTIVE - USE CYTO_GATE_DRAW()
              if(interactive() & cyto_option("CytoExploreR_interactive")) {
                # GATE NEGATIVE POPULATION
                neg_gt <- cyto_gate_draw(
                  x = if(!is.null(neg_events)){
                    neg_events
                  } else {
                    pos_events
                  },
                  alias = paste0(pd_grp[w, "peak"], "-"),
                  channels = pd_grp[w, "peak"],
                  type = "interval",
                  plot = FALSE
                )
                # GATE POSITIVE POPULATION
                pos_gt <- cyto_gate_draw(
                  x = pos_events,
                  alias = paste0(pd_grp[w, "peak"], "+"),
                  channels = pd_grp[w, "peak"],
                  type = "interval",
                  plot = FALSE
                )
              # NON-INTERACTIVE - USE QUANTILES
              } else {
                # NEGATIVE QUANTILE GATE
                neg_gt <- cyto_apply(
                  if(is.null(neg_events)) {
                    pos_events
                  } else {
                    neg_events
                  },
                  channel = pd_grp[w, "peak"],
                  input = "matrix",
                  FUN = "cyto_stat_quantile",
                  probs = c(0.01, 0.4),
                  trans = axes_trans,
                  inverse = FALSE,
                  copy = FALSE
                )
                rownames(neg_gt) <- c("min", "max")
                # CREATE GATE
                neg_gt <- rectangleGate(
                  neg_gt,
                  filterId = paste0(
                    pd_grp[w, "peak"],
                    "-"
                  )
                )
                # NEGATIVE QUANTILE GATE
                pos_gt <- cyto_apply(
                  pos_events,
                  channel = pd_grp[w, "peak"],
                  input = "matrix",
                  FUN = "cyto_stat_quantile",
                  probs = c(0.95, 0.99),
                  trans = axes_trans,
                  inverse = FALSE,
                  copy = FALSE
                )
                rownames(pos_gt) <- c("min", "max")
                # CREATE GATE
                pos_gt <- rectangleGate(
                  pos_gt,
                  filterId = paste0(
                    pd_grp[w, "peak"],
                    "+"
                  )
                )
                # PLOT GATES
                cyto_plot_gate(
                  list(
                    neg_gt,
                    pos_gt
                  )
                )
              }
              # RETURN GATED POSITIVE/NEGATIVE POPULATIONS
              return(
                list(
                  "-" = cyto_gate_apply(
                    if(is.null(neg_events)) {
                      pos_events
                    } else {
                      neg_events
                    },
                    neg_gt
                  )[[1]][[1]], # LIST OF POPULATIONS PER GATE
                  "+" = cyto_gate_apply(
                    pos_events,
                    pos_gt
                  )[[1]][[1]] # LIST OF POPULATIONS PER GATE
                )
              )
            }
          ),
          names = cyto_names(cgs)
        )
        # COMPUTE UNMIXING MATRIX
        do.call(
          "rbind",
          structure(
            lapply(
              seq_along(cs_list),
              function(q) {
                # COMPUTE POSITIVE MEDFI
                pos <- cyto_apply(
                  cs_list[[q]][["+"]],
                  FUN = "cyto_stat_median",
                  channels = channels,
                  input = "matrix",
                  trans = axes_trans,
                  inverse = TRUE,
                  copy = TRUE
                )
                rownames(pos) <- pd_grp[q, "label"]
                # COMPUTE NEGATIVE MEDFI
                neg <- cyto_apply(
                  cs_list[[q]][["-"]],
                  FUN = "cyto_stat_median",
                  channels = channels,
                  input = "matrix",
                  trans = axes_trans,
                  inverse = TRUE,
                  copy = TRUE
                )
                rownames(neg) <- pd_grp[q, "label"]
                # SWEEP OUT NEGATIVE MEDFI & APPEND AUTO
                if(any(auto %in% names(cgs_list)[z]) & q == length(cs_list)) {
                  return(
                    rbind(
                      pos - neg,
                      matrix(
                        apply(
                          neg,
                          2,
                          "mean"
                        ),
                        nrow = 1,
                        ncol = length(channels),
                        dimnames = list(
                          "Autofluorescence",
                          channels
                        )
                      )
                    )
                  )
                } else {
                  return(
                    pos - neg
                  )
                }
              }
            ),
            names = names(cs_list)
          )
        )
      }
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
  
  utils::write.csv(
    unmix_coef,
    "columns.csv"
  )
  
  # REMOVE UNMIXED PARAMETERS (EXTRA PARAMETERS REQUIRED OR CBIND NOT WORK)
  rm <- channels[
    channels %in% colnames(unmix)
  ]
  
  print(rm)
  
  # ALL CHANNELS (CAUSES ISSUES WITH CYTO_CBIND())
  if(length(rm) == length(channels)) {
    # KEEP A CHANNEL FOR CBINDING
    x <- x[, rm[1]]
    # ADD UNMIXED PARAMETERS
    x <- cyto_cbind(
      x,
      unmix_coef
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
    if(length(rm) > 0) {
      x <- x[, -match(rm, colnames(unmix))]
    }
    # ADD UNMIXED PARAMETERS
    x <- cyto_cbind(
      x,
      unmix_coef
    )
    # RETURN UNMIXED CYTOFRAME
    return(x)
  }

}
