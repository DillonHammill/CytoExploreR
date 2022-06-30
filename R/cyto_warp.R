## CYTO_WARP -------------------------------------------------------------------

#' Apply normalisation algorithms to remove batch effects in cytometry data
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to extract for normalisation, set
#'   to the \code{"root"} node by default.
#' @param channels names of the channels or markers to be normalised set to all
#'   channels excluding \code{"Time"}, \code{"Event-ID"} and \code{"Sample-ID"}.
#' @param type indicates the type of normalisation to perform, native options
#'   include \code{"fdaNorm"} or \code{"warpSet"}, \code{"GaussNorm"},
#'   \code{"quantile"} and \code{"CytoNorm"}. Custom normalisation functions can
#'   also be supplied either by name or as a function.
#' @param select selection criteria to restrict the data to a subset of the
#'   samples prior to normalisation, set to \code{NULL} by default to use all
#'   samples.
#' @param group_by grouping arguments passed to \code{cyto_group_by()} to split
#'   the samples into groups prior to normalisation, set to \code{"all"} by
#'   default to normalise all samples together.
#' @param input indicates how the data should be formatted prior to passing it a
#'   custom normalisation function passed through \code{type}.
#' @param target sample selection criteria passed to \code{cyto_select()} within
#'   each group specified by \code{group_by} to identify the reference sample(s)
#'   against which the rest of the samples should be aligned.
#' @param probs vector of quantiles over which quantile normalisation should be
#'   performed, set to every percent by default.
#' @param ... additional arguments passed to the normalisation function
#'   specified by \code{type}.
#'
#' @return a \code{cytoset} or \code{GatingSet} containing the normalised data.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace cytoset
#'
#' @references Finak, G et al. (2013). High throughput flow cytometry data
#'   normalisation for clinical trials. Cytometry A 85(3).
#' @references Van Gassen, S. et al. (2019). CytoNorm: A normalisation algorithm
#'   for cytometry data. Cytometry A 97(3).
#'
#' @export
cyto_warp <- function(x,
                      parent = "root",
                      channels = NULL,
                      type = "quantile",
                      select = NULL,
                      group_by = "all",
                      input = "cytoset",
                      target = NULL,
                      probs = seq(0, 1, 0.01),
                      ...) {
  
  # DEFAULT CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x,
      exclude = c(
        "Time",
        "Event-ID",
        "Sample-ID"
      )
    )
  # PREPARE CHANNELS
  } else {
    channels <- cyto_channels_extract(
      x,
      channels = channels
    )
  }
  
  # EXCESS CHANNELS
  chans <- cyto_channels(x)
  chans <- chans[!chans %in% channels]
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(
      x,
      select
    )
  }
  
  # TARGETS
  if(!is.null(target)) {
    target <- cyto_names(x)[cyto_match(x, target)]
  }
  
  # SPLIT X INTO GROUPS
  x_split <- cyto_group_by(
    x,
    group_by = group_by
  )
  
  # ARGUMENTS
  args <- list(...)
  
  # WARPING FUNCTION - QUANTILE
  if(.grepl("quantile", type, ignore.case = TRUE)) {
    # QUANTILES
    args <- c(
      args, 
      list(probs = probs)
    )
    # INPUT
    input <- "cytoset"
  # WARPING FUNCTION - FDANORM
  } else if(.grepl("fdanorm|warpset", type, ignore.case = TRUE)) {
    # REQUIRE FLOWSTATS
    cyto_require(
      "flowStats",
      source = "BioC"
    )
    # FUN
    type <- cyto_func_match("flowStats::warpSet")
    # INPUT
    input <- "flowSet"
    # ARGUMENTS
    args[["warpFuns"]] <- FALSE
    args[["stains"]] <- channels
  # WARPING FUNCTION - GAUSSNORM
  } else if(.grepl("gaussnorm", type, ignore.case = TRUE)) {
    # REQUIRE FLOWSTATS
    cyto_require(
      "flowStats",
      source = "BioC"
    )
    # FUN
    type <- cyto_func_match("flowStats::gaussNorm")
    # INPUT
    input <- "flowSet"
    # ARGUMENTS
    args[["channel.names"]] <- channels
    # TURN OFF PLOTTING
    if(!"debug" %in% names(args)) {
      args[["debug"]] <- FALSE
    }
  # WARPING FUNCTION - CYTONORM
  } else if(.grepl("cytonorm", type, ignore.case = TRUE)) {
    # REQUIRE CYTONORM
    cyto_require(
      "CytoNorm",
      source = "GitHub",
      repo = "saeyslab/CytoNorm",
      ref = paste0(
        "Van Gassen, S. et al. (2019). CytoNorm: A normalisation algorithm ",
        "for cytometry data. Cytometry A 97(3)."
      )
    )
    # FUN
    type <- "CytoNorm"
  # WARPING FUNCTION - CUSTOM
  } else {
    type <- cyto_func_match(type)
  }
  
  # PERFORM NORMALISATION SEPARATELY ON EACH GROUP
  x_norm <- structure(
    lapply(
      x_split,
      function(z) {
        # TODO: KEEP UN-NORMALISED CHANNELS
        # EXTRACT DATA
        x_data <- cyto_data_extract(
          z,
          parent = parent,
          channels = channels,
          format = input,
          copy = FALSE
        )[[1]]
        # APPLY WARPING FUNCTION
        message(
          paste0(
            "Applying ",
            cyto_func_name(type),
            "() to normalise data..."
          )
        )
        # QUANTILE | CYTONORM
        if(!is.function(type)) {
          # QUANTILE
          if(.grepl("quantile", type, ignore.case = TRUE)) {
            # REFERENCE SAMPLES
            ref <- NULL
            if(!is.null(target)) {
              ref <- which(cyto_names(z) %in% target)
            }
            # NO REFERENCE
            if(length(ref) == 0) {
              ref <- 1:length(z)
            }
            # COMPUTE QUANTILES PER SAMPLE
            q <- cyto_apply(
              x_data,
              input = "matrix",
              FUN = "cyto_stat_quantile",
              probs = probs,
              round = 4,
              simplify = FALSE
            )
            # COMPUTE MEAN QUANTILES FOR REFERENCE GROUP
            rq <- do.call(
              "rbind",
              lapply(
                seq_along(probs),
                function(w) {
                  colMeans(
                    do.call(
                      "rbind",
                      lapply(
                        ref,
                        function(v) {
                          q[[v]][w, ]
                        }
                      )
                    )
                  )
                }
              )
            )
            rownames(rq) <- rownames(q[[1]])
            # FIT QUANTILE MODELS
            qm <- structure(
              lapply(
                seq_along(q),
                function(w) {
                  cnt <- 0
                  apply(
                    q[[w]],
                    2,
                    function(v) {
                      # CHANNEL COUNTER
                      cnt <<- cnt + 1
                      splinefun(
                        q[[w]][, cnt], # current
                        rq[, cnt], # reference
                        method = "monoH.FC"
                      )
                    }
                  )
                }
              ),
              names = cyto_names(z)
            )
            # NORMALISE QUANTILES
            cnt <- 0
            x_data_norm <- cyto_apply(
              x_data,
              input = "cytoframe",
              copy = FALSE,
              FUN = function(cf) {
                # SAMPLE COUNTER
                cnt <<- cnt + 1
                # EXTRACT DATA TO NORMALISE
                exprs <- cyto_exprs(
                  cf, 
                  channels = channels,
                  markers = FALSE,
                  drop = FALSE
                )
                # QUANTILE NORMALISATION
                cnt2 <- 0
                exprs <- apply(
                  exprs,
                  2,
                  function(w) {
                    cnt2 <<- cnt2 + 1
                    qm[[cnt]][[cnt2]](w)
                  }
                )
                # UPDATE CYTOFRAME
                cyto_exprs(cf) <- exprs
                return(cf)
              },
              simplify = FALSE
            )
          # CYTONORM
          } else if(.grepl("cytonorm", type, ignore.case = TRUE)) {
            
          }
        # OTHER WARPING FUNCTIONS
        } else {
          x_data_norm <- cyto_func_call(
            x_data,
            args
          )
        }
        # FORMAT DATA - CYTOSET -> LIST OF CYTOFRAMES
        if(cyto_class(x_data_norm, "flowSet")) {
          x_data_norm <- cyto_apply(
            x_data_norm,
            FUN = function(cf) {
              cyto_convert(cf)
            },
            input = "cytoframe",
            copy = FALSE,
            simplify = FALSE
          )
        # FORMAT DATA - LIST OF MATRICES -> LIST OF CYTOFRAMES
        } else if(cyto_class(x_data_norm, "list", TRUE)) {
          x_data_norm <- structure(
            lapply(
              x_data_norm,
              function(v) {
                if(!cyto_class(v, "cytoframe")) {
                  v <- as(v, "cytoframe")
                }
                return(v)
              }
            ),
            names = names(x_data_norm)
          )
        }
        return(x_data_norm)
      }
    ),
    names = names(x_split)
  )
  
  # CONVERT TO MEGRED NORMALISED CYTOSET
  x_norm <- cytoset(
    do.call(
      "c",
      x_norm
    )
  )

  # NORMALISED CYTOSET
  return(x_norm)
  
}

## CYTO_QUANTILE_NORM ----------------------------------------------------------

#' Apply quantile normalisation to cytometry data
#'
#' @param x object of class cytoset, GatingHierachy or GatingSet.
#' @param parent name of the parent population to normalise when a
#'   GatingHierarchy or GatingSet is supplied. A new GatingSet will be created
#'   post normalisation if the parent population is not set to \code{"root"}.
#' @param channels a vector of channels or markers over which quantile
#'   normalisation should be performed, set to all channels excluding
#'   \code{"Time"}, \code{"Sample|Event-ID"} or any dimension-reduced
#'   parameters.
#' @param select sample selection criteria passed to \code{cyto_select()} to
#'   identify the samples containing the reference quantiles to which the
#'   remaining samples should be normalised. Set to Null by default to use the
#'   average quantiles across all samples as the refernce quantiles for every
#'   sample.
#' @param probs vector of quantiles over which the data is to be normalised, set
#'   to \code{seq(0, 1, 0.01)} by default to obtain quantiles for every percent
#'   of the data.
#' @param ... not in use.
#'
#' @return a \code{cytoset}, \code{GatingHierarchy} or \code{GatingSet}
#'   containing the quantile normalised data.
#'
#' @importFrom flowWorkspace GatingSet gs_cyto_data recompute
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' library(CytoExploreRData)
#' 
#' gs <- GatingSet(Activation)
#' 
#' gs <- cyto_compensate(gs)
#' 
#' gs <- cyto_transform(gs)
#' 
#' gs <- cyto_quantile_norm(gs, channels = c("PE-A", "CD8"))
#'
#' @export
cyto_quantile_norm <- function(x,
                               parent = "root",
                               channels = NULL,
                               select = NULL,
                               probs = seq(0, 1, 0.01),
                               ...) {
  
  # TODO: APPEND NORMALISED CHANNELS TO ORIGINAL
  
  # DEFAULT CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x,
      exclude = c(
        "Event-?ID",
        "Sample-?ID",
        "Time",
        "Original",
        "UMAP",
        "t-?SNE",
        "PCA",
        "EmbedSOM",
        "FIt-?SNE"
      )
    )
  }
  
  # COMPUTE QUANTILES PER SAMPLE
  q <- cyto_apply(
    x,
    parent = parent,
    channels = channels,
    input = "matrix",
    FUN = "cyto_stat_quantile",
    probs = probs,
    round = 4,
    simplify = FALSE
  )
  
  # REFERENCE QUANTILES - MEAN OF ALL SAMPLES
  if(is.null(select)) {
    select <- 1:length(x)
  # REFERENCE QUANTILES - MEAN OF SUBSET OF SAMPLES
  } else {
    select <- cyto_match(
      x,
      select
    )
  }
  
  # COMPUTE MEAN QUANTILES FOR REFERENCE GROUP
  rq <- do.call(
    "rbind",
    lapply(
      seq_along(probs),
      function(z) {
        colMeans(
          do.call(
            "rbind",
            lapply(
              select,
              function(v) {
                q[[v]][z, ]
              }
            )
          )
        )
      }
    )
  )
  rownames(rq) <- rownames(q[[1]])
  
  # FIT QUANTILE MODELS
  qm <- structure(
    lapply(
      seq_along(q),
      function(z) {
        cnt <- 0
        apply(
          q[[z]],
          2,
          function(v) {
            # CHANNEL COUNTER
            cnt <<- cnt + 1
            splinefun(
              q[[z]][, cnt], # current
              rq[, cnt], # reference
              method = "monoH.FC"
            )
          }
        )
      }
    ),
    names = cyto_names(x)
  )
  
  # NORMALISE QUANTILES
  cnt <- 0
  cs <- cyto_apply(
    x,
    parent = parent,
    channels = channels,
    input = "cytoframe",
    copy = FALSE,
    FUN = function(cf) {
      # SAMPLE COUNTER
      cnt <<- cnt + 1
      # EXTRACT DATA TO NORMALISE
      exprs <- cyto_exprs(
        cf, 
        channels = channels,
        markers = FALSE,
        drop = FALSE
      )
      # QUANTILE NORMALISATION
      cnt2 <- 0
      exprs <- apply(
        exprs,
        2,
        function(w) {
          cnt2 <<- cnt2 + 1
          qm[[cnt]][[cnt2]](w)
        }
      )
      # UPDATE CYTOFRAME
      cyto_exprs(cf) <- exprs
      return(cf)
    }
  )
  
  # UPDATE DATA IN GATINGSET
  if(cyto_class(x, "GatingSet")) {
    # REPLACE GATINGSET ROOT NODE
    if(parent == "root") {
      gs_cyto_data(x) <- cs
      suppressMessages(recompute(x))
      return(x)
    # BUILD A NEW GATINGSET
    } else {
      # TRANSFORMATIONS
      trans <- cyto_transformers_extract(x)
      # SPILLOVER
      spill <- cyto_spillover_extract(x)
      # REVERSE TRANSFORMATIONS
      if(!.all_na(trans)) {
        cs <- cyto_transform(
          cs,
          trans = trans,
          inverse = TRUE,
          plot = FALSE,
          quiet = TRUE
        )
      }
      # REVERSE COIMPENSATION
      if(!is.null(spill)) {
        cs <- cyto_compensate(
          cs,
          spillover = spill,
          remove = TRUE
        )
      }
      # BUILD GATINGSET
      gs <- GatingSet(cs)
      # APPLY COMPENSATION
      if(!is.null(spill)) {
        cs <- cyto_compensate(
          cs,
          spillover = spill
        )
      }
      # APPLY TRANSFORMATIONS
      if(!.all_na(trans)) {
        gs <- cyto_transform(
          gs,
          trans = trans,
          inverse = FALSE,
          plot = FALSE,
          quiet = TRUE
        )
      }
      # TRANSFER GATES
      if(length(cyto_nodes(x)) > 1) {
        gs <- cyto_gateTemplate_apply(
          gs,
          x
        )
      }
      return(gs)
    }
  }
  
  # RETURN NORMALISED DATA
  return(x)
  
}