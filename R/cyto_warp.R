## CYTO_WARP -------------------------------------------------------------------

#' Apply signal normalisation algorithms to cytometry data
#' 
#' @noRd
cyto_warp <- function(x,
                      type = "quantile") {
  
  
  
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
  
  # COMPUTE MEAN QUANTILES FOR REFERNCE GROUP
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