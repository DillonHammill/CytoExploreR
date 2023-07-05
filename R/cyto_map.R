## CYTO_MAP --------------------------------------------------------------------

#' Create dimension-reduced maps of cytometry data
#'
#' \code{cyto_map} is a convenient wrapper function to produced
#' dimension-reduced maps o cytometry data using any available
#' dimension-reduction algorithm. \code{cyto_map} comes with native support for
#' PCA, t-SNE, FIt-SNE, UMAP and EmbedSOM. Simply supply the name of the
#' dimension reduction function to the \code{type} argument and additional
#' arguments directly to \code{cyto_map}. \code{cyto_map} takes care of merging
#' the data prior to generating the dimension-reduced maps and automatically
#' splits the data into the original samples for downstream analyses.
#'
#' @param x object of class \code{data.frame}, \code{matrix},
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to extract from the GatingSet for
#'   mapping, set to the \code{"root"} node by default.
#' @param select passed to \code{cyto_select()} to control which samples are
#'   passed to the dimension reduction algorithm, set to all samples by default.
#' @param channels names of the channels or markers to feed into the dimension
#'   reduction algorithm, set to all channels with markers assigned by default.
#' @param type indicates the type of dimension reduction algorithm to apply to
#'   the data, can be either supplied by name \code{e.g. "UMAP"} or as a
#'   function \code{e.g. "uwot::umap"}. Natively supported options include:
#'   \itemize{\item{"PCA"}{ uses rsvd::rpca() to compute principal
#'   components}\item{"t-SNE"}{ uses python openTSNE or R Rtsne::Rtsne() to
#'   compute t-SNE co-ordinates} \item{"FIt-SNE"}{ uses
#'   \url{https://github.com/KlugerLab/FIt-SNE} to compute FIt-SNE co-ordinates
#'   (some additional configuration required)}\item{"UMAP"}{ uses python
#'   umap-learn or R uwot::umap() to compute UMAP co-ordinates}\item{"UMATO"}{
#'   uses UMATO python library to compute dimension-reduced co-ordinates using
#'   hub points and nearest neighbour graphs}\item{"NCViz"}{uses ncvis python
#'   module to compute a noise contrastive dimenion-reduced
#'   co-ordinates}\item{"EmbedSOM"}{ uses EmbedSOM::SOM() and
#'   EmbedSOM::EmbedSOM() to compute EmbedSOM co-ordinates}\item{"PHATE"}{ uses
#'   python phate to compute trajectory embedding co-ordinates}\item{"PaCMAP"}{
#'   uses python pacmap to compute pairwise controlled manifold approximation
#'   co-ordinates}\item{"TriMap"}{ uses python trimap.TRIMAP() to compute
#'   triplet constrained co-ordinates}\item{"IsoMap"}{ uses python
#'   sklearn.manifold.Isomap() to compute non-linear isomapping
#'   co-ordinates}\item{"MDS"}{ uses python sklearn.manifold.MDS() or
#'   \code{stats::cmdscale()} to compute multidimensional scaling
#'   co-ordinates}\item{"kNN"}{ builds a k nearest neighbours (kNN) graph using
#'   RANN::nn2() and uses igraph to compute the layout co-ordinates for the
#'   graph}\item{"sNN"}{ builds a shared nearest neighbours (sNN) graph using
#'   HGC::SNN.construction() and igraph to compute the layout co-ordinates for
#'   the graph}\item{"MST"}{ builds a minimum spanning tree using the igraph 
#'   package}}
#' @param scale optional argument to scale each channel prior to computing
#'   dimension-reduced co-ordinates, options include \code{"global"},
#'   \code{"range"}, \code{"mean"}, \code{"median"} or \code{"zscore"}. Set to
#'   \code{"range"} by default, scaling can be turned off by setting this
#'   argument to FALSE. Global scaling preserves the structure of the data as
#'   expected on the linear scale without having to apply inverse data
#'   transformations to the entire dataset. Scaling is only required when some
#'   channels have been transformed and \code{inverse = FALSE}.
#' @param events number or proportion of events to map, can optionally be
#'   supplied per cytoframe for more fine control over how the data is sampled
#'   prior to mapping. Set to 1 by default to map all events.
#' @param label optional label for new dimension-reduced parameters \code{e.g.
#'   "UMAP"} to which the parameter indices will be appended \code{e.g. "UMAP-1"
#'   and "UMAP-2"}.
#' @param merge_by passed to \code{cyto_merge_by} to control how samples should
#'   be merged prior to performing dimension-reduction, set to \code{"all"} by
#'   default to create a single dimension-reduced map for all samples.
#'   Alternatively, a different dimension-reduced map will be created for each
#'   group of samples merged by \code{cyto-merge_by}.
#' @param nodes a vector of population names representing descendant clustered
#'   populations of the specified \code{parent}. If supplied, a vector of
#'   cluster labels will be generated and passed to \code{Haisu} to provide
#'   improved embedding using any of the natively supported non-linear dimension
#'   reduction algorithms.
#' @param seed numeric passed to \code{cyto_data_extract} to ensure consistent
#'   sampling between runs, set to NULL by default to off this feature.
#' @param inverse logical to indicate whether the data should be inverse
#'   transformed prior to performing dimension reduction, set to FALSE by
#'   default.
#' @param trans object of class \code{transformerList} containing the
#'   transformation definitions applied to the supplied data. Used internally
#'   when \code{inverse} is TRUE to inverse these transformations prior to apply
#'   dimension reduction algorithm.
#' @param plot logical indicating whether a call should be made to
#'   \code{cyto_plot_map} to visualise the produced dimension-reduced maps, set
#'   to TRUE by default.
#' @param ... additional arguments passed to the specified dimension reduction
#'   algorithm
#'
#' @return either a \code{cytoset} or \code{GatingSet} with the
#'   dimension-reduced parameters added as additional channels in each of the
#'   underlying \code{cytoframes}.
#'
#' @importFrom flowWorkspace cytoset gs_cyto_data recompute GatingSet
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_data_extract}}
#' @seealso \code{\link{cyto_save}}
#' @seealso \code{\link{cyto_plot_map}}
#'
#' @export
cyto_map <- function(x,
                     parent = "root",
                     select = NULL,
                     channels = NULL,
                     type = "UMAP",
                     scale = "range",
                     events = 1,
                     label = NULL,
                     merge_by = "all",
                     nodes = NULL,
                     seed = NULL,
                     inverse = FALSE,
                     trans = NA,
                     plot = TRUE,
                     ...) {
  
  # TODO: FIX SAMPLING USING EVENTS ARGUMENT
  
  # CHECKS ---------------------------------------------------------------------
  
  # CLASS
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    x <- as(x, "cytoset")
  }
  
  # SOM CHECK
  som <- cyto_keyword(
    x[[1]],
    "CytoExploreR_SOM"
  )[[1]]
  # SOM GRID DIMENSIONS
  grid <- NULL
  if(length(som) > 0) {
    grid <- as.numeric(
      strsplit(
        som,
        "x"
      )[[1]]
    )
  }
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x,
      exclude = c(
        "Time",
        "Original",
        "Sample",
        "Event",
        "PCA",
        "t-?SNE",
        "FIt-?SNE",
        "UMAP",
        "EmbedSOM"
      )
    )
    channels <- channels[
      channels %in% names(cyto_markers(x))
    ]
  }
  
  # CONVERT CHANNELS
  channels <- cyto_channels_extract(
    x,
    channels = channels,
    plot = FALSE
  )
  
  # TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # SELECT
  x <- cyto_select(
    x,
    select
  )
  
  # COPY - REQUIRED?
  x <- cyto_copy(x)
  
  # NAMES - ALL SELECTED SAMPLES (MAY BE EMPTY)
  x_names <- cyto_names(x)
  
  # COUNTS PER SAMPLE (NEED TO REMOVE ZERO EVENT SAMPLES)
  x_counts <- cyto_apply(
    x,
    parent = parent,
    FUN = "nrow",
    input = "matrix",
    channels = channels[1], # DON'T LOAD ENTIRE MATRIX INTO MEMORY
    copy = FALSE
  )[, 1]
  
  # EXCLUDE ZERO EVENT SAMPLES
  x_data <- cyto_select(
    x,
    names(x_counts)[x_counts != 0]
  )
  
  # EMBEDDING PER GROUP
  x_data_groups <- cyto_group_by(
    x_data,
    group_by = merge_by
  )
  
  # NAMES FOR EXPORTED FILES - SPLIT PER GROUP
  names <- structure(
    lapply(
      x_data_groups, 
      "cyto_names"
    ),
    names = names(x_data_groups)
  )
  
  # LABELS
  if(is.null(label)) {
    if(is.function(type)) {
      label <- tail(
        as.character(
          substitute(
            type
          )
        ), 
        1
      )
    } else {
      # NAMESPACED FUNCTION
      if(grepl(":{2,3}", type)) {
        label <- tail(
          unlist(
            strsplit(
              type, 
              ":{2,3}"
            )
          ),
          1
        )
      } else {
        label <- type
      }
    }
  }
  
  # NODES
  if(!is.null(nodes) & !cyto_class(x, "GatingSet")) {
    nodes <- NULL
    warning(
      "'nodes' is only supported for GatingSet objects!"
    )
  }
  
  # MAP EACH GROUP
  x_map_chans <- c()
  x_data_map <- structure(
    lapply(
      seq_along(x_data_groups), 
      function(z) {
        # CYTO_DATA TO MAP
        cyto_data <- x_data_groups[[z]]
        # CYTOFRAME CONTAINER
        cs <- cyto_data_extract(
          cyto_data,
          parent = parent,
          format = "cytoset",
          coerce = FALSE,
          events = events,
          barcode = TRUE,
          overwrite = TRUE,
          seed = seed,
          trans = trans,
          inverse = inverse,
          copy = FALSE
        )[[1]]
        # COERCE - SAMPLING ABOVE
        cf <- cyto_coerce(
          # SOM - USE FIRST SAMPLE ONLY
          if(length(som) > 0) {
            cs[1]
          } else {
            cs
          },
          format = "cytoframe"
        )
        # NODES (HAISU FUTURE SUPPORT) - MATCH SAMPLED EVENTS
        if(!is.null(nodes)) {
          # PYTHON UNAVAILABLE
          if(!cyto_func_call("reticulate::py_available")) {
            nodes <- NULL
            warning(
              paste0(
                "python must be available to use Haisu for dimensionality ",
                "reduction!"
              )
            )
          # HAISU PYTHON
          } else {
            # CLUSTER LABELS
            nodes <- unlist(
              cyto_gate_indices(
                cyto_data,
                parent = parent,
                nodes = nodes,
                labels = TRUE
              )
            )
            # ALL EVENT IDS
            ids <- cyto_data_extract(
              cyto_data,
              parent = parent,
              format = "matrix",
              channels = "Event-ID",
              coerce = TRUE,
              trans = trans,
              inverse = inverse,
              copy = FALSE
            )[[1]][[1]]
            # RESTRICT NODES TO SAMPLED IDS
            nodes <- nodes[
              ids %in% cyto_exprs(
                cf, 
                channels = "Event-ID",
                drop = TRUE
              )
            ]
            rm(ids)
          }
        }
        # RAW DATA
        cf_exprs <- cyto_exprs(
          cf,
          channels = channels,
          drop = FALSE
        )
        # SCALE TRANSFORMED DATA
        if(!scale %in% FALSE) {
          message(
            paste0(
              "Performing ",
              scale,
              " scaling on channels..."
            )
          )
          # GLOBAL SCALING
          if(scale %in% "global") {
            # LINEAR CHANNEL LIMITS
            cnt <- 0
            limits <- apply(
              cf_exprs,
              2,
              function(z) {
                cnt <<- cnt + 1
                rng <- range(z)
                # INVERSE TRANSFORM -> LINEAR SCALE
                if(colnames(cf_exprs)[cnt] %in% names(trans)) {
                  rng <- .cyto_transform(
                    rng,
                    channel = colnames(cf_exprs)[cnt],
                    trans = trans,
                    inverse = TRUE
                  )
                }
                return(rng)
              }
            )
            # GLOBAL LIMITS ON LINEAR SCALE
            limits <- range(limits)
            # APPLY GLOBAL SCALING
            cnt <- 0
            cf_exprs <- apply(
              cf_exprs,
              2,
              function(z) {
                cnt <<- cnt + 1
                rng <- limits
                # GLOBAL LIMITS ON TRANSFORMED SCALE
                if(colnames(cf_exprs)[cnt] %in% names(trans)) {
                  rng <- .cyto_transform(
                    rng,
                    channel = colnames(cf_exprs)[cnt],
                    trans = trans,
                    inverse = FALSE
                  )
                }
                # RESCALE DATA
                cyto_stat_rescale(
                  z,
                  scale = c(0,1),
                  limits = rng
                )
              }
            )
          # CONVENTIONAL SCALING
          } else {
            cf_exprs <- cyto_stat_scale(
              cf_exprs,
              type = scale
            )
          }
        }
        # PERFORM MAPPING
        coords <- .cyto_map(
          cf_exprs,
          type = type,
          seed = seed,
          label = label,
          nodes = nodes,
          grid = grid,
          ...
        )
        x_map_chans <<- colnames(coords)
        # REPLACING AN EXISTING EMBEDDING
        if(all(colnames(coords) %in% cyto_channels(cf))) {
          cf_exprs <- cyto_exprs(cf)
          cf_exprs[, colnames(coords)] <- coords
          cyto_exprs(cf) <- cf_exprs
        # APPEND A NEW EMBEDDING 
        } else {
          cf <- cyto_cbind(
            cf, 
            coords
          )
        }
        # SOM - CYTOSET
        if(length(som) > 0) {
          cs <- cyto_cbind(
            cs,
            do.call(
              "rbind",
              rep(
                list(coords),
                length(cs)
              )
            )
          )
        # SPLIT - CYTOSET
        } else {
          cs <- cyto_split(
            cf,
            names = names[[z]]
          )
        }
        return(cs)
      }
    ),
    names = names(x_data_groups)
  )
  
  # COMBINE CYTOSETS
  if(length(x_data_map) >  1) {
    x_data_map <- cyto_convert(
      cyto_func_call(
        "rbind2",
        x_data_map
      )
    )
  } else {
    x_data_map <- x_data_map[[1]]
  }

  # CHANNELS - INCLUDE MAPPED CHANNELS
  x_data_map_chans <- cyto_channels(x_data_map)
  
  # MAPPED NON-EMPTY FILE NAMES
  x_data_map_names <- cyto_names(x_data_map)
    
  # COMPLETE MAPPED DATA
  x_data_map <- cytoset(
    structure(
      lapply(
        seq_along(x_names),
        function(z) {
          # CYTOFRAME 
          if(x_names[z] %in% x_data_map_names) {
            cf <- cyto_convert(
              x_data_map[[match(x_names[z], x_data_map_names)]]
            )
          # EMPTY CYTOFRAME
          } else {
            cf <- cyto_empty(
              x_names[z],
              x_data_map_chans
            )
          }
          return(cf)
        }
      ),
      names = x_names
    )
  )
  
  # INHERIT EXPERIMENT DETAILS
  cyto_details(x_data_map) <- cyto_details(x)[
    match(
      rownames(cyto_details(x_data_map)),
      rownames(cyto_details(x))
    ), , drop = FALSE
  ]
  
  # REBUILD GATINGSET
  if(cyto_class(x, "GatingSet")) {
    # REBUILD GATINGSET
    x <- cyto_rebuild(
      x = x,
      y = x_data_map,
      trans = trans,
      inverse = inverse
    )
  # CYTOSET
  } else {
    x <- x_data_map
  }
  
  # PLOT 
  if(plot) {
    # TODO: OVERLAY GATED DESCENDANTS (SLOW & POINT_COL REQUIRES EVENTS=1)
    cyto_plot(
      x,
      parent = parent,
      channels = x_map_chans[1:2],
      merge_by = "all",
      point_shape  = if(length(som) > 0) {
        21
      } else {
        "."
      }
    )
  }
  
  # RETURN FORMATTED DATA
  return(x)
  
}

#' Internal function to perform dimension reduction
#' @importFrom utils tail
#' @noRd
.cyto_map <- function(x,
                      type = "UMAP",
                      seed = NULL,
                      label = NULL,
                      nodes = NULL,
                      grid = NULL,
                      ...){
  
  # SEED
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # ARGUMENTS
  args <- .args_list(...)
  args <- args[
    -match(
      c("type",
        "seed",
        "label"),
      names(args)
    )
  ]
  
  # NODES
  if(is.null(nodes)) {
    args <- args[-match("nodes", names(args))]
  # HAISU ONLY SUPPORTED TSNE | UMAP | PHATE
  } else {
    # REMOVE NODES FROM ARGUMENTS
    if(is.character(type)) {
      if(!grepl("^UMAP$|^t-?SNE$|^PHATE$", type, ignore.case = TRUE)) {
        args <- args[-match("nodes", names(args))]
        stop(
          "'nodes' can only be used with 'type' set to UMAP, t-SNE or PHATE!"
        )
      }
    }
  }
  
  # LABEL
  if(is.null(label)) {
    if(is.function(type)) {
      label <- tail(
        as.character(
          substitute(
            type
          )
        ), 
        1
      )
    } else {
      # NAMESPACED FUNCTION
      if(grepl(":{2,3}", type)) {
        label <- tail(
          unlist(
            strsplit(
              type, 
              ":{2,3}"
            )
          ),
          1
        )
      } else {
        label <- type
      }
    }
  }
  
  # CHARACTER FUNCTION
  if(is.character(type)) {
    # PCA
    if(grepl("^PCA$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_pca"
    # FIt-SNE
    } else if(grepl("^FIt-?SNE$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_fitsne"
    # T-SNE
    } else if(grepl("^t-?SNE$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_tsne"
    # UMAP 
    } else if(grepl("^UMAP$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_umap"
    # UMATO 
    } else if(grepl("^UMATO$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_umato"
    # EMBEDSOM
    } else if(grepl("^Embed-?SOM$", type , ignore.case = TRUE)) {
      type <- ".cyto_map_embedsom"
    # PHATE
    } else if(grepl("^PHATE$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_phate"
    # PACMAP
    } else if(grepl("^PaCMAP$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_pacmap"
    # TRIMAP
    } else if(grepl("^TriMap$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_trimap"
    # ISOMAP
    } else if(grepl("^IsoMap$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_isomap"
    # MDS
    } else if(grepl("^MDS$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_mds"
    # KNN
    } else if(grepl("^kNN$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_knn"
    # SNN
    } else if(grepl("^SNN$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_snn"
    # NCVis
    }else if(grepl("^NCVis$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_ncvis"
    # MST
    } else if(grepl("^MST$", type, ignore.case = TRUE)) {
      type <- ".cyto_map_mst"
    }
    # MST SOM GRID ARGUMENT
    if(!type %in% ".cyto_map_mst") {
      args <- args[
        -match(
          "grid",
          names(args)
        )
      ]
    }
  }
  
  # MESSAGE
  message(
    paste0(
      "Computing ",
      paste0(label, "()"),
      " co-ordinates..."
    )
  )
  # CALL MAPPING FUNCTION
  cyto_map_coords <- cyto_func_call(
    type,
    args
  )
  # EXTRACT CO-ORDINATES (ALL DIMENSIONS)
  if(!cyto_class(cyto_map_coords, "matrix")) {
    # FIND CO-ORDINATES
    coords <- NULL
    lapply(
      seq_along(cyto_map_coords), 
      function(z) {
        if(cyto_class(cyto_map_coords[[z]], "matrix")) {
          if(nrow(cyto_map_coords[[z]]) == nrow(x)) {
            coords <<- cbind(coords, cyto_map_coords[[z]])
          }
        }
      }
    )
    cyto_map_coords <- coords
  }
  
  # DIMENSION REDUCTION
  colnames(cyto_map_coords) <- paste0(
    label,
    "-",
    1:ncol(cyto_map_coords)
  )
  
  # RETURN MAPPED DATA
  return(cyto_map_coords)
  
}

# CROSS ENTROPY TEST SUPPORT REMOVED
# cross_entropy <- function(x,
#                           type = "UMAP",
#                           coords = NULL) {
# 
#   # MESSAGE
#   message(
#     paste0(
#       "Computing ",
#       type,
#       " cross-entropies..."
#     )
#   )
# 
#   # REFERENCE
#   message(
#     paste0(
#       "Roca C. et al. (2021) A Cross Entropy test allows quantitative ",
#       "statistical comparison of t-SNE and UMAP representations. ",
#       "arXiv:2112.04172"
#     )
#   )
# 
#   # SIZE
#   n <- nrow(x)
# 
#   # UMAP -----------------------------------------------------------------------
#   if(type == "UMAP") {
# 
#     # FORMAT NEAREST NEIGHBOURS - ORIGINAL SPACE
#     cnt <- 0
#     coords$nn[[1]]$idx <- t(
#       apply(
#         coords$nn[[1]]$idx,
#         1,
#         function(z) {
#           cnt <<- cnt + 1
#           z[z != cnt]
#         }
#       )
#     )
#     coords$nn[[1]]$dist <- t(
#       apply(
#         coords$nn[[1]]$dist,
#         1,
#         function(z) {
#           z[z != 0]
#         }
#       )
#     )
# 
#     # COMPUTE SIGMA
#     sigma <- apply(
#       coords$nn[[1]]$dist,
#       1,
#       function(dd){
#         # UMAP SIGMA ERROR
#         sigma_error <- function(ss, dd) {
#           p <- exp(-pmax(0, dd - min(dd)) / ss)
#           sum(p) - log2(length(p))
#         }
#         dd_ascen <- sort(dd[dd > 0])
#         ss_lower <- dd_ascen[1]
#         dd_descen <- sort(dd[!is.infinite(dd)], decreasing = TRUE)
#         ss_upper <- dd_descen[1]
#         while(sigma_error(ss_upper, dd) < 0){
#           ss_lower <- ss_upper
#           ss_upper <- 2 * ss_upper
#         }
#         while(sigma_error(ss_lower, dd) > 0) {
#           ss_upper <- ss_lower
#           ss_lower <- ss_lower/2
#         }
#         uniroot(
#           sigma_error,
#           dd,
#           interval = c(ss_lower, ss_upper),
#           tol = (ss_upper - ss_lower) * .Machine$double.eps ^ 0.25
#         )$root
#       }
#     )
# 
#     # COMPUTE PROBABILITIES - ORIGINAL SPACE
#     cnt <- 0
#     orig_p <- t(
#       apply(
#         coords$nn[[1]]$dist,
#         1,
#         function(z) {
#           cnt <<- cnt + 1
#           exp(
#             -pmax(
#               0,
#               z - min(z)
#             ) /
#               sigma[cnt]
#           )
#         }
#       )
#     )
# 
#     # SYMMETRIZE PROBABILITIES - ORIGINAL SPACE
#     for(i in seq_len(nrow(coords$nn[[1]]$idx))) {
#       for(j2 in seq_len(length(coords$nn[[1]]$idx[i, ]))) {
#         j <- coords$nn[[1]]$idx[i, j2]
#         i2 <- match(i, coords$nn[[1]]$idx[j, ])
#         if(!is.na(i2)) {
#           if(j > i) {
#             p_sym <- orig_p[i, j2] + orig_p[j, i2] -
#               orig_p[i, j2] * orig_p[j, i2]
#             orig_p[i, j2] <- p_sym
#             orig_p[j, i2] <- p_sym
#           }
#         } else {
#           # orig_p[i, j2] <- orig_p[i, j2]/2
#         }
#       }
#     }
# 
#     # COMPUTE DISTANCE IN UMAP SPACE FOR ORIGINAL CLOSE NEIGHBOURS
#     cnt <- 0
#     umap_dist <- t(
#       apply(
#         coords$nn[[1]]$idx,
#         1,
#         function(z) {
#           cnt <<- cnt + 1
#           unlist(
#             lapply(
#               z,
#               function(w) {
#                 sum(
#                   (coords$embedding[cnt, ] - coords$embedding[w, ]) ^ 2
#                 )
#               }
#             )
#           )
#         }
#       )
#     )
# 
#     # COMPUTE PROBABILITIES ASSOCIATED TO UMAP REPRESENTATION
#     umap_p <- t(
#       apply(
#         umap_dist,
#         1,
#         function(z) {
#           1/ (1 + coords$a * z ^ coords$b)
#         }
#       )
#     )
# 
#     # COMPUTE FUZZY CROSS ENTROPIES
#     res <- unlist(
#       lapply(
#         1:nrow(orig_p),
#         function(z) {
#           -sum(
#             orig_p[z, ] * log(umap_p[z, ]) +
#               (1 - orig_p[z, ]) * log(1 - umap_p[z, ])
#           )
#         }
#       )
#     )
# 
#   # t-SNE --------------------------------------------------------------------
#   } else if(type == "t-SNE" ) {
# 
#     # RANN REQUIRED - NN2()
#     cyto_require(
#       "RANN",
#       source = "CRAN"
#     )
# 
#     # COMPUTE NEAREST NEIGHBOURS - ORIGINAL SPACE
#     coords$nn <- cyto_func_call(
#       "RANN::nn2",
#       args = list(
#         x,
#         k = coords$perplexity * 3 + 1
#       )
#     )
#     names(coords$nn) <- c("idx", "dist")
# 
#     # FORMAT NEAREST NEIGHBOURS - ORIGINAL SPACE
#     cnt <- 0
#     coords$nn$idx <- t(
#       apply(
#         coords$nn$idx,
#         1,
#         function(z) {
#           cnt <<- cnt + 1
#           z[z != cnt]
#         }
#       )
#     )
#     coords$nn$dist <- t(
#       apply(
#         coords$nn$dist,
#         1,
#         function(z) {
#           z[z != 0]
#         }
#       )
#     )
# 
#     # COMPUTE SIGMA
#     sigma <- apply(
#       coords$nn$dist,
#       1,
#       function(dd2) {
#         # t-SNE PERPLEXITY ERROR
#         sigma_error <- function(ss, dd2) {
#           p <- exp(-dd2/(2*ss^2))
#           if(sum(p) < .Machine$double.eps) {
#             p <- 1
#           }
#           p <- p/sum(p)
#           p <- p[p > 0]
#           2^(-sum(p*log2(p))) - coords$perplexity
#         }
#         dd2_min_idx <- 1
#         dd2_ascen <- sort(dd2)
#         while(dd2_ascen[dd2_min_idx] == 0) {
#           dd2_min_idx <- dd2_min_idx + 1
#         }
#         ss_lower <- dd2_ascen[dd2_min_idx]
#         dd2_max_idx <- 1
#         dd2_descen <- sort(dd2, decreasing = TRUE)
#         while(is.infinite(dd2_descen[dd2_max_idx])) {
#           dd2_max_idx <- dd2_max_idx + 1
#         }
#         ss_upper <- dd2_descen[dd2_max_idx]
#         while(sigma_error(ss_upper, dd2) < 0) {
#           ss_lower <- ss_upper
#           ss_upper <- 2 * ss_upper
#         }
#         while(sigma_error(ss_lower, dd2) > 0) {
#           ss_upper <- ss_lower
#           ss_lower <- ss_lower/2
#         }
#         uniroot(
#           sigma_error,
#           dd2,
#           interval = c(ss_lower, ss_upper),
#           tol = (ss_upper - ss_lower) * .Machine$double.eps ^ 0.25
#         )$root
#       }
#     )
# 
#     # COMPUTE PROBABILITIES - ORIGINAL SPACE
#     cnt <- 0
#     orig_p <- t(
#       apply(
#         coords$nn$dist,
#         1,
#         function(z) {
#           cnt <<- cnt + 1
#           p <- exp(-z/(2*sigma[cnt]^2))
#           p/sum(p)
#         }
#       )
#     )
# 
#     # SYMMETRIZE PROBABILITIES - ORIGINAL SPACE
#     for(i in seq_len(nrow(coords$nn$idx))) {
#       for(j2 in seq_len(length(coords$nn$idx[i, ]))) {
#         j <- coords$nn$idx[i, j2]
#         i2 <- match(i, coords$nn$idx[j, ])
#         if(!is.na(i2)) {
#           if(j > i) {
#             p_sym <- orig_p[i, j2] + orig_p[j, i2] / 2
#             orig_p[i, j2] <- p_sym
#             orig_p[j, i2] <- p_sym
#           }
#         } else {
#           orig_p[i, j2] <- orig_p[i, j2]/2
#         }
#       }
#     }
#     orig_p <- sweep(
#       orig_p,
#       1,
#       rowSums(orig_p),
#       "/"
#     )
# 
#     # COMPUTE DISTANCE IN t-SNE SPACE FOR ORIGINAL CLOSE NEIGHBOURS
#     cnt <- 0
#     tsne_dist <- t(
#       apply(
#         coords$nn$idx,
#         1,
#         function(z) {
#           cnt <<- cnt + 1
#           unlist(
#             lapply(
#               z,
#               function(w) {
#                 sum(
#                   (coords$Y[cnt, ] - coords$Y[w, ]) ^ 2
#                 )
#               }
#             )
#           )
#         }
#       )
#     )
# 
#     # COMPUTE PROBABILITIES ASSOCIATED TO t-SNE REPRESENTATION
#     tsne_p_factor <- nrow(coords$Y) / (2*sum(1/(1+dist(coords$Y)^2)))
#     tsne_p <- t(
#       apply(
#         tsne_dist,
#         1,
#         function(z) {
#           tsne_p_factor / (1 + z)
#         }
#       )
#     )
# 
#     # COMPUTE CROSS ENTROPIES
#     res <- unlist(
#       lapply(
#         1:nrow(orig_p),
#         function(z) {
#           -sum(
#             orig_p[z, ] * log(tsne_p[z, ])
#           )
#         }
#       )
#     )
# 
#   }
# 
#   # RETURN CROSS ENTROPIES PER CELL
#   return(res)
# 
# }
