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
#'   components}\item{"t-SNE"}{ uses Rtsne::Rtsne() to compute t-SNE
#'   co-ordinates} \item{"FIt-SNE"}{ uses
#'   \url{https://github.com/KlugerLab/FIt-SNE} to compute FIt-SNE co-ordinates
#'   (some additional configuration required)}\item{"UMAP"}{ uses uwot::umap()
#'   to compute UMAP co-ordinates}\item{"EmbedSOM"}{ uses EmbedSOM::SOM() and
#'   EmbedSOM::EmbedSOM() to compute EmbedSOM co-ordinates}\item{"PHATE"}{ uses
#'   phateR::phate() to compute trajectory embedding co-ordinates}}
#' @param scale optional argument to scale each channel prior to computing
#'   dimension-reduced co-ordinates, options include \code{"range"},
#'   \code{"mean"}, \code{"median"} or \code{"zscore"}. Set to \code{"range"} by
#'   default, scaling can be turned off by setting this argument to NULL.
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
#' @param cross_entropy logical indicating whether cross entropies should also
#'   be computed for comparison within \code{\link{cyto_map_compare}}, set to
#'   \code{TRUE} by default. The cross entropies are computed for every cell
#'   using the method described by Roca et al. 2022 and stored as a separate
#'   parameter in each cytoframe. Currently, cross entropies can only be
#'   calculated for UMAP and t-SNE dimension reduction algorithms.
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
                     seed = NULL,
                     inverse = FALSE,
                     trans = NA,
                     plot = TRUE,
                     cross_entropy = TRUE,
                     ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CLASS
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    x <- as(x, "cytoset")
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
  
  # MAP EACH GROUP
  x_data_map <- structure(
    lapply(
      seq_along(x_data_groups), 
      function(z) {
        # CYTO_DATA TO MAP
        cyto_data <- x_data_groups[[z]]
        # CYTOFRAME CONTAINER
        cf <- cyto_data_extract(
          cyto_data,
          parent = parent,
          format = "cytoframe",
          coerce = TRUE,
          events = events,
          barcode = TRUE,
          overwrite = TRUE,
          seed = seed,
          trans = trans,
          inverse = inverse,
          copy = FALSE
        )[[1]][[1]]
        # RAW DATA
        cf_exprs <- cyto_exprs(
          cf,
          channels = channels,
          drop = FALSE
        )
        # RE-SCALE DATA
        if(!is.null(scale)) {
          message(
            paste0(
              "Performing ",
              scale,
              " normalisation on channels..."
            )
          )
          cf_exprs <- cyto_stat_scale(
            cf_exprs,
            type = scale
          )
        }
        # PERFORM MAPPING
        coords <- .cyto_map(
          cf_exprs,
          type = type,
          seed = seed,
          label = label,
          cross_entropy = cross_entropy,
          ...
        )
        # APPEND NEW PARAMETERS
        cf <- cyto_cbind(
          cf, 
          coords
        )
        # SPLIT - CYTOSET
        cyto_split(
          cf,
          names = names[[z]]
        )
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
            cyto_convert(
              x_data_map[[match(x_names[z], x_data_map_names)]]
            )
          # EMPTY CYTOFRAME
          } else {
            cyto_empty(
              x_names[z],
              x_data_map_chans
            )
          }
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
  
  # UPDATE CYTOSET IN GATINGSET
  if(cyto_class(x, "GatingSet")) {
    
    # CANNOT ADD NEW PARAMETERS TO SUBSET OF GATINGSET
    # RECONSTRUCT GATINGSET FROM SCRATCH SO SPILL, TRANS & GATES ARE ATTACHED
    
    # RECONSTRUCTING GATINGSET
    message(
      "Building a new GatingSet..."
    )
    # REVERSE TRANSFORMATIONS
    if(!.all_na(trans) & !inverse) {
      x_data_map <- cyto_transform(
        x_data_map,
        trans = trans,
        inverse = TRUE,
        plot = FALSE,
        quiet = TRUE
      )
    }
    # REVERSE COMPENSATION
    spill <- cyto_spillover_extract(x)
    if(!is.null(spill)) {
      x_data_map <- cyto_compensate(
        x,
        remove = TRUE
      )
    }
    # BUILD NEW GATINGSET
    gs <- GatingSet(x_data_map)
    # RE-APPLY COMPENSATION
    if(!is.null(spill)) {
      message(
        "Compensating for fluorescent spillover..."
      )
      gs <- cyto_compensate(
        gs,
        spillover = spill
      )
    }
    # RE-APPLY TRANSFORMATIONS
    if(!.all_na(trans)) {
      message(
        "Re-applying data transformations..."
      )
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
      message(
        "Recomputing gates..."
      )
      x <- cyto_gateTemplate_apply(
        gs, 
        x
      )
    } else {
      x <- gs
    }

  } else {
    x <- x_data_map
  }
  
  # PLOT 
  if(plot) {
    # cyto_plot_map()
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
                      cross_entropy = TRUE,
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
        "label",
        "cross_entropy"),
      names(args)
    )
  ]
  
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
  
  # SUPPORTED CROSS ENTROPY METHOD
  method <- NULL
  
  # CHARACTER FUNCTION
  if(is.character(type)) {
    # PCA
    if(grepl("^PCA$", type, ignore.case = TRUE)) {
      # MESSAGE
      message(
        "Using rsvd::rpca() to compute PCA co-ordinates..."
      )
      # RSVD
      cyto_require(
        "rsvd",
        source = "CRAN",
        ref = paste0(
          "Erichson NB, Voronin S, Brunton SL, Kutz JN (2019). Randomized ",
          "Matrix Decompositions Using R. Journal of Statistical Software, ",
          "89(11), 1–48. doi: 10.18637/jss.v089.i11."
        )
      )
      # TYPE -> FUNCTION
      type <- cyto_func_match("rsvd::rpca")
    # FIt-SNE
    } else if(grepl("^FIt-?SNE$", type, ignore.case = TRUE)) {
      # MESSAGE
      message(
        "Using FIt-SNE::fftRtsne() to compute FIt-SNE co-ordinates... "
      )
      # TYPE -> FUNCTION
      type <- CytoExploreR::fftRtsne
    # T-SNE
    } else if(grepl("^t-?SNE$", type, ignore.case = TRUE)) {
      # MESSAGE
      message(
        paste0(
          "Using Rtsne::Rtsne() to compute t-SNE co-ordinates..."
        )
      )
      # Rtsne
      cyto_require(
        "Rtsne",
        source = "CRAN",
        ref = paste0(
          "L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional",
          " Data Using t-SNE. Journal of Machine Learning Research",
          " 9(Nov):2579-2605, 2008."
        )
      )
      # TYPE -> FUNCTION
      type <- cyto_func_match(
        "Rtsne::Rtsne"
      )
      # SUPPORTED CROSS ENTROPY METHOD
      method <- "t-SNE"
    # UMAP 
    } else if(grepl("^UMAP$", type, ignore.case = TRUE)) {
      # MESSAGE
      message(
        "Using uwot::umap() to compute UMAP co-ordinates..."
      )
      # UWOT
      cyto_require(
        "uwot",
        source = "CRAN",
        ref = paste0(
          "McInnes L., Healy J. & Melville J. (2018) UMAP: Uniform Manifold ",
          "Approximation and Projection for Dimension Reduction. ",
          "arXiv:1802.03426v3"
        )
      )
      # TYPE -> FUNCTION
      type <- cyto_func_match(
        "uwot::umap"
      )
      # SUPPORTED CROSS ENTROPY METHOD
      method <- "UMAP"
      # CROSS ENTROPY
      if(cross_entropy) {
        args[["ret_nn"]] <- TRUE
        args[["ret_model"]] <- TRUE
      }
    # EMBEDSOM
    } else if(grepl("^Embed-?SOM$", type , ignore.case = TRUE)) {
      # MESSAGE
      message(
        "Using EmbedSOM::EmbedSOM() to compute EmbedSOM co-ordinates..."
      )
      # EMBEDSOM
      cyto_require(
        "EmbedSOM",
        source = "CRAN",
        ref = paste0(
          "Kratchovil M., Koladiya A. & Vondrasek J. (2019) Generalised ",
          "EmbedSOM on quadtree-structured self-organising maps. F1000 ",
          "Research (8:2120)"
        )
      )
      # DATA
      names(args)[1] <- "data"
      # CREATE SOM - FLOWSOM NOT SUPPLIED
      if(!any(c("fsom", "map") %in% names(args))) {
        # DEFAULT SOM GRID SIZE - X
        if(!"xdim" %in% names(args)) {
          args[["xdim"]] <- 24
        }
        # DEFAULT SOM GRID SIZE - Y
        if(!"ydim" %in% names(args)) {
          args[["ydim"]] <- 24
        }
        # SOM
        args[["map"]] <- cyto_func_execute(
          "EmbedSOM::SOM",
          args
        )
      }
      # EMBEDSOM
      cyto_map_coords <- cyto_func_execute(
        "EmbedSOM::EmbedSOM",
        args
      )
      colnames(cyto_map_coords) <- c("EmbedSOM-1", "EmbedSOM-2") 
    # PHATE
    } else if(grepl("^PHATE$", type, ignore.case = TRUE)) {
      # MESSAGE
      message(
        "Using phateR::phate() to compute PHATE co-ordinates..."
      )
      # PHATER
      cyto_require(
        "phateR",
        source = "CRAN",
        ref = paste0(
          "Moon K. et al. (2019) Visualising structure and transitions in ",
          "high dimensional biological data. Nature Biotechnology 37(1482-1492)"
        )
      )
      # TYPE -> FUNCTION
      type <- cyto_func_match(
        "phateR::phate"
      )
    # OTHER
    } else {
      type <- tryCatch(
        cyto_func_match(
          type
        ),
        error = function(e){
          stop(
            "'type'is not the name of a valid dimension reduction function!"
          )
        }
      )
    }
  }
  
  # FUNCTION
  if(is.function(type)) {
    # MESSAGE
    message(
      paste0(
        "Computing ",
        paste0(label, "()"),
        " co-ordinates..."
      )
    )
    # FORMAL ARGUMENTS
    cyto_map_fun_args <- cyto_func_args(type)
    names(args)[match("x", names(args))] <- cyto_map_fun_args[1]
    # CALL MAPPING FUNCTION
    cyto_map_coords <- cyto_func_call(
      type,
      args
    )
    # COMPUTE CROSS ENTROPIES
    if(cross_entropy) {
      # CROSS ENTROPY NOT SUPPORTED
      if(length(method) == 0) {
        warning(
          paste0(
            "Cross entropies can only be computed for UMAP and t-SNE ",
            "dimensionality reduction algorithms!"
          )
        )
      # UMAP
      } else if(method == "UMAP") {
        cyto_map_coords$embedding <- cbind(
          cyto_map_coords$embedding,
          cross_entropy(
            x,
            type = method,
            coords = cyto_map_coords
          )
        )
      # t-SNE
      } else if(method == "t-SNE") {
        cyto_map_coords$Y <- cbind(
          cyto_map_coords$Y,
          cross_entropy(
            x,
            type = method,
            coords = cyto_map_coords
          )
        )
      }
    }
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
    # DIMENSION REDUCTION + CROSS ENTROPY
    if(!is.null(method) & cross_entropy) {
      colnames(cyto_map_coords) <- c(
        paste0(
          label,
          "-",
          1:(ncol(cyto_map_coords) - 1)
        ),
        paste0(
          label,
          "_cross_entropy"
        )
      )
    # DIMENION REDUCTION ONLY
    } else {
      colnames(cyto_map_coords) <- paste0(
        label,
        "-",
        1:ncol(cyto_map_coords)
      )
    }
  }
  
  # RETURN MAPPED DATA
  return(cyto_map_coords)
  
}

#' Internal function to compute cross entropies
#' @noRd
cross_entropy <- function(x,
                          type = "UMAP",
                          coords = NULL) {
  
  # MESSAGE
  message(
    paste0(
      "Computing ",
      type,
      " cross-entropies..."
    )
  )
  
  # REFERENCE
  message(
    paste0(
      "Roca C. et al. (2021) A Cross Entropy test allows quantitative ",
      "statistical comparison of t-SNE and UMAP representations. ",
      "arXiv:2112.04172"
    )
  )
  
  # SIZE
  n <- nrow(x)
  
  # UMAP -----------------------------------------------------------------------
  if(type == "UMAP") {
    
    # FORMAT NEAREST NEIGHBOURS - ORIGINAL SPACE
    cnt <- 0
    coords$nn[[1]]$idx <- t(
      apply(
        coords$nn[[1]]$idx,
        1,
        function(z) {
          cnt <<- cnt + 1
          z[z != cnt]
        }
      )
    )
    coords$nn[[1]]$dist <- t(
      apply(
        coords$nn[[1]]$dist,
        1,
        function(z) {
          z[z != 0]
        }
      )
    )

    # COMPUTE SIGMA
    sigma <- apply(
      coords$nn[[1]]$dist,
      1,
      function(dd){
        # UMAP SIGMA ERROR
        sigma_error <- function(ss, dd) {
          p <- exp(-pmax(0, dd - min(dd)) / ss)
          sum(p) - log2(length(p))
        }
        dd_ascen <- sort(dd[dd > 0])
        ss_lower <- dd_ascen[1]
        dd_descen <- sort(dd[!is.infinite(dd)], decreasing = TRUE)
        ss_upper <- dd_descen[1]
        while(sigma_error(ss_upper, dd) < 0){
          ss_lower <- ss_upper
          ss_upper <- 2 * ss_upper
        }
        while(sigma_error(ss_lower, dd) > 0) {
          ss_upper <- ss_lower
          ss_lower <- ss_lower/2
        }
        uniroot(
          sigma_error,
          dd,
          interval = c(ss_lower, ss_upper),
          tol = (ss_upper - ss_lower) * .Machine$double.eps ^ 0.25
        )$root
      }
    )
    
    # COMPUTE PROBABILITIES - ORIGINAL SPACE
    cnt <- 0
    orig_p <- t(
      apply(
        coords$nn[[1]]$dist,
        1,
        function(z) {
          cnt <<- cnt + 1
          exp(
            -pmax(
              0,
              z - min(z)
            ) /
              sigma[cnt]
          )
        }
      )
    )
    
    # SYMMETRIZE PROBABILITIES - ORIGINAL SPACE
    for(i in seq_len(nrow(coords$nn[[1]]$idx))) {
      for(j2 in seq_len(length(coords$nn[[1]]$idx[i, ]))) {
        j <- coords$nn[[1]]$idx[i, j2]
        i2 <- match(i, coords$nn[[1]]$idx[j, ])
        if(!is.na(i2)) {
          if(j > i) {
            p_sym <- orig_p[i, j2] + orig_p[j, i2] - 
              orig_p[i, j2] * orig_p[j, i2]
            orig_p[i, j2] <- p_sym
            orig_p[j, i2] <- p_sym
          }
        } else {
          # orig_p[i, j2] <- orig_p[i, j2]/2
        }
      }
    }

    # COMPUTE DISTANCE IN UMAP SPACE FOR ORIGINAL CLOSE NEIGHBOURS
    cnt <- 0
    umap_dist <- t(
      apply(
        coords$nn[[1]]$idx,
        1,
        function(z) {
          cnt <<- cnt + 1
          unlist(
            lapply(
              z, 
              function(w) {
                sum(
                  (coords$embedding[cnt, ] - coords$embedding[w, ]) ^ 2
                )
              }
            )
          )
        }
      )
    )

    # COMPUTE PROBABILITIES ASSOCIATED TO UMAP REPRESENTATION
    umap_p <- t(
      apply(
        umap_dist,
        1,
        function(z) {
          1/ (1 + coords$a * z ^ coords$b)
        }
      )
    )

    # COMPUTE FUZZY CROSS ENTROPIES
    res <- unlist(
      lapply(
        1:nrow(orig_p),
        function(z) {
          -sum(
            orig_p[z, ] * log(umap_p[z, ]) +
              (1 - orig_p[z, ]) * log(1 - umap_p[z, ])
          )
        }
      )
    )
    
  # t-SNE ----------------------------------------------------------------------
  } else if(type == "t-SNE" ) {
    
    # RANN REQUIRED - NN2()
    cyto_require(
      "RANN",
      source = "CRAN"
    )
    
    # COMPUTE NEAREST NEIGHBOURS - ORIGINAL SPACE
    coords$nn <- cyto_func_call(
      "RANN::nn2",
      args = list(
        x,
        k = coords$perplexity * 3 + 1
      )
    )
    names(coords$nn) <- c("idx", "dist")
    
    # FORMAT NEAREST NEIGHBOURS - ORIGINAL SPACE
    cnt <- 0
    coords$nn$idx <- t(
      apply(
        coords$nn$idx,
        1,
        function(z) {
          cnt <<- cnt + 1
          z[z != cnt]
        }
      )
    )
    coords$nn$dist <- t(
      apply(
        coords$nn$dist,
        1,
        function(z) {
          z[z != 0]
        }
      )
    )
    
    # COMPUTE SIGMA
    sigma <- apply(
      coords$nn$dist,
      1,
      function(dd2) {
        # t-SNE PERPLEXITY ERROR
        sigma_error <- function(ss, dd2) {
          p <- exp(-dd2/(2*ss^2))
          if(sum(p) < .Machine$double.eps) {
            p <- 1
          }
          p <- p/sum(p)
          p <- p[p > 0]
          2^(-sum(p*log2(p))) - coords$perplexity
        }
        dd2_min_idx <- 1
        dd2_ascen <- sort(dd2)
        while(dd2_ascen[dd2_min_idx] == 0) {
          dd2_min_idx <- dd2_min_idx + 1
        }
        ss_lower <- dd2_ascen[dd2_min_idx]
        dd2_max_idx <- 1
        dd2_descen <- sort(dd2, decreasing = TRUE)
        while(is.infinite(dd2_descen[dd2_max_idx])) {
          dd2_max_idx <- dd2_max_idx + 1
        }
        ss_upper <- dd2_descen[dd2_max_idx]
        while(sigma_error(ss_upper, dd2) < 0) {
          ss_lower <- ss_upper
          ss_upper <- 2 * ss_upper
        }
        while(sigma_error(ss_lower, dd2) > 0) {
          ss_upper <- ss_lower
          ss_lower <- ss_lower/2
        }
        uniroot(
          sigma_error,
          dd2,
          interval = c(ss_lower, ss_upper),
          tol = (ss_upper - ss_lower) * .Machine$double.eps ^ 0.25
        )$root
      }
    )
    
    # COMPUTE PROBABILITIES - ORIGINAL SPACE
    cnt <- 0
    orig_p <- t(
      apply(
        coords$nn$dist,
        1,
        function(z) {
          cnt <<- cnt + 1
          p <- exp(-z/(2*sigma[cnt]^2))
          p/sum(p)
        }
      )
    )
    
    # SYMMETRIZE PROBABILITIES - ORIGINAL SPACE
    for(i in seq_len(nrow(coords$nn$idx))) {
      for(j2 in seq_len(length(coords$nn$idx[i, ]))) {
        j <- coords$nn$idx[i, j2]
        i2 <- match(i, coords$nn$idx[j, ])
        if(!is.na(i2)) {
          if(j > i) {
            p_sym <- orig_p[i, j2] + orig_p[j, i2] / 2
            orig_p[i, j2] <- p_sym
            orig_p[j, i2] <- p_sym
          }
        } else {
          orig_p[i, j2] <- orig_p[i, j2]/2
        }
      }
    }
    orig_p <- sweep(
      orig_p,
      1,
      rowSums(orig_p),
      "/"
    )
    
    # COMPUTE DISTANCE IN t-SNE SPACE FOR ORIGINAL CLOSE NEIGHBOURS
    cnt <- 0
    tsne_dist <- t(
      apply(
        coords$nn$idx,
        1,
        function(z) {
          cnt <<- cnt + 1
          unlist(
            lapply(
              z, 
              function(w) {
                sum(
                  (coords$Y[cnt, ] - coords$Y[w, ]) ^ 2
                )
              }
            )
          )
        }
      )
    )
    
    # COMPUTE PROBABILITIES ASSOCIATED TO t-SNE REPRESENTATION
    tsne_p_factor <- nrow(coords$Y) / (2*sum(1/(1+dist(coords$Y)^2)))
    tsne_p <- t(
      apply(
        tsne_dist,
        1,
        function(z) {
          tsne_p_factor / (1 + z)
        }
      )
    )
    
    # COMPUTE CROSS ENTROPIES
    res <- unlist(
      lapply(
        1:nrow(orig_p),
        function(z) {
          -sum(
            orig_p[z, ] * log(tsne_p[z, ])
          )
        }
      )
    )
    
  }
    
  # RETURN CROSS ENTROPIES PER CELL
  return(res)
  
}
