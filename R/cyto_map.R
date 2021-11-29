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
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
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
#'   \itemize{\item{"PCA"}{uses rsvd::rpca() to compute principal
#'   components}\item{"t-SNE"}{uses Rtsne::Rtsne() to compute t-SNE
#'   co-ordinates} \item{"FIt-SNE"}{uses
#'   \url{https://github.com/KlugerLab/FIt-SNE} to compute FIt-SNE co-ordinates
#'   (some additional configuration required)}\item{"UMAP"}{uses uwot::umap() to
#'   compute UMAP co-ordinates}\item{}{EmbedSOM}{uses EmbedSOM::SOM() and
#'   EmbedSOM::EmbedSOM() to compute EmbedSOM co-ordinates}}
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
#' @param ... additional arguments passed to the specified dimension reduction
#'   algorithm
#'
#' @return either a \code{cytoset} or \code{GatingSet} with the
#'   dimension-reduced parameters added as additional channels in each of the
#'   underlyig \code{cytoframes}.
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
                     ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
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
    lapply(seq_along(x_data_groups), function(z){
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
    }),
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
    # EXTRACT CO-ORDINATES (ALL DIMENSIONS)
    if(!cyto_class(cyto_map_coords, "matrix")) {
      # FIND CO-ORDINATES
      coords <- NULL
      lapply(seq_along(cyto_map_coords), function(z){
        if(cyto_class(cyto_map_coords[[z]], "matrix")) {
          if(nrow(cyto_map_coords[[z]]) == nrow(x)) {
            coords <<- cbind(coords, cyto_map_coords[[z]])
          }
        }
      })
      cyto_map_coords <- coords
    }
    # COLUMN NAMES
    colnames(cyto_map_coords) <- paste0(
      label,
      "-",
      1:ncol(cyto_map_coords)
    )
  }
  
  # RETURN MAPPED DATA
  return(cyto_map_coords)
  
}
