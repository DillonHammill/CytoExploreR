## CYTO_MAP --------------------------------------------------------------------

#' Create dimension-reduced maps of cytometry data
#'
#' \code{cyto_map} is a convenient wrapper to produce dimension-reduced maps of
#' cytometry data using PCA, tSNE, FIt-SNE, UMAP and EmbedSOM. These
#' dimensionality reduction functions are called using the default settings, but
#' can be altered by passing relevant arguments through \code{cyto_map}. To see
#' a full list of customisable parameters refer to the documentation for each of
#' these functions by clicking on the links below.
#'
#' If you use \code{cyto_map} to map your cytometry data, be sure to cite the
#' publication that describes the dimensionality reduction algorithm that you
#' have chosen to use. References to these publications can be found in the
#' references section of this document.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}} or
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}.
#' @param parent name of the parent population to extract from
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} objects for mapping,
#'   set to the \code{"root"} node by default.
#' @param select designates which samples should be used for mapping when a
#'   \code{flowSet} or \code{GatingSet} object is supplied. Filtering steps
#'   should be comma separated and wrapped in a list. Refer to
#'   \code{\link{cyto_select}}.
#' @param channels vector of channels names indicating the channels that should
#'   be used by the dimension reduction algorithm to compute the 2-dimensional
#'   map, set to all channels with assigned markers by default. Restricting the
#'   number of channels can greatly improve processing speed and resolution.
#' @param display total number of events to map, all events in the combined data
#'   are mapped by default.
#' @param merge_by vector of experimental variables to split and merge samples
#'   into groups prior to mapping, set to "all" by default to create a single
#'   consensus map.
#' @param type dimension reduction type to use to generate the map, supported
#'   options include "PCA", "tSNE", "FIt-SNE", "UMAP" and "EmbedSOM". Users can
#'   also supply the name of a function to perform custom mappings.
#' @param split logical indicating whether samples merged using
#'   \code{cyto_merge_by} should be split prior to writing fcs files, set to
#'   FALSE by default.
#' @param names original names of the samples prior to merging using
#'   \code{cyto_merge_by}, only required when split is TRUE. These names will be
#'   re-assigned to each of split flowFrames and included in the file names.
#' @param save_as passed to \code{cyto_save} to indicate a folder where the
#'   mapped FCS files should be saved, set to NULL by default to turn off saving
#'   of FCS files.
#' @param inverse logical indicating whether the data should be inverse
#'   transformed prior to writing FCS files, set to FALSE by default. Inverse
#'   transformations of \code{flowFrame} or \code{flowSet} objects requires
#'   passing of transformers through the \code{trans} argument.
#' @param trans object of class \code{transformerList} containing the
#'   transformation definitions applied to the supplied data. Used internally
#'   when \code{inverse_transform} is TRUE, to inverse the transformations prior
#'   to writing FCS files.
#' @param plot logical indicating whether the constructed map should be plotted
#'   using \code{cyto_plot}.
#' @param seed integer to set seed prior to mapping to ensure more consistent
#'   results between runs.
#' @param ... additional arguments passed to the called dimension reduction
#'   function. Links to the documentation for these functions can be found
#'   below.
#'
#' @return cytoframe, cytoset, GatingHierarchy or GatingSet containing the
#'   mapped projection parameters.
#'
#' @importFrom flowCore exprs keyword write.FCS flowSet
#' @importFrom flowWorkspace GatingSet gs_cyto_data<- flowSet_to_cytoset
#'   recompute flowFrame_to_cytoframe cytoset gs_cyto_data
#' @importFrom rsvd rpca
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom EmbedSOM SOM EmbedSOM
#'
#' @seealso \code{\link[rsvd:rpca]{PCA}}
#' @seealso \code{\link[Rtsne:Rtsne]{tSNE}}
#' @seealso \code{\link{fftRtsne}}
#' @seealso \code{\link[umap:umap]{UMAP}}
#' @seealso \code{\link[EmbedSOM:SOM]{SOM}}
#' @seealso \code{\link[EmbedSOM:EmbedSOM]{EmbedSOM}}
#'
#' @references N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'   Randomized Matrix Decompositions Using {R}. Journal of Statistical
#'   Software, 89(11), 1-48. \url{http://doi.org/10.18637/jss.v089.i11}.
#' @references N. Halko, P. Martinsson, and J. Tropp. "Finding structure with
#'   randomness: probabilistic algorithms for constructing approximate matrix
#'   decompositions" (2009). (available at arXiv
#'   \url{http://arxiv.org/abs/0909.4061}).
#' @references Gabriel K. (1971). The biplot graphical display of matrices with
#'   application to principal component analysis. Biometrika 58, 453–467.
#'   \url{doi:10.1093/biomet/58.3.453}.
#' @references Maaten, L. van der, & Hinton, G. (2008). Visualizing Data using
#'   t-SNE. Journal of Machine Learning Research 9, 2579–2605.
#'   \url{http://www.jmlr.org/papers/volume9/vandermaaten08a/}.
#' @references Linderman, G., Rachh, M., Hoskins, J., Steinerberger, S.,
#'   Kluger., Y. (2019). Fast interpolation-based t-SNE for improved
#'   visualization of single-cell RNA-seq data. Nature Methods.
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6402590/}.
#' @references McInnes, L., & Healy, J. (2018). UMAP: uniform manifold
#'   approximation and projection for dimension reduction. Preprint at
#'   \url{https://arxiv.org/abs/1802.03426}.
#' @references Kratochvíl, M., Koladiya, A., Balounova, J., Novosadova, V.,
#'   Fišer, K., Sedlacek, R., Vondrášek, J., and Drbal, K. (2018). Rapid
#'   single-cell cytometry data visualization with EmbedSOM. Preprint at
#'   \url{https://www.biorxiv.org/content/10.1101/496869v1}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_map
NULL

#' @rdname cyto_map
#' @export
cyto_map <- function(x, ...) {
  UseMethod("cyto_map")
}

#' @rdname cyto_map
#' @export
cyto_map.GatingSet <- function(x,
                               parent = "root",
                               select = NULL,
                               channels = NULL,
                               display = 1,
                               type = "UMAP",
                               merge_by = "all",
                               split = TRUE,
                               names = NULL,
                               save_as = NULL,
                               inverse = FALSE,
                               trans = NULL,
                               plot = TRUE,
                               seed = NULL,
                               ...) {

  # SELECT DATA (VIEW)
  if (!is.null(select)) {
    x <- cyto_select(x, select)
  }

  # CLONE GATINGSET VIEW
  gs_copy <- cyto_copy(x)

  # TRANSFORMERS
  if (is.null(trans)) {
    trans <- cyto_transformers_extract(gs_copy)
  }

  # NAMES - ALL NAMES INCLUDING EMPTY PARENT
  gs_names <- cyto_names(gs_copy)

  # REMOVE EMPTY PARENTS
  gs_counts <- cyto_stats_compute(gs_copy,
    alias = parent,
    stat = "count"
  )
  gs_clone <- gs_copy[which(gs_counts[, ncol(gs_counts)] != 0)]

  # GROUP_BY
  gs_list <- cyto_group_by(gs_clone, group_by = merge_by)

  # NAMES
  if (is.null(names)) {
    names <- lapply(gs_list, "cyto_names")
  } else {
    names <- split(names, rep(seq_along(gs_list), LAPPLY(gs_list, "length")))
  }

  # MAPPING FUNCTION
  if(is.function(type)){
    cyto_map_fun <- as.character(substitute(type))
    cyto_map_fun <- cyto_map_fun[length(cyto_map_fun)]
  }else{
    cyto_map_fun <- type
  }
  
  # LOOP THROUGH GATINGSETS - RETURN LIST OF FLOWSETS
  cyto_data <- lapply(seq_along(gs_list), function(z) {
    # GATINGSET
    gs <- gs_list[[z]]
    # EXTRACT FLOWSET
    fs <- cyto_data_extract(gs, 
                            parent = parent)[[parent]]
    # MERGE TO FLOWFRAME
    fr <- cyto_merge_by(fs, merge_by = "all")[[1]]
    # MAPPING - RETURNS MERGED FLOWFRAME & SAVES FILES
    cyto_data <- cyto_map(fr,
      channels = channels,
      display = display,
      type = type,
      split = split,
      names = names[[z]],
      save_as = save_as,
      inverse = inverse,
      trans = trans,
      plot = FALSE,
      seed = seed, 
      cyto_map_fun = cyto_map_fun, ...
    )
    # SPLIT - LIST OF FLOWFRAMES
    cyto_data <- cyto_split(cyto_data, names = names[[z]])
    # CONVERT FLOWFRAME LIST TO FLOWSET
    return(flowSet_to_cytoset(flowSet(cyto_data)))
  })

  # COMBINE FLOWSETS
  cyto_data <- do.call("rbind2", cyto_data)

  # MAPPED CHANNELS
  cyto_data_channels <- cyto_channels(cyto_data)

  # MAPPED FILE NAMES
  cyto_data_names <- cyto_names(cyto_data)

  # COMPLETE CYTO_DATA
  gs_cyto_data <- as.list(cyto_names(gs_copy))
  names(gs_cyto_data) <- cyto_names(gs_copy)
  lapply(seq_along(gs_cyto_data), function(z) {
    if (names(gs_cyto_data)[z] %in% cyto_data_names) {
      ind <- match(names(gs_cyto_data)[z], cyto_data_names)
      gs_cyto_data[[z]] <<- cyto_data[[ind]]
    } else {
      # EMPTY FILE
      gs_cyto_data[[z]] <<- cyto_empty(
        names(gs_cyto_data)[z],
        cyto_data_channels
      )
    }
  })
  lapply(seq_along(gs_cyto_data), function(z) {
    cf <- gs_cyto_data[[z]]
    if (class(cf) == "flowFrame") {
      cf <- flowFrame_to_cytoframe(gs_cyto_data[[z]])
    }
    cyto_names(cf) <- names(gs_cyto_data)[z]
    gs_cyto_data[[z]] <<- cf
  })
  gs_cyto_data <- cytoset(gs_cyto_data)
  cyto_details(gs_cyto_data) <- cyto_details(gs_copy)
  gs_cyto_data(gs_copy) <- gs_cyto_data

  # RECOMPUTE STATISTICS
  suppressMessages(recompute(gs_copy))

  # UPDATE GROUPING
  gs_list <- cyto_group_by(gs_copy,
    group_by = merge_by
  )
  
  # PLOT MAPPING PER GROUP (ONE PLOT PER PAGE)
  if (plot == TRUE) {
    lapply(seq_along(gs_list), function(z) {
      # GATINGSET
      gs <- gs_list[[z]]
      # OVERLAY
      overlay <- tryCatch(gh_pop_get_descendants(gs[[1]],
        parent,
        path = "auto"
      ),
      error = function(e) {
        NA
      }
      )
      # LEGEND
      if (!.all_na(overlay)) {
        legend <- TRUE
      } else {
        legend <- FALSE
      }
      # TITLE
      if (names(gs_list)[z] == "all") {
        title <- paste0("Combined Events", "\n", cyto_map_fun)
      } else {
        title <- paste0(names(gs_list)[z], "\n", cyto_map_fun)
      }
      # POINT_COL - FADE BASE LAYER (OVERLAY)
      if (!.all_na(overlay)) {
        point_col <- "grey"
      } else {
        point_col <- NA
      }
      # CYTO_PLOT DESCENDANTS
      tryCatch(cyto_plot(gs,
        parent = parent,
        channels = cyto_channels(gs, select = cyto_map_fun),
        overlay = overlay,
        group_by = "all",
        display = display,
        title = title,
        legend = legend,
        point_col = point_col
      ),
      error = function(e) {
        if (e$message == "figure margins too large") {
          message("Insufficient plotting space, data mapped successfully.")
        }
      }
      )
    })
  }
  
  # RETURN SPLIT MAPPED FLOWFRAMES
  return(gs_copy)
}

#' @rdname cyto_map
#' @export
cyto_map.flowSet <- function(x,
                             select = NULL,
                             channels = NULL,
                             display = 1,
                             type = "UMAP",
                             merge_by = "all",
                             split = TRUE,
                             names = NULL,
                             save_as = NULL,
                             inverse = FALSE,
                             trans = NULL,
                             plot = TRUE,
                             seed = NULL,
                             ...) {

  # SELECT SAMPLES
  if (!is.null(select)) {
    x <- cyto_select(x, select)
  }

  # COPY
  fs_copy <- cyto_copy(x)

  # NAMES - ALL NAMES NCLUDING EMPTY
  fs_names <- cyto_names(fs_copy)

  # REMOVE EMPTY FLOWFRAMES
  fs_counts <- cyto_stats_compute(fs_copy,
    stat = "count"
  )
  fs_clone <- fs_copy[which(fs_counts[, ncol(fs_counts)] != 0)]

  # GROUP_BY
  fs_list <- cyto_group_by(x, group_by = merge_by)

  # NAMES
  if (is.null(names)) {
    names <- lapply(fs_list, "cyto_names")
  } else {
    names <- split(names, rep(seq_along(fs_list), LAPPLY(fs_list, "length")))
  }

  # MAPPING FUNCTION
  if(is.function(type)){
    cyto_map_fun <- as.character(substitute(type))
    cyto_map_fun <- cyto_map_fun[length(cyto_map_fun)]
  }else{
    cyto_map_fun <- type
  }
  
  # LOOP THROUGH FLOWSETS - RETURN LIST OF FLOWSETS
  cyto_data <- lapply(seq_along(fs_list), function(z) {
    # FLOWSET
    fs <- fs_list[[z]]
    # MERGE TO FLOWFRAME
    fr <- cyto_merge_by(fs, merge_by = "all")[[1]]
    # MAPPING - RETURNS MERGED FLOWFRAME & SAVES FILES
    cyto_data <- cyto_map(fr,
      channels = channels,
      display = display,
      type = type,
      split = split,
      names = names[[z]],
      save_as = save_as,
      inverse = inverse,
      trans = trans,
      plot = plot,
      seed = seed,
      cyto_map_fun = cyto_map_fun,
      ...
    )
    # SPLIT - LIST OF FLOWFRAMES
    cyto_data <- cyto_split(cyto_data, names = names[[z]])
    # CONVERT FLOWFRAME LIST TO FLOWSET
    return(flowSet_to_cytoset(flowSet(cyto_data)))
  })

  # COMBINE FLOWSETS
  cyto_data <- do.call("rbind2", cyto_data)

  # MAPPED FILE NAMES
  cyto_data_names <- cyto_names(cyto_data)

  # MAPPED CHANNELS
  cyto_data_channels <- cyto_channels(cyto_data)

  # COMPLETE CYTO_DATA
  fs_cyto_data <- as.list(cyto_names(fs_copy))
  names(fs_cyto_data) <- cyto_names(fs_copy)
  lapply(seq_along(fs_cyto_data), function(z) {
    if (names(fs_cyto_data)[z] %in% cyto_data_names) {
      ind <- match(names(fs_cyto_data)[z], cyto_data_names)
      fs_cyto_data[[z]] <<- cyto_data[[ind]]
    } else {
      # EMPTY FILE
      fs_cyto_data[[z]] <<- cyto_empty(
        names(fs_cyto_data)[z],
        cyto_data_channels
      )
    }
  })
  lapply(seq_along(fs_cyto_data), function(z) {
    cf <- fs_cyto_data[[z]]
    if (class(cf) == "flowFrame") {
      cf <- flowFrame_to_cytoframe(fs_cyto_data[[z]])
    }
    cyto_names(cf) <- names(fs_cyto_data)[z]
    fs_cyto_data[[z]] <<- cf
  })
  fs_cyto_data <- cytoset(fs_cyto_data)
  cyto_details(fs_cyto_data) <- cyto_details(fs_copy)
  cyto_names(fs_cyto_data) <- cyto_names(fs_copy)
  
  # RETURN MAPPED DATA
  return(fs_cyto_data)
}

#' @rdname cyto_map
#' @export
cyto_map.flowFrame <- function(x,
                               channels = NULL,
                               display = 1,
                               type = "UMAP",
                               split = TRUE,
                               names = NULL,
                               save_as = NULL,
                               inverse = FALSE,
                               trans = NULL,
                               plot = TRUE,
                               seed = NULL,
                               ...) {
  
  # CHANNELS -------------------------------------------------------------------

  # PREPARE CHANNELS
  if (is.null(channels)) {
    channels <- cyto_channels(x,
      exclude = c(
        "Time",
        "Original",
        "Sample ID",
        "Event ID",
        "PCA",
        "tSNE",
        "FIt-SNE",
        "UMAP",
        "EmbedSOM"
      )
    )
    channels <- channels[channels %in% names(cyto_markers(x))]
  }

  # CONVERT CHANNELS
  channels <- cyto_channels_extract(x,
    channels = channels,
    plot = FALSE
  )

  # PREPARE DATA ---------------------------------------------------------------

  # PREPARE DATA - SAMPLING
  x <- cyto_sample(x,
    display = display,
    seed = 56
  )

  # EXTRACT RAW DATA MATRIX
  fr_exprs <- cyto_data_extract(x, 
                                format = "matrix")[[1]]

  # RESTRICT MATRIX BY CHANNELS
  fr_exprs <- fr_exprs[, channels]

  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)
  
  # MAPPING FUNCTION
  if("cyto_map_fun" %in% names(args)){
    cyto_map_fun <- args[["cyto_map_fun"]]
  }else{
    if(is.function(type)){
      cyto_map_fun <- as.character(substitute(type))
      cyto_map_fun <- cyto_map_fun[length(cyto_map_fun)]
    }else{
      cyto_map_fun <- type
    }
  }
  
  # ADD CYTO_MAP_FUN TO ARGS
  args[["cyto_map_fun"]] <- cyto_map_fun
  
  # REMOVE EXCESS ARGUMENTS
  cyto_map_args <- formalArgs(cyto_map.flowFrame)
  cyto_map_args <- cyto_map_args[!cyto_map_args %in% c("type",
                                                       "seed",
                                                       "cyto_map_fun")]
  args <- args[!names(args) %in% cyto_map_args]
  
  # RENAME FR_EXPRS TO X
  names(args)[match("fr_exprs", names(args))] <- "x"
  
  # MAPPING --------------------------------------------------------------------

  # MAPPPING COORDS
  coords <- do.call(".cyto_map", args)

  # ADD MAPPING COORDS TO FLOWFRAME
  x <- cyto_cbind(x, coords)

  # VISUALISATION --------------------------------------------------------------
  
  # CYTO_PLOT - MAP
  if (plot == TRUE) {
    tryCatch(cyto_plot(x,
      channels = colnames(coords),
      title = paste0("Combined Events", "\n", cyto_map_fun)
    ),
    error = function(e) {
      if (e$message == "figure margins too large") {
        message("Insufficient plotting space, data mapped successfully.")
      }
    }
    )
  }

  # SAVE MAPPED SAMPLES --------------------------------------------------------

  # CYTO_SAVE - INVERSE TRANSFORMS ONLY APPLIED FOR SAVING
  if (!is.null(save_as)) {
    cyto_save(x,
      split = FALSE,
      names = names,
      save_as = save_as,
      inverse = inverse,
      trans = trans
    )
  }

  # RESET CYTO_PLOT_MAP
  options("cyto_plot_map" = NULL)
  
  # RETURN MAPPED FLOWFRAME ----------------------------------------------------
  return(x)
}

## INTERNAL MAPPING FUNCTION ---------------------------------------------------

#' Obtain dimension-reduced co-ordinates
#' @param x matrix containing the data to be mapped.
#' @noRd
.cyto_map <- function(x,
                      type = "UMAP",
                      seed = NULL,
                      cyto_map_fun = NULL,
                      ...) {

  # SET SEED - RETURN SAME MAP WITH EACH RUN
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ARGUMENTS
  args <- .args_list(...)
  args <- args[-match(c("type", "seed", "cyto_map_fun"), names(args))]
  
  # MAPPING FUNCTION
  if(is.null(cyto_map_fun)){
    if(is.function(type)){
      cyto_map_fun <- as.character(substitute(type))
      cyto_map_fun <- cyto_map_fun[length(cyto_map_fun)]
    }else{
      cyto_map_fun <- type
    }
  }
  
  # CHARACTER
  if (is.character(type)) {
    # MESSAGE
    message(paste0("Computing ", type, " co-ordinates..."))
    # PCA
    if (grepl(type, "PCA", ignore.case = TRUE)) {
      # MAPPING
      mp <- rpca(x, ...)
      # MAPPING CO-ORDINATES
      coords <- mp$x[, 1:2, drop = FALSE]
      colnames(coords) <- c("PCA-1", "PCA-2")
      # tSNE
    } else if (grepl(type, "tSNE", ignore.case = TRUE)) {
      # MAPPING
      mp <- Rtsne(x, ...)
      # MAPPING CO-ORDINATES
      coords <- mp$Y
      colnames(coords) <- c("tSNE-1", "tSNE-2")
      # FIt-SNE
    } else if (grepl(type, "FIt-SNE", ignore.case = TRUE) |
      grepl(type, "FItSNE", ignore.case = TRUE)) {
      mp <- fftRtsne(x, ...)
      # MAPPING CO-ORDINATES
      coords <- mp[, 1:2, drop = FALSE]
      colnames(coords) <- c("FIt-SNE-1", "FIt-SNE-2")
      # UMAP
    } else if (grepl(type, "UMAP", ignore.case = TRUE)) {
      # MAPPING
      mp <- umap(x, ...)
      # MAPPING CO-ORDINATES
      coords <- mp$layout
      colnames(coords) <- c("UMAP-1", "UMAP-2")
      # EmbedSOM
    } else if (grepl(type, "EmbedSOM", ignore.case = TRUE)) {
      # DATA
      names(args)[1] <- "data"
      # CREATE SOM - FLOWSOM NOT SUPPLIED (fsom)
      if (!"fsom" %in% names(args)) {
        # SOM
        mp <- do.call(
          "SOM",
          args[names(args) %in% formalArgs(EmbedSOM::SOM)]
        )
        # SOM ARGUMENTS
        args[["map"]] <- mp
      }
      # EMBEDSOM
      mp <- do.call(
        "EmbedSOM",
        args[names(args) %in% formalArgs(EmbedSOM::EmbedSOM)]
      )
      # MAPPING CO-ORDINATES
      coords <- mp
      colnames(coords) <- c("EmbedSOM-1", "EmbedSOM-2")
      # UNSUPPORTED TYPE
    } else {
      stop(paste(type, "is not a supported mapping type."))
    }
    # CUSTOM FUNCTION
  } else if(is.function(type)) {
    # MESSAGE
    message(paste0("Computing ", 
                   cyto_map_fun, 
                   " co-ordinates..."))
    # RENAME FIRST ARGUMENT
    names(args)[1] <- formalArgs(type)[1]
    # MAPPING
    mp <- do.call(type, args)
    # MAPPING COORDS
    if(is(mp, "matrix")){
      coords <- mp[, c(1,2)]
    }else{
      # LOOOK FOR COORDS
      coords <- NULL
      lapply(seq_len(mp), function(z){
        if(class(mp[z]) == "matrix"){
          if(nrow(mp[z]) == nrow(x)){
            coords <<- cbind(coords, mp[z])
          }
        }
      })
      coords <- coords[, c(1,2)]
    }
    # COLNAMES
    colnames(coords) <- c(paste0(cyto_map_fun, "-1"), 
                          paste0(cyto_map_fun, "-2"))
  }
  
  # RETURN MAPPED COORDS
  return(coords)
}
