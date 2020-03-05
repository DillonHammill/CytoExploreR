## CYTO_MAP --------------------------------------------------------------------

#' Create dimension-reduced maps of cytometry data
#'
#' \code{cyto_map} is a convenient wrapper to produce dimension-reduced maps of
#' cytometry data using PCA, tSNE, FIt-SNE, UMAP and EmbedSOM. These
#' dimensionality reduction functions are called using the default settings, but
#' can be altered by passing relvant arguments through \code{cyto_map}. To see a
#' full list of customisable parameters refer to the documentation for each of
#' these functions by clicking on the links below.
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#' @param parent name of the parent population to extract from
#'   \code{GatingHierarchy} or \code{GatingSet} objects for mapping, set to the
#'   \code{"root"} node by default.
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
#' @param type dimension reduction type to use to generate the map, supported
#'   options include "PCA", "tSNE", "FIt-SNE", "UMAP" and "EmbedSOM".
#' @param merge logical indicating whether the samples should be merged prior to
#'   mapping to generate a consensus map, set to TRUE by default.
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
#' @return flowFrame, flowSet, GatingHierarchy or GatingSet containing the
#'   mapped projection parameters.
#'
#' @importFrom flowCore exprs keyword write.FCS flowSet fr_append_cols
#' @importFrom flowWorkspace GatingSet gs_cyto_data<- flowSet_to_cytoset
#'   recompute
#' @importFrom stats prcomp
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom EmbedSOM SOM EmbedSOM
#'
#' @seealso \code{\link[stats:prcomp]{PCA}}
#' @seealso \code{\link[Rtsne:Rtsne]{tSNE}}
#' @seealso \code{\link{fftRtsne}}
#' @seealso \code{\link[umap:umap]{UMAP}}
#' @seealso \code{\link[EmbedSOM:SOM]{SOM}}
#' @seealso \code{\link[EmbedSOM:EmbedSOM]{EmbedSOM}}
#'
#' @references Gabriel K. (1971). The biplot graphical display of matrices with
#'   application to principal component analysis. Biometrika 58, 453–467.
#'   \url{doi:10.1093/biomet/58.3.453}.
#' @references Maaten, L. van der, & Hinton, G. (2008). Visualizing Data using
#'   t-SNE. Journal of Machine Learning Research 9, 2579–2605.
#'   \url{http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf}.
#' @references Linderman, G., Rachh, M., Hoskins, J., Steinerberger, S.,
#'   Kluger., Y. (2019). Fast interpolation-based t-SNE for improved
#'   visualization of single-cell RNA-seq data. Nature Methods.
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
cyto_map <- function(x, ...){
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
                               merge = TRUE,
                               split = TRUE,
                               names = NULL,
                               save_as = NULL,
                               inverse = FALSE,
                               trans = NULL,
                               plot = TRUE,
                               seed = NULL,
                               ...){
  
  # SELECT DATA (VIEW)
  if(!is.null(select)){
    x <- cyto_select(x, select)
  }  
  
  # CLONE GATINGSET VIEW
  gs_clone <- cyto_copy(x)
  
  # EXTRACT DATA
  fs <- cyto_extract(gs_clone, 
                     parent = parent)
  
  # NAMES
  if(is.null(names)){
    names <- cyto_names(fs)
  }  
  
  # MERGE
  if(merge == TRUE){
    fr_list <- cyto_merge_by(fs, 
                             merge_by = "all")
  }else{
    split <- FALSE
    fr_list <- lapply(cyto_names(fs), function(z){
      fs[[z]]
    })
    names(fr_list) <- cyto_names(fs)
  }
  
  # FLOWFRAME METHOD - FLOWFRAME/FLOWSET RETURN (DEPENDS ON MERGE/SPLIT)
  map_data <- lapply(fr_list, function(fr){
    cyto_map(fr,
             channels = channels,
             display = display,
             type = type,
             split = split,
             names = names,
             save_as = save_as,
             inverse = inverse,
             trans = trans,
             plot = FALSE,
             seed = seed, ...)
  })
  
  # FLOWFRAME LIST TO FLOWSET
  if(length(map_data) > 1){
    map_data <- flowSet_to_cytoset(flowSet(map_data))
  }else{
    map_data <- map_data[[1]]
  }
  
  # SPLIT MERGED FLOWFRAME
  if(is(map_data, "flowFrame")){
    map_data <- cyto_split(map_data,
                           names = names)
    map_data <- flowSet_to_cytoset(flowSet(map_data))
  }
  
  # UPDATE GATINGSET - CLONE
  gs_cyto_data(gs_clone) <- map_data
  suppressMessages(recompute(gs_clone))

  # PLOT MAP
  if(plot == TRUE){
    # OVERLAY
    overlay <- tryCatch(gh_pop_get_descendants(gs_clone[[1]], 
                                               parent,
                                               path = "auto"), 
                        error = function(e){NA})
    # LEGEND
    if(!.all_na(overlay)){
      legend <- TRUE
    }else{
      legend <- FALSE
    }
    # CYTO_PLOT DESCENDANTS
    tryCatch(cyto_plot(gs_clone,
                       parent = parent,
                       channels = cyto_channels(gs_clone, select = type),
                       overlay = overlay,
                       group_by = "all",
                       display = display,
                       title = paste0("Combined Events", "\n", type),
                       legend = legend,
                       point_col = "grey"), # fade base layer
             error = function(e){
               message("Insufficient plotting space, data mapped successfully.")
             })
  }
  
  # RETURN SPLIT MAPPED FLOWFRAMES
  return(gs_clone)
  
}

#' @rdname cyto_map
#' @export
cyto_map.flowSet <- function(x,
                             select = NULL,
                             channels = NULL,
                             display = 1,
                             type = "UMAP",
                             merge = TRUE,
                             split = TRUE,
                             names = NULL,
                             save_as = NULL,
                             inverse = FALSE,
                             trans = NULL,
                             plot = TRUE,
                             seed = NULL,
                             ...){
  
  # COPY
  x <- cyto_copy(x)
  
  # SELECT SAMPLES
  if(!is.null(select)){
    x <- cyto_select(x, select)
  }
  
  # MERGE MULTIPLE SAMPLES
  if(merge == TRUE){
    fr_list <- cyto_merge_by(x, merge_by = "all")
  }else{
    # CAN'T SPLIT UNMERGED SAMPLES
    split <- FALSE
    fr_list <- lapply(cyto_names(x), function(z){
      x[[z]]
    })
    names(fr_list) <- cyto_names(x)
  }
  
  # NAMES
  if(is.null(names)){
    names <- cyto_names(x)
  }
  
  # FLOWFRAME METHOD - EACH LIST ELEMENT
  x <- lapply(fr_list, function(fr){
    cyto_map(fr,
             channels = channels,
             display = display,
             type = type,
             split = split,
             names = names,
             save_as = save_as,
             inverse = inverse,
             trans = trans,
             plot = plot,
             seed = seed, ...)
  })
  
  # FLOWFRAME LIST TO FLOWSET
  if(length(x) > 1){
    x <- flowSet_to_cytoset(flowSet(x))
  }else{
    x <- x[[1]]
  }
  
  # RETURN MAPPED DATA - FLOWFRAME/FLOWSET (DEPENDS ON MERGE & SPLIT)
  return(x)
  
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
                               ...){
  
  
  # CHANNELS -------------------------------------------------------------------
  
  # PREPARE CHANNELS
  if(is.null(channels)){
    channels <- cyto_channels(x, 
                              exclude = c("Time",
                                          "Original",
                                          "Sample ID",
                                          "Event ID",
                                          "PCA",
                                          "tSNE",
                                          "FIt-SNE",
                                          "UMAP",
                                          "EmbedSOM"))
    channels <- channels[channels %in% names(cyto_markers(x))]
  }
  
  # CONVERT CHANNELS
  channels <- cyto_channels_extract(x, 
                                    channels = channels, 
                                    plot = FALSE)

  # PREPARE DATA ---------------------------------------------------------------
  
  # PREPARE DATA - SAMPLING
  x <- cyto_sample(x, 
                   display = display, 
                   seed = 56)
  
  # EXTRACT RAW DATA MATRIX
  fr_exprs <- cyto_extract(x, raw = TRUE)[[1]]
    
  # RESTRICT MATRIX BY CHANNELS
  fr_exprs <- fr_exprs[, channels]
  
  # MAPPING --------------------------------------------------------------------
  
  # MAPPPING COORDS
  coords <- .cyto_map(fr_exprs,
                      type = type,
                      seed = seed, 
                      ...)
  
  # ADD MAPPING COORDS TO FLOWFRAME
  x <- fr_append_cols(x, coords)
  
  # VISUALISATION --------------------------------------------------------------
  
  # CYTO_PLOT - MAP
  if(plot == TRUE){
    tryCatch(cyto_plot(x,
              channels = colnames(coords),
              title = paste0("Combined Events", "\n", type)),
             error = function(e){
               message("Insufficient plotting space. Data has been mapped.")
             })
  }
  
  # SPLIT MAPPED FLOWFRAME -----------------------------------------------------
  
  # SPLIT -> FLOWSET (ELSE MERGED FLOWFRAME)
  if(split == TRUE & "Sample ID" %in% cyto_channels(x)){
    x <- cyto_split(x, names = names)
    x <- flowSet_to_cytoset(flowSet(x))
  }
  
  # SAVE MAPPED SAMPLES --------------------------------------------------------
  
  # CYTO_SAVE - INVERSE TRANSFORMS ONLY APPLIED FOR SAVING
  if(!is.null(save_as)){
    cyto_save(x,
              split = FALSE,
              names = names,
              save_as = save_as,
              inverse = inverse,
              trans = trans)
  }
  
  # RETURN MAPPED FLOWFRAME/FLOWSET --------------------------------------------
  return(x)
  
}

## INTERNAL MAPPING FUNCTION ---------------------------------------------------

#' Obtain dimension-reduced co-ordinates
#' @param x matrix containing the data to be mapped.
#' @noRd 
.cyto_map <- function(x,
                      type = "UMAP",
                      seed = NULL, 
                      ...){
  
  # MESSAGE
  message(paste0("Computing ", type, " co-ordinates..."))
  
  # SET SEED - RETURN SAME MAP WITH EACH RUN
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # PCA
  if(grepl(type, "PCA", ignore.case = TRUE)){
    # MAPPING
    mp <- prcomp(x, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$x[, 1:2, drop = FALSE]
    colnames(coords) <- c("PCA-1","PCA-2")
    # tSNE
  }else if(grepl(type, "tSNE", ignore.case = TRUE)){
    # MAPPING 
    mp <- Rtsne(x, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$Y
    colnames(coords) <- c("tSNE-1","tSNE-2")
    # FIt-SNE
  }else if(grepl(type, "FIt-SNE", ignore.case = TRUE) |
           grepl(type, "FItSNE", ignore.case = TRUE)){  
    mp <- fftRtsne(x, ...)
    # MAPPING CO-ORDINATES
    coords <- mp
    colnames(coords) <- c("FIt-SNE-1", "FIt-SNE-2")  
    # UMAP 
  }else if(grepl(type, "UMAP", ignore.case = TRUE)){
    # MAPPING
    mp <- umap(x, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$layout
    colnames(coords) <- c("UMAP-1","UMAP-2")
    # EmbedSOM
  } else if(grepl(type, "EmbedSOM", ignore.case = TRUE)){
    # DATA
    data <- x
    # PULL DOWN ARGUMENTS
    args <- .args_list(...)
    # CREATE SOM - FLOWSOM NOT SUPPLIED (fsom)
    if(!"fsom" %in% names(args)){
      # SOM
      mp <- do.call("SOM", 
                    args[names(args) %in% formalArgs(EmbedSOM::SOM)])
      # SOM ARGUMENTS
      args[["map"]] <- mp
    }
    # EMBEDSOM
    mp <- do.call("EmbedSOM",
                  args[names(args) %in% formalArgs(EmbedSOM::EmbedSOM)])
    # MAPPING CO-ORDINATES
    coords <- mp
    colnames(coords) <- c("EmbedSOM-1", "EmbedSOM-2")
    # UNSUPPORTED TYPE  
  }else{
    stop(paste(type, "is not a supported mapping type."))
  }
  
  # RETURN MAPPED COORDS
  return(coords)
  
}
