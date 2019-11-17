## CYTO_MAP --------------------------------------------------------------------

#' Create dimension-reduced maps of cytometry data
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#' @param channels vector of channels names indicating the channels that should
#'   be used by the dimension reduction algorithm to compute the 2-dimensional
#'   map, set to all channels by default. Restricting the number of channels can
#'   greatly improved processing speed but may result in poorer resolution.
#' @param display number of events to map, set to 50 000 events by default.
#' @param type direction reduction method to use to generate the map, supported
#'   options include "PCA", "tSNE" and "UMAP".
#' @param save logical indicating whether the mapped \code{flowFrame} or
#'   \code{flowSet} should be saved as .fcs file(s) in a folder in the current
#'   working directory.
#' @param split logical indicating whether samples merged using
#'   \code{cyto_merge_by} should be split prior to writing fcs files, set to
#'   FALSE by default.
#' @param names original names of the samples prior to merging using
#'   \code{cyto_merge_by}, only required when split is TRUE. These names will be
#'   re-assigned to each of split flowFrames and included in the file names.
#' @param save_as name of the folder to save the .fcs files to when save is
#'   TRUE, set to "cyto_map" by default.
#' @param trans object of class \code{transformerList} containg the
#'   transformation definitions applied to the supplied data. If transformations
#'   are supplied, the data will be inverse transformed prior to saving to
#'   return the data on the original linear scale.
#' @param plot logical indicating whether the constructed map should be plotted
#'   using \code{cyto_plot}.
#' @param seed integer to set seed prior to mapping to ensure consistent results
#'   between runs.
#' @param ... additional arguments passed to the called dimension reduction
#'   function. Links to the documentation for these functions can be found
#'   below.
#'
#' @return flowFrame or list of split flowFrames containing the mapped
#'   projection parameters.
#'
#' @importFrom flowCore exprs keyword
#' @importFrom stats prcomp
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom flowCore write.FCS
#'
#' @seealso \code{\link[stats:prcomp]{PCA}}
#' @seealso \code{\link[Rtsne:Rtsne]{tSNE}}
#' @seealso \code{\link[umap:umap]{UMAP}}
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
cyto_map.flowFrame <- function(x,
                               channels,
                               display = 50000,
                               method = "UMAP",
                               save = TRUE,
                               split = TRUE,
                               names = NULL,
                               save_as = "cyto_map",
                               trans = NULL,
                               plot = TRUE,
                               seed,
                               ...){
  
  
  # CHANNELS -------------------------------------------------------------------
  
  # PREPARE CHANNELS
  if(missing(channels)){
    channels <- cyto_channels(x, 
                              exclude = c("Time",
                                          "Original",
                                          "Sample ID",
                                          "Event ID"))
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
  fr_exprs <- exprs(x)
    
  # RESTRICT MATRIX BY CHANNELS
  fr_exprs <- fr_exprs[, channels]
  
  # MAPPING --------------------------------------------------------------------
  
  # SET SEED - RETURN SAME MAP WITH EACH RUN
  if(!missing(seed)){
    set.seed(seed)
  }
  
  # PCA
  if(grepl(method, "PCA", ignore.case = TRUE)){
    # MAPPING
    mp <- prcomp(fr_exprs, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$x[, 1:2, drop = FALSE]
    colnames(coords) <- c("PCA-1","PCA-2")
    # ADD MAPPING COORDS TO FLOWFRAME
    x <- cbind(x, coords)
  # tSNE
  }else if(grepl(method, "tSNE", ignore.case = TRUE)){
    # MAPPING 
    mp <- Rtsne(fr_exprs, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$Y
    colnames(coords) <- c("tSNE-1","tSNE-2")

  # UMAP 
  }else if(grepl(method, "UMAP", ignore.case = TRUE)){
    # MAPPING
    mp <- umap(fr_exprs, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$layout
    colnames(coords) <- c("UMAP-1","UMAP-2")
  # UNSUPPORTED METHOD  
  }else{
    stop(paste(method, "is not a supported mapping method."))
  }
  
  # ADD MAPPING COORDS TO FLOWFRAME
  x <- cbind(x, coords)
    
  # VISUALISATION --------------------------------------------------------------
  
  # CYTO_PLOT - MAP
  if(plot == TRUE){
    cyto_plot(x,
              channels = colnames(coords))
  }

  # SAVE MAPPED FLOWFRAME(S) ---------------------------------------------------
  
  # CYTO_SAVE
  if(save == TRUE){
    fr_list <- cyto_save(x,
                       split = split,
                       names = names,
                       save_as = save_as,
                       trans = trans)
  }else{
    fr_list <- x
  }
  
  # RETURN MAPPED FLOWFRAME ----------------------------------------------------
  return(fr_list)
  
}
