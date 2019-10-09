# CYTO_MAP ---------------------------------------------------------------------

# UMAP projections can be applied to new datasets but computation is slow. Keep
# implementation to flowFrame method for now.

#' Create dimension-reduced maps of cytometry data
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#' @param channels vector of channels names indicating the channels that should
#'   be used by the dimension reduction algorithm to compute the 2-dimensional
#'   map, set to all channels by default. Restricting the number of channels can
#'   greatly improved processing speed but may result in poorer resolution.
#' @param sample number of events to map, set to 50 000 events by default.
#' @param method direction reduction method to use to generate the map,
#'   supported options include "PCA", "tSNE" and "UMAP".
#' @param save logical indicating whether the mapped \code{flowFrame} or
#'   \code{flowSet} should be saved as .fcs file(s) in a folder in the  current
#'   working directory.
#' @param save_as name of the folder to save the .fcs files to when save is
#'   TRUE, set to "cyto_map" by default.
#' @param plot logical indicating whether the constructed map should be plotted
#'   using \code{cyto_plot}.
#' @param seed integer to set seed prior to mapping to ensure consistent results
#'   between runs.
#' @param ... additional arguments passed to the called dimension reduction
#'   function. Links to the documentation for these functions can be found
#'   below.
#'
#' @importFrom flowCore exprs
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
                               sample = 50000,
                               method = "tSNE",
                               save = FALSE,
                               save_as = "cyto_map",
                               plot = TRUE,
                               seed,
                               ...){
  
  
  # CHANNELS -------------------------------------------------------------------
  
  # PREPARE CHANNELS
  if(missing(channels)){
    channels <- cyto_channels(x, exclude = "Time")
  }
  
  # CONVERT CHANNELS
  channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)

  # PREPARE DATA ---------------------------------------------------------------
  
    # PREPARE DATA - SAMPLING
    x <- cyto_sample(x, display = sample, seed = 56)
  
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
              channels = channels)
  }

  # SAVE MAPPED FLOWFRAME ------------------------------------------------------
  
  # CREATE NEW FOLDER
  if(save == TRUE & !dir.exists(save_as)){
    dir.create(save_as)
  }
  
  # WRITE FCS FILES
  if(save == TRUE){
    write.FCS(x, paste0(save_as, "/", keyword(x, "$FIL")))
  }
  
  # RETURN MAPPED FLOWFRAME ----------------------------------------------------
  return(x)
  
}
