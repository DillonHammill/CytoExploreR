# CYTO_MAP ---------------------------------------------------------------------

#' Create reduced dimension maps of cytometry data
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
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
                               method = "tSNE",
                               save = FALSE,
                               save_as = "cyto_map",
                               plot = TRUE,
                               seed,
                               ...){
  
  # EXTRACT DATA ---------------------------------------------------------------
  
  # MATRIX
  fr_exprs <- exprs(x)
  
  # SET SEED
  if(!missing(seed)){
    set.seed(seed)
  }
  
  # MAPPING --------------------------------------------------------------------
  
  # PCA
  if(grepl(method, "PCA", ignore.case = TRUE)){
    # MAPPING
    mp <- prcomp(fr_exprs, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$x[, 1:2]
  # tSNE
  }else if(grepl(method, "tSNE", ignore.case = TRUE)){
    # MAPPING 
    mp <- Rtsne(fr_exprs, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$Y
  # UMAP 
  }else if(grepl(method, "UMAP", ignore.case = TRUE)){
    # MAPPING
    mp <- umap(fr_exprs, ...)
    # MAPPING CO-ORDINATES
    coords <- mp$layout
  # UNSUPPORTED METHOD  
  }else{
    stop(paste(method, "is not a supported mapping method."))
  }
  
  # CHANNELS -------------------------------------------------------------------
  
  # APPROPRIATE DIMENSION NAMES
  channels <- paste(method, c(1,2), sep = "-")
  
  # APPEND MAPPING TO FLOWFRAME ------------------------------------------------
  
  # APPROPRIATELY NAME MAPPING PARAMETERS
  colnames(coords) <- channels
  
  # CBIND FLOWFRAME & MAPPING
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


#' @rdname cyto_map
#' @export
cyto_map.flowSet <- function(x,
                             method = "tSNE",
                             save = FALSE,
                             save_as = "cyto_map",
                             plot = TRUE,
                             ...){
  
  # MERGE DATA -----------------------------------------------------------------
  
  # MAPPING ON MERGED DATA
  fr <- as(x, "flowFrame")
  
  # MAPPING --------------------------------------------------------------------
  
  # MAPPED FLOWFRAME
  fr <- cyto_map(fr,
                 method = method,
                 plot = plot,
                 save = FALSE,
                 ...)

  # MATCH MAPPING TO INDIVIDUAL FLOWFRAMES -------------------------------------
  
  
  # SAVE MAPPED FLOWSET --------------------------------------------------------
    
  # CREATE NEW FOLDER
  if(save == TRUE & !dir.exists(save_as)){
    dir.create(save_as)
  }
  
  # WRITE FCS FILES
  if(save == TRUE){
    lapply(x, function(z){
      write.FCS(z, paste0(save_as, "/", keyword(z, "$FIL")))
    })
  }
  
  # RETURN MAPPED FLOWSET ------------------------------------------------------
  
  # MAPPING PERFORMED ON POOLED DATA
  return(x)
  
}