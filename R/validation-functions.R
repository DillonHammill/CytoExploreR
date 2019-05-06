## Function Definitions for Exported Validation Functions ----------------------

# Check channels and return channels given marker names ------------------------

#' Extract channels associated with certain markers
#'
#' \code{cyto_channels_extract} will check whether the supplied channels or
#' marker names are valid for the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowset}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and return a vector of
#' valid channel names. \code{cyto_channels_extract} is particularly useful for
#' determining which channel(s) are associated with particular marker(s).
#'
#' @param x an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel and/or marker names (e.g. c("Alexa Fluor
#'   700-A","CD8")).
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to FALSE by default. If set to TRUE an additional check will be
#'   performed to ensure that only 1 or 2 \code{channels} are supplied.
#'
#' @return  A vector of valid channel names.
#'
#' @importFrom flowWorkspace pData getData
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_channels_extract
#'
#' @export
cyto_channels_extract <- function(x, channels, plot = FALSE) {
  UseMethod("cyto_channels_extract")
}

#' @rdname cyto_channels_extract
#' @export
cyto_channels_extract.default <- function(x, channels, plot = FALSE) {
  warning(
    paste("cyto_channels_extract expects objects of class flowFrame,",
          "flowSet, GatingHierarchy or GatingSet.")
  )
}

#' @rdname cyto_channels_extract
#' @export
cyto_channels_extract.flowFrame <- function(x, channels, plot = FALSE) {
  
  # Incorrect channels length
  if (plot == TRUE) {
    if (!length(channels) %in% c(1, 2)) {
      stop("Invalid number of supplied channels.")
    }
  }
  
  # Available channels
  chans <- BiocGenerics::colnames(x)
  fr.data <- pData(parameters(x))
  
  # Channel Indices supplied
  if (is.numeric(channels)) {
    channels <- chans[channels]
  }
  
  # Check if channels match colnames of flowFrame
  if (all(channels %in% chans)) {
    
    # Supplied channels are valid
  } else if (!all(channels %in% chans)) {
    lapply(channels, function(channel) {
      if (channel %in% chans) {
        
        
      } else if (channel %in% fr.data$desc) {
        channels[channels %in% channel] <<- as.character(
          fr.data$name[match(channel, fr.data$desc)]
        )
      } else if (!channel %in% chans & !channel %in% fr.data$desc) {
        stop(paste(channel, "is not a valid channel/marker."))
      }
    })
  }
  
  return(channels)
  
}

#' @rdname cyto_channels_extract
#' @export
cyto_channels_extract.flowSet <- function(x, channels, plot = FALSE) {
  
  cyto_channels_extract.flowFrame(x[[1]], channels, plot)
  
}

#' @rdname cyto_channels_extract
#' @export
cyto_channels_extract.GatingHierarchy <- function(x, channels, plot = FALSE) {
  
  cyto_channels_extract.flowFrame(getData(x, "root"), channels, plot)
  
}

#' @rdname cyto_channels_extract
#' @export
cyto_channels_extract.GatingSet <- function(x, channels, plot = FALSE) {
  
  cyto_channels_extract.flowFrame(getData(x, "root")[[1]], channels, plot)
  
}

# Check markers and return marker names given channels -------------------------

#' Extract marker names for certain channels
#'
#' \code{cyto_markers_extract} will check whether the supplied channels or
#' marker names are valid for the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowset}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and return a vector of
#' marker names. The name of the channel will be returned if there is no
#' associated marker found.
#'
#' @param x an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel and/or marker names (e.g. c("Alexa Fluor
#'   700-A","CD8")).
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to FALSE by default. If set to TRUE an additional check will be
#'   performed to ensure that only 1 or 2 \code{channels} are supplied.
#'
#' @return  A vector of marker names.
#'
#' @importFrom flowWorkspace pData getData
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_markers_extract
#'
#' @export
cyto_markers_extract <- function(x, markers, plot = FALSE) {
  UseMethod("cyto_markers_extract")
}

#' @rdname cyto_markers_extract
#' @export
cyto_markers_extract.flowFrame <- function(x, channels, plot = FALSE){
  
  # Incorrect channels length
  if (plot == TRUE) {
    if (!length(channels) %in% c(1, 2)) {
      stop("Invalid number of supplied channels.")
    }
  }
  
  # Available channels
  chans <- BiocGenerics::colnames(x)
  fr.data <- pData(parameters(x))
  
  # Channel Indices supplied
  if (is.numeric(channels)) {
    channels <- chans[channels]
  }
  
  # Check if any channels match colnames of flowFrame
  if(any(channels %in% chans)){
    
    # Find indices for valid channels
    ind <- which(channels %in% chans)
    
    # Check if channel has an associated marker
    mrks <- as.vector(fr.data[,"desc"][match(channels[ind], fr.data$name)])
    
    # Update channels with markers
    channels[ind[!is.na(mrks)]] <- mrks[!is.na(mrks)]
    
  # channels probably already includes markers
  }else{
    
    if(!all(channels %in% fr.data$name)) {
      warning("'channels' contains invalid channel or marker names.")
    }
    
  }
  
  return(channels)

}

#' @rdname cyto_markers_extract
#' @export
cyto_markers_extract.flowSet <- function(x, channels, plot = FALSE) {
  
  cyto_markers_extract.flowFrame(x[[1]], channels, plot)
  
}

#' @rdname cyto_markers_extract
#' @export
cyto_markers_extract.GatingHierarchy <- function(x, channels, plot = FALSE) {
  
  cyto_markers_extract(getData(x, "root"), channels, plot)
  
}

#' @rdname cyto_markers_extract
#' @export
cyto_markers_extract.GatingSet <- function(x, channels, plot = FALSE) {
  
  cyto_markers_extract(getData(x, "root")[[1]], channels, plot)
  
}

# Return transformList with Transformations or Inverse Transformations ---------

#' Convert transformerLists to transformLists and/or extract inverse
#' transformations
#'
#' @param x object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the samples.
#' @param inverse logical indicating whether the returned transformList should
#'   contain the inverse transformations.
#'
#' @return A \code{\link[flowCore:transformList-class]{transformList}} containing
#'   the desired transformations.
#'
#' @importFrom flowCore inverseLogicleTransform transformList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Add fs to GatingSet
#' gs <- GatingSet(fs)
#'
#' # Get transformList using estimateLogicle
#' trans <- estimateLogicle(fs[[32]], cyto_fluor_channels(fs))
#'
#' # Get transformList containing inverse transformarions
#' inv <- cyto_trans_check(trans, inverse = TRUE)
#'
#' # Convert transformerList into transformList
#' trans <- estimateLogicle(gs[[32]], cyto_fluor_channels(gs))
#'
#' # Convert transformerList into inverse transformList
#' inv <- cyto_trans_check(trans, inverse = FALSE)
#'
#' @rdname cyto_trans_convert
#'
#' @export
cyto_trans_convert <- function(x, inverse = FALSE) {
  UseMethod("cyto_trans_convert")
}

#' @rdname cyto_trans_convert
#' @noRd
cyto_trans_convert.default <- function(x, inverse = FALSE){
  if(is.null(x)){
    return(x)
  }else{
    warning(
      paste("cyto_trans_check expects objects of class",
            "transformList or transformerList.")
    )
  }
}

#' @rdname cyto_trans_convert
#' @export
cyto_trans_convert.transformList <- function(x, inverse = FALSE){
  
  if(inverse){
    x <- inverseLogicleTransform(x)
  }
  
  return(x)
  
}

#' @rdname cyto_trans_convert
#' @export
cyto_trans_convert.transformerList <- function(x, inverse = FALSE){
  
  if(inverse){
    x <- transformList(names(x), lapply(x, `[[`, "inverse"))
  }else{
    x <- transformList(names(x), lapply(x, `[[`, "transform"))
  }
  
  return(x)
  
}
