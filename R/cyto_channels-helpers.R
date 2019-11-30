## CYTO_CHANNELS ---------------------------------------------------------------

#' Extract channel names
#'
#' Simply a wrapper around \code{colnames} to extract the channels associated
#' with a \code{flowFrame}, \code{flowSet}, \code{GatingHierarchy} or
#' \code{GatingSet}.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select vector of channel names to select.
#' @param exclude vector of channel names to exclude.
#'
#' @return vector of channel names.
#'
#' @importFrom BiocGenerics colnames
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels}}
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#'
#' # flowFrame
#' cyto_channels(fs[[1]])
#'
#' # flowSet
#' cyto_channels(fs)
#'
#' # GatingHierarchy
#' cyto_channels(gs[[1]])
#'
#' # GatingSet - exclude FSC & SSC channels
#' cyto_channels(gs, exclude = c("FSC","SSC"))
#'
#' @export
cyto_channels <- function(x, 
                          select = NULL,
                          exclude = NULL){
  
  # CHANNELS
  channels <- colnames(x)
  
  # SELECT
  if(!is.null(select)){
    lapply(select, function(z){
      channels <<- channels[grepl(z, channels, ignore.case = TRUE)]
    })
  }
  
  # EXCLUDE
  if(!is.null(exclude)){
    lapply(exclude, function(z){
      channels <<- channels[!grepl(z, channels, ignore.case = TRUE)]
    })
  }
  
  # RETURN CHANNELS
  return(channels)
  
}

## CYTO_MARKERS ----------------------------------------------------------------

#' Extract marker names
#' 
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'   
#' @return vector of marker names or NULL if no markers have been assigned.
#' 
#' @importFrom flowCore parameters
#' @importFrom flowWorkspace pData
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_channels}}
#' @seealso \code{\link{cyto_fluor_channels}}
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#'
#' # flowFrame
#' cyto_markers(fs[[1]])
#'
#' # flowSet
#' cyto_markers(fs)
#'
#' # GatingHierarchy
#' cyto_markers(gs[[1]])
#'
#' # GatingSet
#' cyto_markers(gs)
#' 
#' @name cyto_markers
NULL

#' @noRd
#' @export
cyto_markers <- function(x){
  UseMethod("cyto_markers")
}

#' @rdname cyto_markers
#' @export
cyto_markers.GatingSet <- function(x){
  # Extract data
  fs <- cyto_extract(x, "root")
  # flowFrame
  fr <- fs[[1]]
  # flowFrame method call
  cyto_markers(fr)
}

#' @rdname cyto_markers
#' @export
cyto_markers.GatingHierarchy <- function(x){
  # Extract data
  fr <- cyto_extract(x, "root")
  # flowFrame method call
  cyto_markers(fr)
}

#' @rdname cyto_markers
#' @export
cyto_markers.flowSet <- function(x){
  # flowFrame
  fr <- x[[1]]
  # flowFrame method call
  cyto_markers(fr)
}

#' @rdname cyto_markers
#' @export
cyto_markers.flowFrame <- function(x){
  # Extract marker information
  markers <- as.character(pData(parameters(x))$desc)
  # Add channels as names
  names(markers) <- as.character(pData(parameters(x))$name)
  # Remove NA entries
  if(.all_na(markers)){
    return(NULL)
  }else{
    markers <- markers[!is.na(markers)]
    return(markers)
  }
}

## CYTO_FLUOR_CHANNELS ---------------------------------------------------------

#' Extract Fluorescent Channels
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Add samples to GatingSet
#' gs <- GatingSet(fs)
#'
#' # Fluorescent channels of flowFrame
#' cyto_fluor_channels(fs[[1]])
#'
#' # Fluorescent channels for a flowSet
#' cyto_fluor_channels(fs)
#'
#' # Fluorescent channels for GatingHierarchy
#' cyto_fluor_channels(gs[[1]])
#'
#' # Fluorescent channels for GatingSet
#' cyto_fluor_channels(gs)
#' @rdname cyto_fluor_channels
#'
#' @export
cyto_fluor_channels <- function(x){
  cyto_channels(x, exclude = c("FSC",
                               "SSC",
                               "Time",
                               "Original",
                               "Sample ID",
                               "Event ID",
                               "UMAP",
                               "tSNE",
                               "PCA",
                               "EmbedSOM"))
}

## CYTO_CHANNELS_EXTRACT -------------------------------------------------------

#' Extract channels associated with certain markers
#'
#' \code{cyto_channels_extract} will check whether the supplied channels or
#' marker names are valid for the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and return a vector of
#' valid channel names. \code{cyto_channels_extract} is particularly useful for
#' determining which channel(s) are associated with particular marker(s).
#'
#' @param x an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
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
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' 
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Add samples to GatingSet
#' gs <- GatingSet(fs)
#'
#' # Extract channels used for CD4 & CD8
#' cyto_channels_extract(gs, c("CD4", "CD8"))
#'
#' @rdname cyto_channels_extract
#'
#' @export
cyto_channels_extract <- function(x,
                                  channels, 
                                  plot = FALSE) {

  # Incorrect channels length
  if (plot == TRUE) {
    if (!length(channels) %in% c(1, 2)) {
      stop("Invalid number of supplied channels.")
    }
  }
  
  # Extract data
  x <- cyto_extract(x)
  
  # Convert to flowFrame
  if(is(x, "flowSet")){
    x <- x[[1]]
  }
  
  # Available channels
  chans <- cyto_channels(x)
  fr_data <- pData(parameters(x))
  
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
        
        
      } else if (channel %in% fr_data$desc) {
        channels[channels %in% channel] <<- as.character(
          fr_data$name[match(channel, fr_data$desc)]
        )
      } else if (!channel %in% chans & !channel %in% fr_data$desc) {
        stop(paste(channel, "is not a valid channel/marker."))
      }
    })
  }
  
  return(channels)
  
}

## CYTO_MARKERS_EXTRACT --------------------------------------------------------

#' Extract marker names for certain channels
#'
#' \code{cyto_markers_extract} will check whether the supplied channels or
#' marker names are valid for the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and return a vector of
#' marker names. The name of the channel will be returned if there is no
#' associated marker found.
#'
#' @param x an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
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
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#'
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Add samples to GatingSet
#' gs <- GatingSet(fs)
#'
#' # Extract markers used for PE-A and Alexa Fluor 488-A
#' cyto_markers_extract(gs, c("Alexa Fluor 488-A","PE-A"))
#'
#' @rdname cyto_markers_extract
#'
#' @export
cyto_markers_extract <- function(x, 
                                 channels,
                                 plot = FALSE) {
  
   # Incorrect channels length
  if (plot == TRUE) {
    if (!length(channels) %in% c(1, 2)) {
      stop("Invalid number of supplied channels.")
    }
  }

  # Extract data 
  x <- cyto_extract(x)
  
  # Extract flowFrame
  if(is(x, "flowSet")){
    x <- x[[1]]
  }
  
  # Available channels
  chans <- cyto_channels(x)
  fr_data <- pData(parameters(x))
  
  # Channel Indices supplied
  if (is.numeric(channels)) {
    channels <- chans[channels]
  }
      
  # Invalid channels or markers
  if(!all(channels %in% c(fr_data$name, fr_data$desc))) {
    stop("'channels' contains invalid channel or marker names.")
  }
  
  # Check if any channels match colnames of flowFrame
  if(any(channels %in% chans)){
    
    # Find indices for valid channels
    ind <- which(channels %in% chans)
    
    # Check if channel has an associated marker
    mrks <- as.vector(fr_data[,"desc"][match(channels[ind], fr_data$name)])
    
    # Update channels with markers
    channels[ind[!is.na(mrks)]] <- mrks[!is.na(mrks)]
    
  }
  
  return(channels)
}

## CYTO_CHANNEL_SELECT ---------------------------------------------------------

#' Select Fluorescent Channel for Compensation Controls
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#'   compensation controls.
#'
#' @return vector of channels in order of compensation Control samples.
#'
#' @importFrom utils menu
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Compensation
#'
#' # Select a channel for each control from dropdown menu
#' cyto_channel_select(fs)
#' }
#' @export
cyto_channel_select <- function(x){

  # Extract data
  x <- cyto_extract(x)
  
  # Convert to list of flowFrames
  x <- cyto_convert(x, "list of flowFrames")
  
  # Channels
  opts <- c(cyto_fluor_channels(x[[1]]), "Unstained")
  
  # Message
  message("Select a fluorescent channel for the following sample:")
  
  # Run through list of flowFrames
  chans <- LAPPLY(seq_len(length(x)), function(z){
    message(cyto_names(x[[z]]))
    opts[menu(choices = opts, graphics = TRUE)]
  })
  
  return(chans)
}
