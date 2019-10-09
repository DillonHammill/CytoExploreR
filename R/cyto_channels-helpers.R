# CYTO_CHANNELS ----------------------------------------------------------------

#' Extract channel names
#'
#' Simply a wrapper around \code{\link[BiocGenerics:row+colnames]{colnames}} to
#' extract the channels associated with a \code{flowFrame}, \code{flowSet},
#' \code{GatingHierarchy} or \code{GatingSet}.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
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
#' # flowset
#' cyto_channels(fs)
#' 
#' # GatingHierachy
#' cyto_channels(gs[[1]])
#' 
#' # GatingSet - exclude FSC & SSC channels
#' cyto_channels(gs, exclude = c("FSC","SSC"))
#'
#' @name cyto_channels
cyto_channels <- function(x, exclude = NULL){
  
  # CHANNELS
  channels <- colnames(x)
  
  # EXCLUDE
  if(!is.null(exclude)){
    lapply(exclude, function(z){
      channels <<- channels[!grepl(z, channels, ignore.case = TRUE)]
    })
  }
  
  # RETURN CHANNELS
  return(channels)
  
}

# CYTO_CHANNELS_EXTRACT --------------------------------------------------------

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

# CYTO_MARKERS_EXTRACT ---------------------------------------------------------

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

# CYTO_FLUOR_CHANNELS ----------------------------------------------------------

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
#' gs <- GtaingSet(fs)
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
cyto_fluor_channels <- function(x) {
  UseMethod("cyto_fluor_channels")
}

#' @rdname cyto_fluor_channels
#' @export
cyto_fluor_channels.flowFrame <- function(x) {
  channels <- unname(BiocGenerics::colnames(x))
  
  # Remove FSC channels
  channels <- channels[!grepl("FSC", channels, ignore.case = TRUE)]
  # Remove SSC channels
  channels <- channels[!grepl("SSC", channels, ignore.case = TRUE)]
  # Remove Original channel
  channels <- channels[!grepl("Original", channels)]
  # Remove Time channel
  channels <- channels[!grepl("Time", channels, ignore.case = TRUE)]
  
  return(channels)
}

#' @rdname cyto_fluor_channels
#' @export
cyto_fluor_channels.flowSet <- function(x) {
  cyto_fluor_channels.flowFrame(x[[1]])
}

#' @rdname cyto_fluor_channels
#' @export
cyto_fluor_channels.GatingHierarchy <- function(x) {
  cyto_fluor_channels(cyto_extract(x, "root"))
}

#' @rdname cyto_fluor_channels
#' @export
cyto_fluor_channels.GatingSet <- function(x) {
  cyto_fluor_channels.flowFrame(cyto_extract(x, "root")[[1]])
}

# CYTO_CHANNEL_SELECT ----------------------------------------------------------

#' Select Fluorescent Channel for Compensation Controls
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#'   compensation controls.
#'
#' @return vector of channels in order of compensation Control samples.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Compensation
#'
#' # Select a channel for each control from dropdown menu
#' cyto_channel_select(fs)
#' @name cyto_channel_select
NULL

#' @noRd
#' @export
cyto_channel_select <- function(x){
  UseMethod("cyto_channel_select")
}

#' @rdname cyto_channel_select
#' @export
cyto_channel_select.flowFrame <- function(x) {
  
  # Channels
  opts <- cyto_fluor_channels(x)
  
  # Print sample name and select channel
  message(
    paste(
      "Select a fluorescent channel for the following sample:",
      cyto_names(x)
    )
  )
  
  chan <- opts[menu(choices = opts, graphics = TRUE)]
  
  return(chan)
}

#' @rdname cyto_channel_select
#' @export
cyto_channel_select.flowSet <- function(x) {
  
  # Assign x to fs
  fs <- x
  
  opts <- c(cyto_fluor_channels(fs), "Unstained")
  
  # Messgae to select sample
  message("Select a fluorescent channel for the following sample:")
  
  # Print sample name and select channel
  chans <- opts[LAPPLY(cyto_names(fs), function(z) {
    message(z)
    
    menu(choices = opts, graphics = TRUE)
  })]
  
  return(chans)
}

#' @rdname cyto_channel_select
#' @export
cyto_channel_select.GatingHierarchy <- function(x) {
  
  # Assign x to gs
  gh <- x
  
  # Extract flowFrame
  fr <- cyto_extract(gh, "root")
  
  opts <- cyto_fluor_channels(fr)
  
  # Print sample name and select channel
  message(
    paste(
      "Select a fluorescent channel for the following sample:",
      cyto_names(fr)
    )
  )
  
  chan <- opts[menu(choices = opts, graphics = TRUE)]
  
  return(chan)
}

#' @rdname cyto_channel_select
#' @export
cyto_channel_select.GatingSet <- function(x) {
  
  # Assign x to gs
  gs <- x
  
  # Extract flowSet
  fs <- cyto_extract(gs, "root")
  
  # Channel options
  opts <- c(cyto_fluor_channels(fs), "Unstained")
  
  # Print sample name and select channel
  chans <- opts[LAPPLY(cyto_names(fs), function(x) {
    message(
      paste(
        "Select a fluorescent channel for the following sample:",
        x
      )
    )
    
    menu(choices = opts, graphics = TRUE)
  })]
  
  return(chans)
}
