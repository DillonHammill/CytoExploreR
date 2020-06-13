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
#' @param ... additional arguments passed to \code{\link[base:grepl]{grepl}} for
#'   character matching. For exact character string matching to override the
#'   default which ignores character case, set \code{fixed} to TRUE.
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
                          exclude = NULL,
                          ...){
  
  # CHANNELS
  channels <- colnames(x)
  
  # SELECT
  if(!is.null(select)){
    ind <- unique(LAPPLY(select, function(z){
      which(
        suppressWarnings(
          grepl(z, 
              channels, 
              ignore.case = TRUE,
              ...))
        )
    }))
    channels <- channels[ind]
  }
  
  # EXCLUDE
  if(!is.null(exclude)){
    lapply(exclude, function(z){
      channels <<- channels[!suppressWarnings(
        grepl(z, 
              channels, 
              ignore.case = TRUE,
              ...)
        )]
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
#' @param select vector of channels or markers for which the channel/marker
#'   combinations should be returned.
#' @param exclude vector of channels or markers for which the channel/marker
#'   combinations should not be returned.
#' @param ... additional arguments passed to \code{\link[base:grepl]{grepl}} for
#'   character matching. For exact character string matching to override the
#'   default which ignores character case, set \code{fixed} to TRUE.
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
cyto_markers <- function(x, ...){
  UseMethod("cyto_markers")
}

#' @rdname cyto_markers
#' @export
cyto_markers.GatingSet <- function(x,
                                   select = NULL,
                                   exclude = NULL,
                                   ...){
  # Extract data
  fs <- cyto_extract(x, "root")
  # flowFrame
  fr <- fs[[1]]
  # flowFrame method call
  cyto_markers(fr,
               select = select,
               exclude = exclude,
               ...)
}

#' @rdname cyto_markers
#' @export
cyto_markers.GatingHierarchy <- function(x,
                                         select = NULL,
                                         exclude = NULL,
                                         ...){
  # Extract data
  fr <- cyto_extract(x, "root")
  # flowFrame method call
  cyto_markers(fr,
               select = select,
               exclude = exclude,
               ...)
}

#' @rdname cyto_markers
#' @export
cyto_markers.flowSet <- function(x,
                                 select = NULL,
                                 exclude = NULL,
                                 ...){
  # flowFrame
  fr <- x[[1]]
  # flowFrame method call
  cyto_markers(fr,
               select = select,
               exclude = exclude,
               ...)
}

#' @rdname cyto_markers
#' @export
cyto_markers.flowFrame <- function(x,
                                   select = NULL,
                                   exclude = NULL,
                                   ...){
  # Extract marker information
  markers <- as.character(pData(parameters(x))$desc)
  # Add channels as names
  names(markers) <- as.character(pData(parameters(x))$name)
  # Remove NA entries
  if(.all_na(markers)){
    return(NULL)
  }else{
    # MARKERS
    markers <- markers[!is.na(markers)]
    # SELECT
    if(!is.null(select)){
      ind <- unique(LAPPLY(select, function(z){
        which(suppressWarnings(grepl(z,
                    markers,
                    ignore.case = TRUE,
                    ...)) |
              suppressWarnings(grepl(z, 
                    names(markers),
                    ignore.case = TRUE,
                    ...)))
      }))
      markers <- markers[ind]
    }
    # EXCLUDE
    if(!is.null(exclude)){
      # EXCLUDE
      lapply(exclude, function(z){
        markers <<- markers[!(suppressWarnings(grepl(z, 
                                   markers,
                                   ignore.case = TRUE,
                                   ...)) |
                              suppressWarnings(grepl(z,
                                    names(markers),
                                    ignore.case = TRUE,
                                    ...)))]
      })
    }
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
  
  # Extract channels
  LAPPLY(channels, function(z){
    if(.all_na(z)){
      
    }else if(z %in% chans){
      
    } else if (z %in% fr_data$desc) {
      channels[channels %in% z] <<- as.character(
        fr_data$name[match(z, fr_data$desc)]
      )
    } else if (!z %in% chans & !z %in% fr_data$desc) {
      stop(paste(z, "is not a valid channel/marker."))
    }
  })
  
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
  
  # Message
  message("Select a fluorescent channel for each of the samples:")
  
  # Channels
  opts <- c(cyto_fluor_channels(x), "Unstained")
  
  # Samples
  nms <- cyto_names(x)
  
  # Channel selection
  chans <- data_editor(data.frame("name" = nms,
                                  "channel" = NA,
                                  stringsAsFactors = FALSE),
                       title = "Channel Selector",
                       type = "selector",
                       options = opts)
  
  # Missing channels
  lapply(seq_along(chans[, "channel"]), function(z){
    if(is.na(chans[z, "channel"])){
      stop(paste0("No channel selected for ", chans[z, "name"], "."))
    }
  })
  
  # RETURN VECTOR OF CHANNELS
  chans <- chans[, "channel"]
  return(chans)
}

## CYTO_CHANNEL_MATCH ----------------------------------------------------------

#' Table Editor for Channel Match File Construction
#'
#' @param x object of \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param save_as name to use for the saved channel match csv file, set to
#'   \code{"date-Channel-Match.csv"}.
#' @param menu logical indicating whether channels should be selected from a
#'   drop down menu instead of manually typing them in.
#'   
#' @return update \code{cyto_details} of \code{flowSet} or \code{GatingSet},
#'   write channel matching to csv file and return channel matching as a
#'   data.frame.
#'
#' @importFrom utils write.csv edit
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#' # Load in CytoExploreRData to access data
#' library(CytoExploreR)
#'
#' # Generate channel match file for compensation controls
#' cyto_channel_match(Compensation)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_channel_match <- function(x,
                               save_as = NULL,
                               menu = TRUE) {
  
  # Set default name for channel_match file
  if (is.null(save_as)) {
    save_as <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Channel-Match.csv"
    )
    # save_as must contain csv file extension
  } else {
    # File extension missing
    if (file_ext(save_as) != "csv") {
      save_as <- paste0(save_as, ".csv")
    }
  }
  
  # Extract sample names
  nms <- cyto_names(x)
  
  # Select channels from menu
  if(menu == TRUE){
    chans <- cyto_channel_select(x)
    channel_match <- data.frame("name" = nms,
                                "channel" = chans)
    colnames(channel_match) <- c("name", "channel")
    rownames(channel_match) <- NULL
    # Manually type in channels
  }else if(menu == FALSE){
    # Construct data.frame for editing
    channel_match <- data.frame("name" = nms, 
                                "channel" = rep("NA", length(nms)))
    colnames(channel_match) <- c("name", "channel")
    rownames(channel_match) <- NULL
    
    # Edit channel_match
    channel_match <- suppressWarnings(edit(channel_match))
  }
  
  # Convert markers to channels
  channel_match$channel <- LAPPLY(channel_match$channel, function(z){
    if(!grepl(z, "unstained", ignore.case = TRUE)){
      cyto_channels_extract(x, z)
    }else{
      z
    }
  })
  
  # Check that all channels are valid or throw an error
  if (!all(channel_match$channel %in% c("Unstained", "unstained",
                                        cyto_fluor_channels(x)))) {
    stop("Some inputs in the channel column are not valid.")
  }
  
  # Write edited channel match file to csv file
  write.csv(channel_match, save_as, row.names = FALSE)
  
  # Update cyto_details
  cyto_details(x)$channel <- channel_match$channel
  
  # Return edited channel match file
  return(channel_match)
}

## CYTO_CHANNELS_RESTRICT ------------------------------------------------------

#' Remove channels of a flowFrame or flowSet
#'
#' \code{cyto_channels_restrict} removes channels from the flowFrame or flowSet
#' which do not have markers assigned. The FSC, SSC and Time parameters are
#' always retained irrespective of marker assignment.
#'
#' @param x object of class flowFrame or flowSet.
#'
#' @return flowFrame or flowSet object with unused channels removed.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#'
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Channels
#' cyto_channels(fs)
#'
#' # Remove unused channels
#' fs <- cyto_channels_restrict(fs)
#'
#' # Channels removed
#' cyto_channels(fs)
#'
#' @export
cyto_channels_restrict <- function(x){
  
  # CHANNELS
  channels <- cyto_channels(x)
  
  # FLUOR CHANNELS
  fluor_channels <- cyto_fluor_channels(x)
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # IGNORE FSC/SSC/Time
  ignore_channels <- channels[!channels %in% fluor_channels]
  
  # CHANNELS USED
  used_channels <- fluor_channels[which(fluor_channels %in% names(markers))]
  
  # CHANNELS TO KEEP
  keep_channels <- c(ignore_channels, used_channels)
  
  # SORT CHANNELS AS BEFORE
  ind <- match(channels, keep_channels)
  ind <- ind[!is.na(ind)]
  
  # RETURN RESTRICTED FLOWFRAME/FLOWSET
  x <- x[, keep_channels[ind]]
  return(x)
  
}
