## CYTO_CHANNELS ---------------------------------------------------------------

#' Extract channel names
#'
#' Simply a wrapper around \code{colnames} to extract the channels associated
#' with a \code{flowFrame}, \code{flowSet} or \code{GatingSet}.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select vector of channel names to select.
#' @param exclude vector of channel names to exclude.
#' @param ... additional arguments passed to \code{\link[base:grep]{grepl}} for
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
#' @param ... additional arguments passed to \code{\link[base:grep]{grepl}} for
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
#' @param ... additional arguments passed to \code{\link{cyto_channels}}.
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
cyto_fluor_channels <- function(x,
                                ...){
  cyto_channels(x, exclude = c("FSC",
                               "SSC",
                               "Time",
                               "Original",
                               "Sample ID",
                               "Event ID",
                               "UMAP",
                               "tSNE",
                               "PCA",
                               "EmbedSOM"),
                ...)
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
#' @importFrom DataEditR data_edit
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
  chans <- data_edit(data.frame("name" = nms,
                                "channel" = NA,
                                stringsAsFactors = FALSE),
                     title = "Channel Selector",
                     logo = CytoExploreR_logo(),
                     col_edit = FALSE,
                     row_edit = FALSE,
                     col_options = list("channel" = opts),
                     col_names = "channel",
                     col_readonly = "name",
                     quiet = TRUE)
  
  # Missing channels
  lapply(seq_along(chans[, "channel"]), function(z){
    if(is.na(chans[z, "channel"])){
      stop(paste0("No channel selected for ", chans[z, "name"], "."))
    }
  })
  
  # RETURN VECTOR OF CHANNELS
  return(chans[, "channel"])
}

## CYTO_CHANNEL_MATCH ----------------------------------------------------------

#' Create a csv file assigning a channel to each compensation control
#'
#' @param x object of \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param save_as name to use for the saved channel match csv file, set to
#'   \code{"date-Channel-Match.csv"}.
#'   
#' @return update \code{cyto_details} of \code{flowSet} or \code{GatingSet},
#'   write channel matching to csv file and return channel matching as a
#'   data.frame.
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
                               save_as = NULL) {
  
  # Set default name for channel_match file
  if (is.null(save_as)) {
    save_as <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Channel-Match.csv"
    )
  }
  
  # Extract sample names
  nms <- cyto_names(x)
  
  # Channels data to edit
  chans <- cyto_channel_select(x)
  channel_match <- data.frame("name" = nms,
                              "channel" = chans)
  colnames(channel_match) <- c("name", "channel")
  rownames(channel_match) <- NULL
  
  # Write edited channel match file to csv file
  write_to_csv(channel_match, save_as)
  
  # Update cyto_details
  cyto_details(x)$channel <- channel_match$channel
  
  # Return edited channel match file
  return(channel_match)
}

## CYTO_CHANNELS_RESTRICT ------------------------------------------------------

#' Restrict the channels of a cytometry object
#'
#' \code{cyto_channels_restrict} removes any unused channels (channels lacking
#' marker assignments) from a \code{flowFrame}, \code{flowSet},
#' \code{GatingHierarchy} or \code{GatingSet}. By default,
#' \code{cyto_channels_restrict} will always retain any FSC, SSC or Time
#' channels irrespective of marker assignment. Removal of channels that contain
#' marker assignments or channels that are privileged channels (FSC/SSC/Time)
#' can be forced through use of the \code{exclude} argument.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet}or \code{GatingSet}.
#' @param exclude vector of privileged channels or markers to remove in addition
#'   to the channels removed by default, set to NULL by default.
#' @param ... additional arguments passed to \code{\link{cyto_channels}} to
#'   control character matching for the exclude argument.
#'
#' @return an object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet} with unused channels removed.
#'
#' @importFrom flowWorkspace gs_cyto_data
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
#' @name cyto_channels_restrict
NULL

#' @noRd
#' @export
cyto_channels_restrict <- function(x, ...){
  UseMethod("cyto_channels_restrict")
}

#' @rdname cyto_channels_restrict
#' @export
cyto_channels_restrict.flowFrame <- function(x, 
                                             exclude = NULL,
                                             ...){
  
  # CHANNELS
  channels <- cyto_channels(x)
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # PERFORM DEFAULT CHANNEL REMOVAL --------------------------------------------
  
  # PRIVELEGED CHANNELS
  channels_exempt <- cyto_channels(x, 
                                   select = c("FSC",
                                              "SSC",
                                              "Time"))
  
  # CHANNELS WITH MARKERS ASSIGNED
  if(!is.null(markers)){
    channels_to_keep <- names(markers)
    names(channels_to_keep) <- markers
    channels_to_keep <- c(channels_exempt, channels_to_keep)
  }else{
    channels_to_keep <- channels_exempt
  }
  
  # REMOVE DUPLICATED CHANNELS
  channels_to_keep <- channels_to_keep[!duplicated(channels_to_keep)]
  
  # REMOVE PRIVILEGED CHANNELS -------------------------------------------------
  
  # EXCLUDE
  if(!is.null(exclude)){
    
    # MARKERS TO REMOVE
    channels_to_remove <- cyto_markers(x, 
                                       select = exclude,
                                       ...)
    if(length(channels_to_remove) > 0){
      channels_to_remove <- names(channels_to_remove)
    }else{
      channels_to_remove <- c()
    }
    
    # CHANNELS TO REMOVE
    channels_to_remove <- c(channels_to_remove,
                            cyto_channels(x,
                                          select = exclude,
                                          ...))
    channels_to_remove <- unique(channels_to_remove)
    
    # UPDATE CHANNELS TO KEEP
    channels_to_keep <- channels_to_keep[!channels_to_keep %in% 
                                           channels_to_remove]
    
  }
  
  # RESTRICT CHANNELS ----------------------------------------------------------
  
  # SORT CHANNELS AS BEFORE
  ind <- match(channels, channels_to_keep)
  ind <- ind[!is.na(ind)]
  
  # RETURN RESTRICTED FLOWFRAME/FLOWSET
  x <- x[, channels_to_keep[ind]]
  return(x)
  
}

#' @rdname cyto_channels_restrict
#' @export
cyto_channels_restrict.flowSet <- function(x,
                                           exclude = NULL,
                                           ...){
  
  # CHANNELS
  channels <- cyto_channels(x)
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # PERFORM DEFAULT CHANNEL REMOVAL --------------------------------------------
  
  # PRIVELEGED CHANNELS
  channels_exempt <- cyto_channels(x, 
                                   select = c("FSC",
                                              "SSC",
                                              "Time"))
  
  # CHANNELS WITH MARKERS ASSIGNED
  if(!is.null(markers)){
    channels_to_keep <- names(markers)
    names(channels_to_keep) <- markers
    channels_to_keep <- c(channels_exempt, channels_to_keep)
  }else{
    channels_to_keep <- channels_exempt
  }
  
  # REMOVE DUPLICATED CHANNELS
  channels_to_keep <- channels_to_keep[!duplicated(channels_to_keep)]
  
  # REMOVE PRIVILEGED CHANNELS -------------------------------------------------
  
  # EXCLUDE
  if(!is.null(exclude)){
    
    # MARKERS TO REMOVE
    channels_to_remove <- cyto_markers(x, 
                                       select = exclude,
                                       ...)
    if(length(channels_to_remove) > 0){
      channels_to_remove <- names(channels_to_remove)
    }else{
      channels_to_remove <- c()
    }
    
    # CHANNELS TO REMOVE
    channels_to_remove <- c(channels_to_remove,
                            cyto_channels(x,
                                          select = exclude,
                                          ...))
    channels_to_remove <- unique(channels_to_remove)
    
    # UPDATE CHANNELS TO KEEP
    channels_to_keep <- channels_to_keep[!channels_to_keep %in% 
                                           channels_to_remove]
    
  }
  
  # RESTRICT CHANNELS ----------------------------------------------------------
  
  # SORT CHANNELS AS BEFORE
  ind <- match(channels, channels_to_keep)
  ind <- ind[!is.na(ind)]
  
  # RETURN RESTRICTED FLOWFRAME/FLOWSET
  x <- x[, channels_to_keep[ind]]
  return(x)
  
}

#' @rdname cyto_channels_restrict
#' @export
cyto_channels_restrict.GatingSet <- function(x,
                                             exclude = NULL,
                                             ...){
  
  # EXTRACT DATA
  cyto_data <- gs_cyto_data(x)
  
  # RESTRICT CHANNELS
  cyto_data <- cyto_channels_restrict(cyto_data,
                                      exclude = exclude,
                                      ...)
  
  # REPLACE DATA
  gs_cyto_data(x) <- cyto_data
  
  # RETURN GATINSET
  return(x)
  
}
