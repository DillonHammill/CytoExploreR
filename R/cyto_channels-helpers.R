## CYTO_CHANNELS ---------------------------------------------------------------

#' Extract channel names
#'
#' Simply a wrapper around \code{colnames} to extract the channels associated
#' with a \code{cytoframe}, \code{cytoset} \code{GatingHierarchy} or
#' \code{GatingSet}.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select vector of channel names to select.
#' @param exclude vector of channel names to exclude.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param ... additional arguments passed to \code{\link[base:grep]{grepl}} for
#'   character matching. For exact character string matching to override the
#'   default which ignores character case, set \code{fixed} to TRUE.
#'
#' @return vector of channel names.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels}}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # GatingSet
#' cyto_channels(gs)
#'
#' # GatingHierarchy
#' cyto_channels(gs[[1]], select = "Alexa")
#'
#' # cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#' cyto_channels(cs)
#'
#' # cytoframe
#' cyto_channels(cs[[1]], exclude = c("FSC","SSC"))
#'
#' @export
cyto_channels <- function(x, 
                          select = NULL,
                          exclude = NULL,
                          append = FALSE,
                          ...){
  
  # LIST
  if(cyto_class(x, "list", TRUE)) {
    x <- unlist(x)[[1]]
  }
  
  # CHANNELS - FLOWFRAME/FLOWSET
  if(cyto_class(x, c("flowFrame", "flowSet"), TRUE)) {
    channels <- BiocGenerics::colnames(x)
  # CHANNELS - CYTOFRAME/CYTOSET/GATINGSET
  } else if(cyto_class(x, c("cytoframe", "cytoset", "GatingSet"))) {
    channels <- flowWorkspace::colnames(x)
  # CHANNELS - DATA.FRAME/MATRIX
  } else {
    channels <- colnames(x)
  }
  
  # SELECT
  if(!is.null(select)){
    ind <- unique(
      LAPPLY(
        select, 
        function(z){
          which(
            suppressWarnings(
              grepl(
                z, 
                channels, 
                ignore.case = TRUE,
                ...
              )
            )
          )
        }
      )
    )
    channels <- channels[ind]
  }
  
  # EXCLUDE
  if(!is.null(exclude)){
    for(z in exclude){
      channels <- channels[!suppressWarnings(
        grepl(z, 
              channels, 
              ignore.case = TRUE,
              ...)
      )]
    }
  }
  
  # MARKERS
  markers <- cyto_markers(x)
  ind <- match(channels, names(markers))
  ind[!is.na(ind)] <- markers[ind[!is.na(ind)]]
  
  # APPEND
  if(append) {
    # APPEND MARKERS
    channels <- paste0(
      "<",
      ind,
      "> ",
      channels
    )
  # STORE MARKERS IN NAMES
  } else {
    names(channels) <- ind
  }
  
  # RETURN CHANNELS
  return(channels)
  
}

## CYTO_CHANNELS REPLACEMENT METHOD --------------------------------------------

#' Replace channel names
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSey}}.
#' @param value vector of new column names to replace the old ones.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_channels}}
#' @seealso \code{\link{cyto_fluor_channels}}
#'
#' @examples 
#' library(CytoExploreRData)
#' 
#' # Activation GatingSet
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                          package = "CytoExploreRData"))
#' 
#' # GatingSet
#' cyto_channels(gs)
#' 
#' # Update first FSC-A to FSC
#' cyto_channels(gs)[1] <- "FSC"
#' cyto_channels(gs)
#'
#' @export
"cyto_channels<-" <- function(x, value) {
  if(cyto_class(x, c("flowFrame", "flowSet"), TRUE)) {
    BiocGenerics::colnames(x) <- value
  } else if(cyto_class(x, c("cytoframe", "cytoset", "GatingSet"))) {
    flowWorkspace::colnames(x) <- value
  } else {
    colnames(x) <- value
  }
  return(x)
}

## CYTO_MARKERS ----------------------------------------------------------------

#' Extract marker names
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} or a list of these
#'   objects.
#' @param select vector of channels or markers for which the channel/marker
#'   combinations should be returned.
#' @param exclude vector of channels or markers for which the channel/marker
#'   combinations should not be returned.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param ... additional arguments passed to \code{\link[base:grep]{grepl}} for
#'   character matching. For exact character string matching to override the
#'   default which ignores character case, set \code{fixed} to TRUE.
#'
#' @return vector of marker names or NULL if no markers have been assigned.
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
#' # Activation GatingSet
#' gs <- cyto_load(system.file("extdata/Activation-GatingSet", 
#'                 package = "CytoExploreRData"))
#'
#' # GatingHierarchy
#' cyto_markers(gs[[1]])
#'
#' # GatingSet
#' cyto_markers(gs)
#'
#' # cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#' cyto_markers(cs)
#'
#' # cytoframe
#' cyto_markers(cs[[1]])
#'
#' @export
cyto_markers <- function(x,
                         select = NULL,
                         exclude = NULL,
                         append = FALSE,
                         ...) {
  
  # LIST
  if(cyto_class(x, "list", TRUE)) {
    x <- unlist(x)[[1]]
  }
  
  # FLOWFRAME/FLOWSET
  if(cyto_class(x, c("flowFrame", "flowSet"), TRUE)) {
    markers <- flowCore::markernames(x)
  # CYTOFRAME/CYTOSET/GATINGHIERARCHY/GATINGSET
  } else if(cyto_class(x, c("cytoframe", "cytoset", "GatingSet"))) {
    markers <- flowWorkspace::markernames(x)
  # MATRIX/DATA.FRAME
  } else {
    markers <- NULL
  }
  
  # INCONSISTENT MARKERS - ALLOW IF SOME SAMPLES UNANNOTATED
  if(cyto_class(markers, "list")) {
    markers[LAPPLY(markers, "length") == 0] <- NULL
    if(length(markers) == 1){
      markers <- markers[[1]]
    }
  }
  
  # MARKER SELECTION/EXCLUSION
  if(!length(markers) == 0) {
    # SELECT
    if(!is.null(select)){
      ind <- unique(LAPPLY(select, function(z){
        which(
          suppressWarnings(
            grepl(
              z,
              markers,
              ignore.case = TRUE,
              ...
            )
          ) |
          suppressWarnings(
            grepl(
              z, 
              names(markers),
              ignore.case = TRUE,
              ...
            )
          )
        )
      }))
      markers <- markers[ind]
    }
    
    # EXCLUDE
    if(!is.null(exclude)){
      # EXCLUDE
      for(z in exclude) {
        markers <- markers[!(suppressWarnings(grepl(z, 
                                                    markers,
                                                    ignore.case = TRUE,
                                                    ...)) |
                               suppressWarnings(grepl(z,
                                                      names(markers),
                                                      ignore.case = TRUE,
                                                      ...)))]
      }
    }
    # APPEND
    if(append) {
      markers <- paste0(
        "<",
        markers,
        "> ",
        names(markers)
      )
    }
  }
  
  # MARKERS
  return(markers)
  
}

## CYTO_FLUOR_CHANNELS ---------------------------------------------------------

#' Extract Fluorescent Channels
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional arguments passed to \code{\link{cyto_channels}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # GatingSet
#' cyto_fluor_channels(gs)
#'
#' # GatingHierarchy
#' cyto_fluor_channels(gs[[1]])
#'
#' # cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#' cyto_fluor_channels(cs)
#'
#' # cytoframe
#' cyto_fluor_channels(cs[[1]])
#' 
#' @rdname cyto_fluor_channels
#'
#' @export
cyto_fluor_channels <- function(x,
                                ...){
  cyto_channels(
    x, 
    exclude = c(
      "FSC",
      "SSC",
      "Time",
      "Original",
      "Sample",
      "Event",
      "UMAP",
      "t-?SNE",
      "PCA",
      "EmbedSOM",
      "FIt-?SNE"
    ),
  ...)
}

## CYTO_CHANNELS_EXTRACT -------------------------------------------------------

#' Extract channels associated with certain markers
#'
#' \code{cyto_channels_extract} will check whether the supplied channels or
#' marker names are valid for the
#' \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#' \code{\link[flowWorkspace:cytoset]{cytoset}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and return a vector of
#' valid channel names. \code{cyto_channels_extract} is particularly useful for
#' determining which channel(s) are associated with particular marker(s).
#'
#' @param x an object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel and/or marker names (e.g. c("Alexa Fluor
#'   700-A","CD8")).
#' @param skip vector of markers/channels in \code{channels} to bypass when
#'   converting to valid channels, for example \code{"Unstained"} is bypassed
#'   when checking channels in the channel match file.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to FALSE by default. If set to TRUE an additional check will be
#'   performed to ensure that only 1 or 2 \code{channels} are supplied.
#' @param ... additional arguments passed to \code{\link[base:grep]{grepl}} for
#'   character matching. For exact character string matching to override the
#'   default which ignores character case, set \code{fixed} to TRUE.
#'
#' @return  A vector of valid channel names.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Extract channels for CD4 & CD8
#' cyto_channels_extract(gs, c("CD4", "CD8"))
#'
#' @export
cyto_channels_extract <- function(x,
                                  channels, 
                                  skip = NULL,
                                  plot = FALSE,
                                  ...) {
  
  # CHANNELS
  chans <- cyto_channels(x)
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # EXTRACT CHANNELS
  res <- c()
  for(z in seq_along(channels)) {
    # SKIP
    if(any(grepl(channels[z], skip, ignore.case = TRUE, ...))) {
      res <- c(res, channels[z])
    # EXACT MARKER MATCH
    } else if(channels[z] %in% markers) {
      res <- c(res, names(markers)[match(channels[z], markers)])
      # EXACT CHANNEL MATCH  
    } else if(channels[z] %in% chans) {
      res <- c(res, chans[match(channels[z], chans)])
      # PARTIAL OR NO MATCH
    } else {
      # PARTIAL MATCHES
      marker_ind <- suppressWarnings(
        which(
          grepl(channels[z],
                markers,
                ignore.case = TRUE,
                ...)
        )
      )
      channel_ind <- suppressWarnings(
        which(
          grepl(channels[z],
                chans,
                ignore.case = TRUE,
                ...)
        )
      )
      # PARTIAL MARKER MATCH
      if(length(marker_ind) != 0) {
        res <- c(res, names(markers)[marker_ind])
        # PARTIAL CHANNEL MATCH
      } else if(length(channel_ind) != 0) {
        res <- c(res, chans[channel_ind])
      } else {
        stop(
          paste0(channels[z], " is not a valid channel or marker for this ", 
                 cyto_class(x, class = TRUE), "!")
        )
      }
    }
  }
  
  # CHECK
  if (plot == TRUE) {
    # res <- res[1:length(channels)]
    if (!length(res) %in% c(1, 2)) {
      stop("Invalid number of supplied channels.")
    }
  }
  
  # CHANNELS
  return(res)
  
}

## CYTO_MARKERS_EXTRACT --------------------------------------------------------

#' Extract marker names for certain channels
#'
#' \code{cyto_markers_extract} will check whether the supplied channels or
#' marker names are valid for the
#' \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#' \code{\link[flowWorkspace:cytoset]{cytoset}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and return a vector of
#' marker names. The name of the channel will be returned if there is no
#' associated marker found.
#'
#' @param x an object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel and/or marker names (e.g. c("Alexa Fluor
#'   700-A","CD8")).
#' @param skip vector of markers/channels in \code{channels} to bypass when
#'   converting to valid markers, for example \code{"Unstained"} is bypassed
#'   when checking markers in the channel match file.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to FALSE by default. If set to TRUE an additional check will be
#'   performed to ensure that only 1 or 2 \code{channels} are supplied.
#' @param ... additional arguments passed to \code{\link[base:grep]{grepl}} for
#'   character matching. For exact character string matching to override the
#'   default which ignores character case, set \code{fixed} to TRUE.
#'
#' @return  A vector of marker names.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Extract channels for CD4 & CD8
#' cyto_markers_extract(gs, c("Alexa Fluor 700-A", "CD8"))
#'
#' @rdname cyto_markers_extract
#'
#' @export
cyto_markers_extract <- function(x, 
                                 channels,
                                 skip = NULL,
                                 append = FALSE,
                                 plot = FALSE,
                                 ...) {
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # CHANNELS
  chans <- cyto_channels(x)
  
  # EXTRACT MARKERS
  res <- c()
  for(z in seq_along(channels)) {
    # SKIP
    if(any(grepl(channels[z], skip, ignore.case = TRUE, ...))) {
      res <- c(res, channels[z])
      names(res[length(res)]) <- channels[z] # append
    # EXACT MARKER MATCH
    } else if(channels[z] %in% markers) {
      res <- c(res, markers[match(channels[z], markers)])
      # EXACT CHANNEL MATCH  
    } else if(channels[z] %in% names(markers)) {
      res <- c(res, markers[match(channels[z], names(markers))])
      # PARTIAL OR NO MATCH
    } else {
      # PARTIAL MATCHES
      marker_ind <- suppressWarnings(
        which(
          grepl(channels[z],
                markers,
                ignore.case = TRUE,
                ...)
        )
      )
      channel_ind <- suppressWarnings(
        which(
          grepl(channels[z],
                names(markers),
                ignore.case = TRUE,
                ...)
        )
      )
      
      # PARTIAL MARKER MATCH
      if(length(marker_ind) != 0) {
        res <- c(res, markers[marker_ind])
        # PARTIAL CHANNEL MATCH
      } else if(length(channel_ind) != 0) {
        res <- c(res, markers[channel_ind])
      } else {
        # CHANNEL UNASSIGNED MARKER - MATCH CHANNEL
        if(channels[z] %in% chans) {
          res <- c(
            res, 
            structure(
              c(channels[z]),
              names = channels[z]
            )
          )
          # CHANNEL UNASSIGNED MARKER - PARTIAL
        } else if(any(grepl(channels[z], chans, ignore.case = TRUE, ...))) {
          res <- c(
            res,
            structure(
              chans[which(grepl(channels[z],
                                chans,
                                ignore.case = TRUE,
                                ...))],
              names = chans[which(grepl(channels[z],
                                        chans,
                                        ignore.case = TRUE,
                                        ...))])
          )
          # INVALID CHANNEL/MARKER
        } else {
          stop(
            paste0(channels[z], " is not a valid channel or marker for this ", 
                   cyto_class(x, class = TRUE), "!")
          )
        }
      }
    }
  }
  
  # CHECK
  if (plot == TRUE) {
    # res <- res[1:length(channels)]
    if (!length(res) %in% c(1, 2)) {
      stop("Invalid number of supplied channels.")
    }
  }
  
  # APPEND
  if(append) {
    res <- paste0(
      "<",
      res,
      "> ",
      names(res)
    )
  }
  
  # MARKERS
  return(res)
}

## CYTO_CHANNEL_SELECT ---------------------------------------------------------

#' Select Fluorescent Channel for Compensation Controls
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#'   compensation controls.
#'
#' @return vector of channels in order of compensation Control samples.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom DataEditR data_edit
#'
#' @examples
#' if(interactive()) {
#' library(CytoExploreRData)
#'
#' # Compensation Gatingset
#' gs <- load_gs(system.file("extdata/Compensation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Select a channel for each control from dropdown menu
#' cyto_channel_select(gs)
#' }
#' 
#' @export
cyto_channel_select <- function(x){
  
  # MESSAGE
  message("Select a fluorescent channel for each of the samples:")
  
  # CHANNEL OPTIONS
  opts <- c(cyto_fluor_channels(x), "Unstained")
  
  # CHANNEL TEMPLATE
  chans <- data.frame("name" = cyto_names(x),
                      "channel" = NA,
                      stringsAsFactors = FALSE)
  
  # CHANNEL SELECTION
  chans <- data_edit(chans,
                     title = "Channel Selector",
                     logo = CytoExploreR_logo(),
                     col_edit = FALSE,
                     row_edit = FALSE,
                     col_options = list("channel" = unname(opts)),
                     col_names = "channel",
                     col_readonly = "name",
                     hide = TRUE,
                     quiet = TRUE,
                     viewer = "pane")
  
  # MISSING CHANNELS
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
#' @param x object of \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channels, one channel for each cytoframe in
#'   \code{x}. Channels are added a \code{"channel"} column within
#'   \code{cyto_details(x)} either by name or the order in which the channels
#'   are supplied. If no channels are supplied the user will be asked to
#'   interactively enter this information using \code{cyto_channel_select()}.
#' @param save_as name to use for the saved channel match csv file, set to
#'   \code{"date-Channel-Match.csv"} by default. To prevent writing to a new
#'   channel match file, set this argument to \code{NA}.
#'
#' @return update \code{cyto_details} of \code{cytoset} or \code{GatingSet},
#'   write channel matching to csv file and return channel matching as a
#'   data.frame.
#'
#' @examples
#' if(interactive()) {
#' library(CytoExploreRData)
#'
#' # Compensation Gatingset
#' gs <- load_gs(system.file("extdata/Compensation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Select a channel for each control from dropdown menu
#' cyto_channel_match(gs)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_channel_match <- function(x,
                               channels = NULL,
                               save_as = NULL) {
  
  # DEFAULT FILENAME
  if (is.null(save_as)) {
    save_as <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Channel-Match.csv"
    )
  }
  
  # SAMPLE NAMES
  nms <- cyto_names(x)
  
  # CHANNELS - SUPPLIED
  if(!is.null(channels)) {
    # INCORRECT NUMBER OF CHANNELS
    if(length(channels) != length(x)) {
      stop(
        "Incorrect number of channels supplied to cyto_channel_match()!"
      )
    # CORRECT NUMBER OF CHANNELS
    } else {
      channel_match <- data.frame(
        "name" = nms,
        "channel" = channels,
        stringsAsFactors = FALSE
      )
      rownames(channel_match) <- NULL
    }
  # NO CHANNELS SUPPLIED
  } else {
    # INTERACTIVELY EDIT CHANNELS
    chans <- cyto_channel_select(x)
    channel_match <- data.frame("name" = nms,
                                "channel" = chans,
                                stringsAsFactors = FALSE)
    rownames(channel_match) <- NULL
  }
  
  # WRITE CHANNEL_MATCH FILE
  if(!.all_na(save_as)) {
    write_to_csv(
      channel_match, 
      save_as
    )
  }
  
  # Update cyto_details
  cyto_details(x)$channel <- channel_match$channel
  
  # Return edited channel match file
  return(channel_match)
}

## CYTO_CHANNEL_MATCH_CHECK ----------------------------------------------------

#' Check and prepare data stored in channel match files
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent vector of parental populations to extract from each sample in
#'   \code{x}, set to the \code{"root"} node by default.
#' @param channel_match a vector of channels one per sample or the name of a CSV
#'   file which indicates the name of the channel associated with each sample If
#'   \code{channel_match} is not supplied or the channels are not annotated in
#'   \code{cyto_details(x)} the user will be asked to interactively enter this
#'   information usng \code{cyto_channel_match()}.
#' @param save_as passed to \code{cyto_channel_match()} to control whether the
#'   channel annotation should be written to a CSV file.
#'
#' @return data.frame containing the updated experimental details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#' library(CytoExploreRData)
#' 
#' # Compensation GatingSet
#' gs <- GatingSet(Compensation)
#' 
#' # Annotate channel matches in experiment details
#' cyto_channel_match_check(gs)
#' }
#'
#' @export
cyto_channel_match_check <- function(x,
                                     parent = "root",
                                     channel_match = NULL,
                                     save_as = NA){
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # ANNOTATE PARENT(S)
  if(!"parent" %in% colnames(pd)) {
    pd[, "parent"] <- rep(
      parent, 
      length.out = length(x)
    )
  }
  
  # CHANNEL_MATCH FILE SEARCH
  if(is.null(channel_match)) {
    channel_match <- grep(
      "channel-match.csv$",
      list.files(),
      ignore.case = TRUE
    )
    if(length(channel_match) == 0) {
      channel_match <- NULL
    } else {
      channel_match <- list.files()[channel_match]
    }
  }
  
  # ANNOTATE CHANNEL(S)
  if(!"channel" %in% colnames(pd)) {
    # CHANNEL_MATCH SUPPLIED
    if(!is.null(channel_match)) {
      # CHANNEL VECTOR OR FILE
      if(is.null(dim(channel_match))) {
        # CHANNEL_MATCH FILE
        if(all(file_exists(file_ext_append(channel_match, ".csv")))) {
          # FIND VALID CHANNEL_MATCH
          channel_match <- structure(
            lapply(
              seq_along(channel_match),
              function(z){
                # READ CHANNEL_MATCH FILE
                cm <- read_from_csv(
                  channel_match[z]
                )
                # CHECK CHANNEL_MATCH
                if(all(c("name", "channel") %in% colnames(cm))) {
                  if(all(cyto_names(x) %in% cm[, "name"])) {
                    return(cm[
                      cyto_match(
                        x,
                        cm[, "name"]
                      ),
                    ])
                  } else {
                    return(NULL)
                  }
                } else {
                  return(NULL)
                }
              }
            )
          )
          channel_match[LAPPLY(channel_match, "is.null")] <- NULL
          if(length(channel_match) == 0) {
            channel_match <- NULL
          } else {
            channel_match <- channel_match[[1]]
          }
          # UPDATE EXPERIMENT DETAILS
          if(!is.null(channel_match)) {
            # UPDATE CHANNEL(S)
            pd[, "channel"] <- channel_match[, "channel"]
            # UPDATE PARENT(S)
            if("parent" %in% colnames(channel_match)) {
              pd[, "parent"] <- channel_match[, "parent"]
            }
          }
        # CHANNEL VECTOR
        } else {
          # UPDATE CHANNEL(S)
          if(length(channel_match) == length(x)) {
            pd[, "channel"] <- cyto_channels_extract(
              x,
              channels = channel_match,
              skip = "Unstained"
            )
          }
        }
      # ARRAY
      } else {
        # NAME COLUMN REQUIRED
        if("name" %in% colnames(pd)) {
          # UPDATE PARENT(S)
          if("parent" %in% colnames(pd)) {
            pd[, "parent"] <- channel_match[
              cyto_match(
                x,
                channel_match[, "name"]
              ),
              "parent"
            ]
          }
          # UPDATE CHANNEL(S)
          if("channel" %in% colnames(pd)) {
            pd[, "channel"] <- channel_match[
              cyto_match(
                x,
                channel_match[, "name"]
              ),
              "channel"
            ]
          }
        }
      }
    }
  }
  
  # CYTO_CHANNEL_MATCH() - LAST RESORT
  if(!"channel" %in% colnames(pd)) {
    channel_match <- cyto_channel_match(
      x,
      save_as = save_as
    )
    pd[, "channel"] <- channel_match[
      cyto_match(
        x,
        channel_match[, "name"]
      ),
      "channel"
    ]
    pd[, "parent"] <- channel_match[
      cyto_match(
        x,
        channel_match[, "name"]
      ),
      "parent"
    ]
  }
  
  # UPDATE EXPERIMENT DETAILS STORED IN X
  cyto_details(x) <- pd[
    cyto_match(
      x,
      pd[, "name"]
    ),
  ]
  
  # RETURN UPDATED EXPERIMENT DETAILS
  return(pd)
  
}


## CYTO_CHANNELS_RESTRICT ------------------------------------------------------

#' Restrict the channels of a cytometry object
#'
#' \code{cyto_channels_restrict} removes any unused channels (channels lacking
#' marker assignments) from a \code{cytoframe}, \code{cytoset},
#' \code{GatingHierarchy} or \code{GatingSet}. By default,
#' \code{cyto_channels_restrict} will always retain any FSC, SSC or Time
#' channels irrespective of marker assignment. Removal of channels that contain
#' marker assignments or channels that are privileged channels (FSC/SSC/Time)
#' can be forced through use of the \code{exclude} argument.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
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
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#
#' # Channels
#' cyto_channels(gs)
#'
#' # Remove unused channels
#' gs <- cyto_channels_restrict(gs)
#'
#' # Channels removed
#' cyto_channels(gs)
#'
#' @export
cyto_channels_restrict <- function(x, 
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
  if(length(markers) != 0){
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
  
  # RESTRICTED CYTOFRAME/CYTOSET
  if(cyto_class(x, c("flowFrame", "flowSet"))) {
    x <- x[, channels_to_keep[ind]]
    # RESTRICTED GATINGSET
  } else {
    gs_cyto_data(x) <- gs_cyto_data(x)[, channels_to_keep[ind]]
  }
  
  # RESTRICTED DATA
  return(x)
  
}
