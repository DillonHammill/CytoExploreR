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
  
  # TODO: SHOULD WE EXCLUDE HEIGHT/WIDTH PARAMTERS HERE?
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
      res <- c(
        res, 
        structure(
          names(markers)[match(channels[z], markers)],
          names = markers[match(channels[z], markers)]
        )
      )
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
        res <- c(
          res, 
          structure(
            names(markers)[marker_ind],
            names = markers[marker_ind]
          )
        )
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
  lapply(
    seq_along(chans[, "channel"]),
    function(z){
      if(is.na(chans[z, "channel"])){
        stop(paste0("No channel selected for ", chans[z, "name"], "."))
      }
    }
  )
  
  # RETURN VECTOR OF CHANNELS
  return(chans[, "channel"])
}

## CYTO_CHANNEL_MATCH ----------------------------------------------------------

#' Match each compensation control to a fluorescent channel
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels names of the possible channels or markers to be matched to
#'   each control, set to all area fluorescent parameters by default.
#' @param file name of a CSV file from which the channel matching information
#'   should be inherited. If not supplied, \code{cyto_channel_match()} will
#'   automatically search for a file named \code{"Compensation-Details.csv"} or
#'   create interactively create such a file.
#' @param save_as name of a CSV file to which the channel matching should be
#'   written for downstream use, set to \code{"Compensation-Details.csv"}
#'   prefixed with the date by default. Users can set custom file names here,
#'   but the file name should contain \code{"Compensation-Details.csv"} in order
#'   to be automatically detected by CytoExploreR within
#'   \code{cyto_spillover_compute()}, \code{cyto_spillover_edit()},
#'   \code{cyto_spillover_spread()} and \code{cyto_plot_compensation()}.
#' @param strip logical indicating whether overlapping characters in file names
#'   should be stripped prior to matching, set to TRUE by default for more
#'   accurate matching.
#' @param ignore.case logical indicating whether all case insensitive matches
#'   should be found, set to TRUE by default.
#' @param insertions logical indicating whether character insertions are allowed
#'   when matching the marker/channel combinations to the sample names, set to
#'   FALSE by default to only allow deletions.
#' @param ... additional arguments passed to \code{link{grep}} when performing
#'   character matching.
#'
#' @return a data.frame written to a CSV file containing information matching
#'   each file name to a channel. This channel matching information is also
#'   automatically added to the \code{cyto_details()} of the supplied samples
#'   where it can be easily accessed by CytoExploreR downstream.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Compensation GatingSet
#' gs <- GatingSet(Compensation)
#'
#' # Channel matching
#' cyto_channel_match(gs)
#'
#' @export
cyto_channel_match <- function(x,
                               channels = NULL,
                               file = NULL,
                               save_as = NULL,
                               strip = TRUE,
                               ignore.case = TRUE,
                               insertions = FALSE,
                               ...) {
  
  # CYTOFRAMES NOT SUPPORTED
  if(cyto_class(x, "flowFrame")) {
    stop(
      paste0(
        "cyto_channel_match() only supports objects of class cytoset, ",
        "GatingHierarchy or GatingSet!"
      )
    )
  }
  
  # BYPASS CHANNEL MATCHING - INTERNAL USE ONLY!
  # USED IN CYTO_SPILLOVER_EDIT() | CYTO_PLOT_COMPENSATION()
  args <- list(...)
  if("channel_match" %in% names(args)) {
    cyto_details(x) <- channel_match[
      match(
        rownames(cyto_details(x)),
        rownames(channel_match)
      ), , drop = FALSE
    ]
    return(x)
  }
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    # EXCLUDE HEIGHT/WIDTH PARAMETERS
    channels <- channels[!grepl("-H$|-W$", channels, ignore.case = TRUE)]
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # SEARCH FOR FILE TO IMPORT DETAILS
  if(is.null(file)) {
    # FILE SEARCH
    pd <- cyto_file_search(
      "Details.csv$",
      colnames = c("name", "group", "parent", "channel"),
      rownames = rownames(cyto_details(x))
    )
    # NO DETAILS FOUND
    if(length(pd) == 0) {
      pd <- cyto_details(x)
      # CHANNEL MATCH FILE LOCATED
    } else {
      # MULTIPLE CHANNEL MATCH FILES LOCATED
      if(length(pd) > 1) {
        # ENQUIRE
        if(interactive() & cyto_option("CytoExploreR_interactive")) {
          message(
            paste0(
              "Multiple files found with channel matching details for this ",
              cyto_class(x),
              ". Which file would you like to inherit these details from?"
            )
          )
          message(
            paste0(
              paste0(
                seq_along(pd),
                ": ",
                names(pd)
              ),
              sep = "\n"
            )
          )
          opt <- cyto_enquire(NULL)
          opt <- tryCatch(
            as.numeric(opt),
            warning = function(w) {
              match(opt, names(pd))
            }
          )
          file <- names(pd)[opt]
          pd <- pd[[opt]]
          # NON-INTERACTIVE FILE SELECTION
        } else {
          warning(
            paste0(
              "Multiple files found with channel matching details for this ",
              cyto_class(x),
              " - resorting to using ",
              names(pd)[1],
              "..."
            )
          )
          file <- names(pd)[1]
          pd <- pd[[1]]
        }
        # SINGLE CHANNEL MATCH FILE LOCATED
      } else {
        message(
          paste0(
            "Importing saved channel matching details from ",
            names(pd)[1],
            "..."
          )
        )
        file <- names(pd)[1]
        pd <- pd[[1]]
      }
    }
    # FILE MANUALLY SUPPLIED
  } else {
    # FILE EXTENSION
    file <- file_ext_append(file, ".csv")
    # FILE EXISTS
    if(file_exists(file, error = FALSE)) {
      message(
        paste0(
          "Importing experiment details from ",
          file,
          "..."
        )
      )
      # READ FILE
      pd <- read_from_csv(file)
      # CHECK FILE
      if(!"name" %in% colnames(pd) |
         !all(rownames(cyto_details(x)) %in% rownames(pd))) {
        stop(
          paste0(
            file,
            " must have rownames and contain entries for every sample in ",
            "this ",
            cyto_class(x),
            "!"
          )
        )
      }
      # FILE DOES NOT EXIST
    } else {
      warning(
        paste0(
          file,
          " does not exist! Resorting to cyto_details(x) instead..."
        )
      )
      pd <- cyto_details(x)
    }
  }
  
  # SAVE_AS
  if(is.null(save_as)) {
    # FILE MISSING
    if(is.null(file)) {
      save_as <- paste0(
        format(Sys.Date(), "%d%m%y"),
        "-", "Compensation-Details.csv"
      )
      # SAVE TO FILE
    } else {
      save_as <- file
    }
  }
  
  # ADD MISSING PARAMETERS
  vars <- c("group", "parent", "channel")
  vars <- vars[!vars %in% colnames(pd)]
  
  # PREPARE MISSING VARIABLES
  if(length(vars) > 0) {
    vars <- matrix(
      NA,
      nrow = nrow(pd),
      ncol = length(vars),
      dimnames = list(
        rownames(pd),
        vars
      )
    )
    pd <- cbind(
      pd,
      vars
    )
  }
  
  # NAMES OF SAMPLES - CHANNEL_MATCH MAY CONTAIN EXTRA SAMPLES
  nms <- rownames(cyto_details(x))
  x_ind <- which(rownames(pd) %in% rownames(cyto_details(x)))
  
  # GROUPS/PARENTS -------------------------------------------------------------
  
  # GROUPS/PARENTS IGNORED FOR CYTOSETS
  if(cyto_class(x, "GatingSet")) {
    # PARENTS MISSING
    ind <- which(is.na(pd$parent))
    # RESTRICT TO AVAILABLE SAMPLES ONLY
    ind <- ind[ind %in% x_ind]
    if(length(ind) > 0) {
      # TERMINAL NODES
      pops <- cyto_nodes(
        x[ind],
        terminal = TRUE,
        path = "auto"
      )
      # COMPUTE COUNTS FOR EACH TERMINAL NODE
      pop_stats <- cyto_apply(
        x[ind],
        parent = pops,
        channels = channels[1],
        input = "matrix",
        FUN = "cyto_stat_count",
        copy = FALSE
      )
      if(cyto_class(pop_stats, "list", TRUE)) {
        pop_stats <- do.call("cbind", pop_stats)
        dimnames(pop_stats) <- list(rownames(pop_stats), pops)
      }
      # PARENT - TERMINAL NODE MOST EVENTS
      pd$parent[ind] <- pops[
        apply(
          pop_stats,
          1,
          "which.max"
        )
      ]
    }
    # GROUPS MISSING - DEFAULT TO PARENTS
    ind <- which(is.na(pd$group))
    if(length(ind) > 0) {
      pd$group[ind] <- pd$parent[ind]
    }
  }
  
  # CHANNELS -------------------------------------------------------------------
  
  # MISSING CHANNEL ASSIGNMENTS
  ind <- which(is.na(pd$channel))
  ind <- ind[ind %in% x_ind]
  
  # DEFAULT CHANNEL ASSIGNMENTS
  if(length(ind) > 0) {
    # NAMES
    file_names <- rownames(pd)[ind]
    names(file_names) <- file_names
    # LOCATE UNSTAINED CONTROL(S)
    unst_ind <- grep("Unst|NIL", file_names, ignore.case = TRUE)
    # ANNOTATE UNSTAINED CONTROLS - REMOVE FROM FILE NAMES
    if(length(unst_ind) > 0) {
      pd$channel[
        match(file_names[unst_ind], rownames(pd))
      ] <- "Unstained"
      file_names <- file_names[-unst_ind]
    }
    # ADDITIONAL FILE NAMES TO MATCH
    if(length(file_names) > 0) {
      # REPLACE ABNORMAL CHARACTERS
      file_names <- gsub(",2f,", "/", file_names, ignore.case = TRUE)
      # REMOVE WHITESPACE & SPECIAL CHARACTERS
      file_names <- gsub("[^[:alnum:]]", "", file_names)
      # STRIP OVERLAPPING CHARACTERS - NO PADDING REQUIRED
      if(strip) {
        file_names <- .cyto_string_strip(file_names)
      }
      # COMBINE MARKERS & CHANNELS FOR MATCHING - "MARKER CHANNEL"
      channels <- unlist(
        lapply(
          seq_along(channels),
          function(z) {
            # CHANNEL
            channel <- channels[z]
            # MARKER(S) ASSIGNED
            if(!is.na(names(channel)) & names(channel) != "NA") {
              # HANDLE MULTIPLE MARKERS PER CHANNEL
              markers <- strsplit(names(channel), "\\||\\/")[[1]]
              # STORE CHANNEL NAME
              channel <- rep(channel, length.out = length(markers))
              # APPEND CHANNEL NAME
              channel <- structure(
                paste(
                  markers,
                  channel
                ),
                names = channel
              )
              # NO MARKER
            } else {
              names(channel) <- channel # STORE CHANNEL NAME 
            }
            return(channel)
          }
        )
      )
      # MATCH FILENAMES TO MARKER/CHANNEL COMBOS
      file_channel_match <- structure(
        lapply(
          seq_along(file_names),
          function(z) {
            # SPLIT FILE NAME INTO INDIVIDUAL CHARACTERS
            file_name_split <- strsplit(file_names[[z]], "")[[1]]
            # LOOP THROUGH MARKER/CHANNEL COMBOS
            channel_match_opts <- structure(
              lapply(
                seq_along(channels),
                function(q) {
                  # SPLIT ALPHNUMERIC CHARACTERS - NO SUFFIX
                  channel_split <- strsplit(
                    gsub(
                      "[^[:alnum:]]",
                      "",
                      gsub(
                        "-A$|-H$|-W$",
                        "",
                        channels[q]
                      ),
                    ),
                    ""
                  )[[1]]
                  # MATCHES PER CHARACTER
                  channel_split_match <- structure(
                    lapply(
                      channel_split,
                      function(char) {
                        grep(
                          char, 
                          file_name_split, 
                          ignore.case = ignore.case, 
                          ...
                        )
                      }
                    ),
                    names = channel_split
                  )
                  # NUMBER OF CHARACTERS
                  n <- length(channel_split_match)
                  # STORE BEST MATCHES PER STARTING CHARACTER
                  cm <- list()
                  # LOOP THROUGH CHARACTER MATCHES
                  for(i in seq_len(n)) {
                    # BEST MATCH PER STARTING LETTER
                    m <- c()
                    # BREAK LOOP - LONGER MATCH IMPOSSIBLE
                    if(length(cm) > 0) {
                      if((n-i) < max(
                        unlist(
                          lapply(
                            cm,
                            function(g) {
                              length(g[g > 0])
                            }
                          )
                        )
                      )) {
                        break()
                      }
                    }
                    # FIRST CHARACTER MUST MATCH
                    if(length(channel_split_match[[i]])> 0) {
                      # LOOP THROUGH STARTING CHARACTER OPTIONS
                      for(v in channel_split_match[[i]]) {
                        # INITIATE WITH STARTING LETTER
                        m_new <- structure(
                          v,
                          names = names(channel_split_match[i])
                        )
                        # DESCEND TREE - SEARCH FOR VALID INDEX
                        for(q in (i + 1):n) {
                          if(q <= n) {
                            # NO MATCH FOUND
                            if(length(channel_split_match[[q]]) == 0) {
                              # SPACER
                              m_new <- c(
                                m_new,
                                structure(
                                  0,
                                  names = names(channel_split_match[q])
                                )
                              )
                              # MATCHES FOUND
                            } else {
                              # MINIMUM INDEX
                              if(length(m_new) == 0) {
                                m_new_ind <- 1
                              } else {
                                m_new_ind <- max(m_new)
                              }
                              # CHECK FOR VALID MATCHES
                              w <- channel_split_match[[q]][
                                channel_split_match[[q]] > m_new_ind
                              ]
                              # MATCH LOCATED
                              if(length(w) > 0) {
                                # NO INSERTIONS
                                if(!insertions) {
                                  if(!any(w == m_new_ind + 1)) {
                                    break()
                                  } else {
                                    w <- w[w == m_new_ind + 1][1]
                                  }
                                }
                                # UPDATE M
                                m_new <- c(
                                  m_new,
                                  structure(
                                    w[1],
                                    names = names(channel_split_match[q])
                                  )
                                )
                                # NO MATCH FOUND
                              } else {
                                # UPDATE M - EMPTY MATCH
                                m_new <- c(
                                  m_new,
                                  structure(
                                    0,
                                    names = names(channel_split_match[q])
                                  )
                                )
                              }
                            }
                          }
                        }
                        # UPDATE M
                        if(length(m_new[m_new > 0]) > length(m[m > 0])) {
                          m <- m_new
                        }
                      }
                    }
                    # UPDATE CM
                    if(length(m[m > 0]) > 0) {
                      cm <- c(
                        cm, 
                        list(m)
                      )
                    }
                  }
                  # RETURN LONGEST MATCH
                  if(length(cm) == 0) {
                    return(0)
                  } else {
                    return(
                      max(
                        unlist(
                          lapply(
                            cm,
                            function(p) {
                              length(p[p>0])
                            }
                          )
                        )
                      )
                    )
                  }
                }
              ),
              names = channels
            )
            # TODO: MUST HAVE AT LEAST 2 CONSECUTIVE CHARACTERS FOR A MATCH
            # UPDATE CHANNEL_MATCH
            if(sum(unlist(channel_match_opts)) > 0) {
              ind <- which(
                channel_match_opts == max(unlist(channel_match_opts))
              )
              # MULTIPLE CHANNEL MATCHES - CHOOSE SHORTEST CHANNEL OPTION
              if(length(ind) > 0) {
                ind <- ind[which.min(nchar(channels[ind]))]
                # CANNOT ASSIGN CHANNELS WITH LENGTH - AMBIGUOUS
                if(length(ind) == 1) {
                  pd[
                    match(names(file_names)[z], rownames(pd)),
                    "channel"
                  ] <<- names(channels)[ind] # NEED TO STORE ORIGINAL CHANNELS
                }
              }
            }
          }
        ),
        names = file_names
      )
    }
  }
  
  # TODO: DECIDE QUALITY OF MATCH & SET TO NA OTHERWISE
  
  # UNMATCHED CHANNELS
  ind <- which(is.na(pd$channel))
  ind <- ind[ind %in% x_ind]
  
  # MATCH CHANNELS BY FLUORESCENT INTENSITIES
  if(length(ind) > 0) {
    message(
      paste0(
        "The following samples could not be matched to a channel by name:",
        "\n",
        paste0(
          pd$name[ind],
          collapse = "\n"
        ),
        "\n",
        "CytoExploreR will make an educated guess for these samples using ",
        "the intensities in the fluorescent channels."
      )
    )
    # EXTRACT DATA - CURRENT SCALE
    cs_list <- structure(
      lapply(
        ind,
        function(z) {
          cyto_data_extract(
            x[z],
            parent = pd$parent[z],
            format = "cytoset",
            copy = FALSE
          )[[1]]
        }
      ),
      names = rownames(pd)[ind]
    )
    # COMPUTE STATISTICS
    chans <- cyto_apply(
      cs_list,
      input = "matrix",
      channels = names(channels),
      FUN = function(z){
        # COMPUTE MEDFI
        res <- cyto_stat_quantile(
          z,
          probs = 0.95
        )
        return(
          names(res)[which.max(res)]
        )
      },
      copy = FALSE
    )
    pd$channel[ind] <- unlist(chans)
  }
  
  # INTERACTIVE EDITING & EXPORT -----------------------------------------------
  
  # INTERACTIVE CHANNEL MATCHING
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    rn <- rownames(pd)
    rownames(pd) <- NULL
    pd <- data_edit(
      pd,
      logo = CytoExploreR_logo(),
      title = "Channel Match Editor",
      row_edit = FALSE,
      col_readonly = "name",
      col_options = list(
        "parent" = if(cyto_class(x, "GatingSet")) {
          cyto_nodes(x, path = "auto")
        } else {
          "root"
        },
        "channel" = c(unname(channels), "Unstained")
      ),
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane"
    )
    rownames(pd) <- rn
  }
  
  # ROWNAMES MISSING IN FILE
  if(is.null(rownames(pd))) {
    rownames(pd) <- pd[, "name"]
  }
  
  # UPDATE DETAILS IN SAMPLES
  cyto_details(x) <- pd[
    match(
      rownames(cyto_details(x)), 
      rownames(pd)
    ), , drop = FALSE]
  
  # SAVE_AS
  if(!.all_na(save_as)) {
    # WRITE TO CSV
    write_to_csv(
      pd,
      save_as,
      row.names = TRUE
    )
  }
  
  # RETURN CHANNEL MATCHING
  return(pd)
  
}

#' Helper function to remove overlapping string fragments
#' @noRd
.cyto_string_strip <- function(x,
                               pad = FALSE) {
  # MULTIPLE STRINGS REQUIRED
  if(length(x) == 1) {
    return(x)
  }
  # MAX SEARCH DEPTH
  depth <- max(nchar(x))
  # LEFT SEARCH
  for(i in 1:depth) {
    if(length(unique(substring(x, 1, 1))) == 1) {
      x <- gsub("^.", "", x)
    } else {
      break()
    }
  }
  # RIGHT SEARCH
  for(i in 1:depth) {
    frag <- unlist(
      lapply(
        x,
        function(z) {
          substring(
            z,
            nchar(z),
            nchar(z)
          )
        }
      )
    )
    if(length(unique(frag)) == 1) {
      x <- gsub(".$", "", x)
    } else {
      break()
    }
  }
  # PADDING - SAME WIDTH
  if(pad) {
    x <- format(x, width = max(nchar(x)))
  }
  # REMAINDER
  return(x)
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
