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
#' @param select channel selection criteria in the form of character strings or
#'   indices.
#' @param exclude channel exclusion criteria in the form of character strings or
#'   indices.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param escape logical to indicate whether special characters should be
#'   escaped prior to performing partial matching for \code{select} and
#'   \code{exclude}, set to TRUE by default.
#' @param ignore.case logical indicating whether case insensitive channel
#'   matches should be returned, set to TRUE bu default.
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
                          escape = TRUE,
                          ignore.case = TRUE,
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
          # INDEX
          if(is.numeric(z)) {
            z
          # NAME
          } else {
            which(
              suppressWarnings(
                .grepl(
                  z, 
                  channels, 
                  escape = escape,
                  fixed = FALSE,
                  ignore.case = ignore.case,
                  ...
                )
              )
            )
          }
        }
      )
    )
  } else {
    ind <- seq_along(channels)
  }
  
  # EXCLUDE
  if(!is.null(exclude)){
    ind_rm <- unique(
      LAPPLY(
        exclude, 
        function(z){
          # INDEX
          if(is.numeric(z)) {
            z
            # NAME
          } else {
            which(
              suppressWarnings(
                .grepl(
                  z, 
                  channels, 
                  escape = escape,
                  fixed = FALSE,
                  ignore.case = ignore.case,
                  ...
                )
              )
            )
          }
        }
      )
    )
  } else {
    ind_rm <- NULL
  }
  
  # CHANNELS EXCLUDE
  if(length(ind_rm) > 0) {
    ind <- ind[!ind %in% ind_rm]
  }
  
  # SUBSET CHANNELS
  channels <- channels[ind]
  
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
#' @param select vector of channels, markers or indices for which the
#'   channel/marker combinations should be returned.
#' @param exclude vector of channels, markers or indices for which the
#'   channel/marker combinations should not be returned.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param escape logical to indicate whether special characters should be
#'   escaped prior to performing partial matching for \code{select} and
#'   \code{exclude}, set to TRUE by default.
#' @param ignore.case logical indicating whether case insensitive channel
#'   matches should be returned, set to TRUE bu default.
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
                         escape = TRUE,
                         ignore.case = TRUE,
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
      ind <- unique(
        LAPPLY(
          select, 
          function(z) {
            if(is.numeric(z)) {
              z
            } else {
              which(
                suppressWarnings(
                  .grepl(
                    z,
                    markers,
                    escape = escape,
                    ignore.case = ignore.case,
                    ...
                  )
                ) |
                suppressWarnings(
                  .grepl(
                    z, 
                    names(markers),
                    escape = escape,
                    ignore.case = ignore.case,
                    ...
                  )
                )
              )
            }
          }
        )
      )
    } else {
      ind <- seq_along(markers)
    }
    
    # EXCLUDE
    if(!is.null(exclude)){
      ind_rm <- unique(
        LAPPLY(
          exclude, 
          function(z) {
            if(is.numeric(z)) {
              z
            } else {
              which(
                suppressWarnings(
                  .grepl(
                    z,
                    markers,
                    escape = escape,
                    ignore.case = ignore.case,
                    ...
                  )
                ) |
                suppressWarnings(
                  .grepl(
                    z, 
                    names(markers),
                    escape = escape,
                    ignore.case = ignore.case,
                    ...
                  )
                )
              )
            }
          }
        )
      )
    } else {
      ind_rm <- NULL
    }
    
    # MARKERS EXCLUDE
    if(length(ind_rm) > 0) {
      ind <- ind[!ind %in% ind_rm]
    }
    
    # SUBSET MARKERS
    markers <- markers[ind]
    
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

## CYTO_MARKERS REPLACEMENT METHOD ---------------------------------------------

#' Replace marker names
#' 
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSey}}.
#' @param value named vector of channels with their associated markers.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_markers}}
#' @seealso \code{\link{cyto_markers_edit}}
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
#' # Set FSC-A marker to FSC
#' cyto_markers(gs) <- c("FSC" = "FSC-A")
#' 
#' @export
"cyto_markers<-" <- function(x, value) {
  # CHECK VALUE
  if(is.null(names(value))) {
    stop(
      "Values to replace must be a named vector of channels with their markers!"
    )
  }
  # TODO: CHACK VALID CHANNELS HAVE BEEN SUPPLIED
  # PREPARE FLIPPED VALUE
  if(all(value %in% cyto_channels(x))) {
    value <- structure(
      names(value),
      names = value
    )
  }
  # UPDATE MARKER ASSIGNMENTS
  if(cyto_class(x, c("flowFrame", "flowSet"), TRUE)) {
    flowCore::markernames(x) <- value
  } else if(cyto_class(x, c("cytoframe", "cytoset", "GatingSet"))) {
    flowWorkspace::markernames(x) <- value
  } else {
    stop(
      paste0(
        "Cannot replace marker names for objects of class ", cyto_class(x), "!"
      )
    )
  }
  return(x)
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
  args <- list(
    "x" = x,
    ...
  )
  args$exclude <- c(
    args$exclude,
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
  )
  
  # CALL CYTO_CHANNELS()
  cyto_func_call(
    "cyto_channels",
    args
  )

}

## CYTO_UNMIX_CHANNELS ---------------------------------------------------------

#' Extract unmixed channels
#' 
#' @param x object of class cytoframe, cytoset, GatingHierrachy or GatingSet 
#' potentially containing unmixed parameters.
#' 
#' @return the names of the unmixed parameters or NULL.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @examples
#' \dontrun{
#' cyto_unmix_channels(gs)
#' }
#' 
#' 
#' @export
cyto_unmix_channels <- function(x) {
  
  # NOTE: CHECK FOR PNTYPE = UNMIXED FLUORESCENCE
  if(cyto_class(x, c("FlowSet", "cytoset", "GatingSet"), TRUE)) {
    x <- x[[1]]
  }
  kw <- cyto_keyword(x)
  kw <- unlist(kw[grepl("\\$P[0-9]{1,3}TYPE", names(kw))])
  kw <- kw[grep("Unmixed_Fluorescence", kw, ignore.case = TRUE)]
  if(length(kw) == 0) {
    return(NULL)
  }
  # EXTRACT UNMIXED PARAMETERS
  return(
    unname(
      unlist(
        cyto_keyword(
          x
        )[gsub(
          "TYPE",
          "N", 
          names(kw)
        )]
      )
    )
  )
  
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
#'   700-A","CD8")) or indices.
#' @param skip vector of markers/channels in \code{channels} to bypass when
#'   converting to valid channels, for example \code{"Unstained"} is bypassed
#'   when checking channels in the channel match file.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to FALSE by default. If set to TRUE an additional check will be
#'   performed to ensure that only 1 or 2 \code{channels} are supplied.
#' @param escape logical to indicate whether special characters should be
#'   escaped prior to performing partial channel matching, set to TRUE by
#'   default.
#' @param ignore.case logical indicating whether case insensitive channel
#'   matches should be returned, set to TRUE by default.
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
                                  append = FALSE,
                                  plot = FALSE,
                                  escape = TRUE,
                                  ignore.case = TRUE,
                                  ...) {
  
  # CHANNELS
  chans <- cyto_channels(x)
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # EXTRACT CHANNELS
  res <- c()
  for(z in seq_along(channels)) {
    # SKIP
    skip <- any(
      .grepl(
        channels[z],
        skip,
        escape = escape,
        ignore.case = ignore.case,
        ...
      )
    )
    if(skip) {
      res <- c(res, channels[z])
    # INDEX
    } else if(is.numeric(channels[z])) {
      res <- c(res, chans[channels[z]]) 
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
          .grepl(
            channels[z],
            markers,
            escape = escape,
            ignore.case = ignore.case,
            ...
          )
        )
      )
      channel_ind <- suppressWarnings(
        which(
          .grepl(
            channels[z],
            chans,
            escape = escape,
            ignore.case = ignore.case,
            ...
          )
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
          paste0(
            channels[z],
            " is not a valid channel or marker for this ", 
            cyto_class(x, class = TRUE),
            "!"
          )
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
  
  # APPEND
  if(append) {
    res <- paste0(
      "<",
      names(res),
      "> ",
      res
    )
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
#'   700-A","CD8")) or channel indices.
#' @param skip vector of markers/channels in \code{channels} to bypass when
#'   converting to valid markers, for example \code{"Unstained"} is bypassed
#'   when checking markers in the channel match file.
#' @param append logical indicating whether the name of the channel should be
#'   appended to the marker names in the form \code{<marker> channel}, set to
#'   FALSE by default.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to FALSE by default. If set to TRUE an additional check will be
#'   performed to ensure that only 1 or 2 \code{channels} are supplied.
#' @param escape logical to indicate whether special characters should be
#'   escaped prior to performing partial marker matching, set to TRUE by
#'   default.
#' @param ignore.case logical indicating whether case insensitive channel
#'   matches should be returned, set to TRUE bu default.
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
                                 escape = TRUE,
                                 ignore.case = TRUE,
                                 ...) {
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # CHANNELS
  chans <- cyto_channels(x)
  
  # EXTRACT MARKERS
  res <- c()
  for(z in seq_along(channels)) {
    # SKIP
    skip <- any(
      .grepl(
        channels[z],
        skip,
        escape = escape,
        ignore.case = ignore.case,
        ...
      )
    )
    if(skip) {
      res <- c(res, channels[z])
      names(res[length(res)]) <- channels[z] # append
    # INDEX
    } else if(is.numeric(channels[z])) {
      res <- c(res, markers[match(chans[channels[z]], names(markers))])
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
          .grepl(
            channels[z],
            markers,
            escape = escape,
            ignore.case = ignore.case,
            ...
          )
        )
      )
      channel_ind <- suppressWarnings(
        which(
          .grepl(
            channels[z],
            names(markers),
            escape = escape,
            ignore.case = ignore.case,
            ...
          )
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
        } else if(
          any(
            .grepl(
              channels[z],
              chans, 
              escape = escape,
              ignore.case = ignore.case,
              ...
            )
          )
        ) {
          res <- c(
            res,
            structure(
              chans[
                which(
                  .grepl(
                    channels[z],
                    chans,
                    escape = escape,
                    ignore.case = ignore.case,
                    ...
                  )
                )
              ],
              names = chans[
                which(
                  .grepl(
                    channels[z],
                    chans,
                    escape = escape,
                    ignore.case = ignore.case,
                    ...
                  )
                )
              ]
            )
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
  chans <- data.frame(
    "name" = cyto_names(x),
    "channel" = NA,
    stringsAsFactors = FALSE
  )
  
  # CHANNEL SELECTION
  chans <- data_edit(
    chans,
    title = "Channel Selector",
    logo = CytoExploreR_logo(),
    col_edit = FALSE,
    row_edit = FALSE,
    col_options = list("channel" = unname(opts)),
    col_names = "channel",
    col_readonly = "name",
    hide = TRUE,
    quiet = TRUE,
    viewer = "pane"
  )
  
  # MISSING CHANNELS
  lapply(
    seq_along(chans[, "channel"]),
    function(z) {
      if(is.na(chans[z, "channel"])) {
        stop(paste0("No channel selected for ", chans[z, "name"], "."))
      }
    }
  )
  
  # RETURN VECTOR OF CHANNELS
  return(chans[, "channel"])
}

## CYTO_CHANNEL_MATCH ----------------------------------------------------------

#' Match each single colour control to a fluorescent channel
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels names of the possible channels or markers to be matched to
#'   each control, set to all area fluorescent parameters by default.
#' @param peaks logical to indicate whether peak detectors should be
#'   automatically assigned using \code{cyto_peaks_assign()}, set to FALSE by
#'   default.
#' @param file name of a CSV file from which the channel matching information
#'   should be inherited. If not supplied, \code{cyto_channel_match()} will
#'   automatically search for a file named \code{"Control-Details.csv"} or
#'   create interactively create such a file.
#' @param save_as name of a CSV file to which the channel matching should be
#'   written for downstream use, set to \code{"Control-Details.csv"}
#'   prefixed with the date by default. Users can set custom file names here,
#'   but the file name should contain \code{"Control-Details.csv"} in order
#'   to be automatically detected by CytoExploreR within
#'   \code{cyto_spillover_compute()}, \code{cyto_spillover_edit()},
#'   \code{cyto_spillover_spread()} and \code{cyto_plot_compensation()}.
#'
#' @return a data.frame written to a CSV file containing information matching
#'   each file name to a channel. This channel matching information is also
#'   automatically added to the \code{cyto_details()} of the supplied samples
#'   where it can be easily accessed by CytoExploreR downstream.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_peaks_assign}}
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
                               peaks = FALSE,
                               file = NULL,
                               save_as = NULL,
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
    channel_match <- args[["channel_match"]]
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
    channels <- channels[
      !grepl(
        "\\-H$|\\-W$",
        channels,
        ignore.case = TRUE
      )
    ]
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # IMPORT DATA FROM FILE
  pd_new <- cyto_file_search(
    "Details.*\\.csv$|Channel-Match.*\\.csv",
    rownames = rownames(pd),
    ignore.case = TRUE,
    data.table = FALSE,
    type = "channel match details",
    files = file
  )
  
  # DETAILS LOCATED IN FILE
  if(!is.null(pd_new)) {
    # SAVE TO IMPORTED FILE
    if(is.null(file)) {
      file <- names(pd_new)
    }
    pd <- pd_new[[1]]
  }
  
  # DEFAULT FILE NAME
  if(is.null(save_as)) {
    save_as <- file
    if(is.null(save_as)) {
      save_as <- cyto_file_name(
        paste0(
          format(
            Sys.Date(), 
            "%d%m%y"
          ), 
          "-Control-Details.csv"
        )
      )
    }
  }
  
  # ADD MISSING PARAMETERS
  vars <- c("group", "parent", "channel", "marker", "label", "select")
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
  
  # TODO: IMPROVE PARENT ASSIGNMENT SO PEAK DTECTOR ASSIGNMENTS ARE CORRECT
  
  # GROUPS/PARENTS IGNORED FOR CYTOSETS
  if(cyto_class(x, "GatingSet")) {
    # SAMPLES WITH PARENTS MISSING
    ind <- which(is.na(pd$parent[x_ind]))
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
        input = "cytoset",
        FUN = "cyto_stat_count",
        copy = FALSE
      )
      if(cyto_class(pop_stats, "list", TRUE)) {
        pop_stats <- do.call("cbind", pop_stats)
        dimnames(pop_stats) <- list(rownames(pop_stats), pops)
      }
      # PARENT - TERMINAL NODE MOST EVENTS
      pd$parent[x_ind[ind]] <- pops[
        apply(
          pop_stats,
          1,
          "which.max"
        )
      ]
    }
    # GROUPS MISSING - DEFAULT TO PARENTS
    ind <- which(is.na(pd$group[x_ind]))
    if(length(ind) > 0) {
      pd$group[x_ind[ind]] <- pd$parent[x_ind[ind]]
    }
  }
  
  # CHANNELS -------------------------------------------------------------------
  
  # LOCATE UNSTAINED CONTROLS
  unst_idx <- grep(
    "Unstained|NIL",
    cyto_names(x),
    ignore.case = TRUE
  )
  
  # SET UNSTAINED CHANNEL
  pd$channel[unst_idx][pd$channel[unst_idx] %in% c(NA, "NA")] <- "Unstained"
  
  # UPDATE METADATA
  cyto_details(x) <- pd[match(rownames(cyto_details(x)), rownames(pd)), ,]
  
  # ASSIGN PEAK DETECTORS
  if(peaks) {
    x <- cyto_peaks_assign(
      x,
      overwrite = FALSE
    )
  }
  
  # UPDATE METADATA
  pd <- cyto_details(x)
  
  # CONVERT SELECT COLUMN TO LOGICAL FOR RENDERING
  pd$select <- as.logical(pd$select)
  
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
        "channel" = c(unname(channels), "Unstained"),
        "select" = c(TRUE, FALSE)
      ),
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane"
    )
    rownames(pd) <- rn
  }
  
  # SET FALSE OPTIONS FOR SELECT COLUMN
  pd[is.na(pd[, "select"]), "select"] <- FALSE
  
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
  channels_exempt <- cyto_channels(
    x, 
    select = c("FSC",
               "SSC",
               "Time",
               "^Event-?ID$",
               "^Sample-?ID$")
  )
  
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
    channels_to_remove <- cyto_markers(
      x, 
      select = exclude,
      ...
    )
    if(length(channels_to_remove) > 0){
      channels_to_remove <- names(channels_to_remove)
    }else{
      channels_to_remove <- c()
    }
    
    # CHANNELS TO REMOVE
    channels_to_remove <- c(
      channels_to_remove,
      cyto_channels(
        x,
        select = exclude,
        ...
      )
    )
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
