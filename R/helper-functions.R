#' Extract Fluorescent Channels
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels,flowFrame-method}}
#' @seealso \code{\link{cyto_fluor_channels,flowSet-method}}
#' @seealso \code{\link{cyto_fluor_channels,GatingSet-method}}
#'
#' @export
setGeneric(
  name = "cyto_fluor_channels",
  def = function(x) {
    standardGeneric("cyto_fluor_channels")
  }
)

#' Extract Fluorescent Channels - flowFrame Method
#'
#' @param x object \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#'
#' @return vector of fluorescent channels.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels,flowSet-method}}
#' @seealso \code{\link{cyto_fluor_channels,GatingSet-method}}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Get fluorescent channels
#' cyto_fluor_channels(fs[[1]])
#' @export
setMethod(cyto_fluor_channels, signature = "flowFrame", definition = function(x) {
  channels <- unname(BiocGenerics::colnames(x))
  channels <- channels[!channels %in% c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W", "Time", "Original")]

  return(channels)
})

#' Extract Fluorescent Channels - flowSet Method
#'
#' @param x object \code{\link[flowCore:flowSet-class]{flowSet}}.
#'
#' @return vector of fluorescent channels.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels,flowFrame-method}}
#' @seealso \code{\link{cyto_fluor_channels,GatingSet-method}}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # get fluorescent channels
#' cyto_fluor_channels(fs)
#' @export
setMethod(cyto_fluor_channels, signature = "flowSet", definition = function(x) {
  cyto_fluor_channels(x[[1]])
})

#' Extract Fluorescent Channels - GatingSet Method
#'
#' @param x object \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return vector of fluorescent channels.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels,flowFrame-method}}
#' @seealso \code{\link{cyto_fluor_channels,flowSet-method}}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#' 
#' # Get fluorescent channels
#' cyto_fluor_channels(gs)
#' @export
setMethod(cyto_fluor_channels, signature = "GatingSet", definition = function(x) {
  fr <- getData(x[[1]], "root")
  cyto_fluor_channels(fr)
})

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
#' @export
setGeneric(
  name = "cyto_channel_select",
  def = function(x) {
    standardGeneric("cyto_channel_select")
  }
)

#' Select Fluorescent Channel for Compensation Controls - flowFrame Method
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#'
#' @return selected channel associated with the supplied flowFrame.
#'
#' @importFrom utils menu
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Select a channel from dropdown menu
#' cyto_channel_select(fs[[1]])
#' }
#' 
#' @export
setMethod(cyto_channel_select, signature = "flowFrame", definition = function(x) {

  # Assign x to fr
  fr <- x

  opts <- cyto_fluor_channels(fr)

  # Print sample name and select channel
  message(paste("Select a fluorescent channel for the following compensation control:", fr@description$GUID))

  if (getOption("CytoRSuite_interact") == TRUE) {
    channel <- opts[menu(choices = opts, graphics = TRUE)]
  } else {
    # Tests use PE Cy7 Control -
    channel <- opts[5]
  }

  return(channel)
})

#' Select Fluorescent Channel for Compensation Controls - flowSet Method
#'
#' @param x object of class
#'   \code{\link[flowCore:flowSet-class]{flowSet}} containing compensation
#'   controls.
#'
#' @return vector of channels in order of flowSet.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom utils menu
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Select channel for each sample from dropdown menu
#' cyto_channel_select(fs)
#' }
#' 
#' @export
setMethod(cyto_channel_select, signature = "flowSet", definition = function(x) {

  # Assign x to fs
  fs <- x

  opts <- c(cyto_fluor_channels(fs), "Unstained")

  # Print sample name and select channel
  channels <- opts[sapply(pData(fs)$name, function(x) {
    message("Select a fluorescent channel for the following compensation control:")

    print(x)

    if (getOption("CytoRSuite_interact") == TRUE) {
      menu(choices = opts, graphics = TRUE)
    } else {

      # Test channels - 7AAD, AF430, APC Cy7, NIL, PE Cy7, PE
      c(4, 7, 11, 12, 5, 2)[match(x, pData(fs)$name)]
    }
  })]

  return(channels)
})

#' Select Fluorescent Channel for Compensation Controls - GatingSet Method
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#'   compensation controls.
#'
#' @return vector of channels in order of GatingSet.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom utils menu
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#' 
#' # Select channel for each sample from dropdown menu
#' cyto_channel_select(gs)
#' }
#' 
#' @export
setMethod(cyto_channel_select, signature = "GatingSet", definition = function(x) {

  # Assign x to gs
  gs <- x

  opts <- c(cyto_fluor_channels(gs), "Unstained")

  # Print sample name and select channel
  channels <- opts[sapply(pData(gs)$name, function(x) {
    message("Select a fluorescent channel for the following compensation control:")

    print(x)

    if (getOption("CytoRSuite_interact") == TRUE) {
      menu(choices = opts, graphics = TRUE)
    } else {

      # Test channels - 7AAD, AF430, APC Cy7, NIL, PE Cy7, PE
      c(4, 7, 11, 12, 5, 2)[match(x, pData(gs)$name)]
    }
  })]

  return(channels)
})

#' Sample a flowFrame
#'
#' @param fr object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param display numeric [0,1] indicating the percentage of events to plot.
#'
#' @return \code{\link[flowCore:flowFrame-class]{flowFrame}} restricted to
#'   \code{size} events.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter
#' @importFrom flowCore Subset
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Constrict first sample by 50%
#' cyto_sample(fs[[1]], 0.5)
#' @export
cyto_sample <- function(fr, display) {

  # Number of events
  events <- nrow(fr)

  # Size
  size <- display * events

  smp <- sampleFilter(size = size)
  fr <- Subset(fr, smp)

  return(fr)
}

#' Assign marker names to flowFrame or flowSet
#'
#' \code{markers_assign} opens an editable table containing a list of channels
#' and markers for a \code{flowFrame} or \code{flowSet}. Users can edit the
#' \code{markers} column as required and these entries will be updated in the
#' \code{flowFrame} or \code{flowSet} upon closing the window.
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters markernames
#' @importFrom utils edit
#'
#' @return NULL and update marker names of \code{x}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Add marker names to channels - edit table
#' cyto_markers(fs)
#' }
#'
#' @export
cyto_markers <- function(x) {
  if (!any(inherits(x, "flowFrame") | inherits(x, "flowSet"))) {
    stop("Please supply either a flowFrame or flowSet object")
  }

  if (inherits(x, "flowFrame")) {

    # Extract pData of parameters
    pd <- pData(parameters(x))
  } else if (inherits(x, "flowSet")) {

    # Extract pData of parameters
    pd <- pData(parameters(x[[1]]))
  }

  # Make data.frame with channel and marker column
  dt <- pd[, c("name", "desc")]
  colnames(dt) <- c("Channel", "Marker")
  rownames(dt) <- NULL

  # Channels with markers
  chans <- as.vector(dt$Channel[!is.na(dt$Marker)])

  # Edit dt
  dt <- suppressWarnings(edit(dt))

  # Channels with markers added
  tb <- dt[!dt$Channel %in% chans, ]

  chns <- tb$Channel[!is.na(tb$Marker)]

  # Pull out assigned markers
  mrk <- dt$Marker[dt$Channel %in% c(chans, chns)]
  names(mrk) <- dt$Channel[dt$Channel %in% c(chans, chns)]

  # Assign markers to x
  markernames(x) <- mrk

  invisible(NULL)
}

#' Interactively edit pData Information for a flowSet or GatingSet
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#'
#' @return NULL and update pData for the \code{flowSet} or \code{GatingSet}.
#'
#' @importFrom flowWorkspace pData
#' @importFrom utils edit
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Edit pData in table editor
#' cyto_annotate(fs)
#' }
#' 
#' @export
cyto_annotate <- function(x) {

  # x should be a flowSet or GatingSet
  if (!any(inherits(x, "flowSet") | inherits(x, "GatingSet"))) {
    stop("Please supply either a flowSet or a GatingSet")
  }

  # Assign x to cyto
  cyto <- x

  # Extract pData
  pd <- pData(cyto)
  rownames(pd) <- NULL

  # Edit pData
  pd <- edit(pd)
  rownames(pd) <- pd$name

  # Update pData
  pData(cyto) <- pd

  assign(deparse(substitute(x)), cyto, envir = globalenv())

  invisible(NULL)
}
