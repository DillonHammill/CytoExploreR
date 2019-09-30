#' Extract Fluorescent Channels
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_fluor_channels,flowFrame-method}}
#' @seealso \code{\link{cyto_fluor_channels,flowSet-method}}
#' @seealso \code{\link{cyto_fluor_channels,GatingSet-method}}
#'
#' @examples
#' library(CytoRSuiteData)
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
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Compensation
#'
#' # Select a channel for each control from dropdown menu
#' cyto_channel_select(fs)
#' @name cyto_channel_select
NULL


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
