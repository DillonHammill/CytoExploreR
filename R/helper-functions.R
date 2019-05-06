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
setMethod(cyto_fluor_channels,
          signature = "flowFrame",
          definition = function(x) {
            channels <- unname(BiocGenerics::colnames(x))
            channels <- channels[!channels %in% c(
              "FSC-A",
              "FSC-H",
              "FSC-W",
              "SSC-A",
              "SSC-H",
              "SSC-W",
              "Time",
              "Original"
            )]
            
            return(channels)
          }
)

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
setMethod(cyto_fluor_channels,
          signature = "flowSet",
          definition = function(x) {
            cyto_fluor_channels(x[[1]])
          }
)

#' Extract Fluorescent Channels - GatingHierarchy Method
#'
#' @param x object
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}}.
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
#' cyto_fluor_channels(gs[[1]])
#' @export
setMethod(cyto_fluor_channels,
          signature = "GatingHierarchy",
          definition = function(x) {
            fr <- getData(x, "root")
            cyto_fluor_channels(fr)
          }
)

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
setMethod(cyto_fluor_channels,
          signature = "GatingSet",
          definition = function(x) {
            fr <- getData(x[[1]], "root")
            cyto_fluor_channels(fr)
          }
)

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
#' @importFrom flowCore identifier
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
setMethod(cyto_channel_select,
          signature = "flowFrame",
          definition = function(x) {
            
            # Assign x to fr
            fr <- x
            
            opts <- cyto_fluor_channels(fr)
            
            # Print sample name and select channel
            message(
              paste(
                "Select a fluorescent channel for the following sample:",
                identifier(fr)
              )
            )

            chan <- opts[menu(choices = opts, graphics = TRUE)]

            return(chan)
          }
)

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
setMethod(cyto_channel_select,
          signature = "flowSet",
          definition = function(x) {
            
            # Assign x to fs
            fs <- x
            
            opts <- c(cyto_fluor_channels(fs), "Unstained")
            
            # Print sample name and select channel
            chans <- opts[unlist(lapply(pData(fs)$name, function(x) {
              message(
                paste("Select a fluorescent channel for the following sample:",
                      x))
              
              menu(choices = opts, graphics = TRUE)
              
            }))]
            
            return(chans)
          }
)

#' Select Fluorescent Channel for Compensation Controls - GatingHierarchy Method
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}}
#'   containing compensation controls.
#'
#' @return vector of channels in order of GatingSet.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom utils menu
#' @importFrom flowWorkspace getData
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
#' cyto_channel_select(gs[[1]])
#' }
#'
#' @export
setMethod(cyto_channel_select,
          signature = "GatingHierarchy",
          definition = function(x) {
            
            # Assign x to gs
            gh <- x
            
            # Extract flowFrame
            fr <- getData(gh, "root")
            
            opts <- cyto_fluor_channels(fr)
            
            # Print sample name and select channel
            message(
              paste(
                "Select a fluorescent channel for the following sample:",
                identifier(fr)
              )
            )
            
            chan <- opts[menu(choices = opts, graphics = TRUE)]
            
            return(chan)
            
            
          }
)

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
#' @importFrom flowWorkspace getData
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
setMethod(cyto_channel_select,
          signature = "GatingSet",
          definition = function(x) {
            
            # Assign x to gs
            gs <- x
            
            # Extract flowSet
            fs <- getData(gs, "root")
            
            # Channel options
            opts <- c(cyto_fluor_channels(fs), "Unstained")
            
            # Print sample name and select channel
            chans <- opts[unlist(lapply(pData(fs)$name, function(x) {
              message(
                paste("Select a fluorescent channel for the following sample:",
                      x))
              
              menu(choices = opts, graphics = TRUE)
              
            }))]
            
            return(chans)
            
            
          }
)
