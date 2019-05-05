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

#' Sample a flowFrame
#'
#' @param fr object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param display numeric [0,1] indicating the percentage of events to plot.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
#'
#' @return \code{\link[flowCore:flowFrame-class]{flowFrame}} restricted to
#'   \code{display} percentage events.
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
#' 
#' @export
cyto_sample <- function(fr, display = 1, seed = NULL) {
  
  # No sampling
  if(display == 1){
    
  }else{
    # Number of events
    events <- nrow(fr)
  
    # Size
    size <- display * events
    
    # Set seed
    if(!is.null(seed)){
      set.seed(seed)
    }

    # Sample
    smp <- sampleFilter(size = size)
    fr <- Subset(fr, smp)
  }

  return(fr)
}

#' Assign marker names to flowFrame or flowSet
#'
#' \code{cyto_markers} opens an editable table containing a list of channels and
#' markers for a \code{flowFrame} or \code{flowSet}. Users can edit the
#' \code{name} or \code{desc} columns with updated channel names or marker names
#' respectively. These entries will be updated in the \code{flowFrame} or
#' \code{flowSet} upon closing the window and saved to a
#' "Experiment-markers.csv" file for future use.
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#' @param file name of csv file containing columns 'Channel' and 'Marker'.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters markernames markernames<-
#' @importFrom utils edit write.csv read.csv
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
cyto_markers <- function(x, file = NULL) {
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
  
  # file missing
  if (is.null(file)) {
    
    if(length(grep("Experiment-Markers.csv", list.files())) != 0){
      
      message("Experiment-Markers.csv found in working directory.")
      dt <- read.csv(list.files()[grep("Experiment-Markers.csv", list.files())],
                     header = TRUE, stringsAsFactors = FALSE)
      
    }else{
      
      # Make data.frame with channel and marker columns
      dt <- pd[, c("name", "desc")]
      colnames(dt) <- c("Channel", "Marker")
      rownames(dt) <- NULL
      
    }
    
  } else {
    if (getOption("CytoRSuite_wd_check")) {
      if (!.file_wd_check(file)) {
        stop(paste(file, "does not exist in this working directory."))
      }
    }
    dt <- read.csv(file, header = TRUE)
  }
  
  # Channels with markers
  chans <- as.vector(dt$Channel[!is.na(dt$Marker)])
  
  # Edit dt
  dt <- suppressWarnings(edit(dt))
  
  # Update channels
  BiocGenerics::colnames(x) <- dt$Channel
  
  # Write result to csv file
  if (length(grep("Experiment-Markers.csv", c(file, list.files()))) != 0) {
    
    write.csv(dt, c(file, list.files())[grep("Experiment-Markers.csv",
                                             c(file, list.files()))[1]]
              , row.names = FALSE)
    
  } else {
    write.csv(dt, paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Markers.csv"
    ), row.names = FALSE)
  }
  
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
#' @param file name of csv file containing experimental information.
#'
#' @return NULL and update pData for the \code{flowSet} or \code{GatingSet}.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore pData<-
#' @importFrom utils edit write.csv read.csv
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
cyto_annotate <- function(x, file = NULL) {
  
  # x should be a flowSet or GatingSet
  if (!any(inherits(x, "flowSet") | inherits(x, "GatingSet"))) {
    stop("Please supply either a flowSet or a GatingSet")
  }
  
  # Assign x to cyto
  cyto <- x
  
  # File missing
  if (is.null(file)) {
    
    if(length(grep("Experiment-Details.csv", list.files())) != 0){
      
      message("Experiment-Details.csv found in working directory.")
      pd <- read.csv(list.files()[grep("Experiment-Details.csv", list.files())],
                     header = TRUE, stringsAsFactors = FALSE)
      
    }else{
      
      # Extract pData
      pd <- pData(cyto)
      rownames(pd) <- NULL
      
    }
    
  } else {
    if (getOption("CytoRSuite_wd_check")) {
      if (!.file_wd_check(file)) {
        stop(paste(file, "does not exist in this working directory"))
      }
    }
    
    # Read in file
    pd <- read.csv(file, header = TRUE)
    
    # Update pData
    pData(cyto) <- pd
  }
  
  # Edit pData
  pd <- suppressWarnings(edit(pd))
  rownames(pd) <- pd$name
  
  # Update pData
  pData(cyto) <- pd
  
  # Write result to csv file
  if (length(grep("Experiment-Details.csv", c(file, list.files()))) != 0) {
    
    write.csv(pd, c(file, list.files())[grep("Experiment-Details.csv",
                                             c(file, list.files()))[1]]
              , row.names = FALSE)
    
  } else {
    
    write.csv(pData(cyto), paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Details.csv"
    ), row.names = FALSE)
    
  }
  
  # Update globally
  assign(deparse(substitute(x)), cyto, envir = globalenv())
  
  invisible(NULL)
}
