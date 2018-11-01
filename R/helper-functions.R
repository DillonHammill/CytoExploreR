#' Extract Fluorescent Channels
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{getChannels,flowFrame-method}}
#' @seealso \code{\link{getChannels,flowSet-method}}
#' @seealso \code{\link{getChannels,GatingSet-method}}
#'
#' @export
setGeneric(name="getChannels",
           def=function(x){standardGeneric("getChannels")}
)

#' Extract Fluorescent Channels - flowFrame Method
#' 
#' @param x object \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' 
#' @return vector of fluorescent channels.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{getChannels,flowSet-method}}
#' @seealso \code{\link{getChannels,GatingSet-method}}
#' 
#' @export
setMethod(getChannels, signature = "flowFrame", definition = function(x){
  
  colnames(x@description$SPILL)
  
})

#' Extract Fluorescent Channels - flowSet Method
#' 
#' @param x object \code{\link[flowCore:flowSet-class]{flowSet}}.
#' 
#' @return vector of fluorescent channels.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{getChannels,flowFrame-method}}
#' @seealso \code{\link{getChannels,GatingSet-method}}
#' 
#' @export
setMethod(getChannels, signature = "flowSet", definition = function(x){
  
  colnames(x[[1]]@description$SPILL)
  
})

#' Extract Fluorescent Channels - GatingSet Method
#' 
#' @param x object \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' 
#' @return vector of fluorescent channels.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{getChannels,flowFrame-method}}
#' @seealso \code{\link{getChannels,flowSet-method}}
#' 
#' @export
setMethod(getChannels, signature = "GatingSet", definition = function(x){
  
  colnames(x@data[[1]]@description$SPILL)
  
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
setGeneric(name="selectChannels",
           def=function(x){standardGeneric("selectChannels")}
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
#' @export
setMethod(selectChannels, signature = "flowFrame", definition = function(x){
  
  # Assign x to fr
  fr <- x
  
  opts <- getChannels(fr)
  
  # Print sample name and select channel
  message(paste("Select a fluorescent channel for the following compensation control:", fr@description$GUID))
  channel <- opts[menu(choices = opts, graphics = TRUE)]
  
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
#' @export
setMethod(selectChannels, signature = "flowSet", definition = function(x){
  
  # Assign x to fs
  fs <- x
  
  opts <- c(getChannels(fs), "Unstained")
  
  # Print sample name and select channel
  channels <- opts[sapply(pData(fs)$name, function(x){
    
    message("Select a fluorescent channel for the following compensation control:")
    
    print(x)
    
    menu(choices = opts, graphics = TRUE)
    
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
#' @export
setMethod(selectChannels, signature = "GatingSet", definition = function(x){
  
  # Assign x to gs
  gs <- x
  
  opts <- c(getChannels(gs), "Unstained")
  
  # Print sample name and select channel
  channels <- opts[sapply(pData(gs)$name, function(x){
    
    message("Select a fluorescent channel for the following compensation control:")
    
    print(x)
    
    menu(choices = opts, graphics = TRUE)
    
  })]
  
  return(channels)
  
})

#' Sample a flowFrame
#'
#' @param fr object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param size numeric indicating the number of events to keep. If \code{size >
#'   nrow(fr)}, size is set to \code{nrow(fr)}.
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
#' @export
sampleFrame <- function(fr, size = 250000){
  
  # Number of events
  events <- nrow(fr)
  
  if(events < size){
    
    size <- events
  
  }else{
    
  }

  smp <- sampleFilter(size = size)
  fr <- Subset(fr, smp)
  
  return(fr)
}

#' Select flowFrames Based on pData
#'
#' @param fs an object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param pData vector of form \code{c("column","row")} indicating the rows of
#'   \code{pData(fs)} to extract.
#'
#' @return an object of calss \code{\link[flowCore:flowSet-class]{flowSet}}
#'   containing the extracted \code{\link[flowCore:flowFrame-class]{flowFrame}}
#'   objects.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
selectFrames <- function(fs, pData){
  
  # Extract pData from flowSet 
  pd <- pData(fs)
  
  # Find which rows to extract
  rows <- pd[[pData[1]]] == pData[2]
  rows[which(is.na(rows))] <- FALSE
  
  # Get sampleNames of these rows - exclude NA
  sn <- rownames(pd)[rows]
  
  # Extract these samples by name from flowSet
  fs <- fs[sn]
  
  return(fs)
  
}
