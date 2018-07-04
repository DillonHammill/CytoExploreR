#' Extract fluorescent parameters
#' 
#' @export
setGeneric(name="getChannels",
           def=function(x){standardGeneric("getChannels")}
)

#' Extract fluorescent parameters from a flowFrame
#' 
#' @param x object \code{flowFrame}
#' 
#' @return vector of fluorescent channels
#' 
#' @export
setMethod(getChannels, signature = "flowFrame", definition = function(x){
  colnames(x@description$SPILL)
})

#' Extract fluorescent parameters from a flowSet
#' 
#' @param x object \code{flowSet}
#' 
#' @return vector of fluorescent channels
#' 
#' @export
setMethod(getChannels, signature = "flowSet", definition = function(x){
colnames(x[[1]]@description$SPILL)
})

#' Extract fluorescent parameters from a GatingSet
#' 
#' @param x object \code{GatingSet}
#' 
#' @return vector of fluorescent channels
#' 
#' @export
setMethod(getChannels, signature = "GatingSet", definition = function(x){
  colnames(x@data[[1]]@description$SPILL)
})

#' Select flourescent channels for compensation controls
#' 
#' @param x object of class \code{flowSet} or \code{GatingSet} containing compensation controls
#' 
#' @return vector of channels in order of \code{flowSet} or \code{GatingSet} samples.
#' 
#' @export
selectChannels <- function(fs){
  
  opts <- c(getChannels(fs), "Unstained")
  
  # Print sample name and select channel
  channels <- opts[sapply(pData(fs)$name, function(x){
    message("Select a fluorescent channel for the following compensation control:")
    print(x)
    menu(choices = opts, graphics = TRUE)
  })]
  
  return(channels)
}

#' Sample a flowFrame
#' 
#' @param fr object of class \code{flowFrame}.
#' @param size numeric indicating the number of events to keep. If \code{size > nrow(fr)}, size is set to \code{nrow(fr)}.
#' 
#' @return \code{flowFrame} restricted to \code{size} events.
#' 
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter
#' @importFrom flowCore Subset
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

#' Select flowFrames based on pData
#' 
#' @param fs an object of class \code{flowSet}.
#' @param pData vector of form \code{c("column","row")} indicating the rows of \code{pData(fs)} to extract. 
#' 
#' @return an object of calss \code{flowSet} containing the extracted \code{flowFrame} objects.
#' 
#' @export
selectFrames <- function(fs, pData){
  
  # Extract pData from flowSet 
  pd <- pData(fs)
  
  # Find which rows to extract
  rows <- pd[[pData[1]]] == pData[2]
  
  # Get sampleNames of these rows
  sn <- rownames(pd)[rows]
  
  # Extract these samples by name from flowSet
  fs <- fs[sn]
  
  return(fs)
  
}