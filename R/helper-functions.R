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

#' Select flurescent channels for compensation controls
#' 
#' @param fs object of class \code{flowSet} containing compensation controls
#' 
#' @return add selected channels to \code{pData(fs)$channel}
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