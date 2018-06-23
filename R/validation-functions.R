#' Check Supplied Channels.
#' 
#' @param x an object of class \code{flowFrame}, \code{flowSet} or \code{GatingSet}.
#' @param channels vector of length 1 or 2 indicating the names of the channels to use for plotting and gating.
#'
#' @return Stops gating process if channels are incorrect or returns nothing.
checkChannels <- function(x, channels){
  
  chans <- colnames(x)
  
  if(!length(channels) %in% c(1,2)){
    
    stop("Incorrect number of channels supplied. Please supply 1-2 channel names for plotting and gating.")
  
  }else if(!all(channels %in% chans)){
    
    if(length(channels[channels %in% chans == FALSE]) == 2){
    
    stop(paste(paste(channels[channels %in% chans == FALSE], collapse = " & "), "are not valid channels for this", class(x)),"!")
    
    }else{
      
      stop(paste(channels[channels %in% chans == FALSE],"is not a valid channel for this", class(x)),"!")
      
    }
    
  }
  
}

#' Check Supplied Gate Types(s).
#' 
#' @param gate_type vector indicating the types of gates to construct using \code{drawGate}.
#' 
#' @return Stops gating process if gate_type is incorrect or returns nothing.
checkGateType <- function(gate_type){
  
  gts <- c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q")
  
  if(!all(gate_type %in% gts)){
    
  if(length(gate_type[gate_type %in% gts == FALSE]) >= 2){
      
    stop(paste(paste(gate_type[gate_type %in% gts == FALSE], collapse = " & "), "are not valid gate_types for drawGate!"))
      
  }else{
      
    stop(paste(gate_type[gate_type %in% gts == FALSE],"is not a valid gate_type for drawGate!"))
      
  }
    
  }
  
}
