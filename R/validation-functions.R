#' Check Supplied Channels
#' 
#' @param x an object of class \code{flowFrame}, \code{flowSet} or \code{GatingSet}.
#' @param channels vector of length 1 or 2 indicating the names of the channels to use for plotting and gating.
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
    
  }else{
    
    return(channels)
  
  }
  
}
