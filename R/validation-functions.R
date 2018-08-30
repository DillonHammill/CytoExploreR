#' Check Channels Supplied to drawGate.
#' 
#' \code{checkChannels} will check whether the supplied channels are valid for the \code{flowSet} or \code{GatingSet} and that only 1or 2 channels are supplied.
#' 
#' @param x an object of class \code{flowFrame}, \code{flowSet} or \code{GatingSet}.
#' @param channels vector of channel names (e.g. c("PE-A","APC-A")).
#'
#' @return Stop gating process if the supplied channels do ot meet the requirements for drawGate.
#' 
#' @importFrom BiocGenerics colnames
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
checkChannels <- function(x, channels){
  
  chans <- colnames(x)
  
  if(!length(channels) %in% c(1,2)){
    
    stop("Incorrect number of channels supplied. Please supply 1-2 channel names for plotting and gating.")
  
  }
  
  if(length(channels) == 1){
    
    if(channels %in% chans == FALSE){
      
      stop(paste(channels[channels %in% chans == FALSE]),"is not a valid channel for this", class(x),"!")
      
    }else{
      
      
    }
  
  }else if(length(channels) == 2){
    
    if(sum(channels %in% chans) == 0){
      
      stop(paste(paste(channels[channels %in% chans == FALSE], collapse = " & "), "are not valid channels for this", class(x)),"!")
      
    }else if(sum(channels %in% chans) == 1){
      
      stop(paste(channels[channels %in% chans == FALSE],"is not a valid channel for this", class(x)),"!")
      
    }else{
      
      
    }
    
    
  }
  
}

#' Check Gate Types(s) Supplied drawGate.
#' 
#' @param gate_type vector indicating the types of gates to construct using \code{drawGate}.
#' @param alias names of the populations to be gated.
#' 
#' @return Stop gating process if gate_type is incorrect or returns \code{gate_type} as full lower case name(s). If a single gate_type is supplied for multiple populations,
#' the same gate_type will be used for all populations.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
checkGateType <- function(gate_type, alias){
  
  gts <- c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q", "web", "Web", "w","W")
  
  if(!all(gate_type %in% gts)){
    
  if(length(gate_type[gate_type %in% gts == FALSE]) >= 2){
      
    stop(paste(paste(gate_type[gate_type %in% gts == FALSE], collapse = " & "), "are not valid gate_types for drawGate!"))
      
  }else{
      
    stop(paste(gate_type[gate_type %in% gts == FALSE],"is not a valid gate_type for drawGate!"))
      
  }
    
  }
  
  gate_type[gate_type %in% c("polygon", "Polygon", "p", "P")] <- "polygon"
  
  gate_type[gate_type %in% c("rectangle", "Rectangle", "r", "R")] <- "rectangle"
  
  gate_type[gate_type %in% c("interval", "Interval", "i", "I")] <- "interval"
  
  gate_type[gate_type %in% c("threshold", "Threshold", "t", "T")] <- "threshold"
  
  gate_type[gate_type %in% c("boundary", "Boundary", "b", "B")] <- "boundary"
  
  gate_type[gate_type %in% c("ellipse", "Ellipse", "e", "E")] <- "ellipse"
  
  gate_type[gate_type %in% c("quadrant", "Quadrant", "q", "Q")] <- "quadrant"
  
  gate_type[gate_type %in% c("web", "Web", "w", "W")] <- "web"
  
  # Repeat gate_type to equal length of alias
  if(length(gate_type) != length(alias) & gate_type != "quadrant" & gate_type != "web"){
  
      gate_type <- rep(gate_type, length(alias))
      
  }
  
  return(gate_type)
}

#' Check Supplied Alias
#' 
#' @param alias vector indicating the names of the populations to be gated.
#' @param gate_type vector indicating the type(s) of gate(s) to be constructed.
#' 
#' @return Stops the gating process if alias is missing or of the incorrect length given the gate_type.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
checkAlias <- function(alias, gate_type){
  
  if(is.null(alias)){
    
    stop("The name(s) of the population(s) to be gated must be supplied as the alias argument.")
  
  }

  if(gate_type == "quadrant" & length(alias) != 4){
    
    stop("Supply 4 population names to alias argument to construct quadrant gates.")
    
  }
  
  if(length(gate_type) > 1){
    
    if(length(alias) != length(gate_type)){
      
      stop("Length of alias must be the same length as gate_type for multi-gates.")
      
    }
  }
}

#' Check Operating System & Open New Graphics Device
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
checkOSGD <- function(){
  
  if(.Platform$OS.type == "windows"){
    
    windows()
    
  }else if(.Platform$OS.type == "unix"){
    
    if(Sys.info()["sysname"] == "Linux"){
      
      X11()
      
    }else if(Sys.info()["sysname"] == "Darwin"){
      
      quartz()
      
    }
  
  }
  
}

#' Check File Exists in Working Directory
#' 
#' @param name filename including file extension to be checked.
#' 
#' @return TRUE/FALSE if file exists in the current working directory.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
checkFile <- function(name){
  
  if(length(which(list.files() == name)) != 0){
    
    # File exists in working directory
    return(TRUE)
  
  }else if(length(which(list.files() == name)) == 0){
    
    # File does not exist in working directory
    return(FALSE)
  }
}

#' Check gatingTemplate for Existing Entry
#' 
#' @param parent name of the parent population.
#' @param alias name of the population of interest.
#' @param gtfile csv file name of the gatingTemplate.
#' 
#' @return Null if entry does not exist or N/Y based on user input.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
checkTemplate <- function(parent, alias, gtfile){
  
  gt <- read.csv(gtfile, header = TRUE)
  
  # Parent and alias entries match file
  if(any(gt$parent %in% parent & gt$alias %in% alias)){
      
      gt <- gt[gt$parent %in% parent & gt$alias %in% alias,]
    
      message(paste(paste(gt$alias, collapse = " & "),"already exists in",gtfile,"."))
      rsp <- readline("Do you want to override the existing gate?(Y/N)")
    
      if(rsp %in% c("y","Y","YES","Yes","yes")){
      
         rsp <- "Y"
    
      }else if(rsp %in% c("n","N","NO","No","no")){
      
         rsp <- "N"
      
      }
      
    }else{
    
    rsp <- NULL
    
  }
  
  return(rsp)
}