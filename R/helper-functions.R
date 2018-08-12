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

#' Remove gate and edit gatingTemplate .csv file
#' 
#' @param gs an object of class \code{GatingSet}.
#' @param alias name(s) of the gate to remove (e.g. "Single Cells").
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g. "gatingTemplate.csv").
#' @param children logical indicating whether descendant populations should also be removed from the gtfile, set to \code{TRUE} by default.
#' 
#' @return an object of class \code{gatingSet} with gate and children removed, as well as gatingTemplate file with population removed.
#' 
#' @export
removeGate <- function(gs, alias = NULL, gtfile = NULL, children = TRUE){
  
  # Supply alias
  if(is.null(alias)){
    
    stop("Please supply the name of the population to be removed.")
  
  }
  
  # Supply gtfile
  if(is.null(gtfile)){
    
    stop("Please supply the name of the gatingTemplate csv file to remove the gate.")
  
  }
  
  # Read in gatingTemplate
  gt <- read.csv(gtfile, header = TRUE)
  
  # For each alias remove from GatingSet and gatingTemplate
  lapply(alias, function(x){
    
    # Remove node & children from GatingSet
    Rm(x, gs)
    
    # Remove row entry from gatingTemplate
    if(children == TRUE){
    
    # Remove rows with alias or parent = alias
    if(length(which(gt$alias == x)) == 0){
      
      stop("Supplied alias is not a valid name for a gated population.")
      
    }else if(length(which(gt$parent == x)) != 0){
      
    indx <- c(which(gt$alias == x), which(gt$parent == x))
    gt <- gt[-indx,]  
    write.csv(gt, gtfile, row.names = FALSE)
      
    }else if(length(which(gt$parent == x)) == 0){
      
      indx <- which(gt$alias == x)
      gt <- gt[-indx,]  
      write.csv(gt, gtfile, row.names = FALSE)
      
    }
    
    }else if(children == FALSE){
      
      indx <- which(gt$alias == x)
      gt <- gt[-indx,]  
      write.csv(gt, gtfile, row.names = FALSE)
      
    }
  })
  
  return(gs)
  
}

#' Extract saved gate(s) from gatingTemplate.
#' 
#' @param parent name of the parental population.
#' @param alias name of the population for which the gate must be extracted.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g. "gatingTemplate.csv") where the gate(s) are saved.
#' 
#' @export
extractGate <- function(parent, alias, gtfile){
  
  # Load gtfile into gatingTemplate
  gt <- gatingTemplate(gtfile)
  
  # Extract population nodes from gt
  pops <- getNodes(gt, only.names = TRUE)
  
  # Extract gates given parent and child node(s)
  gates <- lapply(alias, function(x){
    
    if(parent == "root"){
    
    parent <- "root"
    alias <- paste("/",x,sep="")
  
  } else {
    
    parent <- paste("/",parent,"")
    alias <- paste("/",parent,"/",x,sep="")
  
  }  
  
  gate <- getGate(gt, parent, alias)
  gate <- eval(parameters(gate)$gate)
    
  })
  
  gates <- filters(gates)
  return(gates)
}

#' Edit existing gate(s).
#' 
#' @param x an object of class \code{GatingSet}.
#' @param pData vector of form \code{c("column","row")} indicating the rows of \code{pData(fs)} to extract. 
#' @param parent name of the parental population.
#' @param alias name(s) of the gate to edit (e.g. "Single Cells").
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied with
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval, quadrant and web} which can be abbreviated as upper or lower case first letters as well. Default \code{gate_type}
#' is \code{"polygon"}.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g. "gatingTemplate.csv") where the gate is saved.
#' 
#' @return an object of calss \code{GatingSet} with edited gate applied, as well as gatingTemplate file with editied gate saved.
#' 
#' @export
editGate <- function(x, pData = NULL, parent = NULL, alias = NULL, gate_type = NULL, gtfile = NULL){
  
  # Rename x to gs
  gs <- x
  
  # Restrict to samples matching pData requirements
  if(!is.null(pData)){
    
    # Extract samples using selectFrames
    fs <- selectFrames(fs, pData)
    
  }
  
  fr <- as(fs, "flowFrame")
  
  # Check alias is supplied correctly
  checkAlias(alias, gate_type)
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type)
  
  # Extract gate(s) from gtfile gating_args
  gates <- extractGate(gtfile = gtfile, parent = parent, alias = alias)
  
  # Extract channels from gates
  channels <- parameters(gates[[1]])
  
  # Plot data for visualisation
  drawPlot(fr = fr, channels = channels)
  
  # Plot existing gates
  plotGates(gates, col = "magenta")
  
  # Make new call to drawGate to get new gates - set plot = FALSE
  new <- drawGate(fr, alias = alias, channels = channels, plot = FALSE, gate_type = gate_type)

  # Find and Edit gatingTemplate entries - each alias and gate separate
  gt <-read.csv(gtfile, header = TRUE)
  gt <- as.data.table(gt)
  
  for(i in 1:length(alias)){
    
    gt[gt$alias == alias[i] & gt$parent == parent, "gating_args"] <- .argDeparser(list(gate = new[[i]]))
    
  }
  
  gt <- data.frame(gt)
  write.csv(gt, gtfile, row.names = FALSE)
  
  # Apply entire gatingTemplate to GatingSet
  gt <- gatingTemplate(gtfile)
  rt <- getData(gs, "root")   # Remove all nodes - extract root
  gs <- GatingSet(rt)         # Add root to GatingSet
  
  gating(gt,gs)
  
  assign(deparse(substitute(x)), gs, envir=globalenv())
    
}