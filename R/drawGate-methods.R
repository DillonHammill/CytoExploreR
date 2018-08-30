#' drawGate
#' 
#' \code{drawGate} is a convenient wrapper for the gating functions shipped with \code{cytoSuite} to facilitate analysis of flow cytometry by gate drawing. Using \code{drawGate} 
#' users can specify the type of gate(s) to be constructed through the \code{gate_type} argument and \code{drawGate} will automatically handle plotting the data and make calls
#' to the relevant gating function(s) to construct the gates around populations of interest. \code{drawGate} has methods for \code{flowFrame}, \code{flowSet} and \code{GatingSet}
#' objects, refer to their respective help pages for more information - ?`drawGate,flowFrame-method`, ?`drawGate,flowSet-method` or ?`drawGate,GatingSet-method`. The flowFrame and
#' flowSet methods simply return the constructed gates as a list of \code{filters}, whilst the GatingSet method automatically applies the constructed gates to the GatingSet and saves
#' the constructed gates in an \code{openCyto} gatingTemplate for future use. See ?editGate and ?removeGate to manipulate constructed gates and modify their entries in the gatingTemplate.
#' 
#' @param x object of class \code{flowFrame}, \code{flowSet} or \code{GatingSet}.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
setGeneric(name = "drawGate",
           def = function(x, ...){standardGeneric("drawGate")}
)

#' drawGate flowFrame Method.
#' 
#' @param x object of class \code{flowFrame}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied with
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval, quadrant and web} which can be abbreviated as upper or lower case first letters as well. Default \code{gate_type}
#' is \code{"polygon"}.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted. \code{subSample} is set to 250000 events by default.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to \code{density()} for 1D plots (defaults to 1.5).
#' @param ... additional arguments for \code{plotDens}.
#' 
#' @importFrom BiocGenerics colnames
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
setMethod(drawGate, signature = "flowFrame", definition = function(x, channels = NULL, alias = NULL, subSample = 250000, gate_type = "polygon", axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  fr <- x
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type, alias = alias)
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, gate_type = gate_type) 
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels = channels)

  # Make one call to drawPlot
  if(plot == TRUE){
    
    plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
   
  }else if (plot == FALSE){
    
    plot <- c(FALSE, rep(FALSE, (length(alias) - 1)))
    
  }
  
  # Construct gates save as filters object
  if(length(gate_type) == 1 & gate_type[1] == "quadrant"){
    
  gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot[1], labs = labs,...)
    
  }else if(length(gate_type) == 1 & gate_type[1] == "web"){

  gates <- drawWeb(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot[1], labs = labs,...)
  
  }else{
  
  gates <- mapply(function(gate_type, alias, plot){
    
    if(gate_type == "polygon"){
      
      drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type == "rectangle"){
      
      drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type == "interval"){
      
      drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
      
    }else if(gate_type == "threshold"){
      
      drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type == "boundary"){ 
      
      drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type == "ellipse"){
      
      drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }
    
  }, gate_type, alias, plot = plot)
  
  }
  
  gates <- filters(gates)
  
  return(gates)
})

#' drawGate flowSet Method
#' 
#' @param x object of class \code{flowSet}.
#' @param pData vector of indicating \code{c(column,row)} of \code{pData(fs)} used extract particular samples for gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied wuth
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval, quadrant and web} which can be abbreviated as upper or lower case first letters as well.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted. \code{subSample} is set to 250000 events by default.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to \code{density()} for 1D plots (defaults to 1.5).
#' @param ... additional arguments for \code{plotDens}.
#' 
#' @return list of gates stored in object of class \code{filters}.
#' 
#' @importFrom BiocGenerics colnames
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
setMethod(drawGate, signature = "flowSet", definition = function(x, pData = NULL, channels = NULL, alias = NULL, subSample = 250000, gate_type = "polygon", axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){

  fs <- x
  
  # Restrict to samples matching pData requirements
  if(!is.null(pData)){
    
    # Extract samples using selectFrames
    fs <- selectFrames(fs, pData)
    
  }
  
  fr <- as(fs, "flowFrame")
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type, alias = alias)
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, gate_type = gate_type) 
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels = channels)
  
  # Make one call to drawPlot
  if(plot == TRUE){
    
    plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
    
  }else if (plot == FALSE){
    
    plot <- c(FALSE, rep(FALSE, (length(alias) - 1)))
    
  }
  
  # Construct gates save as filters object
  if(length(gate_type) == 1 & gate_type[1] == "quadrant"){
    
    gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot[1], labs = labs,...)
    
  }else if(length(gate_type) == 1 & gate_type[1] == "web"){
    
    gates <- drawWeb(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot[1], labs = labs,...)
    
  }else{
    
    gates <- mapply(function(gate_type, alias, plot){
      
      if(gate_type == "polygon"){
        
        drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "rectangle"){
        
        drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "interval"){
        
        drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
        
      }else if(gate_type == "threshold"){
        
        drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "boundary"){ 
        
        drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "ellipse"){
        
        drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }
      
    }, gate_type, alias, plot = plot)
    
  }
  
  gates <- filters(gates)
  return(gates)
  
})

#' drawGate GatingSet Method
#' 
#' @param x object of class \code{GatingSet}.
#' @param pData vector of indicating \code{c(column,row)} of \code{pData(gs)} used extract particular samples for gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param parent name of the \code{parent} population to extract for gating.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied wuth
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval, quadrant and web} which can be abbreviated as upper or lower case first letters as well.
#' @param gtfile name of \code{gatingTemplate} csv file to be saved.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted. \code{subSample} is set to 250000 events by default.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to \code{density()} for 1D plots (defaults to 1.5).
#' @param ... additional arguments for \code{plotDens}.
#'
#' 
#' @return object of class \code{GatingSet} with gates applied. The \code{gatingTemplate} can also be optionally written to .csv file
#' and entries add to the R object called \code{Template} generated use the \code{openCyto::add_pop} API.
#' 
#' @importFrom BiocGenerics colnames 
#' @importFrom flowWorkspace getData
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
setMethod(drawGate, signature = "GatingSet", definition = function(x, pData = NULL, parent = "root", alias = NULL, channels = NULL, gate_type = "polygon", subSample = 250000, gtfile = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){

  # Check whether a gatingTemplate ready exists for this population
  if(!is.null(gtfile)){
    
  # Check whether gate already exists in gtfile
  checkTemplate(parent, alias, gtfile)
  
  }
  
  gs <- x
  fs <- flowWorkspace::getData(gs, parent)
  
  # Restrict to samples matching pData requirements
  if(!is.null(pData)){
    
    # Extract samples using selectFrames
    fs <- selectFrames(fs, pData)

  }
  
  fr <- as(fs, "flowFrame")
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type, alias = alias)
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, gate_type = gate_type) 
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels = channels)
  
  # Make one call to drawPlot
  if(plot == TRUE){
    
    plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
    
  }else if (plot == FALSE){
    
    plot <- c(FALSE, rep(FALSE, (length(alias) - 1)))
    
  }
  
  # Construct gates save as filters object
  if(length(gate_type) == 1 & gate_type[1] == "quadrant"){
    
    gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot[1], labs = labs,...)
    
  }else if(length(gate_type) == 1 & gate_type[1] == "web"){
    
    gates <- drawWeb(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot[1], labs = labs,...)
    
  }else{
    
    gates <- mapply(function(gate_type, alias, plot){
      
      if(gate_type == "polygon"){
        
        drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "rectangle"){
        
        drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "interval"){
        
        drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
        
      }else if(gate_type == "threshold"){
        
        drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "boundary"){ 
        
        drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }else if(gate_type == "ellipse"){
        
        drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
        
      }
      
    }, gate_type, alias, plot = plot)
    
  }
  
  gates <- filters(gates)
  
  # Prepare gatingTemplate entries
  pop <- "+"
  
  # Use add_pop to apply gates to GatingSet and construct gatingTemplate
  if(is.null(gtfile)){
    
    message("No gatingTemplate file name supplied - creating gatingTemplate.csv to store gates.")

    pops <- list()
    for(i in 1:length(alias)){
      
      pops[[i]] <- add_pop(
        gs, alias = alias[i], parent = parent, pop = pop, dims = paste(channels, collapse = ","), gating_method = "manualGate",
        gating_args = list(gate = gates[[i]])
      )
    
    }
    pops <- do.call("rbind", pops)
    
    write.csv(pops, "gatingTemplate.csv", row.names = FALSE)
    
  }else if(checkFile(gtfile) == FALSE){
    
    message(paste("Supplied gtfile does not exist in working directory - writing", paste(gtfile),"."))
    
    pops <- list()
    for(i in 1:length(alias)){
      
      pops[[i]] <- add_pop(
        gs, alias = alias[i], parent = parent, pop = pop, dims = paste(channels, collapse = ","), gating_method = "manualGate",
        gating_args = list(gate = gates[[i]])
      )
    
    }
    pops <- do.call("rbind", pops)
    
    write.csv(pops, gtfile, row.names = FALSE)
    
    
  }else if(checkFile(gtfile) == TRUE){
    
    gt <- read.csv(gtfile, header = TRUE)
    
    pops <- list()
    for(i in 1:length(alias)){
      
      pops[[i]] <- add_pop(
        gs, alias = alias[i], parent = parent, pop = pop, dims = paste(channels, collapse = ","), gating_method = "manualGate",
        gating_args = list(gate = gates[[i]])
      )
    
    }
    pops <- do.call("rbind", pops)
    gt <- rbind(gt, pops)
    
    write.csv(gt, gtfile, row.names = FALSE)
    
  }
  
  return(gs)
  
})
