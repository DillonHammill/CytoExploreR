#' drawGate Root
#' 
#' @param x object of class \code{flowFrame}, \code{flowSet} or \code{GatingSet}.
#' 
#' @export
setGeneric(name="drawGate",
           def=function(x, ...){standardGeneric("drawGate")}
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
#' ellipse, threshold, boundary, interval and quadrant} which can be abbreviated as upper or lower case first letters as well. Default \code{gate_type}
#' is \code{"polygon"}.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted. \code{subSample} is set to 250 000 events by default.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to \code{density()} for 1D plots (defaults to 1.5).
#' @param ... additional arguments for \code{plotDens}.
#' 
#' @importFrom BiocGenerics colnames
#' 
#' @export
setMethod(drawGate, signature = "flowFrame", definition = function(x, channels, alias = NULL, subSample = 250000, gate_type = "polygon", axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  fr <- x
  
  # Check alias is upplied correctly
  checkAlias(alias, gate_type)
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels)
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type)

  # Make one call to drawPlot
  plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
  
  # Construct gates save as filters object
  if(gate_type == "quadrant"){
    
  plot <- TRUE
  gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
    
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
#' ellipse, threshold, boundary, interval and quadrant} which can be abbreviated as upper or lower case first letters as well.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted. \code{subSample} is set to 250 000 events by default.
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
#' @export
setMethod(drawGate, signature = "flowSet", definition = function(x, pData = NULL,channels, alias = NULL, subSample = 250000, gate_type = "polygon", axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){

  fs <- x
  
  # Restrict to samples matching pData requirements
  if(!is.null(pData)){
    
    # Extract samples using selectFrames
    fs <- selectFrames(fs, pData)
    
  }
  
  fr <- as(fs, "flowFrame")
  
  # Check alias is upplied correctly
  checkAlias(alias, gate_type)
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels)
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type)
  
  # Make one call to drawPlot
  plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
  
  # Construct gates save as filters object
  if(gate_type == "quadrant"){
    
    plot <- TRUE
    gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
    
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
#' ellipse, threshold, boundary, interval and quadrant} which can be abbreviated as upper or lower case first letters as well.
#' @param template the name of an existing R object in global environment constructed using \code{openCyto::add_pop()} (e.g. "template").
#' @param gtfile name of \code{gatingTemplate} csv file if \code{template = TRUE}.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted. \code{subSample} is set to 250 000 events by default.
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
#' 
#' @export
setMethod(drawGate, signature = "GatingSet", definition = function(x, pData = NULL, parent = "root", alias = NULL, channels, gate_type = "polygon", subSample = 250000, template = NULL, gtfile = "gatingTemplate.csv", axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){

  gs <- x
  fs <- getData(x, parent)
  
  # Restrict to samples matching pData requirements
  if(!is.null(pData)){
    
    # Extract samples using selectFrames
    fs <- selectFrames(fs, pData)

  }
  
  fr <- as(fs, "flowFrame")
  
  # Check alias is upplied correctly
  checkAlias(alias, gate_type)
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels)
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type)
  
  # Make one call to drawPlot
  plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
  
  # Construct gates save as filters object
  if(gate_type == "quadrant"){
    
    plot <- TRUE
    gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
    
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
  
    # Prepare arguments for gatingTemplate Entry
    if(length(channels) == 2){
      channels <- paste(channels, collapse = ",")
    }
    
    if(length(alias) > 1){
      pop <- "*"
    }else{
      pop = "+"
    }  
    
    if(length(alias) > 1){
      alias <- paste(alias, collapse = ",")
    }
    
    # Check if gatingTemplate exists - if not create one
    if (!is.null(template)) {
      
      Template <- get(template, envir = .GlobalEnv)
      
      Template <- rbind(Template, add_pop(
        gs, alias = alias, parent = parent, pop = pop, dims = channels, gating_method = "manualGate",
        gating_args = list(gate = gates)
      ))
      assign(template, Template, envir = .GlobalEnv)
      
    } else {
      
      Template <- add_pop(
        gs, alias = alias, parent = parent, pop = pop, dims = channels, gating_method = "manualGate",
        gating_args = list(gate = gates)
      )
      Template <<- Template
      
    }
    
    write.csv(Template, gtfile, row.names = FALSE)
  
  return(gs)
  
})
