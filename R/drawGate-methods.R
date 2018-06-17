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
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied wuth
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval and quadrant} which can be abbreviated as upper or lower case first letters as well.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to \code{density()} for 1D plots (defaults to 1.5).
#' @param ... additional arguments for \code{plotDens}.
#' 
#' @export
setMethod(drawGate, signature = "flowFrame", definition = function(x, channels, alias, subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  fr <- x
  
  # Supported gate types
  gate_types <- c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q")
  
  # Check that gate_type has been supplied or default to polygon type
  if(missing(gate_type)){
    
    message("No gate type supplied - gate type set to polygon.")
    gate_type <- "polygon"
    
  }else if(anyNA(match(gate_type, gate_types))){
    
    message("Invalid gate type supplied - gate type set to polygon. Supported gate types include polygon, rectangle, interval, threshold, boundary, ellipse and quadrant")
    gate_type[which(is.na(match(gate_type,gate_types)))] <- "polygon"
    
  }
  
  if(length(alias) > 1){
    plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
  }
  
  gates <- mapply(function(gate_type, alias, plot){
    
    if(gate_type %in% c("polygon", "Polygon", "p", "P")){
      
      drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("rectangle", "Rectangle", "r", "R")){
      
      drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("interval", "Interval", "i", "I")){
      
      drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
      
    }else if(gate_type %in% c("threshold", "Threshold", "t", "T")){
      
      drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("boundary", "Boundary", "b", "B")){ 
      
      drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("ellipse", "Ellipse", "e", "E")){
      
      drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("quadrant", "Quadrant", "q", "Q")){
      
      drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }
    
  }, gate_type, alias, plot = plot)
  
  gates <- filters(gates)
  return(gates)
})

#' drawGate flowSet Method
#' 
#' @param x object of class \code{flowSet}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied wuth
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval and quadrant} which can be abbreviated as upper or lower case first letters as well.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to \code{density()} for 1D plots (defaults to 1.5).
#' @param ... additional arguments for \code{plotDens}.
#' 
#' @return list of gates stored in object of class \code{filters}.
#' 
#' @export
setMethod(drawGate, signature = "flowSet", definition = function(x, channels, alias, subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){

  fr <- as(x, "flowFrame")
  
  # Supported gate types
  gate_types <- c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q")
  
  # Check that gate_type has been supplied or default to polygon type
  if(missing(gate_type)){
    
    message("No gate type supplied - gate type set to polygon.")
    gate_type <- "polygon"
    
  }else if(anyNA(match(gate_type, gate_types))){
    
    message("Invalid gate type supplied - gate type set to polygon. Supported gate types include polygon, rectangle, interval, threshold, boundary, ellipse and quadrant")
    gate_type[which(is.na(match(gate_type,gate_types)))] <- "polygon"
    
  }
  
  if(length(alias) > 1){
    plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
  }
  
  gates <- mapply(function(gate_type, alias, plot){
    
    if(gate_type %in% c("polygon", "Polygon", "p", "P")){
      
      drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("rectangle", "Rectangle", "r", "R")){
      
      drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("interval", "Interval", "i", "I")){
      
      drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
      
    }else if(gate_type %in% c("threshold", "Threshold", "t", "T")){
      
      drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("boundary", "Boundary", "b", "B")){ 
      
      drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("ellipse", "Ellipse", "e", "E")){
      
      drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("quadrant", "Quadrant", "q", "Q")){
      
      drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }
    
  }, gate_type, alias, plot = plot)
  
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
#' @param file name of \code{gatingTemplate} csv file if \code{template = TRUE}.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
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
#' @export
setMethod(drawGate, signature = "GatingSet", definition = function(x, pData = NULL, parent = "root", alias, channels, gate_type = NULL, subSample = NULL, template = NULL, file = "gatingTemplate.csv", axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){

  gs <- x
  fs <- getData(x, parent)
  
  # Restrict to samples matching pData requirements
  if(is.null(pData)){
    
    fr <- as(fs, "flowFrame")
    
  }else{
    
    idx <- pData(fs)[[pData[1]]] == pData[2]
    sn <- rownames(pData(fs))[idx]
    fr <- as(fs[sn], "flowFrame")
    
  }
  
  fr <- as(fs, "flowFrame")
  
  # Supported gate types
  gate_types <- c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q")
  
  # Check that gate_type has been supplied or default to polygon type
  if(missing(gate_type)){
    
    message("No gate type supplied - gate type set to polygon.")
    gate_type <- "polygon"
    
  }else if(anyNA(match(gate_type, gate_types))){
    
    message("Invalid gate type supplied - gate type set to polygon. Supported gate types include polygon, rectangle, interval, threshold, boundary, ellipse and quadrant")
    gate_type[which(is.na(match(gate_type,gate_types)))] <- "polygon"
    
  }
  
  if(length(alias) > 1){
    plot <- c(TRUE, rep(FALSE, (length(alias) - 1)))
  }
  
  gates <- mapply(function(gate_type, alias, plot){
    
    if(gate_type %in% c("polygon", "Polygon", "p", "P")){
      
      drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("rectangle", "Rectangle", "r", "R")){
      
      drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("interval", "Interval", "i", "I")){
      
      drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
      
    }else if(gate_type %in% c("threshold", "Threshold", "t", "T")){
      
      drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("boundary", "Boundary", "b", "B")){ 
      
      drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("ellipse", "Ellipse", "e", "E")){
      
      drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("quadrant", "Quadrant", "q", "Q")){
      
      drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }
    
  }, gate_type, alias, plot = plot)
  
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
    
    write.csv(Template, file, row.names = FALSE)
  
  return(gs)
  
})
