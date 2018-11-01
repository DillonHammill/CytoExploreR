#' drawGate
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' \code{drawGate} is a convenient wrapper for the gating functions shipped with
#' \code{cytoRSuite} to facilitate analysis of flow cytometry by gate drawing.
#' Using \code{drawGate} users can specify the type of gate(s) to be constructed
#' through the \code{type} argument and \code{drawGate} will automatically
#' handle plotting the data and make calls to the relevant gating function(s) to
#' construct the gates around populations of interest. \code{drawGate} has
#' methods for \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}} and
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} objects, refer to
#' their respective help pages for more information. The flowFrame and flowSet
#' methods simply return the constructed gates as a list of
#' \code{\link[flowCore:filters-class]{filters}}, whilst the GatingSet method
#' automatically applies the constructed gates to the GatingSet and saves the
#' constructed gates in an \code{openCyto}
#' \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}for future use.
#' See \code{\link{editGate}} and \code{\link{removeGate}} to manipulate
#' constructed gates and modify their entries in the gatingTemplate.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}, \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{drawGate,flowFrame-method}}
#' @seealso \code{\link{drawGate,flowSet-method}}
#' @seealso \code{\link{drawGate,GatingSet-method}}
#'
#' @export
setGeneric(name = "drawGate",
           def = function(x, ...){standardGeneric("drawGate")}
)

#' drawGate flowFrame Method.
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types include \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param subSample  numeric indicating the number of events to plot, set to all
#'   events by default. Reducing the sample size can significantly increase
#'   plotting speed on less powerful machines.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to
#'   \code{\link[stats:density]{density}} for 1-D plots (defaults to 1.5).
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a \code{\link[flowCore:filters-class]{filters}} list containing the
#'   drawn gate objects.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore filters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate,flowSet-method}}
#' @seealso \code{\link{drawGate,GatingSet-method}}
#'
#' @export
setMethod(drawGate, signature = "flowFrame", definition = function(x, channels = NULL, alias = NULL, subSample = 250000, type = "polygon", axis = "x", adjust = 1.5, labels = TRUE, plot = TRUE, ...){
  
  # Assign x to fr
  fr <- x
  
  # Check type argument is valid
  type <- checkGateType(type = type, alias = alias)
  
  # Set default type for 1D gates to interval
  if(length(channels) == 1 & all(type %in% "polygon")){
    
    type <- rep("interval", length(type))
    
  }
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, type = type) 
  
  # Check supplied channel(s) are valid
  channels <- checkChannels(fr, channels = channels, plot = TRUE)

  # Make one call to drawPlot
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, subSample = subSample, popup = TRUE, legend = FALSE, ...)
    
  }
  
  # Construct gates save as filters object
  if(length(type) == 1 & type[1] == "quadrant"){
    
  gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
    
  }else if(length(type) == 1 & type[1] == "web"){

  gates <- drawWeb(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
  
  }else{
  
  gates <- mapply(function(type, alias){
    
    if(type == "polygon"){
      
      drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
      
    }else if(type == "rectangle"){
      
      drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
      
    }else if(type == "interval"){
      
      drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, axis = axis, labels = labels,...)
      
    }else if(type == "threshold"){
      
      drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
      
    }else if(type == "boundary"){ 
      
      drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
      
    }else if(type == "ellipse"){
      
      drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
      
    }
    
  }, type, alias)
  
  }
  
  gates <- filters(gates)
  
  return(gates)
})

#' drawGate flowSet Method
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param select vector containing the indicies of samples within gs to use for
#'   plotting.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param subSample  numeric indicating the number of events to plot, set to all
#'   events by default. Reducing the sample size can significantly increase
#'   plotting speed on less powerful machines.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to
#'   \code{\link[stats:density]{density}} for 1-D plots (defaults to 1.5).
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowSet-method}}.
#'
#' @return a \code{\link[flowCore:filters-class]{filters}} list containing the
#'   drawn gate objects.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore filters
#' @importFrom methods as
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{plotCyto1d,flowSet-method}}
#' @seealso \code{\link{plotCyto2d,flowSet-method}}
#' @seealso \code{\link{drawGate,flowFrame-method}}
#' @seealso \code{\link{drawGate,GatingSet-method}}
#'
#' @export
setMethod(drawGate, signature = "flowSet", definition = function(x, select = NULL, channels = NULL, alias = NULL, subSample = 250000, type = "polygon", axis = "x", adjust = 1.5, labels = TRUE, plot = TRUE, ...){

  # Assign x to fs
  fs <- x
  
  # Restrict to samples matching pData requirements
  if(!is.null(select)){
    
    if(class(select) != "numeric"){
      
      stop("Vector supplied to select argument should contain the numeric indicies of the samples to select.")
      
    }
    
    # Extract samples using selectFrames
    fs <- fs[select]
    
  }
  fr <- as(fs, "flowFrame")
  
  # Check type argument is valid
  type <- checkGateType(type = type, alias = alias)
  
  # Set default type for 1D gates to interval
  if(length(channels) == 1 & all(type %in% "polygon")){
    
    type <- rep("interval", length(type))
    
  }
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, type = type) 
  
  # Check supplied channel(s) are valid
  channels <- checkChannels(fr, channels = channels, plot = TRUE)

  # Make one call to drawPlot
  if(plot == TRUE){
    
    plotCyto(fs, channels = channels, subSample = subSample, popup = TRUE, legend = FALSE, merge = TRUE, ...)
    
  }
  
  # Construct gates save as filters object
  if(length(type) == 1 & type[1] == "quadrant"){
    
    gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
    
  }else if(length(type) == 1 & type[1] == "web"){
    
    gates <- drawWeb(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
    
  }else{
    
    gates <- mapply(function(type, alias){
      
      if(type == "polygon"){
        
        drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "rectangle"){
        
        drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "interval"){
        
        drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, axis = axis, labels = labels,...)
        
      }else if(type == "threshold"){
        
        drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "boundary"){ 
        
        drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "ellipse"){
        
        drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }
      
    }, type, alias)
    
  }
  
  gates <- filters(gates)
  return(gates)
  
})

#' drawGate GatingSet Method
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select vector containing the indicies of samples within gs to use for
#'   plotting.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param parent name of the \code{parent} population to extract for gating.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param gtfile name of \code{gatingTemplate} csv file to be saved.
#' @param subSample  numeric indicating the number of events to plot, set to all
#'   events by default. Reducing the sample size can significantly increase
#'   plotting speed on less powerful machines.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param adjust smooothing factor passed to
#'   \code{\link[stats:density]{density}} for 1-D plots (defaults to 1.5).
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,GatingSet-method}}.
#'
#' @return drawn gates are applied to the
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and saved to a
#'   \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}.
#'   
#' @importFrom BiocGenerics colnames
#' @importFrom flowWorkspace getData
#' @importFrom openCyto add_pop
#' @importFrom methods as
#' @importFrom utils read.csv write.csv
#' @importFrom flowCore filters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @seealso \code{\link{plotCyto,GatingSet-method}}
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate,flowFrame-method}}
#' @seealso \code{\link{drawGate,flowSet-method}}
#'
#' @export
setMethod(drawGate, signature = "GatingSet", definition = function(x, select = NULL, parent = "root", alias = NULL, channels = NULL, type = "polygon", subSample = 250000, gtfile = NULL, axis = "x", adjust = 1.5, labels = TRUE, plot = TRUE, ...){

  # Assign x to gs
  gs <- x
  
  # Check whether a gatingTemplate ready exists for this population
  if(!is.null(gtfile)){
    
    # Check whether gate already exists in gtfile
    checkTemplate(parent, alias, gtfile)
  
  }
  
  fs <- flowWorkspace::getData(x, parent)
  
  # Restrict to samples matching pData requirements
  if(!is.null(select)){
    
    if(class(select) != "numeric"){
      
      stop("Vector supplied to select argument should contain the numeric indicies of the samples to select.")
      
    }
    
    # Extract samples using selectFrames
    fs <- fs[select]

  }
  
  fr <- as(fs, "flowFrame")
  
  # Check type argument is valid
  type <- checkGateType(type = type, alias = alias)
  
  # Set default type for 1D gates to interval
  if(length(channels) == 1 & all(type %in% "polygon")){
    
    type <- rep("interval", length(type))
    
  }
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, type = type) 
  
  # Check supplied channel(s) are valid
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Make one call to drawPlot
  if(plot == TRUE){
    
    plotCyto(x = gs, parent = parent, channels = channels, subSample = subSample, popup = TRUE, legend = FALSE, merge = TRUE, ...)
    
  }
  
  # Construct gates save as filters object
  if(length(type) == 1 & type[1] == "quadrant"){
    
    gates <- drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
    
  }else if(length(type) == 1 & type[1] == "web"){
    
    gates <- drawWeb(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
    
  }else{
    
    gates <- mapply(function(type, alias){
      
      if(type == "polygon"){
        
        drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "rectangle"){
        
        drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "interval"){
        
        drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, axis = axis, labels = labels,...)
        
      }else if(type == "threshold"){
        
        drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "boundary"){ 
        
        drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }else if(type == "ellipse"){
        
        drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = FALSE, labels = labels,...)
        
      }
      
    }, type, alias)
    
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
        gs = x, alias = alias[i], parent = parent, pop = pop, dims = paste(channels, collapse = ","), gating_method = "manualGate",
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
        gs = x, alias = alias[i], parent = parent, pop = pop, dims = paste(channels, collapse = ","), gating_method = "manualGate",
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
        gs = x, alias = alias[i], parent = parent, pop = pop, dims = paste(channels, collapse = ","), gating_method = "manualGate",
        gating_args = list(gate = gates[[i]])
      )
    
    }
    pops <- do.call("rbind", pops)
    gt <- rbind(gt, pops)
    
    write.csv(gt, gtfile, row.names = FALSE)
    
  }
  
  return(pops)
  
})
