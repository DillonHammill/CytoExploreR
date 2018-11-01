#' plotCyto
#'
#' Explore & visualise flow cytometry data.
#'
#' \code{plotCyto} is a versatile wrapper function built on base graphics
#' capable of plotting a variety of existing flow cytometry classes including
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}},
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} and
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}. The type of plot
#' constructed by \code{plotCyto} is determined by the number of channels
#' supplied to the \code{channels} argument. If a single channel is supplied, a
#' 1-D density distribution will be plotted in the specified channel using
#' \code{\link{plotCyto1d,flowFrame-method}}. Supplying two channels to this
#' argument will result in construction of a 2-D scatterplot with a blue-red
#' density colour scale for points using
#' \code{\link{plotCyto2d,flowFrame-method}}.
#'
#' Some key features of \code{plotCyto} include: \itemize{ \item minimal coding
#' \item stacked density distributions \item 2-D contour lines \item back-gating
#' \item overlays \item population labelling \item grid to visualise multiple
#' samples \item plot gate objects \item plot in pop-up windows \item sampling
#' \item fully customisable \item copy and paste arguments from drawGate }
#'
#' For more method-specific information refer to
#' \code{\link{plotCyto1d,flowFrame-method}} and
#' \code{\link{plotCyto2d,flowFrame-method}}.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet} or
#'   \code{GatingSet}.
#' @param ... additional method-specific arguments for plotCyto.
#'
#' @seealso \code{\link{plotCyto,flowFrame-method}}
#' @seealso \code{\link{plotCyto,flowSet-method}}
#' @seealso \code{\link{plotCyto,GatingHierarchy-method}}
#' @seealso \code{\link{plotCyto,GatingSet-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(name = "plotCyto",
           def = function(x, ...){standardGeneric("plotCyto")}
)

#' plotCyto - flowFrame Method
#'
#' Explore & visualise flow cytometry data for a flowFrame.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1D density distribution or length 2 for 2D scatter plot.
#' @param ... additional arguments passed to
#'   \code{\link{plotCyto1d,flowFrame-method}} or
#'   \code{\link{plotCyto2d,flowFrame-method}}.
#'
#' @examples
#' \dontrun{
#'   fs <- Activation
#'   trans <- estimateLogicle(fs[[2]],"PE-A")
#'   fs <- transform(fs, trans)
#'   plotCyto(fs[[2]], channels = c("PE-A"), transList = trans, overlay = fs[[1]])
#'   plotCyto(fs[[2]], channels = c("FSC-A","PE-A"), transList = trans, overlay = fs[[1]])
#' }
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setMethod(plotCyto, signature = "flowFrame", 
          definition = function(x, channels, ...){

  # Assign x to fr
  fr <- x
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Make call to appropriate plotCyto function
  if(length(channels) == 1){
      
      plotCyto1d(x = fr, channel = channels, ...)
    
  }else if(length(channels) == 2){
      
      plotCyto2d(x = fr, channels = channels, ...)
    
  }
  
})

#' plotCyto - flowSet Method
#' 
#' Explore & visualise flow cytometry data for a flowSet.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1D density distribution or length 2 for 2D scatter plot.
#' @param ... additional arguments passed to
#'   \code{\link{plotCyto1d,flowSet-method}} or
#'   \code{\link{plotCyto2d,flowSet-method}}.
#'
#' @examples
#' \dontrun{
#'   fs <- Activation
#'   trans <- estimateLogicle(fs[[2]],"PE-A")
#'   fs <- transform(fs, trans)
#'   plotCyto(fs, channels = c("PE-A"), transList = trans, overlay = fs[[1]])
#'   plotCyto(fs, channels = c("FSC-A","PE-A"), transList = trans, overlay = fs[[1]])
#' }
#'
#' @seealso \code{\link{plotCyto1d,flowSet-method}}
#' @seealso \code{\link{plotCyto2d,flowSet-method}}
#' 
#' @importFrom flowCore exprs fsApply
#' @importFrom flowWorkspace sampleNames
#' @importFrom BiocGenerics colnames
#' @importFrom methods as
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setMethod(plotCyto, signature = "flowSet", 
          definition = function(x, channels, ...){
  
  # Assign x to fs
  fs <- x
  
  # Check channels
  channels <- checkChannels(fs, channels = channels, plot = TRUE)

  # Make call to appropriate plotCyto function
  if(length(channels) == 1){
      
     plotCyto1d(x = fs, channel = channels, ...)
    
  }else if(length(channels) == 2){
    
     plotCyto2d(x = fs, channels = channels, ...)
      
  }

})

#' plotCyto - GatingSet Method
#'
#' Explore & visualise flow cytometry data for a GatingSet.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population containing the events to plot.
#' @param alias name of the gated population for which the gate should be drawn
#'   on the plot.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1D density histogram or length 2 for 2D density plot. If alias is
#'   supplied the channels must be the same as those use to construct the
#'   gate(s).
#' @param transList object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the supplied flowFrame. This
#'   transList object will be used internally to ensure axes labels of the plot
#'   are appropriately transformed. The transList will NOT be applied to the
#'   flowFrame internally and should be applied to the flowFrame prior to
#'   plotting.
#' @param merge logical indicating whether the flowFrames should be merged prior
#'   to plotting to yield a single plot, set to FALSE by default.
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames}, \code{list of flowSets} or
#'   \code{list of flowFrame lists} containing populations to be overlayed onto
#'   the plot(s).
#' @param stack vector of length 2 indicating offset for samples and number of
#'   samples per plot. Set to \code{c(0,1)} to plot each sample in a separate
#'   panel.
#' @param text.labels vector of names to use in population labels, set to
#'   \code{alias} by default.
#' @param main title to use for the plot, set to the name of the sample and the
#'   name of the parent population by default.
#' @param text.legend vector of labels to use for the legend if legend is TRUE.
#' @param ... additional arguments passed to
#'   \code{\link{plotCyto1d,flowSet-method}} or
#'   \code{\link{plotCyto2d,flowSet-method}}.
#'
#' @examples
#' \dontrun{
#'   fs <- Activation
#'   gs <- GatingSet(fs)
#'   drawGate(gs, parent = "root", alias = "Cells", channels = c("FSC-A","SSC-A"))
#'   plotCyto(gs, parent = "root", alias = "Cells", channels = c("FSC-A","SSC-A"))
#'   drawGate(gs, parent = "Cells", alias = "Single Cells", channels = c("FSC-A","FSC-H"))
#'   plotCyto(gs, parent = "Cells", alias = "Single Cells", channels = c("FSC-A","FSC-H"))
#' }
#'
#' @importFrom flowWorkspace getGate getData prettyAxis sampleNames
#'   getTransformations
#' @importFrom flowCore exprs fsApply filters transformList
#' @importFrom BiocGenerics colnames
#' @importFrom methods as
#'
#' @seealso \code{\link{plotCyto1d,flowSet-method}}
#' @seealso \code{\link{plotCyto2d,flowSet-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setMethod(plotCyto, signature = "GatingSet", 
          definition = function(x, parent, alias, channels, transList = NULL, merge = FALSE, overlay = NULL, stack = c(0,1), text.labels, text.legend, main, ...){
  
  # No Parent Supplied
  if(missing(parent)){
    
    stop("Please supply the name of the parent population to plot.")
    
  }
            
  # Assign x to gs
  gs <- x
  
  # Check supplied alias exists in GatingSet
  if(!missing(alias)){
    
    if(!all(alias %in% basename(getNodes(gs)))){
      
      stop("Supplied alias does not exist in the GatingSet.")
      
    }
    
  }
  
  # Extract alias gate for plotting
  if(!missing(alias)){
    
    gates <- lapply(alias, function(pop) getGate(gs[[1]], pop))
    names(gates) <- alias
    gates <- filters(gates)
  
  }else{
    
    gates <- NULL
    
  }
  
  # transList
  if(is.null(transList)){
    
    trns <- getTransformations(gs[[1]], only.function = TRUE)
    
    if(!length(trns) == 0){
      
      transList <- transformList(names(trns), trns)
      
    }
    
    if(length(trns) == 0){
      
      transList <- NULL
      
    }
    
  }
  
  # Labels?
  if(is.null(gates)){
    
    labels <- FALSE
    
  }
  
  # Text for labels
  if(missing(text.labels) & !missing(alias)){
    
    text.labels <- alias
    
  }

  # Extract population from GatingSet
  fs <- getData(gs, parent)
  
  # Titles
  if(missing(main)){
    
    main <- sapply(1:length(gs), function(gh){
      
      if(parent == "root"){
        
        parent <- "All Events"
        
      }
      
      paste(sampleNames(gs)[gh],"\n",parent)
      
    })
    
    if(length(channels) == 1 & stack[1] != 0){
      
      main <- parent
      
    }
    
    if(merge == TRUE){
      
      main <- parent
      
    }
    
  }
  
  # Legend Text
  if(!is.null(overlay)){
    
    if(class(overlay) == "character"){
      
      if(missing(text.legend)){
        
        text.legend <- c(parent, overlay)
        
      }
      
    }
    
  }
  
  # Overlay
  if(!is.null(overlay)){
    
    # Extract populations to overlay - list of flowSets
    if(class(overlay) == "character"){
    
      overlay <- lapply(overlay, function(overlay){getData(gs, overlay)})
    
      if(merge == TRUE){
        
        overlay <- lapply(overlay, function(x){
          
          fr <- as(x,"flowFrame")
          
          if(is.na(match("Original", BiocGenerics::colnames(fr))) == FALSE){
            
            fr <- suppressWarnings(fr[, -match("Original", BiocGenerics::colnames(fr))])
            
          }
          
          return(fr)
          
          })
        
      }
      
    }
    
  }
  
  # Make call to appropriate plotCyto function
  if(length(channels) == 1){
        
     plotCyto1d(x = fs, channel = channels, overlay = overlay, transList = transList, merge = merge, main = main,  gates = gates, text.legend = text.legend, text.labels = text.labels, stack = stack, ...)
    
  }else if(length(channels) == 2){
        
     plotCyto2d(x = fs, channels = channels, overlay = overlay, transList = transList, merge = merge, main = main, gates = gates, text.legend = text.legend, text.labels = text.labels, ...)
        
  }
  
})

#' plotCyto - GatingHierarchy Method
#'
#' Explore & visualise flow cytometry for a GatingHierarchy.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population containing the events to plot.
#' @param alias name of the gated population for which the gate should be drawn
#'   on the plot.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1D density histogram or length 2 for 2D density plot. If alias is
#'   supplied the channels must be the same as those use to construct the
#'   gate(s).
#' @param transList object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the supplied flowFrame. This
#'   transList object will be used internally to ensure axes labels of the plot
#'   are appropriately transformed. The transList will NOT be applied to the
#'   flowFrame internally and should be applied to the flowFrame prior to
#'   plotting.
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames}, \code{list of flowSets} or
#'   \code{list of flowFrame lists} containing populations to be overlayed onto
#'   the plot(s).
#' @param text.labels vector of names to use in population labels, set to
#'   \code{alias} by default.
#' @param main title to use for the plot, set to the name of the sample and the
#'   name of the parent population by default.
#' @param text.legend vector of labels to use for the legend if legend is TRUE.
#' @param ... additional arguments passed to
#'   \code{\link{plotCyto1d,flowSet-method}} or
#'   \code{\link{plotCyto2d,flowSet-method}}.
#'
#' @examples
#' \dontrun{
#'   fs <- Activation
#'   gs <- GatingSet(fs)
#'   drawGate(gs, parent = "root", alias = "Cells", channels = c("FSC-A","SSC-A"))
#'   plotCyto(gs[[1]], parent = "root", alias = "Cells", channels = c("FSC-A","SSC-A"))
#'   drawGate(gs, parent = "Cells", alias = "Single Cells", channels = c("FSC-A","FSC-H"))
#'   plotCyto(gs[[1]], parent = "Cells", alias = "Single Cells", channels = c("FSC-A","FSC-H"))
#' }
#'
#' @importFrom flowWorkspace getGate getData prettyAxis sampleNames
#'   getTransformations
#' @importFrom flowCore exprs fsApply filters transformList
#' @importFrom BiocGenerics colnames
#' @importFrom methods as
#'
#' @seealso \code{\link{plotCyto1d,flowSet-method}}
#' @seealso \code{\link{plotCyto2d,flowSet-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setMethod(plotCyto, signature = "GatingHierarchy", 
          definition = function(x, parent, alias, channels, transList = NULL, overlay = NULL, text.labels, text.legend, main, ...){
            
  # No Parent Supplied
  if(missing(parent)){
              
    stop("Please supply the name of the parent population to plot.")
              
  }
            
  # Assign x to gh
  gh <- x
            
  # Check supplied alias exists in GatingSet
  if(!missing(alias)){
              
     if(!all(alias %in% basename(getNodes(gh)))){
                
        stop("Supplied alias does not exist in the GatingSet.")
                
      }
              
  }
            
  # Extract alias gate for plotting
  if(!missing(alias)){
              
     gates <- lapply(alias, function(pop) getGate(gh, pop))
     names(gates) <- alias
     gates <- filters(gates)
              
  }else{
              
     gates <- NULL
              
  }
            
  # transList
  if(is.null(transList)){
    
    trns <- getTransformations(gh, only.function = TRUE)
    
    if(!length(trns) == 0){
      
      transList <- transformList(names(trns), trns)
      
    }
    
    if(length(trns) == 0){
      
      transList <- NULL
      
    }
    
  }
            
   # Labels?
   if(is.null(gates)){
              
     labels <- FALSE
              
   }
            
   # Text for labels
   if(missing(text.labels) & !missing(alias)){
              
      text.labels <- alias
              
  }
            
  # Extract population from GatingSet
  fr <- getData(gh, parent)
            
  # Titles
  if(missing(main)){
                
      if(parent == "root"){
                  
         parent <- "All Events"
                  
      }
                
      main <- paste(sampleNames(gh),"\n",parent)
              
  }
            
  # Legend Text
  if(!is.null(overlay)){
              
     if(class(overlay) == "character"){
                
        if(missing(text.legend)){
                  
           text.legend <- c(parent, overlay)
                  
         }
                
     }
              
  }
            
  # Overlay
  if(!is.null(overlay)){
              
    # Extract populations to overlay - list of flowSets
    if(class(overlay) == "character"){
                
        overlay <- lapply(overlay, function(overlay){getData(gh, overlay)})
                
    }
              
  }
            
  # Make call to appropriate plotCyto function
  if(length(channels) == 1){
              
     plotCyto1d(x = fr, channel = channels, overlay = overlay, transList = transList, main = main,  gates = gates, text.legend = text.legend, text.labels = text.labels, ...)
              
  }else if(length(channels) == 2){
              
     plotCyto2d(x = fr, channels = channels, overlay = overlay, transList = transList, main = main, gates = gates, text.legend = text.legend, text.labels = text.labels, ...)
              
  }
            
})