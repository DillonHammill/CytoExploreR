#' Interactive Plots of Flow Cytometry Data.
#'
#' \code{drawPlot} opens an interactive plotting window and constructs either a 1-dimensional (1D) density plot or
#' a 2-dimensional (2D) plot based on teh number of channels supplied.
#' 
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param subSample numeric indicating the number of events to display in 2D plots. 
#' @param cex numeric passed to plotDens to alter the point size (defaults to 3).
#' @param adjust numeric smoothing factor used for 1D density plots (defaults to 1.5).
#' @param ... additional arguments for plotDens.
#'
#' @return \code{NULL}
#'
#' @keywords manual, gating, draw, openCyto, gates, plot, density, flow cytometry
#' @importFrom flowDensity plotDens
#' @importFrom flowCore exprs
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawPlot <- function(fr, channels, subSample = NULL, cex = 3, adjust = 1.5, ...){
  
  # Check supplied channel(s) are valid
  checkChannels(fr, channels)
  
  # Open new plotting window
  if(dev.cur() == 1){  # null device
    dev.new() # RStudioGD
    dev.new() # New window
  }else{
    dev.new()
  }
  
  # 1D Density Plot if length(channels) == 1
  if(length(channels) == 1){
    
    # Calculate Kernel Density
    d <- density(flowCore::exprs(fr)[,channels], adjust = adjust)
    
    # Plot Kernel Density - channel as plot title
    plot(d, main=paste(channels), col.main = "red")
    polygon(d, col="dodgerblue", border="black", lwd = 1)
    
  }else if(length(channels) == 2){
    
    # subSample fr to speed up plotting
    if(is.null(subSample)){
      
    # 2D Plot using flowDensity::plotDens
    flowDensity::plotDens(fr, channels = channels, cex = cex, ...)  
      
    }else{
      
      # Restrict fr to # of subSample events
      fr <- sampleFrame(fr, size = subSample)
      
      # 2D Plot using flowDensity::plotDens
      flowDensity::plotDens(fr, channels = channels, cex = 3, ...)  
      
    }
    
  }
  
}

#' Label Plots with Population Name and Frequency.
#'
#' \code{plotLabels} takes on a \code{flowFrame} object, population name \code{alias}, \code{channels} and a gate object to construct
#' a text label for the plot with the population name and frequency.
#' 
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating.
#' @param alias the name of the gated population.
#' @param gate an object of class \code{polygonGate}, \code{rectangleGate} or \code{ellipsoidGate}.
#' @param ... additional arguments for plotDens.
#'
#' @return add text box label to existing plot with population name and frequency.
#'
#' @keywords manual, gating, draw, openCyto, gates, drawPlot, flow cytometry, labels
#' @importFrom plotrix boxed.labels
#' @importFrom flowCore Subset
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
plotLabels <- function(fr, channels, alias, gate){
  
  # Total events
  events <- nrow(fr)
  
  # Population events
  pop <- nrow(flowCore::Subset(fr, gate))
  
  # Population frequency
  freq <- round(pop/events,4)
  freq <- sprintf("%.2f %%",100*freq)
  
  # Center of gate
  if(class(gate) == "polygonGate"){
    
    center <- c((sum(gate@boundaries[,1])/nrow(gate@boundaries)),(sum(gate@boundaries[,2])/nrow(gate@boundaries)))
  
  }else if(class(gate) == "rectangleGate"){
    
    if("Inf" %in% gate@max){
      
    max <- max(exprs(fr)[,names(gate@max[which(gate@max %in% "Inf")])])
    gate@max[which(gate@max %in% "Inf")] <- max
    
    }
    
    if("-Inf" %in% gate@min){
      
    min <- min(exprs(fr)[,names(gate@min[which(gate@min %in% "-Inf")])])
    gate@min[which(gate@min %in% "-Inf")] <- min
    
    }
    
    if(length(gate@max) == 1){
      
      y.min <- par("usr")[3]
      y.max <- par("usr")[4]
      
      center <- center <- c((sum(gate@min[1],gate@max[1])/2),(sum(y.min,y.max)/2))
      
    }else{
  
    center <- c((sum(gate@min[1],gate@max[1])/2),(sum(gate@min[2],gate@max[2])/2))
    
    }
    
  }else if(class(gate) == "ellipsoidGate"){
    
    center <- gate@mean
  
  }
  
  # Add text annotation to plot with population name & frequency
  plotrix::boxed.labels(x = center[1], y = center[2], labels = paste(alias,freq, sep = "\n"), border = FALSE, font = 2)
}
