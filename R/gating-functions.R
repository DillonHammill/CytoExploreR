#' Draw Polygon Gate(s) for Analysis of Flow Cytometry Data.
#'
#' \code{drawPolygon} constructs an interactive plotting window for user to select the coordinates a polygon gate which then constructed
#' into a \code{polygonGate} object and stored in a list.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the constructed \code{polygonGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, polygonGate, openCyto
#' @importFrom flowDensity plotDens
#' @importFrom flowCore polygonGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawPolygon <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(length(channels) != 2) stop("Two channels are required to construct a polygon gate.")
  
  if(plot == TRUE){
    
  drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  gates <- lapply(alias, function(alias){
    
    message(paste("Select at least 3 points to construct a polygon gate around the",alias,"population. \n"))
    
    # Extract gate coordinates
    coords <- locator(type = "o", lwd = 2, pch = 16, col = "red")
    
    if (length(coords$x) < 3) stop("A minimum of 3 points is required to construct a polygon gate.")
    lines(x = coords$x[c(1, length(coords$x))], y = coords$y[c(1, length(coords$x))], lwd = 2, col = "red")
    
    coords <- as.data.frame(coords)
    coords <- as.matrix(coords)
    colnames(coords) <- channels
    
    gate <- flowCore::polygonGate(.gate = coords, filterId = alias)
    
    if(labs == TRUE){
      
    plotLabels(fr = fr, alias = alias, gate = gate)
      
    }
    
    return(gate)
    
  })

  gates <- filters(gates)
  return(gates)
  
}

#' Draw Rectangle Gate(s) for Analysis of Flow Cytometry Data.
#'
#' \code{drawRectangle} constructs an interactive plotting window for user to select the diagonal coordinates of a rectangle which is constructed
#' into a \code{rectangleGate} object and stored in a list.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the constructed \code{rectangleGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto
#' @importFrom flowDensity plotDens
#' @importFrom flowCore rectangleGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawRectangle <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(length(channels) != 2) stop("Two fluorescent channels are required to construct a rectangle gate.")
  
  if(plot == TRUE){
    
    drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  gates <- lapply(alias, function(alias){
    
    message(paste("Select 2 diagonal points to construct a rectangle gate around the",alias,"population. \n"))
    
    # Extract gate coordinates
    coords <- locator(n = 2, type = "p", lwd = 2, pch = 16, col = "red")
    coords <- data.frame(coords)
    coords <- as.matrix(coords)
    colnames(coords) <- channels
    
    rect(xleft = min(coords[,1]), ybottom = min(coords[,2]), xright = max(coords[,1]), ytop = max(coords[,2]), border = "red", lwd = 2)
  
    gate <- flowCore::rectangleGate(.gate = coords, filterId = alias)
    
    if(labs == TRUE){
      
    plotLabels(fr = fr, alias = alias, gate = gate)
      
    }
    
    return(gate)
    })

  gates <- filters(gates)
  return(gates)
}

#' Draw Interval Gate(s) for Analysis of Flow Cytometry Data.
#'
#' \code{drawInterval} constructs an interactive plotting window for user to select the lower and upper bounds of a population which is constructed
#' into a \code{rectangleGate} object and stored in a list. Both 1-D and 2-D interval gates are supported, for 2-D interval gates an additional argument
#' \code{axis} must be supplied to indicate which axis should be gated.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be gated for 2-D interval gates.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the constructed \code{rectangleGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, interval
#' @importFrom flowDensity plotDens
#' @importFrom flowCore rectangleGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawInterval <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, axis = "x", labs = TRUE,...){
  
  if(!length(channels) %in% c(1,2)) stop("Supply 1-2 channels to construct interval gates.")
  
  if(plot == TRUE){
    
    drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  gates <- lapply(alias, function(alias){
    
    message(paste("Select the lower and upper bounds of the",alias,"population to construct an interval gate. \n"))
    
    # Extract gate coordinates
    coords <- locator(n=2, type = "o", lwd = 2, pch = 16, col = "red")
    coords <- data.frame(coords)
    coords <- as.matrix(coords)
    
    if(length(channels) == 1){
      colnames(coords) <- c(channels[1],"Density")
    }else{
      colnames(coords) <- channels
    }
    
    if(axis == "x"){
      abline(v = coords[,1], lwd = 2, col = "red")
    }else if(axis == "y"){
      abline(h = coords[,2], lwd = 2, col = "red")
    }
    
    
    if(axis == "x"){
    if(length(channels) == 1){
      coords <- data.frame(x = coords[,1])
      coords <- as.matrix(coords)
      colnames(coords) <- channels[1]
    }else if(length(channels) == 2){
      coords <- data.frame(x = coords[,1], y = c(-Inf,Inf))
      coords <- as.matrix(coords)
      colnames(coords) <- channels
    }
      gate <- rectangleGate(.gate = coords, filterId = alias)
      
      if(labs == TRUE){
      plotLabels(fr = fr, alias = alias, gate = gate)
      }
      
    }else if(axis == "y"){
      if(length(channels) == 1) stop("Cannot gate y axis if a single channel is supplied.")
      coords <- data.frame(x = c(-Inf,Inf), y = coords[,2])
      coords <- as.matrix(coords)
      colnames(coords) <- channels
        
      gate <- rectangleGate(.gate = coords, filterId = alias)
      
      if(labs == TRUE){
      plotLabels(fr = fr, alias = alias, gate = gate)
      }
    }
    
    return(gate)
    
  })
  
  gates <- filters(gates)
  return(gates)
    
}

#' Draw Threshold Gate(s) for Analysis of Flow Cytometry Data.
#'
#' \code{drawThreshold} constructs an interactive plotting window for user to select the lower bound of a population which is constructed
#' into a \code{rectangleGate} object and stored in a list. Both 1-D and 2-D threshold gates are supported, for 2-D threshold gates all events
#' above the select x and y coordinates are included in the gate.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the constructed \code{rectangleGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, threshold
#' @importFrom flowDensity plotDens
#' @importFrom flowCore rectangleGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawThreshold <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(!length(channels) %in% c(1,2)) stop("Supply 1-2 channels to construct interval gates.")
  
  if(plot == TRUE){
    
    drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  message(paste("Select the lower bound of the",alias,"population to construct a threshold gate. \n"))
  
  if(length(alias) > 1){
    message("Multiple threhold gates are not supported - a single threshold gate will be returned")
  }
  
  # Extract gate coordinates
  coords <- locator(n=1, type = "p", lwd = 2, pch = 16, col = "red")
  
  if(length(channels) == 1){
    pts <- data.frame(x = c(coords$x,Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    abline(v = coords$x, lwd = 2, col = "red")
  }else if(length(channels) == 2){
    pts <- data.frame(x = c(coords$x,Inf), y = c(coords$y,Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rect(xleft = min(coords$x), ybottom = min(coords$y), xright = max(exprs(fr)[,channels[1]]), ytop = max(exprs(fr)[, channels[2]]), border = "red", lwd = 2)
  }
  
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  if(labs == TRUE){
  plotLabels(fr = fr, alias = alias, gate = gate)
  }
  
  gates <- filters(list(gate))
  
}

#' Draw Boundary Gate(s) for Analysis of Flow Cytometry Data.
#'
#' \code{drawBoundary} constructs an interactive plotting window for user to select the upper bound of a population which is constructed
#' into a \code{rectangleGate} object and stored in a list. Both 1-D and 2-D boundary gates are supported, for 2-D boundary gates all events
#' below the select x and y coordinates are included in the gate.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the constructed \code{rectangleGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, boundary
#' @importFrom flowDensity plotDens
#' @importFrom flowCore rectangleGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawBoundary <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(!length(channels) %in% c(1,2)) stop("Supply 1-2 channels to construct interval gates.")
  
  if(plot == TRUE){
    
    drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  message(paste("Select the upper bound of the",alias,"population to construct a boundary gate. \n"))
  
  if(length(alias) > 1){
    message("Multiple boundary gates are not supported - a single boundary gate will be returned")
  }
  
  # Extract gate coordinates
  coords <- locator(n=1, type = "p", lwd = 2, pch = 16, col = "red")
  
  if(length(channels) == 1){
    pts <- data.frame(x = c(-Inf,coords$x))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    abline(v = coords$x, lwd = 2, col = "red")
  }else if(length(channels) == 2){
    pts <- data.frame(x = c(-Inf,coords$x), y = c(-Inf,coords$y))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rect(xleft = min(exprs(fr)[,channels[1]]), ybottom = min(exprs(fr)[,channels[2]]), xright = max(coords$x), ytop = max(coords$y), border = "red", lwd = 2)
  }
  
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  if(labs == TRUE){
  plotLabels(fr = fr, alias = alias, gate = gate)
  }
  
  gates <- filters(list(gate))
  
}

#' Draw Ellipsoid Gate(s) for Analysis of Flow Cytometry Data.
#'
#' \code{drawEllipse} constructs an interactive plotting window for user to select the limits of a population in 2-D (4 points) which is constructed
#' into a \code{ellipsoidGate} object and stored in a list.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the constructed \code{ellipsoidGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, ellipsoidGate, openCyto, ellipse
#' @importFrom flowDensity plotDens
#' @importFrom flowCore ellipsoidGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawEllipse <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(!length(channels) == 2) stop("Supply 2 channels to construct ellipsoid gates.")
  
  if(plot == TRUE){
    
    drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  gates <- lapply(alias, function(alias){
      
    message(paste("Select 4 points to define the limits of the",alias,"population to construct an ellipsoid gate. \n"))
  
    # Extract gate coordinates
    coords <- locator(n=4, type = "p", lwd = 2, pch = 16, col = "red")
    coords <- data.frame(coords)
    
    # Find which points are on major axis
    dst <- as.matrix(stats::dist(coords))
    mj.pts <- coords[which(dst == max(dst), arr.ind = TRUE)[1,],]
    
    # Find which points are on minor axis
    mr.pts <- coords[!coords$x %in% mj.pts$x & !coords$y %in% mj.pts$y,]
    
    # Find center of the major axis
    center <- c((sum(mj.pts$x)/nrow(mj.pts)), (sum(mj.pts$y)/nrow(mj.pts)))
    
    # Find major point which lies above center
    max.pt <- mj.pts[mj.pts$y > center[2] ,]
    
    # Radius of the major axis
    a <- stats::dist(mj.pts)/2
    
    # Radius of the minor axis
    b <- stats::dist(mr.pts)/2
    
    # Angle between horizontal line through center and max.pt
    if(max.pt[1] > center[1]){           # angle < pi/2
      mj.pt.ct <- cbind(max.pt[1],center[2])
      colnames(mj.pt.ct) <- c("x","y")
      adj <- stats::dist(rbind(center,mj.pt.ct))
      angle <- acos(adj/a)
    }else if(max.pt[1] < center[1]){     # angle > pi/2
      mj.pt.ct <- cbind(center[1], max.pt[2])
      colnames(mj.pt.ct) <- c("x","y")
      opp <- stats::dist(as.matrix(rbind(max.pt,mj.pt.ct)))
      angle <- pi/2 + asin(opp/a)
    }
    
    # Covariance matrix
    cinv <- matrix(c(0,0,0,0), nrow = 2, ncol = 2)
    cinv[1,1] <- (((cos(angle)*cos(angle))/(a*a)) + ((sin(angle)*sin(angle))/(b*b)))
    cinv[2,1] <- sin(angle)*cos(angle)*((1/(a*a))-(1/(b*b)))
    cinv[1,2] <- cinv[2,1]
    cinv[2,2] <- (((sin(angle)*sin(angle))/(a*a)) + ((cos(angle)*cos(angle))/(b*b)))
    
    cvm <- solve(cinv)
    
    dimnames(cvm) <- list(channels,channels)
    
    DescTools::DrawEllipse(x = center[1], y = center[2], radius.x = a, radius.y = b, rot = angle, border = "red", lwd = 2)
    
    gate <- ellipsoidGate(.gate = cvm, mean = center, filterId = alias)
    
    if(labs == TRUE){
    plotLabels(fr = fr, alias = alias, gate = gate)
    }
      
    return(gate)
    
    })
  
  gates <- filters(gates)
  
}

#' Draw Quadrant Gates for Analysis of Flow Cytometry Data.
#'
#' \code{drawQuadrants} constructs an interactive plotting window for user to select the crosshair center of 4 populations which is used to construct
#' 4 \code{rectangleGate} objects and stored in a list. Populations are assigned in the following order: bottom left, bottom right, top right and top left.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param alias the name(s) of the populations to be gated. If multiple population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple 
#' gates will be returned. \code{alias} is \code{NULL} by default which will halt the gating routine.
#' @param subSample numeric indicating the number of events to plot to speed up plotting. If \code{subSample} is greater than the total
#' number of events, all events will be plotted which is the default plotting behaviour.
#' @param plot logical indicating whether the data should be plotted. This feature allows for constructing gates of different types over 
#' existing plots which may already contain a different gate type. For demonstration of this feature refer to the package vignette.
#' @param labs logical indicating whether to include \code{plotLabels} for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{plotDens}.
#'
#' @return a list containing the 4 constructed \code{rectangleGate} objects.
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, quadrants
#' @importFrom flowDensity plotDens
#' @importFrom flowCore rectangleGate
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawQuadrants <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(!length(channels) == 2) stop("Supply 2 channels to construct quadrant gates.")
  
  if(plot == TRUE){
    
    drawPlot(fr, channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  if(!length(alias) == 4){
    stop("Supply 4 population names as the alias argument to construct a set of quadrant gates")
  }
  
  message(paste("Select 1 point designating the center point of the populations to construct quadrant gates. \n"))
  
  # Extract points of drawn gate
  pts <- locator(n=1, type = "o", lwd = 2, pch = 16)
  
  lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
  abline(v = pts$x, h = pts$y, lwd = 2)
  
  pts <- as.data.frame(pts)
  colnames(pts) <- channels
  
  # Construct quadrant gates
  
  # Q1 <- Bottom Left
  q1.gate <- data.frame(x = c(-Inf,pts[1,1]), y = c(-Inf, pts[1,2]))
  q1.gate <- as.matrix(q1.gate)
  colnames(q1.gate) <- channels
  q1 <- rectangleGate(.gate = q1.gate, filterId = alias[1])
  
  if(labs == TRUE){
  plotLabels(fr = fr, alias = alias[1], channels = channels, gate = q1)
  }
  
  # Q2 <- Bottom Right
  q2.gate <- data.frame(x = c(pts[1,1], Inf), y = c(-Inf, pts[1,2]))
  q2.gate <- as.matrix(q2.gate)
  colnames(q2.gate) <- channels
  q2 <- rectangleGate(.gate = q2.gate, filterId = alias[2])
  
  if(labs == TRUE){
  plotLabels(fr = fr, alias = alias[2], channels = channels, gate = q2)
  }
  
  # Q3 <- Top Right
  q3.gate <- data.frame(x = c(pts[1,1], Inf), y = c(pts[1,2], Inf))
  q3.gate <- as.matrix(q3.gate)
  colnames(q3.gate) <- channels
  q3 <- rectangleGate(.gate = q3.gate, filterId = alias[3])
  
  if(labs == TRUE){
  plotLabels(fr = fr, alias = alias[3], channels = channels, gate = q3)
  }
  
  # Q4 <- Top Left
  q4.gate <- data.frame(x = c(-Inf, pts[1,1]), y = c(pts[1,2], Inf))
  q4.gate <- as.matrix(q4.gate)
  colnames(q4.gate) <- channels
  q4 <- rectangleGate(.gate = q4.gate, filterId = alias[4])
  
  if(labs == TRUE){
  plotLabels(fr = fr, alias = alias[4], channels = channels, gate = q4)
  }
  
  gates <- filters(list(q1,q2,q3,q4))
  return(gates)

}
