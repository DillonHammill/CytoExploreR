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
#' @importFrom flowCore exprs
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
  
  lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2, col = "red")
  abline(v = pts$x, h = pts$y, lwd = 2, col = "red")
  
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

#' Draw Web Gate for Analysis of Flow Cytometry Data.
#'
#' \code{drawWeb} is a variation of drawQuadrant which allows more flexibility with gate co-ordinates (angled lines) and supports a number 
#' of gates as indicated by the \code{alias} argument. To construct the gate supply the center point and surrounding points on plot edge.
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
#' @return a filters list containing the constructed \code{polygonGate} object(s).
#'
#' @keywords manual, gating, draw, FlowJo, polygonGate, openCyto, drawWeb
#' @importFrom flowDensity plotDens
#' @importFrom flowCore polygonGate
#' @importFrom flowCore exprs
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
drawWeb <- function(fr, channels, alias = NULL, subSample = NULL, plot = TRUE, labs = TRUE,...){
  
  if(length(channels) != 2) stop("Two channels are required to construct a web gate.")
  
  if(plot == TRUE){
    
    drawPlot(fr = fr, channels = channels, subSample = subSample, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Select center of the web gate
  message("Select the center of the web gate.")
  center <- locator(n = 1, type = "p", lwd = 2, pch = 16, col = "red")
  
  # User Prompt
  message("Select surrounding co-ordinates on plot edges to draw a web gate.")
  
  # Extract data for use later & Calculate min and max values
  vals <- exprs(fr)[,channels]
  xmin <- round(min(vals[,channels[1]]), 4)
  xmax <- round(max(vals[,channels[1]]), 4)
  ymin <- round(min(vals[,channels[2]]), 4)
  ymax <- round(max(vals[,channels[2]]), 4)
  
  # Get all gate co-ordinates - c(center, others)
  coords <- lapply(seq(1:length(alias)), function(x){
    
    pt <- locator(n = 1, type = "p", lwd = 2, pch = 16, col = "red")
    lines(x = c(center$x, pt$x), y = c(center$y,pt$y), lwd = 2, col = "red")
    
    return(c(pt$x,pt$y))
    
  })
  coords <- as.data.frame(do.call(rbind, coords))
  colnames(coords) <- c("x","y")
  coords <- rbind(center,coords)
  
  # Determine which quadrants the points are in bottom left anti-clockwise to top right (relative to center)
  quads <- c(0,rep(NA, length(alias)))
  for(i in 2:length(coords$x)){
    
    # Bottom left Q1
    if(coords[i,]$x < center$x & coords[i,]$y <= center$y){
      
      quads[i] <- 1
      
      # Bottom right Q2  
    }else if(coords[i,]$x >= center$x & coords[i,]$y < center$y){
      
      quads[i] <- 2
      
      # Top right Q3
    }else if(coords[i,]$x > center$x & coords[i,]$y >= center$y){
      
      quads[i] <- 3
      
      # Top left Q4  
    }else if(coords[i,]$x <= center$x & coords[i,]$y > center$y){
      
      quads[i] <- 4
      
    }
  }
  coords[,"Q"] <- quads
  coords <- coords[with(coords, order(coords$Q)),]
  
  # Push points to plot limits (intersection with plot limits)
  
  # Quadrant 1: find limit intercept and modify point co-ordinates
  q1 <- coords[coords$Q == 1, ]
  for(x in 1:length(q1$Q)){
    
    # Check if line through center and point intersects with vertical axis
    vint <- lines.intercept(c(center$x,center$y), c(q1[x,"x"], q1[x,"y"]), c(xmin, center$y), c(xmin, ymin), interior.only = TRUE)
    
    # If a vertical intercept exists modify co-ordinates of point
    if(!length(vint) == 0 & !anyNA(vint)){
      
      q1[x,c("x","y")] <- vint
      
    }else if(anyNA(vint) | length(vint) == 0){
      
      # If intercept does not exist - check if line through center and point intersects with horizontal axis
      hint <- lines.intercept(c(center$x,center$y), c(q1[x,"x"], q1[x,"y"]), c(center$x, ymin), c(xmin,ymin), interior.only = TRUE)
      
      # If a horizontal intercept exists modify co-ordinates of point
      if(!length(hint) == 0 & !anyNA(hint)){
        
        q1[x,c("x","y")] <- hint
        
      }
    }
  }
  coords[coords$Q == 1, ] <- q1
  
  # Quadrant 2: find limit intercept and modify point co-ordinates
  q2 <- coords[coords$Q == 2, ] 
  for(x in 1:length(q2$Q)){
    
    # Check if line through center and point intersects with vertical axis
    vint <- lines.intercept(c(center$x,center$y), c(q2[x,"x"], q2[x,"y"]), c(xmax, center$y), c(xmax, ymin), interior.only = TRUE)
    
    # If a vertical intercept exists modify co-ordinates of point
    if(!length(vint) == 0 & !anyNA(vint)){
      
      q2[x,c("x","y")] <- vint
      
    }else if(anyNA(vint) | length(vint) == 0){
      
      # If intercept does not exist - check if line through center and point intersects with horizontal axis
      hint <- lines.intercept(c(center$x,center$y), c(q2[x,"x"], q2[x,"y"]), c(center$x, ymin), c(xmax, ymin), interior.only = TRUE)
      
      # If a horizontal intercept exists modify co-ordinates of point
      if(!length(hint) == 0 & !anyNA(hint)){
        
        q2[x,c("x","y")] <- hint
        
      }
    }
  }
  coords[coords$Q == 2, ] <- q2
  
  # Quadrant 3: find limit intercept and modify point co-ordinates
  q3 <- coords[coords$Q == 3, ]
  for(x in 1:length(q3$Q)){
    
    # Check if line through center and point intersects with vertical axis
    vint <- lines.intercept(c(center$x,center$y), c(q3[x,"x"], q3[x,"y"]), c(xmax,ymax), c(xmax,center$y), interior.only = TRUE)
    
    # If a vertical intercept exists modify co-ordinates of point
    if(!length(vint) == 0 & !anyNA(vint)){
      
      q3[x,c("x","y")] <- vint
      
    }else if(anyNA(vint) | length(vint) == 0){
      
      # If intercept does not exist - check if line through center and point intersects with horizontal axis
      hint <- lines.intercept(c(center$x,center$y), c(q3[x,"x"], q3[x,"y"]), c(center$x, ymax), c(xmax, ymax), interior.only = TRUE)
      
      # If a horizontal intercept exists modify co-ordinates of point
      if(!length(hint) == 0 & !anyNA(hint)){
        
        q3[x,c("x","y")] <- hint
        
      }
    }
  }
  coords[coords$Q == 3, ] <- q3
  
  # Quadrant 4: find limit intercept and modify point co-ordinates
  q4 <- coords[coords$Q == 4, ]
  for(x in 1:length(q4$Q)){
    
    # Check if line through center and point intersects with vertical axis
    vint <- lines.intercept(c(center$x,center$y), c(q4[x,"x"], q4[x,"y"]), c(xmin, ymax), c(xmin, center$y), interior.only = TRUE)
    
    # If a vertical intercept exists modify co-ordinates of point
    if(!length(vint) == 0 & !anyNA(vint)){
      
      q4[x,c("x","y")] <- vint
      
    }else if(anyNA(vint) | length(vint) == 0){
      
      # If intercept does not exist - check if line through center and point intersects with horizontal axis
      hint <- lines.intercept(c(center$x,center$y), c(q4[x,"x"], q4[x,"y"]), c(xmin, ymax), c(center$x, ymax), interior.only = TRUE)
      
      # If a horizontal intercept exists modify co-ordinates of point
      if(!length(hint) == 0 & !anyNA(hint)){
        
        q4[x,c("x","y")] <- hint
        
      }
    }
  }
  coords[coords$Q == 4, ] <- q4
  
  # If multiple points in same quadrant order anticlockwise Q1-Q4
  if(anyDuplicated(coords$Q) != 0){
    
    if(coords$Q[duplicated(coords$Q)] == 1){
      
      # Multiple points in Q1 - sort by -y then +x
      q1 <- coords[coords$Q == 1,]
      q1 <- q1[with(q1, order(-q1$y,q1$x))]
      coords[coords$Q == 1,c("x","y")] <- q1
      
    }else if(coords$Q[duplicated(coords$Q)] == 2){
      
      # Multiple points in Q2 - sort by +x then +y
      q2 <- coords[coords$Q == 2,]
      q2 <- q2[with(q2, order(q2$x,q2$y))]
      coords[coords$Q == 2, c("x","y")] <- q2
      
    }else if(coords$Q[duplicated(coords$Q)] == 3){
      
      # Multiple points in Q3 - sort by +y then -x
      q3 <- coords[coords$Q == 3,]
      q3 <- q3[with(q3, order(q3$y,-q3$x))]
      coords[coords$Q == 3, c("x","y")] <- q3
      
    }else if(coords$Q[duplicated(coords$Q)] == 4){
      
      # Multiple points in Q4 - sort by -x then -y
      q4 <- coords[coords$Q == 4,]
      q4 <- q4[with(q4, order(-q4$x,-q4$y))]
      coords[coords$Q == 4, c("x","y")] <- q4
      
    }
  }
  
  # Construct gates using input points
  # Duplicate first point after last point
  coords[(length(coords$Q)+1), ] <- coords[2, ]
  coords[] <- lapply(coords,round,4)
  
  # Gate coordinates using input points
  gates <- list()
  for(i in 2:(length(coords$Q)-1)){
    gates[[i-1]] <- rbind(coords[1,], coords[i,], coords[i+1,])
  }
  
  # Check if a corner lies between the points - add as gate co-ordinate
  # Calculate corner points using min & max values
  Q1 <- c(xmin, ymin, 1)
  Q2 <- c(xmax, ymin, 2)
  Q3 <- c(xmax, ymax, 3)
  Q4 <- c(xmin, ymax, 4)
  Q <- matrix(c(Q1,Q2,Q3,Q4),byrow = TRUE, nrow = 4)
  colnames(Q) <- c("x","y","Q")
  Q <- data.frame(Q)
  
  indx <- 1:(length(alias)-1)
  
  # Add corners to appropriate gates step wise
  gates[indx] <- lapply(gates[indx], function(x){
    
    # DUPLICATION - points in same quadrant
    if(any(duplicated(x$Q))){
      
      # Quadrant 1:
      if(x$Q[duplicated(x$Q)] == 1){
        
        if(x[2,"x"] == xmin & x[3,"x"] != xmin){                                
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1,2),], Q1, x[3,])
          Q <<- Q[-match(1,Q[,"Q"]),]
          
        }
        # Quadrant 2:  
      }else if(x$Q[duplicated(x$Q)] == 2){
        
        if(x[3,"y"] == ymin & x[3,"y"] != ymin){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),],Q2,x[3,])
          Q <<- Q[-match(2,Q[,"Q"]),]
          
        }
        # Quadrant 3:  
      }else if(x$Q[duplicated(x$Q)] == 3){
        
        if(x[2,"x"] == xmax &  x[3,"x"] != xmax){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),],Q3,x[3,])
          Q <<- Q[-match(3,Q[,"Q"]),]
          
        }
        # Quadrant 4:
      }else if(x$Q[duplicated(x$Q)] == 4){
        
        if(x[2,"y"] == ymax & x[3,"y"] != ymax){
          
          # Include Q4 corner in gate
          x <- rbind(x[c(1,2),], Q4, x[3,])
          Q <<- Q[-match(4,Q[,"Q"]),]
          
        }
      }
      
      # ADJACENT - points in adjacent quadrants 
    }else if(any(x[3,"Q"] - x[2,"Q"] == c(0,1,-3))){
      
      # Q1-Q2
      if(x[2,"Q"] == 1 & x[3,"Q"] == 2){
        
        if(x[2,"x"] == xmin & x[3,"x"] == xmax){
          
          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1,2),], Q1, Q2, x[3,])
          Q <<- Q[-match(c(1,2),Q[,"Q"]),]
          
        }else if(x[2,"x"] == xmin & x[3,"x"] != xmax){
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1,2),], Q1, x[3,])
          Q <<- Q[-match(1,Q[,"Q"]),]
          
        }else if(x[2,"x"] != xmin & x[3,"x"] == xmax){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),], Q2, x[3,])
          Q <<- Q[-match(2,Q[,"Q"]),]
          
        }
        
        # Q2-Q3   
      }else if(x[2,"Q"] == 2 & x[3,"Q"] == 3){
        
        if(x[2,"y"] == ymin & x[3,"y"] == ymax){
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, x[3,])
          Q <<- Q[-match(c(2,3),Q[,"Q"]),]
          
        }else if(x[2,"y"] == ymin & x[3,"y"] != ymax){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),], Q2, x[3,])
          Q <<- Q[-match(2,Q[,"Q"]),]
          
        }else if(x[2,"y"] != ymin & x[3,"y"] == ymax){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),], Q3, x[3,])
          Q <<- Q[-match(3,Q[,"Q"]),]
          
        }
        
        # Q3-Q4 
      }else if(x[2,"Q"] == 3 & x[3,"Q"] == 4){
        
        if(x[2,"x"] == xmax & x[3,"x"] == xmin){
          
          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1,2),], Q3, Q4, x[3,])
          Q <<- Q[-match(c(3,4),Q[,"Q"]),]
          
        }else if(x[2,"x"] == xmax & x[3,"x"] != xmin){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),], Q3, x[3,])
          Q <<- Q[-match(3,Q[,"Q"]),]
          
        }else if(x[2,"x"] != xmax & x[3,"x"] == xmin){
          
          # Include Q4 corner in gate
          x <- rbind(x[c(1,2),], Q4, x[3,])
          Q <<- Q[-match(4,Q[,"Q"]),]
          
        }
        
        # Q4-Q1 
      }else if(x[2,"Q"] == 4 & x[3,"Q"] == 1){
        
        if(x[2,"y"] == ymax & x[3,"y"] == ymin){
          
          # Include Q4 & Q1 corner in gate
          x <- rbind(x[c(1,2),], Q4, Q1, x[3,])
          Q <<- Q[-match(c(1,4),Q[,"Q"]),]
          
        }else if(x[2,"y"] == ymax & x[3,"y"] != ymin){
          
          # Include Q4  corner in gate
          x <- rbind(x[c(1,2),], Q4, x[3,])
          Q <<- Q[-match(4,Q[,"Q"]), ]
          
        }else if(x[2,"y"] != ymax & x[3,"y"] == ymin){
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1,2),], Q1, x[3,])
          Q <<- Q[-match(1,Q[,"Q"]),]
          
        }
        
      }
      
      # SEPARATED - points separated by a quadrant   
    }else if(x[3,"Q"] - x[2,"Q"] == 2){
      
      # Q1-Q3
      if(x[2,"Q"] == 1 & x[3,"Q"] == 3){
        
        if(x[2,"x"] == xmin & x[3,"y"] == ymax){
          
          # Include Q1, Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q1, Q2, Q3, x[3,])
          Q <<- Q[-match(c(1,2,3),Q[,"Q"]),]
          
        }else if(x[2,"x"] == xmin & x[3,"y"] != ymax){
          
          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1,2),], Q1, Q2, x[3,])
          Q <<- Q[-match(c(1,2),Q[,"Q"]),]
          
        }else if(x[2,"x"] != xmin & x[3,"y"] == ymax){
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, x[3,])
          Q <<- Q[-match(c(2,3),Q[,"Q"]),]
          
        }else if(x[2,"x"] != xmin & x[3,"y"] != ymax){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),], Q2, x[3,])
          Q <<- Q[-match(2,Q[,"Q"]),]
          
        }
        
        # Q2-Q4
      }else if(x[2,"Q"] == 2 & x[3,"Q"] == 4){
        
        if(x[2,"y"] == ymin & x[3,"x"] == xmin){
          
          # Include Q2, Q3 & Q4 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, Q4, x[3,])
          Q <<- Q[-match(c(2,3,4),Q[,"Q"]),]
          
        }else if(x[2,"y"] == ymin & x[3,"x"] != xmin){
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, x[3,])
          Q <<- Q[-match(c(2,3),Q[,"Q"]),]
          
        }else if(x[2,"y"] != ymin & x[3,"x"] == xmin){
          
          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1,2),], Q3, Q4, x[3,])
          Q <<- Q[-match(c(3,4),Q[,"Q"]),]
          
        }else if(x[2,"y"] != ymin & x[3,"x"] != xmin){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),], Q3, x[3,])
          Q <<- Q[-match(3,Q[,"Q"]),]
          
        }
        
      }
      
    }
    
    return(x)
  })
  
  # Last gate inherits remaining corners
  if(nrow(Q) != 0){
    
  if(length(which(Q[,"Q"] >= gates[[length(alias)]][2,"Q"])) != 0){
    
    g <- Q[which(Q[,"Q"] >= gates[[length(alias)]][2,"Q"]),]
    
  }
  
  if(length(which(Q[,"Q"] < gates[[length(alias)]][2,"Q"])) != 0){
    
    r <- Q[which(Q[,"Q"] < gates[[length(alias)]][2,"Q"]),]
    
  }
    
  if(exists("g") & exists("r")){
    
    Q <- rbind(g,r)
    
  }else if(exists("g") & !exists("r")){
    
    Q <- g
    
  }else if(!exists("g") & exists("r")){
    
    Q <- r
  }
    
  gates[[length(alias)]] <- rbind(gates[[length(alias)]][c(1,2),], Q, gates[[length(alias)]][3,])
  
  }
  
  
  # Construct the gates
  gates <- lapply(seq(1,length(gates), 1), function(x){
    
    coords <- as.matrix(gates[[x]])[,-3]
    colnames(coords) <- channels
    gate <- flowCore::polygonGate(.gate = coords, filterId = alias[x])
    
    if(labs == TRUE){
      plotLabels(fr = fr, alias = alias[x], channels = channels, gate = gate)
    }
    
    return(gate)
  })
  
  gates <- filters(gates)
  return(gates)
}