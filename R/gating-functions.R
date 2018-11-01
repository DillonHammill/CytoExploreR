#' Draw Polygon(s) to Gate Flow Cytometry Populations.
#'
#' \code{drawPolygon} constructs an interactive plotting window to allow manual
#' selection of the co-ordinates of a ploygon gate(s) (through mouse click)
#' which are constructed into
#' \code{\link[flowCore:polygonGate]{polygonGate}} objects and stored in a
#' \code{\link[flowCore:filters-class]{filters}} list.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3","CD4")}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:polygonGate]{polygonGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, polygonGate, openCyto
#'
#' @importFrom flowCore polygonGate filters
#' @importFrom graphics locator lines
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawPolygon <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE, ...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias   
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }

  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Construct gates
  gates <- lapply(alias, function(alias){
    
    message(paste("Select at least 3 points to construct a polygon gate around the",alias,"population. \n"))
    
    # Extract gate coordinates
    coords <- locator(type = "o", lwd = 2, pch = 16, col = "red")
    
    if (length(coords$x) < 3) stop("A minimum of 3 points is required to construct a polygon gate.")
    lines(x = coords$x[c(1, length(coords$x))], y = coords$y[c(1, length(coords$x))], lwd = 2.5, col = "red")
    
    coords <- as.data.frame(coords)
    coords <- as.matrix(coords)
    colnames(coords) <- channels
    
    gate <- flowCore::polygonGate(.gate = coords, filterId = alias)
    
    if(labels == TRUE){
      
      plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
      
    }
    
    return(gate)
    
  })
  
  gates <- filters(gates)
  return(gates)
  
}

#' Draw Rectangles to Gate Flow Cytometry Populations.
#'
#' \code{drawRectangle} constructs an interactive plotting window to allow
#' manual selection of the co-ordinates of a rectangle gate(s) (through mouse
#' click) which are constructed into
#' \code{\link[flowCore:rectangleGate]{rectangleGate}} objects and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Simply select 2
#' diagonal co-ordinates to construct the rectangleGate(s).
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate]{rectangleGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, rectangleGate, openCyto
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom flowCore rectangleGate filters
#' @importFrom flowCore exprs
#' @importFrom graphics locator rect
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawRectangle <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE, ...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = TRUE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Construct gates
  gates <- lapply(alias, function(alias){
    
    message(paste("Select 2 diagonal points to construct a rectangle gate around the",alias,"population. \n"))
    
    # Extract gate coordinates
    coords <- locator(n = 2, type = "p", lwd = 2, pch = 16, col = "red")
    coords <- data.frame(coords)
    coords <- as.matrix(coords)
    colnames(coords) <- channels
    
    rect(xleft = min(coords[,1]), ybottom = min(coords[,2]), xright = max(coords[,1]), ytop = max(coords[,2]), border = "red", lwd = 2.5)
  
    gate <- flowCore::rectangleGate(.gate = coords, filterId = alias)
    
    if(labels == TRUE){
      
      plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
      
    }
    
    return(gate)
    
    })

  gates <- filters(gates)
  return(gates)
}

#' Draw Interval(s) to Gate Flow Cytometry Populations.
#'
#' \code{drawInterval} constructs an interactive plotting window for user to
#' select the lower and upper bounds of a population (through mouse click) which
#' is constructed into a
#' \code{\link[flowCore:rectangleGate]{rectangleGate}} object and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Both 1-D and 2-D
#' interval gates are supported, for 2-D interval gates an additional argument
#' \code{axis} must be supplied to indicate which axis should be gated.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate]{rectangleGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, rectangleGate, openCyto, interval
#'
#' @importFrom flowCore rectangleGate filters
#' @importFrom graphics locator abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawInterval <- function(fr, channels, alias = NULL, plot = TRUE, axis = "x", labels = TRUE,...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Construct gates
  gates <- lapply(alias, function(alias){
    
    message(paste("Select the lower and upper bounds of the",alias,"population to construct an interval gate. \n"))
    
    # Extract gate coordinates
    coords <- locator(n=2, type = "o", lwd = 2.5, pch = 16, col = "red")
    coords <- data.frame(coords)
    coords <- as.matrix(coords)
    
    if(length(channels) == 1){
      
      colnames(coords) <- c(channels[1],"Density")
      
    }else{
      
      colnames(coords) <- channels
      
    }
    
    if(axis == "x"){
      
      abline(v = coords[,1], lwd = 2.5, col = "red")
      
    }else if(axis == "y"){
      
      abline(h = coords[,2], lwd = 2.5, col = "red")
      
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
      
      if(labels == TRUE){
        
        plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
      
      }
      
    }else if(axis == "y"){
      
      if(length(channels) == 1) stop("Cannot gate y axis if a single channel is supplied.")
      coords <- data.frame(x = c(-Inf,Inf), y = coords[,2])
      coords <- as.matrix(coords)
      colnames(coords) <- channels
        
      gate <- rectangleGate(.gate = coords, filterId = alias)
      
      if(labels == TRUE){
        
        plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
      
      }
    }
    
    return(gate)
    
  })
  
  gates <- filters(gates)
  return(gates)
    
}

#' Draw Threshold(s) to Gate Flow Cytometry Populations.
#'
#' \code{drawThreshold} constructs an interactive plotting window for user to
#' select the lower bound of a population which is constructed into a
#' \code{\link[flowCore:rectangleGate]{rectangleGate}} object and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Both 1-D and 2-D
#' threshold gates are supported, for 2-D threshold gates all events above the
#' select x and y coordinates are included in the gate. Multiple threshold gates
#' are not currently supported.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. Multiple
#'   \code{threshold} gates are not currently supported. \code{alias} is
#'   \code{NULL} by default which will halt the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate]{rectangleGate}} object.
#'
#' @keywords manual, gating, draw, rectangleGate, openCyto, threshold
#'
#' @importFrom flowCore rectangleGate filters
#' @importFrom graphics locator rect abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawThreshold <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE, ...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  
  # Construct gates
  message(paste("Select the lower bound of the",alias,"population to construct a threshold gate. \n"))
  
  if(length(alias) > 1){
    
    stop("Multiple threhold gates are not supported.")
  
  }
  
  # Extract gate coordinates
  coords <- locator(n=1, type = "p", lwd = 2.5, pch = 16, col = "red")
  
  if(length(channels) == 1){
    
    pts <- data.frame(x = c(coords$x,Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    abline(v = coords$x, lwd = 2.5, col = "red")
    
  }else if(length(channels) == 2){
    
    pts <- data.frame(x = c(coords$x,Inf), y = c(coords$y,Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rect(xleft = min(coords$x), ybottom = min(coords$y), xright = max(exprs(fr)[,channels[1]]), ytop = max(exprs(fr)[, channels[2]]), border = "red", lwd = 2.5)
  
  }
  
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  if(labels == TRUE){
    
    plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
    
  }
  
  gates <- filters(list(gate))
  
}

#' Draw Boundary(s) to Gate Flow Cytometry Populations.
#'
#' \code{drawBoundary} constructs an interactive plotting window for user to
#' select the upper bound of a population which is constructed into a
#' \code{\link[flowCore:rectangleGate]{rectangleGate}} object and stored
#' in a \code{\link[flowCore:filters-class]{filters}} list. Both 1-D and 2-D
#' boundary gates are supported, for 2-D boundary gates all events below the
#' select x and y coordinates are included in the gate. Multiple boundary gates
#' ares not currently supported.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. Multiple boundary
#'   gates ares not currently supported. \code{alias} is \code{NULL} by default
#'   which will halt the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:rectangleGate]{rectangleGate}}
#'   object.
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, boundary
#'
#' @importFrom flowCore rectangleGate filters
#' @importFrom graphics locator rect abline
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawBoundary <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE, ...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Construct gates
  message(paste("Select the upper bound of the",alias,"population to construct a boundary gate. \n"))
  
  if(length(alias) > 1){
    
    stop("Multiple boundary gates are not supported.")
 
  }
  
  # Extract gate coordinates
  coords <- locator(n=1, type = "p", lwd = 2.5, pch = 16, col = "red")
  
  if(length(channels) == 1){
    
    pts <- data.frame(x = c(-Inf,coords$x))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    abline(v = coords$x, lwd = 2.5, col = "red")
    
  }else if(length(channels) == 2){
    
    pts <- data.frame(x = c(-Inf,coords$x), y = c(-Inf,coords$y))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rect(xleft = min(exprs(fr)[,channels[1]]), ybottom = min(exprs(fr)[,channels[2]]), xright = max(coords$x), ytop = max(coords$y), border = "red", lwd = 2.5)
  
  }
  
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  if(labels == TRUE){
    
    plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
    
  }
  
  gates <- filters(list(gate))
  
}

#' Draw Ellipse(s) to Gate Flow Cytometry Populations.
#'
#' \code{drawEllipse} constructs an interactive plotting window for user to
#' select the limits of a population in 2 dimensions (4 points) which is
#' constructed into \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}
#' object and stored in a \code{\link[flowCore:filters-class]{filters}} list.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3","CD4")}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, ellipsoidGate, openCyto, ellipse
#'
#' @importFrom flowCore ellipsoidGate filters
#' @importFrom graphics locator polygon
#' @importFrom methods as
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawEllipse <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE,...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Construct gates
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
    mj.center <- c((sum(mj.pts$x)/nrow(mj.pts)), (sum(mj.pts$y)/nrow(mj.pts)))
    
    # Find center of all points
    center <- c(sum(c(mj.pts$x,mr.pts$x))/4, sum(c(mj.pts$y,mr.pts$y))/4)
    
    # Adjust mj.pts to fall on center
    adj <- c((mj.center[1] - center[1]),(mj.center[2] - center[2]))
    mj.pts$x <- mj.pts$x - adj[1]
    mj.pts$y <- mj.pts$y - adj[2]
    
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
    
    gate <- ellipsoidGate(.gate = cvm, mean = center, filterId = alias)
    
    polygon(as(gate, "polygonGate")@boundaries[,1], as(gate, "polygonGate")@boundaries[,2], border = "red", lwd = 2.5)
    
    if(labels == TRUE){
      
      plotLabels(x = fr, gates = gate, channels = channels, alias = alias, cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
      
    }
      
    return(gate)
    
    })
  
  gates <- filters(gates)
  
}

#' Draw Quadrants to Gate Multiple Flow Cytometry Populations.
#'
#' \code{drawQuadrants} constructs an interactive plotting window for user to
#' select the crosshair center of 4 populations which is used to construct 4
#' \code{\link[flowCore:rectangleGate]{rectangleGate}} objects which are
#' stored in a\code{\link[flowCore:filters-class]{filters}}  list. Populations
#' are assigned in the following order: bottom left, bottom right, top right and
#' top left.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the 4 populations to be gated. \code{alias} is
#'   \code{NULL} by default which will halt the gating routine.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the 4
#'   constructed \code{\link[flowCore:rectangleGate]{rectangleGate}}
#'   objects.
#'
#' @keywords manual, gating, draw, FlowJo, rectangleGate, openCyto, quadrants
#'
#' @importFrom flowCore rectangleGate filters
#' @importFrom graphics locator lines abline
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawQuadrants <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE, ...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  if(!length(alias) == 4){
    stop("Supply 4 population names as the alias argument to construct a set of quadrant gates")
  }
  
  # Construct gates
  message(paste("Select 1 point designating the center point of the populations to construct quadrant gates. \n"))
  
  # Extract points of drawn gate
  pts <- locator(n=1, type = "o", lwd = 2, pch = 16)
  
  lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2.5, col = "red")
  abline(v = pts$x, h = pts$y, lwd = 2.5, col = "red")
  
  pts <- as.data.frame(pts)
  colnames(pts) <- channels
  
  # Construct quadrant gates
  
  # Q1 <- Bottom Left
  q1.gate <- data.frame(x = c(-Inf,pts[1,1]), y = c(-Inf, pts[1,2]))
  q1.gate <- as.matrix(q1.gate)
  colnames(q1.gate) <- channels
  q1 <- rectangleGate(.gate = q1.gate, filterId = alias[1])
  
  if(labels == TRUE){
    
    plotLabels(x = fr, gates = q1, channels = channels, alias = alias[1], cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
  
  }
  
  # Q2 <- Bottom Right
  q2.gate <- data.frame(x = c(pts[1,1], Inf), y = c(-Inf, pts[1,2]))
  q2.gate <- as.matrix(q2.gate)
  colnames(q2.gate) <- channels
  q2 <- rectangleGate(.gate = q2.gate, filterId = alias[2])
  
  if(labels == TRUE){
    
    plotLabels(x = fr, gates = q2, channels = channels, alias = alias[2], cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
  
  }
  
  # Q3 <- Top Right
  q3.gate <- data.frame(x = c(pts[1,1], Inf), y = c(pts[1,2], Inf))
  q3.gate <- as.matrix(q3.gate)
  colnames(q3.gate) <- channels
  q3 <- rectangleGate(.gate = q3.gate, filterId = alias[3])
  
  if(labels == TRUE){
    
    plotLabels(x = fr, gates = q3, channels = channels, alias = alias[3], cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
    
  }
  
  # Q4 <- Top Left
  q4.gate <- data.frame(x = c(-Inf, pts[1,1]), y = c(pts[1,2], Inf))
  q4.gate <- as.matrix(q4.gate)
  colnames(q4.gate) <- channels
  q4 <- rectangleGate(.gate = q4.gate, filterId = alias[4])
  
  if(labels == TRUE){
    
    plotLabels(x = fr, gates = q4, channels = channels, alias = alias[4], cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
  
  }
  
  gates <- filters(list(q1,q2,q3,q4))
  return(gates)

}

#' Draw Web to Gate Multiple Flow Cytometry Populations.
#'
#' \code{drawWeb} is a variation of drawQuadrant which allows more flexibility
#' with gate co-ordinates (angled lines) and supports any number of gates as
#' indicated by the \code{alias} argument. To construct the gate simply select
#' the center point and surrounding divider points on plot edge. \code{drawWeb}
#' will contruct the \code{\link[flowCore:polygonGate]{polygonGate}}
#' objects and store them in a \code{\link[flowCore:filters-class]{filters}}
#' list.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for plotting and gating.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3","CD4")}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine. Recommended for 3 or more populations.
#' @param plot logical indicating whether the data should be plotted. This
#'   feature allows for constructing gates of different types over existing
#'   plots which may already contain a different gate type.
#' @param labels logical indicating whether to include \code{\link{plotLabels}}
#'   for the gated population(s), \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return a\code{\link[flowCore:filters-class]{filters}} list containing the
#'   constructed \code{\link[flowCore:polygonGate]{polygonGate}}
#'   object(s).
#'
#' @keywords manual, gating, draw, polygonGate, openCyto, drawWeb
#'
#' @importFrom flowCore polygonGate filters
#' @importFrom flowCore exprs
#' @importFrom graphics locator lines
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#' @seealso \code{\link{drawGate}}
#'
#' @export
drawWeb <- function(fr, channels, alias = NULL, plot = TRUE, labels = TRUE, ...){
  
  # Check channels
  channels <- checkChannels(fr, channels = channels, plot = TRUE)
  
  # Check alias
  if(is.null(alias)){
    
    stop("Please supply a name for a the gated population as the alias argument.")
    
  }
  
  # Call new plot?
  if(plot == TRUE){
    
    plotCyto(fr, channels = channels, popup = TRUE, legend = FALSE, labels = FALSE, ...)
    
  }else if(plot == FALSE){
    
  }
  
  # Construct gates
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
    
    pt <- locator(n = 1, type = "p", lwd = 2.5, pch = 16, col = "red")
    lines(x = c(center$x, pt$x), y = c(center$y,pt$y), lwd = 2.5, col = "red")
    
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
  if(1 %in% coords$Q){
  q1 <- coords[coords$Q == 1, ]
  for(x in 1:length(q1$Q)){
    
    # Calculate intersection with horizontal and vertical axes
    vint <- linesIntercept(c(center$x,center$y), c(q1[x,"x"], q1[x,"y"]), c(xmin, center$y), c(xmin, ymin))
    hint <- linesIntercept(c(center$x,center$y), c(q1[x,"x"], q1[x,"y"]), c(center$x, ymin), c(xmin, ymin))

    # Check which axis the point should be pushed onto
    if(vint[2] >= ymin){
      
      q1[x,c("x","y")] <- vint
      
    }else if(vint[2] < ymin){
      
      q1[x,c("x","y")] <- hint
      
    }
  }
  coords[coords$Q == 1, ] <- q1
  }
  
  # Quadrant 2: find limit intercept and modify point co-ordinates
  if(2 %in% coords$Q){
  q2 <- coords[coords$Q == 2, ] 
  for(x in 1:length(q2$Q)){
    
    # Calculate intersection with horizontal and vertical axes
    vint <- linesIntercept(c(center$x,center$y), c(q2[x,"x"], q2[x,"y"]), c(xmax, center$y), c(xmax, ymin))
    hint <- linesIntercept(c(center$x,center$y), c(q2[x,"x"], q2[x,"y"]), c(center$x, ymin), c(xmax, ymin))
    
    # Check which axis the point should be pushed onto
    if(vint[2] >= ymin){
      
      q2[x,c("x","y")] <- vint
      
    }else if(vint[2] < ymin){
      
      q2[x,c("x","y")] <- hint
      
    }
  }
  coords[coords$Q == 2, ] <- q2
  }
  
  # Quadrant 3: find limit intercept and modify point co-ordinates
  if(3 %in% coords$Q){
  q3 <- coords[coords$Q == 3, ]
  for(x in 1:length(q3$Q)){
    
    # Calculate intersection with horizontal and vertical axes
    vint <- linesIntercept(c(center$x,center$y), c(q3[x,"x"], q3[x,"y"]), c(xmax,ymax), c(xmax,center$y))
    hint <- linesIntercept(c(center$x,center$y), c(q3[x,"x"], q3[x,"y"]), c(center$x, ymax), c(xmax, ymax))
    
    # Check which axis the point should be pushed onto
    if(vint[2] >= ymax){
      
      q3[x,c("x","y")] <- hint
      
    }else if(vint[2] < ymax){
      
      q3[x,c("x","y")] <- vint
      
    }
  }
  coords[coords$Q == 3, ] <- q3
  }
  
  # Quadrant 4: find limit intercept and modify point co-ordinates
  if(4 %in% coords$Q){
  q4 <- coords[coords$Q == 4, ]
  for(x in 1:length(q4$Q)){
    
    # Calculate intersection with horizontal and vertical axes
    vint <- linesIntercept(c(center$x,center$y), c(q4[x,"x"], q4[x,"y"]), c(xmin, ymax), c(xmin, center$y))
    hint <- linesIntercept(c(center$x,center$y), c(q4[x,"x"], q4[x,"y"]), c(xmin, ymax), c(center$x, ymax))
    
    # Check which axis the point should be pushed onto
    if(vint[2] >= ymax){
      
      q4[x,c("x","y")] <- hint
      
    }else if(vint[2] < ymax){
      
      q4[x,c("x","y")] <- vint
      
    }
  }
  coords[coords$Q == 4, ] <- q4
  }
  
  # If multiple points in same quadrant order anticlockwise Q1-Q4
  if(anyDuplicated(coords$Q) != 0){
    
    # Quadrant 1
    if(1 %in% coords$Q[duplicated(coords$Q)]){
      
      # Multiple points in Q1 - sort by -y then +x
      q1 <- coords[coords$Q == 1,]
      q1 <- q1[with(q1, order(-q1$y,q1$x)),]
      coords[coords$Q == 1, c("x","y")] <- q1[,c("x","y")]
      
    }
    
    # Quadrant 2
    if(2 %in% coords$Q[duplicated(coords$Q)]){
      
      # Multiple points in Q2 - sort by +x then +y
      q2 <- coords[coords$Q == 2,]
      q2 <- q2[with(q2, order(q2$x,q2$y)),]
      coords[coords$Q == 2, c("x","y")] <- q2[,c("x","y")]
      
    }
    
    # Quadrant 3
    if(3 %in% coords$Q[duplicated(coords$Q)]){
      
      # Multiple points in Q3 - sort by +y then -x
      q3 <- coords[coords$Q == 3,]
      q3 <- q3[with(q3, order(q3$y,-q3$x)),]
      coords[coords$Q == 3, c("x","y")] <- q3[,c("x","y")]
      
    }
    
    # Quadrant 4
    if(4 %in% coords$Q[duplicated(coords$Q)]){
      
      # Multiple points in Q4 - sort by -x then -y
      q4 <- coords[coords$Q == 4,]
      q4 <- q4[with(q4, order(-q4$x,-q4$y)),]
      coords[coords$Q == 4, c("x","y")] <- q4[,c("x","y")]
      
    }
  }
  
  # Construct gates using input points
  # Duplicate first point after last point
  coords[(length(coords$Q) + 1), ] <- coords[2, ]
  coords[] <- lapply(coords, round, 4)
  
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
  
  # Add corners to appropriate gates step-wise
  gates[indx] <- lapply(gates[indx], function(x){
    
    # DUPLICATION - points in same quadrant
    if(any(duplicated(x$Q))){
      
      # Quadrant 1
      if(1 %in% x$Q[duplicated(x$Q)]){
        
        if(x[2,"x"] == xmin & x[3,"x"] != xmin){                                
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1,2),], Q1, x[3,])
          
          # Remove Q1 from Q
          if(1 %in% Q[,"Q"]){
            
          Q <<- Q[-match(1,Q[,"Q"]),]
          
          }
          
        }
      }
      
      # Quadrant 2
      if(2 %in% x$Q[duplicated(x$Q)]){
        
        if(x[2,"y"] == ymin & x[3,"y"] != ymin){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),],Q2,x[3,])
          
          # Remove Q2 from Q
          if(2 %in% Q[,"Q"]){
            
          Q <<- Q[-match(2,Q[,"Q"]),]
          
          }
          
        }
      }
    
      # Quadrant 3
      if(3 %in% x$Q[duplicated(x$Q)]){
        
        if(x[2,"x"] == xmax &  x[3,"x"] != xmax){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),],Q3,x[3,])
          
          # Remove Q3 from Q
          if(3 %in% Q[,"Q"]){
            
          Q <<- Q[-match(3,Q[,"Q"]),]
          
          }
          
        }
      }
    
      # Quadrant 4
      if(4 %in% x$Q[duplicated(x$Q)]){
        
        if(x[2,"y"] == ymax & x[3,"y"] != ymax){
          
          # Include Q4 corner in gate
          x <- rbind(x[c(1,2),], Q4, x[3,])
          
          # Remove Q4 from Q
          if(4 %in% Q[,"Q"]){
            
          Q <<- Q[-match(4,Q[,"Q"]),]
          
          }
          
        }
      }
      
      # ADJACENT - points in adjacent quadrants 
    }else if(any(x[3,"Q"] - x[2,"Q"] == c(0,1))){
      
      # Q1-Q2
      if(x[2,"Q"] == 1 & x[3,"Q"] == 2){
        
        if(x[2,"x"] == xmin & x[3,"x"] == xmax){
          
          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1,2),], Q1, Q2, x[3,])
          
          # Remove Q1 and Q2 from Q
          if(any(c(1,2) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(1,2), Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] == xmin & x[3,"x"] != xmax){
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1,2),], Q1, x[3,])
          
          # Remove Q1 from Q
          if(any(1 %in% Q[,"Q"])){
            
          Q <<- Q[-match(1, Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] != xmin & x[3,"x"] == xmax){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),], Q2, x[3,])
          
          # Remove Q2 from Q
          if(any(2 %in% Q[,"Q"])){
            
          Q <<- Q[-match(2, Q[,"Q"]),]
          
          }
          
        }
        
        # Q2-Q3   
      }else if(x[2,"Q"] == 2 & x[3,"Q"] == 3){
        
        if(x[2,"y"] == ymin & x[3,"y"] == ymax){
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, x[3,])
          
          # Remove Q2 and Q3 from Q
          if(any(c(2,3) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(2,3),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"y"] == ymin & x[3,"y"] != ymax){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),], Q2, x[3,])
          
          # Remove Q2 from Q
          if(any(2 %in% Q[,"Q"])){
            
          Q <<- Q[-match(2,Q[,"Q"]),]
          
          }
          
        }else if(x[2,"y"] != ymin & x[3,"y"] == ymax){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),], Q3, x[3,])
          
          # Remove Q3 from Q
          if(any(3 %in% Q[,"Q"])){
            
          Q <<- Q[-match(3,Q[,"Q"]),]
          
          }
          
        }
        
        # Q3-Q4 
      }else if(x[2,"Q"] == 3 & x[3,"Q"] == 4){
        
        if(x[2,"x"] == xmax & x[3,"x"] == xmin){
          
          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1,2),], Q3, Q4, x[3,])
          
          # Remove Q3 and Q4 from Q
          if(any(c(3,4) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(3,4),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] == xmax & x[3,"x"] != xmin){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),], Q3, x[3,])
          
          # Remove Q3 from Q
          if(any(3 %in% Q[,"Q"])){
            
          Q <<- Q[-match(3,Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] != xmax & x[3,"x"] == xmin){
          
          # Include Q4 corner in gate
          x <- rbind(x[c(1,2),], Q4, x[3,])
          
          # Remove Q4 from Q
          if(any(4 %in% Q[,"Q"])){
            
          Q <<- Q[-match(4,Q[,"Q"]),]
          
          }
          
        }
        
      }
      
      # SEPARATED - points separated by a quadrant   
    }else if(x[3,"Q"] - x[2,"Q"] == 2){
      
      # Q1-Q3
      if(x[2,"Q"] == 1 & x[3,"Q"] == 3){
        
        if(x[2,"x"] == xmin & x[3,"y"] == ymax){
          
          # Include Q1, Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q1, Q2, Q3, x[3,])
          
          # Remove Q1, Q2 and Q3 from Q
          if(any(c(1,2,3) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(1,2,3),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] == xmin & x[3,"y"] != ymax){
          
          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1,2),], Q1, Q2, x[3,])
          
          # Remove Q1 and Q2 from Q
          if(any(c(1,2) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(1,2),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] != xmin & x[3,"y"] == ymax){
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, x[3,])
          
          # Remove Q2 and Q3 from Q
          if(any(c(2,3) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(2,3),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"x"] != xmin & x[3,"y"] != ymax){
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1,2),], Q2, x[3,])
          
          # Remove Q2 from Q
          if(any(2 %in% Q[,"Q"])){
            
          Q <<- Q[-match(2,Q[,"Q"]),]
          
          }
          
        }
        
        # Q2-Q4
      }else if(x[2,"Q"] == 2 & x[3,"Q"] == 4){
        
        if(x[2,"y"] == ymin & x[3,"x"] == xmin){
          
          # Include Q2, Q3 & Q4 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, Q4, x[3,])
          
          # Remove Q2, Q3 and Q4 from Q
          if(any(c(2,3,4) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(2,3,4),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"y"] == ymin & x[3,"x"] != xmin){
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1,2),], Q2, Q3, x[3,])
          
          # Remove Q2 and Q3 from Q
          if(any(c(2,3) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(2,3),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"y"] != ymin & x[3,"x"] == xmin){
          
          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1,2),], Q3, Q4, x[3,])
          
          # Remove Q3 and Q4
          if(any(c(3,4) %in% Q[,"Q"])){
            
          Q <<- Q[-match(c(3,4),Q[,"Q"]),]
          
          }
          
        }else if(x[2,"y"] != ymin & x[3,"x"] != xmin){
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1,2),], Q3, x[3,])
          
          # Remove Q3 from Q
          if(any(3 %in% Q[,"Q"])){
            
          Q <<- Q[-match(3,Q[,"Q"]),]
          
          
          }
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
    
    if(labels == TRUE){
      
      plotLabels(x = fr, gates = gate, channels = channels, alias = alias[x], cex.text = 1, format.text = c("alias","percent"), alpha = 0.7)
    
    }
    
    return(gate)
  })
  
  gates <- filters(gates)
  
  plotGates(gates = gates, channels = channels, gate.col = "red", gate.lwd = 2.5)
  
  return(gates)
  
}