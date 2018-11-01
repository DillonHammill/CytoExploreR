#' Plot Gate Objects onto an Existing Plot
#'
#' @param x gate object of class
#'   \code{\link[flowCore:rectangleGate]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}, \code{list} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{plotGates,rectangleGate-method}}
#' @seealso \code{\link{plotGates,polygonGate-method}}
#' @seealso \code{\link{plotGates,ellipsoidGate-method}}
#' @seealso \code{\link{plotGates,list-method}}
#' @seealso \code{\link{plotGates,filters-method}}
#'
#' @export
setGeneric(name = "plotGates",
           def = function(x, ...){standardGeneric("plotGates")}
)

#' Plot rectangleGate Objects onto an Existing plot
#'
#' @param x an object of class
#'   \code{\link[flowCore:rectangleGate]{rectangleGate}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param col.gate indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param lwd.gate numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param lty.gate integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param pts.gate logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param pch.gate integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param cex.gate numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @importFrom flowCore parameters
#' @importFrom graphics par rect lines abline points
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{plotGates,polygonGate-method}}
#' @seealso \code{\link{plotGates,ellipsoidGate-method}}
#' @seealso \code{\link{plotGates,list-method}}
#' @seealso \code{\link{plotGates,filters-method}}
#'
#' @export
setMethod(plotGates, signature = "rectangleGate", definition = function(x, channels, col.gate = "red", lwd.gate = 2.5, lty.gate = 1, pts.gate = FALSE, pch.gate = 16, cex.gate = 1){
  
  # Assign x to gt
  gt <- x
  
  # Check Channels
  if(missing(channels)){
    
    channels <- as.vector(parameters(gt))
    
  }
  
  if(!all(channels %in% as.vector(parameters(gt)))){
    
    stop("Channels used to construct the plot do not match those of the supplied gate.")
    
  }
  
  # 2D rectangleGate
  if(length(gt@min) == 2){
    
    # Replace -Inf x values for plotting
    if(is.infinite(gt@min[1])){
      
      gt@min[1] <- par("usr")[1]
      
    }
    
    # Replace Inf x values for plotting
    if(is.infinite(gt@max[1])){
      
      gt@max[1] <- par("usr")[2]
      
    }
    
    # Replace -Inf y values for plotting
    if(is.infinite(gt@min[2])){
      
      gt@min[2] <- par("usr")[3]
      
    }
    
    # Replace Inf y values for plotting
    if(is.infinite(gt@max[2])){
      
      gt@max[2] <- par("usr")[4]
      
    }
    
    if(pts.gate == TRUE){
      
      points(x = c(gt@min[channels[1]],gt@max[channels[1]]), y = c(gt@min[channels[2]],gt@max[channels[2]]), col = col.gate, pch = pch.gate, cex = cex.gate)
      
    }
    
    rect(xleft = gt@min[channels[1]], ybottom = gt@min[channels[2]], xright = gt@max[channels[1]], ytop = gt@max[channels[2]], border = col.gate, lwd = lwd.gate, lty = lty.gate)
    
  }else if(length(gt@min) == 1){
    
    # Replace -Inf values for plotting
    if(is.infinite(gt@min[1])){
      
      gt@min[1] <- par("usr")[1]
      
    }
    
    # Replace Inf values for plotting
    if(is.infinite(gt@max[1])){
      
      gt@max[1] <- par("usr")[2]
      
    }
    
    # height of horizontal line
    hln <- 0.5*par("usr")[4]
    
    # Add points (x1,hln) and (x2, hln)
    if(pts.gate == TRUE){
      
      points(x = c(gt@min[channels[1]], gt@max[channels[1]]), y = c(hln,hln), col = col.gate, pch = pch.gate, cex = cex.gate)
      
    }
    
    # Add vertical lines through x1 and x2
    abline(v = c(gt@min[channels[1]], gt@max[channels[1]]), lwd = lwd.gate, col = col.gate, lty = lty.gate)
    
    # Add horizontal line
    lines(x = c(gt@min[channels[1]], gt@max[channels[1]]), y = c(hln,hln), col = col.gate, lwd = lwd.gate, lty = lty.gate)
    
  }
  
})

#' Plot polygonGate Objects onto an Existing Plot
#'
#' @param x an object of class \code{\link[flowCore:polygonGate]{polygonGate}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param col.gate indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param lwd.gate numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param lty.gate integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param pts.gate logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param pch.gate integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param cex.gate numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @importFrom flowCore parameters
#' @importFrom graphics par points polygon
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{plotGates,rectangleGate-method}}
#' @seealso \code{\link{plotGates,ellipsoidGate-method}}
#' @seealso \code{\link{plotGates,list-method}}
#' @seealso \code{\link{plotGates,filters-method}}
#'
#' @export
setMethod(plotGates, signature = "polygonGate", definition = function(x, channels, col.gate = "red", lwd.gate = 2.5, lty.gate = 1, pts.gate = FALSE, pch.gate = 16, cex.gate = 1){
  
  # Assign x to gt
  gt <- x
  
  # Check Channels
  if(missing(channels)){
    
    channels <- as.vector(parameters(gt))
    
  }
  
  if(!all(channels %in% as.vector(parameters(gt)))){
    
    stop("Channels used to construct the plot do not match those of the supplied gate.")
    
  }
  
  # Replace Inf values with plot limits
  if(!all(is.finite(gt@boundaries))){
    
    cnt <- 0
    lapply(1:ncol(gt@boundaries), function(x){
      
      cnt <<- cnt + 1
      if(any(!is.finite(gt@boundaries[,x]) & any(gt@boundaries[,x] < 0))){
        
        if(cnt == 1){
          
          gt@boundaries[,x][which(gt@boundaries[,x] < 0)] <<- par("usr")[1]
          
        }else if(cnt == 2){
          
          gt@boundaries[,x][which(gt@boundaries[,x] < 0)] <<- par("usr")[3]
          
        }
        
      }else if(any(!is.finite(gt@boundaries[,x]) & any(!gt@boundaries[,x] < 0))){
        
        if(cnt == 1){
          
          gt@boundaries[,x][which(!is.finite(gt@boundaries[,x]) & !gt@boundaries[,x] < 0)] <<- par("usr")[2]
          
        }else if(cnt == 2){
          
          gt@boundaries[,x][which(!is.finite(gt@boundaries[,x]) & !gt@boundaries[,x] < 0)] <<- par("usr")[4]
          
        }
        
      }
      
    })
    
  }
  
  # Plot Gate
  if(pts.gate == TRUE){
    
    points(x = c(gt@boundaries[,channels[1]]), y = c(gt@boundaries[,channels[2]]), pch = pch.gate, col = col.gate, cex = cex.gate)
    
  } 
  
  polygon(gt@boundaries[,channels[1]], gt@boundaries[,channels[2]], border = col.gate, lwd = lwd.gate, lty = lty.gate)
  
  
})

#' Plot ellipsoidGate Objects onto an Existing Plot
#'
#' @param x an object of class
#'   \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param col.gate indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param lwd.gate numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param lty.gate integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param pts.gate logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param pch.gate integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param cex.gate numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowCore parameters
#' @importFrom methods as
#' @importFrom graphics polygon
#'
#' @seealso \code{\link{plotGates,rectangleGate-method}}
#' @seealso \code{\link{plotGates,polygonGate-method}}
#' @seealso \code{\link{plotGates,list-method}}
#' @seealso \code{\link{plotGates,filters-method}}
#'
#' @export
setMethod(plotGates, signature = "ellipsoidGate", definition = function(x, channels, col.gate = "red", lwd.gate = 2.5, lty.gate = 1, pts.gate = FALSE, pch.gate = 16, cex.gate = 1){
  
  # Assign x to gt
  gt <- x
  
  # Check Channels
  if(missing(channels)){
    
    channels <- as.vector(parameters(gt))
    
  }
  
  if(!all(channels %in% as.vector(parameters(gt)))){
    
    stop("Channels used to construct the plot do not match those of the supplied gate.")
    
  }
  
  # Coerce to polygonGate
  gt <- as(gt, "polygonGate")
  
  # Plot gate
  polygon(gt@boundaries[,channels[1]], gt@boundaries[,channels[2]], border = col.gate, lwd = lwd.gate, lty = lty.gate)
  
})

#' Plot List of Gate Objects onto an Existing Plot
#'
#' @param x an object of class \code{list} containing objects of class
#'   \code{rectangleGate}, \code{polygonGate} or \code{ellipsoidGate}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param col.gate indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param lwd.gate numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param lty.gate integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param pts.gate logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param pch.gate integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param cex.gate numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{plotGates,rectangleGate-method}}
#' @seealso \code{\link{plotGates,polygonGate-method}}
#' @seealso \code{\link{plotGates,ellipsoidGate-method}}
#' @seealso \code{\link{plotGates,filters-method}}
#'
#' @export
setMethod(plotGates, signature = "list", definition = function(x, channels, col.gate = "red", lwd.gate = 2.5, lty.gate = 1, pts.gate = FALSE, pch.gate = 16, cex.gate = 1){
  
  # Assign x to gts
  gts <- x
  
  # Check Channels
  if(missing(channels)){
    
    channels <- as.vector(parameters(gts[[1]]))
    
  }
  
  if(!all(channels %in% as.vector(parameters(gts[[1]])))){
    
    stop("Channels used to construct the plot do not match those of the supplied gate.")
    
  }
  
  # Gate colours
  if(length(col.gate) != length(gts)){
    
    if(length(col.gate) == 1){
      
      col.gate <- rep(col.gate, length(gts))
      
    }else if(length(col.gate) > length(gts)){
      
      col.gate <- col.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Line Type
  if(length(lty.gate) != length(gts)){
    
    if(length(lty.gate) == 1){
      
      lty.gate <- rep(lty.gate, length(gts))
      
    }else if(length(lty.gate) > length(gts)){
      
      lty.gate <- lty.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Line Width
  if(length(lwd.gate) != length(gts)){
    
    if(length(lwd.gate) == 1){
      
      lwd.gate <- rep(lwd.gate, length(gts))
      
    }else if(length(lwd.gate) > length(gts)){
      
      lwd.gate <- lwd.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Point Character
  if(length(pch.gate) != length(gts)){
    
    if(length(pch.gate) == 1){
      
      pch.gate <- rep(pch.gate, length(gts))
      
    }else if(length(pch.gate) > length(gts)){
      
      pch.gate <- pch.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Point Expansion
  if(length(cex.gate) != length(gts)){
    
    if(length(cex.gate) == 1){
      
      cex.gate <- rep(cex.gate, length(gts))
      
    }else if(length(cex.gate) > length(gts)){
      
      cex.gate <- cex.gate[1:length(gts)]
      
    }
    
  }
  
  # Plot Gates
  mapply(function(gt, col.gate, lty.gate, lwd.gate, pch.gate, cex.gate){
    
    plotGates(gt, channels = channels, col.gate = col.gate, lwd.gate = lwd.gate, lty.gate = lty.gate, pts.gate = pts.gate, pch.gate = pch.gate, cex.gate = cex.gate)
    
  }, gts, col.gate, lty.gate, lwd.gate, pch.gate, cex.gate)
  
  
})

#' Plot filters List of Gate Objects onto an Existing Plot
#'
#' @param x an object of class \code{\link[flowCore:filters-class]{filters}}
#'   containing objects of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param col.gate indicates the colour of the gate to be constructed, set to
#'   \code{"red"} by default.
#' @param lwd.gate numeric to adjust line thickness of gates, set to \code{2.5}
#'   by default.
#' @param lty.gate integer [0,6] which controls the line type, set to \code{1}
#'   to draw solid lines by default.
#' @param pts.gate logical indicating whether points should be included when
#'   plotting the gates, set to \code{FALSE} by default.
#' @param pch.gate integer [0,25] passed to pch to control the shape of the
#'   points, set to \code{16} to draw filled circles by default. For other
#'   shapes refer to \code{\link[graphics:points]{?pch}}.
#' @param cex.gate numeric character expansion to control the size of the points
#'   in the drawn gate, set to \code{1} by default.
#'
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{plotGates,rectangleGate-method}}
#' @seealso \code{\link{plotGates,polygonGate-method}}
#' @seealso \code{\link{plotGates,ellipsoidGate-method}}
#' @seealso \code{\link{plotGates,list-method}}
#'
#' @export
setMethod(plotGates, signature = "filters", definition = function(x, channels, col.gate = "red", lwd.gate = 2.5, lty.gate = 1, pts.gate = FALSE, pch.gate = 16, cex.gate = 1){
  
  # Assign x to gts
  gts <- x
  
  # Check Channels
  if(missing(channels)){
    
    channels <- as.vector(parameters(gts[[1]]))
    
  }
  
  if(!all(channels %in% as.vector(parameters(gts[[1]])))){
    
    stop("Channels used to construct the plot do not match those of the supplied gate.")
    
  }
  
  # Gate colours
  if(length(col.gate) != length(gts)){
    
    if(length(col.gate) == 1){
      
      col.gate <- rep(col.gate, length(gts))
      
    }else if(length(col.gate) > length(gts)){
      
      col.gate <- col.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Line Type
  if(length(lty.gate) != length(gts)){
    
    if(length(lty.gate) == 1){
      
      lty.gate <- rep(lty.gate, length(gts))
      
    }else if(length(lty.gate) > length(gts)){
      
      lty.gate <- lty.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Line Width
  if(length(lwd.gate) != length(gts)){
    
    if(length(lwd.gate) == 1){
      
      lwd.gate <- rep(lwd.gate, length(gts))
      
    }else if(length(lwd.gate) > length(gts)){
      
      lwd.gate <- lwd.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Point Character
  if(length(pch.gate) != length(gts)){
    
    if(length(pch.gate) == 1){
      
      pch.gate <- rep(pch.gate, length(gts))
      
    }else if(length(pch.gate) > length(gts)){
      
      pch.gate <- pch.gate[1:length(gts)]
      
    }
    
  }
  
  # Gate Point Expansion
  if(length(cex.gate) != length(gts)){
    
    if(length(cex.gate) == 1){
      
      cex.gate <- rep(cex.gate, length(gts))
      
    }else if(length(cex.gate) > length(gts)){
      
      cex.gate <- cex.gate[1:length(gts)]
      
    }
    
  }
  
  # Plot Gates
  mapply(function(gt, col.gate, lty.gate, lwd.gate, pch.gate, cex.gate){
    
    plotGates(gt, channels = channels, col.gate = col.gate, lwd.gate = lwd.gate, lty.gate = lty.gate, pts.gate = pts.gate, pch.gate = pch.gate, cex.gate = cex.gate)
    
  }, gts, col.gate, lty.gate, lwd.gate, pch.gate, cex.gate)
  
  
})
