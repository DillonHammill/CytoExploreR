#' Plot Flow Cytometry Gating Strategies
#'
#' \code{plotCytoGates} automatically plots the entire gating strategy and has
#' full support for back-gating through \code{overlay}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional method specific arguments.
#' 
#' @seealso \code{\link{plotCytoGates,GatingHierarchy-method}}
#' @seealso \code{\link{plotCytoGates,GatingSet-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(name = "plotCytoGates",
           def = function(x, ...){standardGeneric("plotCytoGates")}
)


#' Plot Flow Cytometry Gating Strategies - GatingSet Method
#'
#' \code{plotCytoGates} automatically plots the entire gating strategy and has
#' full support for back-gating through \code{overlay}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param gtfile name of the gatingTemplate csv file used to gate \code{x}. If
#'   not supplied the gating strategy will be obtained directly from the
#'   GatingSet.
#' @param merge logical indicating whether the flowFrames should be merged prior
#'   to plotting to yield a single plot, set to FALSE by default.
#' @param overlay name(s) of population(s) to overlay onto gating strategy for
#'   back-gating.
#' @param title character string to use as the title for the plot layout, set to
#'   "Gating Strategy" by default.
#' @param main vector of titles to use above each plot.
#' @param mfrow a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param popup logical indicating whether the gating strategy should be plotted
#'   in a pop-up window, set to FALSE by default.
#' @param col colour of points in 2D plots set to NA to use default red-blue
#'   colour scale. Control the fill colour of overlays by supplying multiple
#'   colours to this argument (e.g. c("blue","red")).
#' @param fill fill colour for 1D density distributions. Control the colour of
#'   overlays by supplying multiple colours to this argument (e.g. c(NA,"red")).
#' @param legend logical indicating whether a legend should be included when an
#'   overlay is supplied.
#' @param text.legend vector of character strings to use for legend when an
#'   overlay is supplied.
#' @param cex.legend character expansion for legend text, set to 1.2 by default.
#' @param ... extra arguments passed to plotCyto, see
#'   \code{\link{plotCyto1d,flowFrame-method}} and
#'   \code{\link{plotCyto2d,flowFrame-method}} for more details.
#'
#' @importFrom openCyto templateGen
#' @importFrom flowWorkspace getGate
#' @importFrom grDevices colorRampPalette n2mfrow
#' @importFrom graphics mtext legend plot.new
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{plotCytoGates,GatingHierarchy-method}}
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setMethod(plotCytoGates, signature = "GatingSet", 
          definition = function(x, gtfile = NULL, merge = TRUE, overlay = NULL, title = "Gating Strategy", main, mfrow, popup = FALSE, col = NA, fill = "blue", legend = TRUE, text.legend, cex.legend = 1.2, ...){
  
  # Assign x to gs
  gs <- x
          
  # Merge
  if(length(gs) == 1){
    
    merge <- FALSE
    
  }
    
  # Layout title 
  title <- title
            
  # Overlay colours
  cols <- colorRampPalette(c("magenta","hotpink","brown","chocolate","orange","gold","yellow","springgreen","forestgreen","darkmagenta","darkturquoise","cyan"))
  
  # Colours
  if(length(col) != 1 & length(col) != (length(overlay) + 1)){
    
    stop("Please supply colours for the parent and overlayed populations.")
    
  }
  
  # Use GatingSet directly to get gating strategy
  if(is.null(gtfile)){
    
    gt <- templateGen(gs[[1]])
    gts <- data.frame(parent = basename(gt$parent), alias = gt$alias, xchannel = do.call(rbind, strsplit(gt$dims,',', fixed = TRUE))[,1], ychannel = do.call(rbind, strsplit(gt$dims,',', fixed = TRUE))[,2], stringsAsFactors = FALSE)
    
    # Extract Gates
    gates <- lapply(gts$alias, function(pop){
      
      getGate(gs[[1]], pop)
      
    })
    names(gates) <- gts$alias
    
    # Extract unique parents for plotting
    parents <- c()
    lapply(unique(as.vector(gts$parent)), function(parent){
      
      parents <<- c(parents, rep(parent, nrow(unique(gts[gts$parent == parent, c("xchannel","ychannel")]))))
      
    })
    
    # Pop-up window?
    if(popup == TRUE){
      
      checkOSGD()
      
    }
    
    # Calculate mfrow parameters based on number of parents
    if(missing(mfrow)){
      
      mfrow <- c((n2mfrow(length(parents) + 1))[2], n2mfrow((length(parents) + 1))[1])
      par(mfrow = mfrow)
      
    }
    
    if(!is.null(title)){
      
      par(oma = c(0,0,3,0))
      
    }
    
    # Titles
    if(missing(main)){
      
      prnts <- parents
      if("root" %in% prnts){
        
        prnts[prnts %in% "root"] <- "All Events"
        
      }
      main <- prnts
      
    }
    
    mapply(function(parents, main){
      
      parent <- as.character(parents)
      alias <- as.vector(gts[gts$parent == parents, "alias"])
      xchannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "xchannel"])
      ychannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "ychannel"])
      channels <- c(xchannel,ychannel)
      
      if(channels[1] == channels[2]){
        
        channels <- channels[1]
        
      }
      
      if(length(gs) == 1){
        
        if(!is.null(overlay)){
          
          if(length(channels) == 1 & length(fill) == 1){
            
            fill <<- c(fill, cols(length(overlay)))
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
            
          }else if(length(channels) == 1 & length(fill) == (length(overlay) + 1)){
            
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
            
          }else if(length(channels) == 2 & length(col) == 1){
            
            col <<- c(col, cols(length(overlay)))
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, col = col, ...)
            
          }else if(length(channels) == 2 & length(col) == (length(overlay) + 1)){
            
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, col = col, ...)
            
          }
          
        }else{
          
           plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
        
        }
        
      }else if(length(gs) > 1){
        
        if(!is.null(overlay)){
          
          if(length(channels) == 1 & length(fill) == 1){
            
            fill <<- c(fill, cols(length(overlay)))
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
           
          }else if(length(channels) == 1 & length(fill) == (length(overlay) + 1)){
            
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...) 
            
          }else if(length(channels) == 2 & length(col) == 1){
            
            col <<- c(col, cols(length(overlay)))
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, col = col, ...)
            
          }else if(length(channels) == 2 & length(col) == (length(overlay) + 1)){
            
            plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, col = col, ...)
            
          }
          
        }else{
          
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
          
        }
        
      }
      
    }, parents, main)
    
    # Title
    if(!is.null(title)){
      
      mtext(title, outer = TRUE, cex = 1, font = 2)
      
    }
    
    # Legend
    if(!is.null(overlay)){
      
      if(legend == TRUE){
        
        # Legend Text
        if(missing(text.legend)){
          
          if(class(overlay) == "character"){
            
            text.legend <- overlay
            
          }else{
            
            stop("Please supply vector of names to use in the legend.")
            
          }
          
        }
        
        # Add dummy plot
        plot.new()
      
        # Add legend
        legend(x = 0.1, y = 0.8, legend = text.legend, fill = col[2:length(col)], bty = "n", cex = cex.legend, x.intersp = 0.5)
      
      }
      
    }
   
  # Extract gating strategy from gatingTemplate   
  }else if(!is.null(gtfile)){
    
    if(checkFile(gtfile) == FALSE){
      
      stop("Supplied gatingTemplate file does not exist in current working directory.")
      
    }
    
    # Use gtfile to get gating strategy
    gt <- read.csv(gtfile, header = TRUE)
    pops <- gt$alias
    parents <- gt$parent
  
    gates <- mapply(function(pops, parents){
      
      extractGate(gtfile = gtfile, parent = as.character(parents), alias = as.character(pops))  
      
    }, pops, parents)
    names(gates) <- as.character(pops)
    
    channels <- sapply(gates, function(gate){
      
         channels <- as.vector(parameters(gate))

         if(length(channels) == 1){
           
           channels <- c(channels[1],NA)
           
         }
         
         return(channels)
      
      })
    gts <- data.frame(parent = parents, alias = pops, xchannel = channels[1,], ychannel = channels[2,])
  
  # Extract unique parents for plotting
  parents <- c()
  lapply(unique(as.vector(gts$parent)), function(parent){
    
    parents <<- c(parents, rep(parent, nrow(unique(gts[gts$parent == parent, c("xchannel","ychannel")]))))
    
  })
  
  # Pop-up window?
  if(popup == TRUE){
    
    checkOSGD()
    
  }
  
  # Calculate mfrow parameters based on number of parents
  if(missing(mfrow)){
    
    mfrow <- c((n2mfrow(length(parents) + 1))[2], n2mfrow((length(parents) + 1))[1])
    par(mfrow = mfrow)
    
  }
  
  if(!is.null(title)){
    
    par(oma = c(0,0,3,0))
    
  }
  
  # Titles
  if(missing(main)){
    
    prnts <- parents
    if("root" %in% prnts){
      
      prnts[prnts %in% "root"] <- "All Events"
      
    }
    main <- prnts
    
  }
  
  mapply(function(parents, main){
    
    parent <- as.character(parents)
    alias <- as.vector(gts[gts$parent == parents, "alias"])
    xchannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "xchannel"])
    ychannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "ychannel"])
    channels <- c(xchannel,ychannel)
    
    if(channels[1] == channels[2]){
      
      channels <- channels[1]
      
    }
    
    if(length(gs) == 1){
      
      if(!is.null(overlay)){
        
        if(length(channels) == 1 & length(fill) == 1){
          
          fill <<- c(fill, cols(length(overlay)))
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = FALSE, main = main, text.labels = alias, fill = fill, ...)
         
        }else if(length(channels) == 1 & length(fill) == (length(overlay) + 1)){
           
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = FALSE, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 2 & length(col) == 1){
          
          col <<- c(col, cols(length(overlay)))
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = FALSE, main = main, text.labels = alias, col = col, ...)
          
        }else if(length(channels) == 2 & length(col) == (length(overlay) + 1)){
          
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = FALSE, main = main, text.labels = alias, col = col, ...)
          
        }
        
      }else{
        
        plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = FALSE, main = main, text.labels = alias, fill = fill, ...)
        
      }
      
    }else if(length(gs) > 1){
      
      if(!is.null(overlay)){
        
        if(length(channels) == 1 & length(fill) == 1){
          
          fill <<- c(fill, cols(length(overlay)))
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 1 & length(fill) == (length(overlay) + 1)){
          
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 2 & length(col) == 1){
          
          col <<- c(col, cols(length(overlay)))
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, col = col, ...)
          
        }else if(length(channels) == 2 & length(col) == (length(overlay) + 1)){
          
          plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, col = col, ...)
          
        }
        
      }else{
        
        plotCyto(gs, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, merge = merge, main = main, text.labels = alias, fill = fill, ...)
        
      }
      
    }
    
  }, parents, main)
  
  # Title
  if(!is.null(title)){
    
    mtext(title, outer = TRUE, cex = 1, font = 2)
    
  }
  
  # Legend
  if(!is.null(overlay)){
    
    if(legend == TRUE){
      
      # Legend Text
      if(missing(text.legend)){
        
        if(class(overlay) == "character"){
          
          text.legend <- overlay
          
        }else{
          
          stop("Please supply vector of names to use in the legend.")
          
        }
        
      }
      
      # Add dummy plot
      plot.new()
      
      # Add legend
      legend(x = 0.1, y = 0.8, legend = text.legend, fill = col[2:length(col)], bty = "n", cex = cex.legend, x.intersp = 0.5)
      
    }
    
  }
  
  }
  
  # Return default plot layout
  par(mfrow = c(1,1))
  par(oma = c(0,0,0,0))
  
})

#' Plot Flow Cytometry Gating Strategies - GatingHierarchy Method
#'
#' \code{plotCytoGates} automatically plots the entire gating strategy and has
#' full support for back-gating through \code{overlay}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}}.
#' @param gtfile name of the gatingTemplate csv file used to gate \code{x}. If
#'   not supplied the gating strategy will be obtained directly from the
#'   GatingSet.
#' @param overlay name(s) of population(s) to overlay onto gating strategy for
#'   back-gating.
#' @param title character string to use as the title for the plot layout, set to
#'   "Gating Strategy" by default.
#' @param main vector of titles to use above each plot.
#' @param mfrow a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param popup logical indicating whether the gating strategy should be plotted
#'   in a pop-up window, set to FALSE by default.
#' @param col colour of points in 2D plots set to NA to use default red-blue
#'   colour scale. Control the fill colour of overlays by supplying multiple
#'   colours to this argument (e.g. c("blue","red")).
#' @param fill fill colour for 1D density distributions. Control the colour of
#'   overlays by supplying multiple colours to this argument (e.g. c(NA,"red")).
#' @param legend logical indicating whether a legend should be included when an
#'   overlay is supplied.
#' @param text.legend vector of character strings to use for legend when an
#'   overlay is supplied.
#' @param cex.legend character expansion for legend text, set to 1.2 by default.
#' @param ... extra arguments passed to plotCyto, see
#'   \code{\link{plotCyto1d,flowFrame-method}} and
#'   \code{\link{plotCyto2d,flowFrame-method}} for more details.
#'
#' @importFrom openCyto templateGen
#' @importFrom flowWorkspace getGate
#' @importFrom grDevices colorRampPalette n2mfrow
#' @importFrom graphics mtext legend plot.new
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{plotCytoGates,GatingHierarchy-method}}
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setMethod(plotCytoGates, signature = "GatingHierarchy", 
          definition = function(x, gtfile = NULL, overlay = NULL, title = "Gating Strategy", main, popup = FALSE, mfrow, col = NA, fill = "blue", legend = TRUE, text.legend, cex.legend = 1.2, ...){
            
  # Assign x to gh
  gh <- x
            
  # Layout title 
  title <- title
            
  # Overlay colours
  cols <- colorRampPalette(c("magenta","hotpink","brown","chocolate","orange","gold","yellow","springgreen","forestgreen","darkmagenta","darkturquoise","cyan"))
  
  # Colours
  if(length(col) != 1 & length(col) != (length(overlay) + 1)){
    
    stop("Please supply colours for the parent and overlayed populations.")
    
  }
  
  # Use GatingHierarchy directly to get gating strategy
  if(is.null(gtfile)){
    
    gt <- templateGen(gh)
    gts <- data.frame(parent = basename(gt$parent), alias = gt$alias, xchannel = do.call(rbind, strsplit(gt$dims,',', fixed = TRUE))[,1], ychannel = do.call(rbind, strsplit(gt$dims,',', fixed = TRUE))[,2], stringsAsFactors = FALSE)
    
    # Extract Gates
    gates <- lapply(gts$alias, function(pop){
      
      getGate(gh, pop)
      
    })
    names(gates) <- gts$alias
    
    # Extract unique parents for plotting
    parents <- c()
    lapply(unique(as.vector(gts$parent)), function(parent){
      
      parents <<- c(parents, rep(parent, nrow(unique(gts[gts$parent == parent, c("xchannel","ychannel")]))))
      
    })
    
    # Pop-up window?
    if(popup == TRUE){
      
      checkOSGD()
      
    }
    
    # Calculate mfrow parameters based on number of parents
    if(missing(mfrow)){
      
      mfrow <- c((n2mfrow(length(parents) + 1))[2], n2mfrow((length(parents) + 1))[1])
      par(mfrow = mfrow)
      
    }
    
    if(!is.null(title)){
      
      par(oma = c(0,0,3,0))
      
    }
    
    # Titles
    if(missing(main)){
      
      prnts <- parents
      if("root" %in% prnts){
        
        prnts[prnts %in% "root"] <- "All Events"
        
      }
      main <- prnts
      
    }
    
    mapply(function(parents, main){
      
      parent <- as.character(parents)
      alias <- as.vector(gts[gts$parent == parents, "alias"])
      xchannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "xchannel"])
      ychannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "ychannel"])
      channels <- c(xchannel,ychannel)
        
      if(channels[1] == channels[2]){
        
        channels <- channels[1]
        
      }
      
      
      if(!is.null(overlay)){
        
        if(length(channels) == 1 & length(fill) == 1){
          
          fill <<- c(fill, cols(length(overlay)))
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 1 & length(fill) == (length(overlay) + 1)){
          
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 2 & length(col) == 1){
          
          col <<- c(col, cols(length(overlay)))
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, col = col, ...)
          
        }else if(length(channels) == 2 & length(col) == (length(overlay) + 1)){
          
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, col = col, ...)
          
        }
        
      }else{
        
        plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, fill = fill, ...)
        
      }

    }, parents, main)
    
    # Title
    if(!is.null(title)){
      
      mtext(title, outer = TRUE, cex = 1, font = 2)
      
    }
    
    # Legend
    if(!is.null(overlay)){
      
      if(legend == TRUE){
        
        # Legend Text
        if(missing(text.legend)){
          
          if(class(overlay) == "character"){
            
            text.legend <- overlay
            
          }else{
            
            stop("Please supply vector of names to use in the legend.")
            
          }
          
        }
        
        # Add dummy plot
        plot.new()
        
        # Add legend
        legend(x = 0.1, y = 0.8, legend = text.legend, fill = col[2:length(col)], bty = "n", cex = cex.legend, x.intersp = 0.5)
        
      }
      
    }
    
  # Extract gating strategy from gatingTemplate   
  }else if(!is.null(gtfile)){
              
    if(checkFile(gtfile) == FALSE){
                
      stop("Supplied gatingTemplate file does not exist in current working directory.")
                
    }
              
    # Use gtfile to get gating strategy
    gt <- read.csv(gtfile, header = TRUE)
    pops <- gt$alias
    parents <- gt$parent
              
    gates <- mapply(function(pops, parents){
                
      extractGate(gtfile = gtfile, parent = as.character(parents), alias = as.character(pops))  
                
    }, pops, parents)
    names(gates) <- as.character(pops)
              
    channels <- sapply(gates, function(gate){
      
      channels <- as.vector(parameters(gate))
      
      if(length(channels) == 1){
        
        channels <- c(channels[1],NA)
        
      }
      
      return(channels)
      
    })
    gts <- data.frame(parent = parents, alias = pops, xchannel = channels[1,], ychannel = channels[2,])
            
    # Extract unique parents for plotting
    parents <- c()
    lapply(unique(as.vector(gts$parent)), function(parent){
              
      parents <<- c(parents, rep(parent, nrow(unique(gts[gts$parent == parent, c("xchannel","ychannel")]))))
              
    })
        
    # Pop-up window?
    if(popup == TRUE){
      
      checkOSGD()
      
    }
        
    # Calculate mfrow parameters based on number of parents
    if(missing(mfrow)){
              
      mfrow <- c((n2mfrow(length(parents) + 1))[2], n2mfrow((length(parents) + 1))[1])
      par(mfrow = mfrow)
              
    }
    
    if(!is.null(title)){
      
      par(oma = c(0,0,3,0))
      
    }
            
    # Titles
    if(missing(main)){
              
      prnts <- parents
      if("root" %in% prnts){
                
        prnts[prnts %in% "root"] <- "All Events"
                
      }
      main <- prnts
              
    }
            
    mapply(function(parents, main){
              
      parent <- as.character(parents)
      alias <- as.vector(gts[gts$parent == parents, "alias"])
      xchannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "xchannel"])
      ychannel <- as.character(gts[gts$parent == parent & gts$alias == alias[1], "ychannel"])
      channels <- c(xchannel,ychannel)
           
      if(is.na(channels[2])){
        
        channels <- channels[1]
        
      }
           
      if(!is.null(overlay)){
        
        if(length(channels) == 1 & length(fill) == 1){
          
          fill <<- c(fill, cols(length(overlay)))
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 1 & length(fill) == (length(overlay) + 1)){
          
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, fill = fill, ...)
          
        }else if(length(channels) == 2 & length(col) == 1){
          
          col <<- c(col, cols(length(overlay)))
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, col = col, ...)
          
        }else if(length(channels) == 2 & length(col) == (length(overlay) + 1)){
          
          plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, col = col, ...)
          
        }
        
      }else{
        
        plotCyto(gh, parent = parent, alias = alias, overlay = overlay, channels = channels, legend = FALSE, mfrow = FALSE, main = main, text.labels = alias, fill = fill, ...)
        
      }
      
    }, parents, main)
          
    # Title
    if(!is.null(title)){
      
      mtext(title, outer = TRUE, cex = 1, font = 2)
      
    }
            
    # Legend
    if(!is.null(overlay)){
      
      if(legend == TRUE){
        
        # Legend Text
        if(missing(text.legend)){
          
          if(class(overlay) == "character"){
            
            text.legend <- overlay
            
          }else{
            
            stop("Please supply vector of names to use in the legend.")
            
          }
          
        }
        
        # Add dummy plot
        plot.new()
        
        # Add legend
        legend(x = 0.1, y = 0.8, legend = text.legend, fill = col[2:length(col)], bty = "n", cex = cex.legend, x.intersp = 0.5)
        
      }
      
    }
            
  } 
  
    # Return default plot layout
    par(mfrow = c(1,1))
    par(oma = c(0,0,0,0))
    
})
