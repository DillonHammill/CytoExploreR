#' Plot Labels for Gated Populations
#'
#' @param x object of class \code{"flowFrame"}.
#' @param gates object of class
#'   \code{\link[flowCore:rectangleGate]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}, \code{"list"} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param ... additional method-specific arguments.
#'
#' @return add label to existing plot with information about gated population.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom flowCore Subset parameters
#'
#' @seealso \code{\link{plotLabels,flowFrame,rectangleGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,polygonGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,ellipsoidGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,list-method}}
#' @seealso \code{\link{plotLabels,flowFrame,filters-method}}
#'
#' @export
setGeneric(name = "plotLabels",
           def = function(x, gates, ...){standardGeneric("plotLabels")}
)

#' Plot Labels for Gated Populations - rectangleGate Method
#'
#' \code{plotLabels} takes on a \code{flowFrame} object, population name
#' \code{alias}, \code{channels} and a gate object to construct a text label for
#' the plot with the population name and frequency.
#'
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}} gated in the
#'   existing plot.
#' @param gates an object of class \code{\link[flowCore:rectangleGate]{rectangleGate}}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param alias the name of the gated population.
#' @param format.text indicates the type of text to include in the label, can be
#'   either \code{"alias"}, \code{"percent"}, \code{"count"},
#'   \code{c("alias","percent")} or \code{c("alias","count")}. Set to
#'   \code{c("alias","percent")} by default.
#' @param x.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param y.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param font.text integer [1,2,3,4] passed to \code{text} to alter the font,
#'   set to \code{2} by default for a bold font.
#' @param col.text specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param cex.text numeric character expansion used to control the size of the
#'   text in the labels, set to \code{0.8} by default. See \code{?text} for
#'   details.
#' @param alpha numeric [0,1] controls the transparency of the background, set
#'   to \code{0.6} by default.
#'
#' @return add label to existing plot with information about gated population.
#'
#' @importFrom flowCore Subset parameters
#' @importFrom graphics par
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotLabels,flowFrame,polygonGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,ellipsoidGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,list-method}}
#' @seealso \code{\link{plotLabels,flowFrame,filters-method}}
#'
#' @export
setMethod(plotLabels, signature = c("flowFrame", "rectangleGate"), 
          definition = function(x, gates, channels, alias = NULL, format.text = NULL, x.text = NULL, y.text = NULL, font.text = 2, col.text = "black", cex.text = 0.8, alpha = 0.6){
 
  # Assign x to fr
  fr <- x
            
  # Channels needed to position label
  if(missing(channels)){
    
    stop("Supply the name(s) of the fluorescent channels used to construct the plot.")
    
  }
  
  # Alias
  if(is.null(alias) & is.null(format.text)){
    
    # No population name supplied only show percent on label
    format.text <- "percent"
    
  }else if(is.null(alias) & !is.null(format.text)){
    
    if(format.text == c("alias","percent") | format.text == c("alias","count")){
      
      message(paste("No population named supplied to alias - labels will show", format.text[2], "only."))
      format.text <- format.text[2]
      
    }
    
  }

  # Total events
  events <- nrow(fr)
            
  # Population events - percentages
  cnt <- BiocGenerics::nrow(flowCore::Subset(fr, gates))
  prcnt <- round(cnt/events,4)
  prcnt <- sprintf("%.2f %%",100*prcnt)
  
  # Check supplied channels & gate channels
  chans <- parameters(gates)  
  
  # 1D gate supplied
  if(length(chans) == 1){
    
    if(length(channels) == 1){
      
      if(!channels == chans){
        
        stop("Supplied channel does not match that of the supplied gate.")
        
      }else if(channels == chans){
        
        if(is.null(x.text) & is.null(y.text)){
          
          xmin <- gates@min[channels[1]]
          xmax <- gates@max[channels[1]]
          
          if(is.infinite(xmin)){
            
            xmin <- par("usr")[1]
            
          }
          
          if(is.infinite(xmax)){
            
            xmax <- par("usr")[2]
            
          }
          
          ymin <- par("usr")[3]
          ymax <- par("usr")[4]
          
        }
        
      }
      
    }else if(length(channels) == 2){
      
      if(!chans %in% channels){
        
        stop("Supplied channels do not match that of the supplied gate.")
        
      }else if(chans %in% channels){
        
        if(is.null(x.text) & is.null(y.text)){
          
          ind <- match(chans, channels) 
          
          if(ind == 1){
            
            xmin <- gates@min[channels[1]]
            xmax <- gates@max[channels[1]]
          
          if(is.infinite(xmin)){
            
            xmin <- par("usr")[1]
            
          }
          
          if(is.infinite(xmax)){
            
            xmax <- par("usr")[2]
            
          }
          
          ymin <- par("usr")[3]
          ymax <- par("usr")[4]
          
          }else if(ind == 2){
          
            ymin <- gates@min[channels[2]]
            ymax <- gates@max[channels[2]]
            
            if(is.infinite(xmin)){
              
              ymin <- par("usr")[3]
              
            }
            
            if(is.infinite(xmax)){
              
              ymax <- par("usr")[4]
              
            }
            
            xmin <- par("usr")[1]
            xmax <- par("usr")[2]
            
        }
        
      }
        
     }
      
   }
    
  # 2D gate supplied
  }else if(length(chans) == 2){
    
    if(!all(chans %in% channels)){
      
      stop("Supplied channels do not match that of the supplied gate.")
      
    }else if(all(chans %in% channels)){
      
      if(is.null(x.text) & is.null(y.text)){
        
        xmin <- gates@min[channels[1]]
        ymin <- gates@min[channels[2]]
        xmax <- gates@max[channels[1]]
        ymax <- gates@max[channels[2]]
    
        if(is.infinite(xmin)){
      
          xmin <- par("usr")[1]
      
        }
    
        if(is.infinite(ymin)){
      
          ymin <- par("usr")[3]
      
        }
    
        if(is.infinite(xmax)){
      
           xmax <- par("usr")[2]
      
        }
    
        if(is.infinite(ymax)){
      
          ymax <- par("usr")[4]
      
        }
        
      }
      
    }
    
  } 
  
  # Label position
  if(is.null(x.text) & is.null(y.text)){
    
    x.text <- (xmin + xmax)/2
    y.text <- (ymin + ymax)/2
    
  }
   
  # Add labels
  if(all(format.text %in% "alias")){
    
    boxed.labels(x = x.text, y = y.text, labels = alias, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
    
  }else if(all(format.text %in% "percent")){
    
    boxed.labels(x = x.text, y = y.text, labels = prcnt, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
    
  }else if(all(format.text %in% "count")){
    
    boxed.labels(x = x.text, y = y.text, labels = cnt, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
    
  }else if(all(format.text %in% c("alias","percent"))){
    
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, prcnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
    
  }else if(all(format.text %in% c("alias","count"))){
    
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, cnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
    
  }else{
    
    message("Supplied format.text is not valid, reverting to default format.")
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, prcnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
    
  }               
                       
})

#' Plot Labels for Gated Populations - polygonGate Method
#'
#' \code{plotLabels} takes on a \code{flowFrame} object, population name
#' \code{alias}, \code{channels} and a gate object to construct a text label for
#' the plot with the population name and frequency.
#'
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}} gated in the
#'   existing plot.
#' @param gates an object of class
#'   \code{\link[flowCore:polygonGate]{polygonGate}}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param alias the name of the gated population.
#' @param format.text indicates the type of text to include in the label, can be
#'   either \code{"alias"}, \code{"percent"}, \code{"count"},
#'   \code{c("alias","percent")} or \code{c("alias","count")}. Set to
#'   \code{c("alias","percent")} by default.
#' @param x.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param y.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param font.text integer [1,2,3,4] passed to \code{text} to alter the font,
#'   set to \code{2} by default for a bold font.
#' @param col.text specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param cex.text numeric character expansion used to control the size of the
#'   text in the labels, set to \code{0.8} by default. See \code{?text} for
#'   details.
#' @param alpha numeric [0,1] controls the transparency of the background, set
#'   to \code{0.6} by default.
#'
#' @return add label to existing plot with information about gated population.
#'
#' @importFrom flowCore Subset parameters
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotLabels,flowFrame,rectangleGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,ellipsoidGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,list-method}}
#' @seealso \code{\link{plotLabels,flowFrame,filters-method}}
#'
#' @export
setMethod(plotLabels, signature = c("flowFrame", "polygonGate"), 
          definition = function(x, gates, channels, alias = NULL, format.text = NULL, x.text = NULL, y.text = NULL, font.text = 2, col.text = "black", cex.text = 0.8, alpha = 0.6){
            
  # Assign x to fr
  fr <- x
            
  # Channels needed to position label
  if(missing(channels)){
              
    stop("Supply the name(s) of the fluorescent channels used to construct the plot.")
              
  }
            
  # Alias
  if(is.null(alias) & is.null(format.text)){
              
    # No population name supplied only show percent on label
    format.text <- "percent"
              
  }else if(is.null(alias) & !is.null(format.text)){
              
    if(format.text == c("alias","percent") | format.text == c("alias","count")){
                
      message(paste("No population named supplied to alias - labels will show", format.text[2], "only."))
      format.text <- format.text[2]
                
    }
              
  }
            
  # Total events
  events <- nrow(fr)
            
  # Population events - percentages
  cnt <- BiocGenerics::nrow(flowCore::Subset(fr, gates))
  prcnt <- round(cnt/events,4)
  prcnt <- sprintf("%.2f %%",100*prcnt)
  
  # Check supplied channels & gate channels
  chans <- parameters(gates)  
            
  if(!all(chans %in% channels)){
    
    stop("Supplied channels do not match that of the supplied gate.")
    
  }else{
    
    if(is.null(x.text) & is.null(y.text)){
      
      x.text <- sum(gates@boundaries[, channels[1]])/nrow(gates@boundaries)
      y.text <- sum(gates@boundaries[, channels[2]])/nrow(gates@boundaries)
      
    }

  }
            
  # Add labels
  if(all(format.text %in% "alias")){
              
    boxed.labels(x = x.text, y = y.text, labels = alias, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% "percent")){
              
    boxed.labels(x = x.text, y = y.text, labels = prcnt, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% "count")){
              
    boxed.labels(x = x.text, y = y.text, labels = cnt, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% c("alias","percent"))){
              
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, prcnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% c("alias","count"))){
              
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, cnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else{
              
    message("Supplied format.text is not valid, reverting to default format.")
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, prcnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }               
            
})

#' Plot Labels for Gated Populations - ellipsoidGate Method
#'
#' \code{plotLabels} takes on a \code{flowFrame} object, population name
#' \code{alias}, \code{channels} and a gate object to construct a text label for
#' the plot with the population name and frequency.
#'
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}} gated in the
#'   existing plot.
#' @param gates an object of class
#'   \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param alias the name of the gated population.
#' @param format.text indicates the type of text to include in the label, can be
#'   either \code{"alias"}, \code{"percent"}, \code{"count"},
#'   \code{c("alias","percent")} or \code{c("alias","count")}. Set to
#'   \code{c("alias","percent")} by default.
#' @param x.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param y.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param font.text integer [1,2,3,4] passed to \code{text} to alter the font,
#'   set to \code{2} by default for a bold font.
#' @param col.text specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param cex.text numeric character expansion used to control the size of the
#'   text in the labels, set to \code{0.8} by default. See \code{?text} for
#'   details.
#' @param alpha numeric [0,1] controls the transparency of the background, set
#'   to \code{0.6} by default.
#'
#' @return add label to existing plot with information about gated population.
#'
#' @importFrom flowCore Subset parameters
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotLabels,flowFrame,rectangleGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,polygonGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,list-method}}
#' @seealso \code{\link{plotLabels,flowFrame,filters-method}}
#'
#' @export
setMethod(plotLabels, signature = c("flowFrame", "ellipsoidGate"), 
          definition = function(x, gates, channels, alias = NULL, format.text = NULL, x.text = NULL, y.text = NULL, font.text = 2, col.text = "black", cex.text = 0.8, alpha = 0.6){
            
  # Assign x to fr
  fr <- x
            
  # Channels needed to position label
  if(missing(channels)){
              
    stop("Supply the name(s) of the fluorescent channels used to construct the plot.")
              
  }
            
  # Alias
  if(is.null(alias) & is.null(format.text)){
              
    # No population name supplied only show percent on label
    format.text <- "percent"
              
  }else if(is.null(alias) & !is.null(format.text)){
              
    if(format.text == c("alias","percent") | format.text == c("alias","count")){
                
      message(paste("No population named supplied to alias - labels will show", format.text[2], "only."))
      format.text <- format.text[2]
                
    }
              
  }
            
  # Total events
  events <- nrow(fr)
            
  # Population events - percentages
  cnt <- BiocGenerics::nrow(flowCore::Subset(fr, gates))
  prcnt <- round(cnt/events,4)
  prcnt <- sprintf("%.2f %%",100*prcnt)
            
  # Check supplied channels & gate channels
  chans <- parameters(gates)  
            
  if(!all(chans %in% channels)){
              
    stop("Supplied channels do not match that of the supplied gate.")
              
  }else{
              
    if(is.null(x.text) & is.null(y.text)){
                
      x.text <- gates@mean[channels[1]]
      y.text <- gates@mean[channels[2]] 
                
    }
              
  }
            
  # Add labels
  if(all(format.text %in% "alias")){
              
    boxed.labels(x = x.text, y = y.text, labels = alias, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% "percent")){
              
    boxed.labels(x = x.text, y = y.text, labels = prcnt, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% "count")){
              
    boxed.labels(x = x.text, y = y.text, labels = cnt, border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% c("alias","percent"))){
              
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, prcnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else if(all(format.text %in% c("alias","count"))){
              
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, cnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }else{
              
    message("Supplied format.text is not valid, reverting to default format.")
    boxed.labels(x = x.text, y = y.text, labels = paste(alias, prcnt, sep = "\n"), border = FALSE, font = font.text, col = col.text, alpha.bg = alpha, cex = cex.text)
              
  }               
            
})

#' Plot Labels for Gated Populations - list Method
#'
#' \code{plotLabels} takes on a \code{flowFrame} object, population name
#' \code{alias}, \code{channels} and a gate object to construct a text label for
#' the plot with the population name and frequency.
#'
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}} gated in the
#'   existing plot.
#' @param gates an object of class \code{"list"} containing objects of class
#'   \code{\link[flowCore:rectangleGate]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate]{polygonGate}} or
#'   \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param alias the name of the gated population.
#' @param format.text indicates the type of text to include in the label, can be
#'   either \code{"alias"}, \code{"percent"}, \code{"count"},
#'   \code{c("alias","percent")} or \code{c("alias","count")}. Set to
#'   \code{c("alias","percent")} by default.
#' @param x.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param y.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param font.text integer [1,2,3,4] passed to \code{text} to alter the font,
#'   set to \code{2} by default for a bold font.
#' @param col.text specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param cex.text numeric character expansion used to control the size of the
#'   text in the labels, set to \code{0.8} by default. See \code{?text} for
#'   details.
#' @param alpha numeric [0,1] controls the transparency of the background, set
#'   to \code{0.6} by default.
#'
#' @return add label to existing plot with information about gated population.
#'
#' @importFrom flowCore Subset parameters
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotLabels,flowFrame,rectangleGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,polygonGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,ellipsoidGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,filters-method}}
#'
#' @export
setMethod(plotLabels, signature = c("flowFrame", "list"), 
          definition = function(x, gates, channels, alias = NULL, format.text = NULL, x.text = NULL, y.text = NULL, font.text = 2, col.text = "black", cex.text = 0.8, alpha = 0.6){
            
            
  # Make calls to plotLabels
  mapply(function(gate, alias, font.text, col.text, cex.text, alpha){
            
    plotLabels(x = x, gates = gate, channels = channels, alias = alias, format.text = format.text, font.text = font.text, col.text = col.text, cex.text = cex.text, alpha = alpha)
    
  }, gates, alias, font.text, col.text, cex.text, alpha)
            
})

#' Plot Labels for Gated Populations - filters Method
#'
#' \code{plotLabels} takes on a \code{flowFrame} object, population name
#' \code{alias}, \code{channels} and a gate object to construct a text label for
#' the plot with the population name and frequency.
#'
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}} gated in the
#'   existing plot.
#' @param gates an object of class \code{\link[flowCore:filters-class]{filters}}
#'   containing objects of class
#'   \code{\link[flowCore:rectangleGate]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate]{polygonGate}} or
#'   \code{\link[flowCore:ellipsoidGate]{ellipsoidGate}}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param alias the name of the gated population.
#' @param format.text indicates the type of text to include in the label, can be
#'   either \code{"alias"}, \code{"percent"}, \code{"count"},
#'   \code{c("alias","percent")} or \code{c("alias","count")}. Set to
#'   \code{c("alias","percent")} by default.
#' @param x.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param y.text vector containing the x co-ordiantes for the plot labels. Set
#'   to \code{NULL} by default to place labels in the center of the gates.
#' @param font.text integer [1,2,3,4] passed to \code{text} to alter the font,
#'   set to \code{2} by default for a bold font.
#' @param col.text specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param cex.text numeric character expansion used to control the size of the
#'   text in the labels, set to \code{0.8} by default. See \code{?text} for
#'   details.
#' @param alpha numeric [0,1] controls the transparency of the background, set
#'   to \code{0.6} by default.
#'
#' @return add label to existing plot with information about gated population.
#'
#' @importFrom flowCore Subset parameters
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{plotLabels,flowFrame,rectangleGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,polygonGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,ellipsoidGate-method}}
#' @seealso \code{\link{plotLabels,flowFrame,list-method}}
#'
#' @export
setMethod(plotLabels, signature = c("flowFrame", "filters"), 
          definition = function(x, gates, channels, alias = NULL, format.text = NULL, x.text = NULL, y.text = NULL, font.text = 2, col.text = "black", cex.text = 0.8, alpha = 0.6){
            
            
  # Make calls to plotLabels
  mapply(function(gate, alias, font.text, col.text, cex.text, alpha){
              
    plotLabels(x = x, gates = gate, channels = channels, alias = alias, format.text = format.text, font.text = font.text, col.text = col.text, cex.text = cex.text, alpha = alpha)
              
  }, gates, alias, font.text, col.text, cex.text, alpha)
            
})
