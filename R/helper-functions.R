#' Extract fluorescent parameters
#' 
#' @export
setGeneric(name="getChannels",
           def=function(x){standardGeneric("getChannels")}
)

#' Extract fluorescent parameters from a flowFrame
#' 
#' @param x object \code{flowFrame}
#' 
#' @return vector of fluorescent channels
#' 
#' @export
setMethod(getChannels, signature = "flowFrame", definition = function(x){
  colnames(x@description$SPILL)
})

#' Extract fluorescent parameters from a flowSet
#' 
#' @param x object \code{flowSet}
#' 
#' @return vector of fluorescent channels
#' 
#' @export
setMethod(getChannels, signature = "flowSet", definition = function(x){
colnames(x[[1]]@description$SPILL)
})

#' Extract fluorescent parameters from a GatingSet
#' 
#' @param x object \code{GatingSet}
#' 
#' @return vector of fluorescent channels
#' 
#' @export
setMethod(getChannels, signature = "GatingSet", definition = function(x){
  colnames(x@data[[1]]@description$SPILL)
})

#' Select flourescent channels for compensation controls
#' 
#' @param x object of class \code{flowSet} or \code{GatingSet} containing compensation controls
#' 
#' @return vector of channels in order of \code{flowSet} or \code{GatingSet} samples.
#' 
#' @export
selectChannels <- function(fs){
  
  opts <- c(getChannels(fs), "Unstained")
  
  # Print sample name and select channel
  channels <- opts[sapply(pData(fs)$name, function(x){
    message("Select a fluorescent channel for the following compensation control:")
    print(x)
    menu(choices = opts, graphics = TRUE)
  })]
  
  return(channels)
}

#' Sample a flowFrame
#' 
#' @param fr object of class \code{flowFrame}.
#' @param size numeric indicating the number of events to keep. If \code{size > nrow(fr)}, size is set to \code{nrow(fr)}.
#' 
#' @return \code{flowFrame} restricted to \code{size} events.
#' 
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter
#' @importFrom flowCore Subset
#' 
#' @export
sampleFrame <- function(fr, size = 250000){
  
  # Number of events
  events <- nrow(fr)
  
  if(events < size){
    
    size <- events
  
  }else{
    
  }

  smp <- sampleFilter(size = size)
  fr <- Subset(fr, smp)
  
  return(fr)
}

#' Select flowFrames based on pData
#' 
#' @param fs an object of class \code{flowSet}.
#' @param pData vector of form \code{c("column","row")} indicating the rows of \code{pData(fs)} to extract. 
#' 
#' @return an object of calss \code{flowSet} containing the extracted \code{flowFrame} objects.
#' 
#' @export
selectFrames <- function(fs, pData){
  
  # Extract pData from flowSet 
  pd <- pData(fs)
  
  # Find which rows to extract
  rows <- pd[[pData[1]]] == pData[2]
  
  # Get sampleNames of these rows
  sn <- rownames(pd)[rows]
  
  # Extract these samples by name from flowSet
  fs <- fs[sn]
  
  return(fs)
  
}

#' Remove gate and edit gatingTemplate .csv file
#' 
#' @param gs an object of class \code{GatingSet}.
#' @param alias name(s) of the gate to remove (e.g. "Single Cells").
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g. "gatingTemplate.csv").
#' @param children logical indicating whether descendant populations should also be removed from the gtfile, set to \code{TRUE} by default.
#' 
#' @return an object of class \code{gatingSet} with gate and children removed, as well as gatingTemplate file with population removed.
#' 
#' @export
removeGate <- function(gs, alias = NULL, gtfile = NULL, children = TRUE){
  
  # Supply alias
  if(is.null(alias)){
    
    stop("Please supply the name of the population to be removed.")
  
  }
  
  # Supply gtfile
  if(is.null(gtfile)){
    
    stop("Please supply the name of the gatingTemplate csv file to remove the gate.")
  
  }
  
  # Read in gatingTemplate
  gt <- read.csv(gtfile, header = TRUE)
  
  # For each alias remove from GatingSet and gatingTemplate
  lapply(alias, function(x){
    
    # Remove node & children from GatingSet
    Rm(x, gs)
    
    # Remove row entry from gatingTemplate
    if(children == TRUE){
    
    # Remove rows with alias or parent = alias
    if(length(which(gt$alias == x)) == 0){
      
      stop("Supplied alias is not a valid name for a gated population.")
      
    }else if(length(which(gt$parent == x)) != 0){
      
    indx <- c(which(gt$alias == x), which(gt$parent == x))
    gt <- gt[-indx,]  
    write.csv(gt, gtfile, row.names = FALSE)
      
    }else if(length(which(gt$parent == x)) == 0){
      
      indx <- which(gt$alias == x)
      gt <- gt[-indx,]  
      write.csv(gt, gtfile, row.names = FALSE)
      
    }
    
    }else if(children == FALSE){
      
      indx <- which(gt$alias == x)
      gt <- gt[-indx,]  
      write.csv(gt, gtfile, row.names = FALSE)
      
    }
  })
  
  return(gs)
  
}

#' Extract saved gate(s) from gatingTemplate.
#' 
#' @param parent name of the parental population.
#' @param alias name of the population for which the gate must be extracted.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g. "gatingTemplate.csv") where the gate(s) are saved.
#' 
#' @export
extractGate <- function(parent, alias, gtfile){
  
  # Load gtfile into gatingTemplate
  gt <- gatingTemplate(gtfile)
  
  # Extract population nodes from gt
  pops <- getNodes(gt, only.names = TRUE)
  
  # Extract gates given parent and child node(s)
  gates <- lapply(alias, function(x){
    
    if(parent == "root"){
    
    parent <- "root"
    alias <- paste("/",x,sep="")
  
  } else {
    
    parent <- paste("/",parent,"")
    alias <- paste("/",parent,"/",x,sep="")
  
  }  
  
  gate <- getGate(gt, parent, alias)
  gate <- eval(parameters(gate)$gate)
    
  })
  
  gates <- filters(gates)
  return(gates)
}

#' Edit existing gate(s).
#' 
#' @param x an object of class \code{GatingSet}.
#' @param pData vector of form \code{c("column","row")} indicating the rows of \code{pData(fs)} to extract. 
#' @param parent name of the parental population.
#' @param alias name(s) of the gate to edit (e.g. "Single Cells").
#' @param gate_type vector of gate type names used to construct the gates. Multiple \code{gate_types} are supported but should be accompanied with
#' an \code{alias} argument of the same length (i.e. one \code{gate_type} per \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#' ellipse, threshold, boundary, interval, quadrant and web} which can be abbreviated as upper or lower case first letters as well. Default \code{gate_type}
#' is \code{"polygon"}.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g. "gatingTemplate.csv") where the gate is saved.
#' 
#' @return an object of calss \code{GatingSet} with edited gate applied, as well as gatingTemplate file with editied gate saved.
#' 
#' @export
editGate <- function(x, pData = NULL, parent = NULL, alias = NULL, gate_type = NULL, gtfile = NULL){
  
  # Rename x to gs
  gs <- x
  
  # Restrict to samples matching pData requirements
  if(!is.null(pData)){
    
    # Extract samples using selectFrames
    fs <- selectFrames(fs, pData)
    
  }
  
  fr <- as(fs, "flowFrame")
  
  # Extract gate(s) from gtfile gating_args
  gates <- extractGate(gtfile = gtfile, parent = parent, alias = alias)  
  
  # If no gate_type supplied determine using getGateType
  if(is.null(gate_type)){
    
    gate_type <- getGateType(gates)
    
  }
  
  # Check gate_type argument is valid
  gate_type <- checkGateType(gate_type = gate_type, alias = alias)
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, gate_type = gate_type)  

  # Extract channels from gates
  channels <- parameters(gates[[1]])
  
  # Plot data for visualisation
  drawPlot(fr = fr, channels = channels)
  
  # Plot existing gates
  plotGates(gates, col = "magenta")
  
  # 2D Interval gates require axis argument
  if("interval" %in% gate_type){
    
  intvl <- rbind(gates[[match("interval", gate_type)[1]]]@min,gates[[match("interval", gate_type)[1]]]@max)
  
     if(all(is.finite(intvl[1,]))){
  
       axis <- "x"
    
     }else if(all(is.finite(intvl[2,]))){
       
       axis <- "y"
       
     }
  }
  
  # Make new call to drawGate to get new gates - set plot = FALSE
  new <- drawGate(fr, alias = alias, channels = channels, plot = FALSE, gate_type = gate_type, axis = axis)

  # Find and Edit gatingTemplate entries - each alias and gate separate
  gt <-read.csv(gtfile, header = TRUE)
  gt <- as.data.table(gt)
  
  for(i in 1:length(alias)){
    
    gt[gt$alias == alias[i] & gt$parent == parent, "gating_args"] <- .argDeparser(list(gate = new[[i]]))
    
  }
  
  gt <- data.frame(gt)
  write.csv(gt, gtfile, row.names = FALSE)
  
  # Apply entire gatingTemplate to GatingSet
  gt <- gatingTemplate(gtfile)
  rt <- getData(gs, "root")   # Remove all nodes - extract root
  gs <- GatingSet(rt)         # Add root to GatingSet
  
  gating(gt,gs)
  
  assign(deparse(substitute(x)), gs, envir=globalenv())
    
}

#' Get gate_type from from gate.
#' 
#' @param gates an object of class \code{filters} containing the gates from which the \code{gate_type(s)} must be obtained.
#'
#' @return vector of gate_type names for supplied gates.
#'
#' @export
getGateType <- function(gates){
  
  # Combine gate co-ordinates
  pts <- lapply(qt, function(x){rbind(x@min,x@max)})
  pts <- do.call(rbind,pts)
  
  # One gate supplied
  if(length(gates) == 1){
    
    if(class(gates[[1]]) == "ellipsoidGate"){
      
      # gate_type == "ellipse"
      types <- "ellipse"
    
    }else if(class(gates[[1]]) == "rectangleGate"){
      
      # Includes rectangle, interval, threshold and boundary gate_types
      if(length(parameters(gates[[1]])) == 1){
        
        # Gate in One Dimension
        if(is.infinite(gates[[1]]@min)){
          
          types <- "boundary"
          
        }else if(is.infinite(gates[[1]]@max)){
          
          types <- "threshold"
          
        }else if(is.finite(gates[[1]]@min) & is.finte(gates[[1]]@max)){
          
          types = "interval"
          
        }
        
      }else if(length(parameters(gates[[1]])) == 2){
        
        # Gate in 2 Dimensions
        if(is.infinite(gates[[1]]@min[1]) & is.infinite(gates[[1]]@min[2])){
          
          types = "boundary"
          
        }else if(is.infinite(gates[[1]]@max[1]) & is.infinite(gates[[1]]@max[2])){
          
          types <- "threshold"
          
        }else if(all(is.infinite(c(gates[[1]]@min[1], gates[[1]]@max[1]))) | all(is.infinite(c(gates[[1]]@min[2], gates[[1]]@max[2])))){
          
          types <- "interval"
          
        }else{
          
          types <- "rectangle"
          
        }
        
      }
      
    }else if(class(gates[[1]]) == "polygon"){
      
      # gate_type == "polygon"
      types <- "polygon"
      
    }
    
  # Multiple gates supplied
  }else if(length(gates) > 1){
    
    # Get classes of gates
    classes <- sapply(gates, function(x){
     
      class(x)   
     
    })
    
    # All gates are of the same class
    if(all(classes[1] == classes)){
      
      # Gates are all ellipses
      if(classes[1] == "ellipsoidGate"){
        
        types <- rep("ellipse", length(gates))
       
      # Gates are all rectangles 
      }else if(classes[1] == "rectangleGate"){
        
        # if 4 gates are supplied - gate_type may be "quadrant"
        if(length(gates) == 4){
          
          # Quadrant gates should have finite and infinite values in all gates and all finite co-ordinates should be the same
          if(sum(is.finite(pts[,1])) == 4 & sum(is.infinite(pts[,1])) == 4 & sum(duplicated(pts[,1][is.finite(pts[,1])])) == 3){
            
            types <- "quadrant"
           
          # Each gate could be either rectangle, interval, threshold or boundary
          }else{
            
            types <- sapply(gates, function(x){
              
              # Includes rectangle, interval, threshold and boundary gate_types
              if(length(parameters(x)) == 1){
                
                # Gate in One Dimension
                if(is.infinite(x@min)){
                  
                  types <- "boundary"
                  
                }else if(is.infinite(x@max)){
                  
                  types <- "threshold"
                  
                }else if(is.finite(x@min) & is.finte(x@max)){
                  
                  types = "interval"
                  
                }
                
              }else if(length(parameters(x)) == 2){
                
                # Gate in 2 Dimensions
                if(is.infinite(x@min[1]) & is.infinite(x@min[2])){
                  
                  types = "boundary"
                  
                }else if(is.infinite(x@max[1]) & is.infinite(x@max[2])){
                  
                  types <- "threshold"
                  
                }else if(all(is.infinite(c(x@min[1], x@max[1]))) | all(is.infinite(c(x@min[2], x@max[2])))){
                  
                  types <- "interval"
                  
                }else{
                  
                  types <- "rectangle"
                  
                }
                
              }
              
            })
            
          }
          
        }else{
          
            types <- rep("rectangle", length(gates))
            
        }
      
      # Gates are all polygons  
      }else if(classes[1] == "polygonGate"){
        
        # May be gate_type == "web" need to see if any points are conserved
        if(sum(duplicated(pts[,1][is.finite(pts[,1])])) == (length(gates)-1)){
          
          types <- "web"
          
        }else{
          
          types <- rep("polygon", length(gates))
          
        }
        
      }
      
    # Not all supplied gates are of the same class - treat separately
    }else{
      
      types <- sapply(gates, function(x){
        
        if(class(x) == "ellipsoidGate"){
          
          # gate_type == "ellipse"
          types <- "ellipse"
          
        }else if(class(x) == "rectangleGate"){
          
          # Includes rectangle, interval, threshold and boundary gate_types
          if(length(parameters(x)) == 1){
            
            # Gate in One Dimension
            if(is.infinite(x@min)){
              
              types <- "boundary"
              
            }else if(is.infinite(x@max)){
              
              types <- "threshold"
              
            }else if(is.finite(x@min) & is.finte(x@max)){
              
              types = "interval"
              
            }
            
          }else if(length(parameters(x)) == 2){
            
            # Gate in 2 Dimensions
            if(is.infinite(x@min[1]) & is.infinite(x@min[2])){
              
              types = "boundary"
              
            }else if(is.infinite(x@max[1]) & is.infinite(x@max[2])){
              
              types <- "threshold"
              
            }else if(all(is.infinite(c(x@min[1], x@max[1]))) | all(is.infinite(c(x@min[2], x@max[2])))){
              
              types <- "interval"
              
            }else{
              
              types <- "rectangle"
              
            }
            
          }
          
        }else if(class(x) == "polygon"){
          
          # gate_type == "polygon"
          types <- "polygon"
        }
        
      })
      
    }
  
  }
  
  return(types)
  
}