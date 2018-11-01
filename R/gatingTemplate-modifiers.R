#' Remove Gate(s) and Edit gatingTemplate .csv File
#'
#' @param gs an object of class \code{GatingSet}.
#' @param alias name(s) of the population(s) to remove (e.g. "Single Cells"). By
#'   default all descendant populations will be removed as well.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv").
#'
#' @return an object of class \code{gatingSet} with gate and children removed,
#'   as well as gatingTemplate file with population removed.
#'
#' @importFrom flowWorkspace getDescendants Rm
#' @importFrom utils read.csv write.csv
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
removeGate <- function(gs, alias = NULL, gtfile = NULL){
  
  # Supply alias
  if(is.null(alias)){
    
    stop("Please supply the name of the population to be removed.")
    
  }
  
  # Supply gtfile
  if(is.null(gtfile)){
    
    stop("Please supply the name of the gatingTemplate csv file to remove the gate.")
    
  }
  
  # Get children from GatingSet
  chldrn <- sapply(alias, function(x) basename(getDescendants(gs[[1]], x)))
  chldrn <- unlist(chldrn, use.names = FALSE)
  chldrn <- c(alias, unique(chldrn))
  
  # Read in gatingTemplate
  gt <- read.csv(gtfile, header = TRUE)
  
  # Remove all rows with alias = chldrn
  gt <- gt[!gt$alias %in% chldrn, ]
  
  # Remove nodes from GatingSet
  suppressMessages(Rm(alias, gs))
  
  write.csv(gt, gtfile, row.names = FALSE)
  
  assign(deparse(substitute(x)), gs, envir=globalenv())
  
}

#' Extract Saved Gate(s) from gatingTemplate.
#'
#' @param parent name of the parental population.
#' @param alias name of the population for which the gate must be extracted.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the gate(s) are saved.
#'
#' @importFrom flowWorkspace getGate getNodes
#' @importFrom openCyto gating
#' @importFrom flowCore filters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
extractGate <- function(parent, alias, gtfile){
  
  # Load gtfile into gatingTemplate
  gt <- suppressMessages(gatingTemplate(gtfile))
  
  # Extract population nodes from gt
  nds <- getNodes(gt, only.names = TRUE)
  
  #Parent Node
  parent <- names(nds[match(parent,nds)])
  
  # Extract gates given parent and child node(s)
  gates <- lapply(alias, function(x){
    
    # Alias node
    alias <- names(nds[match(x,nds)])
    
    gm <- getGate(gt, parent, alias)
    gate <- eval(parameters(gm)$gate)
    
  })
  
  gates <- filters(gates)
  return(gates)
}

#' Edit Existing Gate(s).
#'
#' @param x an object of class \code{GatingSet}.
#' @param select vector containing the indicies of samples within gs to use for
#'   plotting.
#' @param parent name of the parental population.
#' @param alias name(s) of the gate to edit (e.g. "Single Cells").
#' @param overlay name(s) of the population(s) to overlay onto the plot.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   \code{gate_types} are supported but should be accompanied with an
#'   \code{alias} argument of the same length (i.e. one \code{type} per
#'   \code{alias}). Supported \code{gate_types} are \code{polygon, rectangle,
#'   ellipse, threshold, boundary, interval, quadrant and web} which can be
#'   abbreviated as upper or lower case first letters as well. Default
#'   \code{type} is \code{"polygon"}.
#' @param gtfile name of the \code{gatingTemplate} csv file (e.g.
#'   "gatingTemplate.csv") where the gate is saved.
#' @param ... additional arguments passed to plotCyto, see ?plotCyto for
#'   details.
#'
#' @return an object of calss \code{GatingSet} with edited gate applied, as well
#'   as gatingTemplate file with editied gate saved.
#'
#' @importFrom flowWorkspace getData getTransformations GatingSet
#' @importFrom flowCore parameters
#' @importFrom openCyto gatingTemplate
#' @importFrom data.table as.data.table fread fwrite :=
#' @importFrom methods as
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
editGate <- function(x, select = NULL, parent = NULL, alias = NULL, overlay = NULL, type = NULL, gtfile = NULL, ...){
  
  # Rename x to gs
  gs <- x
  fs <- suppressMessages(getData(gs, parent))
  
  # Overlay
  if(!is.null(overlay)){
    
    # Extract populations to overlay - list of flowSets
    if(class(overlay) == "character"){
      
      overlay <- lapply(overlay, function(overlay){getData(gs, overlay)})
        
      overlay <- lapply(overlay, function(x){
          
        fr <- as(x,"flowFrame")
          
        if(is.na(match("Original", BiocGenerics::colnames(fr))) == FALSE){
            
          fr <- suppressWarnings(fr[, -match("Original", BiocGenerics::colnames(fr))])
            
        }
          
        return(fr)
          
      })
      
    }
    
  }
  
  # Extract transList from gs
  if(length(getTransformations(gs[[1]])) != 0){
    
    transList <- transformList(names(getTransformations(gs[[1]])), getTransformations(gs[[1]]))
    trnsfrmrLst <- transformerList(names(getTransformations(gs[[1]], only.function = FALSE)),getTransformations(gs[[1]], only.function = FALSE))
    
    trans <- TRUE
    
  }else{
    
    transList <- NULL
    
  }
  
  # Restrict to samples matching select requirements
  if(!is.null(select)){
    
    if(class(select) != "numeric"){
      
      stop("Vector supplied to select argument should contain the numeric indicies of the samples to select.")
      
    }
    
    # Extract samples using selectFrames
    fs <- fs[select]
    
  }
  fr <- as(fs, "flowFrame")
  
  if(is.na(match("Original", BiocGenerics::colnames(fr))) == FALSE){
    
    fr <- suppressWarnings(fr[, -match("Original", BiocGenerics::colnames(fr))])
    
  }
  
  # Extract gate(s) from gtfile gating_args
  gates <- extractGate(gtfile = gtfile, parent = parent, alias = alias)  
  
  # If no type supplied determine using getGateType
  if(is.null(type)){
    
    type <- getGateType(gates)
    
  }
  
  # Check type argument is valid
  type <- checkGateType(type = type, alias = alias)
  
  # Check alias is supplied correctly
  checkAlias(alias = alias, type = type)  
  
  # Extract channels from gates
  channels <- parameters(gates[[1]])
  
  # Plot data & existing gates
  plotCyto(fr, channels = channels, overlay = overlay, popup = TRUE, legend = FALSE, gates = gates, col.gate = "magenta", transList = transList, labels = FALSE, main = paste("Combined Events \n", parent), lwd.gate = 2.5, ...)
  
  # 2D Interval gates require axis argument
  if("interval" %in% type){
    
    intvl <- rbind(gates[[match("interval", type)[1]]]@min,gates[[match("interval", type)[1]]]@max)
    
    if(all(is.finite(intvl[,1]))){
      
      axis <- "x"
      
    }else if(all(is.finite(intvl[,2]))){
      
      axis <- "y"
      
    }
    
  }
  
  # Make new call to drawGate to get new gates - set plot = FALSE
  new <- drawGate(fr, alias = alias, channels = channels, type = type, axis = axis, plot = FALSE)
  
  # Re-name parent and pop to be data.table friendly
  prnt <- parent
  als <- alias
  
  # Find and Edit gatingTemplate entries - each alias and gate separate
  gt <- data.table::fread(gtfile)

  for(i in 1:length(alias)){
    
    gt[parent == prnt & alias == als[i], gating_args := .argDeparser(list(gate = new[[i]]))]
    
  }
  
  data.table::fwrite(gt, gtfile)
  
  # Apply entire gatingTemplate to GatingSet
  gt <- gatingTemplate(gtfile)
  rt <- suppressMessages(getData(gs, "root"))   # Remove all nodes - extract root
  
  # Inverse transformations
  if(length(getTransformations(gs[[1]])) != 0){
    
    inv <- transformList(names(trnsfrmrLst), lapply(trnsfrmrLst, `[[`, "inverse"))
    rt <- transform(rt,inv)
    
  }
  
  gs <- suppressMessages(GatingSet(rt))         # Add root to GatingSet
  
  if(trans == TRUE){
    
    gs <- transform(gs, trnsfrmrLst)
    
  }
  
  suppressMessages(gating(gt,gs))
  
  assign(deparse(substitute(x)), gs, envir=globalenv())
  
}

#' Get Gate Type from Saved Gate.
#'
#' @param gates an object of class \code{filters} containing the gates from
#'   which the \code{type(s)} must be obtained.
#'
#' @return vector of gate type names for supplied gates.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
getGateType <- function(gates){
  
  # One gate supplied
  if(length(gates) == 1){
    
    if(class(gates[[1]]) == "ellipsoidGate"){
      
      # type == "ellipse"
      types <- "ellipse"
      
    }else if(class(gates[[1]]) == "rectangleGate"){
      
      # Includes rectangle, interval, threshold and boundary gate_types
      if(length(parameters(gates[[1]])) == 1){
        
        # Gate in One Dimension
        if(is.infinite(gates[[1]]@min)){
          
          types <- "boundary"
          
        }else if(is.infinite(gates[[1]]@max)){
          
          types <- "threshold"
          
        }else if(is.finite(gates[[1]]@min) & is.finite(gates[[1]]@max)){
          
          types <- "interval"
          
        }
        
      }else if(length(parameters(gates[[1]])) == 2){
        
        # Gate in 2 Dimensions
        if(is.infinite(gates[[1]]@min[1]) & is.infinite(gates[[1]]@min[2])){
          
          types <- "boundary"
          
        }else if(is.infinite(gates[[1]]@max[1]) & is.infinite(gates[[1]]@max[2])){
          
          types <- "threshold"
          
        }else if(all(is.infinite(c(gates[[1]]@min[1], gates[[1]]@max[1]))) | all(is.infinite(c(gates[[1]]@min[2], gates[[1]]@max[2])))){
          
          types <- "interval"
          
        }else{
          
          types <- "rectangle"
          
        }
        
      }
      
    }else if(class(gates[[1]]) == "polygonGate"){
      
      # type == "polygon"
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
        
        # Combine gate co-ordinates
        pts <- lapply(gates, function(x){rbind(x@min,x@max)})
        pts <- do.call(rbind,pts)
        
        # if 4 gates are supplied - type may be "quadrant"
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
                  
                }else if(is.finite(x@min) & is.finite(x@max)){
                  
                  types <- "interval"
                  
                }
                
              }else if(length(parameters(x)) == 2){
                
                # Gate in 2 Dimensions
                if(is.infinite(x@min[1]) & is.infinite(x@min[2])){
                  
                  types <- "boundary"
                  
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
        
        # Combine gate co-ordinates
        pts <- lapply(gates, function(x){rbind(x@boundaries)})
        pts <- do.call(rbind,pts)
        dupl <- pts[pts[,1] == pts[which(duplicated(pts))[1],][1],]
        
        # May be type == "web" need to see if any points are conserved
        if(length(dupl[,1]) == 4){
          
          types <- "web"
          
        }else{
          
          types <- rep("polygon", length(gates))
          
        }
        
      }
      
      # Not all supplied gates are of the same class - treat separately
    }else{
      
      types <- sapply(gates, function(x){
        
        if(class(x) == "ellipsoidGate"){
          
          # type == "ellipse"
          types <- "ellipse"
          
        }else if(class(x) == "rectangleGate"){
          
          # Includes rectangle, interval, threshold and boundary gate_types
          if(length(parameters(x)) == 1){
            
            # Gate in One Dimension
            if(is.infinite(x@min)){
              
              types <- "boundary"
              
            }else if(is.infinite(x@max)){
              
              types <- "threshold"
              
            }else if(is.finite(x@min) & is.finite(x@max)){
              
              types <- "interval"
              
            }
            
          }else if(length(parameters(x)) == 2){
            
            # Gate in 2 Dimensions
            if(is.infinite(x@min[1]) & is.infinite(x@min[2])){
              
              types <- "boundary"
              
            }else if(is.infinite(x@max[1]) & is.infinite(x@max[2])){
              
              types <- "threshold"
              
            }else if(all(is.infinite(c(x@min[1], x@max[1]))) | all(is.infinite(c(x@min[2], x@max[2])))){
              
              types <- "interval"
              
            }else{
              
              types <- "rectangle"
              
            }
            
          }
          
        }else if(class(x) == "polygonGate"){
          
          # type == "polygon"
          types <- "polygon"
        }
        
      })
      
    }
    
  }
  
  return(types)
  
}
