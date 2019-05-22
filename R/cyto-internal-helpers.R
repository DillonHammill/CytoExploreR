# CYTO_TRANSFORMED --------------------------------------------------------

#' .cyto_transformed
#'
#' @param x flowFrame, flowSet or GatingSet
#' @param trans transformList or transformerList object
#'
#' @return data which is appropriately transformed
#'
#' @importFrom flowWorkspace pData getData transformerList
#' @importFrom flowCore transform transformList
#'
#' @noRd
.cyto_transformed <- function(x, trans = NULL) {
  
  # Only flowFrame/flowSet/GatingSet
  if (!any(inherits(x, "flowFrame") |
           inherits(x, "flowSet") |
           class(x) == "GatingSet")) {
    stop("'x' must be either a flowFrame, flowSet or GatingSet.")
  }
  
  # Get comlete trans
  trans <- .cyto_transform_complete(x, trans)
  
  # Extract channels which have transformations
  if (inherits(trans, "transformList")) {
    chans <- names(trans@transforms)
  } else if (inherits(trans, "transformerList")) {
    chans <- names(trans)
  }
  
  # Extract summary stats
  if (inherits(x, "flowFrame")) {
    sm <- flowWorkspace::pData(flowCore::parameters(x))
  } else if (inherits(x, "flowSet")) {
    sm <- flowWorkspace::pData(flowCore::parameters(x[[1]]))
  } else if (inherits(x, "GatingSet")) {
    sm <- flowWorkspace::pData(flowCore::parameters(getData(x, "root")[[1]]))
  }
  
  # Extract channels that have been transformed
  chns <- as.vector(sm[, "name"][sm[, "maxRange"] < 6])
  
  # Check all chans have been transformed
  if (length(chns) == 0) {
    
    # No channels transformed
    x <- suppressMessages(flowCore::transform(x, trans))
  } else if (all(chans %in% chns)) {
    
    # All channels have been transformed
  } else {
    
    # Get transformations for untransformed channels
    if (inherits(trans, "transformList")) {
      trans <- transformList(
        chans[!chans %in% chns],
        trans@transforms[chans[!chans %in% chns]][[1]]@f
      )
    } else if (inherits(trans, "transformerList")) {
      trans <- transformerList(
        chans[!chans %in% chns],
        trans[chans[!chans %in% chns]]
      )
    }
    
    # Some channels have been transformed
    x <- suppressMessages(flowCore::transform(x, trans))
  }
  
  return(x)
}

# CYTO_RAW ----------------------------------------------------------------

#' .cyto_raw
#' return data which is untransformed - flowFrame/flowSet/GatingSet
#' GatingSet returns a flowSet of untransformed data at parent node
#' @noRd
.cyto_raw <- function(x, trans = NULL, parent = "root") {
  
  # Only flowFrame/flowSet/GatingSet
  if (!any(inherits(x, "flowFrame") |
           inherits(x, "flowSet") |
           class(x) == "GatingSet")) {
    stop("'x' must be either a flowFrame, flowSet or GatingSet.")
  }
  
  # Data is untransformed
  if (.cyto_transform_check(x) == FALSE) {
    if (inherits(x, "flowFrame") | inherits(x, "flowSet")) {
      return(x)
    } else if (inherits(x, "GatingSet")) {
      return(flowWorkspace::getData(x, parent))
    }
    
    # Data is transformed
  } else {
    if (inherits(x, "flowFrame") | inherits(x, "flowSet")) {
      if (is.null(trans)) {
        stop("Supply a transform object to inverse transformations.")
      }
    }
  }
  
  # Extract transformations from GatingSet
  if (is.null(trans) & inherits(x, "GatingSet")) {
    channels <- colnames(x)
    
    trnsfrms <- lapply(channels, function(channel) {
      getTransformations(x[[1]], channel, only.function = FALSE)
    })
    names(trnsfrms) <- channels
    
    # Remove NULL transforms
    trnsfrms[unlist(lapply(trnsfrms, is.null))] <- NULL
    trans <- transformerList(names(trnsfrms), trnsfrms)
  }
  
  # Get inverse trans
  inv <- cyto_transform_convert(trans, inverse = TRUE)
  
  # Extract channels which have transformations
  if (inherits(trans, "transformList")) {
    chans <- names(trans@transforms)
  } else if (inherits(trans, "transformerList")) {
    chans <- names(trans)
  }
  
  # Extract summary stats
  if (inherits(x, "flowFrame")) {
    sm <- flowWorkspace::pData(flowCore::parameters(x))
  } else if (inherits(x, "flowSet")) {
    sm <- flowWorkspace::pData(flowCore::parameters(x[[1]]))
  } else if (inherits(x, "GatingSet")) {
    sm <- pData(flowCore::parameters(flowWorkspace::getData(x, "root")[[1]]))
  }
  
  # Extract channels that have been transformed - apply inverse transform
  chns <- as.vector(sm[, "name"][sm[, "maxRange"] < 6])
  
  # Extract flowSet from GatingSet
  if (inherits(x, "GatingSet")) {
    x <- flowWorkspace::getData(x, parent)
  }
  
  # Check all chans have been transformed
  if (length(chns) == 0) {
    
    # No channels transformed
  } else if (all(chans %in% chns)) {
    
    # All channels have been transformed
    x <- flowCore::transform(x, inv)
  } else {
    
    # Some channels have been transformed
    trns <- lapply(chans[chans %in% chns], function(x) {
      inv@transforms[[x]]@f
    })
    names(trns) <- chans[chans %in% chns]
    inv <- transformList(names(trns), trns)
    x <- flowCore::transform(x, inv)
  }
  
  return(x)
}

# CYTO_TRANSFORM_CHECK ----------------------------------------------------------

#' .cyto_transform_check
#'
#' Check whether data has been transfomed - return TRUE if
#' any channels transformed
#'
#' @param x flowFrame, flowSet or GatingSet object to check
#'
#' @importFrom flowCore parameters
#' @importFrom flowWorkspace pData getData
#'
#' @noRd
.cyto_transform_check <- function(x) {
  if (inherits(x, "flowFrame")) {
    
    # Extract summary stats
    sm <- pData(parameters(x))
    
    # Check if any maxRange < 6
    if (any(sm[, "maxRange"] < 6)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (inherits(x, "flowSet")) {
    
    # Extract summary stats
    sm <- pData(parameters(x[[1]]))
    
    # Check if any maxRange < 6
    if (any(sm[, "maxRange"] < 6)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (inherits(x, "GatingSet")) {
    
    # Extract root flowSet
    fs <- getData(x, "root")
    
    # Extract summary stats
    sm <- pData(parameters(fs[[1]]))
    
    # Check if any maxRange < 6
    if (any(sm[, "maxRange"] < 6)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

# CYTO_TRANSFORM_COMPLETE ------------------------------------------------------

#' Get a complete transformation object
#'
#' @param x flowFrame, flowSet or GatingSet
#' @param trans transformList or transformerList
#'
#' @return complete transformList or transformerList object for all channels
#'
#' @importFrom flowCore estimateLogicle transformList
#' @importFrom flowWorkspace transformerList GatingSet getTransformations
#' @importFrom methods new
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_transform_complete <- function(x, trans = NA) {
  
  # Check class of trans
  if (!.all_na(trans)) {
    if (!any(inherits(trans, "transformList") |
             inherits(trans, "transformerList"))) {
      stop("'trans' should be a transformList or transformerList object.")
    }
  }
  
  # Extract transformations directly from GatingHierachy
  if(inherits(x, "GatingHierarchy")){
    if(.all_na(trans)){
      trans <- x@transformation
      if(is.null(trans)){
        trans <- NA
      }
    }
  }
  
  # Extract transformations directly from GatingSet
  if(inherits(x, "GatingSet")){
    if(.all_na(trans)){
      trans <- x[[1]]@transformation
      if(is.null(trans)){
        trans <- NA
      }
    }
  }
  
  # Extract fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # If NULL trans get all transformations
  if (.all_na(trans)) {
    if (inherits(x, "flowFrame")) {
      if (.cyto_transform_check(x) == TRUE) {
        stop(paste(
          "Looks like the data is already transformed.",
          "\n",
          "Please supply the transformList/transformerList used."
        ))
      }
      
      trans <- flowCore::estimateLogicle(x, channels)
      return(trans)
    } else if (inherits(x, "flowSet")) {
      if (.cyto_transform_check(x) == TRUE) {
        stop(paste(
          "Looks like the data is already transformed.",
          "\n",
          "Please supply the transformList/transformerList used."
        ))
      }
      
      trans <- flowCore::estimateLogicle(as(x, "flowFrame"), channels)
      return(trans)
    } else if (inherits(x, "GatingSet")) {
      if (.cyto_transform_check(x) == TRUE & length(x@transformation) == 0) {
        stop(paste(
          "Looks like the data is already transformed.",
          "\n",
          "Please supply the transformList/transformerList used."
        ))
      }
      
      # GatingSet is not transformed
      if (length(x@transformation) == 0) {
        
        # GatingSet is not transformed
        fs <- flowWorkspace::getData(x, "root")
        fr <- as(fs, "flowFrame")
        fs <- flowCore::flowSet(fr)
        gs <- suppressMessages(flowWorkspace::GatingSet(fs))
        
        trans <- flowCore::estimateLogicle(gs[[1]], channels)
        return(trans)
        
        # GatingSet contains transformations
      } else if (length(x@transformation) != 0) {
        chans <- names(x@transformation[[1]])
        
        if (any(chans %in% channels)) {
          
          # Extract transformations from GatingSet
          trnsfrms <- lapply(channels[chans %in% channels], function(channel) {
            getTransformations(x[[1]], channel, only.function = FALSE)
          })
          names(trnsfrms) <- channels[chans %in% channels]
          
          # Remove NULL transforms
          trnsfrms[unlist(lapply(trnsfrms, is.null))] <- NULL
          trans <- transformerList(names(trnsfrms), trnsfrms)
          
          if (all(channels %in% names(trans))) {
            
            # GatingSet contains all transformations
            return(trans)
          } else {
            
            # Get remaining transformations with estimateLogicle
            fs <- flowWorkspace::getData(x, "root")
            fr <- as(fs, "flowFrame")
            fs <- flowCore::flowSet(fr)
            gs <- suppressMessages(flowWorkspace::GatingSet(fs))
            
            trnsLst <- estimateLogicle(
              gs[[1]],
              channels[!channels %in% names(trans)]
            )
            trans <- c(trnsLst, trans)
            trans <- flowWorkspace::transformerList(names(trans), trans)
            
            return(trans)
          }
        } else {
          
          # GatingSet does not contain transformations for fluorescent channels
          fs <- flowWorkspace::getData(x, "root")
          fr <- as(fs, "flowFrame")
          fs <- flowCore::flowSet(fr)
          gs <- suppressMessages(flowWorkspace::GatingSet(fs))
          
          trnsLst <- flowCore::estimateLogicle(gs[[1]], channels)
          trans <- c(trnsLst, trans)
          trans <- flowWorkspace::transformerList(names(trans), trans)
          
          return(trans)
        }
      }
    }
  } else if (!.all_na(trans)) {
    
    # flowFrame or flowSet return transformList
    if (inherits(x, "flowFrame") | inherits(x, "flowSet")) {
      
      # Run cyto_transform_convert to get transformList
      trans <- cyto_transform_convert(trans, inverse = FALSE)
      
      # Check which channels have been transformed
      chans <- names(trans@transforms)
      
      # trans contains transformations for all fluorescent channels
      if (all(channels %in% chans)) {
        
        # trans is complete
        return(trans)
        
        # Some fluorescent channels don't have transformations
      } else {
        
        # Convert x to flowSet
        if (inherits(x, "flowFrame")) {
          fs <- flowCore::flowSet(x)
        } else if (inherits(x, "flowSet")) {
          fs <- x
        }
        
        # Generate merged flowFrame for use with estimateLogicle
        fr <- as(fs, "flowFrame")
        
        # Find channels excluded from trans
        excl <- channels[!channels %in% chans]
        
        # Get transformations for these channels using estimateLogicle
        trns <- flowCore::estimateLogicle(fr, excl)
        
        # Combine supplied trans with add transformations
        nms <- c(names(trans@transforms), excl)
        trans <- c(trans, trns)
        names(trans@transforms) <- nms
        
        return(trans)
      }
      
      # GatingSet return transformerList
    } else if (inherits(x, "GatingSet")) {
      
      # Supplied trans is a transformList - convert to transformerList
      if (inherits(trans, "transformList")) {
        chans <- names(trans@transforms)
        
        # Get transform functions
        trans <- lapply(seq_len(length(trans@transforms)), function(x) {
          trans@transforms[[x]]@f
        })
        names(trans) <- chans
        
        # Convert to transform objects
        trans <- lapply(seq_len(length(trans)), function(x) {
          t <- new("transform", .Data = trans[[1]])
          t@transformationId <- names(trans)[x]
          
          return(t)
        })
        
        trans <- lapply(trans, function(t) {
          inv <- flowCore::inverseLogicleTransform(trans = t)
          flowWorkspace::flow_trans("logicle", t@.Data, inv@.Data)
        })
        names(trans) <- chans
        trans <- flowWorkspace::transformerList(names(trans), trans)
      }
      
      # check which channels are covered by trans
      chans <- names(trans)
      
      # transformerList is complete
      if (all(channels %in% chans)) {
        return(trans)
      } else if (!all(channels %in% chans)) {
        
        # GatingSet contains some transformations
        if (length(x@transformation) != 0) {
          trnsfrms <- lapply(channels, function(channel) {
            getTransformations(x[[1]], channel, only.function = FALSE)
          })
          names(trnsfrms) <- channels
          
          # Remove NULL transforms
          trnsfrms[unlist(lapply(trnsfrms, is.null))] <- NULL
          trnsLst <- transformerList(names(trnsfrms), trnsfrms)
          
          # GatingSet contains some transformations
          if (any(channels %in% names(trnsLst))) {
            
            # GatingSet contains all transformations
            if (all(channels %in% names(trnsLst))) {
              return(trnsLst)
            } else {
              
              # GatingSet contains some transformations
              trnsLst <- trnsLst[names(trnsLst) %in% channels]
              
              # See if trans has any additional transformations
              if (any(names(trans) %in%
                      channels[!channels %in% names(trnsLst)])) {
                trans <- transformerList(
                  names(trans[names(trans) %in%
                                channels[!channels %in% names(trnsLst)]]),
                  trans[names(trans) %in%
                          channels[!channels %in% names(trnsLst)]]
                )
                trnsLst <- c(trnsLst, trans)
                trnsLst <- transformerList(names(trnsLst), trnsLst)
              }
              
              # See if all transformations are now present
              if (all(channels %in% names(trnsLst))) {
                return(trnsLst)
              } else {
                
                # Some channels are still missing transformations
                fs <- flowWorkspace::getData(x, "root")
                fr <- as(fs, "flowFrame")
                fs <- flowCore::flowSet(fr)
                gs <- suppressMessages(flowWorkspace::GatingSet(fs))
                
                trans <- estimateLogicle(
                  gs[[1]],
                  channels[!channels %in%
                             names(trnsLst)]
                )
                trans <- c(trnsLst, trans)
                trans <- flowWorkspace::transformerList(names(trans), trans)
                
                return(trans)
              }
            }
          }
          
          # GatingSet has no transformations
        } else if (length(x@transformation) == 0) {
          
          # trans contains all transformations
          if (all(channels %in% chans)) {
            return(trans)
            
            # Get remaining transformations from GatingSet using estimateLogicle
          } else {
            
            # Get remaining transformations with estimateLogicle
            fs <- flowWorkspace::getData(x, "root")
            fr <- as(fs, "flowFrame")
            fs <- flowCore::flowSet(fr)
            gs <- suppressMessages(flowWorkspace::GatingSet(fs))
            
            trnsLst <- estimateLogicle(gs[[1]], channels[!channels %in% chans])
            trans <- c(trnsLst, trans)
            trans <- flowWorkspace::transformerList(names(trans), trans)
            
            return(trans)
          }
        }
      }
    }
  }
}
