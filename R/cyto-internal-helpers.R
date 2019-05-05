#' Check a flowFrame, flowSet, GatingHierarchy or GatingSet has been supplied
#'
#' @param x object of class flowFrame, flowSet, GatingHierarchy or GatingSet to
#'   be checked.
#'   
#' @return TRUE or FALSE if object meets class criteria.   
#'   
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}   
#'   
#' @noRd
.cyto_check <- function(x){
  
  # Check for valid class of object
  if(!any(inherits(x, "flowFrame") |
     inherits(x, "flowSet") |
     inherits(x, "GatingHierarchy") |
     inherits(x, "GatingSet"))){
    
    stop("'x' should be a flowFrame, flowSet, GatingHierarchy or GatingSet.")
    
  }
  
}

#' Extract a valid flowFrame or flowSet
#'
#' @param x object of class flowFrame, flowSet, GatingHierarchy or GatingSet.
#' @param parent name of the parent population to extract from GatingHierachy or
#'   GatingSet objects.
#'   
#' @return either a flowFrame or flowSet object
#' 
#' @importFrom flowWorkspace getData
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}   
#' 
#' @noRd
.cyto_extract <- function(x, parent = "root"){
  
  # Throw error for invalid object
  if(!.cyto_check(x)){
    stop(
      paste("'x' should be either a flowFrame, flowSet, GatingHierarchy",
            "or GatingSet.")
    )
  }
  
  # Extract flowFrame from GatingHierarchy
  if(inherits(x, "GatingHierarchy")){
    x <- getData(x, parent = parent)
  }
  
  # Extract flowSet from GatingSet
  if(inherits(x, "GatingSet")){
    x <- getData(x, parent = parent)
  }
  
  return(x)
}

#' Convert between cytometry objects
#'
#' Automatically removes 'Original' column for coerced objects.
#'
#' @param x flowFrame, flowSet, GatingHierarchy, GatingSet or flowFrame list or
#'   flowSet list.
#' @param return either 'flowFrame', 'flowSet', 'flowFrame list' or 'flowSet
#'   list'.
#' @param parent name of parent population to extract from GtaingHierarchy and
#'   GatingSet objects.
#' @param display percentage of events to include.
#'
#' @return object specified by 'return' argument.
#'
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace getData
#' @importFrom methods as
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_convert <- function(x,
                          return = "flowFrame",
                          parent = "root",
                          display = 1) {
  
  # Convert GatingHierarchy/GatingSet to flowFrame/flowSet objects
  if(inherits(x, "GatingHierarchy") | inherits(x, "GatingSet")){
    x <- getData(x, parent)
  }
  
  if(inherits(x, "flowFrame")){
    
    # Sampling
    if(display != 1){
      x <- cyto_sample(x, display = display)
    }
    
    if(return == "flowFrame list"){
      x <- list(x)
    }else if(return == "flowSet"){
      x <- flowSet(x)
    }else if(return == "flowSet list"){
      x <- list(flowSet(x))
    }
    
  }else if(inherits(x, "flowSet")){
    
    # Sampling
    if(display != 1){
      x <- fsApply(x, function(fr){cyto_sample(fr, display = display)})
    }
    
    if(return == "flowFrame"){
      x <- as(x, "flowFrame")
      if ("Original" %in% BiocGenerics::colnames(x)) {
        x <- suppressWarnings(
          x[, -match("Original", BiocGenerics::colnames(x))]
        )
      }
    }else if(return == "flowFrame list"){
      x <- lapply(seq_len(length(x)), function(y){x[[y]]})
    }else if(return == "flowSet list"){
      x <- list(x)
    }
    
  }else if(all(unlist(lapply(x,"class")) == "flowFrame")){
    
    # Sampling
    if(display != 1){
      x <- lapply(x, function(fr){cyto_sample(fr, display = display)})
    }
    
    if(return == "flowFrame"){
      x <- flowSet(x)
      x <- as(x, "flowFrame")
      if ("Original" %in% BiocGenerics::colnames(x)) {
        x <- suppressWarnings(
          x[, -match("Original", BiocGenerics::colnames(x))]
        )
      }
    }else if(return == "flowSet"){
      x <- flowSet(x)
    }else if(return == "flowSet list"){
      x <- list(flowSet(x))
    }
    
  }else if(all(unlist(lapply(x, function(x){inherits(x,"flowSet")})))){
    
    # Sampling
    if(display != 1){
      x <- lapply(x, function(fs){
        fsApply(fs, function(fr){cyto_sample(fr, display = display)})
      })
    }

    if(return == "flowFrame"){
      x <- lapply(x, function(fs){
        fr <- as(fs, "flowFrame")
        if ("Original" %in% BiocGenerics::colnames(fr)) {
          fr <- suppressWarnings(
            fr[, -match("Original", BiocGenerics::colnames(fr))]
          )
        }
        return(fr)
      })
      x <- flowSet(x)
      x <- as(x,"flowFrame")
      if ("Original" %in% BiocGenerics::colnames(x)) {
        x <- suppressWarnings(
          x[, -match("Original", BiocGenerics::colnames(x))]
        )
      }
      
    }else if(return == "flowFrame list"){
      
      x <- lapply(x, function(fs){
        fr <- as(fs, "flowFrame")
        if ("Original" %in% BiocGenerics::colnames(fr)) {
          fr <- suppressWarnings(
            fr[, -match("Original", BiocGenerics::colnames(fr))]
          )
        }
        return(fr)
      })
      
    }else if(return == "flowSet"){
      
      x <- lapply(x, function(fs){
        fr <- as(fs, "flowFrame")
        if ("Original" %in% BiocGenerics::colnames(fr)) {
          fr <- suppressWarnings(
            fr[, -match("Original", BiocGenerics::colnames(fr))]
          )
        }
        return(fr)
      })
      
      x <- flowSet(x)
      
    }
  }
  return(x)
}
