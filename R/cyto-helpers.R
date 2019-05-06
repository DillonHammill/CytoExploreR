# CYTO_CHECK -------------------------------------------------------------------

#' Check a flowFrame, flowSet, GatingHierarchy or GatingSet has been supplied
#'
#' @param x object of class flowFrame, flowSet, GatingHierarchy or GatingSet to
#'   be checked.
#'   
#' @return TRUE or FALSE if object meets this class criteria.   
#'   
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Valid object
#' cyto_check(Activation)
#' 
#' # Invalid list object
#' cyto_check(list(Activation))
#'   
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}   
#'   
#' @export
cyto_check <- function(x){
  
  # Check for valid class of object
  if(!any(inherits(x, "flowFrame") |
          inherits(x, "flowSet") |
          inherits(x, "GatingHierarchy") |
          inherits(x, "GatingSet"))){
    
    stop("'x' should be a flowFrame, flowSet, GatingHierarchy or GatingSet.")
    
  }
  
}

# CYTO_EXTRACT -----------------------------------------------------------------

#' Extract a valid flowFrame or flowSet
#'
#' \code{cyto_extract} is essentially a wrapper for
#' \code{flowWorkspace::getData()} which also accepts flowFrame or flowSet
#' objects. The \code{parent} population is extracted from GatingHierarchy or
#' GatingSet objects whilst flowFrame or flowSet objects are returned as is.
#'
#' @param x object of class flowFrame, flowSet, GatingHierarchy or GatingSet.
#' @param parent name of the parent population to extract from GatingHierachy or
#'   GatingSet objects.
#'
#' @return either a flowFrame or flowSet object
#'
#' @importFrom flowWorkspace getData
#'
#' @examples
#' library(CytoRSuiteData)
#'
#' # GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Extract root population
#' cyto_extract(gs, "root")
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_extract <- function(x, parent = "root"){
  
  # Throw error for invalid object
  if(!.cyto_check(x)){
    stop(
      paste("'x' should be either a flowFrame, flowSet, GatingHierarchy",
            "or GatingSet.")
    )
  }
  
  # Extract flowFrame from GatingHierarchy
  if(inherits(x, "GatingHierarchy")){
    x <- getData(x, parent)
  }
  
  # Extract flowSet from GatingSet
  if(inherits(x, "GatingSet")){
    x <- getData(x, parent)
  }
  
  return(x)
}

# CYTO_CONVERT -----------------------------------------------------------------

#' Convert between cytometry objects
#'
#' Automatically removes 'Original' column for coerced objects.
#'
#' @param x flowFrame, flowSet, GatingHierarchy, GatingSet or flowFrame list or
#'   flowSet list.
#' @param return either 'flowFrame', 'flowSet', 'flowFrame list' or 'flowSet
#'   list'.
#' @param parent name of parent population to extract from GatingHierarchy and
#'   GatingSet objects.
#' @param display percentage of events to include.
#'
#' @return object specified by 'return' argument.
#'
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace getData GatingSet
#' @importFrom methods as
#'
#' @examples 
#' library(CytoRSuiteData)
#' 
#' # Convert flowSet to list of flowFrames
#' cyto_convert(Activation, "flowFrame list")
#' 
#' # Convert flowSet to flowFrame
#' cyto_convert(Activation, "flowFrame")
#' 
#' # Convert GatingSet to flowFrame
#' cyto_convert(GatingSet(Activation), "flowFrame", parent = "T Cells")
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_convert <- function(x,
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

# CYTO_FILTER ------------------------------------------------------------------

#' Select samples based on experiment variables
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param ...
#'
#' @return \code{flowSet} or \code{GatingSet} restricted to samples which meet
#'   the selection criteria.
#'   
#' @importFrom flowWorkspace pData
#' @importFrom dplyr filter
#' @importFrom tibble as_tibble
#'   
#' @examples
#' library(CytoRSuite)
#' 
#' # Look at experiment details
#' pData(Activation)
#' 
#' # Select Stim-C samples with 0 and 0.5 OVA concentrations
#' fs <- cyto_filter(Activation, Treatment == "Stim-C", OVAConc %in% c(0,0.5))
#'   
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
cyto_filter <- function(x, ...){
  
  # Extract experiment details
  pd <- pData(x)
  
  # Convert pd to tibble for filtering'
  pd <- as_tibble(pd)
  
  # Perform filtering on pd to pull out samples
  pd <- filter(pd, ...)
  
  return(x[pd[,"name"]])
  
}


# CYTO_SAMPLE ------------------------------------------------------------------

#' Sample a flowFrame
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param display numeric [0,1] indicating the percentage of events to keep.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
#'
#' @return \code{\link[flowCore:flowFrame-class]{flowFrame}} restricted to
#'   \code{display} percentage events.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter
#' @importFrom flowCore Subset
#'
#' @examples
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Constrict first sample by 50%
#' cyto_sample(fs[[1]], 0.5)
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
cyto_sample <- function(x, 
                        display = 1, 
                        seed = NULL) {
  
  # No sampling
  if(display == 1){
    
  }else{
    # Number of events
    events <- nrow(x)
    
    # Size
    size <- display * events
    
    # Set seed
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    # Sample
    smp <- sampleFilter(size = size)
    x <- Subset(x, smp)
  }
  
  return(x)
}

# CYTO_MARKERS -----------------------------------------------------------------

#' Assign marker names to flowFrame or flowSet
#'
#' \code{cyto_markers} opens an editable table containing a list of channels and
#' markers for a \code{flowFrame} or \code{flowSet}. Users can edit the
#' \code{name} or \code{desc} columns with updated channel names or marker names
#' respectively. These entries will be updated in the \code{flowFrame} or
#' \code{flowSet} upon closing the window and saved to a
#' "Experiment-markers.csv" file for future use.
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#' @param file name of csv file containing columns 'Channel' and 'Marker'.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters markernames markernames<-
#' @importFrom utils edit write.csv read.csv
#'
#' @return save inputs to "Experiment-Markers.csv" and update marker names of
#'   \code{x}.
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Add marker names to channels - edit table
#' cyto_markers(fs)
#' }
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_markers <- function(x, file = NULL) {
  if (!any(inherits(x, "flowFrame") | inherits(x, "flowSet"))) {
    stop("Please supply either a flowFrame or flowSet object")
  }
  
  if (inherits(x, "flowFrame")) {
    
    # Extract pData of parameters
    pd <- pData(parameters(x))
  } else if (inherits(x, "flowSet")) {
    
    # Extract pData of parameters
    pd <- pData(parameters(x[[1]]))
  }
  
  # file missing
  if (is.null(file)) {
    
    if(length(grep("Experiment-Markers.csv", list.files())) != 0){
      
      message("Experiment-Markers.csv found in working directory.")
      dt <- read.csv(list.files()[grep("Experiment-Markers.csv", list.files())],
                     header = TRUE, stringsAsFactors = FALSE)
      
    }else{
      
      # Make data.frame with channel and marker columns
      dt <- pd[, c("name", "desc")]
      colnames(dt) <- c("Channel", "Marker")
      rownames(dt) <- NULL
      
    }
    
  } else {
    if (getOption("CytoRSuite_wd_check")) {
      if (!.file_wd_check(file)) {
        stop(paste(file, "does not exist in this working directory."))
      }
    }
    dt <- read.csv(file, header = TRUE)
  }
  
  # Channels with markers
  chans <- as.vector(dt$Channel[!is.na(dt$Marker)])
  
  # Edit dt
  dt <- suppressWarnings(edit(dt))
  
  # Update channels
  BiocGenerics::colnames(x) <- dt$Channel
  
  # Write result to csv file
  if (length(grep("Experiment-Markers.csv", c(file, list.files()))) != 0) {
    
    write.csv(dt, c(file, list.files())[grep("Experiment-Markers.csv",
                                             c(file, list.files()))[1]]
              , row.names = FALSE)
    
  } else {
    write.csv(dt, paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Markers.csv"
    ), row.names = FALSE)
  }
  
  # Channels with markers added
  tb <- dt[!dt$Channel %in% chans, ]
  chns <- tb$Channel[!is.na(tb$Marker)]
  
  # Pull out assigned markers
  mrk <- dt$Marker[dt$Channel %in% c(chans, chns)]
  names(mrk) <- dt$Channel[dt$Channel %in% c(chans, chns)]
  
  # Assign markers to x
  markernames(x) <- mrk
  
  invisible(NULL)
}

# CYTO_ANNOTATE ----------------------------------------------------------------

#' Interactively edit pData Information for a flowSet or GatingSet
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param file name of csv file containing experimental information.
#'
#' @return NULL and update pData for the \code{flowSet} or \code{GatingSet}.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore pData<-
#' @importFrom utils edit write.csv read.csv
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Edit pData in table editor
#' cyto_annotate(fs)
#' }
#' 
#' @export
cyto_annotate <- function(x, file = NULL) {
  
  # x should be a flowSet or GatingSet
  if (!any(inherits(x, "flowSet") | inherits(x, "GatingSet"))) {
    stop("Please supply either a flowSet or a GatingSet")
  }
  
  # Assign x to cyto
  cyto <- x
  
  # File missing
  if (is.null(file)) {
    
    if(length(grep("Experiment-Details.csv", list.files())) != 0){
      
      message("Experiment-Details.csv found in working directory.")
      pd <- read.csv(list.files()[grep("Experiment-Details.csv", list.files())],
                     header = TRUE, stringsAsFactors = FALSE)
      
    }else{
      
      # Extract pData
      pd <- pData(cyto)
      rownames(pd) <- NULL
      
    }
    
  } else {
    if (getOption("CytoRSuite_wd_check")) {
      if (!.file_wd_check(file)) {
        stop(paste(file, "does not exist in this working directory"))
      }
    }
    
    # Read in file
    pd <- read.csv(file, header = TRUE)
    
    # Update pData
    pData(cyto) <- pd
  }
  
  # Edit pData
  pd <- suppressWarnings(edit(pd))
  rownames(pd) <- pd$name
  
  # Update pData
  pData(cyto) <- pd
  
  # Write result to csv file
  if (length(grep("Experiment-Details.csv", c(file, list.files()))) != 0) {
    
    write.csv(pd, c(file, list.files())[grep("Experiment-Details.csv",
                                             c(file, list.files()))[1]]
              , row.names = FALSE)
    
  } else {
    
    write.csv(pData(cyto), paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Details.csv"
    ), row.names = FALSE)
    
  }
  
  # Update globally
  assign(deparse(substitute(x)), cyto, envir = globalenv())
  
  invisible(NULL)
}
