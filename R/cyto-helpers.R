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
  
  return(TRUE)
  
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
  if(!cyto_check(x)){
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

#' Filter samples based on experiment variables
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param ... tidyverse-style subsetting using comma separated logical
#'   predicates based on experimental variables stored in \code{pData(x)}. See
#'   examples below for demonstration.
#'
#' @return \code{flowSet} or \code{GatingSet} restricted to samples which meet
#'   the filtering criteria.
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
#' fs <- cyto_filter(Activation,
#'                   Treatment == "Stim-C",
#'                   OVAConc %in% c(0,0.5))
#'
#' # Select Stim-A and Stim-C treatment groups
#' fs <- cyto_filter(Activation,
#'                   Treatment %in% c("Stim-A","Stim-C"))
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_filter <- function(x, ...){
  
  # Check class of x
  if(!any(inherits(x, "flowSet") |
     inherits(x, "GatingSet"))){
    stop("'x' should be an object of class flowSet or GatingSet.")
  }
  
  # Extract experiment details
  pd <- pData(x)
  
  # Perform filtering on pd to pull out samples
  pd_filter <- filter(pd, ...)
  
  # Get indices for selected samples
  if(nrow(pd_filter) == 0){
    message("No samples match the filtering criteria. Returning all samples.")
    ind <- seq_len(length(x))
  }else{
    ind <- match(pd_filter[,"name"], pd[,"name"])
  }
  
  return(x[ind])
  
}

# CYTO_SELECT ------------------------------------------------------------------

# Similar to cyto_filter but acts in a non-tidyverse way.

#' Select samples based on experiment variables
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param ... named list containing experimental variables to be used to select
#'   samples or named arguments containing the levels of the variables to
#'   select. See below examples for use cases.
#'
#' @return \code{flowSet} or \code{GatingSet} restricted to samples which meet
#'   the selection criteria.
#'
#' @examples
#' library(CytoRSuite)
#'
#' # Look at experiment details
#' pData(Activation)
#'
#' # Select Stim-C samples with 0 and 0.5 OVA concentrations
#' fs <- cyto_select(Activation,
#'                   Treatment = "Stim-C",
#'                   OVAConc = c(0,0.5))
#'
#' # Select Stim-A and Stim-C treatment groups
#' fs <- cyto_select(Activation,
#'                   list("Treatment" = c("Stim-A","Stim-C")))
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_select <- function(x, ...) {
  
  # Check class of x
  if(!any(inherits(x, "flowSet") |
          inherits(x, "GatingSet"))){
    stop("'x' should be an object of class flowSet or GatingSet.")
  }
  
  # Pull down ... arguments to list
  args <- list(...)
  
  # ... is already a named list of arguments
  if(class(args[[1]]) == "list"){
    args <- args[[1]]
  }
  
  # Extract experiment details
  pd <- pData(x)
  
  # Check that all varaibles are valid
  if(!all(names(args) %in% colnames(pd))){
    lapply(names(args), function(y){
      if(!y %in% names(pd)){
        stop(paste(y, "is not a valid variable in pData(x)."))
      }
    })
  }
  
  # Check that all variable levels at least exist in pd
  lapply(names(args), function(z){
    var <- factor(pd[,z], exclude = NULL) # keep <NA> as factor level
    lvls <- levels(var)
    # some variable levels do not exist in pd
    if(!all(args[[z]] %in% lvls)){
      lapply(args[[z]], function(v){
        if(!v %in% lvls){
          stop(paste0(v, " is not a valid level for ", z,"!"))
        }
      })
    }
  })
  
  # Get filtered pd
  pd_filter <- pd
  lapply(names(args), function(y){
    ind <- which(pd_filter[,y] %in% args[[y]])
    # No filtering if variable level is missing
    if(length(ind) != 0){
      pd_filter <<- pd_filter[ind, ]
    }
  })
  
  # Get indices for selected samples
  ind <- match(pd_filter[,"name"], pd[,"name"])
  
  return(x[ind])
  
}

# CYTO_GROUP_BY ----------------------------------------------------------------

#' Group a flowSet or GatingSet by experiment variables
#'
#' @param x an object of class \code{flowSet} or \code{GatingSet}.
#' @param parent name of the parent population to extract from GatingSet object.
#' @param group_by names of pData variables to use for merging. Set to "all" to
#'   merge all samples in \code{x}.
#'
#' @return a named list of \code{flowSet} or \code{GatingSet} objects
#'   respectively.
#'
#' @importFrom flowWorkspace pData sampleNames
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_group_by
#'
#' @export
cyto_group_by <- function(x,
                          group_by = "all"){
  
  # Check class of x
  if(!any(inherits(x,"flowSet")|inherits(x,"GatingSet"))){
    stop("'x' should be an object of class flowSet or GatingSet.")
  }
  
  # Extract experiment details
  pd <- pData(x)
  
  # Repace any NA with "NA" to avoid missing rows
  if(any(is.na(pd))){
    pd[is.na(pd)] <- "NA"
  }
  
  # Extract sample names
  nms <- sampleNames(x)
  
  # Check group_by
  if(!all(group_by %in% colnames(pd))){
    lapply(group_by, function(y){
      if(!y %in% colnames(pd)){
        stop(paste0(y, " is not a valid variable for this ", class(x),"."))
      }
    })
  }
  
  # Split pd based on group_by into a named list
  if(group_by == "all"){
    pd_split <- list("all" = pd)
  }else{
    pd_split <- split(pd, pd[,group_by], sep = " ")
  }
  
  # Replace each element of pd_split with matching samples
  x_list <- lapply(seq_len(length(pd_split)), function(z){
    
    x[pd_split[[z]][,"name"]]
    
  })
  
  return(x_list)
  
}

#' Merge samples by pData
#'
#' @param x flowSet or GatingSet object
#' @param parent name of the parent population to extract from GatingSet object.
#' @param group_by names of pData variables to use for merging. Set to "all" to
#'   merge all samples in the flowSet.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#'
#' @return list containing merged flowFrames, named with group.
#'
#' @importFrom flowWorkspace pData getData
#' @importFrom flowCore sampleFilter Subset
#'
#' @noRd
.cyto_merge <- function(x,
                        parent = "root",
                        group_by = "all",
                        display = NULL) {
  
  # check x
  if (inherits(x, "flowFrame") | inherits(x, "GatingHierarchy")) {
    stop("x must be either a flowSet or a GtaingSet object.")
  }
  
  # check group_by
  if (all(!group_by %in% c("all", colnames(pData(x))))) {
    stop("group_by should be the name of pData variables or 'all'.")
  }
  
  # Extract pData information
  pd <- pData(x)
  
  # Sort pd by group_by colnames
  if (!is.null(group_by)) {
    if (group_by[1] != "all") {
      pd <- pd[do.call("order", pd[group_by]), ]
    }
  }
  
  # flowSet for merging
  if (inherits(x, "GatingSet")) {
    fs <- getData(x, parent)
  } else {
    fs <- x
  }
  
  # group_by all samples
  if (length(group_by) == 1 & group_by[1] == "all") {
    pd$group_by <- rep("all", length(x))
    
    fr <- as(fs, "flowFrame")
    
    if ("Original" %in% BiocGenerics::colnames(fr)) {
      fr <- suppressWarnings(
        fr[, -match("Original", BiocGenerics::colnames(fr))]
      )
    }
    
    if (!is.null(display)) {
      fr <- Subset(fr, sampleFilter(size = display * BiocGenerics::nrow(fr)))
    }
    
    fr.lst <- list(fr)
    
    # group_by by one variable
  } else if (length(group_by) == 1) {
    pd$group_by <- pd[, group_by]
    
    fr.lst <- lapply(unique(pd$group_by), function(x) {
      fr <- as(fs[pd$name[pd$group_by == x]], "flowFrame")
      
      if ("Original" %in% BiocGenerics::colnames(fr)) {
        fr <- suppressWarnings(
          fr[, -match("Original", BiocGenerics::colnames(fr))]
        )
      }
      
      if (!is.null(display)) {
        fr <- Subset(fr, sampleFilter(size = display * BiocGenerics::nrow(fr)))
      }
      
      return(fr)
    })
    
    # group_by by multiple variables
  } else {
    pd$group_by <- do.call("paste", pd[, group_by])
    
    fr.lst <- lapply(unique(pd$group_by), function(x) {
      fr <- as(fs[pd$name[pd$group_by == x]], "flowFrame")
      
      if ("Original" %in% BiocGenerics::colnames(fr)) {
        fr <- suppressWarnings(
          fr[, -match("Original", BiocGenerics::colnames(fr))]
        )
      }
      
      if (!is.null(display)) {
        fr <- Subset(fr, sampleFilter(size = display * BiocGenerics::nrow(fr)))
      }
      
      return(fr)
    })
  }
  names(fr.lst) <- unique(pd$group_by)
  
  return(fr.lst)
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

# CYTO_COMPENSATE --------------------------------------------------------------

#' Apply fluorescence compensation to samples
#'
#' \code{cyto_compensate} will apply a saved spillover matrix csv file to all
#' samples. The csv file must contain the fluorescent channels as both the
#' colnames (i.e. first row) and rownames (i.e. first column). If no
#' \code{spillover} csv file is supplied, the spillover matrix will be extracted
#' directly from the first element of \code{x}. To select a different sample for
#' spillover matrix extraction supply the index or name of the sample to the
#' \code{select} argument.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet} or
#'   \code{GatingSet}.
#' @param spillover name of the spillover matrix csv file (e.g.
#'   "Spillover-Matrix.csv") saved in the current working directory.
#' @param select index or name of the sample from which the spillover matrix
#'   should be extracted when no spillover matrix file is supplied to
#'   \code{spillover}. To compensate each sample individually using their stored
#'   spillover matrix file, set \code{select} to NULL.
#'
#' @return a compensated \code{flowFrame}, \code{flowSet} or \code{GatingSet}
#'   object.
#'
#' @importFrom utils read.csv
#' @importFrom flowWorkspace sampleNames getData
#'
#' @examples
#' library(CytoRSuiteData)
#'
#' # Apply stored spillover matrix to flowSet
#' cyto_compensate(Activation)
#'
#' # Save spillover matrix in correct format
#' spill <- Activation[[1]]@description$SPILL
#' rownames(spill) <- colnames(spill)
#' write.csv(spill,
#'           "Spillover-Matrix.csv")
#'           
#' # Apply saved spillover matrix csv file to flowSet
#' cyto_compensate(Activation, "Spillover-Matrix.csv")
#'
#' # Apply stored spillover matrix to GatingSet
#' gs <- GatingSet(Activation)
#' cyto_compensate(gs)
#'
#' # Apply saved spillover matrix csv file to GatingSet
#' cyto_compensate(gs, "Spillover-Matrix.csv")
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_compensate
#'
#' @export
cyto_compensate <- function(x, ...){
  UseMethod("cyto_compensate")
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowFrame <- function(x,
                                      spillover = NULL,
                                      select = 1){
  
  # Read in spillover matrix file
  if(!is.null(spillover)){
    if(getOption("CytoRSuite_wd_check")){
      if(!file_wd_check(spillover)){
        stop(paste(spillover, "does not exist in this working directory."))
      }
    }
    spill <- read.csv(spillover, header = TRUE, row.names = 1)
    colnames(spill) <- rownames(spill)
  }
  
  # Extract spillover directly from x
  if(is.null(spillover)){
    spill <- x@description$SPILL
  }
  
  # Apply compensation
  flowCore::compensate(x, spill)
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowSet <- function(x, 
                                    spillover = NULL,
                                    select = 1){
  
  # Read in spillover matrix file
  if(!is.null(spillover)){
    if(getOption("CytoRSuite_wd_check")){
      if(!file_wd_check(spillover)){
        stop(paste(spillover, "does not exist in this working directory."))
      }
    }
    spill <- read.csv(spillover, header = TRUE, row.names = 1)
    colnames(spill) <- rownames(spill)
    
    # Convert spill into a named list
    spill <- rep(list(spill), length(x))
    names(spill) <- sampleNames(x)
  }
  
  # Extract spillover directly from x
  if(is.null(spillover)){
    if(!is.null(select)){
      spill <- rep(list(x[[select]]@description$SPILL), length(x))
      names(spill) <- sampleNames(x)
    }else{
      spill <- lapply(sampleNames(x), function(y){
        x[[y]]@description$SPILL
      })
    }

  }

  # Apply compensation
  flowCore::compensate(x, spill)
  
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.GatingSet <- function(x, 
                                      spillover = NULL,
                                      select = 1){
  
  # Extract flowSet
  fs <- getData(gs, "root")
  
  # Read in spillover matrix file
  if(!is.null(spillover)){
    if(getOption("CytoRSuite_wd_check")){
      if(!file_wd_check(spillover)){
        stop(paste(spillover, "does not exist in this working directory."))
      }
    }
    spill <- read.csv(spillover, header = TRUE, row.names = 1)
    colnames(spill) <- rownames(spill)
    
    # Convert spill into a named list
    spill <- rep(list(spill), length(x))
    names(spill) <- sampleNames(x)
  }
  
  # Extract spillover directly from x
  if(is.null(spillover)){
    if(!is.null(select)){
      spill <- rep(list(fs[[select]]@description$SPILL), length(x))
      names(spill) <- sampleNames(x)
    }else{
      spill <- lapply(sampleNames(x), function(y){
        fs[[y]]@description$SPILL
      })
    }
    
  }
  
  # Apply compensation
  flowWorkspace::compensate(x, spill)
  
}

# CYTO_TRANSFORM ---------------------------------------------------------------

# Possibly assign transform object to global option that can be accessed
# internally by cyto_plot. 

# Challenges:
# - which samples do you use with estimateLogicle? Maybe sample a percentage of 
#   samples and pool? 
# - needs to support adjustment of transformation parameters on a per channel 
#   basis.

# CYTO_INVERSE_TRANSFORM -------------------------------------------------------

# Extract transform object from global option or supplied. Get inverse
# transformations for designated channels and apply to samples.

# These transform functions can then be used by cyto_plot to make in-line
# transformations. Need to come up with a new argument/s to implement this.
