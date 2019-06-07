# CYTO_LOAD --------------------------------------------------------------------

#' Load .fcs files into ncdfFlowSet
#'
#' \code{cyto_load} is a convenient wrapper around
#' \code{\link[base:list.files]{list.files}} and
#' \code{\link[ncdfFlow:read.ncdfFlowSet]{read.ncdfFlowSet}} which makes it easy
#' to load .fcs files into a ncdfFlowSet.
#'
#' @param path points to the location of the .fcs files to read in.
#' @param ... additional arguments passed to read.ncdfFlowSet.
#'
#' @return \code{ncdfFlowSet}
#'
#' @importFrom ncdfFlow read.ncdfFlowSet
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Get path to Activation .fcs files in CytoRSuiteData
#' datadir <- system.file("extdata", package = "CytoRSuiteData")
#' path <- paste0(datadir, "/Activation")
#' 
#' # Load in .fcs files into ncdfFlowSet
#' fs <- cyto_load(path)
#' 
#' # fs is a ncdfFlowSet
#' class(fs)
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_load <- function(path = ".", ...) {

  # Get paths to files
  files <- list.files(path, full.names = TRUE, pattern = ".fcs")

  # Read into a ncdfFlowSet
  read.ncdfFlowSet(files = files, ...)
}

# CYTO_SETUP -------------------------------------------------------------------

#' Load.fcs files into GatingSet and annotate with experiment details
#'
#' \code{cyto_setup} takes care of all the data loading and annotation steps to
#' prepare your cytometry data for downstream analyses. The .fcs files are first
#' read into a \code{ncdfFlowSet} using \code{cyto_load} which is then added to
#' a \code{GatingSet}. Calls are then made to \code{cyto_markers} and
#' \code{cyto_annotate} to update the GtaingSet with the details of the
#' experiment. These details can be modified later with additional calls to
#' \code{cyto_markers} and/or \code{cyto_annotate}. Users are also asked to
#' provide a name for a gatingTemplate csv file which will be created if
#' necessary and assigned as the active gatingTemplate.
#'
#' @param path points to the location of the .fcs files to read in.
#' @param gatingTemplate name of a gatingTemplate csv file to be used for gate
#'   saving.
#' @param ... additional arguments passed to read.ncdfFlowSet.
#'
#' @return \code{GatingSet}
#'
#' @importFrom flowWorkspace GatingSet
#' @importFrom tools file_ext
#'
#' @examples
#' library(CytoRSuiteData)
#'
#' # Get path to Activation .fcs files in CytoRSuiteData
#' datadir <- system.file("extdata", package = "CytoRSuiteData")
#' path <- paste0(datadir, "/Activation")
#'
#' # Load in .fcs files into an annotated GatingSet
#' fs <- cyto_setup(path)
#'
#' # Markers have been assigned
#' cyto_extract(gs, "root")[[1]]
#'
#' # Experiment details have been updated
#' cyto_details(gs)
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_setup <- function(path = ".",
                       gatingTemplate = NULL, ...) {

  # Load in .fsc files to ncdfFlowSet
  fs <- cyto_load(path, ...)

  # Add fs to GatingSet
  message("Adding samples to a GatingSet.")
  gs <- GatingSet(fs)

  # Markers
  message("Associate markers with their respective channels.")
  cyto_markers(gs)

  # Annotate
  message("Annotate samples with experiment details.")
  cyto_annotate(gs)
  
  # Assign gatingTemplate
  if(is.null(gatingTemplate)){
    
    if(interctive()){
      # User prompt
      gatingTemplate <- readline("Provide a name for the gatingTemplate:")
    }
    
  }
  
  # No gatingTemplate specified in user prompt
  if(!is.null(gatingTemplate)){
    if(.empty(gatingTemplate)){
       gatingTemplate <- NULL
    }
  }
  
  # Check gatingTemplate
  if(!is.null(gatingTemplate)){
    
    # Add file extension if missing
    if(.empty(file_ext(gatingTemplate))){
      gatingTemplate <- paste0(gatingTemplate,".csv")
    }
    
    # Assign globally
    message(paste("Setting", gatingTemplate, "as the active gatingTemplate."))
    cyto_gatingTemplate_select(gatingTemplate)
    
    # Write new csv file if not in current directory
    if(!any(grepl(gatingTemplate, list.files()))){
      message(paste("Creating", gatingTemplate, "."))
      cyto_gatingTemplate_create(gatingTemplate)
    }
    
  }

  return(gs)
}

# CYTO_DETAILS -----------------------------------------------------------------

#' Extract experiment details from flowSet or GatingSet
#'
#' Simply an autocomplete-friendly wrapper around
#' \code{\link[flowWorkspace:pData]{pData}}.
#'
#' @param object of class \code{flowSet} or \code{GatingSet}.
#'
#' @importFrom flowWorkspace pData
#'
#' @return experiment details as data.frame.
#'
#' @export
cyto_details <- function(x) {
  pData(x)
}

# CYTO_NAMES -------------------------------------------------------------------

#' Extract sample names
#'
#' Simply a convenient and autocomplete-friendly wrapper around
#' \code{\link[flowCore:identifier]{identifier}}
#' \code{\link[flowWorkspace:sampleNames]{sampleNames}} to extract the sample
#' names from flowFrame, flowSet GatingHierarchy or GatingSet.
#'
#' @param object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#'
#' @importFrom flowCore identifier
#' @importFrom flowWorkspace sampleNames
#'
#' @return names associated with the supplied object.
#'
#' @rdname cyto_names
#'
#' @export
cyto_names <- function(x) {
  UseMethod("cyto_names")
}

#' @rdname cyto_names
#' @export
cyto_names.flowFrame <- function(x){
  identifier(x)
}

#' @rdname cyto_names
#' @export
cyto_names.flowSet <- function(x){
  sampleNames(x)
}

#' @rdname cyto_names
#' @export
cyto_names.GatingHierarchy <- function(x){
  x@name
}

#' @rdname cyto_names
#' @export
cyto_names.GatingSet <- function(x){
  sampleNames(x)
}

#' @rdname cyto_names
#' @export
cyto_names.list <- function(x){
  unlist(lapply(x,"cyto_names"))
}

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
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_check <- function(x) {

  # Check for valid class of object
  if (!any(inherits(x, "flowFrame") |
    inherits(x, "flowSet") |
    inherits(x, "GatingHierarchy") |
    inherits(x, "GatingSet"))) {
    stop("'x' should be a flowFrame, flowSet, GatingHierarchy or GatingSet.")
  }

  return(TRUE)
}

# CYTO_TRANSFORM ---------------------------------------------------------------

#' Apply logicle transformation to channels of flowSet or GatingSet
#'
#' A convenient wrapper for
#' \code{\link[flowCore:estimateLogicle]{estimateLogicle}} and
#' \code{\link[flowCore:transform]{transform}} which coerces samples prior to
#' logicle parameter estimation. For large \code{flowSet} or \code{GatingSet}
#' objects it is recommended that users select a subset of samples for parameter
#' estimation using the \code{select} argument.
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param channels name(s) of the channels to transform. A call will be made to
#'   \code{\link{cyto_fluor_channels}} to get the names of the fluorescent
#'   channels if no channels are supplied.
#' @param select a named list of experiment variables passed to
#'   \code{\link{cyto_select}} to select particular samples to use for logicle
#'   parameter estimation.
#' @param ... additional arguments passed to
#'   \code{\link[flowCore:estimateLogicle]{estimateLogicle}}.
#'
#' @return transformed \code{flowSet} or \code{GatingSet}.
#'
#' @importFrom flowCore estimateLogicle transform flowSet
#' @importFrom flowWorkspace GatingSet
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_transform
#'
#' @export
cyto_transform <- function(x, ...) {
  UseMethod("cyto_transform")
}

#' @rdname cyto_transform
#' @export
cyto_transform.flowSet <- function(x,
                                   channels,
                                   select = NULL,
                                   ...) {

  # Backwards compatibility - transformList will not be available!

  # Select a subset of samples for parameter estimation
  if (!is.null(select)) {
    x <- cyto_select(x, select)
  }

  # Missing channels - use fluorescent channels
  if (missing(channels)) {
    channels <- cyto_fluor_channels(x)
  }

  # Convert channels argument to valid channel names
  channels <- cyto_channels_extract(x, channels)

  # Coerce flowSet to flowFrame for parameter estimation
  fr <- as(x, "flowFrame")

  # Logicle parameter estimates using estimateLogicle
  trans <- estimateLogicle(fr, channels, ...)

  # Apply transformation using flowCore:;transform
  x <- transform(x, trans)

  return(x)
}

#' @rdname cyto_transform
#' @export
cyto_transform.GatingSet <- function(x,
                                     channels,
                                     select = NULL,
                                     ...) {

  # Select a subset of samples for parameter estimation
  if (!is.null(select)) {
    x <- cyto_select(x, select)
  }

  # Missing channels - use fluorescent channels
  if (missing(channels)) {
    channels <- cyto_fluor_channels(x)
  }

  # Convert channels argument to valid channel names
  channels <- cyto_channels_extract(x, channels)

  # Extract data from GatingSet
  fs <- cyto_extract(x, "root")

  # Coerce flowSet to flowFrame
  fr <- cyto_convert(fs, "flowFrame")

  # Add collapsed flowFrame to GatingSet
  gs <- GatingSet(flowSet(fr))

  # Logicle parameter estimates using estimateLogicle
  trans <- estimateLogicle(gs[[1]], channels, ...)

  # Apply transformation using flowCore::transform
  x <- transform(x, trans)

  return(x)
}

# CYTO_TRANSFORM_INVERSE -------------------------------------------------------



# CYTO_TRANSFORM_CONVERT -------------------------------------------------------

#' Convert transformerLists to transformLists and/or extract inverse
#' transformations
#'
#' @param x object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the samples.
#' @param inverse logical indicating whether the returned transformList should
#'   contain the inverse transformations.
#'
#' @return A \code{\link[flowCore:transformList-class]{transformList}} containing
#'   the desired transformations.
#'
#' @importFrom flowCore inverseLogicleTransform transformList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples to flowSet
#' fs <- Activation
#' 
#' # Add fs to GatingSet
#' gs <- GatingSet(fs)
#' 
#' # Get transformList using estimateLogicle
#' trans <- estimateLogicle(fs[[32]], cyto_fluor_channels(fs))
#' 
#' # Get transformList containing inverse transformarions
#' inv <- cyto_transform_convert(trans, inverse = TRUE)
#' 
#' # Convert transformerList into transformList
#' trans <- estimateLogicle(gs[[32]], cyto_fluor_channels(gs))
#' 
#' # Convert transformerList into inverse transformList
#' inv <- cyto_transform_convert(trans, inverse = FALSE)
#' @rdname cyto_transform_convert
#'
#' @export
cyto_transform_convert <- function(x, inverse = FALSE) {
  UseMethod("cyto_transform_convert")
}

#' @rdname cyto_transform_convert
#' @noRd
cyto_transform_convert.default <- function(x, inverse = FALSE) {
  if (is.null(x)) {
    return(x)
  } else if (.all_na(x)) {
    return(x)
  } else {
    warning(
      paste(
        "cyto_transform_convert expects objects of class",
        "transformList or transformerList."
      )
    )
  }
}

#' @rdname cyto_transform_convert
#' @export
cyto_transform_convert.transformList <- function(x, inverse = FALSE) {
  if (inverse) {
    x <- inverseLogicleTransform(x)
  }

  return(x)
}

#' @rdname cyto_transform_convert
#' @export
cyto_transform_convert.transformerList <- function(x, inverse = FALSE) {
  if (inverse) {
    x <- transformList(names(x), lapply(x, `[[`, "inverse"))
  } else {
    x <- transformList(names(x), lapply(x, `[[`, "transform"))
  }

  return(x)
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
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_extract <- function(x, parent = "root", ...) {

  # Throw error for invalid object
  if (!cyto_check(x)) {
    stop(
      paste(
        "'x' should be either a flowFrame, flowSet, GatingHierarchy",
        "or GatingSet."
      )
    )
  }

  # Extract data from GatingHierarchy or GatingSet
  if (any(inherits(x, "GatingHierarchy") |
    inherits(x, "GatingSet"))) {
    x <- getData(x, parent, ...)
  }

  return(x)
}

# CYTO_CONVERT -----------------------------------------------------------------

#' Convert between cytometry objects
#'
#' Automatically removes 'Original' column for coerced objects.
#'
#' @param x \code{flowFrame}, \code{flowSet}, \code{GatingHierarchy},
#'   \code{GatingSet}.
#' @param return either 'flowFrame', 'flowSet', 'GatingHierarchy', 'GatingSet',
#'   coerced 'flowFrame list' or coreced 'flowSet list'. GatingSet and flowSet
#'   objects can also be converted to a 'list of flowFrames'.
#' @param parent name of parent population to extract from GatingHierarchy and
#'   GatingSet objects.
#'
#' @return object specified by 'return' argument.
#'
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace GatingSet
#' @importFrom methods as
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Convert flowSet to 'list of flowFrames'
#' cyto_convert(Activation, "list of flowFrames")
#' 
#' # Convert flowSet to 'flowFrame'
#' cyto_convert(Activation, "flowFrame")
#' 
#' # Convert GatingSet to flowFrame
#' cyto_convert(GatingSet(Activation), "flowFrame", parent = "T Cells")
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_convert
#'
#' @export
cyto_convert <- function(x, ...) {
  UseMethod("cyto_convert")
}

#' @rdname cyto_convert
#' @export
cyto_convert.flowFrame <- function(x,
                                   return = "flowFrame") {
  if (return == "list of flowFrames") {
    return <- "flowFrame list"
  }

  if (return == "flowFrame") {

  } else if (return == "flowFrame list") {
    x <- list(x)
  } else if (return == "flowSet") {
    x <- flowSet(x)
  } else if (return == "flowSet list") {
    x <- list(flowSet(x))
  } else if (return == "GatingSet") {
    x <- flowSet(x)
    x <- GatingSet(x)
  } else if (return == "GatingHierarchy") {
    x <- flowSet(x)
    x <- GatingSet(x)[[1]]
  }

  return(x)
}

#' @rdname cyto_convert
#' @export
cyto_convert.flowSet <- function(x,
                                 return = "flowSet") {
  if (return == "flowSet") {

  } else if (return == "flowFrame") {
    x <- as(x, "flowFrame")
    if ("Original" %in% BiocGenerics::colnames(x)) {
      x <- suppressWarnings(
        x[, -match("Original", BiocGenerics::colnames(x))]
      )
    }
  } else if (return == "flowFrame list") {
    x <- as(x, "flowFrame")
    if ("Original" %in% BiocGenerics::colnames(x)) {
      x <- suppressWarnings(
        x[, -match("Original", BiocGenerics::colnames(x))]
      )
    }
    x <- list(x)
  } else if (return == "list of flowFrames") {
    x <- lapply(seq_len(length(x)), function(y) {
      x[[y]]
    })
    names(x) <- cyto_names(x)
  } else if (return == "flowSet list") {
    x <- list(x)
  } else if (return == "GatingSet") {
    x <- GatingSet(x)
  } else if (return == "GatingHierarchy") {
    x <- as(x, "flowFrame")
    if ("Original" %in% BiocGenerics::colnames(x)) {
      x <- suppressWarnings(
        x[, -match("Original", BiocGenerics::colnames(x))]
      )
    }
    x <- flowSet(x)
    x <- GatingSet(x)[[1]]
  }

  return(x)
}

#' @rdname cyto_convert
#' @export
cyto_convert.GatingHierarchy <- function(x,
                                         parent = "root",
                                         return = "GatingHierarchy") {
  if (return == "GatingHierarchy") {

  } else if (return == "flowFrame") {
    x <- cyto_extract(x, parent)
  } else if (return == "flowFrame list") {
    x <- list(cyto_extract(x, parent))
  } else if (return == "flowSet") {
    x <- flowSet(cyto_extract(x, parent))
  } else if (return == "flowSet list") {
    x <- list(flowSet(cyto_extract(x, parent)))
  }

  return(x)
}

#' @rdname cyto_convert
#' @export
cyto_convert.GatingSet <- function(x,
                                   parent = "root",
                                   return = "GatingSet") {
  if (return == "GatingSet") {

  } else if (return == "flowFrame") {
    x <- cyto_convert(cyto_extract(x, parent), "flowFrame")
  } else if (return == "flowFrame list") {
    x <- list(cyto_convert(cyto_extract(x, parent), "flowFrame"))
  } else if (return == "list of flowFrames") {
    x <- lapply(seq(1, length(x)), function(z) {
      cyto_extract(x[[z]], parent)
    })
    names(x) <- cyto_names(x)
  } else if (return == "flowSet") {
    x <- cyto_extract(x, parent)
  } else if (return == "flowSet list") {
    x <- list(cyto_extract(x, parent))
  } else if (return == "GatingHierachy") {
    x <- flowSet(cyto_convert(cyto_extract(x, parent), "flowFrame"))
    x <- GatingSet(x)[[1]]
  }

  return(x)
}

# CYTO_FILTER ------------------------------------------------------------------

#' Filter samples based on experiment variables
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param ... tidyverse-style subsetting using comma separated logical
#'   predicates based on experimental variables stored in \code{cyto_details(x)}. See
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
#' cyto_details(Activation)
#' 
#' # Select Stim-C samples with 0 and 0.5 OVA concentrations
#' fs <- cyto_filter(
#'   Activation,
#'   Treatment == "Stim-C",
#'   OVAConc %in% c(0, 0.5)
#' )
#' 
#' # Select Stim-A and Stim-C treatment groups
#' fs <- cyto_filter(
#'   Activation,
#'   Treatment %in% c("Stim-A", "Stim-C")
#' )
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_filter <- function(x, ...) {

  # Check class of x
  if (!any(inherits(x, "flowSet") |
    inherits(x, "GatingSet"))) {
    stop("'x' should be an object of class flowSet or GatingSet.")
  }

  # Extract experiment details
  pd <- cyto_details(x)

  # Perform filtering on pd to pull out samples
  pd_filter <- filter(pd, ...)

  # Get indices for selected samples
  if (nrow(pd_filter) == 0) {
    message("No samples match the filtering criteria. Returning all samples.")
    ind <- seq_len(length(x))
  } else {
    ind <- match(pd_filter[, "name"], pd[, "name"])
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
#' cyto_details(Activation)
#' 
#' # Select Stim-C samples with 0 and 0.5 OVA concentrations
#' fs <- cyto_select(Activation,
#'   Treatment = "Stim-C",
#'   OVAConc = c(0, 0.5)
#' )
#' 
#' # Select Stim-A and Stim-C treatment groups
#' fs <- cyto_select(
#'   Activation,
#'   list("Treatment" = c("Stim-A", "Stim-C"))
#' )
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_select <- function(x, ...) {

  # Check class of x
  if (!any(inherits(x, "flowSet") |
    inherits(x, "GatingSet"))) {
    stop("'x' should be an object of class flowSet or GatingSet.")
  }

  # Pull down ... arguments to list
  args <- list(...)

  # ... is already a named list of arguments
  if (class(args[[1]]) == "list") {
    args <- args[[1]]
  }

  # Extract experiment details
  pd <- cyto_details(x)

  # Check that all varaibles are valid
  if (!all(names(args) %in% colnames(pd))) {
    lapply(names(args), function(y) {
      if (!y %in% names(pd)) {
        stop(paste(y, "is not a valid variable in cyto_details(x)."))
      }
    })
  }

  # Check that all variable levels at least exist in pd
  lapply(names(args), function(z) {
    var <- factor(pd[, z], exclude = NULL) # keep <NA> as factor level
    lvls <- levels(var)
    # some variable levels do not exist in pd
    if (!all(args[[z]] %in% lvls)) {
      lapply(args[[z]], function(v) {
        if (!v %in% lvls) {
          stop(paste0(v, " is not a valid level for ", z, "!"))
        }
      })
    }
  })

  # Get filtered pd
  pd_filter <- pd
  lapply(names(args), function(y) {
    ind <- which(pd_filter[, y] %in% args[[y]])
    # No filtering if variable level is missing
    if (length(ind) != 0) {
      pd_filter <<- pd_filter[ind, ]
    }
  })

  # Get indices for selected samples
  ind <- match(pd_filter[, "name"], pd[, "name"])

  return(x[ind])
}

# CYTO_GROUP_BY ----------------------------------------------------------------

#' Group a flowSet or GatingSet by experiment variables
#'
#' @param x an object of class \code{flowSet} or \code{GatingSet}.
#' @param parent name of the parent population to extract from GatingSet object.
#' @param group_by names of cyto_details variables to use for merging. Set to "all" to
#'   merge all samples in \code{x}.
#'
#' @return a named list of \code{flowSet} or \code{GatingSet} objects
#'   respectively.
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Group flowSet by Treatment
#' cyto_group_by(Activation, "Treatment")
#' 
#' # Group GatingSet by Treatment and OVAConc
#' gs <- GatingSet(Activation)
#' cyto_group_by(gs, c("Treatment", "OVAConc"))
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_group_by
#'
#' @export
cyto_group_by <- function(x,
                          group_by = "all") {

  # Check class of x
  if (!any(inherits(x, "flowSet") | inherits(x, "GatingSet"))) {
    stop("'x' should be an object of class flowSet or GatingSet.")
  }

  # Extract experiment details
  pd <- cyto_details(x)

  # Repace any NA with "NA" to avoid missing rows
  if (any(is.na(pd))) {
    pd[is.na(pd)] <- "NA"
  }

  # Extract sample names
  nms <- cyto_names(x)

  # Check group_by
  if (group_by[1] != "all" & !all(group_by %in% colnames(pd))) {
    lapply(group_by, function(y) {
      if (!y %in% colnames(pd)) {
        stop(paste0(y, " is not a valid variable for this ", class(x), "."))
      }
    })
  }

  # Split pd based on group_by into a named list
  if (group_by[1] == "all") {
    pd_split <- list("all" = pd)
  } else {
    pd_split <- split(pd, pd[, group_by],
      sep = " ",
      lex.order = TRUE,
      drop = TRUE
    )
  }

  # Replace each element of pd_split with matching samples
  x_list <- lapply(seq_len(length(pd_split)), function(z) {
    ind <- match(pd_split[[z]][, "name"], cyto_names(x))
    x[ind]
  })
  names(x_list) <- names(pd_split)

  return(x_list)
}

# CYTO_SAMPLE ------------------------------------------------------------------

#' Sample a flowFrame or flowSet
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#'   \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param display numeric [0,1] indicating the percentage of events to keep.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
#'
#' @return \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#'   \code{\link[flowCore:flowSet-class]{flowSet}} restricted to \code{display}
#'   percentage events.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter Subset fsApply
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Constrict first sample by 50%
#' cyto_sample(fs[[1]], 0.5)
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_sample
#'
#' @export
cyto_sample <- function(x, ...) {
  UseMethod("cyto_sample")
}

#' @rdname cyto_sample
#' @export
cyto_sample.flowFrame <- function(x,
                                  display = 1,
                                  seed = NULL) {

  # Do nothing if no sampling required
  if (display != 1) {
    # Number of events
    events <- nrow(x)

    # Size
    size <- display * events

    # Set seed
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Sample
    smp <- sampleFilter(size = size)
    x <- Subset(x, smp)
  }

  return(x)
}

#' @rdname cyto_sample
#' @export
cyto_sample.flowSet <- function(x,
                                display = 1,
                                seed = NULL) {
  fsApply(x, cyto_sample, display = display, seed = seed)
}

#' @rdname cyto_sample
#' @export
cyto_sample.list <- function(x,
                             display = 1,
                             seed = NULL) {
  lapply(x, function(z) {
    cyto_sample(z, display, seed)
  })
}

# CYTO_MARKERS -----------------------------------------------------------------

#' Assign marker names to flowFrame or flowSet
#'
#' \code{cyto_markers} opens an editable table containing a list of channels and
#' markers for a \code{flowFrame}, \code{flowSet}, \code{GatingHierarchy} or
#' \code{GatingSet}. Users can edit the \code{name} or \code{desc} columns with
#' updated channel names or marker names respectively. These entries will be
#' updated in the \code{x} upon closing the window and saved to a
#' "Experiment-markers.csv" file for future use.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
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

  # check class of x
  cyto_check(x)

  # flowFrame
  if (inherits(x, "flowFrame")) {

    # Extract details of parameters
    pd <- cyto_details(parameters(x))

    # flowSet
  } else if (inherits(x, "flowSet")) {

    # Extract details of parameters
    pd <- cyto_details(parameters(x[[1]]))

    # GatingHierachy
  } else if (inherits(x, "GatingHierarchy")) {
    fr <- cyto_extract(x, "root")
    pd <- cyto_details(parameters(fr))

    # GatingSet
  } else if (inherits(x, "GatingSet")) {
    fr <- cyto_extract(x, "root")[[1]]
    pd <- cyto_details(parameters(fr))
  }

  # file missing
  if (is.null(file)) {
    if (length(grep("Experiment-Markers.csv", list.files())) != 0) {
      message("Experiment-Markers.csv found in working directory.")
      dt <- read.csv(list.files()[grep("Experiment-Markers.csv", list.files())],
        header = TRUE, stringsAsFactors = FALSE
      )
    } else {

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
    write.csv(dt,
      c(file, list.files())[grep(
        "Experiment-Markers.csv",
        c(file, list.files())
      )[1]],
      row.names = FALSE
    )
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

#' Interactively edit cyto_details for a flowSet or GatingSet
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param file name of csv file containing experimental information.
#'
#' @return NULL and update cyto_details for the \code{flowSet} or \code{GatingSet}.
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
#' # Edit cyto_details in table editor
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
    if (length(grep("Experiment-Details.csv", list.files())) != 0) {
      message("Experiment-Details.csv found in working directory.")
      pd <- read.csv(list.files()[grep("Experiment-Details.csv", list.files())],
        header = TRUE, stringsAsFactors = FALSE
      )
    } else {

      # Extract cyto_details
      pd <- cyto_details(cyto)
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

    # Update cyto_details
    cyto_details(cyto) <- pd
  }

  # Edit cyto_details
  pd <- suppressWarnings(edit(pd))
  rownames(pd) <- pd$name

  # Update cyto_details
  cyto_details(cyto) <- pd

  # Write result to csv file
  if (length(grep("Experiment-Details.csv", c(file, list.files()))) != 0) {
    write.csv(pd,
      c(file, list.files())[grep(
        "Experiment-Details.csv",
        c(file, list.files())
      )[1]],
      row.names = FALSE
    )
  } else {
    write.csv(cyto_details(cyto), paste0(
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
#' @importFrom tools file_ext
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
#' write.csv(
#'   spill,
#'   "Spillover-Matrix.csv"
#' )
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
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_compensate
#'
#' @export
cyto_compensate <- function(x, ...) {
  UseMethod("cyto_compensate")
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowFrame <- function(x,
                                      spillover = NULL,
                                      select = 1) {

  # Spillover matrix supplied - matrix, data.frame or csv file
  if(!is.null(spillover)){
    
    # spillover is a character string containing name of csv file
    if(inherits(spillover, "character")){
      # No file extension
      if(file_ext(spillover) == ""){
        spillover <- paste0(spillover, ".csv")
      }
      
      # Check working directory for csv file
      if(getOption("CytoRSuite_wd_check")){
        if(!file_wd_check(spillover)){
          stop(paste(spillover, "does not exist in this working directory."))
        }
      }
      spill <- read.csv(spillover, header = TRUE, row.names = 1)
      colnames(spill) <- rownames(spill)
    
      # column/row names must be valid channels
      if(!all(rownames(spill) %in% BiocGenerics::colnames(x)) |
         !all(rownames(spill) %in% BiocGenerics::colnames(x))){
        stop(
          paste("'spillover' must have valid fluorescent channels as rownames",
                "and colnames.")
        )
      }
      
    # spillover is a matrix/data.frame
    }else if(inherits(spillover, "matrix") |
             inherits(spillover, "data.frame")){
      
      # column names must be valid channels (rownames not essential)
      if(!all(colnames(spillover) %in% BiocGenerics::colnames(x))){
        stop("'spillover' must have valid fluorescent channels as colnames.")
      }else{
        spill <- spillover
      }
    
    }
  
  # Extract spillover matrix directly from x  
  }else if (is.null(spillover)) {
    spill <- x@description$SPILL
  }

  # Apply compensation
  flowCore::compensate(x, spill)
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowSet <- function(x,
                                    spillover = NULL,
                                    select = 1) {

  # Spillover matrix supplied - matrix, data.frame or csv file
  if(!is.null(spillover)){
    
    # spillover is a character string containing name of csv file
    if(inherits(spillover, "character")){
      # No file extension
      if(file_ext(spillover) == ""){
        spillover <- paste0(spillover, ".csv")
      }
      
      # Check working directory for csv file
      if(getOption("CytoRSuite_wd_check")){
        if(!file_wd_check(spillover)){
          stop(paste(spillover, "does not exist in this working directory."))
        }
      }
      spill <- read.csv(spillover, header = TRUE, row.names = 1)
      colnames(spill) <- rownames(spill)
      
      # column/row names must be valid channels
      if(!all(rownames(spill) %in% BiocGenerics::colnames(x)) |
         !all(rownames(spill) %in% BiocGenerics::colnames(x))){
        stop(
          paste("'spillover' must have valid fluorescent channels as rownames",
                "and colnames.")
        )
      }
      
      # Convert spill into a named list
      spill <- rep(list(spill), length(x))
      names(spill) <- cyto_names(x)
      
    # spillover is a matrix/data.frame
    }else if(inherits(spillover, "matrix") |
             inherits(spillover, "data.frame")){
      
      # column names must be valid channels (rownames not essential)
      if(!all(colnames(spillover) %in% BiocGenerics::colnames(x))){
        stop("'spillover' must have valid fluorescent channels as colnames.")
      }else{
        spill <- spillover
      }
      
      # Convert spill into a named list
      spill <- rep(list(spill), length(x))
      names(spill) <- cyto_names(x)
      
    }
    
    # Extract spillover matrix directly from x  
  }else if (is.null(spillover)) {
    if (!is.null(select)) {
      spill <- rep(list(x[[select]]@description$SPILL), length(x))
      names(spill) <- cyto_names(x)
    } else {
      spill <- lapply(cyto_names(x), function(y) {
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
                                      select = 1) {
  
  # Extract flowSet
  fs <- cyto_extract(x, parent = "root")
  
  # Spillover matrix supplied - matrix, data.frame or csv file
  if(!is.null(spillover)){
    
    # spillover is a character string containing name of csv file
    if(inherits(spillover, "character")){
      # No file extension
      if(file_ext(spillover) == ""){
        spillover <- paste0(spillover, ".csv")
      }
      
      # Check working directory for csv file
      if(getOption("CytoRSuite_wd_check")){
        if(!file_wd_check(spillover)){
          stop(paste(spillover, "does not exist in this working directory."))
        }
      }
      spill <- read.csv(spillover, header = TRUE, row.names = 1)
      colnames(spill) <- rownames(spill)
      
      # column/row names must be valid channels
      if(!all(rownames(spill) %in% BiocGenerics::colnames(fs)) |
         !all(rownames(spill) %in% BiocGenerics::colnames(fs))){
        stop(
          paste("'spillover' must have valid fluorescent channels as rownames",
                "and colnames.")
        )
      }
      
      # Convert spill into a named list
      spill <- rep(list(spill), length(fs))
      names(spill) <- cyto_names(fs)
      
      # spillover is a matrix/data.frame
    }else if(inherits(spillover, "matrix") |
             inherits(spillover, "data.frame")){
      
      # column names must be valid channels (rownames not essential)
      if(!all(colnames(spillover) %in% BiocGenerics::colnames(fs))){
        stop("'spillover' must have valid fluorescent channels as colnames.")
      }else{
        spill <- spillover
      }
      
      # Convert spill into a named list
      spill <- rep(list(spill), length(fs))
      names(spill) <- cyto_names(fs)
      
    }
    
    # Extract spillover matrix directly from fs  
  }else if (is.null(spillover)) {
    if (!is.null(select)) {
      spill <- rep(list(fs[[select]]@description$SPILL), length(fs))
      names(spill) <- cyto_names(fs)
    } else {
      spill <- lapply(cyto_names(fs), function(y) {
        fs[[y]]@description$SPILL
      })
    }
  }

  # Apply compensation
  flowWorkspace::compensate(x, spill)
}
