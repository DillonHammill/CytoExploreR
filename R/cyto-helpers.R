## CYTO_LOAD -------------------------------------------------------------------

#' Load .fcs files into ncdfFlowSet
#'
#' \code{cyto_load} is a convenient wrapper around
#' \code{\link[base:list.files]{list.files}} and
#' \code{\link[flowWorkspace:load_cytoset_from_fcs]{load_cytoset_from_fcs}}
#' which makes it easy to load .fcs files into a cytoset. \code{cyto_load} is
#' also a wrapper around \code{\link[flowWorkspace:save_gs]{load_gs}} to load
#' saved GatingSet objects.
#'
#' @param path points to the location of the .fcs files to read in. Preferably
#'   the name of folder in current working directory.
#' @param select vector of file names to select when loading files, set to NULL
#'   be default to select all files in the specified directory.
#' @param exclude vector of file names to exclude when loading files, set to
#'   NULL by default to load all files in the specified directory.
#' @param sort logical indicating whether attempts should be made to sort the
#'   files by name prior to loading, set to \code{TRUE} by default.
#' @param barcode logical indicating whether the flowFrames should be barcoded
#'   using \code{cyto_barcode}, set to FALSE by default.
#' @param restrict logical indicating whether unassigned channels should be
#'   dropped from the returned cytoset, set to FALSE by default. See
#'   \code{\link{cyto_channels_restrict}}.
#' @param ... additional arguments passed to \code{load_cytoset_from_fcs}.
#'
#' @return object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @importFrom flowCore identifier identifier<- parameters
#' @importFrom flowWorkspace load_gs load_cytoset_from_fcs pData
#' @importFrom gtools mixedsort
#' @importFrom tools file_ext
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Get path to Activation .fcs files in CytoExploreRData
#' datadir <- system.file("extdata", package = "CytoExploreRData")
#' path <- paste0(datadir, "/Activation")
#'
#' # Load in .fcs files into cytoset
#' cs <- cyto_load(path)
#'
#' # cs is a cytoset
#' class(cs)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_load <- function(path = ".",
                      select = NULL,
                      exclude = NULL,
                      sort = TRUE,
                      barcode = FALSE,
                      restrict = FALSE, ...) {

  # FILE PATHS
  files <- list.files(path, full.names = TRUE)

  # SELECT
  if(!is.null(select)){
    lapply(select, function(z){
      if(any(grepl(z, files, ignore.case = TRUE))){
        files <<- files[which(grepl(z, files, ignore.case = TRUE))]
      }
    })
  }
  
  # EXCLUDE
  if(!is.null(exclude)){
    lapply(exclude, function(z){
      if(any(grepl(z, files, ignore.case = TRUE))){
        files <<- files[-which(grepl(z, files, ignore.case = TRUE))]
      }
    })
  }
  
  # SORTED FILE PATHS
  if (sort) {
    files <- mixedsort(files)
  }

  # SAVED GATINGSET
  if (all(c("pb", "rds") %in% file_ext(files))) {
    # LOAD GATINGSET
    x <- load_gs(path = path)
    # FCS FILES
  } else {
    # NCDFFLOWSET
    x <- load_cytoset_from_fcs(files = files, ...)

    # CORRECT GUID SLOTS
    nms <- cyto_names(x)
    lapply(seq_len(length(nms)), function(z) {
      suppressMessages(identifier(x[[z]]) <<- nms[z])
    })

    # BARCODING
    if (barcode) {
      x <- cyto_barcode(x)
    }
    
    # CHANNEL RESTRICTION
    if(restrict){
      x <- cyto_channels_restrict(x)
    }
    
  }

  # RETURN NCDFFLOWSET
  return(x)
}

## CYTO_CLEAN ------------------------------------------------------------------

#' Apply flowAI anomaly detection to clean cytometry data
#'
#' @param x object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}. The \code{root} node extracted
#'   when a \code{GatingSet} or \code{GatingHierachy} is supplied.
#' @param ... additional arguments passed to
#'   \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}.
#'
#' @importFrom flowAI flow_auto_qc
#' @importFrom flowWorkspace gs_cyto_data cytoset_to_flowSet flowSet_to_cytoset
#' @importFrom utils capture.output
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Clean Activation flowSet
#' fs <- cyto_clean(fs)
#'
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#'
#' # Clean Activation GatingSet
#' gs <- cyto_clean(gs)
#'
#' @references Monaco,G. et al. (2016) flowAI: automatic and interactive anomaly
#'   discerning tools for flow cytometry data. Bioinformatics. 2016 Aug
#'   15;32(16):2473-80.
#'   
#' @seealso \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}
#'
#' @export
cyto_clean <- function(x, ...){
  
  # GATINGSET/GATINGHIERARCHY
  if(is(x, "GatingSet") | is(x, "GatingHierarchy")){
    # PARENT
    parent <- cyto_nodes(x, path = "auto")[1]
    # EXTRACT DATA
    cyto_data <- cyto_extract(x, parent)
    # flowAI REQUIRES FLOWSET
    if(is(cyto_data, "cytoset")){
      cyto_data <- cytoset_to_flowSet(cyto_data) # REMOVE
    }
    # CLEAN DATA
    invisible(capture.output(cyto_data <-  flow_auto_qc(cyto_data,
                                                        html_report = FALSE,
                                                        mini_report = FALSE,
                                                        fcs_QC = FALSE,
                                                        folder_results = FALSE,
                                                        ...)))
    # RETURN CYTOSET
    if(is(cyto_data, "flowSet")){
      cyto_data <- flowSet_to_cytoset(cyto_data)
    }
    # REPLACE DATA
    gs_cyto_data(x)<- cyto_data
  }else{
    # FLOWSET REQUIRED
    if(is(x, "cytoset")){
      x <- cytoset_to_flowSet(x)
    }
    invisible(capture.output(x <- flow_auto_qc(x, 
                                               html_report = FALSE,
                                               mini_report = FALSE,
                                               fcs_QC = FALSE,
                                               folder_results = FALSE,
                                               ...)))
    # RETURN CYTOSET
    if(is(x, "flowSet")){
      x <- flowSet_to_cytoset(x)
    }
  }
  return(x)
}

## CYTO_SETUP ------------------------------------------------------------------

#' Load.fcs files into GatingSet and annotate with experiment details
#'
#' \code{cyto_setup} takes care of all the data loading and annotation steps to
#' prepare your cytometry data for downstream analyses. The .fcs files are first
#' read into a \code{\link[flowWorkspace:cytoset]{cytoset}} using
#' \code{\link{cyto_load}} which is then added to a
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' Calls are then made to \code{\link{cyto_markers_edit}} and
#' \code{\link{cyto_details_edit}} to update the GatingSet with the details of
#' the experiment. These details can be modified later with additional calls to
#' \code{\link{cyto_markers_edit}} and/or \code{\link{cyto_details_edit}}.
#'
#' Through the \code{clean} argument, the data can then be optionally cleaned
#' using \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}} to automatically remove
#' anomalies in the recorded data. 
#' 
#' Users can optionally provide a name for a
#' gatingTemplate csv file which will be created if necessary and assigned as
#' the active gatingTemplate.
#'
#' @param path points to the location of the .fcs files to read in (e.g. name of
#'   a folder in current working directory).
#' @param gatingTemplate name of a gatingTemplate csv file to be used for gate
#'   saving.
#' @param restrict logical indicating whether unassigned channels should be
#'   dropped from the returned cytoset, set to FALSE by default.  See
#'   \code{\link{cyto_channels_restrict}}.
#' @param clean logical indicating whether the loaded data should be cleaned
#'   using \code{cyto_clean}, set to FALSE by default.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:load_cytoset_from_fcs]{load_cytoset_from_fcs}}.
#'
#' @return object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @importFrom flowWorkspace GatingSet
#' @importFrom tools file_ext
#' @importFrom methods is
#'
#' @examples
#'
#' \dontrun{
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Get path to Activation .fcs files in CytoExploreRData
#' datadir <- system.file("extdata", package = "CytoExploreRData")
#' path <- paste0(datadir, "/Activation")
#'
#' # Load in .fcs files into an annotated GatingSet
#' gs <- cyto_setup(path)
#'
#' # Markers have been assigned
#' cyto_extract(gs, "root")[[1]]
#'
#' # Experiment details have been updated
#' cyto_details(gs)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' \code{\link{cyto_load}}
#' \code{\link{cyto_markers_edit}}
#' \code{\link{cyto_details_edit}}
#' \code{\link{cyto_channels_restrict}}
#' \code{\link{cyto_clean}}
#' \code{\link{cyto_gatingTemplate_select}}
#' \code{\link{cyto_gatingTemplate_create}}
#'
#' @export
cyto_setup <- function(path = ".",
                       gatingTemplate = NULL,
                       restrict = FALSE,
                       clean = FALSE, ...) {

  # CYTOSET/GATINGSET
  message("Loading FCS files into a GatingSet...")
  x <- cyto_load(path = path, restrict = FALSE, ...)

  # MARKERS
  message("Assigning markers to channels...")
  x <- cyto_markers_edit(x)

  # EXPERIMENT DETAILS
  message("Updating experiment details...")
  x <- cyto_details_edit(x)
  
  # FLOWSET LOADED
  if (is(x, "flowSet")) {
    # RESTRICT CHANNELS
    if(restrict){
      message("Removing unassigned channels...")
      x <- cyto_channels_restrict(x)
    }
    # CLEAN DATA
    if(clean){
      message("Cleaning data to remove anomalies...")
      x <- cyto_clean(x)
    }
    # GATINGSET
    x <- GatingSet(x)
  }
  
  # gatingtemplate
  if (!is.null(gatingTemplate)) {

    # FILE EXTENSION
    if (.empty(file_ext(gatingTemplate))) {
      gatingTemplate <- paste0(gatingTemplate, ".csv")
    }

    # ACTIVE GATINGTEMPLATE
    message(paste("Setting", gatingTemplate, "as the active gatingTemplate..."))
    cyto_gatingTemplate_select(gatingTemplate)

    # CREATE GATINGTEMPLATE
    if (.all_na(match(gatingTemplate, list.files()))) {
      message(paste("Creating", gatingTemplate, "..."))
      cyto_gatingTemplate_create(gatingTemplate)
    }
  }

  # SETUP COMPLETE
  message("Done!")
  
  # RETURN GATINGSET
  return(x)
  
}

## CYTO_DETAILS ----------------------------------------------------------------

#' Extract experiment details
#'
#' Simply an autocomplete-friendly wrapper around
#' \code{\link[flowWorkspace:pData-methods]{pData}}. A call is made to
#' \code{\link{cyto_names}} if a
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} is supplied.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return experiment details as data.frame.
#'
#' @importFrom flowWorkspace pData
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Experiment details
#' cyto_details(fs)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_details <- function(x) {

  # Return identifier for flowFrame
  if (is(x, "flowFrame")) {
    return(cyto_names(x))
    # Return experiment details for other objects
  } else {
    # Fix AsIs for name column
    pd <- pData(x)
    #pd$name <- factor(pd$name, levels = pd$name)
    return(pd)
  }
}

#' @importFrom flowWorkspace pData<-
#' @noRd
"cyto_details<-" <- `pData<-`

#' @noRd
#' @export
"cyto_details<-"

## CYTO_NAMES ------------------------------------------------------------------

#' Extract sample names
#'
#' Simply a convenient and autocomplete-friendly wrapper around
#' \code{\link[flowCore:identifier-methods]{identifier}}
#' \code{\link[flowWorkspace:sampleNames]{sampleNames}} to extract the sample
#' names from flowFrame, flowSet GatingHierarchy or GatingSet. Anonymous
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} identifiers will be
#' converted to \code{"Combined Events"}.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return names associated with the supplied object.
#'
#' @importFrom flowCore identifier
#' @importFrom flowWorkspace sampleNames
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#'
#' # flowFrame
#' cyto_names(fs[[1]])
#'
#' # flowSet
#' cyto_names(fs)
#'
#' # GatingHierarchy
#' cyto_names(gs[[1]])
#'
#' # GatingSet
#' cyto_names(gs)
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_names
#'
#' @export
cyto_names <- function(x) {
  UseMethod("cyto_names")
}

#' @rdname cyto_names
#' @export
cyto_names.flowFrame <- function(x) {
  nm <- identifier(x)
  if (nm == "anonymous") {
    nm <- "Combined Events"
  }
  return(nm)
}

#' @rdname cyto_names
#' @export
cyto_names.flowSet <- function(x) {
  sampleNames(x)
}

#' @rdname cyto_names
#' @export
cyto_names.GatingHierarchy <- function(x) {
  sampleNames(x)
}

#' @rdname cyto_names
#' @export
cyto_names.GatingSet <- function(x) {
  sampleNames(x)
}

#' @rdname cyto_names
#' @export
cyto_names.list <- function(x) {
  LAPPLY(x, "cyto_names")
}

## CYTO_CHECK ------------------------------------------------------------------

#' Check a flowFrame, flowSet, GatingHierarchy or GatingSet has been supplied
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} to be checked.
#'
#' @return TRUE or FALSE if object meets this class criteria.
#'
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Valid object
#' cyto_check(Activation)
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_check <- function(x) {

  # Check for valid class of object
  if (!any(is(x, "flowFrame") |
    is(x, "flowSet") |
    is(x, "GatingHierarchy") |
    is(x, "GatingSet"))) {
    stop("'x' should be a flowFrame, flowSet, GatingHierarchy or GatingSet.")
  }

  return(TRUE)
}

## CYTO_TRANSFORM --------------------------------------------------------------

#' Apply Transformations to Cytometry Data
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} containing the
#'   transformation definitions to apply to \code{x}.
#' @param type type of transformation to apply when no trans object is supplied,
#'   options include \code{"log"}, \code{"arcsinh"}, \code{"biex"} and
#'   \code{"logicle"}.
#' @param channels names of the channels to transform, set to all the
#'   fluorescent channels by default. Only required when no \code{trans} object
#'   is supplied.
#' @param parent name of the parent population of \code{GatingHierarchy} or
#'   \code{GatingSet} objects used to visualise the transformations.
#' @param select list of selection criteria passed to \code{\link{cyto_select}}
#'   to select a subset of samples for visualising the transformations.
#' @param inverse logical indicating whether the inverse transformations should
#'   be applied. Currently only supported for \code{flowFrame} and
#'   \code{flowSet} objects.
#' @param plot logical indicating whether the result of the transformations
#'   should be plotted using \code{\link{cyto_plot}}.
#' @param popup logical indicating whether plots should be constructed in a
#'   popup window, set to FALSE by default.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control
#'   over axes limits can be obtained by altering the \code{xlim} and
#'   \code{ylim} arguments.
#' @param ... additional arguments passed to \code{\link{cyto_transformer_log}},
#'   \code{\link{cyto_transformer_arcsinh}}, \code{\link{cyto_transformer_biex}}
#'   or \code{\link{cyto_transformer_logicle}}, when no \code{trans} object is
#'   supplied.
#'
#' @return object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet} with transformations applied.
#'
#' @importFrom flowCore transform
#' @importFrom grDevices n2mfrow
#' @importFrom graphics par
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Automatically transform flowSet
#' fs_trans <- cyto_transform(fs, type = "arcsinh")
#'
#' # Manually construct & apply transformations
#' trans <- cyto_transformer_biex(fs)
#' fs_trans <- cyto_transform(fs, trans)
#'
#' # Add fs to GatingSet
#' gs <- GatingSet(fs)
#'
#' # Automatically transform GatingSet
#' gs_trans <- cyto_transform(gs, type = "logicle")
#'
#' # Manually construct & apply transformations
#' trans <- cyto_transformer_logicle(gs)
#' gs_trans <- cyto_transform(gs, trans)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_transformer_log}}
#' @seealso \code{\link{cyto_transformer_arcsinh}}
#' @seealso \code{\link{cyto_transformer_biex}}
#' @seealso \code{\link{cyto_transformer_logicle}}
#' @seealso \code{\link{cyto_transformer_combine}}
#'
#' @rdname cyto_transform
#'
#' @export
cyto_transform <- function(x, trans = NULL, ...) {
  UseMethod("cyto_transform", trans)
}

#' @rdname cyto_transform
#' @export
cyto_transform.default <- function(x,
                                   trans = NULL,
                                   type = "logicle",
                                   channels = NULL,
                                   parent = "root",
                                   select = NULL,
                                   inverse = FALSE,
                                   plot = TRUE,
                                   popup = FALSE,
                                   axes_limits = "machine",
                                   ...) {

  # No transformations supplied - automatically obtain transform definitions
  if (is.null(trans) | .all_na(trans)) {

    # Message not recommended to auto-transform flowFrame/flowSet objects
    if (is(x, "flowFrame") | is(x, "flowSet")) {
      message(paste(
        "Automatically transforming flowFrame/flowSet objects",
        "is not recommended as transformation definitions will be lost."
      ))
    }

    # Dispatch based on type argument to get TransformerList
    if (type == "log") {
      transformer_list <- cyto_transformer_log(x,
        channels = channels,
        parent = parent,
        select = select,
        plot = FALSE, 
        ...
      )
    } else if (type == "arcsinh") {
      transformer_list <- cyto_transformer_arcsinh(x,
        channels = channels,
        parent = parent,
        select = select,
        plot = FALSE, ...
      )
    } else if (type == "biex") {
      transformer_list <- cyto_transformer_biex(x,
        channels = channels,
        parent = parent,
        select = select,
        plot = FALSE, ...
      )
    } else if (type == "logicle") {
      transformer_list <- cyto_transformer_logicle(x,
        channels = channels,
        parent = parent,
        select = select,
        plot = FALSE, ...
      )
    }
  }

  # Apply transformations
  if (is(x, "flowFrame") |
    is(x, "flowSet")) {

    # Extract transformations from transformerList to transformList
    transform_list <- cyto_transform_extract(transformer_list,
      inverse = inverse
    )

    # Apply transformations
    x <- suppressMessages(transform(x, transform_list))
  } else if (is(x, "GatingHierarchy") |
    is(x, "GatingSet")) {

    # Inverse transformations not yet supported
    if (inverse == TRUE) {
      stop(paste(
        "Inverse transformations are not yet supported for",
        "GatingHierarchy/GatingSet objects."
      ))
    }

    # Apply transformations
    x <- suppressMessages(transform(x, transformer_list))
  }

  # Construct the plots
  if (plot == TRUE) {

    # Plot if space sufficient space
    tryCatch(
      {

        # Pull out flowFrame/flowSet to plot
        cyto_data <- cyto_extract(x, parent)

        # Convert to flowFrame for plotting
        cyto_data <- cyto_convert(cyto_data, "flowFrame")

        # Channels
        channels <- names(transformer_list)

        # Old graphics parameters
        old_pars <- .par("mfrow")
        on.exit(par(old_pars))

        # Set up plotting area
        cyto_plot_new(popup = popup)
        n <- length(channels)
        cyto_plot_layout(
          n2mfrow(n)[1],
          n2mfrow(n)[2]
        )

        # Generate plot for each channel
        lapply(channels, function(chan) {
          if (inverse == FALSE) {
            cyto_plot(cyto_data,
              channels = chan,
              axes_trans = transformer_list,
              title = NA,
              axes_limits = axes_limits
            )
          } else if (inverse == TRUE) {
            cyto_plot(cyto_data,
              channels = chan,
              title = NA,
              axes_limits = axes_limits
            )
          }
        })
      },
      error = function(e) {
        message("Insufficient plotting space, transformations have been applied.")
      }
    )
  }

  # Return transformed data
  return(x)
}

#' @rdname cyto_transform
#' @export
cyto_transform.transformList <- function(x,
                                         trans = NULL,
                                         plot = TRUE,
                                         popup = FALSE,
                                         axes_limits = "machine",
                                         ...) {

  # Added for backwards compatibility - flowFrame/flowSet objects only
  if (is(x, "GatingHierarchy") |
    is(x, "GatingSet")) {
    stop(paste(
      "GatingHierarchy and GatingSet objects require transformerList",
      "objects to apply transformations."
    ))
  }

  # Apply transformations to flowFrame/flowSet
  if (is(x, "flowFrame") |
    is(x, "flowSet")) {

    # Transformations applied as is - allow for inverse transformList
    x <- suppressMessages(transform(x, trans))
  }

  # Construct plots
  if (plot == TRUE) {
    # Plot if sufficient space
    tryCatch(
      {

        # Pull out flowFrame/flowSet to plot
        cyto_data <- cyto_extract(x)

        # Convert to flowFrame for plotting
        cyto_data <- cyto_convert(cyto_data, "flowFrame")

        # Channels
        channels <- names(trans@transforms)

        # Old graphics parameters
        old_pars <- .par("mfrow")
        on.exit(par(old_pars))

        # Set up plotting area
        cyto_plot_new(popup = popup)
        n <- length(channels)
        cyto_plot_layout(
          n2mfrow(n)[1],
          n2mfrow(n)[2]
        )

        # Generate plot for each channel - axes will not be transformed correctly
        lapply(channels, function(chan) {
          cyto_plot(cyto_data,
            channels = chan,
            title = NA,
            axes_limits = axes_limits
          )
        })
      },
      error = function(e) {
        message("Insufficient plotting space, transformations have been applied.")
      }
    )
  }

  # Return transformed data
  return(x)
}

#' @rdname cyto_transform
#' @export
cyto_transform.transformerList <- function(x,
                                           trans = NULL,
                                           inverse = FALSE,
                                           plot = TRUE,
                                           popup = FALSE,
                                           axes_limits = "machine",
                                           ...) {

  # Apply transformations to flowFrame/flowSet
  if (is(x, "flowFrame") |
    is(x, "flowSet")) {

    # Extract transformations to transformList
    transform_list <- cyto_transform_extract(trans, inverse = inverse)

    # Apply transformations
    x <- suppressMessages(transform(x, transform_list))


    # Apply transformations to GatingHierarchy/GatingSet
  } else if (is(x, "GatingHierarchy") |
    is(x, "GatingSet")) {

    # Inverse transformations not supported
    if (inverse == TRUE) {
      stop(paste(
        "Inverse transformations are not yet supported for",
        "GatingHierarchy/GatingSet objects."
      ))
    }

    # Apply transformations
    x <- suppressMessages(transform(x, trans))
  }

  # Construct plots
  if (plot == TRUE) {
    # Plot if sufficient space
    tryCatch(
      {

        # Extract flowFrame/flowSet for plotting
        cyto_data <- cyto_extract(x)

        # Convert to flowFrame for plotting
        cyto_data <- cyto_convert(cyto_data, "flowFrame")

        # Channels
        channels <- names(trans)

        # Old graphics parameters
        old_pars <- .par("mfrow")
        on.exit(par(old_pars))

        # Set up plotting area
        cyto_plot_new(popup = popup)
        n <- length(channels)
        cyto_plot_layout(
          n2mfrow(n)[1],
          n2mfrow(n)[2]
        )

        # Generate plot for each channel
        lapply(channels, function(chan) {
          if (inverse == FALSE) {
            cyto_plot(cyto_data,
              channels = chan,
              axes_trans = trans,
              title = NA,
              axes_limits = axes_limits
            )
          } else if (inverse == TRUE) {
            cyto_plot(cyto_data,
              channels = chan,
              title = NA,
              axes_limits = axes_limits
            )
          }
        })
      },
      error = function(e) {
        message("Insufficient plotting space, transformations have been applied.")
      }
    )
  }

  # Return transformed data
  return(x)
}

## CYTO_TRANSFORM_EXTRACT ------------------------------------------------------

#' Extract Transformations from TransformerList
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}}.
#' @param inverse logical indicating whether the returned
#'   \code{\link[flowCore:transformList-class]{transformList}} should contain
#'   the inverse transformations.
#'
#' @return A \code{\link[flowCore:transformList-class]{transformList}}
#'   containing the desired transformations.
#'
#' @importFrom flowCore transformList
#' @importFrom methods is
#'
#' @examples
#'
#' # Load CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Load in samples to flowSet
#' fs <- Activation
#'
#' # Add fs to GatingSet
#' gs <- GatingSet(fs)
#'
#' # Convert transformerList into transformList
#' trans <- estimateLogicle(gs[[32]], cyto_fluor_channels(gs))
#' trans_list <- cyto_transform_extract(trans)
#'
#' # Convert transformerList into inverse transformList
#' inv <- cyto_transform_extract(trans, inverse = TRUE)
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_transform_extract <- function(x,
                                   inverse = FALSE) {

  # TransformLists are returned unaltered
  if (is(x, "transformList")) {
    return(x)
    # TransformList extracted from transformerList
  } else if (is(x, "transformerList")) {
    # Extract transformations into transformList
    if (inverse == TRUE) {
      x <- transformList(names(x), lapply(x, `[[`, "inverse"))
    } else {
      x <- transformList(names(x), lapply(x, `[[`, "transform"))
    }
    # Return transformList
    return(x)
  }
}

## CYTO_EXTRACT ----------------------------------------------------------------

#' Extract a valid flowFrame or flowSet
#'
#' \code{cyto_extract} is essentially a wrapper for
#' \code{\link[flowWorkspace:gh_pop_get_data]{gs_pop_get_data}} which also
#' accepts \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#' \code{\link[flowCore:flowSet-class]{flowSet}} objects. The \code{parent}
#' population is extracted from
#' \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} objects whilst
#' \code{flowFrame} or \code{flowSet} objects are returned as is.
#'
#' @param x object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param parent name of the parent population to extract from
#'   \code{GatingHierarchy} or \code{GatingSet} objects.
#' @param copy logical indicating whether a deep copy of the extracted data
#'   should be returned.
#' @param raw logical indicating whether a list of raw data matrices should be
#'   returned instead of a flowFrame or flowSet.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:gh_pop_get_data]{gh_pop_get_data}} or
#'   \code{\link[flowWorkspace:gh_pop_get_data]{gs_pop_get_data}}.
#'
#' @return either a \code{flowFrame} or a \code{cytoset}  by default. A list of
#'   raw data matrices when raw is set to TRUE.
#'
#' @importFrom flowWorkspace gs_pop_get_data gh_pop_get_data realize_view
#' @importFrom flowCore exprs
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Extract flowFrame
#' cyto_extract(gs[[1]], parent = "root")
#'
#' # Extract cytoset
#' cyto_extract(gs, parent = "root")
#' 
#' # Extract raw data matrices
#' cyto_extract(gs, parent = "root", raw = TRUE)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_extract <- function(x, 
                         parent = NULL, 
                         copy = FALSE,
                         raw = FALSE, ...) {

  # DEFAULT PARENT
  if(is.null(parent)){
    parent <- cyto_nodes(x, path = "auto")[1]
  }
  
  # EXTRACT
  if (is(x, "GatingHierarchy")) {
    x <- gh_pop_get_data(x, parent, ...)
  } else if (is(x, "GatingSet")) {
    x <- gs_pop_get_data(x, parent, ...)
  }
  
  # COPY
  if(copy){
    x <- cyto_copy(x)
  }
  
  # RAW DATA MATRICES
  if(raw){
    nms <- cyto_names(x)
    if(is(x, "flowFrame")){
      x <- list(exprs(x))
    }else if(is(x, "flowSet")){
      y <- lapply(seq_along(x), function(z){
        exprs(x[[z]])
      })
      x <- y
    }
    names(x) <- nms
  }
  
  # RETURN EXTRACTED DATA
  return(x)
}

## CYTO_CONVERT ----------------------------------------------------------------

#' Convert between cytometry objects
#'
#' @param x \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param return either 'flowFrame', 'flowSet', 'GatingHierarchy', 'GatingSet',
#'   coerced 'flowFrame list' or coerced 'flowSet list'. GatingSet and flowSet
#'   objects can also be converted to a 'list of flowFrames'.
#' @param parent name of parent population to extract from
#'   \code{GatingHierarchy} and \code{GatingSet} objects.
#' @param ... not in use.
#'
#' @return object specified by 'return' argument.
#'
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace GatingSet sampleNames
#' @importFrom methods as
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Convert flowSet to 'list of flowFrames'
#' cyto_convert(Activation, "list of flowFrames")
#'
#' # Convert flowSet to 'flowFrame'
#' cyto_convert(Activation, "flowFrame")
#'
#' # Convert GatingSet to flowFrame
#' cyto_convert(GatingSet(Activation), "flowFrame", parent = "root")
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
                                   return = "flowFrame",
                                   ...) {

  # NAME
  nm <- cyto_names(x)

  # CONVERSIONS
  if (return == "list of flowFrames") {
    return <- "flowFrame list"
  }

  if (return == "flowFrame") {

  } else if (return %in% c("flowFrame list", "list of flowFrames")) {
    x <- list(x)
  } else if (return == "flowSet") {
    x <- flowSet(x)
    sampleNames(x) <- nm
    x <- as(x, "ncdfFlowSet")
  } else if (return == "flowSet list") {
    x <- list(flowSet(x))
    sampleNames(x[[1]]) <- nm
    x <- list(as(x[[1]], "ncdfFlowSet"))
  } else if (return == "GatingSet") {
    x <- flowSet(x)
    sampleNames(x) <- nm
    x <- as(x, "ncdfFlowSet")
    x <- GatingSet(x)
  } else if (return == "GatingHierarchy") {
    x <- flowSet(x)
    sampleNames(x) <- nm
    x <- as(x, "ncdfFlowSet")
    x <- GatingSet(x)[[1]]
  }

  return(x)
}

#' @rdname cyto_convert
#' @export
cyto_convert.flowSet <- function(x,
                                 return = "flowSet",
                                 ...) {
  if (return == "flowSet") {

  } else if (return == "flowFrame") {
    x <- as(x, "flowFrame")
    # REMOVE ORIGINAL PARAMETER
    if ("Original" %in% cyto_channels(x)) {
      # CANNOT BE EMPTY FLOWFRAME
      if (nrow(x) == 0) {
        # ADD EVENT
        x@exprs <- rbind(rep(0, length(colnames(x))), x@exprs)
        # REMOVE ORIGINAL PARAMETER & ADDED EVENT
        x <- suppressWarnings(
          x[-1, -match("Original", cyto_channels(x))]
        )
      } else {
        x <- suppressWarnings(
          x[, -match("Original", cyto_channels(x))]
        )
      }
    }
  } else if (return == "flowFrame list") {
    x <- as(x, "flowFrame")
    # REMOVE ORIGINAL PARAMETER
    if ("Original" %in% cyto_channels(x)) {
      # CANNOT BE EMPTY FLOWFRAME
      if (nrow(x) == 0) {
        # ADD EVENT
        x@exprs <- rbind(rep(0, length(colnames(x))), x@exprs)
        # REMOVE ORIGINAL PARAMETER & ADDED EVENT
        x <- suppressWarnings(
          x[-1, -match("Original", cyto_channels(x))]
        )
      } else {
        x <- suppressWarnings(
          x[, -match("Original", cyto_channels(x))]
        )
      }
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
    # REMOVE ORIGINAL PARAMETER
    if ("Original" %in% cyto_channels(x)) {
      # CANNOT BE EMPTY FLOWFRAME
      if (nrow(x) == 0) {
        # ADD EVENT
        x@exprs <- rbind(rep(0, length(colnames(x))), x@exprs)
        # REMOVE ORIGINAL PARAMETER & ADDED EVENT
        x <- suppressWarnings(
          x[-1, -match("Original", cyto_channels(x))]
        )
      } else {
        x <- suppressWarnings(
          x[, -match("Original", cyto_channels(x))]
        )
      }
    }
    x <- as(flowSet(x), "ncdfFlowSet")
    x <- GatingSet(x)[[1]]
  }

  return(x)
}

#' @rdname cyto_convert
#' @export
cyto_convert.GatingHierarchy <- function(x,
                                         parent = "root",
                                         return = "GatingHierarchy",
                                         ...) {

  # NAME
  nm <- cyto_names(x)

  if (return == "GatingHierarchy") {

  } else if (return == "flowFrame") {
    x <- cyto_extract(x, parent)
  } else if (return %in% c("flowFrame list", "list of flowFrames")) {
    x <- list(cyto_extract(x, parent))
  } else if (return == "flowSet") {
    x <- flowSet(cyto_extract(x, parent))
    sampleNames(x) <- nm
    x <- as(x, "ncdfFlowSet")
  } else if (return == "flowSet list") {
    x <- list(flowSet(cyto_extract(x, parent)))
    sampleNames(x[[1]]) <- nm
    x <- list(as(x[[1]], "ncdfFlowSet"))
  }

  return(x)
}

#' @rdname cyto_convert
#' @export
cyto_convert.GatingSet <- function(x,
                                   parent = "root",
                                   return = "GatingSet",
                                   ...) {
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
  } else if (return == "GatingHierarchy") {
    x <- as(
      flowSet(cyto_convert(cyto_extract(x, parent), "flowFrame")),
      "ncdfFlowSet"
    )
    x <- GatingSet(x)[[1]]
  }

  return(x)
}

## CYTO_FILTER -----------------------------------------------------------------

#' Filter samples based on experiment variables
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... tidyverse-style subsetting using comma separated logical
#'   predicates based on experimental variables stored in
#'   \code{cyto_details(x)}. See examples below for demonstration.
#'
#' @return \code{flowSet} or \code{GatingSet} restricted to samples which meet
#'   the filtering criteria.
#'
#' @importFrom dplyr filter
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
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
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_filter <- function(x, ...) {

  # Check class of x
  if (!any(is(x, "flowSet") |
    is(x, "GatingSet"))) {
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

## CYTO_SELECT -----------------------------------------------------------------

# Similar to cyto_filter but acts in a non-tidyverse way.

#' Select samples based on experiment variables
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... named list containing experimental variables to be used to select
#'   samples or named arguments containing the levels of the variables to
#'   select. See below examples for use cases. Selected samples can be excluded
#'   by setting \code{exclude} to TRUE.
#'
#' @return \code{flowSet} or \code{GatingSet} restricted to samples which meet
#'   the designated selection criteria.
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
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
#'
#' # Exclude Stim-D treatment group
#' fs <- cyto_select(Activation,
#'   Treatment = "Stim-D",
#'   exclude = TRUE
#' )
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_select <- function(x,
                        ...) {

  # Check class of x
  if (!any(is(x, "flowSet") |
    is(x, "GatingSet"))) {
    stop("'x' should be an object of class flowSet or GatingSet.")
  }

  # Pull down ... arguments to list
  args <- list(...)

  # ... is already a named list of arguments
  if (class(args[[1]]) == "list") {
    args <- args[[1]]
  }

  # Exclude
  if (any(grepl("exclude", names(args)))) {
    exclude <- args[[which(grepl("exclude", names(args)))]]
    args <- args[-which(grepl("exclude", names(args)))]
  } else {
    exclude <- FALSE
  }

  # Extract experiment details
  pd <- cyto_details(x)

  # Check that all variables are valid
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

  # Exclude
  if (exclude == TRUE) {
    return(x[-ind])
  } else {
    return(x[ind])
  }
  
}

## CYTO_GROUP_BY ---------------------------------------------------------------

#' Group a flowSet or GatingSet by experiment variables
#'
#' @param x an object of class \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param group_by names of cyto_details variables to use for merging. Set to
#'   "all" to merge all samples in \code{x}.
#'
#' @return a named list of \code{flowSet} or \code{GatingSet} objects
#'   respectively.
#'
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
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
  if (!any(is(x, "flowSet") | is(x, "GatingSet"))) {
    stop("'x' should be an object of class flowSet or GatingSet.")
  }

  # Extract experiment details
  pd <- cyto_details(x)

  # Replace any NA with "NA" to avoid missing rows
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

## CYTO_MERGE_BY ---------------------------------------------------------------

#' Merge a flowSet by experiment variables
#'
#' \code{cyto_merge_by} makes a call to \code{cyto_group_by} to split samples
#' into groups based on experiment variables. The resulting groups are then
#' converted to flowFrames using \code{cyto_convert}. \code{cyto_merge_by} is
#' the preferred way to merge samples in CytoExploreR as it will ensure
#' appropriate sampling in \code{cyto_plot}.
#'
#' @param x object of class \code{flowSet}.
#' @param parent name of the parent population to merge when a \code{GatingSet}
#'   object is supplied, set to the \code{"root"} node by default.
#' @param merge_by vector of \code{\link{cyto_details}} column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to merging.
#' @param select selection critieria passed to \code{\link{cyto_select}} which
#'   indicates which samples in each group to retain prior to merging, set to
#'   NULL by default to merge all samples in each group. Filtering steps should
#'   be comma separated and wrapped in a list. Refer to
#'   \code{\link{cyto_select}} for more details.
#' @param barcode logical indicating whether a call should be made to
#'   \code{\link{cyto_barcode}} prior to grouping and merging samples, set to
#'   TRUE by default. Barcoding helps \code{\link{cyto_sample}} to appropriately
#'   sample events based on the number of merged samples.
#' @param ... additional arguments passed to \code{\link{cyto_barcode}}.
#'
#' @return list of flowFrames merged by the grouping variables specified by
#'   \code{merge_by}.
#'
#' @importFrom flowCore `identifier<-`
#' @importFrom methods is
#'
#' @examples
#'
#' # Load CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#'
#' # Experiment details
#' cyto_details(fs)
#'
#' # Merge samples by 'Treatment'
#' fr_list <- cyto_merge_by(fs, "Treatment")
#'
#' # Merge samples by 'OVAConc'
#' fr_list <- cyto_merge_by(fs, "OVAConc")
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_group_by}}
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_barcode}}
#' @seealso \code{\link{cyto_sample}}
#'
#' @name cyto_merge_by
NULL

#' @noRd
#' @export
cyto_merge_by <- function(x, ...) {
  UseMethod("cyto_merge_by")
}

#' @rdname cyto_merge_by
#' @export
cyto_merge_by.GatingSet <- function(x,
                                    parent = "root",
                                    merge_by = "all",
                                    select = NULL,
                                    barcode = TRUE,
                                    ...) {

  # EXTRACT POPULATON
  fs <- cyto_extract(x, parent)

  # CALL FLOWSET METHOD
  cyto_merge_by(fs,
    merge_by = merge_by,
    select = select,
    barcode = barcode,
    ...
  )
}

#' @rdname cyto_merge_by
#' @export
cyto_merge_by.flowSet <- function(x,
                                  merge_by = "all",
                                  select = NULL,
                                  barcode = TRUE,
                                  ...) {

  # BARCODING ------------------------------------------------------------------

  # SAMPLE ID
  x <- cyto_barcode(x, ...)

  # GROUPING -------------------------------------------------------------------

  # CYTO_GROUP_BY
  fs_list <- cyto_group_by(x, group_by = merge_by)

  # GROUPS
  grps <- names(fs_list)
  
  # COMBINED EVENTS
  if("all" %in% grps){
    grps[which("all" %in% grps)] <- "Combined Events"
  }
  
  # SELECTION ------------------------------------------------------------------

  # ATTEMPT SELECTION OR RETURN ALL SAMPLES
  if (!is.null(select)) {
    fs_list <- lapply(fs_list, function(z) {
      # Select or return all samples if criteria not met
      tryCatch(cyto_select(z, select), error = function(e) {
        z
      })
    })
  }

  # MERGING --------------------------------------------------------------------

  # CONVERT EACH GROUP TO FLOWFRAME
  fr_list <- lapply(fs_list, function(fs) {
    if (!is(fs, "flowFrame")) {
      cyto_convert(fs, "flowFrame")
    } else {
      fs
    }
  })
  names(fr_list) <- grps

  # REPLACE SAMPLENAMES WITH GROUPS
  if (!all(cyto_names(fr_list) %in% grps)) {
    lapply(seq_len(length(fr_list)), function(z) {
      identifier(fr_list[[z]]) <<- grps[z]
    })
  }

  # RETURN PREPARED FLOWFRAME LIST
  return(fr_list)
}

## CYTO_SPLIT ------------------------------------------------------------------

#' Split samples merged with cyto_merge
#'
#' Extract individual samples merged using \code{cyto_merge()} based on
#' \code{"Sample ID"} column created by \code{cyto_barcode()}.
#'
#' @param x object of class \code{flowFrame}.
#' @param names vector of names to assign to each of the extracted flowFrames
#'   when saving the split files. Name should be supplied in the order used
#'   prior to merging.
#'
#' @return list of split flowFrames.
#'
#' @importFrom flowCore flowFrame exprs identifier keyword split
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#'
#' # Load CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Merge samples
#' fr <- cyto_merge_by(fs, "all")[[1]]
#'
#' # Split merged samples
#' fr_list <- cyto_split(fr, names = cyto_names(fs))
#' @seealso \code{\link{cyto_merge_by}}
#' @seealso \code{\link{cyto_barcode}}
#'
#' @export
cyto_split <- function(x,
                       names = NULL) {

  # CHECKS ---------------------------------------------------------------------

  # FLOWFRAME
  if (!is(x, "flowFrame")) {
    stop("cyto_split() expects a flowFrame object.")
  }

  # SAMPLE ID
  if (!"Sample ID" %in% cyto_channels(x)) {
    stop("Merged samples must be barcoded in cyto_merge().")
  }

  # SPLIT INTO FLOWFRAMES ------------------------------------------------------

  # EXTRACT DATA
  fr_exprs <- exprs(x)

  # SAMPLE IDs
  sample_id <- unique(fr_exprs[, "Sample ID"])
  samples <- length(sample_id)

  # SPLIT BY SAMPLE ID
  fr_list <- split(x, factor(fr_exprs[, "Sample ID"], levels = sample_id))

  # NAMES
  if (!is.null(names)) {
    # INSUFFICIENT NAMES
    if (length(names) != length(sample_id)) {
      stop("Supply a name for each file.")
    }
    # NO NAMES
  } else {
    names <- paste0("Sample-", sample_id)
  }

  # NAME SPLIT FILES
  lapply(seq_len(samples), function(z) {
    identifier(fr_list[[z]]) <<- names[z]
  })
  names(fr_list) <- names

  # RETURN SPLIT FLOWFRAMES
  return(fr_list)
}

## CYTO_SAVE -------------------------------------------------------------------

#' Write samples to FCS files in new folder or save GatingSet
#'
#' @param x object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param parent name of the parent population to extract when a
#'   \code{GatingHierarchy} or \code{GatingSet} object is supplied. If the name
#'   of the parent is supplied the samples will be written to FCS files in the
#'   specified \code{save_as} directory. Otherwise the entire \code{GatingSet}
#'   or \code{GatingHierarchy} will be saved to the specified \code{save_as}
#'   directory.
#' @param split logical indicating whether samples merged using
#'   \code{cyto_merge_by} should be split prior to writing FCS files, set to
#'   FALSE by default.
#' @param names original names of the samples prior to merging using
#'   \code{cyto_merge_by}, only required when split is TRUE. These names will be
#'   re-assigned to each of split flowFrames.
#' @param save_as name of the folder to which the written FCS files should be
#'   saved, set to NULL by default to save the files to the current working
#'   directory. To prevent files being overwritten, it is recommended that
#'   \code{save_as} directory not be manually created before running
#'   \code{cyto_save}.
#' @param inverse_transform logical indicating whether the data should be
#'   inverse transformed prior to writing FCS files, set to FALSE by default.
#'   Inverse transformations of \code{flowFrame} or \code{flowSet} objects
#'   requires passing of transformers through the \code{trans} argument.
#' @param trans object of class \code{transformerList} containing the
#'   transformation definitions applied to the supplied data. Used internally
#'   when \code{inverse_transform} is TRUE, to inverse the transformations prior
#'   to writing FCS files.
#' @param ... not in use.
#'
#' @return list of flowFrames containing the data that was saved to the FCS
#'   files.
#'
#' @importFrom flowCore write.FCS exprs
#' @importFrom flowWorkspace save_gs
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#'
#' \dontrun{
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Save each flowFrame to file
#' cyto_save(fs)
#' 
#' }
#'
#' @seealso \code{\link{cyto_split}}
#'
#' @rdname cyto_save
#'
#' @export
cyto_save <- function(x, ...) {
  UseMethod("cyto_save")
}

#' @rdname cyto_save
#' @export
cyto_save.GatingSet <- function(x,
                                parent = NULL,
                                split = FALSE,
                                names = NULL,
                                save_as = NULL,
                                inverse_transform = FALSE,
                                trans = NULL, ...) {

  # SAVE GATINGSET
  if (is.null(parent)) {
    # SAVE GATINGSET
    suppressMessages(save_gs(x, save_as))
    # RETURN GATINGSET
    invisible(x)
    # SAVE FCS FILES
  } else {
    # EXTRACT DATA
    message(paste("Extracting the ", parent, " node from the GatingSet."))
    fs <- cyto_extract(x, 
                       parent = parent,
                       copy = TRUE)
    # TRANSFORMATIONS
    trans <- cyto_transformer_extract(x)
    # FLOWSET METHOD
    fr_list <- cyto_save(
      x = fs,
      split = split,
      names = names,
      save_as = save_as,
      inverse_transform = inverse_transform,
      trans = trans
    )

    # RETURN DATA
    invisible(fr_list)
  }
}

#' @rdname cyto_save
#' @export
cyto_save.GatingHierarchy <- function(x,
                                      parent = NULL,
                                      split = FALSE,
                                      names = NULL,
                                      save_as = NULL,
                                      inverse_transform = FALSE,
                                      trans = NULL,
                                      ...) {

  # SAVE GATINGHIERARCHY
  if (is.null(parent)) {
    # SAVE GATINGHIERARCHY
    suppressMessages(save_gs(x, save_as))
    # RETURN GATINGHIERARCHY
    invisible(x)
    # SAVE FCS FILES
  } else {
    # EXTRACT DATA
    message(paste("Extracting the ", parent, " node from the GatingHierarchy."))
    fr <- cyto_extract(x, 
                       parent = parent,
                       copy = TRUE)
    # TRANSFORMATIONS
    trans <- cyto_transformer_extract(x)
    # FLOWSET METHOD
    fr_list <- cyto_save(
      x = fr,
      split = split,
      names = names,
      save_as = save_as,
      inverse_transform = inverse_transform,
      trans = trans
    )

    # RETURN DATA
    invisible(fr_list)
  }
}

#' @rdname cyto_save
#' @export
cyto_save.flowSet <- function(x,
                              split = FALSE,
                              names = NULL,
                              save_as = NULL,
                              inverse_transform = FALSE,
                              trans = NULL,
                              ...) {

  # COPY
  x <- cyto_copy(x)
  
  # LIST OF FLOWFRAMES
  fr_list <- cyto_convert(x, "list of flowFrames")
  
  # LIST OF SPLIT FLOWFRAMES
  if(split == TRUE){
    # NAMES SUPPLIED - CHECK LENGTH
    if(!is.null(names)){
      # SAMPLES PER FLOWFRAME
      samples_per_file <- lapply(fr_list, function(z){
        length(unique(exprs(z)[, "Sample ID"]))
      })
      # SPLIT NAMES SUPPLIED
      samples <- sum(unlist(samples_per_file))
      if(length(names) != samples){
        stop(paste("Expecting", samples, "names for the split files."))
      }
      # PREPARE NAMES
      ind <- LAPPLY(seq_along(fr_list), function(z){
        rep(z, samples_per_file[[z]])
      })
      names <- split(names, ind)
      # SPLIT FR_LIST
      fr_list <- mapply(function(z, name){
        cyto_split(z, names = name)
      }, fr_list, names)
      fr_list <- unlist(fr_list)
    # NO NAMES SUPPLIED  
    }else{
      fr_list <- LAPPLY(fr_list, function(z){
        cyto_split(z)
      })
    }
  }
  
  # DIRECTORY CHECK
  if (!is.null(save_as) & dir.exists(save_as)) {
    # FILES WILL BE OVERWRITTEN
    if (any(list.files(save_as) %in% cyto_names(fr_list))) {
      message(paste("Files will be overwritten in", save_as, "."))
      opt <- readline("Do you want to continue? (Y/N)")
      if (grepl("n", opt, ignore.case = TRUE)) {
        invisible(NULL)
      }
    }
  }

  # MESSAGE
  if (is.null(save_as)) {
    location <- "current working directory."
  } else {
    location <- save_as
  }
  message(paste("Writing FCS files to", location, "..."))

  # WRITE FCS FILES
  fr_list <- lapply(fr_list, function(z) {
    # INVERSE TRANSFORM
    if (inverse_transform == TRUE){
      # TRANSFORMERS REQUIRED
      if(is.null(trans)){
        stop("Supply transformerList to 'trans' to inverse transformations.")
      }
      # INVERSE TRANSFORM
      z <- cyto_transform(z,
        trans = trans,
        inverse = TRUE,
        plot = FALSE
      )
    }
    # Message
    message(paste0(cyto_names(z), "..."))
    # NO DIRECTORY SPECIFIED
    if (is.null(save_as)) {
      write.FCS(
        z,
        cyto_names(z)
      )
      # DIRECTORY SPECIFIED
    } else {
      # CREATE DIRECTORY
      if (!dir.exists(save_as)) {
        dir.create(save_as)
      }
      write.FCS(
        z,
        paste0(save_as, "/", cyto_names(z))
      )
    }
    return(z)
  })

  # RETURN DATA
  invisible(unlist(fr_list))
}

#' @rdname cyto_save
#' @export
cyto_save.flowFrame <- function(x,
                                split = FALSE,
                                names = NULL,
                                save_as = NULL,
                                inverse_transform = FALSE,
                                trans = NULL,
                                ...) {

  # COPY
  x <- cyto_copy(x)
  
  # SPLIT
  if (split == TRUE) {
    fr_list <- cyto_split(x,
      names = names
    )
  } else {
    fr_list <- list(x)
  }

  # DIRECTORY CHECK
  if (!is.null(save_as) & dir.exists(save_as)) {
    # FILES WILL BE OVERWRITTEN
    if (any(list.files(save_as) %in% cyto_names(fr_list))) {
      message(paste("Files will be overwritten in", save_as, "."))
      opt <- readline("Do you want to continue? (Y/N)")
      if (grepl("n", opt, ignore.case = TRUE)) {
        invisible(NULL)
      }
    }
  }

  # MESSAGE
  if (is.null(save_as)) {
    location <- "current working directory."
  } else {
    location <- save_as
  }
  message(paste("Writing FCS files to", location, "..."))

  # WRITE FCS FILES
  fr_list <- lapply(fr_list, function(z) {
    # INVERSE TRANSFORM
    if (inverse_transform == TRUE){
      #TRANSFORMERS REQUIRED
      if(is.null(trans)){
        stop("Supply transformerList to 'trans' to inverse transformations.")
      }
      # INVERSE TRANSFORM
      z <- cyto_transform(z,
        trans = trans,
        inverse = TRUE,
        plot = FALSE
      )
    }
    # Message
    message(paste0(cyto_names(z), "..."))
    # NO DIRECTORY SPECIFIED
    if (is.null(save_as)) {
      write.FCS(
        z,
        cyto_names(z)
      )
      # DIRECTORY SPECIFIED
    } else {
      # CREATE DIRECTORY
      if (!dir.exists(save_as)) {
        dir.create(save_as)
      }
      write.FCS(
        z,
        paste0(save_as, "/", cyto_names(z))
      )
    }
    return(z)
  })

  # RETURN DATA
  invisible(unlist(fr_list))
}

## CYTO_SAMPLE -----------------------------------------------------------------

#' Sample a flowFrame or flowSet
#'
#' \code{cyto_sample} allows restriction of a flowFrame or flowSet by indicating
#' the percentage or number of events to retain.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#'   \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param display can be either a numeric [0,1] or integer to indicate the
#'   percentage or number of events to keep respectively.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
#' @param plot logical required for lists to indicate whether sampling should be
#'   scaled per flowFrame to retain original ratios, as used in cyto_plot.
#' @param ... not in use.
#'
#' @return \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#'   \code{\link[flowCore:flowSet-class]{flowSet}} restricted to \code{display}
#'   percentage events.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter Subset fsApply exprs
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access files
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Restrict first sample by 50%
#' cyto_sample(fs[[1]], 0.5)
#'
#' # Restrict first sample to 10000 events
#' cyto_sample(fs[[1]], 10000)
#' 
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
                                  seed = NULL,
                                  ...) {

  # NO SAMPLING - EMPTY FLOWFRAME
  if (nrow(x@exprs) == 0) {
    return(x)
  }

  # Do nothing if no sampling required
  if (display != 1) {

    # Number of events
    events <- nrow(x)

    # display is the number of events to keep
    if (display > 1) {

      # display is too large - retain all events
      if (display > events) {
        return(x)
        # display is sample of x
      } else {
        size <- display
      }

      # display is a proportion of events to keep
    } else {

      # Size
      size <- display * events
    }

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
                                seed = NULL,
                                ...) {
  # Apply same degree of sampling to each flowFrame in the flowSet
  fsApply(x, cyto_sample, display = display, seed = seed)
}

#' @rdname cyto_sample
#' @export
cyto_sample.list <- function(x,
                             display = 1,
                             seed = NULL,
                             plot = FALSE,
                             ...) {

  # CYTO_PLOT SAMPLING - RATIOS
  if (plot == TRUE) {

    # LIST OF FLOWSETS
    if (all(LAPPLY(x, function(z) {
      is(z, "flowSet")
    }))) {
      # Same sampling applied to all samples
      x <- lapply(x, function(z) {
        cyto_sample(z, display = display, seed = seed)
      })
      # LIST OF FLOWFRAMES
    } else if (all(LAPPLY(x, function(z) {
      is(z, "flowFrame")
    }))) {
      # BARCODED FLOWFRAMES
      if (any(LAPPLY(x, function(z) {
        "Sample ID" %in% cyto_channels(z)
      }))) {
        # SAMPLES PER FLOWFRAME - ASSUME FLOWFRAME IF NO BARCODE
        samples <- LAPPLY(x, function(z) {
          tryCatch(length(unique(exprs(z)[, "Sample ID"])),
            error = function(e) {
              1
            }
          )
        })
        # EVENTS PER SAMPLE - EACH LAYER
        events_per_sample <- LAPPLY(seq_len(length(x)), function(z) {
          nrow(x[[z]]) / samples[z]
        })
        # EVENTS RATIO - BASE LAYER REFERENCE
        events_ratio <- events_per_sample / events_per_sample[1]
        # SCALE SAMPLING EVENTS USING BASE LAYER AS REFERENCE
        if (display <= 1) {
          events_to_sample <- rep(display * nrow(x[[1]]), length(x)) *
            events_ratio
        } else {
          events_to_sample <- rep(display, length(x)) *
            events_ratio
        }
        # SAMPLING
        lapply(seq_len(length(x)), function(z) {
          x[[z]] <<- cyto_sample(x[[z]],
            display = events_to_sample[z],
            seed = seed
          )
        })
        # NO BARCODING - RELY ON IDENTIFIERS
      } else {
        # Same percentage sampling applied to each flowFrame
        if (display <= 1) {
          # Same sampling applied to all samples
          x <- lapply(x, function(z) {
            cyto_sample(z, display = display, seed = seed)
          })
          # Sampling by event number is more complex
        } else if (display > 1) {
          # Identifiers
          nms <- LAPPLY(x, function(z) {
            cyto_names(z)
          })
          ind <- seq_len(length(nms))
          # Sampling
          x <- lapply(ind, function(z) {
            # Base layer sampled as per usual
            if (z == 1) {
              cyto_sample(x[[z]], display = display, seed = seed)
            } else {
              # Identifier matches base - sample size decreased
              if (nms[z] == nms[1]) {
                # Number of events in base layer
                base_events <- nrow(x[[1]])
                # Number of events in overlay
                overlay_events <- nrow(x[[z]])
                # Proportion of overlay relative to base
                prop <- overlay_events / base_events
                # Update display prop * display
                display <- ceiling(prop * display)
                # Sampling
                cyto_sample(x[[z]], display = display, seed = seed)
                # Identifiers don't match - separate samples - same sampling
              } else {
                cyto_sample(x[[z]], display = display, seed = seed)
              }
            }
          })
        }
      }
    }

    # GENERIC SAMPLING
  } else {

    # ALLOW DIFFERENT SAMPLING PER ELEMENT
    display <- rep(display, length(x))
    x <- mapply(function(x, display) {
      cyto_sample(x, display)
    }, x, display)
    
  }

  # Return sampled list
  return(x)
}

## CYTO_BEADS_SAMPLE -----------------------------------------------------------

#' Sample GatingSet to bead count
#'
#' \code{cyto_beads_sample} can be used on samples containing beads, to
#' downsamples each sample in a GatingSet based on a specific \code{bead_count}.
#' For example, if Sample A contains 100 beads and 75000 events, and Sample B
#' contains 50 beads and 5000 events. \code{cyto_beads_sample} will downsample
#' each sample to contain 50 beads and therefore 5000 events for Sample A and
#' 3750 events for Sample B.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param beads name of the gated bead population to use for the calculation. If
#'   not supplied internal checks will be made for populations named "Single
#'   Beads" or "Beads".
#' @param bead_count minimal bead count to down sample to, set to the minimum
#'   bead count in samples by default. The bead count must be less than or equal
#'   to the minimum bead count among the samples.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_beads_sample <- function(x,
                              beads = NULL,
                              bead_count = NULL){
  
  # GATINGSET
  if(!is(x, "GatingSet")){
    stop("'x' must be a Gatingset object.")
  }
  
  # CLONE GATINGSET
  gs_clone <- cyto_copy(x)
  
  # NODES
  nodes <- cyto_nodes(gs_clone, path = "auto")
  
  # BEADS
  if(is.null(beads)){
    # SINGLE BEADS
    if(any(grepl("single beads", nodes, ignore.case = TRUE))){
      beads <- nodes[which(grepl("single beads", nodes, ignore.case = TRUE))[1]]
    # BEADS
    }else if(any(grepl("beads", nodes, ignore.case = TRUE))){
      beads <- nodes[which(grepl("beads", nodes, ignore.case = TRUE))[1]]
    # BEADS MISSING
    }else{
      stop("Supply the name of the 'beads' population.")
    }
  }
  
  # BEAD COUNTS
  bead_pops <- cyto_extract(gs_clone, beads)
  bead_counts <- suppressMessages(
    cyto_stats_compute(bead_pops,
                       stat = "count",
                       format = "wide")
    )
  
  # BEAD COUNT
  if(is.null(bead_count)){
    bead_count <- min(bead_counts[, ncol(bead_counts)])
  }else{
    if(!bead_count <= min(bead_counts[, ncol(bead_counts)])){
      bead_count <- min(bead_counts[, ncol(bead_counts)])
    }
  }
  
  # BEAD RATIOS
  bead_ratios <- 1/(bead_counts[, ncol(bead_counts)]/bead_count)
  
  # BEADS MISSING
  if(any(is.infinite(bead_ratios))){
    stop(paste0("The following samples do not contain any beads:",
                paste0("\n", cyto_names(gs_clone)[is.infinite(bead_ratios)])))
  }
  
  # SAMPLING - ROOT POPULATION
  pops <- list()
  lapply(seq_along(bead_pops), function(z){
    pops[[z]] <<- cyto_sample(cyto_extract(gs_clone[[z]],
                                           copy = TRUE),
                              bead_ratios[[z]])
  })
  names(pops) <- cyto_names(pops)
  pops <- flowSet(pops)
  
  # REPLACE DATA IN GATINGSET
  gs_cyto_data(gs_clone) <- pops
  
  # RETURN SAMPLED GATINGSET
  return(gs_clone)
  
}

## CYTO_BARCODE ----------------------------------------------------------------

#' Barcode each file in a flowSet with a sample ID
#'
#' Adds a new parameter to each of the flowFrames in the flowSet called
#' \code{"Sample ID"} to barcode events from each flowFrame.
#'
#' @param x object of class \code{flowSet} to be barcoded.
#' @param type indicates whether the \code{"samples"}, \code{"events"} or
#'   \code{"both"} should be barcoded, set \code{"samples"} by default.
#'
#' @return barcoded flowSet with \code{"Sample ID"} and/or \code{"Event ID"}
#'   column added and annotated.
#'
#' @importFrom methods is as
#' @importFrom flowWorkspace `sampleNames<-`
#' @importFrom flowCore fsApply
#'
#' @examples
#'
#' # Load in CytoExploreRData to access files
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Barcode
#' cyto_barcode(fs)
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_barcode <- function(x,
                         type = "samples") {

  # CHECKS ---------------------------------------------------------------------

  # FLOWSET
  if (!is(x, "flowSet")) {
    stop("'cyto_barcode' expects objects of class flowSet.")
  }

  # TYPE
  if (type == "both") {
    type <- c("samples", "events")
  }

  # PREPARE DATA ---------------------------------------------------------------

  # SAMPLENAMES
  nms <- cyto_names(x)

  # SAMPLE ID COLUMN - ONLY IF NOT PRESENT
  if ("samples" %in% type &
    !"Sample ID" %in% cyto_channels(x)) {
    x <- fsApply(x, function(fr) {
      ind <- match(cyto_names(fr), nms)
      mat <- matrix(rep(ind, .cyto_count(fr)),
        ncol = 1
      )
      colnames(mat) <- "Sample ID"
      suppressWarnings(cbind(fr, mat))
    })
  }

  # EVENT ID COLUMN - ONLY IF NOT PRESENT
  if ("events" %in% type &
    !"Event ID" %in% cyto_channels(x)) {
    total_events <- fsApply(x, "nrow")
    total_events <- split(
      seq_len(sum(total_events)),
      rep(seq_len(length(x)),
        times = total_events
      )
    )
    names(total_events) <- cyto_names(x)
    x <- fsApply(x, function(fr) {
      events <- total_events[[cyto_names(fr)]]
      mat <- matrix(events,
        ncol = 1
      )
      colnames(mat) <- "Event ID"
      suppressWarnings(cbind(fr, mat))
    })
  }

  # RETURN BARCODED FLOWSET
  return(x)
}

## CYTO_MARKERS_EDIT -----------------------------------------------------------

#' Assign marker names to flowFrame or flowSet
#'
#' \code{cyto_markers_edit} opens an editable table containing a list of
#' channels and markers for a \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}},
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}. Users can edit the
#' \code{name} or \code{desc} columns with updated channel names or marker names
#' respectively. These entries will be updated in the \code{x} upon closing the
#' window and saved to a "Experiment-markers.csv" file for future use.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{flowSet},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingSet}} or
#'   \code{GatingSet}.
#' @param file name of csv file containing columns 'Channel' and 'Marker'.
#'
#' @return save inputs to "Experiment-Markers.csv" and returns updated samples.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters markernames markernames<-
#' @importFrom utils edit write.csv read.csv
#' @importFrom tools file_ext
#' @importFrom methods is
#'
#' @examples
#'
#' \dontrun{
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Add marker names to channels - edit table
#' fs <- cyto_markers_edit(fs)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_markers_edit <- function(x, file = NULL) {

  # check class of x
  cyto_check(x)

  # flowFrame
  if (is(x, "flowFrame")) {

    # Extract details of parameters
    pd <- cyto_details(parameters(x))

    # flowSet
  } else if (is(x, "flowSet")) {

    # Extract details of parameters
    pd <- cyto_details(parameters(x[[1]]))

    # GatingHierarchy
  } else if (is(x, "GatingHierarchy")) {
    fr <- cyto_extract(x, "root")
    pd <- cyto_details(parameters(fr))

    # GatingSet
  } else if (is(x, "GatingSet")) {
    fr <- cyto_extract(x, "root")[[1]]
    pd <- cyto_details(parameters(fr))
  }

  # file missing
  if (is.null(file)) {

    # Check if file already exists
    if (length(grep("Experiment-Markers.csv", list.files())) != 0) {
      message("Experiment-Markers.csv found in working directory.")

      # Could be multiple files - check for matching channels
      found_files <- list.files()[grep("Experiment-Markers.csv", list.files())]
      n <- length(grep("Experiment-Markers.csv", list.files()))

      # Run through each file and check channels match samples
      dt <- lapply(seq_len(n), function(z) {
        mrks <- read.csv(found_files[z],
          header = TRUE,
          stringsAsFactors = FALSE
        )
        rownames(mrks) <- NULL
        # SampleNames match those of x
        if (identical(mrks$Channel, colnames(x))) {
          return(mrks)
        } else {
          return(NULL)
        }
      })
      names(dt) <- found_files

      # Files found but don't match
      if (all(LAPPLY(dt, "is.null"))) {

        # Make data.frame with channel and marker columns
        dt <- pd[, c("name", "desc")]
        colnames(dt) <- c("Channel", "Marker")
        rownames(dt) <- NULL
      } else {

        # Remove NULL entries from list - result should be of length 1
        dt[LAPPLY(dt, "is.null")] <- NULL
        file <- names(dt)[1]
        dt <- dt[[1]]
      }
    } else {

      # Make data.frame with channel and marker columns
      dt <- pd[, c("name", "desc")]
      colnames(dt) <- c("Channel", "Marker")
      rownames(dt) <- NULL
    }

    # File manually supplied
  } else {

    # File extension missing
    if (file_ext(file) == "") {
      file <- paste0(file, ".csv")
    }

    # File already exists
    if (length(grep(file, list.files())) != 0) {
      message(file, "found in working directory.")
      dt <- read.csv(file,
        header = TRUE,
        stringsAsFactors = FALSE
      )

      # File does not exist (yet)
    } else {

      # Make data.frame with channel and marker columns
      dt <- dt[, c("name", "desc")]
      colnames(dt) <- c("Channel", "Marker")
      rownames(dt) <- NULL
    }
  }

  # Edit dt - data.frame must have no rownames for editor
  dt <- suppressWarnings(edit(dt))
  
  # Update channels
  BiocGenerics::colnames(x) <- as.character(dt$Channel)
  
  # Write result to csv file -  file name manually supplied
  if (!is.null(file)) {
    write.csv(dt, file, row.names = FALSE)

    # Write result to csv file - no file name supplied
  } else {

    # Save file with date and experiment markers
    write.csv(dt, paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Markers.csv"
    ), row.names = FALSE)
  }

  # Markers and channels
  cyto_channels <- dt$Channel
  cyto_markers <- dt$Marker
  names(cyto_markers) <- cyto_channels

  # Only modify markers if supplied
  if (!all(is.na(cyto_markers))) {
    markernames(x) <- cyto_markers[!is.na(cyto_markers)]
  }

  return(x)
}

## CYTO_DETAILS_EDIT -----------------------------------------------------------

#' Interactively edit cyto_details for a flowSet or GatingSet
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param file name of csv file containing experimental information.
#'
#' @return NULL and return \code{flowSet} or \code{GatingSet} with updated
#'   experimental details.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore pData<-
#' @importFrom utils edit write.csv read.csv
#' @importFrom tools file_ext
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#'
#' # Edit cyto_details in table editor
#' cyto_details_edit(fs)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @export
cyto_details_edit <- function(x, file = NULL) {

  # x should be a flowSet or GatingSet
  if (!any(is(x, "flowSet") | is(x, "GatingSet"))) {
    stop("Please supply either a flowSet or a GatingSet")
  }

  # File missing
  if (is.null(file)) {
    if (length(grep("Experiment-Details.csv", list.files())) != 0) {
      message("Experiment-Details.csv found in working directory.")

      # Could be multiple files - check for matching channels
      found_files <- list.files()[grep("Experiment-Details.csv", list.files())]
      n <- length(grep("Experiment-Details.csv", list.files()))

      # Run through each file and check channels match samples
      pd <- lapply(seq_len(n), function(z) {
        pdata <- read.csv(found_files[z],
          header = TRUE,
          stringsAsFactors = FALSE
        )
        # SampleNames match those of x
        if (all(cyto_names(x) %in% as.vector(pdata$name))) {
          return(pdata)
        } else {
          return(NULL)
        }
      })
      names(pd) <- found_files

      # Files found but don't match
      if (all(LAPPLY(pd, "is.null"))) {

        # Extract cyto_details
        pd <- cyto_details(x)
        rownames(pd) <- NULL
      } else {

        # Remove NULL entries from list - result should be of length 1
        pd[LAPPLY(pd, "is.null")] <- NULL
        file <- names(pd)[1]
        pd <- pd[[1]]
      }
    } else {

      # Extract cyto_details
      pd <- cyto_details(x)
      rownames(pd) <- NULL
    }

    # File name supplied manually
  } else {

    # File name lacks csv extension
    if (file_ext(file) == "") {
      file <- paste0(file, ".csv")
    }

    # File already exists
    if (length(grep(file, list.files())) != 0) {
      message(paste(file, "found in working directory."))
      pd <- read.csv(file, header = TRUE)

      # File does not exist
    } else {

      # Extract cyto_details
      pd <- cyto_details(x)
      rownames(pd) <- NULL
    }
  }

  # Edit cyto_details
  pd <- suppressWarnings(edit(pd))
  rownames(pd) <- pd$name
  
  # Update cyto_details
  cyto_details(x) <- pd

  # Write result to csv file - file name manually supplied
  if (!is.null(file)) {
    write.csv(cyto_details(x), file, row.names = FALSE)

    # Write result to csv file - no file name supplied
  } else {

    # Save file with date and experiment markers
    write.csv(cyto_details(x), paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Details.csv"
    ), row.names = FALSE)
  }

  return(x)
}

## CYTO_COMPENSATE -------------------------------------------------------------

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
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param spillover name of the spillover matrix csv file (e.g.
#'   "Spillover-Matrix.csv") saved in the current working directory.
#' @param select index or name of the sample from which the spillover matrix
#'   should be extracted when no spillover matrix file is supplied to
#'   \code{spillover}. To compensate each sample individually using their stored
#'   spillover matrix file, set \code{select} to NULL.
#' @param ... not in use.
#'
#' @return a compensated \code{flowFrame}, \code{flowSet} or \code{GatingSet}
#'   object.
#'
#' @importFrom utils read.csv
#' @importFrom tools file_ext
#' @importFrom methods is
#'
#' @examples
#'
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
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
#' 
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
cyto_compensate.GatingSet <- function(x,
                                      spillover = NULL,
                                      select = 1,
                                      ...) {

  # Compensation already applied
  if(!is.null(cyto_spillover_extract(x))){
    stop("Looks like the data has already been compensated...")
  }
  
  # Extract flowSet
  fs <- cyto_extract(x, parent = "root")

  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (is(spillover, "character")) {
      # No file extension
      if (file_ext(spillover) == "") {
        spillover <- paste0(spillover, ".csv")
      }
      # Check working directory for csv file
      if (getOption("CytoExploreR_wd_check")) {
        if (!file_wd_check(spillover)) {
          stop(paste(spillover, "does not exist in this working directory."))
        }
      }
      spill <- read.csv(spillover, header = TRUE, row.names = 1)
      colnames(spill) <- rownames(spill)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% BiocGenerics::colnames(fs)) |
        !all(rownames(spill) %in% BiocGenerics::colnames(fs))) {
        stop(
          paste(
            "'spillover' must have valid fluorescent channels as rownames",
            "and colnames."
          )
        )
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(fs))
      names(spill) <- cyto_names(fs)
      # spillover is a matrix/data.frame
    } else if (is(spillover, "matrix") |
      is(spillover, "data.frame")) {
      # column names must be valid channels (rownames not essential)
      if (!all(colnames(spillover) %in% BiocGenerics::colnames(fs))) {
        stop("'spillover' must have valid fluorescent channels as colnames.")
      } else {
        spill <- spillover
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(fs))
      names(spill) <- cyto_names(fs)
    }
    # Extract spillover matrix directly from fs
  } else if (is.null(spillover)) {
    if (!is.null(select)) {
      spill <- cyto_spillover_extract(fs[[select]])
      if(is.null(spill)){
        stop("Unable to extract spillover matrix from selected sample.")
      }
      spill <- rep(spill, length(fs))
      names(spill) <- cyto_names(fs)
    } else {
      
      spill <- lapply(cyto_names(fs), function(y) {
        sm <- cyto_spillover_extract(fs[[y]])[[1]]
        if(is.null(sm)){
          stop(paste0("Unable to extract spillover matrix from ",
                      cyto_names(fs[[y]]), "."))
        }
        return(sm)
      })
    }
  }

  print(spill)
  
  # Channels
  fluor_channels <- cyto_fluor_channels(fs)
  
  # Spillover may contain more channels than in samples
  spill <- lapply(spill, function(z){
    # Select rows - square matrix
    if(nrow(z) == ncol(z)){
      z <- z[match(fluor_channels, colnames(z)), ]
    }
    # Select columns
    z <- z[, fluor_channels]
    return(z)
  })
  
  # Apply compensation
  flowWorkspace::compensate(x, spill)
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowSet <- function(x,
                                    spillover = NULL,
                                    select = 1,
                                    ...) {

  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (is(spillover, "character")) {
      # No file extension
      if (file_ext(spillover) == "") {
        spillover <- paste0(spillover, ".csv")
      }
      # Check working directory for csv file
      if (getOption("CytoExploreR_wd_check")) {
        if (!file_wd_check(spillover)) {
          stop(paste(spillover, "does not exist in this working directory."))
        }
      }
      spill <- read.csv(spillover, header = TRUE, row.names = 1)
      colnames(spill) <- rownames(spill)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% BiocGenerics::colnames(x)) |
        !all(rownames(spill) %in% BiocGenerics::colnames(x))) {
        stop(
          paste(
            "'spillover' must have valid fluorescent channels as rownames",
            "and colnames."
          )
        )
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(x))
      names(spill) <- cyto_names(x)
      # spillover is a matrix/data.frame
    } else if (is(spillover, "matrix") |
      is(spillover, "data.frame")) {
      # column names must be valid channels (rownames not essential)
      if (!all(colnames(spillover) %in% BiocGenerics::colnames(x))) {
        stop("'spillover' must have valid fluorescent channels as colnames.")
      } else {
        spill <- spillover
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(x))
      names(spill) <- cyto_names(x)
    }
    # Extract spillover matrix directly from x
  } else if (is.null(spillover)) {
    if (!is.null(select)) {
      spill <- cyto_spillover_extract(x[[select]])
      if(is.null(spill)){
        stop("Unable to extract spillover matrix from selected sample.")
      }
      spill <- rep(spill, length(x))
      names(spill) <- cyto_names(x)
    } else {
      spill <- lapply(cyto_names(x), function(y) {
        sm <- cyto_spillover_extract(x[[y]])[[1]]
        if(is.null(sm)){
          stop(paste0("Unable to extract spillover matrix from ",
                      cyto_names(x[[y]]), "."))
        }
        return(sm)
      })
    }
  }

  print(spill)
  
  # Channels
  fluor_channels <- cyto_fluor_channels(x)
  
  # Spillover may contain more channels than in samples
  spill <- lapply(spill, function(z){
    # Select rows - square matrix
    if(nrow(z) == ncol(z)){
      z <- z[match(fluor_channels, colnames(z)), ]
    }
    # Select columns
    z <- z[, fluor_channels]
    return(z)
  })
  
  # Apply compensation
  flowCore::compensate(x, spill)
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowFrame <- function(x,
                                      spillover = NULL,
                                      select = 1,
                                      ...) {

  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (is(spillover, "character")) {
      # No file extension
      if (file_ext(spillover) == "") {
        spillover <- paste0(spillover, ".csv")
      }
      # Check working directory for csv file
      if (getOption("CytoExploreR_wd_check")) {
        if (!file_wd_check(spillover)) {
          stop(paste(spillover, "does not exist in this working directory."))
        }
      }
      spill <- read.csv(spillover, header = TRUE, row.names = 1)
      colnames(spill) <- rownames(spill)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% BiocGenerics::colnames(x)) |
        !all(rownames(spill) %in% BiocGenerics::colnames(x))) {
        stop(
          paste(
            "'spillover' must have valid fluorescent channels as rownames",
            "and colnames."
          )
        )
      }
      # spillover is a matrix/data.frame
    } else if (is(spillover, "matrix") |
      is(spillover, "data.frame")) {
      # column names must be valid channels (rownames not essential)
      if (!all(colnames(spillover) %in% BiocGenerics::colnames(x))) {
        stop("'spillover' must have valid fluorescent channels as colnames.")
      } else {
        spill <- spillover
      }
    }
    # Extract spillover matrix directly from x
  } else if (is.null(spillover)) {
    spill <- cyto_spillover_extract(x)[[1]]
    if(is.null(spill)){
      stop(paste0("Unable to extract spillover matrix from", 
                  cyto_names(x), "."))
    }
  }
  
  # Channels
  fluor_channels <- cyto_fluor_channels(x)
  
  # Select rows - square matrix
  if(nrow(spill) == ncol(spill)){
    spill <- spill[match(fluor_channels, colnames(spill)), ]
  }
  # Select columns
  spill <- spill[, fluor_channels]
  
  # Apply compensation
  flowCore::compensate(x, spill)
}

## CYTO_NODES ------------------------------------------------------------------

#' Extract Names of Gated Populations in GatingHierarchy or GatingSet
#'
#' \code{cyto_nodes} is simply an autocomplete-friendly wrapper for
#' \code{\link[flowWorkspace:gs_get_pop_paths]{gs_get_pop_paths}}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:gs_get_pop_paths]{gs_get_pop_paths}}.
#'
#' @return character vector of gated node/population names.
#'
#' @importFrom flowWorkspace gh_get_pop_paths gs_get_pop_paths
#' @importFrom methods is
#'
#' @export
cyto_nodes <- function(x, ...) {

  # Make call to gh/gs_get_pop_paths
  if (is(x, "GatingHierarchy")) {
    gh_get_pop_paths(x, ...)
  } else if (is(x, "GatingSet")) {
    gs_get_pop_paths(x, ...)
  }
}

## CYTO_CHANNEL_MATCH ----------------------------------------------------------

#' Table Editor for Channel Match File Construction
#'
#' @param x object of \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param save_as name to use for the saved channel match csv file, set to
#'   \code{"date-Channel-Match.csv"}.
#'
#' @return save constructed channel_match csv file.
#'
#' @importFrom utils write.csv edit
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#' # Load in CytoExploreRData to access data
#' library(CytoExploreR)
#'
#' # Generate channel match file for compensation controls
#' cyto_channel_match(Compensation)
#' }
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_channel_match <- function(x,
                               save_as = NULL) {

  # Set default name for channel_match file
  if (is.null(save_as)) {
    save_as <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Channel-Match.csv"
    )
    # save_as must contain csv file extension
  } else {
    # File extension missing
    if (file_ext(save_as) != "csv") {
      save_as <- paste0(save_as, ".csv")
    }
  }

  # Extract sample names
  nms <- cyto_names(x)

  # Construct data.frame for editing
  channel_match <- data.frame("name" = nms, "channel" = rep("NA", length(nms)))
  colnames(channel_match) <- c("name", "channel")
  rownames(channel_match) <- NULL

  # Edit channel_match
  channel_match <- suppressWarnings(edit(channel_match))

  # Check that all channels are valid or throw an error
  if (!all(channel_match$channel %in% c("Unstained", cyto_fluor_channels(x)))) {
    stop("Some inputs in the channel column are not valid.")
  }

  # Write edited channel match file to csv file
  write.csv(channel_match, save_as, row.names = FALSE)

  # Return edited channel match file
  return(channel_match)
}

## CYTO_EMPTY ------------------------------------------------------------------

#' Construct an empty flowFrame
#'
#' @param name name to add to the constructed flowFrame.
#' @param channels channels to include in the constructed flowFrame.
#' @param ... additional arguments passed to
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#'
#' @importFrom flowCore identifier<- flowFrame
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#'
#' # Construct empty flowFrame
#' cyto_empty(name = "Test.csv", channels = c("FSC-A", "SSC-A", "PE-A"))
#' 
#' @export
cyto_empty <- function(name = NULL,
                       channels = NULL, ...) {

  # CHANNELS
  if (is.null(channels)) {
    stop("Supply the names of the channels to include in the flowFrame.")
  }

  # CONSTRUCT EMPTY FLOWFRAME
  empty_flowFrame <- matrix(0,
    ncol = length(channels),
    nrow = 1,
    byrow = TRUE
  )
  colnames(empty_flowFrame) <- channels
  empty_flowFrame <- flowFrame(empty_flowFrame, ...)
  empty_flowFrame <- empty_flowFrame[-1, ]

  # NAME
  if (!is.null(name)) {
    identifier(empty_flowFrame) <- name
  }

  # RETURN EMPTY FLOWFRAME
  return(empty_flowFrame)
}

## CYTO_COPY -------------------------------------------------------------------

#' Copy a cytoset, cytoframe or GatingSet
#' 
#' @param x cytoframe, cytoset or GatingSet to be copied.
#' 
#' @return copied cytoframe, cytoset or GatingSet.
#' 
#' @importFrom flowWorkspace gs_clone realize_view
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
cyto_copy <- function(x){
  
  # GatingSet
  if(is(x, "GatingSet")){
    x <- gs_clone(x)
  }
  
  # cytoset
  if(is(x, "cytoset") | is(x, "cyto_frame")){
    x <- realize_view(x)
  }
  
  # RETURN COPY
  return(x)
  
}

## CYTO_SPILLOVER_EXTRACT ------------------------------------------------------

#' Extract spillover matrix from cytometry object
#' 
#' @param x object of class flowFrame, flowSet, GatingHierachy or GatingSet.
#' 
#' @return list of spillover matrices or NULL.
#' 
#' @importFrom methods is
#' @importFrom flowWorkspace gh_get_compensations gs_get_compensations
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @examples 
#' library(CytoExploreRData)
#' 
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#' 
#' # Apply compensation
#' gs <- cyto_compensate(gs)
#' 
#' # Extract spillover matrices
#' spill <- cyto_spillover_extract(gs)
#' 
#' @export
cyto_spillover_extract <- function(x){
  
  # GATINGSET
  if(is(x, "GatingSet")){
    spill <- gs_get_compensations(x)
    if(all(LAPPLY(spill, "is.null"))){
      spill <- NULL
    }else{
      spill <- LAPPLY(spill, function(z){z@spillover})
    }
  # GATINGHIERACHY  
  }else if(is(x, "GatingHierarchy")){
    spill <- gh_get_compensations(x)
    if(!is.null(spill)){
      spill <- list(spill@spillover)
      names(spill) <- cyto_names(x)
    }
  # FLOWSET  
  }else if(is(x, "flowSet")){
    spill <- lapply(seq_along(x), function(z){
      # CyTOF lacks spill slot (just in case)
      tryCatch(x[[z]]@description$SPILL, error = function(e){NULL})
    })
    names(spill) <- cyto_names(x)
    if(all(LAPPLY(spill, "is.null"))){
      spill <- NULL
    }
  # FLOWFRAME  
  }else if(is(x, "flowFrame")){
    spill <- tryCatch(x@description$SPILL, error = function(e){NULL})
    if(!is.null(spill)){
      spill <- list(spill)
      names(spill) <- cyto_names(x)
    }
  }
  
  # RETURN LIST OF SPILLOVER MATRICES
  return(spill)
}
