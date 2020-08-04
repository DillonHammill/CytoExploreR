## CYTO_IMPORT -----------------------------------------------------------------

#' Import cytometry data analyzed using other platforms
#' @param path path to xml or wsp file and associated fcs files.
#' @noRd
cyto_import <- function(path = ".",
                        type = "flowJo",
                        select = NULL,
                        exclude = NULL,
                        sort = TRUE,
                        barcode = FALSE,
                        restrict = FALSE,
                        markers = TRUE,
                        details = TRUE,
                        gatingTemplate = NULL, ...) {
  
  # CytoML EXISTS
  if (requireNamespace("CytoML")) {
    
    # FILES IN DIRECTORY
    file_paths <- list.files(path,
                             full.names = TRUE)
    
    # FILE NAMES
    file_names <- basename(file_paths)
    
    # FILE EXTENSIONS
    file_ext <- file_ext(file_names)
    
    # FCS FILES
    fcs_files <- file_paths[!file_ext %in% c("", "fcs", "FCS")]
    
    # IMPORT CYTOBANK TO GATINGSET
    if(grepl("cytobank", type, ignore.case = TRUE)){
      # ACS
      if(any(grepl("acs", file_ext, ignore.case = TRUE))) {
        acs_file <- file_paths[which(file_ext == "acs")]
        cytobank_exp <- CytoML::open_cytobank_experiment(acs_file)
        gs <- CytoML::cytobank_to_gatingset(cytobank_exp)
      # XML 
      }else if(any(grepl("xml", file_ext, ignore.case = TRUE))){
        xml_file <- file_paths[which(grepl(file_ext, 
                                           "xml", 
                                           ignore.case = TRUE))]
        gs <- CytoML::cytobank_to_gatingset(xml_file, fcs_files)
      }
    # IMPORT DIVA TO GATINGSET
    }else if(grepl("diva", type, ignore.case = TRUE)){
      xml_file <- file_paths[which(grepl(file_ext, 
                                         "xml", 
                                         ignore.case = TRUE))]
      diva_ws <- CytoML::open_diva_xml(xml_file)
      gs <- CytoML::diva_to_gatingset(diva_ws, fcs_files)
    # IMPORT FLOWJO TO GATINGSET
    }else if(grepl("flowjo", type, ignore.case = TRUE)){
      # WPS
      if(any(grepl("wsp", file_ext, ignore.case = TRUE))){
        wps_file <- file_paths[which(grepl(file_ext, 
                                           "wps", 
                                           ignore.case = TRUE))]
        gs <- CytoML::flowjo_to_gatingset(wps_file, path = path)
      # XML
      }else if(any(grepl("xml", file_ext, ignore.case = TRUE))){
        xml_file <- file_paths[which(grepl(file_ext, 
                                           "xml", 
                                           ignore.case = TRUE))]
        flowjo_ws <- CytoML::open_flowjo_xml(xml_file)
        gs <- CytoML::flowjo_to_gatingset(flowjo_ws, fcs_files)
      }
    }
    
    # MARKERS
    if (markers != FALSE) {
      message("Assigning markers to channels...")
      # DEFAULT FILE NAME
      if (markers == TRUE) {
        gs <- cyto_markers_edit(gs)
      } else {
        gs <- cyto_markers_edit(gs,
                                file = markers
        )
      }
    }
    
    # EXPERIMENT DETAILS
    if (details != FALSE) {
      message("Updating experiment details...")
      if (details == TRUE) {
        gs <- cyto_details_edit(gs)
      } else {
        gs <- cyto_details_edit(gs,
                                file = details
        )
      }
    }
    
    # GENERATE GATINGTEMPLATE
    if(!is.null(gatingTemplate)){
      # FILE EXTENSION
      if(file_ext(gatingTemplate) != "csv"){
        gatingTemplate <- paste0(gatingTemplate, ".csv")
      }
      # CREATE GATINGTEMPLATE
      message("Creating CytoExploreR gatingTemplate...")
      cyto_gatingTemplate_generate(gs, gatingTemplate)
    }
    
    # RETURN GATINGSET
    return(gs)
    
  # CytoML MISSING
  } else {
    stop(paste0(
      "cyto_import requires the CytoML package ",
      "- BiocManager::install('CytoML')"
    ))
  }
}

## CYTO_EXPORT -----------------------------------------------------------------

#' Export cytometry data for use in flowJo or Cytobank
#'
#' Simply a wrapper around CytoML \code{gatingset_to_cytobank} and
#' \code{gatingset_to _flowjo} to export your CytoExploreR analyses in a format
#' accepted by these platforms. Users will need to have \code{docker} running on
#' their computers to export flowJo workspace files.
#'
#' @param x object of class GatingSet to export.
#' @param save_as the name of a flowJo workspace file with .wsp extension or a
#'   Cytobank xml file with .xml extension.
#' @param ... additional arguments passed to \code{gatingset_to_cytobank} or
#'   \code{gatingset_to_flowjo}.
#'
#' @importFrom tools file_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#' library(CytoML)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transformations
#' gs <- cyto_transform(gs)
#'
#' # Gating
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#'
#' # Export to cytobank xml
#' cyto_export(gs, "cytobank.xml")
#'
#' # Export to flowjo workspace
#' cyto_export(gs, "flowjo.wsp")
#' }
#'
#' @export
cyto_export <- function(x,
                        save_as = NULL,
                        ...) {
  
  # CytoML
  if (requireNamespace("CytoML")) {
    
    # FILE NAME
    if (is.null(save_as)) {
      stop("Supply either a wsp or xml file name to 'save_as'.")
    } else {
      # FLOWJO EXPORT BY DEFAULT
      if (.empty(file_ext(save_as))) {
        save_as <- paste0(save_as, ".wsp")
      }
    }
    
    # SAVE
    if (file_ext(save_as) == "xml") {
      message("Saving GatingSet to Cytobank XML file...")
      CytoML::gatingset_to_cytobank(x, 
                                    save_as,
                                    ...)
    } else if (file_ext(save_as) == "wsp") {
      message("Saving GatingSet to flowJo workspace file...")
      CytoML::gatingset_to_flowjo(x, 
                                  save_as,
                                  ...)
    }
    
    # CytoML MISSING
  } else {
    stop(paste0(
      "cyto_export requires the CytoML package ",
      "- BiocManager::install('CytoML')"
    ))
  }
}

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
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_load <- function(path = ".",
                      select = NULL,
                      exclude = NULL,
                      sort = TRUE,
                      barcode = FALSE,
                      restrict = FALSE, ...) {

  # DIRECTORY NOT FOUND
  if (!dir.exists(path)) {
    stop("Specified path does not exist.")
  }

  # FILE PATHS
  files <- list.files(path, full.names = TRUE)

  # VALID FILES
  files_ext <- file_ext(files)
  files_ind <- which(!files_ext %in% c("", "fcs", "FCS"))
  
  # NO VALID FILES
  if(length(files_ind) == length(files)){
    stop(paste0(path,
                " does not contain any valid FCS files."))
  }
  
  # EXCLUDE IRRELEVANT FILES
  if(length(files_ind) > 0){
    files <- files[-files_ind]
  }
  
  # SELECT
  if (!is.null(select)) {
    file_ind <- c()
    lapply(select, function(z) {
      if (any(grepl(z, files, ignore.case = TRUE))) {
        file_ind <<- c(file_ind,
                       which(grepl(z, files, ignore.case = TRUE)))
      }
    })
    files <- files[unique(file_ind)]
  }

  # EXCLUDE
  if (!is.null(exclude)) {
    file_ind <- c()
    lapply(exclude, function(z) {
      if (any(grepl(z, files, ignore.case = TRUE))) {
        file_ind <<- c(file_ind,
                       which(grepl(z, files, ignore.case = TRUE)))
      }
    })
    files <- files[-unique(file_ind)]
  }

  # SORTED FILE PATHS
  if (sort) {
    files <- mixedsort(files)
  }

  # SAVED GATINGSET
  if ("pb" %in% file_ext(files)) {
    # LOAD GATINGSET
    x <- load_gs(path = path)
    # FCS FILES
  } else {
    # CYTOSET
    x <- load_cytoset_from_fcs(files = normalizePath(files), ...)

    # CORRECT GUID SLOTS - NECESSARY?
    nms <- cyto_names(x)
    lapply(seq_len(length(nms)), function(z) {
      suppressMessages(identifier(x[[z]]) <<- nms[z])
    })

    # BARCODING
    if (barcode) {
      x <- cyto_barcode(x)
    }

    # CHANNEL RESTRICTION
    if (restrict) {
      x <- cyto_channels_restrict(x)
    }
  }

  # RETURN CYTOSET
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
#' @references Monaco,G. et al. (2016) flowAI: automatic and interactive anomaly
#'   discerning tools for flow cytometry data. Bioinformatics. 2016 Aug
#'   15;32(16):2473-80.
#'   \url{https://academic.oup.com/bioinformatics/article/32/16/2473/2240408}
#' @seealso \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}
#'
#' @export
cyto_clean <- function(x, ...) {

  # GATINGSET/GATINGHIERARCHY
  if (is(x, "GatingSet") | is(x, "GatingHierarchy")) {
    # PARENT
    parent <- cyto_nodes(x, path = "auto")[1]
    # EXTRACT DATA
    cyto_data <- cyto_extract(x, parent)
    # flowAI REQUIRES FLOWSET
    if (is(cyto_data, "cytoset")) {
      cyto_data <- cytoset_to_flowSet(cyto_data) # REMOVE
    }
    # flowAI messes with experiment details :(
    if (is(cyto_data, "flowSet")) {
      pd <- cyto_details(cyto_data)
    }
    # CLEAN DATA
    invisible(capture.output(cyto_data <- flow_auto_qc(cyto_data,
      html_report = FALSE,
      mini_report = FALSE,
      fcs_QC = FALSE,
      folder_results = FALSE,
      ...
    )))

    # RETURN CYTOSET
    if (is(cyto_data, "flowSet")) {
      cyto_data <- flowSet_to_cytoset(cyto_data)
      cyto_details(cyto_data) <- pd
    }
    # REPLACE DATA
    gs_cyto_data(x) <- cyto_data
  } else {
    # FLOWSET REQUIRED
    if (is(x, "cytoset")) {
      x <- cytoset_to_flowSet(x)
    }
    # flowAI messes with experiment details :(
    if (is(x, "flowSet")) {
      pd <- cyto_details(x)
    }
    invisible(capture.output(x <- flow_auto_qc(x,
      html_report = FALSE,
      mini_report = FALSE,
      fcs_QC = FALSE,
      folder_results = FALSE,
      ...
    )))
    # RETURN CYTOSET
    if (is(x, "flowSet")) {
      x <- flowSet_to_cytoset(x)
      cyto_details(x) <- pd
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
#' Users can optionally provide a name for a gatingTemplate csv file which will
#' be created if necessary and assigned as the active gatingTemplate.
#'
#' @param path points to the location of the .fcs files to read in (e.g. name of
#'   a folder in current working directory).
#' @param gatingTemplate name of a gatingTemplate csv file to be used for gate
#'   saving.
#' @param restrict logical indicating whether unassigned channels should be
#'   dropped from the returned cytoset, set to FALSE by default. Alternatively,
#'   users can supply a vector of channels to remove, see
#'   \code{\link{cyto_channels_restrict}} for details.
#' @param clean logical indicating whether the loaded data should be cleaned
#'   using \code{cyto_clean}, set to FALSE by default. Alternatively, users can
#'   indicate which types of anomalies should be checked as expected by
#'   \code{remove_from} in \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}.
#' @param markers logical indicating whether a call should be made to
#'   \code{cyto_markers_edit} to update the markers associated with channels in
#'   the loaded sampes, set to TRUE by default. The name of the csv to which
#'   these details will be supplied can also be passed to this argument.
#' @param parse_names logical indicating whether the file names should be parsed
#'   into experiment details using \code{cyto_names_parse}, set to FALSE by
#'   default. If you need to parse the names using a different dlimiter, supply
#'   the delimiter to this argument instead of TRUE.
#' @param details logical indicating whether a call should be made to
#'   \code{cyto_details_edit} to update the experimental details associated with
#'   the loaded samples, set to TRUE by default. The name of the csv to which
#'   these details will be supplied can also be passed to this argument.
#' @param sample logical indicating whether all samples should be downsampled to
#'   the minimum number of events in a sample, set to FALSE by default.
#'   Alternatively, users can supply a numeric to indicate the desired number of
#'   events to keep in each sample.
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
#' @seealso \code{\link{cyto_load}}
#' @seealso \code{\link{cyto_markers_edit}}
#' @seealso \code{\link{cyto_names_parse}}
#' @seealso \code{\link{cyto_details_edit}}
#' @seealso \code{\link{cyto_channels_restrict}}
#' @seealso \code{\link{cyto_clean}}
#' @seealso \code{\link{cyto_gatingTemplate_select}}
#' @seealso \code{\link{cyto_gatingTemplate_create}}
#'
#' @export
cyto_setup <- function(path = ".",
                       gatingTemplate = NULL,
                       restrict = FALSE,
                       clean = FALSE,
                       markers = TRUE,
                       parse_names = FALSE,
                       details = TRUE,
                       sample = FALSE, ...) {
  
  # CYTOSET/GATINGSET
  message("Loading FCS files into a GatingSet...")
  x <- cyto_load(path = path, restrict = FALSE, ...)
  
  # MARKERS
  if (markers != FALSE) {
    message("Assigning markers to channels...")
    # DEFAULT FILE NAME
    if (markers == TRUE) {
      x <- cyto_markers_edit(x)
    } else {
      x <- cyto_markers_edit(x,
                             file = markers
      )
    }
  }
  
  # PARSE NAMES
  if(parse_names != FALSE){
    message("Parsing file names into experiment details...")
    if(parse_names == TRUE){
      x <- cyto_names_parse(x)
    }else{
      x <- cyto_names_parse(x,
                            split = parse_names)
    }
  }
  
  # EXPERIMENT DETAILS
  if (details != FALSE) {
    message("Updating experiment details...")
    if (details == TRUE) {
      x <- cyto_details_edit(x)
    } else {
      x <- cyto_details_edit(x,
                             file = details
      )
    }
  }
  
  # FLOWSET LOADED
  if (is(x, "flowSet")) {
    # RESTRICT CHANNELS
    if (restrict != FALSE) {
      if(restrict == TRUE){
        restrict <- NULL
      }
      message("Removing unassigned channels...")
      x <- cyto_channels_restrict(x,
                                  exclude = restrict)
    }
    # CLEAN DATA
    if (clean != FALSE) {
      if(clean == TRUE){
        clean <- "all"
      }
      message("Cleaning data to remove anomalies...")
      x <- cyto_clean(x,
                      remove_from = clean)
    }
    # SAMPLING
    if(sample != FALSE){
      if(sample == TRUE){
        sample <- min(
          cyto_apply(x, 
                     .cyto_count)
        )
      }
      message(
        paste("Downsampling each sample to",
              sample, "events.")
      )
      x <- cyto_sample(x, 
                       display = sample,
                       seed = 56)
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
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
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
    # pd$name <- factor(pd$name, levels = pd$name)
    return(pd)
  }
}

# CYTO_DETAILS REPLACEMENT METHOD ----------------------------------------------

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

# CYTO_NAMES REPLACEMENT METHOD ------------------------------------------------

#' Replacement method for cyto_names
#'
#' @param x object of class flowFrame, flowSet, GatingHierarchy or GatingSet.
#' @param value vector of replacement names.
#'
#' @importFrom flowWorkspace sampleNames<-
#' @importFrom flowCore identifier<-
#'
#' @examples
#'
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Sample names
#' cyto_names(fs)
#'
#' # Change first sample name
#' cyto_names(fs)[1] <- "first_sample"
#'
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#'
#' # Change last sample name
#' cyto_names(gs)[length(gs)] <- "last_sample"
#'
#' # Updated sample names
#' cyto_names(gs)
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
"cyto_names<-" <- function(x, value) {
  if (is(x, "flowSet") |
    is(x, "GatingHierarchy") |
    is(x, "GatingSet")) {
    sampleNames(x) <- value
  } else if (is(x, "flowFrame")) {
    identifier(x) <- value
  }
  return(x)
}

## CYTO_NAMES_PARSE ------------------------------------------------------------

#' Parse file names into experiment details
#'
#' Reward users that use consistent names for files by parsing file names into
#' new variables and adding them to \code{cyto_details}. Variable construction
#' may not be perfect if certain files are named differently, but it will give a
#' starting point for editing these values using \code{cyto_details_edit}
#' without having to manually enter a lot of these details.
#'
#' @param x object of class flowSet or GatingSet.
#' @param vars vector containing the names of the variables to be added to
#'   \code{cyto_details}, set to var_1, var_2 ... by default.
#' @param split delimiter to split the file name into fragments, set to "_" by
#'   default.
#' @param exclude vector of indices indicating which text fragments should not be
#'   included in the experiment details, set to NULL by default.
#' @param ... additional arguments passed to
#'   \code{\link[base:strsplit]{strsplit}}.
#'
#' @return flowSet or GatingSet with \code{cyto_details} updated with new
#'   variables extracted from the file names.
#'
#' @importFrom tools file_path_sans_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Parse file names to variables
#' fs <- cyto_names_parse(Activation,
#' vars = c("sample_type", "sample_id"),
#' split = "_")
#'
#' # Updated experiment details
#' cyto_details(fs)
#'
#' @export
cyto_names_parse <- function(x,
                            vars = NULL,
                            split = "_",
                            exclude = NULL,
                            ...) {
  
  # NAMES
  cyto_names <- cyto_names(x)
  
  # STRIP EXTENSION
  cyto_names <- file_path_sans_ext(cyto_names)
  
  # SPLIT
  cyto_names_split <- strsplit(cyto_names, 
                               split,
                               ...)
  
  # SKIP
  if(!is.null(exclude)){
    cyto_names_split <- lapply(cyto_names_split, function(z){
      z[-exclude]
    })
  }
  
  # REQUIRED VARIABLE LENGTH
  var_length <- max(LAPPLY(cyto_names_split, "length"))
  
  # VARIABLES
  if(is.null(vars)){
    vars <- paste0("var_", seq_len(var_length))
  }else{
    if(length(vars) != var_length){
      stop(paste0("Require ",
                  var_length, 
                  " variable names to 'vars' to update cyto_details."))
    }
  }  
  
  # FILL WITH NA
  cyto_names_split <- lapply(cyto_names_split, function(z){
    if(length(z) < var_length){
      z <- rep(c(z, rep(NA, var_length)), length.out = var_length)
      return(z)
    }else{
      return(z)
    }
  })
  cyto_names_split <- do.call("rbind", cyto_names_split)
  colnames(cyto_names_split) <- vars
  
  # UPDATE CYTO_DETAILS
  cyto_details(x) <- cbind(cyto_details(x), cyto_names_split)
  
  # RETURN FLOWSET/GATINGSET
  return(x)
  
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
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param ... additional arguments passed to \code{\link{cyto_transformer_log}},
#'   \code{\link{cyto_transformer_arcsinh}}, \code{\link{cyto_transformer_biex}}
#'   or \code{\link{cyto_transformer_logicle}}, when no \code{trans} object is
#'   supplied.
#'
#' @return object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet} with transformations applied.
#'
#' @importFrom flowWorkspace recompute gh_get_pop_paths gh_pop_get_parent
#'   gh_pop_get_gate gh_pop_set_gate gs_pop_set_gate
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

  # TRANSFORM FLOWFRAME OR FLOWSET
  if (is(x, "flowFrame") | is(x, "flowSet")) {

    # Extract transformations from transformerList to transformList
    transform_list <- cyto_transform_extract(transformer_list,
      inverse = inverse
    )

    # Apply transformations
    x <- suppressMessages(transform(x, transform_list))

    # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (is(x, "GatingHierarchy") | is(x, "GatingSet")) {

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
          c(
            n2mfrow(n)[1],
            n2mfrow(n)[2]
          )
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

  # TRANSFORM FLOWFRAME OR FLOWSET
  if (is(x, "flowFrame") | is(x, "flowSet")) {

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
          c(
            n2mfrow(n)[1],
            n2mfrow(n)[2]
          )
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

  # TRANSFORM FLOWFRAME OR FLOWSET
  if (is(x, "flowFrame") | is(x, "flowSet")) {

    # Extract transformations to transformList
    transform_list <- cyto_transform_extract(trans, inverse = inverse)

    # Apply transformations
    x <- suppressMessages(transform(x, transform_list))


    # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (is(x, "GatingHierarchy") | is(x, "GatingSet")) {

    # Inverse transformations not yet supported
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
          c(
            n2mfrow(n)[1],
            n2mfrow(n)[2]
          )
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
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}}.
#' @param copy logical indicating whether a deep copy of the extracted data
#'   should be returned.
#' @param raw logical indicating whether a list of raw data matrices should be
#'   returned instead of a flowFrame or flowSet.
#' @param channels names of the markers or channels for which data should be
#'   extracted, set to all channels by default.
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
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_extract <- function(x,
                         parent = NULL,
                         select = NULL,
                         copy = FALSE,
                         raw = FALSE,
                         channels = NULL,
                         ...) {
  
  # DEFAULT PARENT
  if (is.null(parent)) {
    parent <- cyto_nodes(x, path = "auto")[1]
  }

  # EXTRACT
  if (is(x, "GatingHierarchy")) {
    x <- gh_pop_get_data(x, parent, ...)
  } else if (is(x, "GatingSet")) {
    x <- gs_pop_get_data(x, parent, ...)
  }

  # COPY
  if (copy) {
    x <- cyto_copy(x)
  }

  # SELECT
  if(!is.null(select)){
    x <- cyto_select(x, select)
  }
  
  # RESTRICT
  if(!is.null(channels)){
    channels <- cyto_channels_extract(x, channels = channels)
    x <- x[, channels]
  }
  
  # RAW DATA MATRICES
  if (raw) {
    nms <- cyto_names(x)
    if (is(x, "flowFrame")) {
      x <- list(exprs(x))
    } else if (is(x, "flowSet")) {
      y <- lapply(seq_along(x), function(z) {
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
#' # Select Stim-C samples with 0 and 500 nM OVA concentrations
#' fs <- cyto_filter(
#'   Activation,
#'   Treatment == "Stim-C",
#'   OVAConc %in% c(0, 500)
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
#' # Select Stim-C samples with 0 and 500 nM OVA concentrations
#' fs <- cyto_select(Activation,
#'   Treatment = "Stim-C",
#'   OVAConc = c(0, 500)
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
cyto_select <- function(x, ...) {

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

  # INDICES SUPPLIED
  if (length(args) == 1 &
    (is.null(names(args)) | .empty(names(args))) &
    is.numeric(unlist(args))) {
    # INDICES TO SELECT
    ind <- unlist(args)
  } else {
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
        pd_filter <<- pd_filter[ind, , drop = FALSE]
      }
    })

    # Get indices for selected samples
    ind <- match(pd_filter[, "name"], pd[, "name"])
  }

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
#' @param group_by names of cyto_details variables to use for merging, set to
#'   "all" to group all samples in \code{x}. The order of the grouping can be
#'   controlled by specifying the factor levels in a list (e.g. list(Treatment =
#'   c("Stim-A","Stim-C","Stim-B", "Stim-D"))).
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

  # group_by is a list with factor levels - should not be "all"
  if (is(group_by, "list")) {
    # Check variables and factor levels
    lapply(seq_along(group_by), function(z) {
      # Variable
      var <- names(group_by)[z]
      # Expected variable levels
      var_levels <- unique(pd[, var])
      # Watch out for NA
      if (any(LAPPLY(group_by[[z]], "is.na"))) {
        ind <- which(LAPPLY(group_by[[z]], "is.na"))
        group_by[[z]][ind] <- "NA"
      }
      # Check variables
      if (!var %in% colnames(pd)) {
        stop(paste0(
          var,
          " is not a valid variable for this ",
          class(x), "."
        ))
      }
      # Incorrect factor levels
      if (!all(group_by[[z]] %in% unique(pd[, var]))) {
        lapply(group_by[[z]], function(y) {
          if (!y %in% unique(pd[, var])) {
            stop(paste0(
              y, " is not a valid factor level for ",
              group_by[[z]],
              "."
            ))
          }
        })
      }
      # Update factor levels in pd
      if (!all(var_levels %in% group_by[[z]])) {
        missing_levels <- as.vector(var_levels[!var_levels %in% group_by[[z]]])
        group_by[[z]] <<- c(
          group_by[[z]],
          missing_levels
        )
      }
      # Convert pd variable to factor and set levels
      pd[, var] <<- factor(pd[, var], levels = group_by[[z]])
    })
    # Convert group_by to vector
    group_by <- names(group_by)
    # group_by is a vector of variable names
  } else {
    # Check variables
    if (group_by[1] != "all" & !all(group_by %in% colnames(pd))) {
      lapply(group_by, function(y) {
        if (!y %in% colnames(pd)) {
          stop(paste0(y, " is not a valid variable for this ", class(x), "."))
        }
      })
    }
  }

  # Split pd based on group_by into a named list
  if (length(group_by) == 1) {
    if (group_by == "all") {
      pd_split <- list("all" = pd)
    } else if (group_by == "name") {
      pd_split <- lapply(nms, function(z) {
        pd[pd$name == z, , drop = FALSE]
      })
      names(pd_split) <- nms
    } else {
      pd_split <- split(pd, pd[, group_by],
        sep = " ",
        lex.order = TRUE,
        drop = TRUE
      )
    }
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
  if ("all" %in% grps) {
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
#' @param inverse logical indicating whether the data should be
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
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Save each flowFrame to file
#' cyto_save(fs, save_as = "Samples")
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
                                inverse = FALSE,
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
      copy = TRUE
    )
    # TRANSFORMATIONS
    trans <- cyto_transformer_extract(x)
    # FLOWSET METHOD
    fr_list <- cyto_save(
      x = fs,
      split = split,
      names = names,
      save_as = save_as,
      inverse = inverse,
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
                                      inverse = FALSE,
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
      copy = TRUE
    )
    # TRANSFORMATIONS
    trans <- cyto_transformer_extract(x)
    # FLOWSET METHOD
    fr_list <- cyto_save(
      x = fr,
      split = split,
      names = names,
      save_as = save_as,
      inverse = inverse,
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
                              inverse = FALSE,
                              trans = NULL,
                              ...) {

  # COPY
  x <- cyto_copy(x)

  # LIST OF FLOWFRAMES
  fr_list <- cyto_convert(x, "list of flowFrames")

  # LIST OF SPLIT FLOWFRAMES
  if (split == TRUE) {
    # NAMES SUPPLIED - CHECK LENGTH
    if (!is.null(names)) {
      # SAMPLES PER FLOWFRAME
      samples_per_file <- lapply(fr_list, function(z) {
        length(unique(exprs(z)[, "Sample ID"]))
      })
      # SPLIT NAMES SUPPLIED
      samples <- sum(unlist(samples_per_file))
      if (length(names) != samples) {
        stop(paste("Expecting", samples, "names for the split files."))
      }
      # PREPARE NAMES
      ind <- LAPPLY(seq_along(fr_list), function(z) {
        rep(z, samples_per_file[[z]])
      })
      names <- split(names, ind)
      # SPLIT FR_LIST
      fr_list <- mapply(function(z, name) {
        cyto_split(z, names = name)
      }, fr_list, names)
      fr_list <- unlist(fr_list)
      # NO NAMES SUPPLIED
    } else {
      fr_list <- LAPPLY(fr_list, function(z) {
        cyto_split(z)
      })
    }
  }

  # DIRECTORY CHECK
  if (!is.null(save_as) & dir.exists(save_as)) {
    # FILES WILL BE OVERWRITTEN
    if (any(list.files(save_as) %in% cyto_names(fr_list))) {
      message(paste0("Files will be overwritten in ", save_as, "."))
      opt <- readline("Do you want to continue? (Y/N)")
      if (grepl("n", opt, ignore.case = TRUE)) {
        return(NULL)
      }
    }
  }

  # MESSAGE
  if (is.null(save_as)) {
    location <- "current working directory."
  } else {
    location <- save_as
  }
  message(paste0("Writing FCS files to ", location, "..."))

  # WRITE FCS FILES
  fr_list <- lapply(fr_list, function(z) {
    # INVERSE TRANSFORM
    if (inverse == TRUE) {
      # TRANSFORMERS REQUIRED
      if (is.null(trans) | .all_na(trans)) {
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
                                inverse = FALSE,
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
      message(paste0("Files will be overwritten in ", save_as, "."))
      opt <- readline("Do you want to continue? (Y/N)")
      if (grepl("n", opt, ignore.case = TRUE)) {
        return(NULL)
      }
    }
  }

  # MESSAGE
  if (is.null(save_as)) {
    location <- "current working directory."
  } else {
    location <- save_as
  }
  message(paste0("Writing FCS files to ", location, "..."))

  # WRITE FCS FILES
  fr_list <- lapply(fr_list, function(z) {
    # INVERSE TRANSFORM
    if (inverse == TRUE) {
      # TRANSFORMERS REQUIRED
      if (is.null(trans) | .all_na(trans)) {
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
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param display can be either a numeric [0,1] or integer to indicate the
#'   percentage or number of events to keep respectively.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
#' @param plot logical required for lists to indicate whether sampling should be
#'   scaled per flowFrame to retain original ratios, as used in cyto_plot.
#' @param ... not in use.
#'
#' @return object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} sampled to
#'   \code{display} events.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowCore sampleFilter Subset flowSet exprs
#' @importFrom methods is
#' @importFrom flowWorkspace gs_cyto_data gs_cyto_data<- cytoset
#'   flowFrame_to_cytoframe flowSet_to_cytoset recompute
#'
#' @examples
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
cyto_sample.GatingHierarchy <- function(x,
                                        display = 1,
                                        seed = NULL,
                                        ...) {
  
  # EXTRACT CYTOSET
  cs <- gs_cyto_data(x)
  
  # SAMPLING
  cs <- cyto_sample(cs,
                    display = display,
                    seed = seed,
                    ...)
  
  # REPLACE DATA
  gs_cyto_data(x) <- cs
  
  # RECOMPUTE
  suppressMessages(recompute(x))
  
  # RETURN GATINGHIERACHY
  return(x)
  
}

#' @rdname cyto_sample
#' @export
cyto_sample.GatingSet <- function(x,
                                  display = 1,
                                  seed = NULL,
                                  ...){
  
  # EXTRACT CYTOSET
  cs <- gs_cyto_data(x)
  
  # SAMPLING
  cs <- cyto_sample(cs,
                    display = display,
                    seed = seed,
                    ...)
  
  # REPLACE DATA
  gs_cyto_data(x) <- cs
  
  # RECOMPUTE
  suppressMessages(recompute(x))
  
  # RETURN GATINGSET
  return(x)
  
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
  
  # FLOWSET/CYTOSET
  cs <- cyto_apply(x, 
                   cyto_sample, 
                   display = display, 
                   seed = seed, 
                   ...)
  
  # RETURN FLOWSET/CYTOSET
  return(cs)
  
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
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace flowSet_to_cytoset gs_cyto_data
#'
#' @export
cyto_beads_sample <- function(x,
                              beads = NULL,
                              bead_count = NULL) {

  # GATINGSET
  if (!is(x, "GatingSet")) {
    stop("'x' must be a Gatingset object.")
  }

  # CLONE GATINGSET
  gs_clone <- cyto_copy(x)

  # NODES
  nodes <- cyto_nodes(gs_clone, path = "auto")

  # BEADS
  if (is.null(beads)) {
    # SINGLE BEADS
    if (any(grepl("single beads", nodes, ignore.case = TRUE))) {
      beads <- nodes[which(grepl("single beads", nodes, ignore.case = TRUE))[1]]
      # BEADS
    } else if (any(grepl("beads", nodes, ignore.case = TRUE))) {
      beads <- nodes[which(grepl("beads", nodes, ignore.case = TRUE))[1]]
      # BEADS MISSING
    } else {
      stop("Supply the name of the 'beads' population.")
    }
  }

  # BEAD COUNTS
  bead_pops <- cyto_extract(gs_clone, beads)
  bead_counts <- suppressMessages(
    cyto_stats_compute(bead_pops,
      stat = "count",
      format = "wide"
    )
  )
  bead_counts <- bead_counts[, ncol(bead_counts), drop = TRUE]

  # BEADS MISSING
  if (any(bead_counts == 0)) {
    ind <- which(bead_counts == 0)
    stop(paste0(
      "The following samples do not contain any beads:",
      paste0("\n", cyto_names(gs_clone)[ind])
    ))
  }

  # BEAD COUNT
  if (is.null(bead_count)) {
    bead_count <- min(bead_counts)
  } else {
    if (!bead_count <= min(bead_counts)) {
      bead_count <- min(bead_counts)
    }
  }

  # BEAD RATIOS
  bead_ratios <- lapply(seq_len(length(bead_counts)), function(z) {
    1 / (bead_counts[z] / bead_count)
  })

  # SAMPLING - ROOT POPULATION
  pops <- list()
  lapply(seq_along(bead_pops), function(z) {
    pops[[z]] <<- cyto_sample(
      cyto_extract(gs_clone[[z]],
        copy = TRUE
      ),
      bead_ratios[[z]]
    )
  })
  names(pops) <- cyto_names(pops)
  pops <- flowSet_to_cytoset(flowSet(pops))

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
#' @importFrom flowCore fsApply fr_append_cols
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
      suppressWarnings(fr_append_cols(fr, mat))
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
      suppressWarnings(fr_append_cols(fr, mat))
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
#' @param file name of csv file containing columns 'channel' and 'marker'.
#' @param ... additional arguments passed to \code{data_editor}.
#'
#' @return save inputs to "Experiment-Markers.csv" and returns updated samples.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters
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
cyto_markers_edit <- function(x,
                              file = NULL,
                              ...) {

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
        # Channels must match
        if (all(cyto_channels(x) %in% mrks$channel)) {
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
        colnames(dt) <- c("channel", "marker")
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
      colnames(dt) <- c("channel", "marker")
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
      dt <- pd[, c("name", "desc")]
      colnames(dt) <- c("channel", "marker")
      rownames(dt) <- NULL
    }
  }

  # File name not supplied
  if (is.null(file)) {
    file <- paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Markers.csv")
  }

  # Edit dt using data_editor
  dt <- data_editor(dt,
    title = "Experiment Markers Editor",
    save_as = file,
    ...
  )
  
  # Update channels
  BiocGenerics::colnames(x) <- as.character(dt$channel)

  # # TODO - CHECK IF FILE EXISTS WITH DIFFERENT CHANNELS - NEED NEW FILE NAME
  # if(length(grep(file, list.files())) != 0){
  #   # Check if file contains the same
  #
  # }

  # Markers and channels
  cyto_channels <- dt$channel
  cyto_markers <- dt$marker
  names(cyto_markers) <- cyto_channels
  
  # Only modify markers if supplied
  if (!all(is.na(cyto_markers))) {
    if(class(x) %in% c("flowFrame", "flowSet")){
      flowCore::markernames(x) <- cyto_markers
    }else{
      flowWorkspace::markernames(x) <- cyto_markers
    }
  }
  
  # Return updated samples
  return(x)
}

## CYTO_DETAILS_EDIT -----------------------------------------------------------

#' Interactively edit cyto_details for a flowSet or GatingSet
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param file name of csv file containing experimental information.
#' @param ... additional arguments passed to \code{data_editor}.
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
cyto_details_edit <- function(x,
                              file = NULL,
                              ...) {

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

  # FILE
  if (is.null(file)) {
    file <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Details.csv"
    )
  }

  # Edit cyto_details
  pd <- data_editor(pd,
    title = "Experiment Details Editor",
    save_as = file,
    ...
  )
  rownames(pd) <- pd$name

  # Update cyto_details
  cyto_details(x) <- pd

  # Return updated samples
  return(x)
}

## CYTO_DETAILS_SAVE -----------------------------------------------------------

#' Save experiment details to csv file
#'
#' @param x object of class \code{flowSet} or \code{GatingSet} annotated with
#'   experiment details.
#' @param save_as name of csv file to which the experiment details shuld be
#'   saved.
#'
#' @return write experiment details to named csv file.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Modify experiment details manually
#' cyto_details(gs)$Treatment <- c(
#'   rep("A", 8),
#'   rep("B", 8),
#'   rep("C", 8),
#'   rep("D", 8),
#'   NA
#' )
#'
#' # Save experiment details to file
#' cyto_details_save(gs)
#' }
#'
#' @importFrom utils write.csv
#'
#' @export
cyto_details_save <- function(x,
                              save_as = NULL) {

  # SAVE AS
  if (is.null(save_as)) {
    save_as <- paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Details.csv")
  }

  # WRITE CSV FILE
  pd <- cyto_details(x)
  write.csv(pd,
    save_as,
    row.names = FALSE
  )

  return(pd)
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
#' \code{select} argument. Compensation applied to samples can be easily removed
#' by setting the \code{remove} argument to TRUE.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param spillover name of the spillover matrix csv file (e.g.
#'   "Spillover-Matrix.csv") saved in the current working directory or spillover
#'   matrix object of class \code{matrix} or \code{data.frame} with channel
#'   names assigned to each column. If supplied, \code{spillover} will be
#'   applied to all samples.
#' @param select index or name of the sample from which the spillover matrix
#'   should be extracted when no spillover matrix file is supplied to
#'   \code{spillover}. To compensate each sample individually using their stored
#'   spillover matrix file, set \code{select} to NULL.
#' @param remove logical indicating whether applied compensation should be
#'   removed from the supplied samples, set to FALSE by default.
#' @param ... not in use.
#'
#' @return a compensated \code{flowFrame}, \code{flowSet} or \code{GatingSet}
#'   object.
#'
#' @importFrom utils read.csv
#' @importFrom tools file_ext
#' @importFrom methods is
#' @importFrom flowCore decompensate
#' @importFrom flowWorkspace cytoset flowFrame_to_cytoframe gs_cyto_data
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
#' spill <- cyto_spillover_extract(Activation)[[1]]
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
#' # Remove applied compensation
#' cyto_compensate(gs, remove = TRUE)
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
cyto_compensate.GatingSet <- function(x,
                                      spillover = NULL,
                                      select = 1,
                                      remove = FALSE,
                                      ...) {

  # Extract flowSet
  fs <- cyto_extract(x, parent = "root")

  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (is(spillover, "character")) {
      # file extension
      spillover <- file_ext_append(spillover, ".csv")
      # file does not exist
      if(!file.exists(spillover)){
        stop(
          paste(spillover, 
                "does not exist or does not have required permissions.")
          )
      }
      # read in spillover matrix
      spill <- read.csv(spillover, 
                        header = TRUE, 
                        row.names = 1)
      colnames(spill) <- rownames(spill)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% cyto_channels(fs)) |
        !all(rownames(spill) %in% cyto_channels(fs))) {
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
      if (!all(colnames(spillover) %in% cyto_channels(fs))) {
        stop("'spillover' must have valid fluorescent channels as colnames.")
      } else {
        spill <- spillover
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(fs))
      names(spill) <- cyto_names(fs)
      # spillover is a list
    } else {
      spill <- spillover
    }
    # Extract spillover matrix directly from fs
  } else if (is.null(spillover)) {
    if (!is.null(select)) {
      spill <- cyto_spillover_extract(fs[[select]])
      if (is.null(spill)) {
        stop("Unable to extract spillover matrix from selected sample.")
      }
      spill <- rep(spill, length(fs))
      names(spill) <- cyto_names(fs)
    } else {
      spill <- lapply(cyto_names(fs), function(y) {
        sm <- cyto_spillover_extract(fs[[y]])[[1]]
        if (is.null(sm)) {
          stop(paste0(
            "Unable to extract spillover matrix from ",
            cyto_names(fs[[y]]), "."
          ))
        }
        return(sm)
      })
    }
  }

  # Channels
  fluor_channels <- cyto_fluor_channels(fs)

  # Spillover may contain more channels than in samples
  spill <- lapply(spill, function(z) {
    # Select rows - square matrix
    if (nrow(z) == ncol(z)) {
      z <- z[match(fluor_channels, colnames(z)), ]
    }
    # Select columns
    z <- z[, fluor_channels]
    return(z)
  })

  # REMOVE COMPENSATION
  if (remove == TRUE) {
    cf_list <- lapply(seq_along(fs), function(z) {
      cf <- decompensate(fs[[z]], spill[[z]])
      if(!is(cf, "cytoframe")){
        cf <- flowFrame_to_cytoframe(cf)
      }
      return(cf)
    })
    names(cf_list) <- cyto_names(fs)
    cs <- cytoset(cf_list)
    gs_cyto_data(x) <- cs
    # APPLY COMPENSATION
  } else if (remove == FALSE) {
    flowWorkspace::compensate(x, spill)
  }

  # RETURN GATINGSET
  return(x)
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowSet <- function(x,
                                    spillover = NULL,
                                    select = 1,
                                    remove = FALSE,
                                    ...) {

  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (is(spillover, "character")) {
      # file extension
      spillover <- file_ext_append(spillover, ".csv")
      # file does not exist
      if(!file.exists(spillover)){
        stop(
          paste(spillover, 
                "does not exist or does not have the required permissions.")
          )
      }
      # read in spillover matrix
      spill <- read.csv(spillover, 
                        header = TRUE, 
                        row.names = 1)
      colnames(spill) <- rownames(spill)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% cyto_channels(x)) |
        !all(rownames(spill) %in% cyto_channels(x))) {
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
      if (!all(colnames(spillover) %in% cyto_channels(x))) {
        stop("'spillover' must have valid fluorescent channels as colnames.")
      } else {
        spill <- spillover
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(x))
      names(spill) <- cyto_names(x)
      # SPILLOVER IS A LIST
    } else {
      spill <- spillover
    }
    # Extract spillover matrix directly from x
  } else if (is.null(spillover)) {
    if (!is.null(select)) {
      spill <- cyto_spillover_extract(x[[select]])
      if (is.null(spill)) {
        stop("Unable to extract spillover matrix from selected sample.")
      }
      spill <- rep(spill, length(x))
      names(spill) <- cyto_names(x)
    } else {
      spill <- lapply(cyto_names(x), function(y) {
        sm <- cyto_spillover_extract(x[[y]])[[1]]
        if (is.null(sm)) {
          stop(paste0(
            "Unable to extract spillover matrix from ",
            cyto_names(x[[y]]), "."
          ))
        }
        return(sm)
      })
    }
  }

  # Channels
  fluor_channels <- cyto_fluor_channels(x)

  # Spillover may contain more channels than in samples
  spill <- lapply(spill, function(z) {
    # Select rows - square matrix
    if (nrow(z) == ncol(z)) {
      z <- z[match(fluor_channels, colnames(z)), ]
    }
    # Select columns
    z <- z[, fluor_channels]
    return(z)
  })

  # REMOVE COMPENSATION
  if (remove == TRUE) {
    cf_list <- lapply(seq_along(x), function(z) {
      cf <- decompensate(x[[z]], spill[[z]])
      if(!is(cf, "cytoset")){
        cf <- flowFrame_to_cytoframe(cf)
      }
      return(cf)
    })
    names(cf_list) <- cyto_names(x)
    x <- cytoset(cf_list)
    return(x)
    # APPLY COMPENSATION
  } else if (remove == FALSE) {
    flowCore::compensate(x, spill)
  }
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowFrame <- function(x,
                                      spillover = NULL,
                                      select = 1,
                                      remove = FALSE,
                                      ...) {

  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (is(spillover, "character")) {
      # file extension
      spillover <- file_ext_append(spillover, ".csv")
      # file does not exist
      if(!file.exists(spillover)){
        stop(
          paste(spillover,
                "does not exist or does not have the required permissions.")
        )
      }
      # read in spillover matrix
      spill <- read.csv(spillover, 
                        header = TRUE, 
                        row.names = 1)
      colnames(spill) <- rownames(spill)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% cyto_channels(x)) |
        !all(rownames(spill) %in% cyto_channels(x))) {
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
      if (!all(colnames(spillover) %in% cyto_channels(x))) {
        stop("'spillover' must have valid fluorescent channels as colnames.")
      } else {
        spill <- spillover
      }
      # spillover is a list
    } else {
      spill <- spillover[[1]]
    }
    # Extract spillover matrix directly from x
  } else if (is.null(spillover)) {
    spill <- cyto_spillover_extract(x)[[1]]
    if (is.null(spill)) {
      stop(paste0(
        "Unable to extract spillover matrix from",
        cyto_names(x), "."
      ))
    }
  }

  # Channels
  fluor_channels <- cyto_fluor_channels(x)

  # Select rows - square matrix
  if (nrow(spill) == ncol(spill)) {
    spill <- spill[match(fluor_channels, colnames(spill)), ]
  }
  # Select columns
  spill <- spill[, fluor_channels]

  # REMOVE COMPENSATION
  if (remove == TRUE) {
    decompensate(x, spill)
    # APPLY COMPENSATION
  } else if (remove == FALSE) {
    flowCore::compensate(x, spill)
  }
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

## CYTO_NODES_CONVERT ----------------------------------------------------------

#' Convert to unique node paths
#'
#' \code{cyto_nodes_convert} checks whether the supplied nodes exist in the
#' supplied \code{GatingHierarchy} or \code{GatingSet} and returns unique node
#' paths for each of the supplied nodes. In the case of ambiguous nodes,
#' \code{cyto_nodes_convert} will attempt to anchor the node path to a known
#' unique parental node.
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param nodes vectors of nodes to check and convert.
#' @param anchor unique path to parental node to use as an anchor for ambiguous
#'   nodes.
#' @param path specifies whether the returned nodes should be in the "full" or
#'   "auto" format, set to "auto" by default
#'
#' @return vector of unique paths for each of the supplied nodes.
#'
#' @author @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_nodes_convert <- function(x,
                               nodes = NULL,
                               anchor = NULL,
                               path = "auto") {

  # PATHS
  nodes_full <- cyto_nodes(x, path = "full")
  nodes_auto <- cyto_nodes(x, path = "auto")
  nodes_terminal <- basename(nodes_full)

  # STRIP REFERENCE TO ROOT
  nodes <- LAPPLY(nodes, function(node){
    if(grepl("root/", node)){
      node <- gsub("root/", "/", node)
      return(node)
    }
    return(node)
  })
  
  # ANCHOR
  if (!is.null(anchor)) {
    # INVALID ANCHOR
    if (!anchor %in% c(nodes_full, nodes_auto, nodes_terminal)) {
      stop(paste0(anchor, " does not exist in this ", class(x), "!"))
    }
    # TERMINAL ANCHOR
    if (anchor %in% nodes_terminal) {
      anchor_match <- which(LAPPLY(nodes_terminal, function(node_terminal) {
        anchor == node_terminal
      }))
      # AUTO NODE
    } else if (anchor %in% nodes_auto) {
      anchor_match <- which(LAPPLY(nodes_auto, function(node_auto) {
        anchor == node_auto
      }))
      # FULL NODE
    } else if (anchor %in% nodes_full) {
      anchor_match <- which(LAPPLY(nodes_full, function(node_full) {
        anchor == node_full
      }))
    }
    # ANCHOR MUST BE UNIQUE
    if (length(anchor_match) > 1) {
      stop(paste0(anchor, " is not unique within this ", class(x), "!"))
    }
    # AUTO ANCHOR
    anchor <- nodes_auto[anchor_match]
  }

  # CONVERT NODES
  nodes <- LAPPLY(nodes, function(node) {
    # TERMINAL NODE
    if (node %in% nodes_terminal) {
      nodes_match <- which(LAPPLY(nodes_terminal, function(node_terminal) {
        node == node_terminal
      }))
      # AUTO NODE
    } else if (node %in% nodes_auto) {
      nodes_match <- which(LAPPLY(nodes_auto, function(node_auto) {
        node == node_auto
      }))
      # FULL NODE
    } else if (node %in% nodes_full) {
      nodes_match <- which(LAPPLY(nodes_full, function(node_full) {
        node == node_full
      }))
      # NODE DOES NOT EXIST
    } else {
      stop(paste0(node, " does not exist in this ", class(x), "!"))
    }
    # RETURN NODE
    if (length(nodes_match) == 1) {
      nodes <- cyto_nodes(x, path = path)[nodes_match]
      return(nodes)
      # AMBIGUOUS NODE - ANCHOR TO PARENTAL NODE
    } else if (length(nodes_match) > 1) {
      # NO ANCHOR
      if (is.null(anchor)) {
        stop(paste0(node, " is ambiguous in this ", class(x), "!"))
        # ANCHOR
      } else if (!is.null(anchor)) {
        # SPLIT NODE
        node_split <- .cyto_nodes_split(node)
        # SPLIT NODES FULL
        nodes_full_split <- .cyto_nodes_split(nodes_full)
        # SPLIT ANCHOR
        anchor_split <- .cyto_nodes_split(anchor)
        # ROOT ANCHOR
        if (anchor_split == "root") {
          anchor_split <- NULL
        }
        # UNIQUE NODE EXISTS
        ind <- which(LAPPLY(nodes_full_split, function(z) {
          all(unique(c(node_split, anchor_split)) %in% z)
        }))
        # CANNOT FIND UNIQUE NODE PATH
        if (length(ind) == 0) {
          stop(paste0(
            "Failed to generate unique path for ", node,
            " relative to ", anchor, "."
          ))
          # UNIQUE NODE PATH EXISTS
        } else {
          if (length(ind) == 1) {
            node <- cyto_nodes(x, path = path)[ind]
            return(node)
          } else if (length(ind) > 1) {
            # USE SHORTEST PATH (REMOVE DESCENDANTS)
            nodes_unique <- nodes_full_split[ind]
            nodes_lengths <- LAPPLY(nodes_unique, "length")
            nodes_length_min <- min(nodes_lengths)
            if (length(nodes_lengths[nodes_lengths == nodes_length_min]) > 1) {
              stop(paste0(
                "Failed to generate unique path for ", node,
                " relative to ", anchor, "."
              ))
            }
            ind <- which(nodes_lengths == nodes_length_min)
            node <- cyto_nodes(x, path = path)[ind]
            return(node)
          }
        }
      }
    }
  })

  # RETURN UNIQUE NODES
  return(nodes)
}

#' Split node paths into fragments
#' @noRd
.cyto_nodes_split <- function(nodes = NULL) {
  nodes_split <- lapply(nodes, function(node) {
    node <- unlist(strsplit(node, "\\/"))
    node <- node[!LAPPLY(node, ".empty")]
    return(node)
  })
  names(nodes_split) <- nodes

  return(nodes_split)
}

## CYTO_NODES_ANCESTOR ---------------------------------------------------------

#' Get name of most recent common ancestor shared by nodes
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param nodes vector of nodes for which the most recent common ancestor should
#'   be returned.
#' @param ... additional arguments passed to \code{\link{cyto_nodes_convert}} to
#'   control the format of the returned ancestral node.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' library(CytoExploreRData)
#' 
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#' 
#' # Compensation
#' gs <- cyto_compensate(gs)
#' 
#' # Transformations
#' gs <- cyto_transform(gs)
#' 
#' # Gating
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#' 
#' # Ancestral node of CD4 T Cells & CD8 T Cells
#' cyto_nodes_ancestor(gs, nodes = c("CD4 T Cells", "CD8 T Cells"))
#'
#' @export
cyto_nodes_ancestor <- function(x,
                                nodes = NULL, 
                                ...){
  
  # GATINHIERARCHY OR GATINGSET
  if(!is(x, "GatingHierarchy") &
     !is(x, "GatingSet")) {
    stop("'x' must be a GatingHierarchy or GatingSet object.")
  }
  
  # MISSING NODES
  if(is.null(nodes)){
    stop("Supply the name of the nodes to 'nodes'.")
  }
  
  # GET FULL NODES
  nodes_full <- cyto_nodes_convert(x,
                                   nodes = nodes,
                                   path = "full")
  
  # GET SPLIT NODES
  nodes_split <- .cyto_nodes_split(nodes_full)
  
  # GET COMMON ANCESTOR - SHORTEST SHARED PATH WITH FIRST NODE
  ancestor <- c()
  for(i in rev(seq_len(length(nodes_split[[1]])))){
    if(all(LAPPLY(nodes_split[-1], function(z){
      all(nodes_split[[1]][seq_len(i)] %in% z)
    }))){
      ancestor <- c(ancestor, 
                    paste0("/", paste0(nodes_split[[1]][seq_len(i)],
                                       collapse = "/")))
      break()
    }
  }
  
  # NO COMMON ANCESTOR
  if(length(ancestor) == 0){
    ancestor  <- "root"
  }

  # CONVERT ANCESTRAL NODE
  ancestor <- cyto_nodes_convert(x,
                                 nodes = ancestor,
                                 ...)
  
  # RETURN COMMON ANCESTOR
  return(ancestor)
  
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
cyto_copy <- function(x) {

  # GatingSet
  if (is(x, "GatingSet")) {
    x <- gs_clone(x)
  }

  # cytoset
  if (is(x, "cytoset") | is(x, "cytoframe")) {
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
#' @export
cyto_spillover_extract <- function(x) {

  # GATINGSET
  if (is(x, "GatingSet")) {
    spill <- gs_get_compensations(x)
    if (all(LAPPLY(spill, "is.null"))) {
      spill <- NULL
    } else {
      spill <- lapply(spill, function(z) {
        z@spillover
      })
    }
    # GATINGHIERARCHY
  } else if (is(x, "GatingHierarchy")) {
    spill <- gh_get_compensations(x)
    if (!is.null(spill)) {
      spill <- list(spill@spillover)
      names(spill) <- cyto_names(x)
    }
    # FLOWSET
  } else if (is(x, "flowSet")) {
    spill <- lapply(seq_along(x), function(z) {
      # CyTOF lacks spill slot (just in case)
      tryCatch(x[[z]]@description$SPILL, error = function(e) {
        NULL
      })
    })
    names(spill) <- cyto_names(x)
    if (all(LAPPLY(spill, "is.null"))) {
      spill <- NULL
    }
    # FLOWFRAME
  } else if (is(x, "flowFrame")) {
    spill <- tryCatch(x@description$SPILL, error = function(e) {
      NULL
    })
    if (!is.null(spill)) {
      spill <- list(spill)
      names(spill) <- cyto_names(x)
    }
  }

  # RETURN LIST OF SPILLOVER MATRICES
  return(spill)
}

## CYTO_CALIBRATE --------------------------------------------------------------

#' Calibrate channel ranges for cyto_plot
#'
#' The density colour scale for points in \code{cyto_plot} can now be made to
#' represent fluorescent intensity of a third channel, instead of local density
#' in the 2D plotting space. In order for this feature to work properly, the
#' full range for the channel/marker of interest is required to appropriately
#' set the colour scale. \code{cyto_calibrate} performs this range calibration,
#' so that this information can be used in all downstream \code{cyto_plot}
#' calls.
#'
#' @param x object of class \code{cytoframe}, \code{cytoset},
#'   \code{GatingHierarchy} or \code{GatingSet} to use for the calibration. For
#'   the best calibration is recommended that users supply samples containing
#'   both negative and positive events in each channel.
#' @param parent name of the parent population to use for channel calibration
#'   when a \code{GatingHierarchy} or \code{GatingSet} is supplied, set to the
#'   \code{"root"} node by default.
#' @param type indicates the type of calibration to perform, options include
#'   \code{"range"} or \code{"quantile"}. Range calibration simply uses the full
#'   range of values across samples for the calibration. Quantile calibration
#'   computes an lower and upper quantile for each channel, values falling
#'   outside the calibration range are assigned the bottom or top colour.
#' @param probs vector of lower and upper probabilities passed to
#'   \code{stats:quantile} to compute quantiles, set to \code{c(0.01, 0.99)} by
#'   default.
#' @param ... not in use.
#'
#' @return saves calibration settings for use by \code{\link{cyto_plot}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation flowSet
#' fs <- Activation
#'
#' # Calibration
#' cyto_calibrate(fs)
#'
#' # Colour based on Hoechst-405 staining
#' cyto_plot(fs[1],
#' channels = c("FSC-A", "SSC-A"),
#' point_col = "Hoechst-405")
#'
#' @export
cyto_calibrate <- function(x,
                           parent = "root",
                           type = "range",
                           probs = c(0.01, 0.99),
                           ...){
  
  # COMPUTE CHANNEL RANGES
  if(grepl("^r", type, ignore.case = TRUE)){
    cyto_cal <- .cyto_range(x,
                            parent = parent,
                            axes_limits = "data",
                            buffer = 0,
                            ...)
  # COMPUTE CHANNEL QUANTILES
  }else if(grepl("^q", type, ignore.case = TRUE)){
    cyto_cal<- .cyto_quantile(x,
                              parent = parent,
                              probs = probs,
                              ...)
  }
  
  # SAVE RDS TO TEMPFILE
  tempfile <- paste0(tempdir(),
                     .Platform$file.sep,
                     "cyto_calibrate.rds")
  saveRDS(cyto_cal,
          tempfile)
  
}

## CYTO_CALIBRATE_RESET --------------------------------------------------------

#' Remove current calibration settings
#' 
#' @return remove current calibration settings created by \code{cyto_calibrate}.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
cyto_calibrate_reset <- function(){
  
  # TEMP FILES
  temp_dir <- tempdir()
  temp_files <- list.files(temp_dir)
  
  # REMOVE CALIBRATION SETTINGS
  if(any(grepl("cyto_calibrate.rds", temp_files, ignore.case = TRUE))){
    files <- temp_files[grepl("cyto_calibrate.rds", 
                              temp_files,
                              ignore.case = TRUE)]
    file.remove(paste0(temp_dir,
                       .Platform$file.sep,
                       files))
  }
  
}

## CYTO_CALIBRATE_RECALL -------------------------------------------------------

#' Recall saved calibration settings
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
.cyto_calibrate_recall <- function(){
  
  # TEMPFILES
  temp_dir <- tempdir()
  temp_files <- list.files(temp_dir)
  
  # CALIBRATION
  if(any(grepl("cyto_calibrate.rds", temp_files))){
    return(readRDS(paste0(temp_dir,
                          .Platform$file.sep,
                          "cyto_calibrate.rds")))
  }else{
    return(NULL)
  }
  
}

## CYTO_APPLY ------------------------------------------------------------------

#' Apply a function to elements of cytoframe, cytoset, GatingHierarchy or
#' GatingSet
#'
#' \code{cyto_apply} is convenient wrapper around \code{lapply} and \code{apply}
#' to apply a function to a \code{cytoframe}, \code{GatingHierarchy} or over
#' elements of a \code{cytoset}, \code{GatingSet} or list of \code{cytoframes}.
#' \code{cyto_apply} is extremely flexible by supporting functions that accept
#' the data in either \code{cytoframe}, \code{matrix} or \code{vector} formats.
#' All the data processing steps are handled internally prior to passing the
#' data to the specified function. It is important that the arguments in
#' \code{FUN} do not conflict with the arguments of \code{cyto_apply} and they
#' should be supplied to \code{cyto_apply} by name.
#'
#' @param x object of class \code{cytoframe}, \code{cytoset},
#'   \code{GatingHierarchy}, \code{GatingSet} or a list of these objects.
#' @param FUN name of a function to apply to each element of \code{x}.
#' @param ... additional arguments passed to \code{FUN}. Multiple arguments are
#'   supported but must be named, see examples below.
#' @param simplify logical indicating whether attempts should be made to coerce
#'   the output to a \code{cytoset} or \code{matrix}, set to TRUE by default. A
#'   \code{cytoframe} or \code{cytoset} will be returned for
#'   \code{GatingHierarchy} and \code{GatingSet} objects respectively.
#' @param input indicates the data input format as required by \code{FUN} can be
#'   either 1 - "cytoframe", 2 - "matrix", 3 - "column" or 4 - "row", set to
#'   "cytoframe" by default. \code{cyto_apply} will take care of all the data
#'   formatting prior to passing it \code{FUN}. The \code{"column"} and
#'   \code{"row"} options are for functions that expect vectors as the input.
#' @param parent name of the parent population to extract from
#'   \code{GatingHierarchy} or \code{GatingSet} objects, set to \code{"root"} by
#'   default.
#' @param copy logical indicating whether the indicated function should be
#'   applied to a copy of each \code{cytoframe}, set to TRUE by default. Apply
#'   the function to a copy of each \code{cytoframe} ensures that the underlying
#'   data remains unchanged.
#' @param channels vector of channels which should be included in the data
#'   passed to \code{FUN}, set to all channels by default.
#' @param trans object of class \code{transformerList} containing the
#'   definitions of the transformers applied to the supplied data. These
#'   transformers are required to inverse transformations when \code{inverse =
#'   TRUE}.
#' @param inverse logical indicating whether each \code{cytoframe} should be
#'   inverse transformed prior to applying \code{FUN}.
#'
#' @importFrom methods is
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace cytoset
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#' library(flowWorkspace)
#'
#' # Activation cytoset
#' fs <- Activation
#' cs <- flowSet_to_cytoset(fs)
#'
#' # Compute quantiles
#' cyto_apply(cs,
#' cyto_quantile,
#' probs = c(0.05, 0.95),
#' na.rm = TRUE,
#' input = 3)
#'
#' @name cyto_apply
NULL

#' @noRd
#' @export
cyto_apply <- function(x, 
                       ...){
  UseMethod("cyto_apply")
}

#' @rdname cyto_apply
#' @export
cyto_apply.list <- function(x,
                            FUN,
                            ...,
                            simplify = TRUE,
                            input = "cytoframe",
                            parent = NULL,
                            copy = TRUE,
                            channels = NULL,
                            trans = NA,
                            inverse = FALSE){
  
  # FUNCTION MISSING
  if(missing(FUN)) {
    stop("Supply a function to 'FUN'.")
  }
  
  # NAMED LIST
  if(is.null(names(x))){
    names(x) <- LAPPLY(x, cyto_names)
  }
  
  # APPLY FUNCTION
  res <- lapply(x, function(z) {
    cyto_apply(z,
               FUN = FUN,
               ...,
               simplify = simplify,
               input = input,
               parent = parent,
               copy = copy,
               channels = channels,
               trans = trans,
               inverse = inverse)
  })
  names(res) <- names(x)
  
  # SIMPLIFY
  if(simplify == TRUE){
    # DON'T COERCE TO CYTOSET
    if(!all(LAPPLY(res, is, "flowFrame"))){
      res <- do.call("rbind", res)
    }
  }
  
  # RETURN OUTPUT
  return(res)
  
}

#' @rdname cyto_apply
#' @export
cyto_apply.GatingSet <- function(x,
                                 FUN,
                                 ..., 
                                 simplify = TRUE,
                                 input = "cytoframe", 
                                 parent = NULL,
                                 copy = TRUE,
                                 channels = NULL,
                                 trans = NA,
                                 inverse = FALSE) {
  
  # FUNCTION MISSING
  if(missing(FUN)) {
    stop("Supply a function to 'FUN'.")
  }
  
  # EXTRACT POPULATION
  cs <- cyto_extract(x,
                     parent = parent)
  
  # EXTRACT TRANSFORMERS
  trans <- cyto_transformer_extract(x)
  
  # CYTOSET METHOD
  res <- cyto_apply(cs,
                    FUN = FUN,
                    simplify = simplify,
                    input = input,
                    copy = copy,
                    channels = channels,
                    trans = trans,
                    inverse = inverse,
                    ...)
  
  # RETURN OUTPUT
  return(res)
  
}

#' @rdname cyto_apply
#' @export
cyto_apply.GatingHierarchy <- function(x,
                                       FUN,
                                       ..., 
                                       simplify = TRUE,
                                       input = "cytoframe", 
                                       parent = NULL,
                                       copy = TRUE,
                                       channels = NULL,
                                       trans = NA,
                                       inverse = FALSE) {
  
  # FUNCTION MISSING
  if(missing(FUN)) {
    stop("Supply a function to 'FUN'.")
  }
  
  # EXTRACT POPULATION
  cf <- cyto_extract(x,
                     parent = parent)
  
  # EXTRACT TRANSFORMERS
  trans <- cyto_transformer_extract(x)
  
  # CYTOSET METHOD
  res <- cyto_apply(cf,
                    FUN = FUN,
                    simplify = simplify,
                    input = input,
                    copy = copy,
                    channels = channels,
                    trans = trans,
                    inverse = inverse,
                    ...)
  
  # RETURN OUTPUT
  return(res)
  
}

#' @rdname cyto_apply
#' @export
cyto_apply.flowSet <- function(x,
                               FUN,
                               ..., 
                               simplify = TRUE,
                               input = "cytoframe", 
                               copy = TRUE,
                               channels = NULL,
                               trans = NA,
                               inverse = FALSE){
  
  # FUNCTION MISSING
  if(missing(FUN)) {
    stop("Supply a function to 'FUN'.")
  }
  
  # NAMES
  nms <- cyto_names(x)
  
  # APPLY FUNCTION
  res <- structure(lapply(nms, function(z) {
    # CYTOFRAME
    y <- x[[z]]
    # APPLY FUNCTION
    cyto_apply(y,
               FUN = FUN,
               ...,
               simplify = simplify,
               input = input,
               copy = copy,
               channels = channels,
               trans = trans,
               inverse = inverse)
  }), names = nms)
  
  # SIMPLIFY
  if(simplify){
    # CYTOFRAMES TO CYTOSET
    if(all(LAPPLY(res, is, "flowFrame"))){
      if(is(res[[1]], "cytoframe")){
        res <- cytoset(res)
      }else{
        res <- flowSet(res)
      }
      # LIST OF MATRICES TO MATRIX
    }else{
      res <- do.call("rbind", res)
    }
  }
  
  # RETURN OUTPUT
  return(res)
  
}

#' @rdname cyto_apply
#' @export
cyto_apply.flowFrame <- function(x,
                                 FUN,
                                 ...,
                                 simplify = TRUE,
                                 input = "cytoframe",
                                 copy = TRUE,
                                 channels = NULL,
                                 trans = NA,
                                 inverse = FALSE){
  
  # FUNCTION MISSING
  if(missing(FUN)) {
    stop("Supply a function to 'FUN'.")
  }
  
  # FUNCTION NAME
  if(is.function(FUN)) {
    FUN_NAME <- deparse(substitute(FUN))
    # UNNAMED FUNCTION
    if(length(FUN_NAME) > 1) {
      FUN_NAME <- "FUN"
    }
  }else {
    FUN_NAME <- FUN
  }
  
  # PREPARE FUNCTION
  FUN <- match.fun(FUN)
  
  # CYTOFRAME INPUT
  if(input == 1 | grepl("^cy", input) | grepl("^cf", input)) {
    input <- "cytoframe"
    # MATRIX INPUT
  } else if(input == 2 | grepl("^m", input)) {
    input <- "matrix"
    # COLUMN/CHANNEL INPUT
  } else if(input == 3 | grepl("^co", input) | grepl("^ch", input)) {
    input <- "column"
    # ROW/CELL INPUT
  } else if(input == 4 | grepl("^r", input) | grepl("^ce", input)) {
    input <- "row"
  }
  
  # COPY
  if(copy){
    cf <- cyto_copy(x)
  }else{
    cf <- x
  }
  
  # CHANNELS
  if(!is.null(channels)){
    cf <- cf[, cyto_channels_extract(cf, channels)]
  }
  
  # INVERSE
  if(inverse & !.all_na(trans)){
    cf <- cyto_transform(cf,
                         trans = trans,
                         inverse = TRUE,
                         plot = FALSE)
  }
  
  # CYTOFRAME INPUT
  if(input == "cytoframe"){
    res <- FUN(cf, ...)
    # MATRIX/VECTOR INPUTS
  }else{
    # MATRIX
    cf_raw <- cyto_extract(cf,
                           raw = TRUE)[[1]]
    # ENTIRE MATRIX
    if(input == "matrix"){
      res <- FUN(cf_raw, ...)
      # COLUMN-WISE
    }else if(input == "column"){
      res <- apply(cf_raw, 2, FUN, ...)
      # ROW-WISE
    }else if(input == "row"){
      res <- apply(cf_raw, 1, FUN, ...)
    }
  }
  
  # VECTOR OUTPUT - INHERIT FUNCTION NAME
  if(is.null(dim(res))){
    if(is.null(names(res))){
      if(length(res) == 1){
        names(res) <- FUN_NAME
      }else{
        names(res) <- paste0(FUN_NAME, "-", seq_len(length(res)))
      }
    }
  }
  
  # SIMPLIFY
  if(simplify == TRUE){
    # VECTOR OUTPUT TO MATRIX
    if(is.null(dim(res))){
      res <- matrix(res,
                    ncol = length(res),
                    dimnames = list(cyto_names(x),
                                    names(res)))
    }
    # ROW NAMES
    if(any(c("matrix", "data.frame") %in% is(res))) {
      if(is.null(rownames(res))){
        if(nrow(res) > 1){
          rownames(res) <- paste(cyto_names(x), 1:nrow(res), sep = "|")
        }else{
          rownames(res) <- cyto_names(x)
        }
      }else{
        if(!all(rownames(res) == cyto_names(x))) {
          rownames(res) <- paste(cyto_names(x), rownames(res), sep = "|")
        }
      }
    }
  }
  
  # RETURN OUTPUT
  return(res)
  
}

## CYTO_CBIND ------------------------------------------------------------------

#' Bind new columns to cytoframe or cytoset
#'
#' @param x object of class \code{cytoframe}, \code{cytoset} or a list of
#'   \code{cytoframes}.
#' @param cols matrix of columns to be added to \code{x} can be supplied in a
#'   list for \code{cytoset} and \code{list} methods.
#'
#' @importFrom flowWorkspace cf_append_cols
#' @importFrom flowCore fr_append_cols
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_cbind
NULL

#' @noRd
#' @export
cyto_cbind <- function(x,
                       ...){
  UseMethod("cyto_cbind")
}

#' @rdname cyto_cbind
#' @export
cyto_cbind.flowFrame <- function(x,
                                 cols = NULL){
  
  # MATRIX
  if(!is.matrix(cols)){
    stop("cols must be a matrix")
  }
  
  # CYTOFRAME
  if(is(x, "cytoframe")){
    cf_append_cols(x, cols)
    # FLOWFRAME  
  }else{
    fr_append_cols(x, cols)
  }
  
}

#' @rdname cyto_cbind
#' @export
cyto_cbind.flowSet <- function(x, 
                               cols = NULL){
  
  # MATRIX
  if(is.matrix(cols)){
    # COUNTS
    cyto_counts <- cyto_apply(x, 
                              "nrow", 
                              input = 2,
                              copy = FALSE)
    # SAME NUMBER OF EVENTS
    if(nrow(cols) != sum(cyto_counts)){
      stop(
        paste0("'cols' does not contain the same number of events as this", 
               class(x))
      )
      # SPLIT MATRIX INTO LIST
    }else{
      cols <- lapply(seq_along(cyto_counts), function(z){
        if(z == 1){
          cols[1:cyto_counts[z], ]
        }else{
          start <- sum(cyto_counts[1:(z-1)]) + 1
          end <- start + cyto_counts[z] - 1
          cols[start:end, ]
        }
      })
      names(cols) <- cyto_names(x)
    }
  }
  
  # BIND COLUMNS
  cf_list <- lapply(cyto_names(x), function(z){
    cyto_cbind(x[[z]], cols[[z]])
  })
  names(cf_list) <- cyto_names(x)
  
  # RETURN FLOWSET/CYTOSET
  if(is(cf_list[[1]], "cytoframe")){
    return(cytoset(cf_list))
  }else{
    return(flowSet(cf_list))
  }
  
}

#' @noRd
#' @export
cyto_cbind.list <- function(x,
                            cols = NULL){
  
  # MATRIX
  if(is.matrix(cols)){
    # COUNTS
    cyto_counts <- cyto_apply(x, 
                              "nrow", 
                              input = 2,
                              copy = FALSE)
    # SAME NUMBER OF EVENTS
    if(nrow(cols) != sum(cyto_counts)){
      stop(
        paste0("'cols' does not contain the same number of events as this", 
               class(x))
      )
      # SPLIT MATRIX INTO LIST
    }else{
      cols <- lapply(seq_along(cyto_counts), function(z){
        if(z == 1){
          cols[1:cyto_counts[z], ]
        }else{
          start <- sum(cyto_counts[1:(z-1)]) + 1
          end <- start + cyto_counts[z] - 1
          cols[start:end, ]
        }
      })
      names(cols) <- cyto_names(x)
    }
  }
  
  # BIND COLUMNS
  cf_list <- lapply(seq_along(x), function(z){
    cyto_cbind(x[[z]], cols[[z]])
  })
  names(cf_list) <- cyto_apply(x, 
                               "cyto_names", 
                               input = "cytoframe",
                               copy = FALSE)
  
  # RETURN LIST OF APPENDED CYTOFRAMES
  return(cf_list)
  
}

## CYTO_LIST -------------------------------------------------------------------

#' Convert cytoframe or cytoset to named cytoframe list
#' 
#' @param x object of class \code{cytoframe} or \code{cytoset} to be listed.
#' 
#' @return list of \code{cytoframe} objects.
#' 
#' @importFrom methods is
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @examples 
#' library(CytoExploreRData)
#' 
#' # Activation dataset
#' fs <- Activation
#'
#' # Convert to list
#' cf_list <- cyto_list(fs)
#' 
#' @export
cyto_list <- function(x){
  if(is(x, "flowFrame")) {
    structure(list(x), names = cyto_names(x))
  } else if(is(x, "flowSet")) {
    cyto_apply(x, function(z){return(z)}, simplify = FALSE)
  }
}
