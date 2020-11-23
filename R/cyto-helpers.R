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
      # APPEND FILE EXTENSION
      gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
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
      save_as <- file_ext_append(save_as, ".wsp")
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
#' @param ... additional arguments passed to \code{\link{cyto_load}}.
#'
#' @return object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @importFrom flowCore identifier identifier<- parameters
#' @importFrom flowWorkspace load_gs load_cytoset_from_fcs pData
#'
#' @examples
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
  
  # DIRECTORY NOT FOUND
  if (!dir.exists(path)) {
    stop("Specified path does not exist.")
  }
  
  # FILE PATHS
  files <- list.files(path, full.names = TRUE)
  
  # REMOVE DIRECTORIES
  files <- files[!dir.exists(files)]
  
  # SAVED GATINGSET
  if ("pb" %in% file_ext(files)) {
    # LOAD GATINGSET
    x <- load_gs(path = path,
                 backend_readonly = FALSE)
    # FCS FILES
  } else {
    
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
      files <- file_sort(files)
    }
    
    # CYTOSET
    x <- load_cytoset_from_fcs(files = normalizePath(files), ...)
    
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
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{CatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}. The \code{root}
#'   node extracted when a \code{GatingSet} or \code{GatingHierachy} is
#'   supplied.
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
#'   \url{https://academic.oup.com/bioinformatics/article/32/16/2473/2240408}
#'   
#' @seealso \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}
#'
#' @export
cyto_clean <- function(x, ...) {
  
  # GATINGSET/GATINGHIERARCHY
  if (cyto_class(x, c("GatingHierarchy", "GatingSet"), TRUE)) {
    # PARENT
    parent <- cyto_nodes(x, path = "auto")[1]
    # EXTRACT DATA
    cyto_data <- cyto_extract(x, parent)
    # flowAI REQUIRES FLOWSET
    if (cyto_class(cyto_data, "cytoset", TRUE)) {
      cyto_data <- cytoset_to_flowSet(cyto_data) # REMOVE
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
    if (cyto_class(cyto_data, "flowSet", TRUE)) {
      cyto_data <- flowSet_to_cytoset(cyto_data)
    }
    # REPLACE DATA
    gs_cyto_data(x) <- cyto_data
  } else {
    # FLOWSET REQUIRED
    if (cyto_class(x, "cytoset", TRUE)) {
      x <- cytoset_to_flowSet(x)
    }
    # CLEAN DATA
    invisible(capture.output(x <- flow_auto_qc(x,
                                               html_report = FALSE,
                                               mini_report = FALSE,
                                               fcs_QC = FALSE,
                                               folder_results = FALSE,
                                               ...
    )))
    # RETURN CYTOSET
    if (cyto_class(x, "flowSet", TRUE)) {
      x <- flowSet_to_cytoset(x)
    }
  }
  return(x)
}

## CYTO_SETUP ------------------------------------------------------------------

#' Load, annotate and prepare FCS files into a GatingSet
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
#' @param ... additional arguments passed to \code{\link{cyto_load}}.
#'
#' @return object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @importFrom flowWorkspace GatingSet
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
  if (cyto_class(x, "flowSet")) {
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
                     "nrow",
                     input = 2,
                     copy = FALSE)
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
  
  # gatingTemplate
  if (!is.null(gatingTemplate)) {
    
    # FILE EXTENSION APPEND
    gatingTemplate <- file_ext_append(gatingTemplate, ".csv")
    
    # ACTIVE GATINGTEMPLATE
    message(paste("Setting", gatingTemplate, "as the active gatingTemplate..."))
    cyto_gatingTemplate_select(gatingTemplate)
    
    # CREATE GATINGTEMPLATE
    if(!file_exists(gatingTemplate)){
      message(paste("Creating", gatingTemplate, "..."))
      cyto_gatingTemplate_create(gatingTemplate)
    }
  }
  
  # SETUP COMPLETE
  message("Done!")
  
  # RETURN GATINGSET
  return(x)
  
}

## CYTO_CLASS ------------------------------------------------------------------

#' Get class of a cytometry object
#'
#' @param x either a \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} for which the class
#'   should be checked or returned.
#' @param expect a class of object supplied as a character string, if supplied a
#'   logical value will be returned indicating whether \code{x} is of the
#'   expected class.
#' @param class indicates whether the class of \code{x} should be matched based
#'   on the objects class (exact class match) or based on its reference class.
#'   This argument is set to TRUE by default if no class expectation is
#'   supplied, otherwise it is set to FALSE if not specified manually. See
#'   examples for more details.
#'
#' @return class of \code{x} or a logical indicating whether \code{x} matches
#'   the expected class.
#'
#' @importFrom methods is
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # flowSet
#' fs <- Activation
#' cyto_class(fs)
#' cyto_class(fs, expect = "flowSet")
#'
#' # GatingSet
#' gs <- GatingSet(fs)
#' cyto_class(gs)
#' cyto_class(gs, expect = "GatingSet")
#'
#' # cytoset
#' cs <- cyto_extract(gs, "root")
#' cyto_class(cs)
#' cyto_class(cs, expect = "flowSet", class = FALSE) # reference class
#'
#' # cytoframe
#' cf <- cs[[1]]
#' cyto_class(cf)
#' cyto_class(cf, expect = "flowFrame", class = FALSE) # reference class
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_class <- function(x,
                       expect = NULL,
                       class) {
  
  # EXACT CLASS
  if(missing(class)) {
    if(is.null(expect)) {
      class <- TRUE
    } else {
      class <- FALSE
    }
  }
  
  # CLASS
  if(is.null(expect)) {
    if(class){
      return(class(x))
    } else {
      return(is(x))
    }
    # CLASS CHECK
  } else {
    if(class) {
      return(any(class(x) %in% expect))
    } else {
      return(any(is(x) %in% expect))
    }
  }
  
}

## CYTO_DETAILS ----------------------------------------------------------------

#' Extract experiment details
#'
#' Simply an autocomplete-friendly wrapper around
#' \code{\link[flowWorkspace:pData-methods]{pData}} with some additional
#' formatting options.
#'
#' @param x object of class \code{\link[flowWorkspace-cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param convert logical indicating whether
#'   \code{\link[utils:type.convert]{type.convert}} should be called on each
#'   column.
#' @param drop logical indicating whether row names should be dropped, set to
#'   FALSE by default.
#' @param factor logical indicating whether characters should be converted to
#'   factors when convert is TRUE, set to FALSE by default.
#' @param ... additional arguments passed to
#'   \code{\link[utils:type.convert]{type.convert}}.
#'
#' @return experiment details as data.frame.
#'
#' @importFrom flowWorkspace pData
#' @importFrom utils type.convert
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'                           
#' # Experiment variables
#' cyto_details(gs)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_details <- function(x,
                         convert = FALSE,
                         drop = FALSE,
                         factor = FALSE, 
                         ...) {
  
  # CYTOFRAME
  if(cyto_class(x, "flowFrame")) {
    stop(
      paste0("'cyto_details' cannot be extracted from ",
             cyto_class(x, class = TRUE), " objects!")
    )
  }
  
  # EXPERIMENT VARIABLES
  pd <- pData(x)
  # ROWNAMES
  if(drop) {
    rownames(pd) <- NULL
  }
  # DATA.FRAME
  if(convert){
    pd <- type.convert(pd, 
                       as.is = !factor,
                       ...)
  }
  
  # EXPERIMENT DETAILS
  return(pd)
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
#' \code{\link[flowWorkspace:sampleNames]{sampleNames}} to extract the sample
#' names from a cytoset, GatingHierarchy or GatingSet.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return names associated with the supplied object.
#'
#' @importFrom flowWorkspace sampleNames
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                          package = "CytoExploreRData"))
#'
#' # GatingSet
#' cyto_names(gs)
#'
#' # GatingHierarchy
#' cyto_names(gs[[1]])
#'
#' # cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#' cyto_names(cs)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_names
#'
#' @export
cyto_names <- function(x) {
  
  if(cyto_class(x, "list", TRUE)) {
    LAPPLY(x, function(z){
      if(cyto_class(x, "flowFrame")) {
        stop(
          paste0("'cyto_names' cannot be extracted from ", 
                 cyto_class(z, class = TRUE), " objects!")
        )
      }
      sampleNames(z)
    })
  }else {
    if(cyto_class(x, "flowFrame")) {
      stop(
        paste0("'cyto_names' cannot be extracted from ", 
               cyto_class(x, class = TRUE), " objects!")
      )
    }
    sampleNames(x)
  }
  
}

# CYTO_NAMES REPLACEMENT METHOD ------------------------------------------------

#' Replacement method for cyto_names
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param value vector of replacement names.
#'
#' @importFrom flowWorkspace sampleNames<-
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                          package = "CytoExploreRData"))
#'
#' # GatingSet
#' cyto_names(gs)[1] <- "Activation_001.fcs"
#' cyto_names(gs)
#'
#' # cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#' cyto_names(cs)[1:2] <- c("Activation_1.fcs", "Activation_2.fcs")
#' cyto_names(cs)
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
"cyto_names<-" <- function(x, value) {
  if (cyto_class(x, c("flowSet", "GatingHierarchy", "GatingSet"))) {
    sampleNames(x) <- value
  } else {
    stop(
      paste0("'cyto_names' cannot be replaced for objects of class ", 
             cyto_class(x, class = TRUE), "!")
    )
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
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'  \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param vars vector containing the names of the variables to be added to
#'   \code{cyto_details}, set to var_1, var_2 ... by default.
#' @param split delimiter to split the file name into fragments, set to "_" by
#'   default.
#' @param exclude vector of indices indicating which text fragments should not
#'   be included in the experiment details, set to NULL by default.
#' @param ... additional arguments passed to
#'   \code{\link[base:strsplit]{strsplit}}.
#'
#' @return cytoset or GatingSet with \code{cyto_details} updated with new
#'   variables extracted from the file names.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'                           
#' # Parse file names to variables
#' gs <- cyto_names_parse(gs,
#' vars = c("sample_type", "sample_id"),
#' split = "_")
#'
#' # Updated experiment details
#' cyto_details(gs)
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
  cyto_names <- file_ext_remove(cyto_names)
  
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

## CYTO_TRANSFORM --------------------------------------------------------------

#' Apply Transformations to Cytometry Data
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
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
#' @importFrom flowCore transform transformList
#' @importFrom grDevices n2mfrow
#' @importFrom graphics par
#' @importFrom methods as
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
    if (cyto_class(x, c("flowFrame", "flowSet"))) {
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
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    
    # Extract transformations from transformerList to transformList
    transform_list <- cyto_transform_extract(transformer_list,
                                             inverse = inverse
    )
    
    # Apply transformations
    x <- suppressMessages(transform(x, transform_list))
    
    # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (cyto_class(x, "GatingSet")) {
    
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
  if (cyto_class(x, "GatingSet")) {
    stop(paste(
      "GatingHierarchy and GatingSet objects require transformerList",
      "objects to apply transformations."
    ))
  }
  
  # TRANSFORM FLOWFRAME OR FLOWSET
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    
    # Restrict (subsetted data may lack channels)
    trans <- trans@transforms[names(trans) %in% cyto_channels(x)]
    trans <- transformList(names(trans),
                           lapply(trans, function(z){
                             z@f
                           }))
    
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
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    
    # Restrict transformers
    trans <- cyto_transformer_combine(trans[names(trans) %in% cyto_channels(x)])
    
    # Extract transformations to transformList
    transform_list <- cyto_transform_extract(trans, inverse = inverse)
    
    # Apply transformations
    x <- suppressMessages(transform(x, transform_list))
    
    
    # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (cyto_class(x, "GatingSet")) {
    
    # Inverse transformations not yet supported
    if (inverse == TRUE) {
      stop(paste(
        "Inverse transformations are not yet supported for",
        "GatingHierarchy/GatingSet objects."
      ))
    }
    
    # Restrict transformers - subsetted data may lack channels
    trans <- cyto_transformer_combine(trans[names(trans) %in% cyto_channels(x)])
    
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

## .CYTO_TRANSFORM -------------------------------------------------------------

#' Transform vector of values given transformers
#' 
#' @param x values to transform.
#' @param trans transformerList containing transformation definitions
#' @param channel which transformer to apply.
#' @param inverse logical indicating whether inverse transformation should be
#'   applied.
#'   
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'   
#' @noRd
.cyto_transform <- function(x,
                            trans = NULL,
                            channel = NULL,
                            inverse = FALSE,
                            ...){
  
  # NO TRANSFORMERS
  if(is.null(trans)){
    return(x)
    # NO TRANSFORMERS
  }else if(.all_na(trans)){
    return(x)
    # TRANSFORMER TO USE
  }else{
    if(is.null(channel)){
      stop("Indicate which transformer to apply through 'channel' argument.")
    }
  }
  
  # TRANSFORMER
  if(channel %in% names(trans)){
    transformer <- trans[[channel]]
  }else{
    return(x)
  }
  
  # SKIP NA
  ind <- which(!is.na(x))
  
  # TRANSFORM
  if(inverse == FALSE){
    x[ind] <- transformer$transform(x[ind])
    # INVERSE TRANSFORM  
  }else{
    x[ind] <- transformer$inverse(x[ind])
  }
  
  # RETURN TRANSFORMED VALUES
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
#'
#' @examples
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
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_transform_extract <- function(x,
                                   inverse = FALSE) {
  
  # TransformLists are returned unaltered
  if (cyto_class(x, "transformList")) {
    return(x)
    # TransformList extracted from transformerList
  } else if (cyto_class(x, "transformerList")) {
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

## CYTO_DATA_EXTRACT -----------------------------------------------------------

#' Extract data from cytosets, GatingHierarchies or GatingSets
#'
#' \code{cyto_data_extract} will extract populations specified by \code{parent}
#' and return a list of cytosets. \code{cyto_data_extract} can also optionally
#' perform additional tasks on the extracted cytosets, such subsetting channels,
#' inversing transformations and returning raw data matrices.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to extract from
#'   \code{GatingHierarchy} or \code{GatingSet} objects.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}}.
#' @param copy logical indicating whether a deep copy of the extracted data
#'   should be returned, set to TRUE by default to ensure that any changes to
#'   the data will not affect the supplied data.
#' @param format can be either \code{"cytoframe"}, \code{"cytoset"} or
#'   \code{"matrix"} to indicate in what format the data should be returned in
#'   the lists, set to \code{"cytoset"} by default.
#' @param channels names of the markers or channels for which data should be
#'   extracted, set to all channels by default.
#' @param markers logical indicating whether the column names of the raw data
#'   matrices should be replaced with the marker names where possible, set to
#'   FALSE by default.
#' @param split logical indicating whether cytosets should be split into
#'   cytosets containing a single sample each, set to FALSE by default.
#' @param trans object of class transformerList required to inverse
#'   transformations applied to the data when \code{inverse = TRUE}.
#' @param inverse logical indicating whether transformations applied to the data
#'   should be inversed, set to FALSE by default.
#'
#' @return either a list of cytosets, list of cytoframe lists or a list of raw
#'   data matrix lists.
#'
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#'
#' @examples
#' # Load in CytoExploreRData to access data
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(system.file("extdata/Activation-GatingSet",
#'                 package = "CytoExploreRData"))
#'
#' # Extract cytosets & inverse transform
#' cyto_data_extract(gs[[1]],
#' parent = c("CD4 T Cells", "CD8 T Cells"),
#' inverse = TRUE)
#'
#' # Extract raw data matrices
#' cyto_data_extract(gs,
#' parent = "root",
#' format = "matrix",
#' channels = c("CD4", "CD8"))
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_data_extract <- function(x,
                              parent = "root",
                              select = NULL,
                              copy = TRUE,
                              format = "cytoset",
                              channels = NULL,
                              markers = FALSE,
                              split = FALSE,
                              trans = NA,
                              inverse = FALSE) {
  
  # TRANSFORMERS
  trans <- cyto_transformer_extract(x)
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(x, select)
  }
  
  # CYTOSET
  if(cyto_class(x, "flowSet")) {
    cs_list <- list(x)
    # GATINGHIERARCHY/GATINGSET
  } else {
    cs_list <- lapply(parent, function(z){
      gs_pop_get_data(x, z)
    })
    names(cs_list) <- parent
  }
  
  # PREPARE CYTOSET LIST
  res <- lapply(cs_list, function(cs) {
    # COPY
    if(copy) {
      cs <- cyto_copy(cs)
    }
    # RESTRICT
    if(!is.null(channels)) {
      channels <- cyto_channels_extract(cs, channels)
      cs <- cs[, channels, drop = FALSE]
    }
    # TRANSFORM
    if(inverse & !.all_na(trans)) {
      cs <- cyto_transform(cs,
                           trans = trans,
                           inverse = inverse,
                           plot = FALSE)
    }
    # CYTOFRAME
    if(format == "cytoframe") {
      structure(
        lapply(cyto_names(cs), function(z){
          cs[[z]]
        }), names = cyto_names(cs)
      )
      # CYTOSET
    } else if(format == "cytoset") {
      if(split) {
        structure(lapply(seq_along(cs), function(z){
          cs[z]
        }), names = cyto_names(cs))
      } else {
        cs
      }
      # MATRIX
    } else if(format == "matrix") {
      structure(
        lapply(cyto_names(cs), function(z){
          m <- exprs(cs[[z]])
          # MARKERS
          if(markers){
            colnames(m) <- unname(cyto_markers_extract(cs[[z]], colnames(m)))
          }
          return(m)
        }), names = cyto_names(cs)
      )
    }
  })
  names(res) <- names(cs_list)
  
  # EXTRACTED DATA
  return(res)
}

#' @export
cyto_extract <- function(...){
  .Defunct("cyto_data_extract")
}

## CYTO_FILTER -----------------------------------------------------------------

#' Filter samples based on experiment variables
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... tidyverse-style subsetting using comma separated logical
#'   predicates based on experimental variables stored in
#'   \code{cyto_details(x)}. See examples below for demonstration.
#'
#' @return \code{cytoset} or \code{GatingSet} restricted to samples which meet
#'   the filtering criteria.
#'
#' @importFrom dplyr filter
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
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
  if (!cyto_class(x, "flowSet") & !cyto_class(x, "GatingSet", TRUE)) {
    stop("'x' should be an object of class cytoset or GatingSet.")
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
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... named list containing experimental variables to be used to select
#'   samples or named arguments containing the levels of the variables to
#'   select. See below examples for use cases. Selected samples can be excluded
#'   by setting \code{exclude} to TRUE.
#'
#' @return \code{cytoset} or \code{GatingSet} restricted to samples which meet
#'   the designated selection criteria.
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Look at experiment details
#' cyto_details(gs)
#'
#' # Select Stim-C samples with 0 and 500 nM OVA concentrations
#' gs_stim_c <- cyto_select(gs,
#'   Treatment = "Stim-C",
#'   OVAConc = c(0, 500)
#' )
#'
#' # Select Stim-A and Stim-C treatment groups
#' gs_stim_ac <- cyto_select(
#'   gs,
#'   list("Treatment" = c("Stim-A", "Stim-C"))
#' )
#'
#' # Exclude Stim-D treatment group
#' gs_stim_abc <- cyto_select(gs,
#'   Treatment = "Stim-D",
#'   exclude = TRUE
#' )
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_select <- function(x, ...) {
  
  # Check class of x
  if (!cyto_class(x, "flowSet") & !cyto_class(x, "GatingSet", TRUE)) {
    stop("'x' should be an object of class cytoset or GatingSet.")
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

## CYTO_GROUPS -----------------------------------------------------------------

#' Retrieve details about experimental groups
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param group_by names of cyto_details variables to use for grouping For more
#'   control over the group order, specify the factor levels for each variable
#'   in a list (e.g. list(Treatment = c("Stim-A","Stim-C","Stim-B", "Stim-D"))).
#' @param details logical indicating whether the split experimental details
#'   should be returned instead of the group names, set to FALSE by default.
#'
#' @return names of experimental groups or a list of experiment details per
#'   experimental group.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Treatment & OVAConc group names
#' cyto_groups(gs,
#' group_by = list(Treatment = c("Stim-A", "Stim-B", "Stim-C", "Stim-D"),
#' OVAConc = c(0, 5, 50, 500)))
#'
#' # Treatment & OVAConc group details
#' cyto_groups(gs,
#' group_by = list(Treatment = c("Stim-A", "Stim-B", "Stim-C", "Stim-D"),
#' OVAConc = c(0, 5, 50, 500)),
#' details = TRUE)
#'
#' @export
cyto_groups <- function(x, 
                        group_by = "all",
                        details = FALSE){
  
  # Check class of x
  if (!cyto_class(x, "flowSet") & !cyto_class(x, "GatingSet", TRUE)) {
    stop("'x' should be an object of class cytoset or GatingSet.")
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
  
  # RETURN SPLIT DETAILS
  if(details == TRUE){
    return(pd_split)
    # RETURN GROUP NAMES
  }else{
    return(names(pd_split))
  }
  
}

## CYTO_SORT_BY ----------------------------------------------------------------

#' Sort a flowSet or GatingSet by experiment variables
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param sort_by names of cyto_details variables to use for sorting. For more
#'   control over the sorting order, specify the factor levels for each variable
#'   in a list (e.g. list(Treatment = c("Stim-A","Stim-C","Stim-B", "Stim-D"))).
#'
#' @return \code{cytoset} or \code{GatingSet} object sorted based on experiment
#'   variables.
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Sort by Treatment and OVAConc
#' gs <- cyto_sort_by(gs,
#' list(Treatment = c("Stim-A", "Stim-B", "Stim-C", "Stim-D"),
#' OVAConc = c(0, 5, 50, 500)))
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_sort_by <- function(x,
                         sort_by = NULL) {
  
  # Experiment details per group
  pd_split <- cyto_groups(x,
                          group_by = sort_by,
                          details = TRUE)
  
  # Combine pd
  pd <- do.call("rbind", pd_split)
  
  # Sorting indices
  ind <- match_ind(cyto_names(x), pd[, "name"])
  
  # Return sorted object
  return(x[ind])
  
}

## CYTO_GROUP_BY ---------------------------------------------------------------

#' Group a flowSet or GatingSet by experiment variables
#'
#' @param x an object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param group_by names of cyto_details variables to use for merging, set to
#'   "all" to group all samples in \code{x}. The order of the grouping can be
#'   controlled by specifying the factor levels in a list (e.g. list(Treatment =
#'   c("Stim-A","Stim-C","Stim-B", "Stim-D"))).
#'
#' @return a named list of \code{cytoset} or \code{GatingSet} objects
#'   respectively.
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Group by Treatment
#' cyto_group_by(gs, "Treatment")
#'
#' # Group GatingSet by Treatment and OVAConc
#' cyto_group_by(gs, c("Treatment", "OVAConc"))
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_group_by
#'
#' @export
cyto_group_by <- function(x,
                          group_by = "all") {
  
  # Experiment details per group
  pd_split <- cyto_groups(x, 
                          group_by = group_by,
                          details = TRUE)
  
  # Replace each element of pd_split with matching samples
  x_list <- lapply(seq_len(length(pd_split)), function(z) {
    ind <- match(pd_split[[z]][, "name"], cyto_names(x))
    x[ind]
  })
  names(x_list) <- names(pd_split)
  
  return(x_list)
}

## CYTO_MERGE_BY ---------------------------------------------------------------

#' Merge a cytoset by experiment variables
#'
#' \code{cyto_merge_by} makes a call to \code{cyto_group_by} to split samples
#' into groups based on experiment variables. The resulting groups are then
#' coerced in a \code{cytoset} containing a single \code{cytoframe} with the
#' merged data. \code{cyto_merge_by} is the preferred way to merge samples in
#' CytoExploreR as it will ensure appropriate sampling in \code{cyto_plot}.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to merge when a
#'   \code{GatingHierarchy} or \code{GatingSet} object is supplied, set to the
#'   \code{"root"} node by default.
#' @param merge_by vector of \code{\link{cyto_details}} column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to merging.
#' @param select selection criteria passed to \code{\link{cyto_select}} which
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
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Experiment details
#' cyto_details(gs)
#'
#' # Merge samples by 'Treatment'
#' cs_list <- cyto_merge_by(gs, "Treatment")
#'
#' # Merge samples by 'OVAConc'
#' cs_list <- cyto_merge_by(gs, "OVAConc")
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_group_by}}
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_barcode}}
#' @seealso \code{\link{cyto_sample}}
#'
#' @export
cyto_merge_by <- function(x,
                          parent = "root",
                          merge_by = "all",
                          select = NULL,
                          barcode = TRUE,
                          ...) {
  
  # EXTRACT DATA ---------------------------------------------------------------
  
  # CYTOSET
  x <- cyto_data_extract(x,
                         parent = parent,
                         copy = FALSE)[[1]]
  
  # BARCODING ------------------------------------------------------------------
  
  # SAMPLE ID
  x <- cyto_barcode(x, ...)
  
  # GROUPING -------------------------------------------------------------------
  
  # CYTO_GROUP_BY
  cs_list <- cyto_group_by(x, 
                           group_by = merge_by)
  
  # GROUPS
  grps <- names(cs_list)
  
  # COMBINED EVENTS
  if ("all" %in% names(cs_list)) {
    names(cs_list)[which("all" %in% names(cs_list))] <- "Combined Events"
  }
  
  # SELECTION ------------------------------------------------------------------
  
  # ATTEMPT SELECTION OR RETURN ALL SAMPLES
  if (!is.null(select)) {
    cs_list <- lapply(cs_list, function(z) {
      # Select or return all samples if criteria not met
      tryCatch(cyto_select(z, select), error = function(e) {
        z
      })
    })
  }
  
  # MERGING --------------------------------------------------------------------
  
  # CONVERT EACH GROUP TO MERGED CYTOSET
  structure(lapply(seq_along(cs_list), function(z){
    # CONVERT TO MERGED CYTOSET
    cytoset(structure(list(flowFrame_to_cytoframe(as(z, "flowFrame"))),
                      names = names(cs_list)[z]))
  }), names = names(cs_list))
  
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
  if (!cyto_class(x, "flowFrame")) {
    stop("cyto_split() requires a cytoframe object.")
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

#' Sample a cytoframe, cytoset, GatingHierarchy or GatingSet
#'
#' \code{cyto_sample} allows restriction of a \code{cytoframe}, \code{cytoset},
#' \code{GatingHierarchy} or \code{GatingSet} by indicating the percentage or
#' number of events to retain.
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
#' @importFrom flowCore sampleFilter Subset flowSet exprs
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
  if (nrow(exprs(x)) == 0) {
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
                   "cyto_sample", 
                   display = display, 
                   seed = seed, 
                   copy = FALSE,
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
      cyto_class(z, "flowSet")
    }))) {
      # Same sampling applied to all samples
      x <- lapply(x, function(z) {
        cyto_sample(z, display = display, seed = seed)
      })
      # LIST OF FLOWFRAMES
    } else if (all(LAPPLY(x, function(z) {
      cyto_class(z, "flowFrame")
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

## CYTO_SAMPLE_TO_NODE ---------------------------------------------------------

#' Sample GatingSet to node
#'
#' \code{cyto_sample_to_node} can be used to downsample the root node of a
#' GatingSet based on the number of events in a particular node. A classic
#' example is downsampling to a bead population to give an indication of the
#' relative number of events per sample. For example, if Sample A contains 100
#' beads and 75000 events, and Sample B contains 50 beads and 5000 events.
#' \code{cyto_sample_to_node} will downsample each sample to contain 50 beads
#' and therefore 5000 events for Sample A and 3750 events for Sample B.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param node name of the gated bead population to use for the calculation. If
#'   not supplied internal checks will be made for populations named "Single
#'   Beads" or "Beads".
#' @param count minimum number of events to downsample to, set to the minimum
#'   count for the specified node across samples by default. \code{count} must
#'   be less than or equal to the minimum count across the samples.
#' @param ... additional arguments passed to \code{\link{cyto_sample}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace cytoset flowSet_to_cytoset gs_cyto_data recompute
#'
#' @export
cyto_sample_to_node <- function(x,
                                node = NULL,
                                count = NULL,
                                ...) {
  
  # GATINGSET
  if (!cyto_class(x, "GatingSet", TRUE)) {
    stop("'x' must be a Gatingset object.")
  }
  
  # CLONE GATINGSET
  gs_clone <- cyto_copy(x)
  
  # NODES
  nodes <- cyto_nodes(gs_clone, path = "auto")
  
  # BEADS
  if (is.null(node)) {
    # SINGLE BEADS
    if (any(grepl("single beads", nodes, ignore.case = TRUE))) {
      node <- nodes[which(grepl("single beads", nodes, ignore.case = TRUE))[1]]
      # BEADS
    } else if (any(grepl("beads", nodes, ignore.case = TRUE))) {
      node <- nodes[which(grepl("beads", nodes, ignore.case = TRUE))[1]]
      # BEADS MISSING
    } else {
      stop("Supply the name of the 'node' to downsample the GatingSet.")
    }
  }
  
  # NODE COUNTS
  node_pops <- cyto_data_extract(gs_clone, node)[[node]]
  node_counts <- suppressMessages(
    cyto_stats_compute(node_pops,
                       stat = "count",
                       format = "wide",
                       details = FALSE)
  )
  node_counts <- node_counts[, ncol(node_counts), drop = TRUE]
  
  # NODE MISSING EVENTS
  if (any(node_counts == 0)) {
    ind <- which(node_counts == 0)
    stop(paste0(
      "The following samples do not contain any events in the specified node:",
      paste0("\n", cyto_names(gs_clone)[ind])
    ))
  }
  
  # NODE COUNT
  if (is.null(count)) {
    count <- min(node_counts)
  } else {
    if (!count <= min(node_counts)) {
      count <- min(node_counts)
    }
  }
  
  # NODE RATIOS
  node_ratios <- lapply(seq_len(length(node_counts)), function(z) {
    1 / (node_counts[z] / count)
  })
  
  # SAMPLING - ROOT POPULATION
  pops <- lapply(seq_along(node_pops), function(z) {
    cyto_sample(
      cyto_data_extract(gs_clone[[z]])[["root"]],
      node_ratios[[z]],
      ...
    )
  })
  names(pops) <- cyto_names(node_pops)
  
  # CYTOSET
  if(cyto_class(pops[[1]], "cytoframe", TRUE)) {
    pops <- cytoset(pops)
  } else {
    pops <- flowSet_to_cytoset(flowSet(pops))
  }
  
  # REPLACE DATA IN GATINGSET
  gs_cyto_data(gs_clone) <- pops
  recompute(gs_clone)
  
  # RETURN SAMPLED GATINGSET
  return(gs_clone)
}

#' @export
cyto_beads_sample <- function(...){
  .Defunct("cyto_sample_to_node")
}

## CYTO_BARCODE ----------------------------------------------------------------

#' Barcode each sample or event in a cytoset or GatingSet
#'
#' Adds a new parameter to each of the cytoframes in the cytoset called
#' \code{"Sample ID"} or \code{"Event ID"} to barcode the samples or events from
#' each cytoframe.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} to be barcoded.
#' @param type indicates whether the \code{"samples"}, \code{"events"} or
#'   \code{"both"} should be barcoded, set \code{"samples"} by default.
#'
#' @return barcoded cytoset or GatingSey with \code{"Sample ID"} and/or
#'   \code{"Event ID"} column added and annotated.
#'
#' @importFrom flowWorkspace gs_cyto_data
#'
#' @examples
#' library(CytoExploreRData)
#' 
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Barcode
#' gs <- cyto_barcode(gs)
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_barcode <- function(x,
                         type = "samples") {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CYTOSET
  if(cyto_class(x, "GatingSet", TRUE)) {
    cs <- cyto_data_extract(x,
                            parent = "root")[["root"]]
  } else {
    cs <- x
  }
  
  # TYPE
  if (type == "both") {
    type <- c("samples", "events")
  }
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # SAMPLE ID COLUMN - ONLY IF NOT PRESENT
  if ("samples" %in% type &
      !"Sample ID" %in% cyto_channels(cs)) {
    lapply(seq_along(cs), function(z) {
      suppressWarnings(
        cyto_cbind(cs[[z]],
                   matrix(rep(z, cyto_stat_count(cs[[z]])),
                          ncol = 1,
                          dimnames = list(NULL, "Sample ID")))
      )
    })
  }
  
  # EVENT ID COLUMN - ONLY IF NOT PRESENT
  if ("events" %in% type &
      !"Event ID" %in% cyto_channels(cs)) {
    total_events <- cyto_apply(cs, "cyto_stat_count")
    total_events <- split(
      seq_len(sum(total_events)),
      rep(seq_len(length(cs)),
          times = total_events
      )
    )
    lapply(seq_along(cs), function(z) {
      suppressWarnings(
        cyto_cbind(cs[[z]],
                   matrix(total_events[[z]],
                          ncol = 1,
                          dimnames = list(NULL, "Event ID")))
      )
    })
  }
  
  # RETURN GATINGSET
  if(cyto_class(x, "GatingSet", TRUE)) {
    gs_cyto_data(x) <- cs
    return(x)
  }
  
  # RETURN BARCODED CYTOSET
  return(cs)
  
}

## CYTO_MARKERS_EDIT -----------------------------------------------------------

#' Assign marker names to cytoframe or cytoset
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
#' @param ... not in use.
#'
#' @return save inputs to "Experiment-Markers.csv" and returns updated samples.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters
#' @importFrom DataEditR data_edit
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
  
  # flowFrame
  if (cyto_class(x, "flowFrame")) {
    
    # Extract details of parameters
    pd <- cyto_details(parameters(x))
    
    # flowSet
  } else if (cyto_class(x, "flowSet")) {
    
    # Extract details of parameters
    pd <- cyto_details(parameters(x[[1]]))
    
    # GatingHierarchy
  } else if (cyto_class(x, "GatingHierarchy", TRUE)) {
    fr <- cyto_extract(x, "root")
    pd <- cyto_details(parameters(fr))
    
    # GatingSet
  } else if (cyto_class(x, "GatingSet", TRUE)) {
    fr <- cyto_extract(x, "root")[[1]]
    pd <- cyto_details(parameters(fr))
  } else {
    stop("'x' must be a valid cytometry object.")
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
        mrks <- read_from_csv(found_files[z])
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
    
    # Append file extension
    file <- file_ext_append(file, ".csv")
    
    # File already exists
    if (file_exists(file)) {
      message(
        paste0("Editing data in ", file, "...")
      )
      dt <- read_from_csv(file)
      
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
  
  # Edit dt using data_edit
  if(interactive()) {
    dt <- data_edit(dt,
                    logo = CytoExploreR_logo(),
                    title = "Experiment Markers Editor",
                    row_edit = FALSE, # cannot add/remove rows
                    col_edit = FALSE,
                    col_names = c("channel", "marker"),
                    save_as = file,
                    write_fun = "write_to_csv",
                    quiet = TRUE, 
                    hide = TRUE,
                    ...)
  }
  
  # Update channels
  if(!all(as.character(dt$channel) %in% cyto_channels(x))) {
    cyto_channels(x) <- as.character(dt$channel)
  }
  
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
    if(cyto_class(x, c("flowFrame", "flowSet"), TRUE)){
      flowCore::markernames(x) <- cyto_markers
    }else{
      flowWorkspace::markernames(x) <- cyto_markers
    }
  }
  
  # Return updated samples
  return(x)
}

## CYTO_DETAILS_EDIT -----------------------------------------------------------

#' Interactively edit cyto_details for a cytoset or GatingSet
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param file name of csv file containing experimental information.
#' @param ... not in use.
#'
#' @return NULL and return \code{flowSet} or \code{GatingSet} with updated
#'   experimental details.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore pData<-
#' @importFrom DataEditR data_edit
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
  if (!cyto_class(x, "flowSet") & !cyto_class(x, "GatingSet", TRUE)) {
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
        pdata <- read_from_csv(found_files[z])
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
    
    # File append extension
    file <- file_ext_append(file, ".csv")
    
    # File already exists
    if (length(grep(file, list.files())) != 0) {
      message(
        paste0("Editing data in ", file, "...")
      )
      pd <- read_from_csv(file)
      
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
  
  # REMOVE ROW NAMES
  cyto_names <- pd$name
  rownames(pd) <- NULL
  
  # Edit cyto_details
  if(interactive()) {
    pd <- data_edit(pd,
                    logo = CytoExploreR_logo(),
                    title = "Experiment Details Editor",
                    row_edit = FALSE,
                    col_readonly = "name",
                    save_as = file,
                    write_fun = "write_to_csv",
                    quiet = TRUE,
                    hide = TRUE,
                    ...)
  }
  
  # Update cyto_details
  cyto_details(x) <- pd
  
  # Return updated samples
  return(x)
}

## CYTO_DETAILS_SAVE -----------------------------------------------------------

#' Save experiment details to csv file
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} annotated with
#'   experiment details.
#' @param save_as name of csv file to which the experiment details should be
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
#' @export
cyto_details_save <- function(x,
                              save_as = NULL) {
  
  # SAVE AS
  if (is.null(save_as)) {
    save_as <- paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Details.csv")
  }
  
  # WRITE CSV FILE
  pd <- cyto_details(x)
  write_to_csv(pd,
               save_as)
  
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
    if (cyto_class(spillover, "character", TRUE)) {
      # read in spillover matrix
      spill <- read_from_csv(spilllover)
      # column names must be valid channels
      if (!all(colnames(spill) %in% cyto_channels(fs))) {
        stop(
          paste(
            "'spillover' must have valid fluorescent channels as colnames."
          )
        )
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(fs))
      names(spill) <- cyto_names(fs)
      # spillover is a matrix/data.frame
    } else if (cyto_class(spillover, c("matrix", "data.frame"))) {
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
      if(!cyto_class(cf, "cytoframe", TRUE)){
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
    if (cyto_class(spillover, "character", TRUE)) {
      # read in spillover matrix
      spill <- read_from_csv(spillover)
      # column/row names must be valid channels
      if (!all(colnames(spill) %in% cyto_channels(x))) {
        stop(
          paste(
            "'spillover' must have valid fluorescent channels as colnames."
          )
        )
      }
      # Convert spill into a named list
      spill <- rep(list(spill), length(x))
      names(spill) <- cyto_names(x)
      # spillover is a matrix/data.frame
    } else if (cyto_class(spillover, c("matrix", "data.frame"))) {
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
      if(!cyto_class(cf, "cytoset", TRUE)){
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
    if (cyto_class(spillover, "character", TRUE)) {
      # read in spillover matrix
      spill <- read_from_csv(spillover)
      # column/row names must be valid channels
      if (!all(rownames(spill) %in% cyto_channels(x)) |
          !all(rownames(spill) %in% cyto_channels(x))) {
        stop(
          paste(
            "'spillover' must have valid fluorescent channels as colnames."
          )
        )
      }
      # spillover is a matrix/data.frame
    } else if (cyto_class(spillover, c("matrix", "data.frame"))) {
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
#'
#' @export
cyto_nodes <- function(x, ...) {
  
  # Make call to gh/gs_get_pop_paths
  if (cyto_class(x, "GatingHierarchy", TRUE)) {
    gh_get_pop_paths(x, ...)
  } else if (cyto_class(x, "GatingSet", TRUE)) {
    gs_get_pop_paths(x, ...)
  }
}

## CYTO_NODES_CHECK ------------------------------------------------------------

#' Check if nodes are unique and exist in GatingSet or GatingHierarchy
#'
#' @param x object of class \code{GatingHierarchy} or \code{GatingSet}.
#' @param nodes vector of node paths to check.
#'
#' @return supplied nodes or throw an error if any nodes do not exist or are not
#'   unique within the \code{GatingHierarchy} or \code{GatingSet}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_nodes_check <- function(x,
                             nodes = NULL){
  
  # NO NODES SUPPLIED
  if(is.null(nodes)){
    stop("Supply a vector of nodes to check to the 'nodes' argument.")
  }
  
  # REFERENCE NODE PATHS
  nodes_auto <- cyto_nodes(x, path = "auto")
  nodes_auto_split <- .cyto_nodes_split(nodes_auto)
  
  # SPLIT NODES
  nodes_split <- .cyto_nodes_split(nodes)
  
  # CHECK NODES
  lapply(seq_along(nodes), function(z){
    node <- nodes[z]
    node_split <- nodes_split[[z]]
    # CHECK AGAINST AUTO PATHS
    ind <- which(LAPPLY(nodes_auto_split, function(y){
      # SHORT NODE PATH
      if(length(node_split) < length(y)){
        # PARTIAL
        if(any(node_split %in% y)){
          return(TRUE)
        }else{
          return(FALSE)
        }
        # LONGER NODE PATH
      }else{
        check <- rev(rev(seq_len(length(node_split)))[seq_len(length(y))])
        # PARTIAL OR EXACT
        if(any(node_split[check] %in% y)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }))
    # PARTIAL MATCHES
    if(length(ind) > 1){
      stop(paste0(node, " is not unique in this ", class(x), ".",
                  " Use either ", paste0(nodes_auto[ind], 
                                         collapse = " or "),
                  "."))
    }
    # NO MATCH
    if(length(ind) == 0){
      stop(paste0(node, " does not exist in this ", class(x), "."))
    }
  })
  
  # RETURN NODES
  return(nodes)
  
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
  if(!cyto_class(x, c("GatingHierarchy", "GatingSet"))) {
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
#' @examples 
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' gs_copy <- cyto_copy(gs)
#'
#' @export
cyto_copy <- function(x) {
  
  # GatingSet
  if (cyto_class(x, "GatingSet", TRUE)) {
    x <- gs_clone(x)
  }
  
  # cytoset
  if (cyto_class(x, c("cytoframe", "cytoset"), TRUE)) {
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
#' @importFrom flowCore keyword
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
      sp <- tryCatch(keyword(x[[z]], "SPILL"), error = function(e) {
        NULL
      })
      if(!is.null(sp)) {
        sp <- sp[[1]]
      }
    })
    names(spill) <- cyto_names(x)
    if (all(LAPPLY(spill, "is.null"))) {
      spill <- NULL
    }
    # FLOWFRAME
  } else if (is(x, "flowFrame")) {
    spill <- tryCatch(keyword(x, "SPILL"), error = function(e) {
      NULL
    })
    if (!is.null(spill)) {
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

# Modified version of flowCore fsApply which simplifies tibbles and
# appropriately handles cytosets.

#' Apply a function over cytometry objects
#'
#' @param x object of class \code{flowSet} or \code{cytoset}.
#' @param FUN function to apply to each \code{flowFrame} or \code{cytoframe} in
#'   the supplied \code{flowSet} or \code{cytoset}.
#' @param ... optional arguments to \code{FUN}.
#' @param simplify logical indicating whether the resulting list should be
#'   simplified into a flowSet or bound into a matrix or data.frame using
#'   \code{rbind}, set to TRUE by default.
#' @param use.exprs logical indicating whether the function should be applied to
#'   the raw data in each flowFrame, set to FALSE by default.
#'
#' @importFrom methods is
#' @importFrom flowCore fsApply phenoData phenoData<-
#' @importFrom flowWorkspace flowSet_to_cytoset
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation dataset
#' fs <- Activation
#'
#' # Sample each flowFrame
#' fs <- cyto_apply(fs, cyto_sample)
#'
#' @seealso \code{\link[flowCore:fsApply]{fsApply}}
#'
#' @export
cyto_apply <- function(x, 
                       FUN,
                       ...,
                       simplify = TRUE,
                       use.exprs = FALSE){
  
  # FUNCTION MISSING
  if(missing(FUN)) {
    stop("Supply a function to 'FUN'.")
  }
  
  # PREPARE FUNCTION
  FUN <- match.fun(FUN)
  
  # INVALID FUNCTION
  if(!is.function(FUN)){
    stop("'FUN' is not a function.")
  }
  
  # APPLY FUNCTION
  res <- structure(lapply(cyto_names(x), function(z) {
    y <- x[[z]]
    FUN(if(use.exprs){
      cyto_extract(y, raw = TRUE)[[1]]
    }else{
      y
    }, ...)
  }), names = cyto_names(x))
  
  # SIMPLIFY OUTPUT
  if(simplify) {
    if(all(LAPPLY(res, is, "flowFrame"))) {
      res <- as(res,"flowSet")
      phenoData(res) <- phenoData(x)[cyto_names(x),]
    # CHECK LENGTHS
    } else if(diff(range(LAPPLY(res, length))) == 0) {
      res <- do.call(rbind, res)
    }
  }
  
  # FLOWSET RETURNED
  if(is(res, "flowSet")){
    if(is(x, "cytoset")){
      return(flowSet_to_cytoset(res))
    }else{
      return(res)
    }
  }else{
    return(res)
  }
  
}
