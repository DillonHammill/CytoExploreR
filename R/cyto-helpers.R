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
  files <- list.files(path, 
                      full.names = TRUE,
                      recursive = TRUE) # search internal directories
  
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
      if(length(file_ind) > 0){
        files <- files[unique(file_ind)]
      }
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
      if(length(file_ind) > 0){
        files <- files[-unique(file_ind)]
      }
    }
    
    # SORTED FILE PATHS
    if (sort & length(files) > 1) {
      files <- file_sort(files)
    }
    
    # CYTOSET
    x <- load_cytoset_from_fcs(files = normalizePath(files), ...)
    
    # BARCODE EVENTS - REQUIRED FOR CYTO_PLOT
    x <- cyto_barcode(x, "events")
    
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
#' @param type the method to use when cleaning the data, options include
#'   \code{"flowAI"}, \code{"flowClean"} or \code{"flowCut"}.
#' @param ... additional arguments passed to \code{flowAI::flow_auto_qc},
#'   \code{flowClean::flowClean}, or \code{flowCut::flowCut}.
#'
#' @importFrom flowWorkspace gs_cyto_data flowSet_to_cytoset
#' @importFrom utils capture.output
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
#'
#' # Clean Activation GatingSet
#' gs <- cyto_clean(gs)
#'
#' @references Monaco,G. et al. (2016) flowAI: automatic and interactive anomaly
#'   discerning tools for flow cytometry data. Bioinformatics. 2016 Aug
#'   15;32(16):2473-80.
#'   \url{https://academic.oup.com/bioinformatics/article/32/16/2473/2240408}
#'
#' @references Fletez-Brant K, Spidlen J, Brinkman R, Roederer M, Chattopadhyay
#'   P (2016). flowClean: Automated identification and removal of fluorescence
#'   anaomalies in flow cytometry data. Cytometry A 89(5).
#'   \url{https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22837}
#'
#' @references Meskas J, Wang S, Brinkman R (2021). flowCut --- An R Package for
#'   precise and accurate automated removal of outlier events and flagging of
#'   files based on time versus fluorescence analysis. bioRxiv.
#'   \url{https://doi.org/10.1101/2020.04.23.058545}
#'
#' @seealso \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}
#'
#' @export
cyto_clean <- function(x,
                       type = "flowAI",
                       ...) {
  
  # EXTRACT DATA
  if(cyto_class(x, "GatingSet")) {
    cyto_data <- cyto_data_extract(x,
                                   parent = "root",
                                   copy = FALSE)[["root"]]
  } else {
    cyto_data <- x
  }
  
  # FLOWAI
  if(grepl("ai", type, ignore.case = TRUE)) {
    cyto_require("flowAI",
                 source = "BioC",
                 repo = "giannimonaco/flowAI",
                 ref = paste0("Monaco G, Chen H, Poidinger M, Chen J,",
                              " de Magalhaes J, Larbi A (2016). flowAI:",
                              " automatic and interactive anomaly discerning",
                              " tools for flow cytometry data. Bioinformatics,",
                              " 32(16)."))
    invisible(
      capture.output(
        cyto_data <- cyto_apply(
          cyto_data,
          "flowAI::flow_auto_qc",
          input = "cytoframe",
          copy = FALSE,
          html_report = FALSE,
          mini_report = FALSE,
          fcs_QC = FALSE,
          folder_results = FALSE,
          ...
        )
        # cyto_data <- flow_auto_qc(
        #   cyto_data,
        #   html_report = FALSE,
        #   mini_report = FALSE,
        #   fcs_QC = FALSE,
        #   folder_results = FALSE,
        #   ...
        # )
      )
    )
  # FLOWCLEAN
  } else if(grepl("clean", type, ignore.case = TRUE)) {
    cyto_require("flowClean",
                 source = "BioC",
                 repo = "cafletezbrant/flowClean",
                 ref = paste0(
                   "Fletez-Brant K, Spidlen J, Brinkman R, Roederer M,",
                   " Chattopadhyay P (2016). flowClean: Automated",
                   " identification and removal of fluorescence anaomalies",
                   " in flow cytometry data. Cytometry A 89(5)"
                 ))
    cyto_data <- cyto_apply(
      cyto_data,
      "flowClean::clean",
      input = "flowFrame",
      slot = "frame",
      copy = FALSE,
      diagnostic = FALSE, 
      ...
    )
    
  # FLOWCUT
  } else {
    cyto_require("flowCut",
                 source = "BioC",
                 repo = "jmeskas/flowCut",
                 ref = paste0(
                   "Meskas J, Wang S, Brinkman R (2021). flowCut --- An R",
                   " Package for precise and accurate automated removal of",
                   " outlier events and flagging of files based on time",
                   " versus fluorescence analysis. bioRxiv"))
    cyto_data <- cyto_apply(
      cyto_data,
      "flowCut::flowCut",
      input = "flowFrame",
      slot = "frame",
      copy = FALSE,
      Plot = "None", 
      ...
    )
  }
  
  # CYTOSET
  if(cyto_class(cyto_data, "flowSet", TRUE)) {
    cyto_data <- flowSet_to_cytoset(cyto_data)
  }
  
  # RETURN CLEAN DATA
  if(cyto_class(x, "GatingSet")) {
    gs_cyto_data(x) <- cyto_data
    return(x)
  } else {
    return(cyto_data)
  }
  
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
#'   default. If you need to parse the names using a different delimiter, supply
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
#' cyto_data_extract(gs, "root")[[1]]
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
#' @seealso \code{\link{cyto_gatingTemplate_active}}
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
                     input = "matrix",
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
    cyto_gatingTemplate_active(gatingTemplate)
    
    # CREATE GATINGTEMPLATE
    if(!file_exists(gatingTemplate)){
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
#' cs <- cyto_data_extract(gs, "root")[["root"]]
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
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
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
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
#'  
#' # Activation cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#'  
#' # GatingSet
#' cyto_names(gs)
#'
#' # cytoset
#' cyto_names(cs)
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
  # DEFUNCT
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
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param value vector of replacement names.
#'
#' @importFrom flowWorkspace sampleNames<-
#' @importFrom flowCore identifier
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
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
  if(cyto_class(x, "flowFrame")) {
    identifier(x) <- value
  } else if (cyto_class(x, c("flowSet", "GatingHierarchy", "GatingSet"))) {
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
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
#'                           
#' # Parse file names to variables
#' gs <- cyto_names_parse(
#'   gs,
#'   vars = c("sample_type", "sample_id"),
#'   split = "_"
#'   )
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
#' @param quiet logical which prevents the prointing of messages when set to
#'   TRUE.
#' @param ... additional arguments passed to
#'   \code{\link{cyto_transformers_define}}, when no \code{trans} object is
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
#' @seealso \code{\link{cyto_transformers_define}}
#' @seealso \code{\link{cyto_transformers_combine}}
#' @seealso \code{\link{cyto_transformers_extract}}
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
                                   quiet = FALSE,
                                   ...) {
  
  # No transformations supplied - automatically obtain transform definitions
  if (is.null(trans) | .all_na(trans)) {
    
    # Message not recommended to auto-transform flowFrame/flowSet objects
    if (cyto_class(x, c("flowFrame", "flowSet"))) {
      if(!quiet) {
        message(paste(
          "Automatically transforming cytoframe/cytoset objects",
          "is not recommended as transformation definitions will be lost."
        ))
      }
    }
    
    # Dispatch based on type argument to get TransformerList
    transformer_list <- cyto_transformers_define(x,
                                                 channels = channels,
                                                 parent = parent,
                                                 type = type,
                                                 select = select,
                                                 plot = FALSE, 
                                                 ...)
    
  }
  
  # MESSAGE
  if(!quiet)(
    message(
      paste0("Applying ", 
             ifelse(inverse, "inverse ", ""),
             "data transformations...")
    )
  )
  
  # TRANSFORM FLOWFRAME OR FLOWSET
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    
    # Extract transformations from transformerList to transformList
    transform_list <- cyto_transform_extract(transformer_list,
                                             inverse = inverse
    )
    
    # APPLY TRANSFORMATIONS
    x <- suppressMessages(transform(x, transform_list))
    
    # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (cyto_class(x, "GatingSet")) {
    
    # INVERSE TRANSFORM NOT SUPPORTED
    if (inverse == TRUE) {
      stop(paste(
        "Inverse transformations are not yet supported for",
        "GatingHierarchy/GatingSet objects."
      ))
    }
    
    # APPLY TRANSFORMATIONS
    x <- suppressMessages(transform(x, transformer_list))
  }
  
  # COMPLETE
  if(!quiet){
    message(
      "DONE!"
    )
  }
  
  # VISUALISE TRANSFORMATIONS
  if(plot == TRUE) {
    # INVERSE
    if(inverse == TRUE) {
      transformer_list <- NA
    }
    # PLOT DATA TRANSFORMATIONS
    tryCatch(
      cyto_plot_profile(x,
                        channels = channels,
                        select = select,
                        axes_trans = transformer_list,
                        axes_limits = axes_limits,
                        merge_by = "all",
                        display = 1), # sampling performed above
      error = function(e){
        message("Insufficient plotting space to display transformations!")
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
                                         quiet = FALSE,
                                         ...) {
  
  # Added for backwards compatibility - flowFrame/flowSet objects only
  if (cyto_class(x, "GatingSet")) {
    stop(paste(
      "GatingHierarchy and GatingSet objects require transformerList",
      "objects to apply transformations."
    ))
  }
  
  # MESSAGE
  if(!quiet){
    message(
      paste0("Applying ", 
             ifelse(inverse, "inverse ", ""),
             "data transformations...")
    )
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
  
  # COMPLETE
  if(!quiet){
    message(
      "DONE!"
    )
  }
  
  # VISUALISE TRANSFORMATIONS
  if(plot == TRUE) {
    # PLOT DATA TRANSFORMATIONS
    tryCatch(
      cyto_plot_profile(x,
                        channels = names(trans),
                        axes_trans = trans,
                        axes_limits = axes_limits,
                        merge_by = "all",
                        display = 1, 
                        ...), # sampling performed above
      error = function(e){
        message("Insufficient plotting space to display transformations!")
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
                                           quiet = FALSE,
                                           ...) {
  
  # MESSAGE
  if(!quiet){
    message(
      paste0("Applying ", 
             ifelse(inverse, "inverse ", ""),
             "data transformations...")
    )
  }
  # TRANSFORM FLOWFRAME OR FLOWSET
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    
    # Restrict transformers
    trans <- cyto_transformers_combine(trans[names(trans) %in% cyto_channels(x)])
    
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
    trans <- cyto_transformers_combine(trans[names(trans) %in% cyto_channels(x)])
    
    # Apply transformations
    x <- suppressMessages(transform(x, trans))
  }
  
  # COMPLETE
  if(!quiet){
    message(
      "DONE!"
    )
  }
  
  # VISUALISE TRANSFORMATIONS
  if(plot == TRUE) {
    # INVERSE
    if(inverse == TRUE) {
      transformer_list <- NA
    }
    # PLOT DATA TRANSFORMATIONS
    tryCatch(
      cyto_plot_profile(x,
                        parent = "root",
                        channels = names(trans),
                        axes_trans = trans,
                        axes_limits = axes_limits,
                        merge_by = "all",
                        display = 1,
                        ...), # sampling performed above
      error = function(e){
        message("Insufficient plotting space to display transformations!")
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
#' @param sample passed to \code{\link{cyto_sample}} to extract a subset of
#'   events from each sample. \code{sample} indicates the number or proportion
#'   of events to extract from each sample.
#' @param coerce logical to indicate whether the data should be merged into a
#'   single object prior to returning the data, set to FALSE by default. If
#'   \code{coerce} is TRUE, \code{sample} controls the proportion or number of
#'   events to keep in the merged object.
#' @param barcode additional argument passed to \code{cyto_coerce} when
#'   \code{coerce} is TRUE to control whether samples should be barcoded prior
#'   to coercion, set to FALSE by default.
#' @param seed values used to \code{set.seed()} internally to return the same
#'   extract the same subset of events on each run when \code{sample} is
#'   supplied.
#' @param path can be either \code{"auto"} or \code{"full"} to control whether
#'   the extracted list should be labelled with the shortest or complete path to
#'   each population in parent.
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
#' mt <- cyto_data_extract(gs,
#' parent = "root",
#' format = "matrix",
#' channels = c("CD4", "CD8"))
#' mt[["Activation_1.fcs"]]
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
                              inverse = FALSE,
                              sample = NULL,
                              coerce = FALSE,
                              barcode = FALSE,
                              seed = NULL,
                              path = "auto") {
  
  # TODO: ADD COERCE ARGUMENT? TO MERGE MATRICES ETC.
  
  # PARENT - CYTO_STATS_COMPUTE ALIAS NULL
  if(is.null(parent)) {
    parent = "root"
  }
  
  # PARENTAL PATH
  if(cyto_class(x, "GatingSet")) {
    parent <- cyto_nodes_convert(x,
                                 nodes = parent,
                                 path = path)
  }
  
  # EXTRACT TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # SELECT - CANNOT WORK FOR CYTOFRAMES
  if(!is.null(select)) {
    x <- cyto_select(x, select)
  }
  
  # CYTOFRAME/CYTOSET
  if(cyto_class(x, c("flowFrame", "flowSet"))) {
    cs_list <- list(x)
    # GATINGHIERARCHY/GATINGSET
  } else {
    cs_list <- structure(
      lapply(parent, function(z){
        gs_pop_get_data(x, z)
      }),
      names = parent)
  }
  
  # PREPARE CYTOSET LIST
  res <- lapply(seq_along(cs_list), function(id) {
    # CYTOSET
    cs <- cs_list[[id]]
    # COPY
    if(copy) {
      cs <- cyto_copy(cs)
    }
    # RESTRICT
    if(!is.null(channels)) {
      channels <- cyto_channels_extract(cs, channels)
      cs <- cs[, channels, drop = FALSE]
    }
    # COERCE - NAME WITH PARENT
    if(coerce) {
      cs <- cyto_coerce(
        cs,
        format = "cytoset",
        display = ifelse(is.null(sample), 1, sample),
        name = paste0(names(cs_list)[id], "-merge"),
        barcode = barcode,
        seed = seed
      )
      # SAMPLE
    }else if(!is.null(sample)) {
      cs <- cyto_sample(
        cs,
        display = sample,
        seed = seed
      )
    }
    # TRANSFORM
    if(inverse & !.all_na(trans)) {
      cs <- cyto_transform(cs,
                           trans = trans,
                           inverse = inverse,
                           plot = FALSE,
                           quiet = TRUE)
    }
    # CYTOFRAME
    if(format == "cytoframe") {
      # CYTOFRAME -> CYTOFRAME
      if(cyto_class(cs, "flowFrame")) {
        list(cs) # CANNOT CALL CYTO_NAMES HERE
        # CYTOSET -> CYTOFRAME
      } else {
        structure(
          lapply(cyto_names(cs), function(z){
            cs[[z]]
          }), names = cyto_names(cs)
        )
      }
      # CYTOSET
    } else if(format == "cytoset") {
      # CYTOFRAME -> CYTOSET
      if(cyto_class(cs, "flowFrame")) {
        cytoset(
          structure(
            list(cs), names = "cf"
          )
        )
        # CYTOSET -> CYTOSET
      } else {
        if(split) {
          structure(lapply(seq_along(cs), function(z){
            cs[z]
          }), names = cyto_names(cs))
        } else {
          cs
        }
      }
      # MATRIX
    } else if(format == "matrix") {
      # CYTOFRAME -> MATRIX
      if(cyto_class(cs, "flowFrame")) {
        structure(
          list(
            cyto_exprs(cs, 
                       markers = markers,
                       drop = FALSE)
          ), names = "cf_raw"
        )
        # CYTOSET -> MATRIX
      } else {
        structure(
          lapply(cyto_names(cs), function(z){
            cyto_exprs(cs[[z]],
                       markers = markers,
                       drop = FALSE)
          }), names = cyto_names(cs)
        )
      }
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

## CYTO_EXPRS ------------------------------------------------------------------

#' Extract raw data matrix from cytoframe or cytoset
#'
#' A convenient wrapper around \code{exprs} which simplifies extraction of raw
#' data from certain parameters. \code{cyto_exprs} is primarily for use within
#' CytoExploreR and users should instead use the \code{cyto_data_extract} API.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}} or
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}.
#' @param channels vector of channel or marker names for which raw data should
#'   be extracted.
#' @param markers logical indicating whether the column names of the extracted
#'   matrix should be converted to markers where possible, set to FALSE by
#'   default.
#' @param drop logical to control whether single channel extractions should be
#'   converted to vectors by default, set to TRUE by default.
#' @param ... additional arguments passed to data extraction from matrix.
#'
#' @return matrix or list of matrices containing raw data.
#'
#' @importFrom flowCore exprs
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' #' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # cytoset
#' cs <- cyto_data_extract(gs, "T Cells")[["T Cells"]]
#'
#' # matrices
#' cyto_exprs(cs, c("CD44", "CD69"))
#'
#' @seealso \code{\link{cyto_data_extract}}
#'
#' @name cyto_exprs
NULL

#' @rdname cyto_exprs
#' @export
cyto_exprs <- function(x,...) {
  UseMethod("cyto_exprs")
}

#' @rdname cyto_exprs
#' @export
cyto_exprs.flowFrame <- function(x, 
                                 channels = NULL,
                                 markers = FALSE,
                                 drop = TRUE,
                               ...) {
  # CHANNELS
  if(is.null(channels)) {
    mt <- exprs(x)[, 
                   , 
                   drop = drop,
                   ...]
  } else {
    mt <- exprs(x)[, 
                   cyto_channels_extract(x, channels), 
                   drop = drop, 
                   ...]
  }

  # MARKERS
  if(markers) {
    colnames(mt) <- unname(cyto_markers_extract(x, colnames(mt)))
  }
  
  return(mt)
  
}

#' @rdname cyto_exprs
#' @export
cyto_exprs.flowSet <- function(x,
                               channels = NULL,
                               markers = FALSE,
                               drop = TRUE,
                               ...) {
  structure(
    lapply(seq_along(x), function(z){
      cyto_exprs(x[[z]],
                 channels = channels,
                 markers = markers,
                 drop = drop,
                 ...)
    }),
    names = cyto_names(x)
  )
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
#' # Look at experiment details
#' cyto_details(gs)
#'
#' # Select Stim-C samples with 0 and 500 nM OVA concentrations
#' fs <- cyto_filter(
#'   gs,
#'   Treatment == "Stim-C",
#'   OVAConc %in% c(0, 500)
#' )
#'
#' # Select Stim-A and Stim-C treatment groups
#' fs <- cyto_filter(
#'   gs,
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
    ind <- match(rownames(pd_filter), rownames(pd))
  }
  
  return(x[ind])
}


## CYTO_MATCH ------------------------------------------------------------------

#' Get index of samples matching experimental criteria
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... named list containing experimental variables to be used to select
#'   samples or named arguments containing the levels of the variables to
#'   select. See below examples for use cases. Neagtive exclusion indices can be
#'   returned by setting \code{exclude = TRUE}.
#'
#' @return vector of indices that representat the location of the selected
#'   sample(s) within the supplied \code{cytoset} or \code{GatingSet}.
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
#' # Look at experiment details
#' cyto_details(gs)
#'
#' # Indices of Stim-C samples with 0 and 500 nM OVA concentrations
#' cyto_match(gs,
#'   Treatment = "Stim-C",
#'   OVAConc = c(0, 500)
#' )
#'
  #' # Indices of Stim-A and Stim-C treatment groups
#' cyto_match(
#'   gs,
#'   list("Treatment" = c("Stim-A", "Stim-C"))
#' )
#'
#' # Exclusion indices of Stim-D treatment group
#' cyto_match(gs,
#'   Treatment = "Stim-D",
#'   exclude = TRUE
#' )
#' 
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_filter}}
#' 
#' @export
cyto_match <- function(x, 
                       ...) {
  
  # CYTOSET/GATINGSET
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    stop(
      "'x' should be an object of class cytoset or GatingSet."
    )
  }
  
  # ARGUMENTS
  args <- list(...)
  
  # ... NAMED LIST OF ARGUMENTS
  if(cyto_class(args[[1]], "list")) {
    args <- args[[1]]
  }
  
  # EXCLUDE
  if(any(grepl("exclude", names(args)))) {
    exclude <- args[[which(grepl("exclude", names(args)))]]
    args <- args[-which(grepl("exclude", names(args)))]
  } else {
    exclude <- FALSE
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # INDICES/NAMES
  if(length(args) == 1 & (is.null(names(args)) | .empty(names(args)))) {
    # NULL
    if(is.null(unlist(args))) {
      ind <- seq_along(x)
    # INDICES SUPPLIED
    } else if(is.numeric(unlist(args))) {
      ind <- unlist(args)
    # UNNAMED CHARACTERS
    } else {
      # MATCH NAMES - ROWNAMES/NAME/PARTIAL
      ind <- LAPPLY(unlist(args), function(z){
        # ROWNAMES - EXACT MATCH
        if(z %in% rownames(pd)) {
          match(z, rownames(pd))
        # ROWNAMES - PARTIAL MATCH
        } else if(any(grepl(z, rownames(pd), ignore.case = TRUE))) {
          which(grepl(z, rownames(pd), ignore.case = TRUE))
        # NAME COLUMN - EXACT MATCH
        } else if(z %in% pd[, "name"]) {
          match(z, pd[, "name"])
        # NAME COLUMN - PARTIAL MATCH
        } else if(any(grepl(z, pd[, "name"], ignore.case = TRUE))) {
          which(grepl(z, pd[, "name"], ignore.case = TRUE))
        # NO MATCH
        } else {
          NULL
        }
      })
    }
  # EXPERIMENTAL VARIABLES
  } else {
    # VALID VARIABLES
    if (!all(names(args) %in% colnames(pd))) {
      lapply(names(args), function(y) {
        if (!y %in% names(pd)) {
          stop(paste(y, "is not a valid variable in cyto_details(x)."))
        }
      })
    }
    
    # VARIABLE LEVELS EXIST
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
    
    # FILTERED EXPERIMENT DETAILS
    pd_filter <- pd
    lapply(names(args), function(y) {
      ind <- which(pd_filter[, y] %in% args[[y]])
      # SKIP FILTER IF VARIABLE LEVEL MISSING
      if (length(ind) != 0) {
        pd_filter <<- pd_filter[ind, , drop = FALSE]
      }
    })
    
    # INDICES
    ind <- match(rownames(pd_filter), rownames(pd))
  }
  
  # NEGATIVE EXCLUSION INDICES
  if (exclude == TRUE) {
    return(-ind)
  # POSITIVE SELECTION INDICES
  } else {
    return(ind)
  }
  
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
#' @seealso \code{\link{cyto_match}}
#' @seealso \code{\link{cyto_filter}}
#'
#' @export
cyto_select <- function(x, 
                        ...) {
  
  # MATCH - INDICES
  ind <- cyto_match(x,
                    ...)
  
  # RETURN SELECTED SAMPLES
  return(x[ind])
  
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
#' @param select vector of indices or named list containing experimental
#'   variables passed to \code{\link{cyto_select}} to be used to select samples
#'   prior to retrieving group-wise experimental details.
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
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_group_by}}
#' @seealso \code{\link{cyto_merge_by}}
#' @seealso \code{\link{cyto_sort_by}}
#'
#' @export
cyto_groups <- function(x, 
                        group_by = "all",
                        details = FALSE,
                        select = NULL){
  
  # Check class of x
  if (!cyto_class(x, c("flowSet", "GatingSet"))) {
    stop("'x' should be an object of class cytoset or GatingSet.")
  }
  
  # Extract experiment details
  pd <- cyto_details(x)
  
  # Replace any NA with "NA" to avoid missing rows
  if (any(is.na(pd))) {
    pd[is.na(pd)] <- "NA"
  }
  
  # group_by is a list with factor levels - should not be "all"
  if (cyto_class(group_by, "list")) {
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
      pd_split <- lapply(cyto_names(x), function(z) {
        pd[rownames(pd) == z, , drop = FALSE] # name column may not match
      })
      names(pd_split) <- cyto_names(x)
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
  
  # CONVERT ALL -> COMBINED EVENTS
  names(pd_split)[names(pd_split) == "all"] <- "Combined Events"
  
  # RETURN SPLIT DETAILS
  if(details == TRUE){
    return(pd_split)
    # RETURN GROUP NAMES
  }else{
    return(names(pd_split))
  }
  
}

## CYTO_SORT_BY ----------------------------------------------------------------

#' Sort a cytoset or GatingSet by experiment variables
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param sort_by names of cyto_details variables to use for sorting. For more
#'   control over the sorting order, specify the factor levels for each variable
#'   in a list (e.g. list(Treatment = c("Stim-A","Stim-C","Stim-B", "Stim-D"))).
#' @param select vector of indices or named list containing experimental
#'   variables passed to \code{\link{cyto_select}} to be used to select samples
#'   prior to sorting.
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
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_groups}}
#' @seealso \code{\link{cyto_group_by}}
#' @seealso \code{\link{cyto_merge_by}}
#'
#' @export
cyto_sort_by <- function(x,
                         sort_by = NULL,
                         select = NULL) {
  
  # Experiment details per group
  pd_split <- cyto_groups(x,
                          select = select,
                          group_by = sort_by,
                          details = TRUE)
  
  # Combine pd
  pd <- do.call("rbind", pd_split)
  
  # Sorting indices
  ind <- match_ind(rownames(cyto_details(x)), rownames(pd))
  
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
#' @param select vector of indices or named list containing experimental
#'   variables passed to \code{\link{cyto_select}} to be used to select samples
#'   prior to grouping.
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
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_groups}}
#' @seealso \code{\link{cyto_sort_by}}
#' @seealso \code{\link{cyto_merge_by}}
#'
#' @export
cyto_group_by <- function(x,
                          group_by = "all",
                          select = NULL) {
  
  # Experiment details per group
  pd_split <- cyto_groups(x, 
                          select = select,
                          group_by = group_by,
                          details = TRUE)
  
  # Replace each element of pd_split with matching samples
  x_list <- lapply(seq_len(length(pd_split)), function(z) {
    ind <- match(rownames(pd_split[[z]]), rownames(cyto_details(x)))
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
#' @param format either \code{"cytoframe"} or \code{"cytoset"} to indicate the
#'   desired format for merged groups, set to \code{"cytoset"} by default.
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
#' @seealso \code{\link{cyto_groups}}
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
                          format = "cytoset",
                          ...) {
  
  # EXTRACT DATA ---------------------------------------------------------------
  
  # CYTOSET
  x <- cyto_data_extract(x,
                         parent = parent,
                         copy = FALSE)[[1]]
  
  # BARCODING ------------------------------------------------------------------
  
  # SAMPLE ID
  if(barcode) {
    x <- cyto_barcode(x, ...)
  }
  
  # GROUPING -------------------------------------------------------------------
  
  # CYTO_GROUP_BY
  cs_list <- cyto_group_by(x, 
                           group_by = merge_by)
  
  # COMBINED EVENTS
  if ("all" %in% names(cs_list)) {
    names(cs_list)[which("all" %in% names(cs_list))] <- "Combined Events"
  }
  
  # SELECTION ------------------------------------------------------------------
  
  # ATTEMPT SELECTION OR RETURN ALL SAMPLES
  if (!is.null(select)) {
    cs_list <- lapply(cs_list, function(z) {
      tryCatch(cyto_select(z, select), error = function(e) {
        z
      })
    })
  }
  
  # MERGING --------------------------------------------------------------------
  
  # CONVERT EACH GROUP TO MERGED CYTOSET/FLOWSET
  if(grepl("s", format, ignore.case = TRUE)) {
    structure(
      lapply(seq_along(cs_list), function(z){
        do.call(
          cyto_class(cs_list[[z]]), # flowSet() / cytoset()
          list(
            structure(
              list(
                as(cs_list[[z]], 
                   cyto_class(cs_list[[z]][[1]])) # flowFrame / cytoframe
              ),
              names = names(cs_list)[z]
            )
          )
        )
      }),
      names = names(cs_list)
    )
  # CONVERT EACH GROUP TO CYTOFRAME
  } else {
    structure(
      lapply(seq_along(cs_list), function(z){
        as(cs_list[[z]], cyto_class(cs_list[[z]][[1]])) # flowFrame / cytoframe
      }),
      names = names(cs_list)
      )
  }
  
}

## CYTOFRAME COERCION METHODS --------------------------------------------------

# These are wrappers for cytoframe coercion methods that are currently 
# missing in flowWorkspace (not exported to avoid conflict with flowWorkspace).

#' Coerce cytoset to cytoframe
#' 
#' @importFrom methods setAs
#' @importFrom flowWorkspace flowFrame_to_cytoframe
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs("cytoset",
      "cytoframe",
      function(from){
        flowFrame_to_cytoframe(
          as(from, "flowFrame")
        )
      })

#' Coerce flowSet to cytoframe
#' 
#' @importFrom methods setAs
#' @importFrom flowWorkspace flowFrame_to_cytoframe
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs("flowSet",
      "cytoframe",
      function(from){
        flowFrame_to_cytoframe(
          as(from, "flowFrame")
        )
      })

## CYTO_SPLIT ------------------------------------------------------------------

#' Split samples merged with cyto_merge
#'
#' Extract individual samples merged using \code{cyto_merge_by()} or
#' \code{cyto_coerce()} based on \code{"Sample-ID"} column created by
#' \code{cyto_barcode()}.
#'
#' @param x object of class \code{cytoframe} or \code{cytoset}.
#' @param id can be either the name of a channel or marker containing the IDs
#'   for samples or a vector containing an index for each event within \code{x},
#'   set to \code{Sample-ID} channel by default.
#' @param names vector of names to assign to each of the extracted cytoframes
#'   when saving the split files. Names should be supplied in the order used
#'   prior to merging (i.e. in Sample-ID order).
#'
#' @return a cytoset containing the split cytoframes.
#'
#' @importFrom flowCore identifier split
#' @importFrom utils type.convert
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
                       id = "\\Sample.*ID$",
                       names = NULL) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # FLOWFRAME
  if (!cyto_class(x, c("flowFrame", "flowSet"))) {
    stop("cyto_split() requires a cytoframe or cytoset object.")
  }
  
  # ID - CHANNEL NAME
  if(is.character(id)) {
    id <- cyto_channels_extract(x, id)[1] # FIRST MATCH
  # ID - VECTOR
  } else {
    # EVENTS
    if(cyto_class(x, "flowFrame")) {
      cnts <- matrix(
        nrow(x),
        dimnames = list(NULL, "nrow")
      )
    } else {
      cnts <- cyto_apply(
        x,
        FUN = "nrow",
        input = "matrix",
        copy = FALSE
      )
    }
    # INSUFFICIENT IDS
    if(length(id) != sum(cnts)) {
      stop(
        "'id' must contain an ID for every event in 'x'!"
      )
    # NUMERIC/FACTOR IDS
    } else {
      id <- split(
        id,
        LAPPLY(nrow(cnts), function(z){
          rep(z, cnts[z, 1])
        }),
        drop = FALSE
      )
    }
  }
  
  # SPLIT INTO CYTOFRAMES ------------------------------------------------------
  
  # CYTOFRAME -> CYTOSET
  if(cyto_class(x, "flowFrame")) {
    x <- cytoset(
      list("cytoframe" = x)
    )
  }
  
  # LOOP THROUGH CYTOSET
  cf_list <- lapply(seq_along(x), function(z){
    # CYTOFRAME
    cf <- cyto_data_extract(
      x[z],
      format = "cytoframe",
      copy = FALSE
    )[[1]][[1]]
    # IDS
    if(cyto_class(id, "list", TRUE)) {
      ids <- id[[z]]
      if(!cyto_class(ids, "factor")) {
        ids <- factor(ids)
      }
    # IDS - STORED IN A CHANNEL
    } else {
      ids <- factor(
        cyto_exprs(
          cf, 
          channels = id,
          drop = TRUE
        )
      )
    }
    # BYPASS EMPTY SAMPLES
    cf_list <- split(
      cf,
      ids
    )
    # NAMES
    if (!is.null(names)) {
      names(cf_list) <- tryCatch(
        names[seq(1, length(cf_list))],
        error = function(e){
          stop(
            "Insufficient sample names passed to 'names'!"
          )
        }
      )
      # REMOVE USED NAMES
      names <<- names[-seq( 1, length(cf_list))]
      # NO NAMES
    } else {
      # FACTOR LEVELS - MAY BE ASSIGNED
      if(!is.character(type.convert(names(cf_list), as.is = TRUE))) {
        names(cf_list) <- paste0("Sample-", names(cf_list))
      }
    }
    # IDENTIFIERS - NOT REQUIRED
    cf_list <- structure(
      lapply(seq_along(cf_list), function(v) {
        identifier(cf_list[[v]]) <<- names(cf_list)[v]
        return(cf_list[[v]])
      }),
      names = names(cf_list)
    )
    return(cf_list)
  })
  
  # FLATTEN CYTOFRAME LIST
  cf_list <- unlist(cf_list)
  
  # RETURN CYTOSET
  return(cytoset(cf_list))
}

## CYTO_SAVE -------------------------------------------------------------------

#' Write samples to FCS files in new folder or save GatingSet
#'
#' @param x object of class \code{flowSet}, \code{GatingHierarchy} or
#'   \code{GatingSet}.
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
#' @param inverse logical indicating whether the data should be inverse
#'   transformed prior to writing FCS files, set to FALSE by default. Inverse
#'   transformations of \code{flowFrame} or \code{flowSet} objects requires
#'   passing of transformers through the \code{trans} argument.
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
    message(paste("Saving GatingSet to ", save_as, "..."))
    suppressMessages(save_gs(x, save_as))
    # RETURN GATINGSET
    invisible(x)
    # SAVE FCS FILES
  } else {
    # EXTRACT DATA
    message(paste("Extracting the ", parent, " node from the GatingSet."))
    cs <- cyto_data_extract(
      x,
      parent = parent,
      copy = TRUE
    )[[1]]
    # TRANSFORMATIONS
    trans <- cyto_transformers_extract(x)
    # FLOWSET METHOD
    cs <- cyto_save(
      x = cs,
      split = split,
      names = names,
      save_as = save_as,
      inverse = inverse,
      trans = trans
    )
    # RETURN DATA
    invisible(cs)
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
    cs <- cyto_data_extract(
      x,
      parent = parent,
      copy = TRUE
    )[[1]]
    # TRANSFORMATIONS
    trans <- cyto_transformers_extract(x)
    # FLOWSET METHOD
    cs <- cyto_save(
      x = cs,
      split = split,
      names = names,
      save_as = save_as,
      inverse = inverse,
      trans = trans
    )
    
    # RETURN DATA
    invisible(cs)
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
  
  # CYTOSET OF SPLIT CYTOFRAMES
  if (split == TRUE) {
    x <- cyto_split(x, 
                    names = names)
  }
  
  # DIRECTORY CHECK
  if (!is.null(save_as) & dir.exists(save_as)) {
    # FILES WILL BE OVERWRITTEN
    if (any(list.files(save_as) %in% cyto_names(x))) {
      message(paste0("Files will be overwritten in ", save_as, "."))
      if(!cyto_enquire(
        "Do you want to continue? (Y/N)",
        options = c("T", "Y")
      )) {
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
  
  # INVERSE TRANSFORM?
  if(inverse == TRUE) {
    # TRANSFORMERS REQUIRED
    if (is.null(trans) | .all_na(trans)) {
      stop("Supply transformerList to 'trans' to inverse transformations.")
    }
    # INVERSE TRANSFORM
    x <- cyto_transform(cyto_copy(x),
                        trans = trans,
                        inverse = inverse,
                        plot = FALSE)
  }
  
  # WRITE FCS FILES
  lapply(seq_along(x), function(z){
    # Message
    message(paste0(cyto_names(x)[z], "..."))
    # NO DIRECTORY SPECIFIED
    if (is.null(save_as)) {
      write.FCS(
        x[[z]],
        cyto_names(x)[z]
      )
      # DIRECTORY SPECIFIED
    } else {
      # CREATE DIRECTORY
      if (!dir.exists(save_as)) {
        dir.create(save_as)
      }
      write.FCS(
        x[[z]],
        paste0(save_as, "/", cyto_names(x)[z])
      )
    }
  })
  
  # RETURN DATA
  invisible(x)
}

## CYTO_LIST -------------------------------------------------------------------

#' Convert cytoset or GatingSet to a list of its elements
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param extract logical indicating whether the data should be extracted from
#'   the cytoset or GatingSet container to a cytoframe or GatingHierarchy
#'   respectively, set to TRUE by default.
#'
#' @return a list of cytoset or GatingSet elements.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_list <- function(x,
                      extract = TRUE) {
  
  # CYTOFRAME SUPPLIED
  if(cyto_class(x, "flowFrame")) {
    return(list(x)) # no names
  }
  
  # CYTOFRAME
  if(extract) {
    structure(
      lapply(seq_along(x), function(z){
        x[[z]]
      }), 
      names = cyto_names(x)
    )
  # CYTOSET
  } else {
    structure(
      lapply(seq_along(x), function(z){
        x[z]
      }),
      names = cyto_names(x)
    )
  }
  
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
#'   percentage or number of events to keep respectively. For \code{cytoset} or
#'   \code{GatingSet} objects, the same degree of sampling is applied to each
#'   \code{cytoframe} by default. However, cytoframes will be separately sampled
#'   if a vector the same length as the \code{cytoset} or \code{GatingSet} is
#'   supplied.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
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
  
  # PREPARE DISPLAY - INDIVIDUAL CYTOFRAMES
  display <- rep(display, 
                 length.out = length(x))
  
  # SAPLE EACH CYTOFRAME
  cytoset(
    structure(
      lapply(seq_along(x), function(z){
        cyto_sample(x[[z]],
                    display = display[z],
                    seed = seed,
                    ...)
      }),
      names = cyto_names(x)
    )
  )

}

#' @rdname cyto_sample
#' @export
cyto_sample.list <- function(x,
                             display = 1,
                             seed = NULL,
                             ...) {
  
  # SAME SAMPLE SIZE PER LAYER
  display <- rep(display, length(x))
  
  # SAMPLING
  x <- mapply(function(x, display) {
    cyto_sample(x, 
                display = display,
                seed = seed)
  }, x, display)
  
  # Return sampled list
  return(x)
}

## CYTO_SAMPLE_N ---------------------------------------------------------------

#' Compute sample number for each cytoframe in a cytoset
#'
#' This function is used within \code{cyto_coerce} and \code{cyto_plot} to
#' compute the number of events to extract from each cytoframe in cytoset prior
#' to coercion.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param display can be either a numeric [0,1] or integer to indicate the
#'   percentage or number of events to keep respectively.
#' @param parent name of the parental population to extract from GatingHierarchy
#'   or GatingSet objects.
#'
#' @return named vector containing the number of events to extract from each
#'   cytoframe within the cytoset.
#'
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
#'
#'  # SAMPLE SIZE
#'  cyto_sample_n(gs,
#'  parent = "T Cells",
#'  display = 50000)
#'
#' @export
cyto_sample_n <- function(x,
                          display = 1,
                          parent = NULL) {
  
  # COUNTS
  counts <- cyto_apply(x,
                       "nrow",
                       parent = parent,
                       input = "matrix",
                       copy = FALSE)
  
  # FORMAT
  counts <- structure(
    counts[, 1, drop = TRUE],
    names = rownames(counts)
    )
  
  # TOTAL COUNT
  total_count <- sum(counts)
  
  # DISPLAY CANNOT BE LARGER THAN TOTAL COUNT
  if(display > 1 & display > total_count) {
    display <- total_count # display may be set to zero here
  }
  
  # DISPLAY - COUNT
  if(display <= 1) {
    display <- display * total_count
  }

  # DISPLAY - FREQ
  if(display == 0) {
    display_freq <- 0 # avoid possible NaN
  } else {
    display_freq <- display / total_count
  }
  
  # APPROX COUNTS
  sample_counts <- ceiling(display_freq * counts)
  sample_total <- sum(sample_counts)
  
  # EXACT COUNTS - REMOVE EXCESS EVENTS
  excess <- sample_total - display
  if(excess > 0) {
    sample_ind <- which(sample_counts != 0)
    for(i in seq_len(excess)) {
      ind <- sample(sample_ind, 1)
      sample_counts[ind] <- sample_counts[ind] - 1
      if(sample_counts[ind] == 0) {
        sample_ind <- sample_ind[-match(names(ind), names(sample_ind))]
      }
    }
  }
  
  # SAMPLE COUNTS
  return(sample_counts)
  
}

## CYTO_COERCE -----------------------------------------------------------------

#' Coerce cytoframes in a cytoset
#'
#' \code{cyto_coerce} is an efficient method for simultaneously coercing and
#' downsampling elements of a \code{cytoset}. Basically, the individual
#' \code{cytoframe} elements are sampled prior to coercion using a combination
#' of \code{cyto_sample_n} and \code{cyto_sample}. \code{format} provides
#' control over the format in which the coerced data should be returned.
#'
#' @param x \code{\link[flowWorkspace:cytoset]{cytoset}}.
#' @param display passed to \code{\link{cyto_sample}} to control the number of
#'   events to retain in the coerced \code{cytoframe}, set to 1 by default to
#'   retain all events.
#' @param seed numeric passed to \code{\link{set.seed}} to return the same
#'   sample with each run, set to NULL by default for random sampling.
#' @param barcode logical indicating whether the samples should be barcoded
#'   prior to coercion, set to FALSE by default.
#' @param format can be either 1 - matrix, 2 - cytoframe or 3 - cytoset to
#'   control the format of the returned data, set to \code{"cytoset"} by
#'   default.
#' @param name character string indicating the name to use for the coerced
#'   object, set to \code{"merge"} by default.
#' @param overwrite logical passed to \code{cyto_barcode} to control whether
#'   existing barcodes should be overwritten, the user will be asked
#'   interactively if not manually supplied.
#' @param ... additional arguments passed to \code{\link{cyto_sample}}.
#'
#' @return a sampled and coerced \code{matrix}, \code{cytoframe} or
#'   \code{cytoset}.
#'
#' @importFrom flowWorkspace flowFrame_to_cytoframe cytoset
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(
#'  system.file(
#'    "extdata/Activation-GatingSet",
#'    package = "CytoExploreRData"
#'    )
#'  )
#'
#' # Activation cytoset
#' cs <- cyto_data_extract(gs, "root")[["root"]]
#'
#' # Coerce & Sample
#' cyto_coerce(cs, display = 50000, seed = 56)
#'
#' @export
cyto_coerce <- function(x,
                        display = 1,
                        seed = NULL,
                        barcode = FALSE,
                        format = "cytoset",
                        name = "merge",
                        overwrite = NULL,
                        ...) {
  # IDENTIFIERS
  ids <- cyto_names(x)
  
  # SAMPLE COUNTS
  cnts <- cyto_sample_n(x,
                        display = display)
  
  # REMOVE EMPTY SAMPLES
  if(any(!names(cnts) %in% ids)) {
    x <- cyto_select(x,
                     list("name" = ids[ids %in% names(cnts)]))
  }

  # SAMPLING
  x <- cytoset(
    structure(
      lapply(
        cyto_names(x), 
        function(z){
          fr <- cyto_sample(x[[z]], 
                            display = cnts[z], 
                            seed = seed)
          if(cyto_class(fr, "flowFrame", TRUE)) {
            fr <- flowFrame_to_cytoframe(fr)
          }
          return(fr)
        }
      ),
    names = cyto_names(x))
  )

  # BARCODE
  if(barcode) {
    cyto_barcode(x,
                 overwrite = overwrite)
  }
  
  # COERCE
  if(length(x) == 1) {
    x <- x[[1]]
  } else {
    x <- flowFrame_to_cytoframe(
      as(x, "flowFrame")
    )
  }
  
  # FORMAT
  if(format == "cytoset") {
    x <- cytoset(
      structure(
        list(x),
        names = name
      )
    )
  }
  
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
  node_pops <- cyto_data_extract(gs_clone, 
                                 parent = node,
                                 copy = FALSE)[[node]]
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
      cyto_data_extract(gs_clone[[z]],
                        parent = "root",
                        copy = FALSE)[["root"]][[1]],  # cytoframe
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
#' @param overwrite logical to indicate whether existing barcodes should be
#'   overwritten, thus providing a non-interactive way to control how existing
#'   barcodes are handled.
#'
#' @return barcoded cytoset or GatingSey with \code{"Sample ID"} and/or
#'   \code{"Event ID"} column added and annotated.
#'
#' @importFrom flowWorkspace gs_cyto_data realize_view cytoset
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
                         type = "samples",
                         overwrite = NULL) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CYTOSET
  cs <- cyto_data_extract(x,
                          parent = "root",
                          copy = FALSE,
                          format = "cytoset")[[1]]
  
  # TYPE
  if (grepl("^b", type, ignore.case = TRUE)) {
    type <- c("samples", "events")
  }
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # BARCODE?
  barcode <- TRUE
  
  # CHECK FOR SAMPLE IDs
  if (any(grepl("^s", type, ignore.case = TRUE))){
    # SAMPLE IDs EXIST - BACKWARDS COMPATIBLE
    if(any(grepl("\\Sample.*\\ID$", cyto_channels(cs), ignore.case = TRUE))) {
      # OVERWRITE
      if(!is.logical(overwrite)) {
        overwrite <- cyto_enquire(
          "Override existing sample IDs? (Y/N): ",
          options = c("T", "Y")
        )
      }
      # OVERWRITE BARCODES
      if(overwrite) {
        # REMOVE SAMPLE ID COLUMN
        cs <- realize_view(cs[, -which(grepl("\\Sample.*\\ID$", 
                                             cyto_channels(cs), 
                                             ignore.case = TRUE))])
      } else {
        barcode <- FALSE
      }
    }
    # BARCODE SAMPLES
    if(barcode){
      cs <- cytoset(
        structure(
          lapply(seq_along(cs), function(z) {
            suppressWarnings(
              cyto_cbind(cs[[z]],
                         matrix(rep(z, cyto_stat_count(cs[[z]])),
                                ncol = 1,
                                dimnames = list(NULL, "Sample-ID")))
            )
          }),
          names = cyto_names(cs)
        )
      )
    }
  }
  
  # CHECK FOR EVENT IDs
  if (any(grepl("^e", type, ignore.case = TRUE))) {
    # EVENT IDs EXIST
    if(any(grepl("\\Event.*\\ID$", cyto_channels(cs), ignore.case = TRUE))) {
      # OVERWRITE
      if(!is.logical(overwrite)) {
        overwrite <- cyto_enquire(
          "Override existing event IDs? (Y/N): ",
          options = c("T", "Y")
        )
      }
      # OVERWRITE BARCODES
      if(overwrite){
        # REMOVE EVENT ID COLUMN
        cs <- realize_view(cs[, -which(grepl("\\Event.*\\ID$", 
                                             cyto_channels(cs), 
                                             ignore.case = TRUE))])
      } else {
        barcode <- FALSE
      }
    }
    # BARCODE EVENTS
    if(barcode) {
      # EVENT IDs
      cnt <- 0
      event_ids <- lapply(seq_along(cs), function(z){
        if(nrow(cs[[z]]) == 0) {
           return(NA)
        } else {
          ids <- seq(cnt + 1, cnt + nrow(cs[[z]]))
          cnt <<- ids[length(ids)]
          return(ids)
        }
      })
      cs <- cytoset(
        structure(
          lapply(seq_along(cs), function(z) {
            suppressWarnings(
              cyto_cbind(cs[[z]],
                         matrix(event_ids[[z]],
                                ncol = 1,
                                dimnames = list(NULL, "Event-ID")))
            )
          }),
          names = cyto_names(cs)
        )
      )
    }
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
  
  # CHANNELS
  chans<- cyto_channels(x)
  
  # MARKERS
  marks<- cyto_markers(x)
  
  # CHANNELS/MARKERS DATA.FRAME
  pd <- data.frame(
    "channel" = chans,
    "marker" = LAPPLY(chans, function(z){
      if(z %in% names(marks)) {
        return(marks[z])
      } else {
        return(NA)
      }
    })
  )
  
  # FILE MISSING
  if (is.null(file)) {
    # FILE EXISTS
    if (length(grep("Experiment-Markers.csv", list.files())) != 0) {
      message("Experiment-Markers.csv found in working directory.")
      # MULTIPLE FILES? - CHECK MATCHING CHANNELS
      found_files <- list.files()[grep("Experiment-Markers.csv", list.files())]
      dt <- lapply(found_files, function(z) {
        mrks <- read_from_csv(z)
        rownames(mrks) <- NULL
        # CHANNELS MUST MATCH
        if (all(cyto_channels(x) %in% mrks$channel)) {
          return(mrks)
        } else {
          return(NULL)
        }
      })
      names(dt) <- found_files
      # FILES FOUND WITHOUT CHANNEL MATCH
      if (all(LAPPLY(dt, "is.null"))) {
        # CHANNEL/MARKER DATA.FRAME
        dt <- pd
      # FILES FOUND WITH CHANNEL MATCH
      } else {
        # REMOVE EMPTY ENTRIES - SELECT FIRST ONE
        dt[LAPPLY(dt, "is.null")] <- NULL
        file <- names(dt)[1]
        dt <- dt[[1]]
      }
    # FILE DOESN'T EXIST
    } else {
      # CHANNEL/MARKER DATA.FRAME
      dt <- pd
    }
  # FILE SUPPLIED
  } else {
    # FILE EXISTS
    if (file_exists(file)) {
      message(
        paste0("Editing data in ", file, "...")
      )
      dt <- read_from_csv(file)
      # DT MUST CONTAIN ALL CHANNELS IN DATA
      if(any(!pd$channel %in% dt$channel)) {
        dt <- rbind(
          dt,
          pd[which(!pd$channel %in% dt$channel), , drop = FALSE]
        )
      }
    # FILE DOESN'T EXIST
    } else {
      # CHANNEL/MARKER DATA.FRAME
      dt <- pd
    }
  }
  
  # DEFAULT FILE NAME
  if (is.null(file)) {
    file <- paste0(format(Sys.Date(), "%d%m%y"), "-Experiment-Markers.csv")
  }
  
  # EDIT
  if(interactive()) {
    dt_edit <- data_edit(dt,
                         logo = CytoExploreR_logo(),
                         title = "Experiment Markers Editor",
                         row_edit = FALSE, # cannot add/remove rows
                         col_edit = FALSE,
                         col_names = c("channel", "marker"),
                         quiet = TRUE,
                         hide = TRUE,
                         viewer = "pane",
                         ...)
    # WRITE_TO_CSV NOT IN SCOPE
    write_to_csv(dt_edit,
                 file,
                 row.names = FALSE)
  }
  
  # ONLY UPDATE CHANNELS/MARKERS RELEVANT TO DATA
  ind <- LAPPLY(dt$channel, match, chans)
  cyto_chans<- dt_edit$channel[ind]
  cyto_marks <- dt_edit$marker[ind]
  names(cyto_marks) <- cyto_chans
  
  # UPDATE CHANNELS
  if(any(!cyto_chans %in% chans)) {
    cyto_channels(x) <- cyto_chans
  }
  
  # EMPTY -> NA
  ind <- which(LAPPLY(cyto_marks, ".empty"))
  if (length(ind) > 0) {
    cyto_marks[ind] <- NA
  }
  
  # UPDATE MARKERS
  if (!.all_na(cyto_marks)) {
    if (cyto_class(x, c("flowFrame", "flowSet"), TRUE)){
      flowCore::markernames(x) <- cyto_marks
    } else {
      flowWorkspace::markernames(x) <- cyto_marks
    }
  }
  
  # UPDATE CHANNELS/MARKERS IN PLACE
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
      
      # Run through each file and check names match samples
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
      } else {
        
        # Remove NULL entries from list - result should be of length 1
        pd[LAPPLY(pd, "is.null")] <- NULL
        file <- names(pd)[1]
        pd <- pd[[1]]
      }
    } else {
      
      # Extract cyto_details
      pd <- cyto_details(x)
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
    }
  }
  
  # FILE
  if (is.null(file)) {
    file <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-Experiment-Details.csv"
    )
  }
  
  # Edit cyto_details - rownames cannot be edited
  if(interactive()) {
    # REMOVE ROW NAMES
    cyto_names <- rownames(pd)
    rownames(pd) <- NULL
    # EDIT
    pd <- data_edit(pd,
                    logo = CytoExploreR_logo(),
                    title = "Experiment Details Editor",
                    row_edit = FALSE,
                    # col_readonly = "name",
                    quiet = TRUE,
                    hide = TRUE,
                    viewer = "pane",
                    ...)
    # REPLACE ROW NAMES
    rownames(pd) <- cyto_names
  }
  
  # Update cyto_details
  cyto_details(x) <- pd
  
  # Save - cannot save above as rownames have been removed
  write_to_csv(pd,
               file = file,
               row.names = TRUE)
  
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
  fs <- cyto_data_extract(x, parent = "root")[["root"]]
  
  # Spillover matrix supplied - matrix, data.frame or csv file
  if (!is.null(spillover)) {
    # spillover is a character string containing name of csv file
    if (cyto_class(spillover, "character", TRUE)) {
      # read in spillover matrix
      spill <- read_from_csv(spillover)
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

#' Names of gated populations in gatingTemplate/GatingHierarchy/GatingSet
#'
#' \code{cyto_nodes} is simply an autocomplete-friendly wrapper for
#' \code{\link[flowWorkspace:gs_get_pop_paths]{gs_get_pop_paths}}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingSet}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}},
#'   \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}, or the name of
#'   a gatingTemplate CSV file.
#' @param path can be either \code{"full"} or \code{"auto"} to return the
#'   complete path or shortest unique path to each node.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:gs_get_pop_paths]{gs_get_pop_paths}} or
#'   \code{\link[openCyto:gt_get_nodes]{gt_get_nodes}}.
#'
#' @return character vector of gated node/population names.
#'
#' @importFrom flowWorkspace gh_get_pop_paths gs_get_pop_paths
#' @importFrom openCyto gatingTemplate gt_get_nodes
#'
#' @export
cyto_nodes <- function(x,
                       path = "full",
                       ...) {
  
  # GATINGHIERARCHY
  if (cyto_class(x, "GatingHierarchy", TRUE)) {
    gh_get_pop_paths(x,
                     path = path,
                     ...)
  # GATINGSET
  } else if (cyto_class(x, "GatingSet", TRUE)) {
    gs_get_pop_paths(x,
                     path = path,
                     ...)
  # GATINGTEMPLATE
  } else {
    # GATINGTEMPLATE NAME
    if(is.character(x)) {
      x <- file_ext_append(x, ".csv")
      x <- suppressMessages(gatingTemplate(x))
    }
    # NODES
    if(path == "full") {
      names(gt_get_nodes(x, only.names = TRUE, ...))
    } else {
      unname(gt_get_nodes(x, only.names = TRUE, ...))
    }
  }
}

## CYTO_NODES_CHECK ------------------------------------------------------------

#' Check if nodes are unique and exist in GatingSet or GatingHierarchy
#'
#' @param x object of class \code{GatingHierarchy}, \code{GatingSet},
#'   \code{gatingTemplate} or the name of a gatingTemplate CSV file.
#' @param nodes vector of node paths to check.
#'
#' @return supplied nodes or throw an error if any nodes do not exist or are not
#'   unique within the \code{GatingHierarchy}, \code{GatingSet} or
#'   \code{gatingTemplate}.
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
#' @param x object of class \code{GatingHierarchy}, \code{GatingSet},
#'   \code{gatingTemplate} or the name of a gatingTemplate CSV file.
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
#' @param x object of class \code{GatingHierarchy},\code{GatingSet},
#'   \code{gatingTemplate} or the name of a gatingTemplate CSV file.
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
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return list of spillover matrices or NULL.
#'
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace gh_get_compensations gs_get_compensations
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
#' # Extract spillover matrices
#' spill <- cyto_spillover_extract(gs)
#' @export
cyto_spillover_extract <- function(x) {
  
  # GATINGSET
  if (cyto_class(x, "GatingSet", TRUE)) {
    spill <- gs_get_compensations(x)
    if (all(LAPPLY(spill, "is.null"))) {
      spill <- NULL
    } else {
      spill <- lapply(spill, function(z) {
        z@spillover
      })
    }
    # GATINGHIERARCHY
  } else if (cyto_class(x, "GatingHierarchy", TRUE)) {
    spill <- gh_get_compensations(x)
    if (!is.null(spill)) {
      spill <- list(spill@spillover)
      names(spill) <- cyto_names(x)
    }
    # CYTOSET
  } else if (cyto_class(x, "flowSet")) {
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
    # CYTOFRAME
  } else if (cyto_class(x, "flowFrame")) {
    # CANNOT SET NAMES
    spill <- tryCatch(keyword(x, "SPILL"), error = function(e) {
      NULL
    })
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
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} to use for the
#'   calibration. For the best calibration is recommended that users supply
#'   samples containing both negative and positive events in each channel.
#' @param parent name of the parent population to use for channel calibration
#'   when a \code{GatingHierarchy} or \code{GatingSet} is supplied, set to the
#'   \code{"root"} node by default.
#' @param type indicates the type of calibration to perform, options include
#'   \code{"range"} or \code{"quantile"}, set to \code{"quantile"} by default.
#'   Range calibration simply uses the full range of values across samples for
#'   the calibration. Quantile calibration computes an lower and upper quantile
#'   for each channel, values falling outside the calibration range are assigned
#'   the bottom or top colour.
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
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Calibration
#' cyto_calibrate(gs)
#'
#' # Colour based on Hoechst-405 staining
#' cyto_plot(gs[1],
#' parent = "root",
#' channels = c("FSC-A", "SSC-A"),
#' point_col = "Hoechst-405")
#'
#' @export
cyto_calibrate <- function(x,
                           parent = "root",
                           type = "quantile",
                           probs = c(0.01, 0.95),
                           ...){
  
  # COMPUTE CHANNEL RANGES
  if(grepl("^r", type, ignore.case = TRUE)){
    cyto_cal <- cyto_apply(x,
                           "cyto_stat_range",
                           parent = parent,
                           input = "matrix",
                           inverse = FALSE,
                           copy = FALSE,
                           ...)
    # COMPUTE CHANNEL QUANTILES
  }else if(grepl("^q", type, ignore.case = TRUE)){
    if(length(probs) != 2) {
      stop("A lower and upper quantile probability must be supplied.")
    }
    cyto_cal <- cyto_apply(x,
                           "cyto_stat_quantile",
                           parent = parent,
                           input = "matrix",
                           inverse = FALSE,
                           copy = FALSE,
                           probs = probs,
                           ...)
  }
  
  # MIN/MAX
  cyto_cal <- apply(cyto_cal, 2, function(x){
    c("min" = min(x), "max" = max(x))
  })
  
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

# flowFrame/flowSet methods require parent argument (not used) otherwise parent
# gets passed to FUN. cyto_apply cannot dispatch to cyto_stats functions unless
# called within the package - users should interact with cyto_stats_compute
# directly.

#' Apply a function to elements of cytoset, GatingHierarchy or GatingSet
#'
#' \code{cyto_apply} is convenient wrapper around \code{lapply} and \code{apply}
#' to apply a function over elements of a \code{cytoset}, \code{GatingHierarchy}
#' or \code{GatingSet}. \code{cyto_apply} is extremely flexible by supporting
#' functions that accept the data in either \code{cytoset}, \code{cytoframe},
#' \code{matrix} or \code{vector} formats. All the data processing steps are
#' handled internally by \code{\link{cyto_data_extract}} prior to passing the
#' data to the specified function. It is important that the arguments in
#' \code{FUN} do not conflict with the arguments of \code{cyto_apply} and they
#' should be supplied to \code{cyto_apply} by name.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param FUN name of a function to apply to each element of \code{x}.
#' @param ... additional arguments passed to \code{FUN}. Multiple arguments are
#'   supported but must be named, see examples below.
#' @param simplify logical indicating whether attempts should be made to coerce
#'   the output to a \code{cytoset} or \code{matrix}, set to TRUE by default.
#' @param input indicates the data input format as required by \code{FUN} can be
#'   either 1 - "cytoset", 2 - "cytoframe", 3 - "matrix", 4 - "column" or 5 -
#'   "row", set to "cytoframe" by default. \code{cyto_apply} will take care of
#'   all the data formatting prior to passing it \code{FUN}. The \code{"column"}
#'   and \code{"row"} options are for functions that expect vectors as the
#'   input.
#' @param parent name of the parent(s) population to extract from
#'   \code{GatingHierarchy} or \code{GatingSet} objects, set to \code{"root"} by
#'   default.
#' @param copy logical indicating whether the data should be copied prior to
#'   preprocessing the data, set to TRUE by default to ensure that the original
#'   data remains unchanged.
#' @param channels vector of channels which should be included in the data
#'   passed to \code{FUN}, set to all channels by default.
#' @param trans object of class \code{transformerList} containing the
#'   definitions of the transformers applied to the supplied data. These
#'   transformers are only required when a cytoset is supplied to inverse
#'   transformations when \code{inverse = TRUE}.
#' @param inverse logical indicating whether the data should be inverse
#'   transformed prior to applying \code{FUN}, set to FALSE by default.
#' @param slot name of the slot to extract from the output if the object
#'   returned by the function is a list or S4 R object, set to NULL by default.
#'   the output of the function
#'
#' @importFrom methods is
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace cytoset flowFrame_to_cytoframe
#'   cytoframe_to_flowFrame flowSet_to_cytoset
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' # Compute CD44 & CD69 quantiles for CD4 T Cells & CD8 T Cells
#' cyto_apply(gs,
#' "quantile",
#' probs = c(0.25, 0.5, 0.75),
#' na.rm = TRUE,
#' input = "column",
#' parent = c("CD4 T Cells", "CD8 T Cells"),
#' channels = c("CD44", "CD69"))
#'
#' @name cyto_apply
NULL

#' @noRd
#' @export
cyto_apply <- function(x,
                       ...) {
  UseMethod("cyto_apply")
}

#' @export
cyto_apply.default <- function(x,
                               FUN,
                               ...,
                               simplify = TRUE,
                               parent = "root",
                               input = "cytoframe",
                               copy = TRUE,
                               channels = NULL,
                               trans = NA,
                               inverse = FALSE,
                               slot = NULL) {
  
  # GATINGHIERARCHY/GATINGSET
  if(cyto_class(x,  "GatingSet")) {
    
    # TRANSFORMERS
    trans <- cyto_transformers_extract(x)
    
    # LIST OF CYTOSETS
    x <- cyto_data_extract(x,
                           parent = parent)
    
    # APPLY FUNCTION
    res <- structure(
      lapply(x, function(z){
        cyto_apply(z,
                   FUN = FUN,
                   simplify = simplify,
                   input = input,
                   copy = copy,
                   channels = channels,
                   trans = trans,
                   inverse = inverse,
                   slot = slot,
                   ...)
      }), names = names(x)
    )
    
    # MATRIX - SINGLE PARENT
    if(length(parent) == 1) {
      res <- res[[1]]
    }
    
    # RETURN
    return(res)
    
    # UNSUPPORTED METHOD  
  } else {
    stop(
      paste0(
        "'cyto_apply' does not accept objects of class ", 
        cyto_class(x, class = TRUE), "!"
      )
    )
  }
  
}

#' @export
cyto_apply.flowSet <- function(x,
                               FUN,
                               ...,
                               simplify = TRUE,
                               parent = "root",
                               input = "cytoframe",
                               copy = TRUE,
                               channels = NULL,
                               trans = NA,
                               inverse = FALSE,
                               slot = NULL) {
  
  # FUNCTION
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
  FUN <- match_fun(FUN) # namespaced function character covered
  
  # COPY
  if(copy) {
    x <- cyto_copy(x)
  }
  
  # CHANNELS
  if(!is.null(channels)) {
    x <- x[, cyto_channels_extract(x, channels), drop = FALSE]
  }
  
  # INVERSE TRANSFORM
  if(!.all_na(trans) & inverse) {
    x <- cyto_transform(x,
                        trans = trans,
                        inverse = inverse,
                        plot = FALSE,
                        quiet = TRUE)
  }
  
  # CYTOSET INPUT 
  if(input == 1 | 
     grepl("^cytoset", input, ignore.case = TRUE) |
     grepl("^cs", input, ignore.case = TRUE)) {
    input <- "cytoset"
    # CYTOFRAME INPUT
  } else if(input == 2 | 
            grepl("^cytoframe", input, ignore.case = TRUE) |
            grepl("^cf", input, ignore.case = TRUE)) {
    input  <- "cytoframe"
    # MATRIX INPUT
  } else if(input == 3 | 
            grepl("^m", input, ignore.case = TRUE)) {
    input <- "matrix"
    # COLUMN/CHANNEL INPUT
  } else if(input == 4 |
            grepl("^co", input, ignore.case = TRUE) |
            grepl("^ch", input, ignore.case = TRUE)) {
    input <- "column"
    # ROW/CELL
  } else if(input == 5 |
            grepl("^r", input, ignore.case = TRUE) |
            grepl("ce", input, ignore.case = TRUE)) {
    input <- "row"
  }
  
  # DISPATCH
  if(input == "cytoset") {
    res <- structure(
      lapply(seq_along(x), function(z){
        output <- FUN(x[z], ...)
        # SLOT
        if(!is.null(slot)) {
          output <- output[[slot]]
        }
        return(cyto_convert(output))
      }), names = cyto_names(x))
  } else if(input == "flowFrame") { # internal use only
    res <- structure(
      lapply(cyto_names(x), function(z){
        output <- FUN(cytoframe_to_flowFrame(x[[z]]), ...)
        # SLOT
        if(!is.null(slot)) {
          output <- output[[slot]]
        }
        return(cyto_convert(output))
      }), names = cyto_names(x))
  } else if(input == "cytoframe") {
    res <- structure(
      lapply(cyto_names(x), function(z){
        output <-FUN(x[[z]], ...)
        # SLOT 
        if(!is.null(slot)) {
          output <- output[[slot]]
        }
        return(cyto_convert(output))
      }), names = cyto_names(x))
  } else if(input == "matrix") {
    res <- structure(
      lapply(cyto_names(x), function(z){
        output <- FUN(exprs(x[[z]]), ...)
        # SLOT
        if(!is.null(slot)) {
          output <- output[[slot]]
        }
        return(cyto_convert(output))
      }), names = cyto_names(x))
  } else if(input == "column") {
    # TODO: Add support for passing channel-specific arguments through here
    # named vector or named list
    res <- structure(
      lapply(cyto_names(x), function(z){
        output <- apply(exprs(x[[z]]), 2, FUN, ...)
        # SLOT
        if(!is.null(slot)) {
          output <- output[[slot]]
        }
        return(cyto_convert(output))
      }), names = cyto_names(x))
  } else if(input =="row") {
    res <- structure(
      lapply(cyto_names(x), function(z){
        output <- apply(exprs(x[[z]]), 1, FUN, ...)
        # SLOT
        if(!is.null(slot)) {
          output <- output[[slot]]
        }
        return(cyto_convert(output))
      }), names = cyto_names(x))
  }
  
  # SIMPLIFY
  if(simplify) {
    # CYTOFRAMES -> CYTOSET
    if(cyto_class(res[[1]], "flowFrame")) {
      res <- cytoset(res)
      # CYTOSET
    }  else if(cyto_class(res[[1]], "flowSet") &
               length(res[[1]]) == 1) {
      res <- cytoset(
        structure(
          lapply(res, `[[`, 1), 
          names = names(res)
        )
      )
      # VECTOR/MATRIX/LIST
    } else {
      # VECTORS TO MATRIX
      if(is.null(dim(res[[1]]))) {
        res <- structure(
          lapply(names(res), function(z){
            matrix(res[[z]],
                   nrow = 1,
                   dimnames = list(z,
                                   names(res[[z]])))
          }), names = names(res))
      }
      # PREPARE MATRICES
      res <- lapply(names(res), function(z){
        # ROWNAMES
        if(is.null(rownames(res[[z]]))) {
          if(nrow(res[[z]]) > 1) {
            rownames(res[[z]]) <- paste(z, 1:nrow(res[[z]]), sep = "|")
          } else {
            rownames(res[[z]]) <- z
          }
        } else {
          if(!all(rownames(res[[z]]) == z)) {
            rownames(res[[z]]) <- paste(z, rownames(res[[z]]), sep = "|")
          }
        }
        # COLNAMES
        if(is.null(colnames(res[[z]]))) {
          if(ncol(res[[z]]) == 1) {
            colnames(res[[z]]) <- FUN_NAME
          } else {
            colnames(res[[z]]) <- paste0(FUN_NAME, "-", 1:nrow(res[[z]]))
          }
        }
        return(res[[z]])
      })
      # BIND MATRICES
      res <- do.call("rbind", res)
    }
  }
  
  # OUTPUT
  return(res)
  
}

#' @rdname cyto_apply
#' @export
cyto_apply.list <- function(x,
                            FUN,
                            ...,
                            simplify = TRUE,
                            parent = "root",
                            input = "cytoframe",
                            copy = TRUE,
                            channels = NULL,
                            trans = NA,
                            inverse = FALSE,
                            slot = NULL) {
  
  structure(
    lapply(x, function(z){
      cyto_apply(z,
                 FUN = FUN,
                 ...,
                 simplify = simplify,
                 parent = parent,
                 input = input,
                 copy = copy,
                 channels = channels,
                 trans = trans,
                 inverse = inverse,
                 slot = slot)
    }), names = names(x)
  )
  
}

## CYTO_CONVERT ----------------------------------------------------------------

#' Convert flowFrame/flowSet objects to cytoframe/cytoset objects
#'
#' @param x object of class \code{flowFrame}, \code{flowSet}, \code{cytoframe}
#'   or \code{cytoset}.
#' @param ... additional arguments passed to
#'
#' @return a cytoframe or cytoset.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace flowFrame_to_cytoframe flowSet_to_cytoset
#'
#' @seealso
#'   \code{\link[flowWorkspace:flowFrame_to_cytoframe]{flowFrame_to_cytoframe}}
#' @seealso
#'   \code{\link[flowWorkspace:flowSet_to_cytoset]{flowSet_to_cytoset}}  
#'
#' @examples 
#' library(CytoExploreRData)
#' 
#' # flowSet to cytoset
#' cs <- cyto_convert(Activation)
#'
#' @export
cyto_convert <- function(x,
                         ...) {

  if(cyto_class(x, "flowFrame", TRUE)) {
    x <- flowFrame_to_cytoframe(x, ...)
  } else if(cyto_class(x, "flowSet", TRUE)) {
    x <- flowSet_to_cytoset(x)
  }
  return(x)
  
}

## CYTO_CBIND ------------------------------------------------------------------

#' Bind new columns to cytoframe or cytoset
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}.
#' @param cols matrix of columns to be added to \code{x} can be supplied in a
#'   list for \code{cytoset} method.
#'
#' @importFrom flowWorkspace cf_append_cols realize_view cf_is_subsetted
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
    stop("'cols' must be a matrix!")
  }
  
  # HANDLE EMPTY FLOWFRAME/CYTOFRAME
  if(cyto_stat_count(x) == 0 & nrow(cols) > 0) {
    cols <- cols[-seq_len(nrow(cols)), ,drop = FALSE]
  }
  
  # CYTOFRAME
  if(cyto_class(x, "cytoframe", TRUE)){
    # REALIZE VIEW
    if(cf_is_subsetted(x)) {
      x <- realize_view(x)
    }
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
                              input = "matrix",
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
          cols[1:cyto_counts[z], , drop = FALSE]
        }else{
          start <- sum(cyto_counts[1:(z-1)]) + 1
          end <- start + cyto_counts[z] - 1
          cols[start:end, , drop = FALSE]
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
  if(cyto_class(cf_list[[1]], "cytoframe", TRUE)){
    return(cytoset(cf_list))
  }else{
    return(flowSet(cf_list))
  }
  
}

## CYTO_ENQUIRE ----------------------------------------------------------------

#' Ask CytoExploreR user a question
#'
#' @param x question to ask user.
#' @param options possible options for user response to be TRUE.
#' @param ignore.case logical indicating whether to ignore case when checking
#'   against options, set to TRUE by default.
#' @param ... additional arguments passed to \code{\link[base:grepl]{grepl}}.
#'
#' @return the response or logical if options are specified.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
cyto_enquire <- function(x,
                         options = NULL,
                         ignore.case = TRUE, 
                         ...) {
  # QUESTION
  answer <- readline(x)
  # OPTIONS
  if(!is.null(options)) {
    if(any(LAPPLY(options, 
                  function(z) {
                    grepl(z, 
                          answer, 
                          ignore.case = ignore.case,
                          ...)
                    }))) {
      answer <- TRUE
    } else {
      answer <- FALSE
    }
  }
  return(answer)
}

## CYTO_REQUIRE ----------------------------------------------------------------

#' Load external packages for use within CytoExploreR
#'
#' @param x name of the required package.
#' @param source where to install the package from if it is not found on the
#'   user's machine, options include \code{CRAN}, \code{BioC} or \code{GitHub}.
#'   Set to \code{CRAN} by default.
#' @param repo name of the GitHub repository from which the package should be
#'   installed when \code{source = GitHub}.
#' @param version minimal version requirement if the package is located on the
#'   user's machine.
#' @param ref citation to print to the console when loading the package for use
#'   within CytoExploreR.
#' @param ... additional arguments passed to \code{install.packages},
#'   \code{BiocManager::install} or \code{remotes::install_github}.
#'
#' @return NULL
#'
#' @importFrom utils install.packages installed.packages
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_require <- function(x,
                         source = "CRAN",
                         repo = NULL,
                         version = NULL,
                         ref = NULL,
                         ...) {
  
  # PACKAGE LOCATED
  pkgs <- as.data.frame(installed.packages())
  
  # PACKAGE FOUND
  if(x %in% pkgs$Package) {
    # CHECK VERSION
    if(!is.null(version)) {
      # PACKAGE VERSION
      pkg_version <- pkgs[pkgs$Package == x, "Version"]
      # VERSION TOO OLD
      if(tryCatch(pkg_version < version, error = function(e){TRUE})) {
        inst <- TRUE
      # VERSION CORRECT
      } else {
        inst <- FALSE
      }
    # VERSION DOES NOT MATTER
    } else {
      inst <- FALSE
    }
  # PACKAGE MISSING 
  } else {
    inst <- TRUE
  }
  
  # INSTALL PACKAGE
  if(inst) {
    # MESSAGE
    message(
      paste0("Installing required package: ",
             x, "...")
    )
    # GITHUB - REPO
    if(!is.null(repo)) {
      if(!"remotes" %in% pkgs$Package) {
        install.packages("remotes")
      }
      if(requireNamespace("remotes")) {
        remotes::install_github(repo, ...)
      }
    # CRAN
    } else if(grepl("^c", source, ignore.case = TRUE)){
      install.packages(x)
    # BIOCONDUCTOR
    } else if (grepl("^b", source, ignore.case = TRUE)) {
      if(!"BiocManager" %in% pkgs$Package) {
        install.packages("BiocManager", ...)
      }
      if(requireNamespace("BiocManager")) {
        BiocManager::install(x, ...)
      }
    }
  }
  
  # LOAD PACKAGE
  requireNamespace(x)
  
  # REFERENCE
  if(!is.null(ref)) {
    message(ref)
  }
  
}
