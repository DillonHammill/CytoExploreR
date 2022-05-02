## CYTO_IMPORT -----------------------------------------------------------------

#' Import data from other platforms using CytoML
#'
#'
#' @param path path to xml or wsp file and associated fcs files.
#' @param type indicates the type of data to import, options include
#'   \code{"flowJo"}, \code{"cytobank"} or \code{"diva"}. Set to \code{"flowJo"}
#'   by default.
#' @param markers logical indicating whether a call should be made to
#'   \code{cyto_markers_edit} to update the markers associated with channels in
#'   the loaded samples, set to TRUE by default. The name of the csv to which
#'   these details will be supplied can also be passed to this argument.
#' @param details logical indicating whether a call should be made to
#'   \code{cyto_details_edit} to update the experimental details associated with
#'   the loaded samples, set to TRUE by default. The name of the csv to which
#'   these details will be supplied can also be passed to this argument.
#' @param gatingTemplate passed to \code{cyto_gatingTemplate_generate()} to
#'   create a CytoExploreR friendly gatingTemplate based on gates applied to the
#'   loaded samples.
#' @param ... additional arguments passed to \code{cytobank_to_gatingset()},
#'   \code{flowjo_to_gatingset} or \code{diva_to_gatingset}. Refer to
#'   documentation in CytoML package for more details.
#'
#' @return a GatingSet and write new gatingTemplate to CSV file.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#' # flowJo
#' gs <- cyto_import(
#'   "flowjo-samples",
#'   gatingTemplate = "gatingTemplate.csv"
#' )
#' }
#'
#' @export
cyto_import <- function(path = ".",
                        type = "flowJo",
                        markers = TRUE,
                        details = TRUE,
                        gatingTemplate = NULL,
                        ...) {
  
  # TODO: ADD SUPPORT FOR SELECT|EXCLUDE|SORT|RESTRICT
  # PRE OR POST GATINGSET?
  
  # REQUIRE CYTOML
  cyto_require(
    "CytoML",
    source = "BioC",
    repo = "RGLab/CytoML",
    ref = paste(
      "Finak G, Jiang W, Gottardo R (2018). CytoML for cross-platform",
      "cytometry data sharing. Cytometry A,",
      "93(12)."
    )
  )
  
  # FILES IN DIRECTORY
  file_paths <- list.files(
    path,
    full.names = TRUE
  )
  
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
      cytobank_exp <- cyto_func_call(
        "CytoML::open_cytobank_experiment",
        list(acs_file)
      )
      gs <- cyto_func_call(
        "CytoML::cytobank_to_gatingset",
        list(cytobank_exp, ...)
      )
    # XML 
    }else if(any(grepl("xml", file_ext, ignore.case = TRUE))){
      xml_file <- file_paths[
        which(
          grepl(
            "xml",
            file_ext,
            ignore.case = TRUE
          )
        )
      ]
      gs <- cyto_func_call(
        "CytoML::cytobank_to_gatingset",
        list(xml_file, fcs_files)
      )
    }
  # IMPORT DIVA TO GATINGSET
  }else if(grepl("diva", type, ignore.case = TRUE)){
    xml_file <- file_paths[
      which(
        grepl(
          "xml",
          file_ext, 
          ignore.case = TRUE
        )
      )
    ]
    diva_ws <- cyto_func_call(
      "CytoML::open_diva_xml",
      list(xml_file)
    )
    gs <- cyto_func_call(
      "CytoML::diva_to_gatingset",
      list(diva_ws, fcs_files, ...)
    )
  # IMPORT FLOWJO TO GATINGSET
  }else if(grepl("flowjo", type, ignore.case = TRUE)){
    # WSP
    if(any(grepl("wsp", file_ext, ignore.case = TRUE))){
      wsp_file <- file_paths[
        which(
          grepl(
            file_ext, 
            "wsp", 
            ignore.case = TRUE
          )
        )
      ]
      wsp <- open_flowjo_xml(wsp_file)
      gs <- cyto_func_call(
        "CytoML::flowjo_to_gatingset",
        list(wsp, path = path, ...)
      )
    # XML
    }else if(any(grepl("xml", file_ext, ignore.case = TRUE))){
      xml_file <- file_paths[
        which(
          grepl(
            "xml",
            file_ext,
            ignore.case = TRUE
          )
        )
      ]
      flowjo_ws <- cyto_func_call(
        "CytoML::open_flowjo_xml",
        list(xml_file)
      )
      gs <- cyto_func_call(
        "CytoML::flowjo_to_gatingset",
        list(flowjo_ws, fcs_files, ...)
      )
    }
  }
  
  # MARKERS
  if (markers != FALSE) {
    message("Assigning markers to channels...")
    # DEFAULT FILE NAME
    if (markers == TRUE) {
      gs <- cyto_markers_edit(
        gs
      )
    } else {
      gs <- cyto_markers_edit(
        gs,
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
      gs <- cyto_details_edit(
        gs,
        file = details
      )
    }
  }
  
  # BARCODE
  gs <- cyto_barcode(
    gs, 
    type = "events"
  )
  
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
      CytoML::gatingset_to_cytobank(
        x, 
        save_as,
        ...
      )
    } else if (file_ext(save_as) == "wsp") {
      message("Saving GatingSet to flowJo workspace file...")
      CytoML::gatingset_to_flowjo(
        x, 
        save_as,
        ...
      )
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
#'   the name of folder in current working directory, but paths to individual
#'   files are also allowed. Alternatively, users can supply a named list of
#'   matrices to be converted into a cytoset for use within CytoExploreR.
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
#' @param fixed logical passed to \code{grepl()} during file name matching in
#'   \code{select} and \code{exclude} to control whether an exact match is
#'   required, set to FALSE by default for more flexible matching.
#' @param ignore.case logical passed to \code{grepl()} during file name matching
#'   in \code{select} and \code{exclude} to control whether case insensitive
#'   matching is required, set to TRUE by default.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:load_cytoset_from_fcs]{load_cytoset_from_fcs()}}.
#'
#' @return object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowWorkspace load_gs load_cytoset_from_fcs
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
                      restrict = FALSE,
                      fixed = FALSE,
                      ignore.case = TRUE,
                      ...) {

  # PATH - NAMED LIST OF MATRICES
  if(cyto_class(path, "list", TRUE)) {
    # NAMES REQUIRED
    if(is.null(names(path)) | any(is.null(LAPPLY(path, "dim")))) {
      stop(
        "'path' must be a named list of matrices!"
      )
    }
    # CONVERT TO CYTOFRAMES & WRAP AS CYTOSET
    x <- cytoset(
      structure(
        lapply(
          path,
          function(z) {
            # DATA.FRAME -> MATRIX
            if(! cyto_class(z, "matrix", TRUE)) {
              z <- data.matrix(z)
            }
            # CONVERT TO CYTOFRAME
            cyto_convert(
              flowFrame(
                z
              )
            )
          }
        ),
        names = names(path)
      )
    )
    
    # RETURN CYTOSET
    return(x)
  }
  
  # PATH - DIRECTORY/FILES
  files <- LAPPLY(
    path, 
    function(z){
      # DIRECTORY
      if(dir.exists(z)){
        return(
          list.files(
            z,
            full.names = TRUE,
            recursive = TRUE
          )
        )
        # FILE
      } else if(file.exists(z)) {
        return(z)
        # DIRECTORY/ FILE DOES NOT EXIST
      } else {
        stop(
          paste0(
            z,
            " is inaccessible or does not exist!"
          )
        )
      }
    }
  )
  
  # SAVED GATINGSET
  if ("pb" %in% file_ext(files)) {
    # LOAD GATINGSET
    x <- load_gs(
      path = path,
      backend_readonly = FALSE
    )
    # BARCODE
    if(!any(.grepl("^Event-?ID$", cyto_channels(x), ignore.case = TRUE))) {
      x <- cyto_barcode(
        x,
        "events"
      )
    }
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
      lapply(
        select, 
        function(z) {
          file_ind <<- c(
            file_ind,
            .grep(
              z, 
              files, 
              ignore.case = ignore.case, 
              fixed = fixed
            )
          )
        }
      )
      if(length(file_ind) > 0){
        files <- files[unique(file_ind)]
      }
    }
    
    # EXCLUDE
    if (!is.null(exclude)) {
      file_ind <- c()
      lapply(
        exclude, 
        function(z) {
          file_ind <<- c(
            file_ind,
            .grep(
              z,
              files,
              ignore.case = ignore.case,
              fixed = fixed
            )
          )
        }
      )
      if(length(file_ind) > 0){
        files <- files[-unique(file_ind)]
      }
    }
    
    # SORTED FILE PATHS
    if (sort & length(files) > 1) {
      files <- file_sort(files)
    }
    
    # CYTOSET
    x <- load_cytoset_from_fcs(
      files = normalizePath(files),
      ...
    )
    
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

#' Apply anomaly detection algorithms to obtain clean cytometry data
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{CatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}. The \code{root}
#'   node extracted when a \code{GatingSet} or \code{GatingHierachy} is
#'   supplied.
#' @param type the method to use when cleaning the data, options include
#'   \code{"flowAI"}, \code{"flowClean"}, \code{"flowCut"} or \code{"PeacoQC"},
#'   set to \code{"flowAI"} by default.
#' @param ... additional arguments passed to \code{flowAI::flow_auto_qc},
#'   \code{flowClean::flowClean}, \code{flowCut::flowCut} or
#'   \code{PeacoQC::PeacoQC}.
#'
#' @importFrom flowWorkspace gs_cyto_data flowSet_to_cytoset GatingSet
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
#' @references Monaco G, et al. (2016) flowAI: automatic and interactive anomaly
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
#' @references Emmaneel A, et al. (2021) PeacoQC: peak-based selection of high
#'   quality cytometry data. Cytometry A.
#'   \url{https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.24501}
#'
#' @seealso \code{\link[flowAI:flow_auto_qc]{flow_auto_qc}}
#'
#' @export
cyto_clean <- function(x,
                       type = "flowAI",
                       ...) {
  
  # EXTRACT DATA
  if(cyto_class(x, "GatingSet")) {
    cyto_data <- cyto_data_extract(
      x,
      parent = "root",
      copy = FALSE
    )[["root"]]
  } else {
    cyto_data <- x
  }
  
  # FLOWAI
  if(.grepl("ai", type)) {
    cyto_require(
      "flowAI",
      source = "BioC",
      repo = "giannimonaco/flowAI",
      ref = paste0("Monaco G, Chen H, Poidinger M, Chen J,",
                   " de Magalhaes J, Larbi A (2016). flowAI:",
                   " automatic and interactive anomaly discerning",
                   " tools for flow cytometry data. Bioinformatics,",
                   " 32(16).")
    )
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
      )
    )
  # FLOWCLEAN
  } else if(.grepl("clean", type)) {
    cyto_require(
      "flowClean",
      source = "BioC",
      repo = "cafletezbrant/flowClean",
      ref = paste0(
        "Fletez-Brant K, Spidlen J, Brinkman R, Roederer M,",
        " Chattopadhyay P (2016). flowClean: Automated",
        " identification and removal of fluorescence anaomalies",
        " in flow cytometry data. Cytometry A 89(5)"
      )
    )
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
  } else if(.grepl("cut", type)) {
    cyto_require(
      "flowCut",
      source = "BioC",
      repo = "jmeskas/flowCut",
      ref = paste0(
        "Meskas J, Wang S, Brinkman R (2021). flowCut --- An R",
        " Package for precise and accurate automated removal of",
        " outlier events and flagging of files based on time",
        " versus fluorescence analysis. bioRxiv"
      )
    )
    cyto_data <- cyto_apply(
      cyto_data,
      "flowCut::flowCut",
      input = "flowFrame",
      slot = "frame",
      copy = FALSE,
      Plot = "None", 
      ...
    )
  # PEACOQC
  } else {
    cyto_require(
      "PeacoQC",
      source = "BioC",
      repo = "saeyslab/PeacoQC",
      ref = paste0(
        "Emmaneel A, et al. (2021) PeacoQC: peak-based selection of high ",
        "quality cytometry data. Cytometry A."
      )
    )
    # PEACOQC
    cyto_data <- cyto_apply(
      cyto_data,
      function(fr) {
        # PREPARE ARGUMENTS
        args <- list(...)
        # TURN OFF PLOTTING
        if(!"plot" %in% names(args)) {
          args[["plot"]] <- FALSE
        }
        # DON'T WRITE FCS FILES
        if(!"save_fcs" %in% names(args)) {
          args[["save_fcs"]] <- FALSE
        }
        # DEFAULT CHANNELS
        if(!"channels" %in% names(args)) {
          args[["channels"]] <- cyto_channels(
            fr, 
            exclude = c("Event", "Sample", "Time")
          )
        }
        # REMOVE MARGINAL EVENTS - RANGES SHOULD BE SUPPLIED MANUALLY
        fr <- cyto_func_execute(
          "PeacoQC::RemoveMargins",
          c(list("ff" = fr), 
            args)
        )
        # RUN PEACOQC
        fr <- cyto_func_execute(
          "PeacoQC::PeacoQC",
          c(list("ff" = fr), 
            args)
        )
        # RETURN HIGH QUALITY FLOWFRAME
        return(fr[["FinalFF"]])
      },
      input = "flowFrame",
      copy = FALSE
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
#'   a folder in current working directory). Alternatively, users can supply a
#'   named list of matrices to be converted into a cytoset for use within
#'   CytoExploreR.
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
#' @param events numeric passed to \code{cyto_sample} to control the number or
#'   proportion of events to retain in each sample, set to 1 by default to keep
#'   all events. Setting \code{events} to 0 will result in downsampling to the
#'   minimum number of events across samples.
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
                       events = 1, 
                       ...) {
  
  # CYTOSET/GATINGSET
  message(
    paste0(
      "Loading ",
      ifelse(
        cyto_class(
          path, 
          "list", 
          TRUE
        ),
        "data",
        "FCS files"
      ), 
      " into a GatingSet..."
    )
  )
  x <- cyto_load(
    path = path, 
    restrict = FALSE, 
    ...
  )
  
  # MARKERS
  if (markers != FALSE) {
    message("Assigning markers to channels...")
    # DEFAULT FILE NAME
    if (markers == TRUE) {
      x <- cyto_markers_edit(x)
    } else {
      x <- cyto_markers_edit(
        x,
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
      x <- cyto_names_parse(
        x,
        split = parse_names
      )
    }
  }
  
  # EXPERIMENT DETAILS
  if (details != FALSE) {
    message("Updating experiment details...")
    if (details == TRUE) {
      x <- cyto_details_edit(x)
    } else {
      x <- cyto_details_edit(
        x,
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
      x <- cyto_channels_restrict(
        x, 
        exclude = restrict
      )
    }
    # CLEAN DATA
    if (clean != FALSE) {
      if(clean == TRUE){
        clean <- "all"
      }
      message("Cleaning data to remove anomalies...")
      x <- cyto_clean(
        x,
        remove_from = clean
      )
    }
    # SAMPLING
    if(all(events != 1)) {
      x <- cyto_sample(
        x, 
        events = events,
        seed = 56
      )
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
    pd <- type.convert(
      pd, 
      as.is = !factor,
      ...
    )
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
#' @importFrom flowWorkspace sampleNames cf_get_uri
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
  # FLOWFRAME -> IDENTIFIER
  if(cyto_class(x, "flowFrame", TRUE)) {
    # IDENTIFIER - CANNOT ACCESS URI - NOT ON DISK
    nm <- identifier(x)
  # CYTOFRAME -> URI
  } else {
    nm <- file_ext_remove(
      basename(
        cf_get_uri(x)
      )
    )
  }
  # COMBINED EVENTS
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
  cyto_names_split <- strsplit(
    cyto_names, 
    split,
    ...
  )
  
  # SKIP
  if(!is.null(exclude)){
    cyto_names_split <- lapply(
      cyto_names_split, 
      function(z){
        z[-exclude]
      }
    )
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
  cyto_names_split <- lapply(
    cyto_names_split,
    function(z){
     if(length(z) < var_length){
        z <- rep(c(z, rep(NA, var_length)), length.out = var_length)
        return(z)
      } else {
        return(z)
      }
    }
  )
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
#' @param copy logical indicating whether the data should be copied prior to
#'   applying the transformations, set to FALSE by default. Setting this
#'   argument to TRUE will leave the supplied data untouched and return a copy
#'   of the data with the transformations applied.
#' @param plot logical indicating whether the result of the transformations
#'   should be plotted using \code{\link{cyto_plot}}.
#' @param popup logical indicating whether plots should be constructed in a
#'   popup window, set to FALSE by default.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param quiet logical which prevents the printing of messages when set to
#'   TRUE.
#' @param ... additional arguments passed to
#'   \code{\link{cyto_transformers_define}}, when no \code{trans} object is
#'   supplied.
#'
#' @return object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet} with transformations applied.
#'
#' @importFrom flowWorkspace gs_cyto_data GatingSet
#' @importFrom flowCore transform transformList
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
                                   copy = FALSE,
                                   plot = TRUE,
                                   popup = FALSE,
                                   axes_limits = "machine",
                                   quiet = FALSE,
                                   ...) {
  
  # TRANSFORMERS MISSING
  if(is.null(trans) | .all_na(trans)) {
    # WARNING TRANSFORMATION OF CYTOFRAME/CYTOSET OBJECTS
    if(!quiet & !cyto_class(x, "GatingSet")) {
      warning(
        paste(
          "Automatically transforming cytoframe/cytoset objects",
          "is not recommended as transformation definitions will be lost."
        )
      )
    }
    # INVERSE TRANSFORM GATINGHIERARCHY/GATINGSET
    if(cyto_class(x, "GatingSet") & inverse) {
      transformer_list <- cyto_transformers_extract(x)
    # TRANSFORMERS FOR CYTOFRAME/CYTOSET/GATINGHIERARCHY/GATINGSET
    } else {
      transformer_list <- cyto_transformers_define(
        x,
        channels = channels,
        parent = parent,
        type = type,
        select = select,
        plot = FALSE, 
        ...
      )
    }
  # TRANS NOT TRANSFORMLIST OR TRANSFORMERLIST
  } else {
    stop(
      "'trans' must be an object of class transformerList!"
    )
  }
  
  # MESSAGE
  if(!quiet)(
    message(
      paste0(
        "Applying ", 
        ifelse(inverse, "inverse ", ""),
        "data transformations..."
      )
    )
  )
  
  # COPY
  if(copy){
    x <- cyto_copy(x)
  }
  
  # TRANSFORM FLOWFRAME OR FLOWSET
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    
    # Extract transformations from transformerList to transformList
    transform_list <- cyto_transform_extract(
      transformer_list,
      inverse = inverse
    )
    
    # APPLY TRANSFORMATIONS
    x <- suppressMessages(
      transform(
        x,
        transform_list
      )
    )
    
    # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (cyto_class(x, "GatingSet")) {
    # NO TRANSFORMERS
    if(!.all_na(transformer_list)) {
      # EXTRACT DATA
      cs <- cyto_data_extract(
        x,
        parent = "root",
        copy = ifelse(copy, FALSE, TRUE) # COPIED ABOVE?
      )[[1]]

      # GATETEMPLATE
      gateTemplate <- cyto_gateTemplate(x)
      
      # INVERSE TRANSFORMATIONS - APPLIED AT CYTOSET LEVEL - NOT ATTACHED
      if(inverse) {
        # TRANSFORMLIST
        transform_list <- cyto_transform_extract(
          transformer_list,
          inverse = inverse
        )
        # APPLY INVERSE TRANSFORMATIONS
        cs <- suppressMessages(
          transform(
            cs,
            transform_list
          )
        )
        # RECONSTRUCT GATINGSET
        x <- GatingSet(cs)
      # DATA TRANSFORMATIONS - APPLIED GATINGSET LEVEL - ATTACHED
      } else {
        # RECONSTRUCT GATINGSET
        x <- GatingSet(cs)
        # APPLY DATA TRANSFORMATIONS
        x <- suppressMessages(
          transform(
            x,
            transformer_list
          )
        )
      }
      # GATINGHIERARCHY - PREPARE GATES
      if(cyto_class(gateTemplate, "gateTemplate")) {
        # GATINGHIERARCHY
        x <- x[[1]]
        # ATTEMPT TO TRANSFORM GATES
        if(length(gateTemplate) > 0) {
          gateTemplate <- structure(
            lapply(
              gateTemplate,
              function(z) {
                if(any(z$channels %in% names(transformer_list))) {
                  z$gate <- tryCatch(
                    cyto_gate_transform(
                      z$gate,
                      trans = transformer_list,
                      inverse = inverse
                    ),
                    error = function(e) {
                      message(
                        paste0(
                          "Transformation of ", cyto_class(z$gate),
                          " objects is not currently supported!"
                        )
                      )
                      return(z$gate)
                    }
                  )
                }
                return(z)
              }
            ),
            names = names(gateTemplate),
            class = "gateTemplate"
          )
        }
      # GATINGSET - PREPARE GATES
      } else {
        # PREPARE NEW GATES
        if(all(LAPPLY(gateTemplate, "length") > 0)) {
          gateTemplate <- structure(
            lapply(
              seq_along(x),
              function(v) {
                structure(
                  lapply(
                    gateTemplate[[v]],
                    function(z) {
                      if(any(z$channels %in% names(transformer_list))) {
                        z$gate <- tryCatch(
                          cyto_gate_transform(
                            z$gate,
                            trans = transformer_list,
                            inverse = inverse
                          ),
                          error = function(e) {
                            message(
                              paste0(
                                "Transformation of ", cyto_class(z$gate),
                                " objects is not currently supported!"
                              )
                            )
                            return(z$gate)
                          }
                        )
                      }
                      return(z)
                    }
                  ),
                  names = names(gateTemplate[[v]]),
                  class = "gateTemplate"
                )
              }
            ),
            names = cyto_names(x)
          )
        }
      }
      # TRANSFER GATES
      x <- cyto_gateTemplate_apply(
        x,
        gateTemplate
      )
    }
  }

  
  # COMPLETE
  if(!quiet){
    message(
      "DONE!"
    )
  }
  
  # VISUALISE TRANSFORMATIONS
  if(plot) {
    # INVERSE
    if(inverse) {
      transformer_list <- NA
    }
    # PLOT DATA TRANSFORMATIONS
    tryCatch(
      cyto_plot_profile(
        x,
        channels = channels,
        select = select,
        axes_trans = transformer_list,
        axes_limits = axes_limits,
        merge_by = "all"
      ),
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
                                         copy = FALSE,
                                         plot = TRUE,
                                         popup = FALSE,
                                         axes_limits = "machine",
                                         quiet = FALSE,
                                         ...) {
  
  # Added for backwards compatibility - flowFrame/flowSet objects only
  if (cyto_class(x, "GatingSet")) {
    stop(
      paste(
        "GatingHierarchy and GatingSet objects require transformerList",
        "objects to apply transformations."
      )
    )
  }
  
  # MESSAGE
  if(!quiet){
    message(
      paste0(
        "Applying ", 
        ifelse(inverse, "inverse ", ""),
        "data transformations..."
      )
    )
  }
  
  # COPY
  if(copy){
    x <- cyto_copy(x)
  }
  
  # TRANSFORM FLOWFRAME OR FLOWSET
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    # VALID TRANSFORMERS ONLY
    if(any(names(trans) %in% cyto_channels(x))) {
      # RESTRICT TRANSFORMERS - SUBSETTED DATA MAY LACK CHANNELS
      trans <- trans@transforms[names(trans) %in% cyto_channels(x)]
      trans <- transformList(
        names(trans),
        lapply(
          trans, 
          function(z){
            z@f
          }
        )
      )
      # Transformations applied as is - allow for inverse transformList
      x <- suppressMessages(
        transform(
          x, 
          trans
        )
      )
    }
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
      cyto_plot_profile(
        x,
        channels = names(trans),
        axes_trans = trans,
        axes_limits = axes_limits,
        merge_by = "all"
      ), 
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
                                           copy = FALSE,
                                           plot = TRUE,
                                           popup = FALSE,
                                           axes_limits = "machine",
                                           quiet = FALSE,
                                           ...) {
  
  # MESSAGE
  if(!quiet){
    message(
      paste0(
        "Applying ", 
        ifelse(inverse, "inverse ", ""),
        "data transformations..."
      )
    )
  }
  
  # COPY
  if(copy) {
    x <- cyto_copy(x)
  }
  
  # TRANSFORM FLOWFRAME OR FLOWSET
  if (cyto_class(x, c("flowFrame", "flowSet"))) {
    # VALID TRANSFORMERS ONLY
    if(any(names(trans) %in% cyto_channels(x))) {
      # RESTRICT TRANSFORMERS - SUBSETTED DATA MAY LACK CHANNELS
      trans <- cyto_transformers_combine(
        trans[names(trans) %in% cyto_channels(x)]
      )
      # TRANSFORMLIST
      transform_list <- cyto_transform_extract(
        trans, 
        inverse = inverse
      )
      # APPLY TRANSFORMATIONS
      x <- suppressMessages(
        transform(
          x, 
          transform_list
        )
      )
    }
  # TRANSFORM GATINGHIERARCHY OR GATINGSET
  } else if (cyto_class(x, "GatingSet")) {
    # VALID TRANSFORMERS ONLY
    if(any(names(trans) %in% cyto_channels(x))) {
      # RESTRICT TRANSFORMERS - SUBSETTED DATA MAY LACK CHANNELS
      trans <- cyto_transformers_combine(
        trans[names(trans) %in% cyto_channels(x)]
      )
      # EXTRACT DATA
      cs <- cyto_data_extract(
        x,
        parent = "root",
        copy = ifelse(copy, FALSE, TRUE) # COPIED ABOVE?
      )[[1]]
      # GATETEMPLATE
      gateTemplate <- cyto_gateTemplate(x)
      # INVERSE TRANSFORMATIONS - APPLIED AT CYTOSET LEVEL - NOT ATTACHED
      if(inverse) {
        # TRANSFORMLIST
        transform_list <- cyto_transform_extract(
          trans,
          inverse = inverse
        )
        # APPLY INVERSE TRANSFORMATIONS
        cs <- suppressMessages(
          transform(
            cs,
            transform_list
          )
        )
        # RECONSTRUCT GATINGSET
        x <- GatingSet(cs)
        # DATA TRANSFORMATIONS - APPLIED GATINGSET LEVEL - ATTACHED
      } else {
        # RECONSTRUCT GATINGSET
        x <- GatingSet(cs)
        # APPLY DATA TRANSFORMATIONS
        x <- suppressMessages(
          transform(
            x,
            trans
          )
        )
      }
      # GATINGHIERARCHY - PREPARE GATES
      if(cyto_class(gateTemplate, "gateTemplate")) {
        # GATINGHIERARCHY
        x <- x[[1]]
        # ATTEMPT TO TRANSFORM GATES
        if(length(gateTemplate) > 0) {
          gateTemplate <- structure(
            lapply(
              gateTemplate,
              function(z) {
                if(any(z$channels %in% names(trans))) {
                  z$gate <- tryCatch(
                    cyto_gate_transform(
                      z$gate,
                      trans = trans,
                      inverse = inverse
                    ),
                    error = function(e) {
                      message(
                        paste0(
                          "Transformation of ", cyto_class(z$gate),
                          " objects is not currently supported!"
                        )
                      )
                      return(z$gate)
                    }
                  )
                }
                return(z)
              }
            ),
            names = names(gateTemplate),
            class = "gateTemplate"
          )
        }
        # GATINGSET - PREPARE GATES
      } else {
        # PREPARE NEW GATES
        if(all(LAPPLY(gateTemplate, "length") > 0)) {
          gateTemplate <- structure(
            lapply(
              seq_along(x),
              function(v) {
                structure(
                  lapply(
                    gateTemplate[[v]],
                    function(z) {
                      if(any(z$channels %in% names(trans))) {
                        z$gate <- tryCatch(
                          cyto_gate_transform(
                            z$gate,
                            trans = trans,
                            inverse = inverse
                          ),
                          error = function(e) {
                            message(
                              paste0(
                                "Transformation of ", cyto_class(z$gate),
                                " objects is not currently supported!"
                              )
                            )
                            return(z$gate)
                          }
                        )
                      }
                      return(z)
                    }
                  ),
                  names = names(gateTemplate[[v]]),
                  class = "gateTemplate"
                )
              }
            ),
            names = cyto_names(x)
          )
        }
      }
      # TRANSFER GATES
      x <- cyto_gateTemplate_apply(
        x,
        gateTemplate
      )
    }
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
      cyto_plot_profile(
        x,
        parent = "root",
        channels = names(trans),
        axes_trans = trans,
        axes_limits = axes_limits,
        merge_by = "all"
      ),
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
#'   should be returned, set to FALSE by default. Users should set this argument
#'   to TRUE if they wish to keep the original data unchanged.
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
#' @param events passed to \code{\link{cyto_sample}} to extract a subset of
#'   events from each sample. \code{events} indicates the number or proportion
#'   of events to extract from each sample, set to \code{1} by default to
#'   extract all events.
#' @param coerce logical to indicate whether the data should be merged into a
#'   single object prior to returning the data, set to FALSE by default. If
#'   \code{coerce} is TRUE, \code{sample} controls the proportion or number of
#'   events to keep in the merged object.
#' @param barcode logical to control whether \code{cyto_barcode} should be
#'   called on the data prior to export. This argument also accepts options that
#'   can be passed to the \code{type} argument of \code{cyto_barcode} to control
#'   the type of barcoding to perform.
#' @param overwrite logical passed to \code{cyto_barcode} to provide a
#'   non-interactive way of controlling how existing barcodes should be handled.
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
#' @importFrom flowWorkspace gs_pop_get_data cf_get_uri cytoframe_to_flowFrame
#'   cytoset_to_flowSet
#' @importFrom flowCore flowSet
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
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_barcode}}
#' @seealso \code{\link{cyto_coerce}}
#' @seealso \code{\link{cyto_sample}}
#' @seealso \code{\link{cyto_transform}}
#'
#' @export
cyto_data_extract <- function(x,
                              parent = "root",
                              select = NULL,
                              copy = FALSE,
                              format = "cytoset",
                              channels = NULL,
                              markers = FALSE,
                              split = FALSE,
                              trans = NA,
                              inverse = FALSE,
                              events = 1,
                              coerce = FALSE,
                              barcode = FALSE,
                              overwrite = NULL,
                              seed = NULL,
                              path = "auto") {
  
  # PARENT - CYTO_STATS_COMPUTE ALIAS NULL
  if(is.null(parent)) {
    parent = "root"
  }
  
  # PARENTAL PATH
  if(cyto_class(x, "GatingSet")) {
    parent <- cyto_nodes_convert(
      x,
      nodes = parent,
      path = path
    )
  }
  
  # EXTRACT TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # SELECT
  if(!is.null(select) & !cyto_class(x, "flowFrame")) {
    x <- cyto_select(x, select)
  }
  
  # CYTOFRAME/CYTOSET LIST
  if(cyto_class(x, c("flowFrame", "flowSet"))) {
    cs_list <- list(
      cyto_convert(x) # SLOW WITH OLD DATA STRUCTURES
    )
  # GATINGHIERARCHY/GATINGSET
  } else {
    cs_list <- structure(
      lapply(
        parent, 
        function(z){
          gs_pop_get_data(x, z)
        }
      ),
      names = parent
    )
  }
  
  # PREPARE CYTOFRAME|CYTOSET LIST
  res <- lapply(
    seq_along(cs_list), 
    function(id) {
      # CYTOFRAME|CYTOSET
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
      # BARCODE
      if(barcode != FALSE) {
        if(barcode == TRUE) {
          barcode <- "samples"
        }
        cs <- cyto_barcode(
          cs,
          type = barcode,
          overwrite = overwrite
        )
      }
      # COERCE - NAME WITH PARENT - BYPASS CYTOFRAME
      if(coerce & cyto_class(cs, "flowSet")) {
        # SAMPLE COERCED CYTOSET
        cs <- cyto_coerce(
          cs,
          format = "cytoset",
          events = events,
          name = ifelse(is.null(names(cs_list)), 
                        paste0("merge-", id),
                        paste0(names(cs_list)[id], "-merge")),
          barcode = FALSE,
          seed = seed
        )
      # SAMPLE
      } else if(.all_na(events) | all(events != 1)) {
        # SAMPLE EACH CYTOFRAME
        cs <- cyto_sample(
          cs,
          events = events,
          seed = seed
        )
      }
      # TRANSFORM
      if(inverse & !.all_na(trans)) {
        cs <- cyto_transform(
          cs,
          trans = trans,
          inverse = inverse,
          plot = FALSE,
          quiet = TRUE
        )
      }
      # FLOWFRAME - INTERNAL BACKWARDS COMPATIBILITY
      if(.grepl("flowFrame", format)) {
        # CYTOFRAME -> FLOWFRAME
        if(cyto_class(cs, "flowFrame")) {
          structure(
            list(
              cytoframe_to_flowFrame(cs)
            ),
            names = cyto_names(cs)
          )
        # CYTOSET -> FLOWFRAME
        } else {
          structure(
            lapply(
              cyto_names(cs), 
              function(z) {
                cs[[z]]
              }
            ), names = cyto_names(cs)
          )
        }
      # FLOWSET - INTERNAL BACKWARDS COMPATIBILITY
      } else if(.grepl("flowSet", format)) {
        # CYTOFRAME -> FLOWSET
        if(cyto_class(cs, "flowFrame")) {
          flowSet(
            cytoframe_to_flowFrame(
              cs
            )
          )
        # CYTOSET -> FLOWSET
        } else {
          cytoset_to_flowSet(
            cs
          )
        }
      # CYTOFRAME
      } else if(.grepl("cytoframe", format)) {
        # CYTOFRAME -> CYTOFRAME
        if(cyto_class(cs, "flowFrame")) {
          structure(
            list(
              cs
            ),
            names = cyto_names(cs)
          )
        # CYTOSET -> CYTOFRAME
        } else {
          structure(
            lapply(
              cyto_names(cs), 
              function(z) {
                cs[[z]]
              }
            ), names = cyto_names(cs)
          )
        }
      # CYTOSET
      } else if(.grepl("cytoset", format)) {
        # CYTOFRAME -> CYTOSET
        if(cyto_class(cs, "flowFrame")) {
          cytoset(
            structure(
              list(cs), 
              names = cyto_names(cs)
            )
          )
        # CYTOSET -> CYTOSET
        } else {
          if(split) {
            structure(
              lapply(
                seq_along(cs), 
                function(z) {
                  cs[z]
                }
              ), 
              names = cyto_names(cs)
            )
          } else {
            cs
          }
        }
      # MATRIX
      } else if(.grepl("matrix", format)) {
        # CYTOFRAME -> MATRIX
        if(cyto_class(cs, "flowFrame")) {
          structure(
            list(
              cyto_exprs(
                cs, 
                markers = markers,
                drop = FALSE
              )
            ), 
            names = cyto_names(cs)
          )
        # CYTOSET -> MATRIX
        } else {
          structure(
            lapply(
              cyto_names(cs), 
              function(z) {
                cyto_exprs(
                  cs[[z]],
                  markers = markers,
                  drop = FALSE
                )
              }
            ), 
            names = cyto_names(cs)
          )
        }
      }
    }
  )
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
  if(!is.null(channels)) {
    x <- x[, cyto_channels_extract(x, channels)]
  }
  
  # EXTRACT DATA
  mt <- exprs(x)[, , drop = drop, ...]

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
  
  # CHANNELS
  if(!is.null(channels)) {
    x <- x[, cyto_channels_extract(x, channels)]
  }
  
  # EXTRACT DATA
  structure(
    lapply(
      seq_along(x), 
      function(z){
        cyto_exprs(
          x[[z]],
          channels = NULL,
          markers = markers,
          drop = drop,
          ...
        )
      }
    ),
    names = cyto_names(x)
  )
}

## CYTO_EXPRS REPLACEMENT METHOD -----------------------------------------------

#' @noRd
#' @export
"cyto_exprs<-" <- function(object, value) {
  UseMethod("cyto_exprs<-")
}

#' @importFrom flowCore exprs<-
#' @noRd
#' @export
"cyto_exprs<-.flowFrame" <- function(object, value) {
  exprs(object) <- value
  return(object)
}

#' @importFrom flowCore exprs<-
#' @noRd
#' @export
"cyto_exprs<-.flowSet" <- function(object, value) {
  
  # TODO: THROWS ERROR WHEN REPLACING SUBSET OF CYTOSET
  # CYTO_EXPRS(CS[1:2]) <- LIST(MAT, MAT) - S4 NOT SUBSETTABLE
  
  # VALUE - LIST OF MATRICES
  if(!all(LAPPLY(value, "cyto_class", "matrix")) |
     length(value) != length(object)) {
    stop(
      paste0(
        "Replacement value must be a list of matrices, one per each ",
        cyto_class(object[[1]])[1], 
        " in this ",
        cyto_class(object)[1],
        "!"
      )
    )
  }
  # ORDER MATRICES
  if(all(cyto_names(object) %in% names(value))) {
    value <- value[cyto_names(object)]
  }
  # REPLACE EXPRESSION DATA
  lapply(
    seq_along(object),
    function(z) {
      exprs(object[[z]]) <- value[[z]]
    }
  )
  return(object)
}

## CYTO_EMPTY ------------------------------------------------------------------

#' Construct an empty cytoframe
#'
#' @param name name to add to the constructed cytoframe
#' @param channels channels to include in the constructed cytoframe
#' @param ... additional arguments passed to
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#'
#' @importFrom flowCore identifier<- flowFrame
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' # Construct empty cytoframe
#' cyto_empty(name = "Test.csv", 
#' channels = c("FSC-A", "SSC-A", "PE-A"))
#' 
#' @export
cyto_empty <- function(name = NULL,
                       channels = NULL,
                       ...) {
  
  # CHANNELS
  if (is.null(channels)) {
    stop(
      "Supply the names of the channels to include in the cytoframe"
    )
  }
  
  # CONSTRUCT EMPTY FLOWFRAME
  fr <- matrix(
    0,
    ncol = length(channels),
    nrow = 1,
    byrow = TRUE
  )
  colnames(fr) <- channels
  fr <- flowFrame(fr, ...)
  fr <- fr[-1, ]
  
  # NAME
  if (!is.null(name)) {
    identifier(fr) <- name
  }
  
  # UPDATE MARKERS
  if(!is.null(names(channels))) {
    channels <- channels[!names(channels) %in% c("NA", NA)]
    if(length(channels) > 0) {
      cyto_markers(fr) <- channels 
    }
  }
  
  # FLOWFRAME -> CYTOFRAME
  return(
    cyto_convert(
      fr
    )
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
cyto_filter <- function(x,
                        ...) {
  
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
#'   select. Supplied lists can optionally include the above \code{exclude} and
#'   \code{exact} arguments 
#' @param exclude logical to indicate whether negative exclusion indices should
#'   be returned for samples matching the supplied experimental criteria, set to
#'   FALSE by default.
#' @param exact logical indicating whether to only search for samples by exact
#'   matching to experimental criteria, set to FALSE by default.
#'
#' @return vector of indices that indicate the location of the selected
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
                       ...,
                       exclude = FALSE,
                       exact = FALSE) {
  
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
  
  # EXACT
  if("exact" %in% names(args)) {
    exact <- args[["exact"]]
    args <- args[!names(args) %in% "exact"]
  }
  
  # EXCLUDE
  if("exclude" %in% names(args)) {
    exclude <- args[["exclude"]]
    args <- args[!names(args) %in% "exclude"]
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
        # ROWNAMES - EXACT
        if(z %in% rownames(pd)) {
          match(z, rownames(pd))
        # NAME - EXACT
        } else if(z %in% pd[, "name"]) {
          match(z, pd[, "name"])
        # PARTIAL MATCH
        } else if(exact == FALSE) {
          # ROWNAMES - PARTIAL
          if(any(.grepl(z, rownames(pd)))) {
            .grep(z, rownames(pd))
          # NAME - PARTIAL
          } else if(any(.grepl(z, pd[, "name"]))) {
            .grep(z, pd[, "name"])
          # NO MATCH
          } else {
            NULL
          }
        # NO MATCH
        } else {
          NULL
        }
      })
      # NO MATCH
      if(length(ind) == 0) {
        stop(
          "Could not match samples by name using supplied matching criteria!"
        )
      }
    }
  # EXPERIMENTAL VARIABLES
  } else {
    # INDICES PER VARIABLE
    ind <- lapply(
      names(args), 
      function(z) {
        # DIRECT MATCH
        if(z %in% colnames(pd)) {
          var_ind <- match_ind(z, colnames(pd))
          # PARTIAL MATCH
        } else if(exact == FALSE) {
          var_ind <- .grep(
            z, 
            colnames(pd)
          )
        }
        # NO MATCH
        if(length(var_ind) == 0) {
          stop(
            paste0(
              z,
              " is not a valid variable in cyto_details(x)."
            )
          )
          # MULTIPLE MATCHES
        } else if(length(var_ind) > 1) {
          stop(
            paste0(
              z,
              " matches multiple variables in cyto_details(x)."
            )
          )
        }
        # ORDER AS SUPPLIED
        LAPPLY(args[[z]], function(w){
          # EXACT MATCH
          if(w %in% pd[, var_ind]) {
            levels <- which(pd[, var_ind] %in% as.character(w))
            # PARTIAL MATCH
          } else if(exact == FALSE) {
            levels <- .grep(
              as.character(w),
              pd[, var_ind]
            )
          }
          # NO MATCH
          if(length(levels) == 0) {
            stop(
              paste0(
                w, 
                " is not a valid level for ",
                z,
                "!"
              )
            )
          }
          return(levels)
        }
      )
    })
    # INTERSECTION
    if(length(ind) > 1) {
      ind <- Reduce("intersect", ind)
    } else {
      ind <- ind[[1]]
    }
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
#'   by setting \code{exclude = TRUE} and exact matches can be performed by
#'   setting \code{exact = TRUE}. See \code{\link{cyto_match}} or more details.
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
  ind <- cyto_match(x, ...)
  
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
#'   Set this argument to \code{"all"} or \code{NA} to place all samples within
#'   the same group.
#' @param details logical indicating whether the split experimental details
#'   should be returned instead of the group names, set to FALSE by default.
#' @param select vector of indices or named list containing experimental
#'   variables passed to \code{\link{cyto_select}} to be used to select samples
#'   prior to retrieving group-wise experimental details.
#' @param sep indicates the separator to use when constructing new group names
#'   from multiple variables, set to \code{" "} by default. For compatibility
#'   with openCyto, we can set \code{sep = ":"}.
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
                        select = NULL,
                        sep = " "){
  
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
    lapply(
      seq_along(group_by), 
      function(z) {
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
          missing_levels <- as.vector(
            var_levels[!var_levels %in% group_by[[z]]]
          )
          group_by[[z]] <<- c(
            group_by[[z]],
            missing_levels
          )
        }
        # Convert pd variable to factor and set levels
        pd[, var] <<- factor(pd[, var], levels = group_by[[z]])
      }
    )
    # Convert group_by to vector
    group_by <- names(group_by)
    # group_by is a vector of variable names
  } else {
    # GROUP_BY = NA
    if(.all_na(group_by)) {
      group_by <- "all"
    }
    # Check variables
    if(!all(group_by %in% "all") & !all(group_by %in% colnames(pd))) {
      lapply(
        group_by, 
        function(y) {
          if (!y %in% colnames(pd)) {
            stop(paste0(y, " is not a valid variable for this ", class(x), "."))
          }
        }
      )
    }
  }
  
  # Split pd based on group_by into a named list
  if (length(group_by) == 1) {
    if (group_by == "all") {
      pd_split <- list("all" = pd)
    } else if (group_by == "name") {
      pd_split <- lapply(
        cyto_names(x), 
        function(z) {
          pd[rownames(pd) == z, , drop = FALSE] # name column may not match
        }
      )
      names(pd_split) <- cyto_names(x)
    } else {
      pd_split <- split(
        pd, pd[, group_by],
        sep = sep,
        lex.order = TRUE,
        drop = TRUE
      )
    }
  } else {
    pd_split <- split(
      pd, pd[, group_by],
      sep = sep,
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
  pd_split <- cyto_groups(
    x,
    select = select,
    group_by = sort_by,
    details = TRUE
  )
  
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
#'   c("Stim-A","Stim-C","Stim-B", "Stim-D"))). Set this argument to
#'   \code{"all"} or \code{NA} to place all samples within the same group.
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
  pd_split <- cyto_groups(
    x, 
    select = select,
    group_by = group_by,
    details = TRUE
  )
  
  # Replace each element of pd_split with matching samples
  x_list <- lapply(
    seq_len(length(pd_split)), 
    function(z) {
      ind <- match(rownames(pd_split[[z]]), rownames(cyto_details(x)))
      x[ind]
    }
  )
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
#'   prior to merging. Set this argument to \code{"all"} or \code{NA} to merge
#'   all samples.
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
  x <- cyto_data_extract(
    x,
    parent = parent,
    copy = FALSE
  )[[1]]
  
  # BARCODING ------------------------------------------------------------------
  
  # SAMPLE ID
  if(barcode) {
    x <- cyto_barcode(x, ...)
  }
  
  # GROUPING -------------------------------------------------------------------
  
  # CYTO_GROUP_BY
  cs_list <- cyto_group_by(
    x, 
    group_by = merge_by
  )
  
  # COMBINED EVENTS
  if ("all" %in% names(cs_list)) {
    names(cs_list)[which("all" %in% names(cs_list))] <- "Combined Events"
  }
  
  # SELECTION ------------------------------------------------------------------
  
  # ATTEMPT SELECTION OR RETURN ALL SAMPLES
  if (!is.null(select)) {
    cs_list <- lapply(
      cs_list, 
      function(z) {
        tryCatch(cyto_select(z, select), error = function(e) {
          z
          }
        )
      }
    )
  }
  
  # MERGING --------------------------------------------------------------------
  
  # CONVERT EACH GROUP TO MERGED CYTOSET/FLOWSET
  if(.grepl("s", format)) {
    structure(
      lapply(seq_along(cs_list), function(z){
        do.call(
          cyto_class(cs_list[[z]]), # flowSet() / cytoset()
          list(
            structure(
              list(
                as(cs_list[[z]], 
                   cyto_class(cs_list[[z]][[1]])) # flowFrame|cytoframe
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
      lapply(
        seq_along(cs_list), 
        function(z) {
          as(cs_list[[z]], cyto_class(cs_list[[z]][[1]])) # flowFrame|cytoframe
        }
      ),
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
setAs(
  "cytoset",
  "cytoframe",
  function(from){
    flowFrame_to_cytoframe(
      as(from, "flowFrame"),
      emptyValue = FALSE
    )
  }
)

#' Coerce flowSet to cytoframe
#' 
#' @importFrom methods setAs
#' @importFrom flowWorkspace flowFrame_to_cytoframe
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "flowSet",
  "cytoframe",
  function(from){
    flowFrame_to_cytoframe(
      as(from, "flowFrame"),
      emptyValue = FALSE
    )
  }
)

#' matrix to cytoframe
#' 
#' @importFrom flowWorkspace flowFrame_to_cytoframe
#' @importFrom flowCore flowFrame
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "matrix",
  "cytoframe",
  function(from) {
    flowFrame_to_cytoframe(
      flowFrame(
        data.matrix(from)
      ),
      emptyValue = FALSE
    )
  }
)

#' data.frame to cytoframe
#' 
#' @importFrom flowWorkspace flowFrame_to_cytoframe
#' @importFrom flowCore flowFrame
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "data.frame",
  "cytoframe",
  function(from) {
    flowFrame_to_cytoframe(
      flowFrame(
        data.matrix(from)
      ),
      emptyValue = FALSE
    )
  }
)

#' data.frame to cytoset
#' 
#' @importFrom flowWorkspace cytoset cf_get_uri
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "data.frame",
  "cytoset",
  function(from) {
    cf <- as(from, "cytoframe")
    cytoset(
      structure(
        list(cf),
        names = file_ext_remove(
          basename(
            cf_get_uri(cf)
          )
        )
      )
    )
  }
)

#' matrix to cytoset
#' 
#' @importFrom flowWorkspace cytoset cf_get_uri
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "matrix",
  "cytoset",
  function(from) {
    cf <- as(from, "cytoframe")
    cytoset(
      structure(
        list(cf),
        names = file_ext_remove(
          basename(
            cf_get_uri(cf)
          )
        )
      )
    )
  }
)

#' cytoframe to cytoset
#' 
#' @importFrom flowWorkspace cytoset cf_get_uri
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "cytoframe",
  "cytoset",
  function(from) {
    cytoset(
      structure(
        list(from),
        names = file_ext_remove(
          basename(
            cf_get_uri(from)
          )
        )
      )
    )
  }
)

#' flowFrame to cytoset
#' 
#' @importFrom flowWorkspace cytoset cf_get_uri
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
setAs(
  "flowFrame",
  "cytoset",
  function(from) {
    cf <- cyto_convert(from)
    cytoset(
      structure(
        list(
          cf
        ),
        names = file_ext_remove(
          basename(
            cf_get_uri(cf)
          )
        )
      )
    )
  }
)

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
                       id = NULL,
                       names = NULL) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # FLOWFRAME
  if (!cyto_class(x, c("flowFrame", "flowSet"))) {
    stop("cyto_split() requires a cytoframe or cytoset object.")
  }
  
  # ID - MISSING
  if(is.null(id)) {
    id <- "\\Sample.*ID$"
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
        LAPPLY(
          nrow(cnts), 
          function(z){
            rep(z, cnts[z, 1])
          }
        ),
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
  cf_list <- lapply(
    seq_along(x), 
    function(z){
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
        # CYTOEXPLORER - SAMPLE-ID
        if(id %in% "Sample-ID") {
          ids <- factor(
            cyto_exprs(
              cf, 
              channels = id,
              drop = TRUE
            )
          )
        # OTHER CHANNEL - FLOWJO COMPATIBILITY
        } else {
          # EXTRACT CHANNEL
          exprs <- cyto_exprs(
            cf,
            channels = id,
            drop = TRUE
          )
          # COMPUTE KERNEL DENSITY
          d <- cyto_stat_density(
            exprs,
            stat = "density"
          )
          # COMPUTE CUT POINTS
          breaks <- c(
            min(d$x),
            d$x[which(diff(diff(-d$y) >= 0) < 0)],
            max(d$x)
          )
          # CONVERT TO FACTOR
          ids <- cut(
            exprs,
            breaks = breaks,
            labels = seq_len(length(breaks) - 1)
          )
        }
      }
      # BYPASS EMPTY SAMPLES
      cf_list <- split(
        cf,
        ids
      )
      # STORE SPLIT DATA IN SEPARATE H5 FILES
      cf_list <- lapply(cf_list, "cyto_copy")
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
        lapply(
          seq_along(cf_list),
          function(v) {
            identifier(cf_list[[v]]) <<- names(cf_list)[v]
            return(cf_list[[v]])
          }
        ),
        names = names(cf_list)
      )
      return(cf_list)
    }
  )
  
  # FLATTEN CYTOFRAME LIST
  cf_list <- unlist(cf_list)
  
  # RETURN CYTOSET
  return(cytoset(cf_list))
}

## CYTO_SAVE -------------------------------------------------------------------

#' Write samples to FCS files in new folder or save an entire GatingSet
#'
#' @param x object of class \code{cytoframe}, \code{cytoset},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param save_as name of the folder to which the written FCS files should be
#'   saved, set to NULL by default to save the files to the current working
#'   directory. To prevent files being overwritten, it is recommended that
#'   \code{save_as} directory not be manually created before running
#'   \code{cyto_save}.
#' @param parent name of the parent population to extract when a
#'   \code{GatingHierarchy} or \code{GatingSet} object is supplied. If the name
#'   of the parent is supplied the samples will be written to FCS files in
#'   separate named folders within the \code{save_as} directory. Otherwise the
#'   entire \code{GatingSet} or \code{GatingHierarchy} will be saved to the
#'   specified \code{save_as} directory.
#' @param alias an alternative way to \code{parent} to specify the names of the
#'   populations to be saved.
#' @param split logical indicating whether samples merged using
#'   \code{cyto_merge_by} should be split prior to writing FCS files, set to
#'   FALSE by default.
#' @param id passed to \code{cyto_split},  can be either the name of a channel
#'   or marker containing the IDs for samples or a vector containing an index
#'   for each event within \code{x}, set to \code{Sample-ID} channel by default.
#' @param names original names of the samples prior to merging using
#'   \code{cyto_merge_by}, only required when split is TRUE. These names will be
#'   re-assigned to each of split cytoframes.
#' @param overwrite logical flag to control how \code{save_as} should be handled
#'   if files already exist in this directory, users will be asked interactively
#'   if \code{overwrite} is not manually supplied here.
#' @param ... additional arguments passed to \code{\link{cyto_data_extract}} to
#'   allow control over how the data is formatted prior to saving. Setting
#'   \code{inverse = TRUE} will export data on the linear scale but be sure to
#'   set \code{copy = TRUE} as well if you only want to inverse data
#'   transformations for export but leave your data unchanged.
#'
#' @return a GatingHierarchy, GatingSet or a list of cytosets.
#'
#' @importFrom flowCore write.FCS
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
#' @seealso \code{\link{cyto_data_extract}}
#'
#' @export
cyto_save <- function(x,
                      save_as = ".",
                      parent = NULL,
                      alias = NULL,
                      split = FALSE,
                      id = NULL,
                      names = NULL,
                      overwrite = NULL,
                      ...) {
  
  # PARENT
  parent <- c(parent, alias)
  
  # GATINGHIERARCHY OR GATINGSET -> SAVE ENTIRE GATINGHIERARCHY/GATINGSET
  if(cyto_class(x, "GatingSet") & is.null(parent)) {
    # SAVE_AS -> EMPTY DIRECTORY REQUIRED
    if(is.null(save_as)) {
      stop(
        paste0(
          "Supply the name of a new/empty directory to 'save_as' ",
          "to save this ",
          cyto_class(x), "."
        )
      )
    }
    # EMPTY DIRECTORY CHECK
    if(length(list.files(save_as)) > 0) {
      # OVERWRITE
      if(is.null(overwrite)) {
        overwrite <- cyto_enquire(
          paste0(
            save_as,
            " is non-empty. Overwrite? (Y/N)"
          ),
          options = c("Y", "T")
        )
      }
      # REPLACE DIRECTORY
      if(overwrite){
        unlink(
          save_as, 
          recursive = TRUE
        )
        # ABORT
      } else {
        invisible(NULL)
      }
    }
    # SAVE GATINGSET
    message(
      paste0(
        "Saving ",
        cyto_class(x), 
        " to ", 
        save_as, 
        "..."
      )
    )
    suppressMessages(
      save_gs(x, save_as)
    )
    # RETURN GATINGSET - INVISIBLE DOES NOT TERMINATE FUNCTION
    invisible(x)
    # WRITE FCS FILES
  } else {
    # EXTRACT DATA - NAMED LIST PER PARENT
    cs_list <- cyto_data_extract(
      x,
      parent = parent,
      format = "cytoset",
      ...
    )
    # IF NAMED LIST CREATE NEW FOLDER PER POPULATION
    cs_list <- structure(
      lapply(
        seq_along(cs_list),
        function(z) {
          # CYTOSET
          cs <- cs_list[[z]]
          # SPLIT DATA
          if(split) {
            cs <- cyto_split(
              cs,
              names = names,
              id = id
            )
          }
          # CREATE SAVE_AS DIRECTORY
          if(!dir.exists(save_as)) {
            dir.create(save_as)
          }
          # PARENTAL DIRECTORIES
          if(!is.null(names(cs_list)[z])) {
            # CREATE NEW PARENTAL DIRECTORY IN SAVE_AS
            dir <- paste0(
              save_as, 
              .Platform$file.sep,
              names(cs_list)[z]
            )
            if(!dir.exists(dir)) {
              dir.create(dir)
            }
            # SAVE_AS DIRECTORY
          } else {
            dir <- save_as
          }
          # DIRECTORY CHECK
          if(dir.exists(dir) & 
             any(list.files(dir) %in% cyto_names(cs))) {
            # OVERWRITE ENQUIRY
            if(is.null(overwrite)) {
              # FILES WILL BE OVERWRITTEN
              message(
                paste0(
                  "The following files will be overwritten in ",
                  dir, 
                  ":\n",
                  paste(
                    cyto_names(cs)[cyto_names(cs) %in% list.files(dir)],
                    collapse = "\n"
                  )
                )
              )
              # OVERWRITE?
              overwrite <- cyto_enquire(
                "Do you want to continue? (Y/N)",
                options = c("Y", "T")
              )
            }
            # ABORT
            if(!overwrite) {
              return(cs)
            }
          }
          # MESSAGE - WRITING FCS FILES
          message(
            paste0(
              "\nWriting FCS files to ",
              dir,
              "..."
            )
          )
          # INDIVIDUAL FCS FILES
          lapply(
            seq_along(cs),
            function(v) {
              # FCS FILE NAME
              message(
                cyto_names(cs)[v]
              )
              # WRITE FCS FILES
              write.FCS(
                cs[[v]],
                paste0(
                  dir,
                  .Platform$file.sep,
                  cyto_names(cs)[v]
                )
              )
            }
          )
        }
      )
    )
    # RETURN (SPLIT) DATA
    invisible(cs_list)
  }
  
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
#' @param events can be either a numeric [0,1] or integer to indicate the
#'   percentage or number of events to keep respectively. For \code{cytoset} or
#'   \code{GatingSet} objects, the same degree of sampling is applied to each
#'   \code{cytoframe} by default. However, cytoframes will be separately sampled
#'   if a vector the same length as the \code{cytoset} or \code{GatingSet} is
#'   supplied. Setting \code{events} to NA will result in downsampling to the
#'   minimum number of events across samples.
#' @param seed value used to \code{set.seed()} internally. Setting a value for
#'   seed will return the same result with each run.
#' @param ... not in use.
#'
#' @return object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} sampled to
#'   \code{events} events.
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
                                        events = 1,
                                        seed = NULL,
                                        ...) {
  
  # TRANSFORMERS
  trans <- cyto_transformers_extract(x)
  
  # COMPENSATION
  spill <- cyto_spillover_extract(x)
  
  # GATES
  gateTemplate <- cyto_gateTemplate(x)
    
  # EXTRACT CYTOSET
  cs <- cyto_copy(
    gs_cyto_data(x)
  )
  
  # SAMPLING
  cs <- cyto_sample(
    cs,
    events = events,
    seed = seed,
    ...
  )
  
  # REVERSE TRANSFORMATIONS
  if(!.all_na(trans)) {
    cs <- cyto_transform(
      cs,
      trans = trans,
      inverse = TRUE,
      plot = FALSE,
      quiet = TRUE,
      copy = FALSE
    )
  }
  
  # REVERSE COMPENSATION
  if(!is.null(spill)) {
    cs <- cyto_compensate(
      cs,
      spillover = spill,
      remove = TRUE,
      copy = FALSE
    )
  }
  
  # RECONSTRUCT GATINGHIERARCHY
  message(
    "Rebuilding the GatingHierarchy..."
  )
  x <- GatingSet(cs)[[1]]
  
  # RE-APPLY COMPENSATION
  if(!is.null(spill)) {
    message(
      "Compensating for fluorescence spillover..."
    )
    x <- cyto_compensate(
      x,
      spillover = spill,
      copy = FALSE
    )
  }
  
  # RE-APPLY TRANSFORMATIONS
  if(!.all_na(trans)) {
    message(
      "Re-applying data transformations..."
    )
    x <- cyto_transform(
      x,
      trans = trans,
      plot = FALSE,
      quiet = TRUE,
      copy = FALSE
    )
  }
  
  # TRANSFER GATES
  message(
    "Transferring gates..."
  )
  x <- cyto_gateTemplate_apply(
    x,
    gateTemplate
  )
  
  # RETURN SAMPLED GATINGHIERACHY
  return(x)
  
}

#' @rdname cyto_sample
#' @export
cyto_sample.GatingSet <- function(x,
                                  events = 1,
                                  seed = NULL,
                                  ...){
  
  # TRANSFORMERS
  trans <- cyto_transformers_extract(x)
  
  # COMPENSATION
  spill <- cyto_spillover_extract(x)
  
  # GATES
  gateTemplate <- cyto_gateTemplate(x)
  
  # EXTRACT CYTOSET
  cs <- cyto_copy(
    gs_cyto_data(x)
  )
  
  # SAMPLING
  cs <- cyto_sample(
    cs,
    events = events,
    seed = seed,
    ...
  )
  
  # REVERSE TRANSFORMATIONS
  if(!.all_na(trans)) {
    cs <- cyto_transform(
      cs,
      trans = trans,
      inverse = TRUE,
      plot = FALSE,
      quiet = TRUE,
      copy = FALSE
    )
  }
  
  # REVERSE COMPENSATION
  if(!is.null(spill)) {
    cs <- cyto_compensate(
      cs,
      spillover = spill,
      remove = TRUE,
      copy = FALSE
    )
  }
  
  # RECONSTRUCT GATINGHIERARCHY
  message(
    "Rebuilding the GatingSet..."
  )
  x <- GatingSet(cs)
  
  # RE-APPLY COMPENSATION
  if(!is.null(spill)) {
    message(
      "Compensating for fluorescence spillover..."
    )
    x <- cyto_compensate(
      x,
      spillover = spill,
      copy = FALSE
    )
  }
  
  # RE-APPLY TRANSFORMATIONS
  if(!.all_na(trans)) {
    message(
      "Re-applying data transformations..."
    )
    x <- cyto_transform(
      x,
      trans = trans,
      plot = FALSE,
      quiet = TRUE,
      copy = FALSE
    )
  }
  
  # TRANSFER GATES
  message(
    "Transferring gates..."
  )
  x <- cyto_gateTemplate_apply(
    x,
    gateTemplate
  )
  
  # RETURN SAMPLED GATINGSET
  return(x)
  
}

#' @rdname cyto_sample
#' @export
cyto_sample.flowFrame <- function(x,
                                  events = 1,
                                  seed = NULL,
                                  ...) {
  
  # NO SAMPLING - EMPTY FLOWFRAME
  if (nrow(exprs(x)) == 0) {
    return(x)
  }
  
  # Do nothing if no sampling required
  if (events != 1) {
    # Number of events
    n <- nrow(x)
    # n is the number of events to keep
    if (events > 1) {
      # n is too large - retain all events
      if (events > n) {
        return(x)
        # n is sample of x
      } else {
        size <- events
      }
      # n is a proportion of events to keep
    } else {
      # Size
      size <- events * n
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
                                events = 1,
                                seed = NULL,
                                ...) {
  
  # SAMPLE TO MINIMUM EVENTS
  if(.all_na(events)) {
    events <- cyto_apply(
        x,
        "nrow",
        input = "matrix",
        copy = FALSE
      )
    # BYPASS ZERO EVENT SAMPLING
    if(min(events) == 0) {
      warning(
        paste0(
          "The following samples contain zero events: \n",
          paste0(
            rownames(events[events == 0, ,drop = FALSE]),
            collapse = "\n"
          ),
          "\n",
          "Downsampling every sample to zero events is illogical so the ",
          "samples have been returned as is."
        )
      )
      return(x)
    }
    # DOWNSAMPLE TO MINIMUM EVENTS
    events <- min(events)
    message(
      paste("Downsampling each sample to",
            events, "events.")
    )
  }
  
  # PREPARE EVENTS - INDIVIDUAL CYTOFRAMES
  events <- rep(events, length.out = length(x))
  
  # SAPLE EACH CYTOFRAME
  cytoset(
    structure(
      lapply(
        seq_along(x), 
        function(z){
          cyto_sample(
            x[[z]],
            events = events[z],
            seed = seed,
            ...
          )
        }
      ),
      names = cyto_names(x)
    )
  )

}

#' @rdname cyto_sample
#' @export
cyto_sample.list <- function(x,
                             events = 1,
                             seed = NULL,
                             ...) {
  
  # SAME SAMPLE SIZE PER LAYER
  events <- rep(events, length(x))
  
  # SAMPLING
  x <- mapply(
    function(x,
             events) {
      cyto_sample(
        x, 
        events = events,
        seed = seed
      )
    }, 
    x, 
    events
  )
  
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
#' @param events can be either a numeric [0,1] or integer to indicate the
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
#'  events = 50000)
#'
#' @export
cyto_sample_n <- function(x,
                          events = 1,
                          parent = NULL) {
  
  # COUNTS
  counts <- cyto_apply(
    x,
    "nrow",
    parent = parent,
    input = "matrix",
    copy = FALSE
  )

  # FORMAT
  counts <- structure(
    counts[, 1, drop = TRUE],
    names = rownames(counts)
  )
  
  # TOTAL COUNT
  total_count <- sum(counts)
  
  # DISPLAY CANNOT BE LARGER THAN TOTAL COUNT
  if(events > 1 & events > total_count) {
    events <- total_count # n may be set to zero here
  }
  
  # DISPLAY - COUNT
  if(events <= 1) {
    events <- events * total_count
  }

  # DISPLAY - FREQ
  if(events == 0) {
    events_freq <- 0 # avoid possible NaN
  } else {
    events_freq <- events / total_count
  }
  
  # APPROX COUNTS
  sample_counts <- ceiling(events_freq * counts)
  sample_total <- sum(sample_counts)
  
  # EXACT COUNTS - REMOVE EXCESS EVENTS
  excess <- sample_total - events
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
#' @param events numeric to control how events should be sampled when coercing
#'   samples. If a single value is supplied it is passed to \code{cyto_sample_n}
#'   to compute the number of events to extract from each sample to obtain
#'   \code{events} in the coerced data. Alternatively, users can manually supply
#'   a \code{events} value for each file for more control over the way in which
#'   samples are sampled prior to merging.
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
#' # Coerce & Sample 50000 total events
#' cyto_coerce(cs, events = 50000, seed = 56)
#'
#' # Sample 5000 events each & coerce
#' cyto_coerce(cs, events = rep(5000, 33))
#'
#' @export
cyto_coerce <- function(x,
                        events = 1,
                        seed = NULL,
                        barcode = FALSE,
                        format = "cytoset",
                        name = "merge",
                        overwrite = NULL,
                        ...) {
  # IDENTIFIERS
  ids <- cyto_names(x)
  
  # SAMPLE REQUIRES CYTO_SAMPLE_N
  if(length(events) == 1) {
    events <- cyto_sample_n(
      x,
      events = events
    )
  # SAMPLE MANUALLY SUPPLIED
  } else if(length(events) == length(x)) {
    if(is.null(names(events))) {
      names(events) <- ids
    }
  # SAMPLE PER CYTOFRAME REQUIRED
  } else {
    stop(
      paste0(
        "'events' must have the same length as 'x' to use custom sampling ",
        "per cytoframe!"
      )
    )
  }
  
  # REMOVE EMPTY SAMPLES
  if(any(!names(events) %in% ids)) {
    x <- cyto_select(
      x,
      list("name" = ids[ids %in% names(events)])
    )
  }

  # SAMPLING
  x <- cytoset(
    structure(
      lapply(
        cyto_names(x), 
        function(z){
          fr <- cyto_sample(
            x[[z]], 
            events = events[z], 
            seed = seed
          )
          if(cyto_class(fr, "flowFrame", TRUE)) {
            fr <- flowFrame_to_cytoframe(
              fr,
              emptyValue = FALSE # REQUIRED?
            )
          }
          return(fr)
        }
      ),
    names = cyto_names(x))
  )

  # BARCODE
  if(barcode) {
    x <- cyto_barcode(
      x,
      overwrite = overwrite
    )
  }
  
  # COERCE
  if(length(x) == 1) {
    x <- x[[1]]
  } else {
    x <- as(x, "flowFrame")
  }
  
  # CYTOFRAME
  if(cyto_class(x, "flowFrame", TRUE)) {
    x <- flowFrame_to_cytoframe(
      x,
      emptyValue = FALSE # REQUIRED?
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
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param node name of the gated bead population to use for the calculation. If
#'   not supplied internal checks will be made for populations named "Single
#'   Beads" or "Beads".
#' @param events minimum number of events to downsample to, set to the minimum
#'   count for the specified node across samples by default. \code{events} must
#'   be less than or equal to the minimum event count across the samples.
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
                                events = NULL,
                                ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGHIERARCHY OR GATINGSET
  if (!cyto_class(x, "GatingSet")) {
    stop("'x' must be an object of class GatingHierarchy or Gatingset.")
  }
  
  # NODES
  nodes <- cyto_nodes(
    x, 
    path = "auto"
  )
  
  # BEADS
  if (is.null(node)) {
    # SINGLE BEADS
    if (any(.grepl("single beads", nodes))) {
      node <- nodes[.grep("single beads", nodes)[1]]
      # BEADS
    } else if (any(.grepl("beads", nodes))) {
      node <- nodes[.grep("beads", nodes)[1]]
      # BEADS MISSING
    } else {
      stop("Supply the name of the 'node' to downsample the GatingSet.")
    }
  }
  
  # TRANSFORMERS
  trans <- cyto_transformers_extract(x)
  
  # COMPENSATION
  spill <- cyto_spillover_extract(x)
  
  # GATES
  gateTemplate <- cyto_gateTemplate(x)
  
  # COMPUTE SAMPLE SIZES -------------------------------------------------------
  
  # NODE COUNTS
  node_counts <- cyto_stats_compute(
    x,
    alias = node,
    stat = "count",
    details = FALSE
  )
  node_counts <- node_counts[, ncol(node_counts), drop = TRUE]
  # node_pops <- cyto_data_extract(gs_clone, 
  #                                parent = node,
  #                                copy = FALSE)[[node]]
  # node_counts <- suppressMessages(
  #   cyto_stats_compute(node_pops,
  #                      stat = "count",
  #                      format = "wide",
  #                      details = FALSE)
  # )
  # node_counts <- node_counts[, ncol(node_counts), drop = TRUE]
  
  # # NODE MISSING EVENTS
  # if (any(node_counts == 0)) {
  #   ind <- which(node_counts == 0)
  #   stop(
  #     paste0(
  #     "The following samples do not contain any events in the specified node:",
  #     paste0("\n", cyto_names(gs_clone)[ind])
  #     )
  #   )
  # }
  
  # DOWNSAMPLE EACH NODE - EXTRACT EVENT-ID FROM ROOT
  
  # NODE COUNT - BYPASS ZERO EVENT SAMPLE NODES
  if (is.null(events)) {
    events <- min(node_counts[node_counts > 0])
  } else {
    if (!events <= min(node_counts[node_counts > 0])) {
      events <- min(node_counts[node_counts > 0])
    }
  }
  
  # NODE RATIOS
  node_ratios <- LAPPLY(
    node_counts, 
    function(z) {
      if(z > 0) {
        return(1 / (z / events))
      # BYPASS ZERO EVENT SAMPLES
      } else {
        return(0)
      }
    }
  )
  
  # SAMPLE ROOT NODE -----------------------------------------------------------
  
  # SAMPLED ROOT CYTOSET
  cs <- cyto_data_extract(
    x,
    parent = "root",
    events = node_ratios, # DOES THIS WORK?
    copy = TRUE
  )[[1]]
  
  # INVERSE TRANSFORMATIONS
  if(!.all_na(trans)) {
    cs <- cyto_transform(
      cs,
      trans = trans,
      inverse = TRUE,
      plot = FALSE,
      quiet = TRUE,
      copy = FALSE
    )
  }

  # REMOVE COMPENSATION
  if(!is.null(spill)) {
    cs <- cyto_compensate(
      cs,
      spillover = spill,
      remove = TRUE,
      copy = FALSE
    )
  }
  
  # REBUILD GATINGHIERARCHY/GATINGSET ------------------------------------------
  
  # RECONSTRUCT GATINGSET
  message(
    paste0(
      "Rebuilding the ",
      cyto_class(x), 
      "..."
    )
  )
  x <- GatingSet(cs)
  
  # RE-APPLY COMPENSATION
  if(!is.null(spill)) {
    message(
      "Compensating for fluorescence spillover..."
    )
    x <- cyto_compensate(
      x,
      spillover = spill,
      copy = FALSE
    )
  }
  
  # RE-APPLY TRANSFORMATIONS
  if(!.all_na(trans)) {
    message(
      "Re-applying data transformations..."
    )
    x <- cyto_transform(
      x,
      trans = trans,
      plot = FALSE,
      quiet = TRUE,
      copy = FALSE
    )
  }
  
  # GATINGHIERARCHY
  if(cyto_class(gateTemplate, "gateTemplate")) {
    x <- x[[1]]
  }
  
  # TRANSFER GATES
  message(
    "Transferring gates..."
  )
  x <- cyto_gateTemplate_apply(
    x,
    gateTemplate
  )
  
  # RETURN SAMPLED GATINGSET
  return(x)
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
#' @return barcoded cytoset or GatingSet with \code{"Sample-ID"} and/or
#'   \code{"Event-ID"} column added and annotated.
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
  cs <- cyto_data_extract(
    x,
    parent = "root",
    copy = FALSE,
    format = "cytoset"
  )[[1]]
  
  # TYPE
  if (.grepl("^b", type)) {
    type <- c("samples", "events")
  }
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # BARCODE?
  barcode <- TRUE
  
  # CHECK FOR SAMPLE IDs
  if (any(.grepl("^s", type))){
    # SAMPLE IDs EXIST - BACKWARDS COMPATIBLE
    if(any(.grepl("^Sample-?ID$", cyto_channels(cs)))) {
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
        cs <- cyto_copy(
          cs[, -which(.grepl("^Sample-?ID$", cyto_channels(cs)))]
        )
      } else {
        barcode <- FALSE
      }
    }
    # BARCODE SAMPLES
    if(barcode){
      cs <- cytoset(
        structure(
          lapply(
            seq_along(cs), 
            function(z) {
              suppressWarnings(
                cyto_cbind(
                  cs[[z]],
                  matrix(
                    rep(z, cyto_stat_count(cs[[z]])),
                    ncol = 1,
                    dimnames = list(NULL, "Sample-ID")
                  )
                )
              )
            }
          ),
          names = cyto_names(cs)
        )
      )
    }
  }
  
  # CHECK FOR EVENT IDs
  if (any(.grepl("^e", type))) {
    # EVENT IDs EXIST
    if(any(.grepl("^Event-?ID$", cyto_channels(cs)))) {
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
        cs <- cyto_copy(
          cs[, -which(.grepl("^Event-?ID$",cyto_channels(cs)))]
        )
      } else {
        barcode <- FALSE
      }
    }
    # BARCODE EVENTS
    if(barcode) {
      # EVENT IDs
      cnt <- 0
      event_ids <- lapply(
        seq_along(cs), 
        function(z){
          if(nrow(cs[[z]]) == 0) {
            return(NA)
          } else {
            ids <- seq(cnt + 1, cnt + nrow(cs[[z]]))
            cnt <<- ids[length(ids)]
            return(ids)
          }
        }
      )
      cs <- cytoset(
        structure(
          lapply(
            seq_along(cs), 
            function(z) {
              suppressWarnings(
                cyto_cbind(
                  cs[[z]],
                  matrix(
                    event_ids[[z]],
                    ncol = 1,
                    dimnames = list(NULL, "Event-ID")
                  )
                )
              )
            }
          ),
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
#' @param save_as name of a CSV to which the edited marker assignments should
#'   be saved, defaults to \code{Experiment-Markers.csv} prefixed with the date
#'   unless \code{file} is specified or a file is found when searching the
#'   current directory for experimental markers. Setting this argument to NA
#'   will prevent writing of edited marker assignments to file. Custom file
#'   names passed to \code{save_as} should end in \code{"Details.csv"} so that
#'   CytoExploreR can automatically find and read in these details when
#'   required.
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
#' @name cyto_markers_edit
NULL

#' @rdname cyto_markers_edit
#' @export
cyto_markers_edit <- function(x,
                              file = NULL,
                              save_as = NULL,
                              ...) {
  
  # CHANNELS
  chans <- cyto_channels(x)
  
  # CHANNELS/MARKERS DATA.FRAME
  pd <- data.frame(
    "channel" = chans,
    "marker" = names(chans)
  )
  
  # FILE MISSING
  if(is.null(file)) {
    # FILE SEARCH - BYPASS SUFFIX
    pd_new <- cyto_file_search(
      "Markers.*\\.csv$",
      colnames = c("channel", "marker"),
      channel = chans
    )
    # SINLGE FILE MARKER ASSIGNMENTS FOUND
    if(length(pd_new) > 0) {
      # SINGLE FILE
      if(length(pd_new) == 1) {
        message(
          paste0(
            "Importing saved marker assignments from ",
            names(pd_new)[1],
            "..."
          )
        )
        file <- names(pd_new)[1]
        pd <- pd_new[[1]]
        # MULTIPLE FILES
      } else {
        # ENQUIRE
        if(interactive() & cyto_option("CytoExploreR_interactive")) {
          message(
            paste0(
              "Multiple files found with marker assignments for this ",
              cyto_class(x),
              ". Which file would you like to import marker assignments from?"
            )
          )
          message(
            paste0(
              paste0(
                1:length(pd_new), 
                ": ",
                names(pd_new)
              ),
              sep = "\n"
            )
          )
          opt <- cyto_enquire(NULL)
          opt <- tryCatch(
            as.numeric(opt),
            warning = function(w){
              return(
                match(opt, names(pd_new))
              )
            }
          )
          file <- names(pd_new)[opt]
          pd <- pd_new[[opt]]
        } else {
          warning(
            paste0(
              "Multiple files found with marker assignments for this ",
              cyto_class(x),
              " - resorting to using ",
              names(pd_new)[1],
              "..."
            )
          )
          file <- names(pd_new)[1]
          pd <- pd_new[[1]]
        }
      }
    } 
    # FILE SUPPLIED
  } else {
    # FILE EXTENSION
    file <- file_ext_append(file, ".csv")
    # FILE EXISTS
    if(file_exists(file, error = FALSE)) {
      message(
        paste0(
          "Importing marker assignments from ",
          file,
          "..."
        )
      )
      # READ FILE
      pd_new <- read_from_csv(file)
      # ALL CHANNELS REQUIRED
      if(any(!pd$channel %in% pd_new$channel)) {
        pd <- rbind(
          pd_new,
          pd[which(!pd$channel %in% pd_new$channel), , drop = FALSE]
        )
      } else {
        pd <- pd_new
      }
    }
  }

  # DEFAULT FILE NAME
  if (is.null(save_as)) {
    if(is.null(file)) {
      save_as <- cyto_file_name(
        paste0(
          format(
            Sys.Date(), 
            "%d%m%y"
          ), 
          "-Experiment-Markers.csv"
        )
      )
    } else {
      save_as <- file
    }
  }
  
  # EDIT - NON-INTERACTIVE UPDATE WITHOUT EDITING
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    pd_new <- data_edit(
      pd,
      logo = CytoExploreR_logo(),
      title = "Experiment Markers Editor",
      row_edit = FALSE, # cannot add/remove rows
      col_edit = FALSE,
      col_names = c("channel", "marker"),
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane",
      ...
    )
  }
  
  # WRITE CSV
  if(!.all_na(save_as)) {
    write_to_csv(
      pd_new,
      save_as,
      row.names = FALSE
    )
  }
  
  # ONLY UPDATE CHANNELS/MARKERS RELEVANT TO DATA
  ind <- match_ind(chans, pd$channel)
  cyto_chans <- pd_new$channel[ind]
  cyto_marks <- pd_new$marker[ind]
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

## CYTO_CHANNELS_EDIT ----------------------------------------------------------

#' @rdname cyto_markers_edit
#' @export
cyto_channels_edit <- function(x,
                               file = NULL,
                               save_as = NULL,
                               ...) {
  
  cyto_markers_edit(
    x,
    file = file,
    save_as = save_as,
    ...
  )
  
}

## CYTO_DETAILS_EDIT -----------------------------------------------------------

#' Interactively edit cyto_details for a cytoset or GatingSet
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param file name of csv file containing the experimental information
#'   associated with each sample in the supplied data. This file must contain a
#'   column called \code{"name"} which contains the names of the samples.
#' @param save_as name of a CSV to which the edited experimental details should
#'   be saved, defaults to \code{Experiment-Details.csv} prefixed with the date
#'   unless \code{file} is specified or a file is found when searching the
#'   current directory for experimental details. Setting this argument to NA
#'   will prevent writing of edited experimental details to file. Custom file
#'   names passed to \code{save_as} should end in \code{"Details.csv"} so that
#'   CytoExploreR can automatically find and read in these details when
#'   required.
#' @param ... additional arguments passed to \code{DataEditR::data_edit()} when
#'   \code{cyto_details_edit()} is run in interactive mode.
#'
#' @return NULL and return \code{flowSet} or \code{GatingSet} with updated
#'   experimental details.
#'
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
#'
#' @export
cyto_details_edit <- function(x,
                              file = NULL,
                              save_as = NULL,
                              ...) {
  
  # CYTOSET/GATINGSET
  if (!cyto_class(x, "flowSet") & !cyto_class(x, "GatingSet", TRUE)) {
    stop("Please supply either a flowSet or a GatingSet")
  }
  
  # FILE MISSING
  if (is.null(file)) {
    # FILE SEARCH - BYPASS SUFFIX
    pd <- cyto_file_search(
      "Details.*\\.csv$", # Compensation-Details | Experiment-Details
      colnames = "name",
      rownames = rownames(cyto_details(x))
    )
    # NO DETAILS FOUND
    if(length(pd) == 0) {
      pd <- cyto_details(x)
    # DETAILS FOUND
    } else {
      # MULTIPLE FILES FOUND
      if(length(pd) > 1) {
        # ENQUIRE
        if(interactive() & cyto_option("CytoExploreR_interactive")) {
          message(
            paste0(
              "Multiple files found with experimental details for this ",
              cyto_class(x),
              ". Which file would you like to import experiment details from?"
            )
          )
          message(
            paste0(
              paste0(
                1:length(pd), 
                ": ",
                names(pd)
              ),
              sep = "\n"
            )
          )
          opt <- cyto_enquire(NULL)
          opt <- tryCatch(
            as.numeric(opt),
            warning = function(w){
              return(
                match(opt, names(pd))
              )
            }
          )
          file <- names(pd)[opt]
          pd <- pd[[opt]]
        } else {
          warning(
            paste0(
              "Multiple files found with experimental details for this ",
              cyto_class(x),
              " - resorting to using ",
              names(pd)[1],
              "..."
            )
          )
          file <- names(pd)[1]
          pd <- pd[[1]]
        }
      # SINGLE FILE FOUND
      } else {
        message(
          paste0(
            "Importing saved experiment details from ",
            names(pd)[1],
            "..."
          )
        )
        file <- names(pd)[1]
        pd <- pd[[1]]
      }
    }
  # FILE MANUALLY SUPPLIED
  } else {
    # FILE EXTENSION
    file <- file_ext_append(file, ".csv")
    # FILE EXISTS
    if(file_exists(file, error = FALSE)) {
      message(
        paste0(
          "Importing experiment details from ",
          file,
          "..."
        )
      )
      # READ FILE
      pd <- read_from_csv(file)
      # CHECK FILE
      if(!"name" %in% colnames(pd) |
         !all(rownames(cyto_details(x)) %in% rownames(pd))) {
        stop(
          paste0(
            file,
            " must have rownames and contain entries for every sample in ",
            "this ",
            cyto_class(x),
            "!"
          )
        )
      }
    } else {
      warning(
        paste0(
          file,
          " does not exist! Resorting to cyto_details(x) instead..."
        )
      )
      pd <- cyto_details(x)
    }
  }
  
  # SAVE_AS
  if(is.null(save_as)) {
    if(is.null(file)) {
      save_as <- cyto_file_name(
        paste0(
          format(Sys.Date(), "%d%m%y"),
          "-Experiment-Details.csv"
        )
      )
    } else {
      save_as <- file
    }
  }
  
  # EDIT DETAILS - ROWNAMES CANNOT BE EDITED
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    # REMOVE ROW NAMES
    cyto_names <- rownames(pd)
    rownames(pd) <- NULL
    # EDIT
    pd <- data_edit(
      pd,
      logo = CytoExploreR_logo(),
      title = "Experiment Details Editor",
      row_edit = FALSE,
      # col_readonly = "name",
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane",
      ...
    )
    # REPLACE ROW NAMES - REQUIRED TO UPDATE CYTO_DETAILS
    rownames(pd) <- cyto_names
  }
  
  # DETAILS REQUIRE ROWNAMES FOR UPDATING - (ROWNAMES MISSING IN FILE)
  if(is.null(rownames(pd))) {
    rownames(pd) <- pd[, "name"]
  }
  
  # UPDATE DETAILS - PD MAY CONTAIN EXTRA INFORMATION
  cyto_details(x) <- pd[match_ind(
    rownames(
      cyto_details(x)
    ), 
    rownames(pd)
  ), , drop = FALSE]
   
  # SAVE UPDATED DETAILS - CANNOT SAVE ABOVE AS ROWNAMES REMOVED
  if(!.all_na(save_as)) {
    write_to_csv(
      pd,
      file = save_as,
      row.names = TRUE
    )
  }
  
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
    save_as <- cyto_file_name(
      paste0(
        format(Sys.Date(), "%d%m%y"),
        "-Experiment-Details.csv"
      )
    )
  }
  
  # WRITE CSV FILE
  pd <- cyto_details(x)
  write_to_csv(pd, save_as)
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
#' @param copy logical to control whether the supplied data should be copied
#'   prior to applying compensation for fluorescent spillover, set to FALSE by
#'   default. Setting this argument to TRUE will leave the supplied data
#'   untouched and return a copy of the data with compensation applied.
#' @param ... not in use.
#'
#' @return a compensated \code{flowFrame}, \code{flowSet} or \code{GatingSet}
#'   object.
#'
#' @importFrom flowCore decompensate
#' @importFrom flowWorkspace cytoset flowFrame_to_cytoframe GatingSet
#'   gs_cyto_data
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
                                      remove = FALSE,
                                      copy = FALSE,
                                      ...) {
  
  # WARNING
  if(!is.null(cyto_spillover_extract(x)) & !remove) {
    message(
      paste0(
        "Compensation has already been applied to this ",
        cyto_class(x, class = TRUE),
        "!"
      )
    )
    # ENQUIRE
    if(interactive() & cyto_option("CytoExploreR_interactive")) {
      opt <- cyto_enquire(
        "Do you want to continue? (Y/N)",
        options = c("Y", "T")
      )
    # ABORT
    } else {
      opt <- FALSE
    }
    # ABORT
    if(!opt) {
      return(x)
    }
  }
  
  # TRANSFORMATIONS
  trans <- cyto_transformers_extract(x)
  
  # CYTOSET
  cs <- cyto_data_extract(
    x, 
    parent = "root",
    copy = copy,
  )[["root"]]
  
  # PREPARE SPILLOVER MATRIX
  spill <- .cyto_spillover_prepare(
    cs,
    spillover = spillover,
    select = select
  )
  
  # INVERSE TRANSFORMATIONS
  if(!.all_na(trans)) {
    message(
      "Applying inverse data transformations..."
    )
    cs <- cyto_transform(
      cs,
      trans = trans,
      inverse = TRUE,
      quiet = TRUE,
      plot = FALSE
    )
  }
  
  # DECOMPENSATE
  if(remove) {
    message(
      "Removing applied compensation..."
    )
    cf_list <- lapply(
      seq_along(cs), 
      function(z) {
        cyto_convert(
          decompensate(
            cs[[z]], 
            spill[[z]]
          )
        )
      })
    names(cf_list) <- cyto_names(cs)
    cs <- cytoset(cf_list)
  }
  
  # CREATE NEW GATINGSET
  gs <- GatingSet(cs)
  
  # COMPENSATE
  if(!remove) {
    message(
      "Compensating for fluorescence spillover..."
    )
    flowWorkspace::compensate(gs, spill)
  }
  
  # RE-APPLY TRANSFORMATIONS
  if(!.all_na(trans)) {
    message(
      "Re-applying data transformations..."
    )
    gs <- cyto_transform(
      gs,
      trans = trans,
      inverse = FALSE,
      quiet = TRUE,
      plot = FALSE
    )
  }
  
  # TRANSFER GATES
  if(length(cyto_nodes(x)) > 1) {
    message(
      "Recomputing gates..."
    )
    x <- cyto_gateTemplate_apply(
      gs,
      x
    )
  } else {
    x <- gs
  }
  
  # RETURN UPDATED GATINGSET
  return(x)
}

#' @rdname cyto_compensate
#' @export
cyto_compensate.flowSet <- function(x,
                                    spillover = NULL,
                                    select = 1,
                                    remove = FALSE,
                                    copy = FALSE,
                                    ...) {
  
  # PREPARE SPILLOVER MATRIX
  spill <- .cyto_spillover_prepare(
    x,
    spillover = spillover,
    select = select
  )
  
  # COPY
  if(copy) {
    x <- cyto_copy(x)
  }
  
  # REMOVE COMPENSATION
  if (remove == TRUE) {
    cf_list <- lapply(
      seq_along(x), 
      function(z) {
        cyto_convert(
          decompensate(
             x[[z]],
             spill[[z]]
          )
        )
      }
    )
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
                                      copy = FALSE,
                                      ...) {
  
  # PREPARE SPILLOVER MATRIX
  spill <- .cyto_spillover_prepare(
    x,
    spillover = spillover
  )[[1]]
  
  # COPY
  if(copy) {
    x <- cyto_copy(x)
  }
  
  # REMOVE COMPENSATION
  if (remove == TRUE) {
    decompensate(x, spill)
    # APPLY COMPENSATION
  } else if (remove == FALSE) {
    flowCore::compensate(x, spill)
  }
}

## .CYTO_SPILLOVER_PREPARE -----------------------------------------------------

#' Prepare a supplied spillover matrix
#'
#' @param x cytoset or GatingSet.
#' @param spillover an array, list or character vector.
#' @param select index or name of the sample from which the spillover matrix
#'   should be extracted when no spillover matrix file is supplied to
#'   \code{spillover}. To compensate each sample individually using their stored
#'   spillover matrix file, set \code{select} to NULL.
#'
#' @return a named list of spillover matrices for each sample in \code{x}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_spillover_prepare <- function(x,
                                    spillover = NULL,
                                    select = NULL) {
  
  # EXTRACT CHANNELS
  chans <- cyto_channels(x)
  
  # EXTRACT SPILLOVER MATRIX
  if(is.null(spillover)) {
    # CYTOFRAME
    if(cyto_class(x, "flowFrame")) {
      spill <- cyto_spillover_extract(x)
      if (is.null(spill)) {
        stop("Unable to extract spillover matrix from selected sample.")
      }
      names(spill) <- cyto_names(x)
      # CYTOSET/GATINGSET
    } else {
      if (!is.null(select)) {
        spill <- cyto_spillover_extract(x[[select]])
        if (is.null(spill)) {
          stop("Unable to extract spillover matrix from selected sample.")
        }
        spill <- rep(spill, length(x))
        names(spill) <- cyto_names(x)
      } else {
        spill <- lapply(
          cyto_names(x), 
          function(y) {
            sm <- cyto_spillover_extract(x[y])[[1]]
            if (is.null(sm)) {
              stop(paste0(
                "Unable to extract spillover matrix from ",
                cyto_names(x[y]), "."
              ))
            }
            return(sm)
          }
        )
      }
    }
  } else {
    spill <- spillover
  }
  
  # READ SPILLOVER MATRIX FROM FILE
  if(cyto_class(spill, "character", TRUE)) {
    spill <- read_from_csv(spill)
  }

  # CREATE LIST OF SPILLOVER MATRICES - BYPASS CHECKING LISTS
  if(!cyto_class(spill, "list", TRUE)) {
    # NON-SQUARE OR UNLABELLED MATRIX
    if(!all(colnames(spill) %in% chans)) {
      # ANY COLUMNS CONTAINING CHANNEL NAMES?
      chans_ind <- LAPPLY(
        seq_len(ncol(spill)),
        function(z){
          if(all(spill[, z] %in% chans)) {
            return(z)
          } else {
            return(NULL)
          }
        }
      )
      # COLUMN CONTAINING SPILLOVER VALUES
      spill_ind <- LAPPLY(
        seq_len(ncol(spill)),
        function(z){
          if(is.numeric(spill[, z])) {
            return(z)
          } else {
            return(NULL)
          }
        }
      )
      # LONG MATRIX - 2 CHANNEL COLUMNS & VALUE COLUMN
      if(length(chans_ind) == 2 & length(spill_ind) >= 1) {
        # SPILLOVER CHANNELS
        spill_chans <- unique(unlist(spill[, chans_ind]))
        # INITIALISE SQUARE SPILLOVER MATRIX
        spill_mat <- diag(
          length(spill_chans)
        )
        colnames(spill_mat) <- spill_chans
        rownames(spill_mat) <- spill_chans
        # FILL MATRIX
        lapply(
          seq_len(nrow(spill)), 
          function(z){
            spill_mat[rownames(spill_mat) %in% spill[z, chans_ind[1]],
                      colnames(spill_mat) %in% spill[z, chans_ind[2]]] <<- 
              spill[z, spill_ind[1]]
          }
        )
        spill <- spill_mat
        # INVALID SPILLOVER MATRIX
      } else {
        stop(
          paste(
            "'spillover' should be a square matrix with channels",
            "as column names!"
          )
        )
      }
    }
    # SPILLOVER LIST
    spill <- rep(list(spill), length(x))
    names(spill) <- cyto_names(x)
  }
  
  # EXCLUDE MISSING PARAMETERS
  spill <- structure(
    lapply(
      spill, 
      function(z){
        # REMOVE EXCESS ROWS/COLUMNS - MATRIX MAY BE NON-SQUARE
        cols_rm <- which(!colnames(z) %in% chans)
        if(length(cols_rm) > 0) {
          z <- z[, -cols_rm]
          row_rm <- LAPPLY(
            cols_rm, 
            function(w) {
              which(z[, w] == 1)
            }
          )
          row_rm <- row_rm[!is.na(row_rm)]
          if(length(row_rm) > 0) {
            z <- z[-row_rm, ]
          }
        }
        # # PERCENTAGES -> DECIMAL - SPILLOVER < 1000%
        # if(any(z >= 10)){
        #   z <- z/100
        # }
        return(z)
      }
    ),
    names = names(spill)
  )
  
  # PREPARED SPILLOVER MATRIX
  return(spill)

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
#' @param terminal logical to indicate whether only the paths to terminal nodes
#'   should be returned, set to FALSE by default.
#' @param hidden logical to control whether the paths to hidden nodes should be
#'   included, set to FALSE by default. \code{hidden} acts in place of the
#'   \code{showHidden} argument for consistency with other CytoExploreR APIs.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:gs_get_pop_paths]{gs_get_pop_paths}} or
#'   \code{\link[openCyto:gt_get_nodes]{gt_get_nodes}}.
#'
#' @return character vector of gated node/population names.
#'
#' @importFrom flowWorkspace gh_get_pop_paths gs_get_pop_paths
#'   gh_pop_get_children gs_pop_get_children
#' @importFrom openCyto gatingTemplate gt_get_nodes
#'
#' @export
cyto_nodes <- function(x,
                       path = "full",
                       terminal = FALSE,
                       hidden = FALSE,
                       ...) {
  
  # GATINGHIERARCHY
  if (cyto_class(x, "GatingHierarchy", TRUE)) {
    nodes <- gh_get_pop_paths(
      x,
      path = path,
      showHidden = hidden,
      ...
    )
    # TERMINAL NODES
    if(terminal) {
      nodes <- LAPPLY(
        nodes,
        function(node) {
          if(length(gh_pop_get_children(x, node)) == 0) {
            return(node)
          } else {
            return(NULL)
          }
        }
      )
    }
  # GATINGSET
  } else if (cyto_class(x, "GatingSet", TRUE)) {
    nodes <- gs_get_pop_paths(
      x,
      path = path,
      showHidden = hidden,
      ...
    )
    # TERMINAL NODES
    if(terminal) {
      nodes <- LAPPLY(
        nodes,
        function(node) {
          if(length(gs_pop_get_children(x, node)) == 0) {
            return(node)
          } else {
            return(NULL)
          }
        }
      )
    }
  # GATINGTEMPLATE
  } else {
    # GATINGTEMPLATE NAME
    if(is.character(x)) {
      x <- file_ext_append(x, ".csv")
      x <- suppressMessages(gatingTemplate(x))
    }
    # NODES
    if(path == "full") {
      nodes <- names(
        gt_get_nodes(
          x, 
          only.names = TRUE,
          ...
        )
      )
      names(nodes) <- nodes
    } else {
      nodes <- gt_get_nodes(
        x, 
        only.names = TRUE,
        ...
      )
    }
    # TERMINAL NODES
    if(terminal) {
      nodes <- LAPPLY(
        seq_along(nodes),
        function(z) {
          # FULL NODE PATH REQUIRED HERE
          if(length(gt_get_children(x, names(nodes)[z])) == 0) {
            return(nodes[z])
          } else {
            return(NULL)
          }
        }
      )
    }
    nodes <- unname(nodes)
  }
  
  # NODE PATHS
  return(nodes)
  
}

## CYTO_NODES_CHECK ------------------------------------------------------------

#' Check if nodes are unique and exist in GatingSet or GatingHierarchy
#'
#' @param x object of class \code{GatingHierarchy}, \code{GatingSet},
#'   \code{gatingTemplate} or the name of a gatingTemplate CSV file.
#' @param nodes vector of node paths to check.
#' @param hidden logical to indicate whether checks should include hidden nodes
#'   as well, set to TRUE by default.
#'
#' @return supplied nodes or throw an error if any nodes do not exist or are not
#'   unique within the \code{GatingHierarchy}, \code{GatingSet} or
#'   \code{gatingTemplate}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_nodes_check <- function(x,
                             nodes = NULL,
                             hidden = TRUE){
  
  # NO NODES SUPPLIED
  if(is.null(nodes)){
    stop("Supply a vector of nodes to check to the 'nodes' argument.")
  }
  
  # REFERENCE NODE PATHS
  nodes_auto <- cyto_nodes(
    x, 
    path = "auto",
    hidden = hidden
  )
  nodes_auto_split <- .cyto_nodes_split(nodes_auto)
  
  # SPLIT NODES
  nodes_split <- .cyto_nodes_split(nodes)
  
  # CHECK NODES
  lapply(
    seq_along(nodes), 
    function(z){
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
        stop(
          paste0(
            node, 
            " is not unique in this ", class(x), ".",
            " Use either ", 
            paste0(nodes_auto[ind], collapse = " or "),
            "."
          )
        )
      }
      # NO MATCH
      if(length(ind) == 0){
        stop(paste0(node, " does not exist in this ", class(x), "."))
      }
    }
  )
  
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
#'   "auto" format, set to "auto" by default.
#' @param hidden logical indicating whether hidden nodes should be included, set
#'   to FALSE by default.
#' @param sort logical indicating whether the returned nodes should be sorted to
#'   match their order in the gating tree, set to FALSE by default.
#'
#' @return vector of unique paths for each of the supplied nodes.
#'
#' @author @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_nodes_convert <- function(x,
                               nodes = NULL,
                               anchor = NULL,
                               path = "auto",
                               hidden = FALSE,
                               sort = FALSE) {
  
  # TODO: ALLOW DIFFERENT ANCHOR PER NODE 
  
  # PATHS
  nodes_full <- cyto_nodes(
    x, 
    path = "full",
    hidden = hidden
  )
  nodes_auto <- cyto_nodes(
    x, 
    path = "auto",
    hidden = hidden
  )
  nodes_terminal <- basename(nodes_full)
  
  # STRIP REFERENCE TO ROOT
  nodes <- LAPPLY(
    nodes, 
    function(node){
      if(grepl("root/", node)){
        node <- gsub("root/", "/", node)
        return(node)
      }
      return(node)
    }
  )
  
  # ANCHOR
  if (!is.null(anchor)) {
    # INVALID ANCHOR
    if (!anchor %in% c(nodes_full, nodes_auto, nodes_terminal)) {
      stop(paste0(anchor, " does not exist in this ", class(x), "!"))
    }
    # TERMINAL ANCHOR
    if (anchor %in% nodes_terminal) {
      anchor_match <- which(
        LAPPLY(
          nodes_terminal,
          function(node_terminal) {
            anchor == node_terminal
          }
        )
      )
    # AUTO NODE
    } else if (anchor %in% nodes_auto) {
      anchor_match <- which(
        LAPPLY(
          nodes_auto, 
          function(node_auto) {
            anchor == node_auto
          }
        )
      )
    # FULL NODE
    } else if (anchor %in% nodes_full) {
      anchor_match <- which(
        LAPPLY(
          nodes_full, 
          function(node_full) {
            anchor == node_full
          }
        )
      )
    }
    # ANCHOR MUST BE UNIQUE
    if (length(anchor_match) > 1) {
      stop(paste0(anchor, " is not unique within this ", class(x), "!"))
    }
    # AUTO ANCHOR
    anchor <- nodes_auto[anchor_match]
  }
  
  # CONVERT NODES
  nodes <- LAPPLY(
    nodes, 
    function(node) {
      # TERMINAL NODE
      if (node %in% nodes_terminal) {
        nodes_match <- which(
          LAPPLY(
            nodes_terminal, 
            function(node_terminal) {
              node == node_terminal
            }
          )
        )
        # AUTO NODE
      } else if (node %in% nodes_auto) {
        nodes_match <- which(
          LAPPLY(
            nodes_auto, 
            function(node_auto) {
              node == node_auto
            }
          )
        )
        # FULL NODE
      } else if (node %in% nodes_full) {
        nodes_match <- which(
          LAPPLY(
            nodes_full,
            function(node_full) {
              node == node_full
            }
          )
        )
        # NODE DOES NOT EXIST
      } else {
        stop(paste0(node, " does not exist in this ", class(x), "!"))
      }
      # RETURN NODE
      if (length(nodes_match) == 1) {
        nodes <- cyto_nodes(
          x,
          path = path,
          hidden = hidden
        )[nodes_match]
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
          ind <- which(
            LAPPLY(
              nodes_full_split, 
              function(z) {
                all(unique(c(node_split, anchor_split)) %in% z)
              }
            )
          )
          # CANNOT FIND UNIQUE NODE PATH
          if (length(ind) == 0) {
            stop(paste0(
              "Failed to generate unique path for ", node,
              " relative to ", anchor, "."
            ))
            # UNIQUE NODE PATH EXISTS
          } else {
            if (length(ind) == 1) {
              node <- cyto_nodes(
                x, 
                path = path,
                hidden = hidden
              )[ind]
              return(node)
            } else if (length(ind) > 1) {
              # USE SHORTEST PATH (REMOVE DESCENDANTS)
              nodes_unique <- nodes_full_split[ind]
              nodes_lengths <- LAPPLY(nodes_unique, "length")
              nodes_length_min <- min(nodes_lengths)
              if (length(nodes_lengths[nodes_lengths == nodes_length_min]) > 1) {
                stop(
                  paste0(
                    "Failed to generate unique path for ", node,
                    " relative to ", anchor, "."
                  )
                )
              }
              ind <- which(nodes_lengths == nodes_length_min)
              node <- cyto_nodes(
                x, 
                path = path,
                hidden = hidden
              )[ind]
              return(node)
            }
          }
        }
      }
    }
  )
  
  # SORT NODES
  if(sort) {
    nodes <- structure(
      LAPPLY(
        nodes,
        function(node){
          if(path == "auto") {
            nodes_auto[match(node, nodes_auto)]
          } else {
            nodes_full[match(node, nodes_full)]
          }
        }
      ),
      names = names(nodes)
    )
  }
  
  # RETURN UNIQUE NODES
  return(nodes)
}

#' Split node paths into fragments
#' @noRd
.cyto_nodes_split <- function(nodes = NULL) {
  nodes_split <- lapply(
    nodes, 
    function(node) {
      node <- unlist(strsplit(node, "\\/"))
      node <- node[!LAPPLY(node, ".empty")]
      return(node)
    }
  )
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
#' @param hidden logical indicating whether hidden nodes should be included in
#'   the search for the most recent ancestral node, set to FALSE by default.
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
                                hidden = FALSE,
                                ...){
  
  # MISSING NODES
  if(is.null(nodes)){
    stop("Supply the name of the nodes to 'nodes'.")
  }
  
  # GET FULL NODES
  nodes_full <- cyto_nodes_convert(
    x,
    nodes = nodes,
    path = "full",
    hidden = hidden
  )
  
  # GET SPLIT NODES
  nodes_split <- .cyto_nodes_split(nodes_full)
  
  # GET COMMON ANCESTOR - SHORTEST SHARED PATH WITH FIRST NODE
  ancestor <- c()
  for(i in rev(seq_len(length(nodes_split[[1]])))){
    if(all(LAPPLY(nodes_split[-1], function(z){
      all(nodes_split[[1]][seq_len(i)] %in% z)
    }))){
      ancestor <- c(
        ancestor, 
        paste0(
          "/", 
          paste0(
            nodes_split[[1]][seq_len(i)],
            collapse = "/"
          )
        )
      )
      break()
    }
  }
  
  # NO COMMON ANCESTOR
  if(length(ancestor) == 0){
    ancestor  <- "root"
  }
  
  # CONVERT ANCESTRAL NODE
  ancestor <- cyto_nodes_convert(
    x,
    nodes = ancestor,
    hidden = hidden,
    ...
  )
  
  # RETURN COMMON ANCESTOR
  return(ancestor)
  
}

## CYTO_NODES_KIN --------------------------------------------------------------

#' Get paths to nodes
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param nodes vector of nodes for which the relative nodes should be returned.
#' @param type options include \code{"children"}, \code{"parent"},
#'   \code{"grandparent"} or \code{"descendants"} to indicate the type of
#'   relationship required for each of the supplied nodes.
#' @param path specifies whether the returned nodes should be in the "full" or
#'   "auto" format, set to "auto" by default.
#' @param hidden logical indicating whether hidden nodes should be included in
#'   the search for relatives, set to TRUE by default.
#' @param terminal logical to indicate whether only the relative with terminal
#'   nodes should be returned, set to FALSE by default.
#' @param ... not in use.
#'
#' @return the paths of the nodes of relation \code{type} to each of the
#'   supplied nodes.
#'
#' @importFrom flowWorkspace gh_pop_get_parent gh_pop_get_children
#'   gh_pop_get_descendants
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#' gs <- cyto_compensate(gs)
#' gs <- cyto_transform(gs)
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#'
#' # T Cells parent
#' cyto_nodes_kin(
#'   "T Cells",
#'   "parent"
#' )
#'
#' # Descendants of Live Cells
#' cyto_nodes_kin(
#'   "Live Cells",
#'   "descendants"
#' )
#'
#' @export
cyto_nodes_kin <- function(x,
                           nodes = NULL,
                           type = "children",
                           path = "auto",
                           hidden = TRUE,
                           terminal = FALSE,
                           ...) {
  
  # NOTE: HIDDEN NOT SUPPORTED FOR PARENTS IN FLOWWORKSPACE
  
  # GATINGHIERARCHY
  gh <- x
  if(cyto_class(x, "GatingSet", TRUE)) {
    gh <- x[[1]]
  }
  
  # RELATIVE NODES
  LAPPLY(
    nodes,
    function(node) {
      # CHILDREN
      if(.grepl("^c", type)) {
        pops <- gh_pop_get_children(
          gh,
          node,
          showHidden = hidden,
          path = path
        )
      # DESCENDANTS
      } else if(.grepl("^d", type)) {
        pops <- gh_pop_get_descendants(
          gh,
          node,
          showHidden = hidden,
          path = path
        )
      # PARENT
      } else if(.grepl("^p", type)) {
        pops <- gh_pop_get_parent(
          gh,
          node,
          path = path
        )
      # GRANDPARENT
      } else if(.grepl("^g", type)) {
        parent <- gh_pop_get_parent(
          gh,
          node,
          path = path
        )
        pops <- gh_pop_get_parent(
          gh,
          parent,
          path = path
        )
      # UNSUPPORTED TYPE
      } else {
        stop(
          paste0(
            "'type' must be either 'parent', 'children', ",
            "'descendants' or 'grandparent'!"
          )
        )
      }
      names(pops) <- rep(node, length(pops))
      # TERMINAL
      if(terminal) {
        pops <- pops[
          LAPPLY(
            pops,
            function(pop) {
              length(
                gh_pop_get_children(
                  gh,
                  unname(pop),
                  showHidden = hidden
                )
              ) == 0
            }
          )
        ]
      }
      return(pops)
    }
  )
  
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
      spill <- lapply(
        spill, 
        function(z) {
          z@spillover
        }
      )
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
    spill <- lapply(
      seq_along(x), 
      function(z) {
        # KEYWORDS
        kw <- keyword(x[[z]])
        # SPILLOVER SLOT
        ind <- grep("SPILL", names(kw), ignore.case = TRUE)
        if(length(ind) > 0) {
          return(kw[[ind]])
        } else {
          return(NULL)
        }
      }
    )
    names(spill) <- cyto_names(x)
    if (all(LAPPLY(spill, "is.null"))) {
      spill <- NULL
    }
  # CYTOFRAME
  } else if (cyto_class(x, "flowFrame")) {
    # KEYWORDS
    kw <- keyword(x)
    ind <- grep("SPILL", names(kw), ignore.case = TRUE)
    if(length(ind) > 0) {
      spill <- structure(kw[ind], names = cyto_names(x))
    } else {
      spill <- NULL
    }
  }
  
  # RETURN LIST OF SPILLOVER MATRICES
  return(spill)
}

## CYTO_APPLY ------------------------------------------------------------------

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
#' @param parent name of the parent(s) population to extract from
#'   \code{GatingHierarchy} or \code{GatingSet} objects, set to \code{"root"} by
#'   default.
#' @param select passed to \code{cyto_select()} to select the samples upon which
#'   the function should be applied, set to NULL by default to use all supplied
#'   samples.
#' @param coerce logical passed to \code{cyto_coerce()} to indicate whether the
#'   data should be coerced prior to applying the specified function, set to
#'   FALSE by default.
#' @param events numeric passed to \code{cyto_sample()} to control the number of
#'   events over which the specified function is to be applied, set to 1 by
#'   default to use all events.
#' @param input indicates the data input format as required by \code{FUN} can be
#'   either 1 - "cytoset", 2 - "cytoframe", 3 - "matrix", 4 - "column" or 5 -
#'   "row", set to "cytoframe" by default. \code{cyto_apply} will take care of
#'   all the data formatting prior to passing it \code{FUN}. The \code{"column"}
#'   and \code{"row"} options are for functions that expect vectors as the
#'   input.
#' @param copy logical indicating whether the data should be copied prior to
#'   preprocessing the data, set to FALSE by default. Users should set this
#'   argument to TRUE when apply inverse transformations to ensure that the
#'   original data remains unchanged.
#' @param channels vector of channels which should be included in the data
#'   passed to \code{FUN}, set to all channels by default.
#' @param trans object of class \code{transformerList} containing the
#'   definitions of the transformers applied to the supplied data. These
#'   transformers are only required when a cytoset is supplied to inverse
#'   transformations when \code{inverse = TRUE}.
#' @param inverse logical indicating whether the data should be inverse
#'   transformed prior to applying \code{FUN}, set to FALSE by default.
#' @param slot allows extraction of particular data from the output of
#'   \code{FUN}, it can be either a be the name of  a function or the name of a
#'   slot to extract from a list or an array (i.e. column name), set to NULL by
#'   default to return the entire object returned by \code{FUN}.
#' @param simplify logical indicating whether attempts should be made to coerce
#'   the output to a \code{cytoset} or \code{matrix}, set to TRUE by default.
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
                               parent = "root",
                               select = NULL,
                               coerce = FALSE,
                               events = 1,
                               input = "cytoframe",
                               copy = FALSE,
                               channels = NULL,
                               trans = NA,
                               inverse = FALSE,
                               slot = NULL,
                               simplify = TRUE) {
  
  # GATINGHIERARCHY/GATINGSET
  if(cyto_class(x,  "GatingSet")) {
    
    # TRANSFORMERS
    trans <- cyto_transformers_extract(x)
    
    # LIST OF CYTOSETS PER PARENT
    x <- cyto_data_extract(
      x,
      parent = parent,
      format = "cytoset"
    )
    
    # APPLY FUNCTION
    res <- structure(
      lapply(x, function(z){
        cyto_apply(z,
                   FUN = FUN,
                   ...,
                   select = select,
                   coerce = coerce,
                   events = events,
                   input = input,
                   copy = copy,
                   channels = channels,
                   trans = trans,
                   inverse = inverse,
                   slot = slot,
                   simplify = simplify)
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
                               parent = "root",
                               select = NULL,
                               coerce = FALSE,
                               events = 1,
                               input = "cytoframe",
                               copy = FALSE,
                               channels = NULL,
                               trans = NA,
                               inverse = FALSE,
                               slot = NULL,
                               simplify = TRUE) {
  
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
  FUN <- cyto_func_match(FUN) # namespaced function character covered
  
  # CHANNELS
  if(!is.null(channels)) {
    x <- x[, cyto_channels_extract(x, channels), drop = FALSE]
  }
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(x, select)
  }
  
  # SAMPLING/COERCION
  if(coerce) {
    x <- cyto_coerce(
      x, 
      events = events
    )
  } else {
    if(!all(events == 1)) {
      x <- cyto_sample(
        x, 
        events = events
      )
    }
  }
  
  # COPY
  if(copy) {
    x <- cyto_copy(x)
  }
  
  # INVERSE TRANSFORM
  if(!.all_na(trans) & inverse) {
    x <- cyto_transform(
      x,
      trans = trans,
      inverse = inverse,
      plot = FALSE,
      quiet = TRUE
    )
  }
  
  # CYTOSET INPUT 
  if(input == 1 | .grepl("^cytoset|^cs", input)) {
    input <- "cytoset"
    # CYTOFRAME INPUT
  } else if(input == 2 | .grepl("^cytoframe|^cf", input)) {
    input  <- "cytoframe"
    # MATRIX INPUT
  } else if(input == 3 | .grepl("^m", input)) {
    input <- "matrix"
    # COLUMN/CHANNEL INPUT
  } else if(input == 4 | .grepl("^co|^ch", input)) {
    input <- "column"
    # ROW/CELL
  } else if(input == 5 | .grepl("^r|ce", input)) {
    input <- "row"
  }
  
  # DISPATCH
  if(input == "cytoset") {
    res <- structure(
      lapply(
        seq_along(x),
        function(z){
          output <- cyto_slot(
            FUN(
              x[z],
              ...
            ),
            slot = slot
          )
          return(cyto_convert(output))
        }
      ), 
      names = cyto_names(x)
    )
  } else if(input == "flowFrame") { # internal use only
    res <- structure(
      lapply(
        cyto_names(x), 
        function(z){
          output <- cyto_slot(
            FUN(
              cytoframe_to_flowFrame(x[[z]]), 
              ...
            ),
            slot = slot
          )
          return(cyto_convert(output))
        }
      ), 
      names = cyto_names(x)
    )
  } else if(input == "cytoframe") {
    res <- structure(
      lapply(
        cyto_names(x),
        function(z){
          output <- cyto_slot(
            FUN(
              x[[z]], 
              ...
            ),
            slot = slot
          )
          return(cyto_convert(output))
        }
      ), 
      names = cyto_names(x)
    )
  } else if(input == "matrix") {
    res <- structure(
      lapply(
        cyto_names(x), 
        function(z){
          output <- cyto_slot(
            FUN(
              cyto_exprs(x[[z]], drop = FALSE),
              ...
            ),
            slot = slot
          )
          return(cyto_convert(output))
        }
      ), 
      names = cyto_names(x)
    )
  } else if(input == "column") {
    # TODO: Add support for passing channel-specific arguments through here
    # named vector or named list
    res <- structure(
      lapply(
        cyto_names(x), 
        function(z){
          output <- cyto_slot(
            apply(
              cyto_exprs(x[[z]], drop = FALSE),
              2, 
              FUN, 
              ...
            ),
            slot = slot
          )
          return(cyto_convert(output))
        }
      ), 
      names = cyto_names(x)
    )
  } else if(input =="row") {
    res <- structure(
      lapply(
        cyto_names(x), 
        function(z){
          output <- cyto_slot(
            apply(
              cyto_exprs(x[[z]], drop = FALSE),
              1, 
              FUN,
              ...
            ),
            slot = slot
          )
          return(cyto_convert(output))
        }
      ), 
      names = cyto_names(x)
    )
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
      # VECTORS TO MATRIX - CANNOT USE DIM HERE - S4 OBJECTS NO DIMS
      if(cyto_class(res[[1]], 
                    c("integer", "numeric", "logical", "character"), TRUE)) {
        # DONT MERGE CHANNELS INTO SINGLE COLUMN
        if(!all(names(res[[1]]) %in% cyto_channels(x))) {
          # ATTEMPT MATRIX CONVERSION - VECTORS DIFFERENT OF SIZES
          res <- structure(
            lapply(
              names(res),
              function(z){
                tryCatch(
                  matrix(
                    res[[z]],
                    nrow = length(res[[z]]),
                    ncol = 1,
                    dimnames = list(
                      paste0(
                        names(res[z]), 
                        if(!is.null(names(res[[z]])) & length(res[[z]]) > 1) {
                          paste0("|", names(res[[z]]))
                        } else {
                          ""
                        }
                      ),
                      FUN_NAME
                    )
                  ),
                  error = function(e) {
                    return(res[[z]])
                  }
                )
              }
            ),
            names = names(res)
          )
        } else {
          res <- structure(
            lapply(
              names(res),
              function(z){
                matrix(
                  res[[z]],
                  ncol = length(res[[z]]),
                  nrow = 1,
                  dimnames = list(
                    names(res[z]),
                    names(res[[z]])
                  )
                )
              }
            ),
            names = names(res)
          )
        }
      }
      # LIST OF MATRICES
      if(all(!is.null(LAPPLY(res, "dim")))) {
        # PREPARE & FORMAT MATRICES
        res <- lapply(
          names(res), 
          function(z){
            # ROWNAMES
            if(is.null(rownames(res[[z]]))) {
              if(nrow(res[[z]]) > 1) {
                rownames(res[[z]]) <- paste(
                  z, 
                  1:nrow(res[[z]]), 
                  sep = "|"
                )
              } else {
                rownames(res[[z]]) <- z
              }
            } else {
              # CYTO_NAMES MISSING IN OUTPUT - CHECK FIRST ENTRY ONLY
              if(!grepl(z, rownames(res[[z]])[1], fixed = TRUE)) {
                rownames(res[[z]]) <- paste(
                  z, 
                  rownames(res[[z]]), 
                  sep = "|"
                )
              }
            }
            # COLNAMES
            if(is.null(colnames(res[[z]]))) {
              if(ncol(res[[z]]) == 1) {
                colnames(res[[z]]) <- FUN_NAME
              } else {
                colnames(res[[z]]) <- paste0(FUN_NAME, "-", 1:ncol(res[[z]]))
              }
            }
            return(res[[z]])
          })
        # RBIND MATRICES - SAME DIMENSIONS
        if(length(unique(LAPPLY(res, "ncol"))) == 1) {
          res <- do.call("rbind", unname(res))
        }
      }
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
                            parent = "root",
                            select = NULL,
                            coerce = FALSE,
                            events = 1,
                            input = "cytoframe",
                            copy = FALSE,
                            channels = NULL,
                            trans = NA,
                            inverse = FALSE,
                            slot = NULL,
                            simplify = TRUE) {
  
  structure(
    lapply(
      x, 
      function(z) {
        cyto_apply(
          z,
          FUN = FUN,
          ...,
          parent = parent,
          select = select,
          coerce = coerce,
          events = events,
          input = input,
          copy = copy,
          channels = channels,
          trans = trans,
          inverse = inverse,
          slot = slot,
          simplify = simplify
        )
      }
    ), 
    names = names(x)
  )
  
}

## CYTO_SLOT -------------------------------------------------------------------

#' Accessory function to extract elements from objects
#'
#' @param x an object from which element(s) are to be extracted.
#' @param slot can be either the name of a slot to extract or a function that
#'   converts \code{x} to the desired output.
#'
#' @return the element extracted or formatted from \code{x}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#' library(CytoExploreRData)
#' 
#' # Create named list
#' res <- list("A" = 3, "B" = 4)
#' 
#' # Extract B from list
#' cyto_slot(res, "B")
#' }
#'
#' @export
cyto_slot <- function(x,
                      slot = NULL) {
  
  # NO SLOT
  if(is.null(slot)) {
    return(x)
  }
  
  # TRY MATCH A FUNCTION 
  slot <- tryCatch(
    cyto_func_match(
      slot
    ),
    error = function(e) {
      return(slot)
    }
  )
  
  # EXTRACT SLOT USING FUNCTION
  if(cyto_class(slot, "function")) {
    return(
      cyto_func_call(
        slot,
        list(x)
      )
    )
  # EXTRACT SLOT BY NAME
  } else {
    # LIST OR VECTOR
    if(is.null(dim(x))) {
      if(slot %in% names(x)) {
        return(
          x[[slot]]
        )
      } else {
        warning(
          paste0(
            slot, 
            " does not exist in this ", 
            cyto_class(x)[1], 
            "!"
          )
        )
        return(x)
      }
    # ARRAY
    } else {
      if(slot %in% colnames(x)) {
        return(
          x[, slot]
        )
      } else {
        warning(
          paste0(
            slot, 
            " does not exist in the column names of this ", 
            cyto_class(x)[1], 
            "!"
          )
        )
        return(x)
      }
    }
  }
  
}

## CYTO_CONVERT ----------------------------------------------------------------

#' Convert flowFrame/flowSet objects to cytoframe/cytoset objects
#'
#' @param x object of class \code{flowFrame}, \code{flowSet}, \code{cytoframe}
#'   or \code{cytoset}.
#' @param ... additional arguments passed to the appropriate
#'   \code{\link[flowWorkspace:convert]{conversion}} function.
#'
#' @return a cytoframe or cytoset.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace load_cytoframe_from_fcs flowSet_to_cytoset
#'
#' @seealso \code{\link[flowWorkspace:convert]{flowFrame_to_cytoframe}}
#' @seealso \code{\link[flowWorkspace:convert]{flowSet_to_cytoset}}
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
    # # FLOWFRAME_TO_CYTOFRAME() BUT RETAIN ORIGINAL FILENAME
    # tmp <- paste0(tempdir, .Platform$file.sep, cyto_names(x))
    # write.FCS(x, tmp)
    # x <- load_cytoframe_from_fcs(tmp, ...)
    x <- flowFrame_to_cytoframe(x, emptyValue = FALSE, ...)
  } else if(cyto_class(x, "flowSet", TRUE)) {
    x <- flowSet_to_cytoset(x, ...)
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
                       cols = NULL){
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
    cyto_counts <- cyto_apply(
      x, 
      "nrow", 
      input = "matrix",
      copy = FALSE
    )
    # SAME NUMBER OF EVENTS
    if(nrow(cols) != sum(cyto_counts)){
      stop(
        paste0("'cols' does not contain the same number of events as this", 
               class(x))
      )
      # SPLIT MATRIX INTO LIST
    }else{
      cols <- lapply(
        seq_along(cyto_counts), 
        function(z){
          if(z == 1){
            cols[1:cyto_counts[z], , drop = FALSE]
          }else{
            start <- sum(cyto_counts[1:(z-1)]) + 1
            end <- start + cyto_counts[z] - 1
            cols[start:end, , drop = FALSE]
          }
        }
      )
      names(cols) <- cyto_names(x)
    }
  } else {
    stop("'cols' must be a matrix!")
  }
  
  # BIND COLUMNS
  cf_list <- lapply(
    cyto_names(x),
    function(z){
      cyto_cbind(x[[z]], cols[[z]])
    }
  )
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
        cyto_func_call("remotes::install_github", list(repo, ...))
        # remotes::install_github(repo, ...)
    # CRAN
    } else if(grepl("^c", source, ignore.case = TRUE)){
      install.packages(x, ...)
    # BIOCONDUCTOR
    } else if (grepl("^b", source, ignore.case = TRUE)) {
      if(!"BiocManager" %in% pkgs$Package) {
        install.packages("BiocManager", ...)
      }
        cyto_func_call("BiocManager::install", list(x, ...))
        # BiocManager::install(x, ...)
    }
  }
  
  # LOAD PACKAGE
  requireNamespace(x)
  
  # REFERENCE
  if(!is.null(ref)) {
    message(ref)
  }
  
}

## CYTO_PROGRESS ---------------------------------------------------------------

#' CytoExploreR progress bar
#'
#' @param pb an existing progress bar to increment.
#' @param label text to include to the left of the progress bar.
#' @param total the maximum number of iterations for the progress bar.
#' @param clear logical indicating whether the progress bar should be removed
#'   from the console once complete, set to FALSE by default.
#' @param ... additional arguments passed to \code{progress_bar}.
#'
#' @return progress bar object.
#'
#' @importFrom progress progress_bar
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#'   pb <- cyto_progress(total = 10)
#'   lapply(1:10, function(z) {
#'     cyto_progress(pb)
#'     Sys.sleep(0.2)
#'   })
#' }
#' 
#' @export
cyto_progress <- function(pb = NULL,
                          label = NULL,
                          total = 100,
                          clear = FALSE,
                          ...) {
  
  # ITERATE EXISTING PROGRESS BAR
  if(!is.null(pb)) {
    # INCREMENT PROGRESS BAR
    if(!pb$finished) {
      pb$tick()
    }
    # PROGRESS BAR COMPLETE
    if(pb$finished) {
      cyto_option("CytoExploreR_progress", NULL)
    }
  } else {
    # CREATE NEW PROGRESS BAR - NO ACTIVE PROGRESS BARS
    if(is.null(cyto_option("CytoExploreR_progress"))) {
      pb <- progress_bar$new(
        format = paste0(
          "",
          label,
          " [:bar] :percent ETA::eta"
        ),
        total = total,
        clear = clear,
        show_after = 0,
        ...
      )
      pb$tick(0)
      cyto_option(
        "CytoExploreR_progress", 
        structure(
          list(pb), 
          names = label
        )
      ) 
    }
  }
  return(pb)
  
}
