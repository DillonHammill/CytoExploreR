## CYTO_GATE_CLEAN -------------------------------------------------------------

#' Apply anomaly detection algorithms to gate high quality cytometry data
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parental population to sample, set to the
#'   \code{"root"} node by default.
#' @param alias name of the population containing high quality events to be
#'   created in the GatingSet, set to \code{"clean"} by default.
#' @param type the method to use when cleaning the data, options include
#'   \code{"flowAI"}, \code{"flowClean"}, \code{"flowCut"} or \code{"PeacoQC"},
#'   set to \code{"flowAI"} by default. Custom functions are also supported
#'   through \code{type} if they return a logical vector indicating the events
#'   to gate. The \code{slot} argument provides flexibility over how the
#'   resulting logical vector is obtained from the object returned by the custom
#'   gating function.
#' @param channels vector of channels/markers over which anomaly detection
#'   algorithms should be applied to gate high quality events, set to all
#'   channels except \code{"Event-ID"} and \code{"Sample-ID"} by default. The
#'   \code{"Time"} channel should be included here when performing flow rate and
#'   signal acquisition checks.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'   \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'   saved, set to \code{cyto_gatingTemplate_active()} by default.
#' @param input allows flexibility over how the data should be formatted prior
#'   to passing it a custom gating function, set to \code{flowFrame} by default
#'   to match the current behaviour in \code{openCyto}.
#' @param slot provides flexibility over how the resulting logical vector is
#'   obtained from the object returned by the custom gating function.
#' @param trans object of class transformerList containing the definitions of
#'   transformers already applied to the supplied data, only required for
#'   \code{cytoset} objects.
#' @param inverse logical indicating whether inverse data transformations should
#'   be applied to the data prior to cleaning, set to FALSE by default.
#' @param ... additional arguments to passed to the gating function.
#'
#' @importFrom openCyto gs_add_gating_method
#'
#' @return a \code{GatingHierarchy} or \code{GatingSet} with new node containing
#'   high quality events and an updated \code{gatingTemplate}.
#'
#' @seealso \code{\link{cyto_slot}}
#' @seealso \code{\link{cyto_clean}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Prepare Activation GatingSet
#' gs <- GatingSet(Activation)
#' gs <- cyto_compensate(gs)
#' gs <- cyto_transform(gs)
#'
#' # Write gatingTemplate to file
#' cyto_gatingTemplate_create("gatingTemplate.csv", active = TRUE)
#'
#' # Gate high quality events
#' cyto_gate_clean(
#'   gs,
#'   parent = "root",
#'   alias = "clean",
#'   type = "PeacoQC",
#'   gatingTemplate = "gatingTemplate.csv"
#' )
#' }
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
#' @export
cyto_gate_clean <- function(x,
                            parent = "root",
                            alias = "clean",
                            type = "PeacoQC",
                            channels = NULL,
                            gatingTemplate = NULL,
                            input = "flowFrame",
                            slot = NULL,
                            trans = NA,
                            inverse = FALSE,
                            ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET REQUIRED
  if(!cyto_class(x, "GatingSet")) {
    stop(
      paste0(
        "'cyto_gate_sample()' only supports GatingHierarchy and GatingSet ",
        "objects!"
      )
    )
  }
  
  # ALIAS
  if(length(alias) != 1) {
    stop(
      "Supply a name for the sampled population to 'alias'."
    )
  }
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x,
      exclude = c("Event", "Sample")
    )
  } else {
    channels <- cyto_channels_extract(
      x, 
      channels
    )
  }
  
  # TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # GATINGTEMPLATE CHECKS
  if(cyto_class(x, "GatingSet")) {
    # ACTIVE GATINGTEMPLATE
    if (is.null(gatingTemplate)) {
      gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
    }
    # CHECK EXISTING ENTRIES IN GATINGTEMPLATE
    gt <- .cyto_gatingTemplate_check(
      parent, 
      alias, 
      gatingTemplate
    )
    # CREATE GATINGTEMPLATE
    if (is.null(gt)) {
      message(
        paste("Creating", gatingTemplate, "to save the constructed gate(s).")
      )
      cyto_gatingTemplate_create(gatingTemplate, active = TRUE)
      gt <- cyto_gatingTemplate_read(gatingTemplate, data.table = TRUE)
    }
  }
  
  # ALGORITHM REFERNCES --------------------------------------------------------
  
  # GATING MESSAGE
  message(
    paste0(
      "Applying ", cyto_func_name(type), "() ", 
      "anomaly detection algorithm to gate high ",
      "quality events in the ", parent, " population...\n"
    )
  )
  
  # DEFAULT ALGORITHM TYPES
  if(is.character(type)) {
    # FLOWAI
    if(grepl("^FlowAI$", type, ignore.case = TRUE)) {
      # LOAD FLOWAI
      cyto_require(
        "flowAI",
        source = "BioC",
        repo = "giannimonaco/flowAI",
        ref = paste0(
          "Monaco G, Chen H, Poidinger M, Chen J,",
          " de Magalhaes J, Larbi A (2016). flowAI:",
          " automatic and interactive anomaly discerning",
          " tools for flow cytometry data. Bioinformatics,",
          " 32(16)."
        ),
        version = "1.27.3"
      )
      input <- "flowFrame"
      inverse <- TRUE
    # FLOWCUT
    } else if(grepl("^FlowCut$", type, ignore.case = TRUE)) {
      # LOAD FLOWCUT
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
      input <- "flowFrame"
      inverse <- FALSE
    # FLOWCLEAN
    } else if(grepl("^FlowClean$", type, ignore.case = TRUE)) {
      # LOAD FLOWCLEAN
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
      input <- "flowFrame"
      inverse <- TRUE
    # PEACOQC
    } else if(grepl("^PeacoQC$", type, ignore.case = TRUE)) {
      # LOAD PEACOQC
      cyto_require(
        "PeacoQC",
        source = "BioC",
        repo = "saeyslab/PeacoQC",
        ref = paste0(
          "Emmaneel A, et al. (2021) PeacoQC: peak-based selection of high ",
          "quality cytometry data. Cytometry A."
        )
      )
      input <- "flowFrame"
      inverse <- FALSE
    }
  }
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # APPLY GATES TO GATINGSET & CREATE GATINGTEMPLATE ENTRY
  pop <- suppressWarnings(
    suppressMessages(
      gs_add_gating_method(
        gs = x,
        alias = alias,
        parent = parent,
        pop = "*",
        dims = "", # EMPTY DIM WARNING
        gating_method = "cyto_gate_clean",
        gating_args = list(
          "type" = type,
          "input" = input,
          "params" = channels,
          "slot" = slot,
          "inverse" = inverse,
          "openCyto.minEvents" = -1,
          ...
        ),
        groupBy = NA,
        collapseDataForGating = FALSE,
        preprocessing_method = "pp_cyto_gate_clean",
        preprocessing_args =  NA
      )
    )
  )
  
  # ADD POPULATIONS TO GATINGTEMPLATE
  gt <- rbind(gt, pop)
  
  # WRITING NEW GATINGTEMPLATE ENTRIES
  message(paste("Re-writing", gatingTemplate, "with new gating entries..."))
  
  # SAVE UPDATED GATINGTEMPLATE
  cyto_gatingTemplate_write(gt, gatingTemplate)
  
  # RETURN GATINGSET
  return(x)
  
}

#' Preprocessing method to pass transformers to gating function
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.pp_cyto_gate_clean <- function(fs,
                                gs,
                                gm,
                                channels,
                                groupBy = NA,
                                isCollapse = NA,
                                ...) {
  
  return(
    list(
      cyto_transformers_extract(gs)
    )
  )
  
}

#' Anomaly detection algorithm plugin for openCyto
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.cyto_gate_clean <- function(fr,
                             pp_res,
                             channels,
                             type = "PeacoQC",
                             input = "flowFrame",
                             params = NULL,
                             slot = NULL,
                             inverse = FALSE,
                             ...) {
  
  # BYPASS EMPTY CYTOFRAMES
  if(nrow(fr) == 0) {
    return(logical(0))
  }
  
  # CHANNELS FOR GATING PASSED THROUGH PARAMS ARGUMENT
  if(is.null(params)) {
    params <- cyto_channels(
      fr,
      exclude = c(
        "Event",
        "Sample"
      )
    )
  }
  
  # TODO: MORE FLEXIBLE MATCHING TO INPUT
  # CONVERT DATA TO REQUIRED FORMAT - DATA RESTRICTED TO CHANNELS
  if(!cyto_class(fr, input, TRUE)) {
    fr <- cyto_data_extract(
      fr,
      format = input,
      channels = params,
      trans = pp_res$trans,
      inverse = inverse,
      copy = TRUE
    )[[1]][[1]]
  }
  
  # DISPATCH BASED ON TYPE - DEFAULT METHODS
  if(is.character(type)) {
    # FLOWAI
    if(grepl("^FlowAI$", type, ignore.case = TRUE)) {
      # FLOWAI ARGUMENTS - DEFAULTS
      args_default <- list(
        output = 3,
        html_report = FALSE,
        mini_report = FALSE,
        fcs_QC = FALSE,
        folder_results = FALSE,
        emptyValue = FALSE
      )
      # COMBINE INCOMING ARGUMENTS
      args <- list(...)
      args <- c(
        list(fr),
        args,
        args_default[!names(args_default) %in% names(args)]
      )
      # RUN FLOWAI - REMOVAL INDICES
      remove <- suppressPrint(
        cyto_func_call(
          "flowAI::flow_auto_qc",
          args
        )[[1]]
      )
      # LOGICAL VECTOR
      gate <- rep(TRUE, nrow(fr))
      if(length(remove) > 0) {
        gate[remove] <- FALSE
      }
      return(gate)
    # FLOWCUT
    } else if(grepl("^FlowCut$", type, ignore.case = TRUE)) {
      # FLOWCUT ARGUMENTS - DEFAULTS
      args_default <- list(
        Plot = "None"
      )
      # COMBINE INCOMING ARGUMENTS
      args <- list(...)
      args <- c(
        list(fr),
        args,
        args_default[!names(args_default) %in% names(args)]
      )
      # RUN FLOWAI - REMOVAL INDICES
      remove <- invisible(
        capture.output(
          cyto_func_call(
            "flowCut::flowCut",
            args
          )[["ind"]]
        )
      )
      # LOGICAL VECTOR
      gate <- rep(TRUE, nrow(fr))
      if(length(remove) > 0) {
        gate[remove] <- FALSE
      }
      return(gate)
    # FLOWCLEAN
    } else if(grepl("^FlowClean", type, ignore.case = TRUE)) {
      # FLOWCUT ARGUMENTS - DEFAULTS
      args_default <- list(
        filePrefixWithDir = cyto_names(fr),
        ext = ".fcs",
        diagnostic = FALSE
      )
      # COMBINE INCOMING ARGUMENTS
      args <- list(...)
      args <- c(
        list(fr),
        args,
        args_default[!names(args_default) %in% names(args)]
      )
      # RUN FLOWCLEAN
      gate <- invisible(
        capture.output(
          cyto_exprs(
            cyto_func_call(
              "flowClean::clean",
              args
            ),
            "GoodVsBad",
            drop = TRUE
          ) < 10000
        )
      )
      # LOGICAL VECTOR
      return(gate)
    # PEACOQC
    } else if(grepl("^PeacoQC", type, ignore.case = TRUE)) {
      # PEACOQC ARGUMENTS - DEFAULTS
      args_default <- list(
        output = "full",  # MARGIN INDICES
        channels = cyto_channels(fr),
        plot = FALSE,
        save_fcs = FALSE
      )
      # COMBINE INCOMING ARGUMENTS
      args <- list(...)
      args <- c(
        args,
        args_default[!names(args_default) %in% names(args)]
      )
      # EVENT INDICES
      ind <- seq_len(nrow(fr))
      # REMOVE MARGINS - LIST(FR, INDICES)
      res <- cyto_func_execute(
        "PeacoQC::RemoveMargins",
        args = c(
          list("ff" = fr),
          args
        )
      )
      # REMOVE MARGIN EVENTS
      if(length(res[[2]]) > 0) {
        ind <- ind[-c(res[[2]])]
      }
      # RUN PEACOQC
      keep <- cyto_func_execute(
        "PeacoQC::PeacoQC",
        args = c(
          list("ff" = res[[1]]),
          args
        )
      )[["GoodCells"]]
      # INDICES TO KEEP
      if(length(keep) > 0) {
        ind <- ind[keep]
      }
      gate <- rep(FALSE, nrow(fr))
      if(length(ind) > 0) {
        gate[ind] <- TRUE
      }
      return(gate)
    }
  }
  
  # NON-DEFAULT TYPE
  gate <- cyto_slot(
    cyto_func_call(
      type,
      list(fr, ...)
    ),
    slot = slot
  )
  
  # RETURN LOGICAL VECTOR
  return(gate)
  
}
