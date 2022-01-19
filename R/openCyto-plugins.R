## CYTO_GATE_DRAW --------------------------------------------------------------

#' Manual Gate Drawing Plugin for openCyto
#'
#' \code{gate_manual} if a wrapper for gate_draw which allows the user to manual
#' draw gates using the \code{openCyto} gating pipeline. This plugin lacks the
#' ability to save drawn gates to the gatingTemplate, this feature is however
#' included in \code{cyto_gate_draw} which invisibly returns these
#' gatingTemplate entries but operates independently of the openCyto gating
#' pipeline.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for gating.
#' @param pp_res output of preprocessing function.
#' @param channels name(s) of the channel(s) to be used for plotting.
#' @param alias name of the population to be gated. This is not inherited from
#'   the gatingTemplate and must be supplied manually to the gating_args.
#' @param ... additional arguments passsed to \code{\link{cyto_plot}}.
#'
#' @return a \code{filters} object containing the contructed gate objects which
#'   is passed onto the \code{openCyto} gating pipeline to apply these gates to
#'   the GatingSet.
#'
#' @importFrom openCyto registerPlugins
#'
#' @seealso \code{\link{cyto_plot}}
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#' fs <- Activation # load in .fcs files
#'
#' gs <- GatingSet(fs) # add flowSet to GatingSet
#'
#' template <- add_pop(
#'   gs,
#'   alias = "Lymphocytes", pop = "+", parent = "root",
#'   dims = "FSC-A,SSC-A", gating_method = "cyto_gate_manual",
#'   gating_args = "events=0.5,alias='Lymphocytes',type='ellipse'",
#'   collapseDataForGating = TRUE, groupBy = 2
#' )
#'
#' # gating window will open to construct gate left click vertices on plot
#' cyto_plot(gs[[1]],
#'   parent = "root",
#'   alias = "Lymphocytes",
#'   channels = c("FSC-A", "SSC-A")
#' )
#' }
#'
#' @noRd
.cyto_gate_manual <- function(fr,
                              pp_res,
                              channels,
                              alias, 
                              ...) {
  
  # CALL CYTO_GATE_DRAW
  cyto_gate_draw(
    x = fr, 
    channels = channels, 
    alias = alias,
    ...
  )
  
}

## CYTO_GATE_MATCH -------------------------------------------------------------

#' Plugin function to extract gates stored in openCyto gatingTemplates
#'
#' @param fr object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}}
#'   containing data for a group.
#' @param pp_res return value from a preprocessing function, in this case it
#'   should be paired with \code{.pp_cyto_gate_match()} which returns the index
#'   or name of the set of gates to extract from the gatingTemplate.
#' @param gate passed manually as a list of filters to the openCyto
#'   gatingTemplate. openCyto will passed the prepared gates to this function
#'   for further processing.
#'   
#' @return a list of filters to be applied to this group.
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_gate_match <- function(fr,
                             pp_res,
                             channels,
                             gate) {
  
  # PP_RES = NULL - NO GROUPING
  if(is.null(pp_res)) {
    ind <- 1
  }
  
  # GATE INDEX
  if(is.numeric(pp_res)) {
    ind <- pp_res
  } else if(is.character(pp_res)) {
    ind <- match(
      pp_res,
      names(gate)
    )
  }
  
  # EXTRACT GATE
  return(gate[[ind]])
  
}

#' Preprocessing openCyto plugin pass group index/name to gating function
#'
#' @param fs an object of class
#'   \code{\link[flowWorkspace:cytoset-class]{cytoset}} containing the samples
#'   for this group.
#' @param gs an object of class \code{\link[flowWorkspace:GatingSet-class]}
#'   containing data for all samples, including samples within \code{fs}.
#' @param gm name of the gating method to use, not in use.
#' @param groupBy specifies how the samples have been grouped based experiment
#'   variables. Variables should be separated by a colon and index grouping is
#'   not supported due to inconsistencies in gating.
#' @param isCollapse used by openCyto as a means of passing
#'   \code{CollpaseDataForGating} flag to the gating method.
#' @param ... not in use.
#' 
#' @return index or name of the group to be gated.
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.pp_cyto_gate_match <- function(fs,
                                gs,
                                gm,
                                channels,
                                groupBy = NA,
                                isCollapse = NA,
                                ...) {
  
  # GROUPBY
  if(is.character(groupBy) & nchar(groupBy) == 0) {
    groupBy <- NA
  }
  
  # GROUPING VARIABLES SEPARATED BY COLON
  if(is.character(groupBy)) {
    groupBy <- unlist(strsplit(groupBy, ":"))
  }
  
  # EXPERIMENT GROUPS
  if(!.all_na(groupBy)) {
    # VARIABLE NAMES
    if(all(grepl("^[A-Za-z]+$", groupBy, ignore.case = TRUE))) {
      # MATCH NAMES FS TO GROUPS IN GS
      grps <- cyto_groups(
        gs,
        group_by = groupBy,
        details = TRUE
      )
      # WE COULD ALSO PASS GROUP NAMES HERE
      ind <- names(grps)[
        LAPPLY(
          grps,
          function(grp) {
            all(cyto_names(fs) %in% c(rownames(grp), grp$name))
          }
        )
      ]
    # INDICES NOT SUPPORTED - GATE ORDER INCONSISTENT
    } else {
      if(groupBy == length(gs)) {
        ind <- 1
      } else {
        message(
          "Numeric groupBy is not supported. Grouping all samples."
        )
        ind <- 1
      }
    }
  }
  
  # NO GROUPING - USE FIRST SET OF GATES
  if(.all_na(groupBy)) {
    ind <- 1
  }
  
  # RETURN INDEX OF GATE OR GROUP NAME (PREFERRED)
  return(ind)
  
}

## CYTO_GATE_DRAW --------------------------------------------------------------

#' Manual gating openCyto plugin
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.cyto_gate_draw <- function(fr,
                            pp_res,
                            channels,
                            gate) {
  
  cyto_func_call(
    ".cyto_gate_match",
    .args_list()
  )

}

#' Manual gating preprocessing openCyto plugin
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.pp_cyto_gate_draw <- function(fs,
                               gs,
                               gm,
                               channels,
                               groupBy = NA,
                               isCollapse = NA,
                               ...) {
  
  cyto_func_call(
    ".pp_cyto_gate_match",
    .args_list(...)
  )
  
}

## CYTO_GATE_SAMPLE ------------------------------------------------------------

#' Sampling plugin gating function for openCyto
#'
#' \code{.pp_cyto_gate_sample()} computes the number of events to extract from
#' each sample and returns a logical vector for these events in the merged
#' cytoframe passed to \code{.cyto_gate_sample()}.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_gate_sample <- function(fr,
                              pp_res,
                              channels) {
  
  # OPENCYTO WILL PASS LOGICAL VECTOR FROM PREPROCESSING
  return(pp_res)

}

#' Sampling preprocessing function for openCyto
#'
#' See above for details.
#'
#' @importFrom flowWorkspace gh_pop_get_indices
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
#' @noRd
.pp_cyto_gate_sample <- function(fs,
                                 gs,
                                 gm,
                                 channels,
                                 groupBy = NA,
                                 isCollapse = NA,
                                 events = 1,
                                 seed = NULL,
                                 ...) {
    
  # GROUPBY
  if(is.character(groupBy) & nchar(groupBy) == 0) {
    groupBy <- NA
  }
  
  # GROUPING VARIABLES SEPARATED BY COLON
  if(is.character(groupBy)) {
    groupBy <- unlist(strsplit(groupBy, ":"))
  }
  
  # NO GROUPBY
  if(.all_na(groupBy)) {
    groupBy <- "name"
  }
  
  # GATING GROUPS
  grps <- cyto_groups(
    gs,
    group_by = groupBy,
    details = TRUE
  )
  
  # EVENTS PER GROUP
  events <- rep(events, length.out = length(grps))
  names(events) <- names(grps)
  
  # WE COULD ALSO PASS GROUP NAMES HERE
  grp_ind <- names(grps)[
    LAPPLY(
      grps,
      function(grp) {
        all(cyto_names(fs) %in% c(rownames(grp), grp$name))
      }
    )
  ]
  
  # COMPUTE EVENTS - THIS GROUP
  grp_events <- cyto_sample_n(
    fs,
    events = events[grp_ind]
  )
  
  # SEED
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # LOGICAL FILTER
  gate <- cyto_apply(
    fs, 
    input = "cytoset",
    FUN = function(cs) {
      # LOGICAL SAMPLE FILTER
      if(nrow(cs[[1]]) == 0) {
        return(logical(0))
      } else {
        keep <- rep(FALSE, nrow(cs[[1]]))
        keep[
          sample(
            seq_len(nrow(cs[[1]])),
            grp_events[cyto_names(cs)]
          )
        ] <- TRUE
        return(keep)
      }
    },
    copy = FALSE,
    simplify = FALSE
  )
  
  # RETURN LOGICAL SAMPLE FILTER
  return(gate)
  
}

## CYTO_GATE_CLEAN -------------------------------------------------------------

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
                             ...) {
  
  # BYPASS EMPTY CYTOFRAMES
  if(nrow(fr) == 0) {
    return(logical(0))
  }
  
  # CHANNELS FOR GATING PASSED THROUGH PARAMS ARGUMENT
  if(is.null(params)) {
    params <- cyto_channels(
      fr,
      exclude = c("Event", "Sample", "Time")
    )
  }
  
  # CONVERT DATA TO REQUIRED FORMAT - DATA RESTRICTED TO CHANNELS
  if(!cyto_class(fr, input, TRUE)) {
    fr <- cyto_data_extract(
      fr,
      format = input,
      channels = params,
      copy = FALSE
    )[[1]][[1]]
  }
  
  # DISPATCH BASED ON TYPE - DEFAULT METHODS
  if(is.character(type)) {
    # FLOWAI
    if(grepl("^FlowAI$", type, ignore.case = TRUE)) {
      # LOAD FLOWAI
      cyto_require(
        "flowAI",
        source = "BioC",
        repo = "giannimonaco/flowAI",
        ref = NULL # DROP REFERENCE - PRINTED MULTIPLE TIMES
      )
      # FLOWAI ARGUMENTS - DEFAULTS
      args_default <- list(
        output = 3,
        html_report = FALSE,
        mini_report = FALSE,
        fcs_QC = FALSE,
        folder_results = FALSE
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
            "flowAI::flow_auto_qc",
            args
          )[[1]]
        )
      )
      # LOGICAL VECTOR
      gate <- rep(TRUE, nrow(fr))
      if(length(remove) > 0) {
        gate[remove] <- FALSE
      }
      return(gate)
      # FLOWCUT
    } else if(grepl("^FlowCut$", type, ignore.case = TRUE)) {
      # LOAD FLOWCUT
      cyto_require(
        "flowCut",
        source = "BioC",
        repo = "jmeskas/flowCut",
        ref = NULL
      )
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
      # LOAD FLOWCLEAN
      cyto_require(
        "flowClean",
        source = "BioC",
        repo = "cafletezbrant/flowClean",
        ref = NULL
      )
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
      # LOAD PEACOQC
      cyto_require(
        "PeacoQC",
        source = "BioC",
        repo = "saeyslab/PeacoQC",
        ref = NULL
      )
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
        c(list("ff" = fr),
          args)
      )
      # REMOVE MARGIN EVENTS
      if(length(res[[2]]) > 0) {
        ind <- ind[-c(res[[2]])]
      }
      # RUN PEACOQC
      keep <- cyto_func_execute(
        "PeacoQC::PeacoQC",
        c(list("ff" = res[[1]]),
          args)
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
