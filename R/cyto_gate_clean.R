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
#'   channels except \code{"Time"}, \code{"Event-ID"} and \code{"Sample-ID"} by
#'   default.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'   \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'   saved, set to \code{cyto_gatingTemplate_active()} by default.
#' @param input allows flexibility over how the data should be formatted prior
#'   to passing it a custom gating function, set to \code{flowFrame} by default
#'   to match the current behaviour in \code{openCyto}.
#' @param slot provides flexibility over how the resulting logical vector is
#'   obtained from the object returned by the custom gating function.
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
                            ...) {
  
  # TODO: DISPLAY REFERENCES HERE?
  
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
      exclude = c("Event", "Sample", "Time")
    )
  } else {
    channels <- cyto_channels_extract(
      x, 
      channels
    )
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
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # MESSAGE
  message(
    "Applying anomaly detection algorithms to gate high quality events..."
  )

  # APPLY GATES TO GATINGSET & CREATE GATINGTEMPLATE ENTRY
  pop <- suppressWarnings(
    suppressMessages(
      gs_add_gating_method(
        gs = x,
        alias = alias,
        parent = parent,
        pop = "+",
        dims = NA, # EMPTY DIM WARNING
        gating_method = "cyto_gate_clean",
        gating_args = list(
          type = type,
          input = input,
          params = channels,
          slot = slot,
          openCyto.minEvents = -1,
          ...
        ),
        groupBy = NA,
        collapseDataForGating = FALSE,
        preprocessing_method = NA,
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
