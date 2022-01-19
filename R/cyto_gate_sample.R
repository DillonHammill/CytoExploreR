## CYTO_GATE_SAMPLE ------------------------------------------------------------

#' Gate a sample of a node in a GatingHierarchy or GatingSet
#'
#' \code{cyto_gate_sample()} provides a way to downsample events within nodes in
#' a GatingHierarchy or GatingSet and store these events as a new node.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parental population to sample, set to the
#'   \code{"root"} node by default.
#' @param alias name of the sampled population to be created in the GatingSet.
#' @param merge_by vector of \code{\link{cyto_details}} column names (e.g.
#'   c("Treatment","Concentration")) indicating how the samples should be grouped
#'   prior to sampling, set to "name" by default to downsample each sample
#'   separately.
#' @param events numeric to control the number of events to downsample each
#'   group in \code{merge_by}, set to 1 by default to include all events.
#'   \code{events} can be supplied per group otherwise the same degree of
#'   sampling will be applied to all groups. \code{cyto_gate_sample()} calls
#'   \code{cyto_sample_n()} on each group to determine the proportionate number
#'   of events to gate in each sample.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to which the
#'   \code{gatingTemplate} entries for the \code{GatingSet} method should be
#'   saved, set to \code{cyto_gatingTemplate_active()} by default.
#' @param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'   sampling is applied with each call to \code{cyto_gate_sample()}, set to an
#'   arbitrary numeric by default. This behaviour can be turned off by setting
#'   this argument to NULL.
#' @param group_by same as \code{merge_by} included for backwards compatibility
#'   with older versions of CytoExploreR, \code{group_by} was renamed to
#'   \code{merge_by} in CytoExploreR v2.0.0 to maintain consistency with
#'   \code{cyto_plot}. Users should therefore use \code{merge_by} as support for
#'   \code{group_by} will be ended in the future.
#' @param ... not in use.
#'
#' @importFrom openCyto gs_add_gating_method
#'
#' @return a \code{GatingHierarchy} or \code{GatingSet} with new sampled node
#'   added and an updated \code{gatingTemplate}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_sample}}
#' @seealso \code{\link{cyto_gatingTemplate_active}}
#' @seealso \code{\link{cyto_gatingTemplate_edit}}
#' @seealso \code{\link{cyto_gate_remove}}
#' @seealso \code{\link{cyto_gate_rename}}
#'
#' @examples 
#' \dontrun{
#' library(CytoExploreRData)
#' 
#' # Prepare Activation GatingSet
#' gs <- GatingSet(Activation)
#' gs <- cyto_compensate(gs)
#' gs <- cyto_transform(gs)
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#' 
#' # Write gatingTemplate to file
#' cyto_gatingTemplate_write(Activation_gatingTemplate, "gatingTemplate.csv")
#' 
#' # Add sampled T Cells node - 500 events per Treatment group
#' cyto_gate_sample(
#'   gs,
#'   parent = "T Cells",
#'   alias = "T Cells Sample",
#'   merge_by = "Treatment",
#'   events = 500,
#'   gatingTemplate = "gatingTemplate.csv"
#' )
#' }
#'
#' @export
cyto_gate_sample <- function(x,
                             parent = "root",
                             alias = NULL,
                             merge_by = "name",
                             events = 1,
                             gatingTemplate = NULL,
                             seed = 56,
                             group_by = NULL,
                             ...) {
  
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
  
  # ACTIVE GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active(ask = TRUE)
  }
  # CHECK EXISTING ENTRIES IN GATINGTEMPLATE
  gt <- .cyto_gatingTemplate_check(parent, 
                                   alias, 
                                   gatingTemplate)
  # CREATE GATINGTEMPLATE
  if (is.null(gt)) {
    message(
      paste("Creating", gatingTemplate, "to save the constructed gate(s).")
    )
    cyto_gatingTemplate_create(gatingTemplate, active = TRUE)
    gt <- cyto_gatingTemplate_read(gatingTemplate, data.table = TRUE)
  }
  
  # GATINGTEMPLATE ENTRIES -----------------------------------------------------
  
  # GROUP_BY -> MERGE_BY
  if(length(group_by) > 0) {
    message(
      "Support for 'group_by' is ending. Please use 'merge_by' instead."
    )
    merge_by <- group_by
  }
  
  # PREPARE GROUP_BY
  if(all(is.character(merge_by))) {
    if(all(merge_by == "all")) {
      group_by <- "NA"
    } else {
      group_by <- paste(merge_by, collapse = ":")
    }
  } else {
    group_by <- "NA"
  }
  
  # gs_add_gating_method - GATINGTEMPLATE ENTRY & APPLY TO GATINGSET
  message(paste("Adding new sample filters to", cyto_class(x), "..."))
  
  # ADD POPULATION TO GATINGSET
  pop <- suppressWarnings(
    suppressMessages(
      gs_add_gating_method(
        gs = x,
        alias = alias,
        parent = parent,
        pop = "+",
        dims = NA, # EMPTY DIM WARNING
        gating_method = "cyto_gate_sample",
        gating_args = list(
          openCyto.minEvents = -1
        ),
        groupBy = group_by,
        collapseDataForGating = FALSE, # REQUIRED
        preprocessing_method = "pp_cyto_gate_sample",
        preprocessing_args =  list(
          events = events,
          seed = seed
        )
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
