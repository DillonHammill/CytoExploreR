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
#'   gating_args = "display=0.5,alias='Lymphocytes',type='ellipse'",
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
                              alias, ...) {
  
  # Determine vertices of polygon using gate_draw
  gates <- cyto_gate_draw(x = fr, 
                          channels = channels, 
                          alias = alias, ...)
  
  return(gates)
}

#' Apply Manually Drawn Gates to GatingSet
#'
#' This \code{openCyto} plugin for \code{cyto_gate_draw} is required to extract
#' gates stored in the gatingTemplate and appropraitely apply them to the
#' GatingSet.
#'
#' @param fr a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#'   containing the flow cytometry data for gating.
#' @param pp_res output of preprocessing function.
#' @param channels fluorescent channel(s) to use for gating.
#' @param gate stored \code{cyto_gate_draw} gate in csv file for passing to
#'   openCyto.
#'
#' @return pass saved gate to openCyto to apply to all samples.
#'
#' @importFrom openCyto registerPlugins
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_gate_draw <- function(fr,
                            pp_res,
                            channels,
                            gate) {
  
  # pp_res is NULL - no grouping
  if (is.null(pp_res)) {
    gt <- 1
  }
  
  # Index of gate to use
  if (is.numeric(pp_res)) {
    gt <- pp_res
  } else if (is.character(pp_res)) {
    
    # names of gates must contain merged groupBy info
    gt <- match(pp_res, names(gate))
  }
  
  return(gate[[gt]])
}

#' cyto_gate_draw Preprocessing Function
#'
#' This preprocessing function passes on the index of the gate to apply to the
#' samples.
#'
#' @param fs a \code{\link[flowCore:flowSet-class]{flowSet}} object containing
#'   samples in group for gating.
#' @param gs a \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#'   all the samples to be gated.
#' @param gm gating method to use for gating, not used.
#' @param channels, names of the channels used to construct the gate, not used.
#' @param groupBy grouping variable to partition data, can either a single
#'   numeric or vector of pData column names.
#' @param isCollapse logical indicating when the data should be collapsed prior
#'   to gating, not used.
#' @param ... not in use.
#'
#' @return index of gate to apply to samples.
#'
#' @importFrom openCyto registerPlugins
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.pp_cyto_gate_draw <- function(fs,
                               gs,
                               gm,
                               channels = NA,
                               groupBy = NA,
                               isCollapse = NA, ...) {
  
  # Samples
  smp <- length(gs)
  
  # Extract pData information
  pd <- cyto_details(gs)
  
  # groupBy
  if (is.character(groupBy) & nchar(groupBy) == 0) {
    groupBy <- NA
  }
  
  # split groupBy
  if (is.character(groupBy)) {
    groupBy <- unlist(strsplit(groupBy, ":"))
  }
  
  # Add groupBy info to pData gs
  if (!all(is.na(groupBy)) & all(grepl("^[A-Za-z]+$", groupBy))) {
    
    # groupBy is composed of characters
    pd$groupby <- do.call(paste, pd[, groupBy, drop = FALSE])
    
    grpby <- pd[, "groupby"][match(cyto_details(fs[1])[, "name"], pd[, "name"])]
  } else if (!all(is.na(groupBy)) & !all(grepl("^[A-Za-z]+$", groupBy))) {
    if (groupBy == smp) {
      grpby <- 1
    } else {
      message("Numeric groupBy is not supported. Grouping all samples.")
      grpby <- 1
    }
  }
  
  # No grouping select first gate - only 1 gate expected
  if (all(is.na(groupBy))) {
    grpby <- 1
  }
  
  return(grpby)
}
