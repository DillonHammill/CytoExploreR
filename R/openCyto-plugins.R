#' drawGate Plugin for openCyto
#'
#' \code{drawGate} allows the user to draw polygon gates directly onto plots of flow cytometry data.
#' Simply left click to gate and right click to close the gate.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for gating.
#' @param pp_res output of preprocessing function.
#' @param channels a vector of length 2 indicating the channels used to construct the 2D plot.
#' @param filterId gate name assigned by openCyto from the \code{gatingTemplate}.
#' @param gate_range range in which gate should be constructed (only needed for autogating functions).
#' @param gate_type type of gate to be constructed, supported types include 
#' \code{c("polygon", "rectangle", "interval", "threshold", "ellipse", "quadrant")}.
#' @param min argument passed to \code{truncate_flowFrame} to restrict data to values > \code{min}.
#' @param max argument passed to \code{truncate_flowFrame} to restrict data to values < \code{max}.
#' @param N an integer indicating the number of gates to construct, set to 1 by default.
#' @param axis indicates the axis to use for gating for \code{gate_type="interval"} when 2 fluorescent channel are supplied.
#' @param adjust numeric smoothing factor used for 1D density plots.
#' @param ... additional arguments passsed to \code{DrawGate}.
#'
#' @return a \code{polygonGate} constructed from coordinates supplied by \code{DrawGate}.
#'
#' @keywords manual, gating, polygon, polygonGate
#' @import flowDensity openCyto
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @examples{
#' library(openCyto)
#' library(cytoTools)
#' 
#' fs <- read.flowSet(path = "Samples", pattern = ".fcs") # load in .fcs files
#' 
#' gs <- GatingSet(fs) # add flowSet to GatingSet
#' 
#' template <- add_pop(
#' gs, alias = "Lymphocytes", pop = "+", parent = "root", dims = "FSC-A,SSC-A", gating_method = "drawGate", gating_args = "gate_type='ellipse'",
#' gating_args = "subSample=10000", collapseDataForGating = TRUE, groupBy = 2
#' )
#' 
#' # gating window will open to construct gate left click vertices on plot and close gate by right click and selecting "stop".
#' 
#' ggcyto(gs[[1]], subset = "root", aes(x = "FSC-A",y = "SSC-A")) + geom_hex(bins = 100) + geom_stats()
#' 
#' }
gate_draw <- function(fr, pp_res, channels, alias, subSample, gate_range = NULL, min = NULL, max = NULL, gate_type = c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q"), axis = c("x","y"), adjust = 1.5, plot = TRUE, labs = TRUE,...){
  
  gate_type <- match.arg(gate_type)
  
  if(missing(axis)){
    axis <- "x"
  }
  
  axis <- match.arg(axis)
  
  # Two fluorescent channels must be supplied
  if(missing(channels) | !length(channels) %in% c(1,2)){
    stop("Two fluorescent channels must be specified to draw polygon gate.")
  }
  
  # Truncate flowFrame if min and max arguments are supplied
  if(!(is.null(min) && is.null(max))){
    fr <- openCyto::.truncate_flowframe(fr, channels = channels, min = min, max = max)
  }
  
  # Determine vertices of polygon using DrawGate
  gates <- drawGate(x, channels, alias, subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...)
  
  return(gates)
}
registerPlugins(gate_draw, "drawGate")

#' Apply Stored Gates to GatingSet
#'
#' \code{drawGate} allows the user to draw gates directly onto plots of flow cytometry data.This plugin provides
#' a way of passing saved gates to \code{openCyto}.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for gating.
#' @param pp_res output of preprocessing function.
#' @param channels fluorescent channel(s) to use for gating.
#' @param gate stored \code{drawGate} gate in csv file for passing to openCyto.
#'
#' @return pass saved gate to openCyto to apply to all samples.
#'
#' @keywords manual, gating, polygon, polygonGate
#' @import flowDensity openCyto
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
gate_manual <- function(fr, pp_res, channels, gate){
  
  return(gate)
}
registerPlugins(gate_manual, "manualGate")