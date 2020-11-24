## CYTO_PLOT_HEATMAP -----------------------------------------------------------

#' Construct a heatmap of cytometry data
#'
#' The \code{HeatmapR} package is required for heatmap construction, it can be
#' downloaded by running
#' \code{devtools::install_github("DillonHammill/HeatmapR")} in the console.
#'
#' @param x obejct of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param alias names of the populations for which the statistics should be
#'   calculated, set to NULL by default.
#' @param parent names of the parent populations for which frequency statistics
#'   should be calculated.
#' @param channels names of the channels or markers for which statistics should
#'   be calculated, set all channels with assigned markers by default.
#' @param stat name of the statistic passed to \code{\link{cyto_stats_compute}}
#'   to perform the calculation, set to \code{median} by default to compute
#'   median fluorescent intensity in the selected channels.
#' @param merge_by vector of experiment variable names to merge the samples into
#'   groups prior to constructing the heatmap, set to NULL by default.
#' @param ... additional arguments passed to \code{heat_map()} to construct the
#'   heatmap.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @name cyto_plot_heatmap
NULL

#' @noRd
#' @export
cyto_plot_heatmap <- function(x, 
                              ...){
  UseMethod("cyto_plot_heatmap")
}

#' @rdname cyto_plot_heatmap
#' @export
cyto_plot_heatmap.flowSet <- function(x, 
                                      channels = NULL,
                                      stat = "median",
                                      merge_by = NULL,
                                      ...){
  
  # CHANNELS
  if(is.null(channels)){
    channels <- names(cyto_markers(x))
  }
  
  # MERGE_BY
  if(!is.null(merge_by)){
    fr_list <- cyto_merge_by(x,
                             merge_by = merge_by)
  }else{
    fr_list <- cyto_convert(x, "list of flowFrames")
  }
  
  # COMPUTE STATISTICS
  cyto_stats <- lapply(fr_list, function(fr){
    cyto_stats_compute(fr,
                       stat = stat,
                       channels = channels,
                       format = "wide")
  })
  cyto_stats <- do.call("rbind", cyto_stats)
  
  # EXPERIMENT DETAILS
  if(is.null(merge_by)){
    cyto_details <- cyto_details(x)
  }else{
    cyto_details <- cyto_details(x)[, merge_by]
  }
  cyto_stats <- cbind(cyto_details, cyto_stats)
  
  # HEATMAP
  if(requireNamespace("HeatmapR")){
    cyto_heat <- HeatmapR::heat_map(cyto_stats, 
                                    ...)
    return(cyto_heat)
  }else{
    stop(paste0("The HeatmapR package is required to construct heatmaps:",
                " devtools::install_github('DillonHammill/HeatmapR')"))
  }
  
}

#' @rdname cyto_plot_heatmap
#' @export
cyto_plot_heatmap.flowFrame <- function(x, 
                                        channels = NULL,
                                        stat = "median",
                                        ...){
  
  # CHANNELS
  if(is.null(channels)){
    channels <- names(cyto_markers(x))
  }
  
  # COMPUTE STATISTICS
  cyto_stats <- cyto_stats_compute(x,
                                   channels = channels,
                                   stat = stat,
                                   format = "wide")
  
  # HEATMAP
  if(requireNamespace("HeatmapR")){
    cyto_heat <- HeatmapR::heat_map(cyto_stats, 
                                    ...)
    return(cyto_heat)
  }else{
    stop(paste0("The HeatmapR package is required to construct heatmaps:",
                " devtools::install_github('DillonHammill/HeatmapR')"))
  }
  
}