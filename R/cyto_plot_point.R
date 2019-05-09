# CYTO_PLOT_POINT --------------------------------------------------------------

#' Add points to empty plots
#'
#' @param x object of class \code{flowFrame} or a list of \code{flowFrame}
#'   onjects.
#' @param channels names of the channels used to construct the plot.
#' @param overlay optional argument if x is a flowFrame to overlay a list of
#'   flowFrame objects.
#' @param point_shape shape(s) to use for points in 2-D scatterplots, set to
#'   \code{"."} by default to maximise plotting speed.  See
#'   \code{\link[graphics:par]{pch}} for alternatives.
#' @param point_size numeric to control the size of points in 2-D scatter plots
#'   set to 2 by default.
#' @param point_col_scale vector of ordered colours to use for the density
#'   colour gradient of points.
#' @param point_cols vector colours to draw from when selecting colours for
#'   points if none are supplied to point_col.
#' @param point_col colour(s) to use for points in 2-D scatter plots, set to NA
#'   by default to use a blue-red density colour scale.
#' @param point_col_alpha numeric [0,1] to control point colour transparency in
#'   2-D scatter plots, set to 1 by default to use solid colours.
#'
#' @importFrom flowCore exprs
#' @importFrom graphics points
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_plot_point
#'
#' @export
cyto_plot_point <- function(x, ...){
  UseMethod("cyto_plot_point")
}

#' @rdname cyto_plot_point
#' @export
cyto_plot_point.flowFrame <- function(x,
                                      channels,
                                      overlay,
                                      point_shape = ".",
                                      point_size = 2,
                                      point_col_scale = NA,
                                      point_cols = NA,
                                      point_col = NA,
                                      point_col_alpha = NA, ...){
 
  # Combine x and overlay into a list
  if(!.all_na(overlay)){
    fr_list <- c(list(x), .cyto_convert(overlay, "flowFrame list"))
  }else{
    fr_list <- list(x)
  }
   
  # Pull down arguments to named list
  args <- args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Get point colours
  args[["point_col"]] <- .cyto_plot_point_col(args[["fr_list"]],
                                              args[["channels"]],
                                              args[["point_col_scale"]],
                                              args[["point_cols"]],
                                              args[["point_col"]],
                                              args[["point_col_alpha"]])
  
  # Extract data and plot points
  mapply(function(x,
                  point_shape,
                  point_size,
                  point_col){
    
    # Extract data for plotting
    fr_exprs <- exprs(x)[,channels]
    
    # Plot points
    points(x = fr_exprs[,channels[1]],
           y = fr_exprs[,channels[2]],
           pch = point_shape,
           cex = point_size,
           col = point_col)
    
  }, args[["fr_list"]],
  args[["point_shape"]],
  args[["point_size"]],
  args[["point_col"]])
  
}

#' @rdname cyto_plot_point
#' @export
cyto_plot_point.list <- function(x,
                                 channels,
                                 point_shape = ".",
                                 point_size = 2,
                                 point_col_scale = NA,
                                 point_cols = NA,
                                 point_col = NA,
                                 point_col_alpha = NA, ...){
  
  # Must
  
  # Pull down arguments to named list
  args <- args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Get point colours
  args[["point_col"]] <- .cyto_plot_point_col(args[["x"]],
                                              args[["channels"]],
                                              args[["point_col_scale"]],
                                              args[["point_cols"]],
                                              args[["point_col"]],
                                              args[["point_col_alpha"]])
  
  # Extract data and plot points
  mapply(function(x,
                  point_shape,
                  point_size,
                  point_col){
    
    # Extract data for plotting
    fr_exprs <- exprs(x)[,channels]
    
    # Plot points
    points(x = fr_exprs[,channels[1]],
           y = fr_exprs[,channels[2]],
           pch = point_shape,
           cex = point_size,
           col = point_col)
    
  }, args[["x"]],
  args[["point_shape"]],
  args[["point_size"]],
  args[["point_col"]])
  
  
}