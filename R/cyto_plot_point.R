# CYTO_PLOT_POINT --------------------------------------------------------------

#' Add points and contour lines to empty cyto_plot
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
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 0 by default to exclude contour lines.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 2 by default.
#' @param contour_line_col colour to use for contour lines, set to "black" by
#'   default.
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
                                      overlay = NA,
                                      point_shape = ".",
                                      point_size = 2,
                                      point_col_scale = NA,
                                      point_cols = NA,
                                      point_col = NA,
                                      point_col_alpha = 1,
                                      contour_lines = 0,
                                      contour_line_type = 1,
                                      contour_line_width = 1,
                                      contour_line_col = "black", ...){
 
  # PREPARE DATA ---------------------------------------------------------------
  
  # Combine x and overlay into a list
  if(!.all_na(overlay)){
    fr_list <- c(list(x), .cyto_convert(overlay, "flowFrame list"))
  }else{
    fr_list <- list(x)
  }
   
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE
  .args_update(args)
  
  # POINT_COL
  point_col <- .cyto_plot_point_col(fr_list,
                                    channels = channels,
                                    point_col_scale = point_col_scale,
                                    point_cols = point_cols,
                                    point_col = point_col,
                                    point_col_alpha = point_col_alpha)
  
  # SPLIT ARGUMENTS ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args,
                                channels = channels,
                                n = length(fr_list),
                                plots = 1,
                                layers = length(fr_list),
                                gates = 0)
  
  # UNLIST ARGUMENTS - SINGLE PLOT
  args <- lapply(args, function(z){z[[1]]})
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # POINT & CONTOUR LINES LAYERS -----------------------------------------------
  
  # LAYERS
  lapply(seq_len(length(fr_list)), function(z){
    
    # EXTRACT DATA
    fr_exprs <- exprs(fr_list[[z]])[, channels]
    
    # POINTS - SKIP NO EVENTS
    if(!is.null(nrow(fr_exprs))){
      
      # POINTS
      if(nrow(fr_exprs) != 0){
        points(x = fr_exprs[,channels[1]],
               y = fr_exprs[,channels[2]],
               pch = point_shape[z],
               cex = point_size[z],
               col = point_col[z])
      }
      
    }
    
    # CONTOUR_LINES
    if(contour_lines[z] != 0){
      cyto_plot_contour(fr_list[[z]],
                        channels = channels,
                        contour_lines = contour_lines[z],
                        contour_line_type = contour_line_type[z],
                        contour_line_width = contour_line_width[z],
                        contour_line_col = contour_line_col[z])
      
    }
    
  })
  
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
                                 point_col_alpha = 1,
                                 contour_lines = 0,
                                 contour_line_type = 1,
                                 contour_line_width = 1,
                                 contour_line_col = "black", ...){
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Update arguments
  .args_update(args)
  
  # POINT_COL
  point_col <- .cyto_plot_point_col(x,
                                    channels = channels,
                                    point_col_scale = point_col_scale,
                                    point_cols = point_cols,
                                    point_col = point_col,
                                    point_col_alpha = point_col_alpha)

  # SPLIT ARGUMENTS ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args,
                                channels = channels,
                                n = length(x),
                                plots = 1,
                                layers = length(x),
                                gates = 0)
  
  # UNLIST ARGUMENTS - SINGLE PLOT
  args <- lapply(args, function(z){z[[1]]})
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # POINTS & CONTOUR LINES LAYERS ----------------------------------------------
  
  # LAYERS
  lapply(seq_len(length(x)), function(z){
    
    # EXTRACT DATA
    fr_exprs <- exprs(x[[z]])[, channels]
    
    # POINTS - SKIP NO EVENTS
    if(!is.null(nrow(fr_exprs))){
      
      # POINTS
      if(nrow(fr_exprs) != 0){
        points(x = fr_exprs[,channels[1]],
               y = fr_exprs[,channels[2]],
               pch = point_shape[z],
               cex = point_size[z],
               col = point_col[z])
      }
      
    }
    
    # CONTOUR_LINES
    if(contour_lines[z] != 0){
      cyto_plot_contour(x[[z]],
                        channels = channels,
                        contour_lines = contour_lines[z],
                        contour_line_type = contour_line_type[z],
                        contour_line_width = contour_line_width[z],
                        contour_line_col = contour_line_col[z])
      
    }
    
  })
  
}