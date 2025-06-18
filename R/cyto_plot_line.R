## CYTO_PLOT_LINE --------------------------------------------------------------

#' Add lines with points to empty cyto_plot
#'
#' @param x a list of matrices or data.frames containing the \code{"x"} and
#'   \code{"y"} co-ordinates of the lines.
#' @param line_type integer [1,4] to control the type of each line, set to 1 by
#'   default for solid lines.
#' @param line_width numeric to control the width of each line, set to 1 by
#'   default.
#' @param line_col colour(s) to use for each line, set to "red" by default.
#' @param line_col_alpha numeric [0, 1] to control the transparency of line
#'   colours, set to 1 by default for solid colours.
#' @param line_fill colour(s) to use to fill the polygons below each line, set
#'   to NA by default to remove the polygons.
#' @param line_fill_alpha numeric [0, 1] to control the transparency of the
#'   polygon fill colour(s), set to 1 by default for solid colour(s).
#' @param line_point_shape shape(s) to use for points, set to \code{NA} by
#'   default to remove points from the plot. See \code{\link[graphics:par]{pch}}
#'   for alternatives.
#' @param line_point_size numeric to control the size of points per line, set to
#'   1 by default.
#' @param line_point_col colour(s) to use for the points in each line, set to
#'   "red" by default.
#' @param line_point_col_alpha numeric [0, 1] to control the transparency of
#'   point colours, set to 1 by default for solid colours.
#' @param ... not in use.
#'
#' @return NULL
#'
#' @importFrom grDevices adjustcolor
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_plot_line <- function(x,
                           line_type = 1,
                           line_width = 1,
                           line_col = "red",
                           line_col_alpha = 1,
                           line_fill = NA,
                           line_fill_alpha = 1,
                           line_point_shape = NA,
                           line_point_size = 1,
                           line_point_col = NA,
                           line_point_col_alpha = NA,
                           ...) {
  
  # TODO: INHERIT CYTO_PLOT() DEFAULT COLOURS
  
  # ADD LINES TO PLOT - MAPPLY REPEAT ARGUMENTS
  mapply(
    function(
      x,
      line_type,
      line_width,
      line_col,
      line_col_alpha,
      line_point_shape,
      line_point_size,
      line_point_col,
      line_point_col_alpha
    ) {
      
      # X - MATRIX OR DATA.FRAME OF XY COORDS
      if(!cyto_class(x, c("matrix", "data.frame"))){
        stop(
          "'x' must be a list of XY coordinates in a matrix or data.frame!"
        )
      }
      
      # XY COLUMNS REQUIRED
      if(!all(c("x", "y") %in% colnames(x))) {
        stop(
          "'x' must contain columns 'x' and 'y'!"
        )
      }
        
      # ONLY ADD IF COORDS EXIST
      if(nrow(x) > 0) {
        # ADD LINE FILL
        if(!is.na(line_fill)) {
          polygon(
            x = c(min(x[, "x"], x[, "x"], max(x[, "x"]))),
            y = c(min(c(0, x[, "y"])), x[, "y"], min(c(0, x[, "y"]))),
            border = NA,
            col = adjustcolor(
              line_fill,
              line_fill_alpha
            )
          )
        }
        # ADD LINE TO PLOT
        lines(
          x = x[, "x"],
          y = x[, "y"],
          type = "l",
          lty = line_type,
          lwd = line_width,
          col = adjustcolor(line_col, line_col_alpha)
        )
        # ADD POINTS TO PLOT
        if(!is.na(line_point_shape)) {
          points(
            x = x[, "x"],
            y = x[, "y"],
            pch = line_point_shape,
            cex = line_point_size,
            col = adjustcolor(
              line_point_col,
              line_point_col_alpha
            )
          )
        }
      }
    },
    x,
    line_type,
    line_width,
    line_col,
    line_col_alpha,
    line_point_shape,
    line_point_size,
    line_point_col,
    line_point_col_alpha
  )
  
  # INVISIBLE NULL RETURN
  invisible(NULL)
  
}