## CYTO_PLOT_GATE --------------------------------------------------------------

#' Plot Gate Objects onto an existing cyto_plot
#'
#' @param gate gate object of class
#'   \code{\link[flowCore:rectangleGate-class]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate-class]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}},
#'   \code{multiRangeGate} \code{list} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to
#'   \code{1} to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to
#'   \code{2.5} by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set
#'   to \code{"red"} by default.
#' @param gate_line_col_alpha numeric to control the transparency of gate lines,
#'   set to 1 by default.
#' @param gate_point_shape point shape to use for gate vertices, set to NA by
#'   default to plot lines only.
#' @param gate_point_size numeric to control the size of gate vertices, set to 1
#'   by default.
#' @param gate_point_col colour(s) to use for gate vertices, set to \code{"red"}
#'   by default.
#' @param gate_point_col_alpha numeric [0, 1] to control the transparency of
#'   gate vertices, set to 1 by default for solid colours.
#' @param gate_fill fill colour(s) to use for gates, set to "white" by default.
#' @param gate_fill_alpha numeric [0,1] to control gate fill colour
#'   transparency, set to 0 by default to make gate fill colour completely
#'   transparent.
#' @param ... not in use.
#'
#' @importFrom graphics par rect polygon abline
#' @importFrom grDevices adjustcolor
#' @importFrom flowCore parameters
#' @importFrom methods is
#'
#' @return invisibly return modified gate objects with dimensions appropriate
#'   for the constructed plot.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_plot_gate
NULL

#' @noRd
#' @export
cyto_plot_gate <- function(gate, ...){
  UseMethod("cyto_plot_gate")
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.default <- function(gate, ...){
  stop(paste("cyto_plot does not support objects of class", class(gate),"."))
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.rectangleGate <- function(gate,
                                         channels = NULL,
                                         gate_line_type = 1,
                                         gate_line_width = 2.5,
                                         gate_line_col = "red",
                                         gate_line_col_alpha = 1,
                                         gate_fill = "white",
                                         gate_fill_alpha = 0,
                                         gate_point_shape = NA,
                                         gate_point_size = 1,
                                         gate_point_col = "red",
                                         gate_point_col_alpha = 1,
                                         ...){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04)
  xmin <- xmin + 0.5 * xpad # 2% BUFFER
  xmax <- xmax - 0.5 * xpad # 2% BUFFER
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04)
  ymin <- ymin + 0.5 * ypad # 2% BUFFER
  ymax <- ymax - 0.5 * ypad # 2% BUFFER
  yrng <- ymax - ymin
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!is.null(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    channels <- parameters(gate)
  }
  
  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # REPLACE INFINITE X COORDS
    if(is.infinite(gate@min)){
      gate@min <- xmin
    }
    if(is.infinite(gate@max)){
      gate@max <- xmax
    }
    # PLOT GATE
    rect(
      xleft = gate@min,
      xright = gate@max,
      ybottom = ymin,
      ytop = ymax,
      lty = gate_line_type,
      lwd = gate_line_width,
      border = adjustcolor(gate_line_col, gate_line_col_alpha),
      col = adjustcolor(gate_fill, gate_fill_alpha)
    )
    # GATE VERTICES
    if(!is.na(gate_point_shape)) {
      points(
        x = c(
          gate@min,
          gate@min, 
          gate@max, 
          gate@max
        ),
        y = c(
          ymin, 
          ymax, 
          ymax,
          ymin
        ),
        pch = gate_point_shape,
        cex = gate_point_size,
        col = adjustcolor(
          gate_point_col,
          gate_point_col_alpha
        )
      )
    }
  # 2D PLOT  
  }else if(length(channels) == 2){
    # REPLACE INFINITE X COORDS
    if(is.infinite(gate@min[channels[1]])){
      # Quadrant
      if(any(grepl("quad", names(attributes(gate))))){
        gate@min[channels[1]] <- lims[1] - 0.025*(lims[2] - lims[1]) # off plot
      }else{
        gate@min[channels[1]] <- xmin
      }
    }
    if(is.infinite(gate@max[channels[1]])){
      # Quadrant
      if(any(grepl("quad", names(attributes(gate))))){
        gate@max[channels[1]] <- lims[2] + 0.025*(lims[2] - lims[1]) # off plot
      }else{
        gate@max[channels[1]] <- xmax
      }
    }
    # REPLACE INFINITE Y COORDS
    if(is.infinite(gate@min[channels[2]])){
      # Quadrant
      if(any(grepl("quad", names(attributes(gate))))){
        gate@min[channels[2]] <- lims[3] - 0.025*(lims[4] - lims[3]) # off plot
      }else{
        gate@min[channels[2]] <- ymin
      }
    }
    if(is.infinite(gate@max[channels[2]])){
      # Quadrant
      if(any(grepl("quad", names(attributes(gate))))){
        gate@max[channels[2]] <- lims[4] + 0.025*(lims[4] - lims[3]) # off plot
      }else{
        gate@max[channels[2]] <- ymax
      }
    }
    
    # PLOT GATE
    rect(
      xleft = gate@min[channels[1]],
      xright = gate@max[channels[1]],
      ybottom = gate@min[channels[2]],
      ytop = gate@max[channels[2]],
      lty = gate_line_type,
      lwd = gate_line_width,
      border = adjustcolor(gate_line_col, gate_line_col_alpha),
      col = adjustcolor(gate_fill, gate_fill_alpha)
    )
    # GATE VERTICES
    if(!is.na(gate_point_shape)) {
      points(
        x = c(
          gate@min[channels[1]],
          gate@min[channels[1]],
          gate@max[channels[1]], 
          gate@max[channels[1]]
        ),
        y = c(
          gate@min[channels[2]],
          gate@max[channels[2]],
          gate@max[channels[2]], 
          gate@min[channels[2]]
        ),
        pch = gate_point_shape,
        cex = gate_point_size,
        col = adjustcolor(
          gate_point_col,
          gate_point_col_alpha
        )
      )
    }
  }
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.polygonGate <- function(gate,
                                       channels = NULL,
                                       gate_line_type = 1,
                                       gate_line_width = 2.5,
                                       gate_line_col = "red",
                                       gate_line_col_alpha = 1,
                                       gate_fill = "white",
                                       gate_fill_alpha = 0,
                                       gate_point_shape = NA,
                                       gate_point_size = 1,
                                       gate_point_col = "red",
                                       gate_point_col_alpha = 1,
                                       ...){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04)
  xmin <- xmin + 0.5 * xpad # 2% BUFFER
  xmax <- xmax - 0.5 * xpad # 2% BUFFER
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04)
  ymin <- ymin + 0.5 * ypad # 2% BUFFER
  ymax <- ymax - 0.5 * ypad # 2% BUFFER
  yrng <- ymax - ymin
  
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!is.null(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    channels <- parameters(gate)
  }
  
  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT - RECTANGLE
  if( length(channels) == 1){
    # REPLACE INFINITE X COORDS
    if(is.infinite(gate@min[channels[1]])){
      gate@min[channels[1]] <- xmin
    }
    if(is.infinite(gate@max[channels[1]])){
      gate@max[channels[1]] <- xmax
    }
    # PLOT GATE
    rect(
      xleft = gate@min[channels[1]],
      xright = gate@max[channels[1]],
      ybottom = ymin,
      ytop = ymax,
      lty = gate_line_type,
      lwd = gate_line_width,
      border = adjustcolor(gate_line_col, gate_line_col_alpha),
      col = adjustcolor(gate_fill, gate_fill_alpha)
    )
    # GATE VERTICES
    if(!is.na(gate_point_shape)) {
      points(
        x = c(
          gate@min[channels[1]], 
          gate@min[channels[1]],
          gate@max[channels[1]], 
          gate@max[channels[1]]
        ),
        y = c(
          ymin, 
          ymax, 
          ymax,
          ymin
        ),
        pch = gate_point_shape,
        cex = gate_point_size,
        col = adjustcolor(
          gate_point_col,
          gate_point_col_alpha
        )
      )
    }
    # 2D PLOT - POLYGON OR RECTANGLE
  }else if(length(channels) == 2){
    # 2D RECTANGLE GATE
    if(is(gate, "rectangleGate")){
      # REPLACE INFINITE X COORDS
      if(is.infinite(gate@min[channels[1]])){
        gate@min[channels[1]] <- xmin
      }
      if(is.infinite(gate@max[channels[1]])){
        gate@max[channels[1]] <- xmax
      }
      # REPLACE INFINITE Y COORDS
      if(is.infinite(gate@min[channels[2]])){
        gate@min[channels[2]] <- ymin
      }
      if(is.infinite(gate@max[channels[2]])){
        gate@max[channels[2]] <- ymax
      }
      # PLOT GATE
      rect(
        xleft = gate@min[channels[1]],
        xright = gate@max[channels[1]],
        ybottom = gate@min[channels[2]],
        ytop = gate@max[channels[2]],
        lty = gate_line_type,
        lwd = gate_line_width,
        border = adjustcolor(gate_line_col, gate_line_col_alpha),
        col = adjustcolor(gate_fill, gate_fill_alpha)
      )
      # GATE VERTICES
      if(!is.na(gate_point_shape)) {
        points(
          x = c(
            gate@min[channels[1]], 
            gate@min[channels[1]],
            gate@max[channels[1]], 
            gate@max[channels[1]]
          ),
          y = c(
            gate@min[channels[2]], 
            gate@max[channels[2]], 
            gate@max[channels[2]],
            gate@min[channels[2]]
          ),
          pch = gate_point_shape,
          cex = gate_point_size,
          col = adjustcolor(
            gate_point_col,
            gate_point_col_alpha
          )
        )
      }
    # 2D POLYGONGATE  
    }else if(is(gate, "polygonGate")){
      # CO-ORDINATES
      x_coords <- gate@boundaries[, channels[1]]
      y_coords <- gate@boundaries[, channels[2]]
      # REPLACE INFINITE X COORDS
      if(any(is.infinite(x_coords))){
        ind <- which(is.infinite(x_coords))
        lapply(ind, function(z){
          if(x_coords[ind] < 0){
            x_coords[z] <<- xmin
          }else if(x_coords[ind] > 0){
            x_coords[z] <<- xmax
          }
        })
      }
      # REPLACE INFINITE Y COORDS
      if(any(is.infinite(y_coords))){
        ind <- which(is.infinite(y_coords))
        lapply(ind, function(z){
          if(y_coords[ind] < 0){
            y_coords[z] <<- ymin
          }else if(y_coords[ind] > 0){
            y_coords[z] <<- ymax
          }
        })
      }
      # PLOT GATE
      polygon(
        x_coords,
        y_coords,
        lty = gate_line_type,
        lwd = gate_line_width,
        border = adjustcolor(gate_line_col, gate_line_col_alpha),
        col = adjustcolor(gate_fill, gate_fill_alpha)
      )
      # GATE VERTICES
      if(!is.na(gate_point_shape)) {
        points(
          x = x_coords,
          y = y_coords,
          pch = gate_point_shape,
          cex = gate_point_size,
          col = adjustcolor(
            gate_point_col,
            gate_point_col_alpha
          )
        )
      }
    }
    
  }
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate
#' @export 
cyto_plot_gate.spectralGate <- function(gate,
                                        channels = NULL,
                                        gate_line_type = 1,
                                        gate_line_width = 2.5,
                                        gate_line_col = "red",
                                        gate_line_col_alpha = 1,
                                        gate_fill = "white",
                                        gate_fill_alpha = 0,
                                        gate_point_shape = NA,
                                        gate_point_size = 1,
                                        gate_point_col = "red",
                                        gate_point_col_alpha = 1,
                                        ...) {
  
  # NOTE: CHANNELS = ALL DETECTORS ON X AXIS
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04)
  xmin <- xmin + 0.5 * xpad # 2% BUFFER
  xmax <- xmax - 0.5 * xpad # 2% BUFFER
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04)
  ymin <- ymin + 0.5 * ypad # 2% BUFFER
  ymax <- ymax - 0.5 * ypad # 2% BUFFER
  yrng <- ymax - ymin
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # PLOT GATE ------------------------------------------------------------------
  
  # EXTRACT CO-ORDINATES
  coords <- gate@coords
  
  # GATE DIMENSIONS
  chans <- rownames(coords)
  
  # SUBSET DETECTORS
  coords <- coords[match(chans, rownames(coords)), ]
  
  # CO-ORDINATES
  x_coords <- coords[, 1]
  y_coords <- coords[, 2]
  
  # REPLACE INFINITE X COORDS
  if(any(is.infinite(x_coords))){
    ind <- which(is.infinite(x_coords))
    lapply(
      ind, 
      function(z) {
        if(x_coords[ind] < 0){
          x_coords[z] <<- xmin
        }else if(x_coords[ind] > 0){
          x_coords[z] <<- xmax
        }
      }
    )
  }
  # REPLACE INFINITE Y COORDS
  if(any(is.infinite(y_coords))){
    ind <- which(is.infinite(y_coords))
    lapply(
      ind, 
      function(z) {
        if(y_coords[ind] < 0) {
          y_coords[z] <<- ymin
        } else if(y_coords[ind] > 0) {
          y_coords[z] <<- ymax
        }
      }
    )
  }
  # PLOT FILL POLYGON
  polygon(
    c(1, x_coords, max(x_coords)),
    c(0, y_coords, 0),
    lty = gate_line_type,
    lwd = gate_line_width,
    border = NA,
    col = adjustcolor(
      gate_fill, 
      gate_fill_alpha
    )
  )
  # PLOT LINE
  lines(
    x_coords,
    y_coords,
    lty = gate_line_type,
    lwd = gate_line_width,
    col = adjustcolor(
      gate_line_col,
      gate_line_col_alpha
    )
  )
  # PLOT GATE VERTICES
  if(!is.na(gate_point_shape)) {
    points(
      x = x_coords,
      y = y_coords,
      pch = gate_point_shape,
      cex = gate_point_size,
      col = adjustcolor(
        gate_point_col,
        gate_point_col_alpha
      )
    )
  }
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
  
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.ellipsoidGate <- function(gate,
                                         channels = NULL,
                                         gate_line_type = 1,
                                         gate_line_width = 2.5,
                                         gate_line_col = "red",
                                         gate_line_col_alpha = 1,
                                         gate_fill = "white",
                                         gate_fill_alpha = 0,
                                         gate_point_shape = NA,
                                         gate_point_size = 1,
                                         gate_point_col = "red",
                                         gate_point_col_alpha = 1,
                                         ...){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04)
  xmin <- xmin + 0.5 * xpad # 2% BUFFER
  xmax <- xmax - 0.5 * xpad # 2% BUFFER
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04)
  ymin <- ymin + 0.5 * ypad # 2% BUFFER
  ymax <- ymax - 0.5 * ypad # 2% BUFFER
  yrng <- ymax - ymin
  
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!is.null(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    channels <- parameters(gate)
  }
  
  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # REPLACE INFINITE X COORDS
    if(is.infinite(gate@min[channels[1]])){
      gate@min[channels[1]] <- xmin
    }
    if(is.infinite(gate@max[channels[1]])){
      gate@max[channels[2]] <- xmax
    }
    # PLOT GATE
    rect(
      xleft = gate@min[channels[1]],
      xright = gate@max[channels[1]],
      ybottom = ymin,
      ytop = ymax,
      lty = gate_line_type,
      lwd = gate_line_width,
      border = adjustcolor(gate_line_col, gate_line_col_alpha),
      col = adjustcolor(
        gate_fill, 
        gate_fill_alpha
      )
    )
    # GATE VERTICES
    if(!is.na(gate_point_shape)) {
      points(
        x = c(
          gate@min[channels[1]], 
          gate@min[channels[1]],
          gate@max[channels[1]], 
          gate@max[channels[1]]
        ),
        y = c(
          ymin, 
          ymax, 
          ymax,
          ymin
        ),
        pch = gate_point_shape,
        cex = gate_point_size,
        col = adjustcolor(
          gate_point_col,
          gate_point_col_alpha
        )
      )
    }
    # 2D PLOT
  }else if(length(channels) == 2){
    # 2D RECTANGLEGATE
    if(is(gate, "rectangleGate")){
      # REPLACE INFINITE X COORDS
      if(is.infinite(gate@min[channels[1]])){
        gate@min[channels[1]] <- xmin
      }
      if(is.infinite(gate@max[channels[1]])){
        gate@max[channels[1]] <- xmax
      }
      # REPLACE INFINITE Y COORDS
      if(is.infinite(gate@min[channels[2]])){
        gate@min[channels[2]] <- ymin
      }
      if(is.infinite(gate@max[channels[2]])){
        gate@max[channels[2]] <- ymax
      }
      # PLOT GATE
      rect(
        xleft = gate@min[channels[1]],
        xright = gate@max[channels[1]],
        ybottom = gate@min[channels[2]],
        ytop = gate@max[channels[2]],
        lty = gate_line_type,
        lwd = gate_line_width,
        border = adjustcolor(gate_line_col, gate_line_col_alpha),
        col = adjustcolor(gate_fill, gate_fill_alpha)
      )
      # INTERPOLATED GATE VERTICES
      if(!is.na(gate_point_shape)) {
        points(
          x = c(
            gate@min[channels[1]], 
            gate@min[channels[1]],
            gate@max[channels[1]], 
            gate@max[channels[1]]
          ),
          y = c(
            gate@min[channels[2]], 
            gate@max[channels[2]], 
            gate@max[channels[2]],
            gate@min[channels[2]]
          ),
          pch = gate_point_shape,
          cex = gate_point_size,
          col = adjustcolor(
            gate_point_col,
            gate_point_col_alpha
          )
        )
      }
      # 2D ELLIPSOIDGATE  
    }else if(is(gate, "ellipsoidGate")){
      # COERCE TO POLYGONGATE
      gate <- as(gate, "polygonGate")
      # CO-ORDINATES
      x_coords <- gate@boundaries[,channels[1]]
      y_coords <- gate@boundaries[,channels[2]]
      # REPLACE INFINITE X COORDS
      if(any(is.infinite(x_coords))){
        ind <- which(is.infinite(x_coords))
        lapply(ind, function(z){
          if(x_coords[z] < 0){
            x_coords[z] <<- xmin
          }else if(x_coords[z] > 0){
            x_coords[z] <<- xmax
          }
        })
      }
      # REPLACE INFINITE Y COORDS
      if(any(is.infinite(y_coords))){
        ind <- which(is.infinite(y_coords))
        lapply(ind, function(z){
          if(y_coords[z] < 0){
            y_coords[z] <<- ymin
          }else if(y_coords[z] > 0){
            y_coords[z] <<- ymax
          }
        })
      }
      # PLOT GATE
      polygon(
        x_coords,
        y_coords,
        lty = gate_line_type,
        lwd = gate_line_width,
        border = adjustcolor(gate_line_col, gate_line_col_alpha),
        col = adjustcolor(gate_fill, gate_fill_alpha)
      )
      # GATE VERTICES
      if(!is.na(gate_point_shape)) {
        points(
          x = x_coords,
          y = y_coords,
          pch = gate_point_shape,
          cex = gate_point_size,
          col = adjustcolor(
            gate_point_col,
            gate_point_col_alpha
          )
        )
      }
    }
    
  }
  
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.quadGate <- function(gate,
                                    channels = NULL,
                                    gate_line_type = 1,
                                    gate_line_width = 2.5,
                                    gate_line_col = "red",
                                    gate_line_col_alpha = 1,
                                    gate_fill = "white",
                                    gate_fill_alpha = 0,
                                    gate_point_shape = NA,
                                    gate_point_size = 1,
                                    gate_point_col = "red",
                                    gate_point_col_alpha = 1,
                                    ...){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # GATE_FILL
  gate_fill <- rep_len(gate_fill, 4)
  gate_fill_alpha <- rep_len(gate_fill_alpha, 4)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!is.null(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    channels <- parameters(gate)
  }
  
  # PLOT GATE ------------------------------------------------------------------
  
  # GATE FILL - PER QUADRANT
  if(!all(gate_fill == "white")){
    lapply(seq_len(4), function(z){
      # QUADRANT 1 - TOP RIGHT
      if(z == 1){
        xleft <- gate@boundary[channels[1]]
        xright <- xmax
        ybottom <- gate@boundary[channels[2]]
        ytop <- ymax
        # QUADRANT 2 - TOP LEFT
      }else if(z == 2){
        xleft <- xmin
        xright <- gate@boundary[channels[1]]
        ybottom <- gate@boundary[channels[2]]
        ytop <- ymax
        # QUADRANT 3 - BOTTOM RIGHT
      }else if(z == 3){
        xleft <- gate@boundary[channels[1]]
        xright <- xmax
        ybottom <- ymin
        ytop <- gate@boundary[channels[2]]
        # QUADRANT 4 - BOTTOM LEFT
      }else if(z == 4){
        xleft <- xmin
        xright <- gate@boundary[channels[1]]
        ybottom <- ymin
        ytop <- gate@boundary[channels[2]]
      }
      # PLOT GATE FILL
      rect(
        xleft = xleft,
        xright = xright,
        ybottom = ybottom,
        ytop = ytop,
        col = adjustcolor(gate_fill[z], gate_fill_alpha[z]),
        border = NA
      )
    })
  }
  
  # 2D PLOT GATE
  abline(
    v = gate@boundary[channels[1]],
    h = gate@boundary[channels[2]],
    lty = gate_line_type,
    lwd = gate_line_width, 
    col = adjustcolor(gate_line_col, gate_line_col_alpha)
  )
  
  # GATE VERTICES
  if(!is.na(gate_point_shape)) {
    points(
      x = gate@boundary[channels[1]],
      y = gate@boundary[channels[2]],
      pch = gate_point_shape,
      cex = gate_point_size,
      col = adjustcolor(
        gate_point_col,
        gate_point_col_alpha
      )
    )
  }
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.multiRangeGate <- function(gate,
                                          channels = NULL,
                                          gate_line_type = 1,
                                          gate_line_width = 2.5,
                                          gate_line_col = "red",
                                          gate_line_col_alpha = 1,
                                          gate_fill = "white",
                                          gate_fill_alpha = 0,
                                          gate_point_shape = NA,
                                          gate_point_size = 1,
                                          gate_point_col = "red",
                                          gate_point_col_alpha = 1,
                                          ...) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04)
  xmin <- xmin + 0.5 * xpad # 2% BUFFER
  xmax <- xmax - 0.5 * xpad # 2% BUFFER
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04)
  ymin <- ymin + 0.5 * ypad # 2% BUFFER
  ymax <- ymax - 0.5 * ypad # 2% BUFFER
  yrng <- ymax - ymin
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!is.null(channels)){
    gate <- cyto_gate_convert(
      gate,
      channels = channels
    )
  }else{
    channels <- parameters(gate)
  }
  
  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    mapply(
      function(z,
               gate_line_type,
               gate_line_width,
               gate_line_col,
               gate_line_col_alpha,
               gate_fill,
               gate_fill_alpha,
               gate_point_shape,
               gate_point_size,
               gate_point_col,
               gate_point_col_alpha) {
        # GATE COORDINATES
        coords <- c(
          gate@ranges[["min"]][z],
          gate@ranges[["max"]][z]
        )
        # REPLACE INFINITE X COORDS
        if(is.infinite(min(coords))){
          coords[coords == min(coords)] <- xmin
        }
        if(is.infinite(max(coords))){
          coords[coords == max(coords)] <- xmax
        }
        # PLOT GATE(S)
        rect(
          xleft = coords[1],
          xright = coords[2],
          ybottom = ymin,
          ytop = ymax,
          lty = gate_line_type,
          lwd = gate_line_width,
          border = adjustcolor(gate_line_col, gate_line_col_alpha),
          col = adjustcolor(gate_fill, gate_fill_alpha)
        )
        # GATE VERTICES
        if(!is.na(gate_point_shape)) {
          points(
            x = c(
              coords[1], 
              coords[1],
              coords[2], 
              coords[2]
            ),
            y = c(
              ymin, 
              ymax, 
              ymax,
              ymin
            ),
            pch = gate_point_shape,
            cex = gate_point_size,
            col = adjustcolor(
              gate_point_col,
              gate_point_col_alpha
            )
          )
        }
      },
      seq_along(gate@ranges[[1]]),
      gate_line_type,
      gate_line_width,
      gate_line_col,
      gate_line_col_alpha,
      gate_fill,
      gate_fill_alpha,
      gate_point_shape,
      gate_point_size,
      gate_point_col,
      gate_point_col_alpha
    )
  # 2D PLOT  
  }else if(length(channels) == 2) {
    mapply(
      function(z,
               gate_line_type,
               gate_line_width,
               gate_line_col,
               gate_line_col_alpha,
               gate_fill,
               gate_fill_alpha,
               gate_point_shape,
               gate_point_size,
               gate_point_col,
               gate_point_col_alpha) {
        # GATE CHANNELS
        chan <- parameters(gate)
        # GATE CO-ORDINATES
        coords <- c(
          gate@ranges[["min"]][z],
          gate@ranges[["max"]][z]
        )
        # REPLACE INFINITE X COORDS
        if(is.infinite(min(coords))){
          coords[coords == min(coords)] <- if(chan %in% channels[1]) {
            xmin
          } else {
            ymin
          }
        }
        if(is.infinite(max(coords))){
          coords[coords == max(coords)] <- if(chan %in% channels[1]) {
            xmax
          } else {
            ymax
          }
        }
        # GATE LIMITS
        if(chan %in% channels[1]) {
          xmin <- min(coords)
          xmax <- max(coords)
        }
        if(!chan %in% channels[1]) {
          ymin <- min(coords)
          ymax <- max(coords)
        }
        # PLOT GATE(S)
        rect(
          xleft = xmin,
          xright = xmax,
          ybottom = ymin,
          ytop = ymax,
          lty = gate_line_type,
          lwd = gate_line_width,
          border = adjustcolor(gate_line_col, gate_line_col_alpha),
          col = adjustcolor(gate_fill, gate_fill_alpha)
        )
        # GATE VERTICES
        if(!is.na(gate_point_shape)) {
          points(
            x = c(
              xmin, 
              xmin,
              xmax, 
              xmax
            ),
            y = c(
              ymin, 
              ymax, 
              ymax,
              ymin
            ),
            pch = gate_point_shape,
            cex = gate_point_size,
            col = adjustcolor(
              gate_point_col,
              gate_point_col_alpha
            )
          )
        }
      },
      seq_along(gate@ranges[[1]]),
      gate_line_type,
      gate_line_width,
      gate_line_col,
      gate_line_col_alpha,
      gate_fill,
      gate_fill_alpha,
      gate_point_shape,
      gate_point_size,
      gate_point_col,
      gate_point_col_alpha
    )
  }
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.filters <- function(gate,
                                   channels = NULL,
                                   gate_line_type = 1,
                                   gate_line_width = 2.5,
                                   gate_line_col = "red",
                                   gate_line_col_alpha = 1,
                                   gate_fill = "white",
                                   gate_fill_alpha = 0,
                                   gate_point_shape = NA,
                                   gate_point_size = 1,
                                   gate_point_col = "red",
                                   gate_point_col_alpha = 1,
                                   ...){
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # LIST OF GATES
  gate <- unlist(gate)
  
  # LIST GATE METHOD -----------------------------------------------------------
  
  # PLOT GATES
  gate <- cyto_plot_gate(
    gate = gate,
    channels = channels,
    gate_line_type = gate_line_type,
    gate_line_width = gate_line_width,
    gate_line_col = gate_line_col,
    gate_line_col_alpha = gate_line_col_alpha,
    gate_fill = gate_fill,
    gate_fill_alpha = gate_fill_alpha,
    gate_point_shape = gate_point_shape,
    gate_point_size = gate_point_size,
    gate_point_col = gate_point_col,
    gate_point_col_alpha = gate_point_col_alpha
  )
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE(S) WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate
#' @export
cyto_plot_gate.list <- function(gate,
                                channels = NULL,
                                gate_line_type = 1,
                                gate_line_width = 2.5,
                                gate_line_col = "red",
                                gate_line_col_alpha = 1,
                                gate_fill = "white",
                                gate_fill_alpha = 0,
                                gate_point_shape = NA,
                                gate_point_size = 1,
                                gate_point_col = "red",
                                gate_point_col_alpha = 1,
                                ...){
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # LIST OF GATE OBJECTS - WATCH OUT FOR FILTERS & DUPLICATES
  gate <- unique(unlist(gate))
  
  # GATE & POPULATION COUNTS ---------------------------------------------------
  
  # GATE COUNT
  gate_count <- length(gate)
  
  # POPULATION COUNT - GATE_FILL ARGUMENTS
  pop_count <- c()
  lapply(
    gate, 
    function(z) {
      if(cyto_class(z, "quadGate")) {
        pop_count <<- c(pop_count, 4)
      }else{
        pop_count <<- c(pop_count, 1)
      }
    }
  )
  
  # PREPARE GATE_FILL ARGUMENTS ------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # REPEAT ARGUMENTS
  lapply(
    c("gate_fill","gate_fill_alpha"),
    function(z) {
      res <- rep(args[[z]], length.out = sum(pop_count))
      args[[z]] <<- split(
        res, 
        rep(
          seq_len(
            gate_count
          ),
          times = pop_count
        )
      )
    }
  )
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CALL METHOD ----------------------------------------------------------------
  
  # LOOP THROUGH GATE
  gate <- mapply(
    function(
      gate,
      gate_line_type,
      gate_line_width,
      gate_line_col,
      gate_line_col_alpha,
      gate_fill,
      gate_fill_alpha,
      gate_point_shape,
      gate_point_size,
      gate_point_col,
      gate_point_col_alpha
    ) {
    
    cyto_plot_gate(
      gate,
      channels = channels,
      gate_line_type = gate_line_type,
      gate_line_width = gate_line_width,
      gate_line_col = gate_line_col,
      gate_line_col_alpha = gate_line_col_alpha,
      gate_fill = gate_fill,
      gate_fill_alpha = gate_fill_alpha,
      gate_point_shape = gate_point_shape,
      gate_point_size = gate_point_size,
      gate_point_col = gate_point_col,
      gate_point_col_alpha = gate_point_col_alpha
    )
    
  }, 
  gate,
  gate_line_type,
  gate_line_width,
  gate_line_col,
  gate_line_col_alpha,
  gate_fill,
  gate_fill_alpha, 
  gate_point_shape,
  gate_point_size,
  gate_point_col,
  gate_point_col_alpha,
  SIMPLIFY = FALSE)
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE(S) WITH CORRECT DIMENSIONS
  invisible(gate)
  
}
