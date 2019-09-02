# CYTO_PLOT_GATE ---------------------------------------------------------------

#' Plot Gate Objects onto an Existing cyto_plot
#'
#' @param x gate object of class
#'   \code{\link[flowCore:rectangleGate-class]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate-class]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}}, \code{list} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param channels fluorescent channels to used to construct the plot.
#' @param gate_line_type integer [0,6] which controls the line type, set to
#'   \code{1} to draw solid lines by default.
#' @param gate_line_width numeric to adjust line thickness of gates, set to
#'   \code{2.5} by default.
#' @param gate_line_col indicates the colour of the gate to be constructed, set
#'   to \code{"red"} by default.
#' @param gate_fill fill colour(s) to use for gates, set to "white" by default.
#' @param gate_fill_alpha numeric [0,1] to control gate fill colour
#'   transparency, set to 0 by default to make gate fill colour completely
#'   transparent.
#'
#' @importFrom graphics par rect polygon abline
#' @importFrom grDevices adjustcolor
#' @importFrom flowCore parameters
#' @importFrom tools file_ext
#'
#' @return invisibly return modified gate objects with dimensions appropriate
#'   for the constructed plot.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_plot_gate3
NULL

#' @noRd
#' @export
cyto_plot_gate3 <- function(x, ...){
  UseMethod("cyto_plot_gate3")
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.default <- function(x, ...){
  stop(paste("cyto_plot does not support objects of class", class(x),"."))
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.rectangleGate <- function(gate,
                                          channels,
                                          gate_line_type = 1,
                                          gate_line_width = 2.5,
                                          gate_line_col = "red",
                                          gate_fill = "white",
                                          gate_fill_alpha = 0){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # LIMITS
  lims <- par("usr")
  ymin <- lims[3]
  ymax <- lims[4]
  ypad <- abs(ymin)/2
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
  if(!missing(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    chanenls <- parameters(gate)
  }

  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    rect(xleft = gate@min[channels[1]],
         xright = gate@max[channels[1]],
         ybottom = ymin + ypad,
         ytop = ymax - ypad,
         lty = gate_line_type,
         lwd = gate_line_width,
         border = gate_line_col,
         col = adjustcolor(gate_fill, gate_fill_alpha))
  # 2D PLOT  
  }else if(length(channels) == 2){
    rect(xleft = gate@min[channels[1]],
         xright = gate@max[channels[1]],
         ybottom = gate@min[channels[2]],
         ytop = gate@max[channels[2]],
         lty = gate_line_type,
         lwd = gate_line_width,
         border = gate_line_col,
         col = adjustcolor(gate_fill, gate_fill_alpha))
  }
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.polygonGate <- function(gate,
                                        channels,
                                        gate_line_type = 1,
                                        gate_line_width = 2.5,
                                        gate_line_col = "red",
                                        gate_fill = "white",
                                        gate_fill_alpha = 0){
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!missing(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    channels <- parameters(gate)
  }

  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT - rectangle
  if( length(channels) == 1){
    rect(xleft = gate@min[channels[1]],
         xright = gate@max[channels[1]],
         ybottom = ymin + ypad,
         ytop = ymax - ypad,
         lty = gate_line_type,
         lwd = gate_line_width,
         border = gate_line_col,
         col = adjustcolor(gate_fill, gate_fill_alpha))
  # 2D PLOT - polygon
  }else if(length(channels) == 2){
    polygon(gate@boundaries[, channels[1]],
            gate@boundaries[, channels[2]],
            lty = gate_line_type,
            lwd = gate_line_width,
            border = gate_line_col,
            col = adjustcolor(gate_fill, gate_fill_alpha)
    )
  }

  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.ellipsoidGate <- function(gate,
                                          channels,
                                          gate_line_type = 1,
                                          gate_line_width = 2.5,
                                          gate_line_col = "red",
                                          gate_fill = "white",
                                          gate_fill_alpha = 0){
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CONVERT GATE ---------------------------------------------------------------
  
  # GATE DIMENSIONS
  if(!missing(channels)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }else{
    channels <- parameters(gate)
  }
  
  # POLYGONGATE
  gate <- as(gate, "polygonGate")
  
  # PLOT GATE ------------------------------------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    rect(xleft = gate@min[channels[1]],
         xright = gate@max[channels[1]],
         ybottom = ymin + ypad,
         ytop = ymax - ypad,
         lty = gate_line_type,
         lwd = gate_line_width,
         border = gate_line_col,
         col = adjustcolor(gate_fill, gate_fill_alpha))
  # 2D PLOT
  }else if(length(channels) == 2){
    polygon(gate@boundaries[, channels[1]],
            gate@boundaries[, channels[2]],
            lty = gate_line_type,
            lwd = gate_line_width,
            border = gate_line_col,
            col = adjustcolor(gate_fill, gate_fill_alpha)
    )
  }

  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.quadGate <- function(gate,
                                     channels,
                                     gate_line_type = 1,
                                     gate_line_width = 2.5,
                                     gate_line_col = "red",
                                     gate_fill = "white",
                                     gate_fill_alpha = 0){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  lims <- par("usr")
  xmin <- lims[1]
  xmax <- lims[2]
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
  if(!missing(channels)){
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
      rect(
        xleft = xleft,
        xright = xright,
        ybottom = ybottom,
        ytop = ytop,
        col = adjustcolor(gate_fill[z], gate_fill_alpha[z])
      )
    })
  }
    
  # 2D PLOT GATE
  abline(
    v = gate@boundary[channels[1]],
    h = gate@boundary[channels[2]],
    lty = gate_line_type,
    lwd = gate_line_width, 
    col = gate_line_col
  )
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.filters <- function(gate,
                                    channels,
                                    gate_line_type = 1,
                                    gate_line_width = 2.5,
                                    gate_line_col = "red",
                                    gate_fill = "white",
                                    gate_fill_alpha = 0){
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # LIST OF GATES
  gate <- unlist(gate)
  
  # LIST GATE METHOD -----------------------------------------------------------
  
  # PLOT GATES
  gate <- cyto_plot_gate3(gate = gate,
                          channels = channels,
                          gate_line_type = gate_line_type,
                          gate_line_width = gate_line_width,
                          gate_line_col = gate_line_col,
                          gate_fill = gate_fill,
                          gate_fill_alpha = gate_fill_alpha)
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE(S) WITH CORRECT DIMENSIONS
  invisible(gate)
  
}

#' @rdname cyto_plot_gate3
#' @export
cyto_plot_gate3.list <- function(gate,
                                 channels,
                                 gate_line_type = 1,
                                 gate_line_width = 2.5,
                                 gate_line_col = "red",
                                 gate_fill = "white",
                                 gate_fill_alpha = 0){
   
  # PREPARE GATE ---------------------------------------------------------------
  
  # LIST OF GATE OBJECTS - WATCH OUT FOR FILTERS
  gate <- unlist(gate)
  
  # UNIQUE GATE ----------------------------------------------------------------
  
  # REPLICATE GATES - FIRST GATE ONLY
  if(length(unique(gate)) != 1){
    gate <- gate[1]
  }
  
  # GATE & POPULATION COUNTS ---------------------------------------------------
  
  # GATE COUNT
  gate_count <- length(gate)
  
  # POPULATION COUNT - GATE_FILL ARGUMENTS
  pop_count <- c()
  pop_count <- LAPPLY(gate, function(z){
    if(class(z) == "quadGate"){
      P <<- c(P, 4)
    }else{
      P <<- c(P, 1)
    }
  })
  
  # PREPARE GATE_FILL ARGUMENTS ------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # REPEAT ARGUMENTS
  lapply(c("gate_fill","gate_fill_alpha"), function(z){
    res <- rep(args[[z]], pop_count)
    args[[z]] <- spilt(res, rep(seq_len(gate_count),
                                times = pop_count))
  })
    
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CALL METHOD ----------------------------------------------------------------
  
  # LOOP THROUGH GATE
  gate <- mapply(function(gate,
                          gate_line_type,
                          gate_line_width,
                          gate_line_col,
                          gate_fill,
                          gate_fill_alpha){
    
    cyto_plot_gate3(gate,
                    channels = channels,
                    gate_line_type = gate_line_type,
                    gate_line_width = gate_line_width,
                    gate_line_col = gate_line_col,
                    gate_fill = gate_fill,
                    gate_fill_alpha = gate_fill_alpha)
    
  }, gate,
  gate_line_type,
  gate_line_width,
  gate_line_col,
  gate_fill,
  gate_fill_alpha, SIMPLIFY = FALSE)
  
  # RETURN GATE ----------------------------------------------------------------
  
  # GATE(S) WITH CORRECT DIMENSIONS
  invisible(gate)
  
}
