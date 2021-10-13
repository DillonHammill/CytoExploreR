## CYTO_GATE INTERNAL HELPERS --------------------------------------------------

# A collection of internal functions to compute gate centers and counts.

## .CYTO_GATE_CENTER -----------------------------------------------------------

#' Compute gate centers
#'
#' @param x gate object.
#' @param channels channels used to construct the plot.
#' @param text_x x co-ordinate for label.
#' @param text_y y co-ordinate for label.
#'
#' @return matrix with "x" and "y" colnames containing the gate center
#'   co-ordinates.
#'
#' @importFrom graphics par
#' @importFrom flowCore rectangleGate
#' @importFrom methods is
#'
#' @noRd
.cyto_gate_center <- function(x, ...) {
  UseMethod(".cyto_gate_center")
}

#' @noRd
.cyto_gate_center.default <- function(x,
                                      channels,
                                      text_x = NA,
                                      text_y = NA) {
  
  # CANNOT COMPUTE GATE CENTER - FILTER
  matrix(c(text_x, text_y),
         ncol = 2,
         byrow = FALSE,
         dimnames = list(NULL, c("x","y")))
  
}

#' @noRd
.cyto_gate_center.rectangleGate <- function(x,
                                            channels,
                                            text_x = NA,
                                            text_y = NA) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng/1.04)/2 # 2% BUFFER EITHER SIDE
  xmin <- xmin + 0.5 * xpad # 1% BUFFER
  xmax <- xmax - 0.5 * xpad # 1% BUFFER
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng/1.04)/2 # 2% BUFFER EITHER SIDE
  ymin <- ymin + 0.5 * ypad # 1% BUFFER
  ymax <- ymax - 0.5 * ypad # 1% BUFFER
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # CORRECT DIMENSIONS
  x <- cyto_gate_convert(x, channels = channels)
  
  # GATE CENTER ----------------------------------------------------------------
  
  # 1D GATE - 1D PLOT
  if(length(channels) == 1){
    # X COORD
    if(.all_na(text_x)){
      gate_xmin <- x@min
      gate_xmax <- x@max
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_xmin)){
        gate_xmin <- xmin
      }
      if(is.infinite(gate_xmax)){
        gate_xmax <- xmax
      }
    }
    # Y COORD
    if(.all_na(text_y)){
      gate_ymin <- ymin
      gate_ymax <- ymax
    }
    
    # 2D GATE - 2D PLOT
  }else if(length(channels) == 2){
    # X COORD
    if(.all_na(text_x)){
      gate_xmin <- x@min[channels[1]]
      gate_xmax <- x@max[channels[1]]
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_xmin)){
        gate_xmin <- xmin
      }
      if(is.infinite(gate_xmax)){
        gate_xmax <- xmax
      }
    }
    # Y COORD
    if(.all_na(text_y)){
      gate_ymin <- x@min[channels[2]]
      gate_ymax <- x@max[channels[2]]
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_ymin)){
        gate_ymin <- ymin
      }
      if(is.infinite(gate_ymax)){
        gate_ymax <- ymax
      }
    }
  }
  
  # GATE CENTER ----------------------------------------------------------------
  
  # GATE CENTER -  X COORD
  if(.all_na(text_x)){
    text_x <- (gate_xmin + gate_xmax) / 2
  }
  
  # GATE CENTER - Y COORD
  if(.all_na(text_y)){
    text_y <- (gate_ymin + gate_ymax) / 2
  }
  
  # RETURN GATE CENTER ---------------------------------------------------------
  
  # GATE CENTER MATRIX
  text_xy <- matrix(c(text_x, text_y),
                    ncol = 2,
                    byrow = FALSE)
  colnames(text_xy) <- c("x", "y")
  
  # RETURN GATE CENTER MATRIX
  return(text_xy)
}

#' @noRd
.cyto_gate_center.polygonGate <- function(x,
                                          channels,
                                          text_x = NA,
                                          text_y = NA) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng/1.04)/2 # 2% BUFFER EITHER SIDE
  xmin <- xmin + 0.5 * xpad # 1% BUFFER
  xmax <- xmax - 0.5 * xpad # 1% BUFFER
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng/1.04)/2 # 2% BUFFER EITHER SIDE
  ymin <- ymin + 0.5 * ypad # 1% BUFFER
  ymax <- ymax - 0.5 * ypad # 1% BUFFER
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # CORRECT DIMENSIONS
  x <- cyto_gate_convert(x, channels = channels)
  
  # GATE CENTER ----------------------------------------------------------------
  
  # GATE CENTER X COORD
  if (.all_na(text_x)) {
    # 2D RECTANGLEGATE
    if(is(x, "rectangleGate")){
      gate_xmin <- x@min[channels[1]]
      gate_xmax <- x@max[channels[1]]
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_xmin)){
        gate_xmin <- xmin
      }
      if(is.infinite(gate_xmax)){
        gate_xmax <- xmax
      }
      text_x <- (gate_xmin + gate_xmax)/2
      # 2D POLYGONGATE  
    }else{
      # COORDS
      coords <- x@boundaries[, channels[1]]
      # REPLACE INFINITE COORDS
      if(any(is.infinite(coords))){
        # -Inf
        if(coords[is.infinite(coords)] < 0){
          coords[is.infinite(coords)] <- xmin
        }
        # Inf
        if(coords[is.infinite(coords)] > 0){
          coords[is.infinite(coords)] <- xmax
        }
      }
      text_x <- sum(coords) / length(coords)
    }
  }
  
  # GATE CENTER Y COORD
  if (.all_na(text_y)) {
    # 2D RECTANGLEGATE
    if(is(x, "rectangleGate")){
      # 1D PLOT 
      if(length(channels) == 1){
        gate_ymin <- ymin
        gate_ymax <- ymax
      }else{
        gate_ymin <- x@min[channels[2]]
        gate_ymax <- x@max[channels[2]]
      }
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_ymin)){
        gate_ymin <- ymin
      }
      if(is.infinite(gate_ymax)){
        gate_ymax <- ymax
      }
      text_y <- (gate_ymin + gate_ymax)/2
      # 2D POLYGONGATE  
    }else{
      # COORDS
      coords <- x@boundaries[, channels[2]]
      # REPLACE INFINITE COORDS
      if(any(is.infinite(coords))){
        # -Inf
        if(coords[is.infinite(coords)] < 0){
          coords[is.infinite(coords)] <- ymin
        }
        # Inf
        if(coords[is.infinite(coords)] > 0){
          coords[is.infinite(coords)] <- ymax
        }
      }
      text_y <- sum(coords) / length(coords)
    }
  }
  
  # RETURN GATE CENTER ---------------------------------------------------------
  
  # GATE CENTER MATRIX
  text_xy <- matrix(c(text_x, text_y),
                    ncol = 2,
                    byrow = FALSE
  )
  colnames(text_xy) <- c("x", "y")
  
  # RETURN GATE CENTER MATRIX
  return(text_xy)
}

#' @noRd
.cyto_gate_center.ellipsoidGate <- function(x,
                                            channels,
                                            text_x = NA,
                                            text_y = NA) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng/1.04)/2 # 2% BUFFER EITHER SIDE
  xmin <- xmin + 0.5 * xpad # 1% BUFFER
  xmax <- xmax - 0.5 * xpad # 1% BUFFER
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng/1.04)/2 # 2% BUFFER EITHER SIDE
  ymin <- ymin + 0.5 * ypad # 1% BUFFER
  ymax <- ymax - 0.5 * ypad # 1% BUFFER
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # CORRECT DIMENSIONS
  x <- cyto_gate_convert(x, channels = channels)
  
  # GATE CENTER ----------------------------------------------------------------
  
  # GATE CENTER X COORD
  if(.all_na(text_x)){
    # 2D RECTANGLEGATE
    if(is(x, "rectangleGate")){
      gate_xmin <- x@min[channels[1]]
      gate_xmax <- x@max[channels[1]]
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_xmin)){
        gate_xmin <- xmin
      }
      if(is.infinite(gate_xmax)){
        gate_xmax <- xmax
      }
      text_x <- (gate_xmin + gate_xmax)/2
      # ELLIPSOIDGATE
    }else{
      text_x <- x@mean[channels[1]]
    }
  }
  
  # GATE CENTER Y COORD
  if(.all_na(text_y)){
    # 2D RECTANGLEGATE
    if(is(x, "rectangleGate")){
      # 1D PLOT
      if(length(channels) == 1){
        gate_ymin <- ymin
        gate_ymax <- ymax
      }else{
        gate_ymin <- x@min[channels[2]]
        gate_ymax <- x@max[channels[2]]
      }
      # REPLACE INFINITE COORDS
      if(is.infinite(gate_ymin)){
        gate_ymin <- ymin
      }
      if(is.infinite(gate_ymax)){
        gate_ymax <- ymax
      }
      text_y <- (gate_ymin + gate_ymax)/2
      # ELLIPSOIDGATE
    }else{
      text_y <- x@mean[channels[2]]
    }
  }
  
  # RETURN GATE CENTER ---------------------------------------------------------
  
  # GATE CENTER MATRIX
  text_xy <- matrix(c(text_x, text_y),
                    ncol = 2,
                    byrow = FALSE,
  )
  colnames(text_xy) <- c("x", "y")
  
  # RETURN GATE CENTER MATRIX
  return(text_xy)
}

#' @noRd
.cyto_gate_center.quadGate <- function(x,
                                       channels,
                                       text_x = NA,
                                       text_y = NA){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng/1.04)/2 # 2% BUFFER EITHER SIDE
  xmin <- xmin + 0.5 * xpad # 1% BUFFER
  xmax <- xmax - 0.5 * xpad # 1% BUFFER
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng/1.04)/2 # 2% BUFFER EITHER SIDE
  ymin <- ymin + 0.5 * ypad # 1% BUFFER
  ymax <- ymax - 0.5 * ypad # 1% BUFFER
  
  # PREPARE GATE ---------------------------------------------------------------
  
  # CORRECT DIMENSIONS
  x <- cyto_gate_convert(x, channels = channels)
  
  # REPEAT ARGUMENTS -----------------------------------------------------------
  
  # QUADGATE - 4 CENTERS
  text_x <- rep(text_x, length.out = 4)
  text_y <- rep(text_y, length.out = 4)
  
  # QUADGATE CO-ORDINATES ------------------------------------------------------
  
  x_coord <- x@boundary[channels[1]]
  y_coord <- x@boundary[channels[2]]
  
  # LABEL CENTERS --------------------------------------------------------------
  
  # UPDATE MISSING LABEL CO-ORDINATES
  lapply(seq_len(4), function(z){
    # QUADRANT 1 - TOP LEFT
    if(z == 1){
      # X COORD 
      if(.all_na(text_x[z])){
        text_x[z] <<- mean(c(x_coord, xmin))
      }
      # Y COORD
      if(.all_na(text_y[z])){
        text_y[z] <<- mean(c(y_coord, ymax))
      }
      # QUADRANT 2 - TOP RIGHT
    }else if(z == 2){
      # X COORD 
      if(.all_na(text_x[z])){
        text_x[z] <<- mean(c(x_coord, xmax))
      }
      # Y COORD
      if(.all_na(text_y[z])){
        text_y[z] <<- mean(c(y_coord, ymax))
      }
      # QUADRANT 3 - BOTTOM RIGHT
    }else if(z == 3){
      # X COORD 
      if(.all_na(text_x[z])){
        text_x[z] <<- mean(c(x_coord, xmax))
      }
      # Y COORD
      if(.all_na(text_y[z])){
        text_y[z] <<- mean(c(y_coord, ymin))
      }
      # QUADRANT 4 - BOTTOM LEFT
    }else if(z == 4){
      # X COORD 
      if(.all_na(text_x[z])){
        text_x[z] <<- mean(c(x_coord, xmin))
      }
      # Y COORD
      if(.all_na(text_y[z])){
        text_y[z] <<- mean(c(y_coord, ymin))
      }
    }
  })
  
  # RETURN GATE CENTERS --------------------------------------------------------
  
  # GATE CENTERS MATRIX
  text_xy <- as.matrix(cbind(text_x, text_y))
  colnames(text_xy) <- c("x","y")
  
  # RETURN GATES CENTERS MATRIX
  return(text_xy)
  
}

#' @noRd
.cyto_gate_center.filters <- function(x,
                                      channels,
                                      text_x = NA,
                                      text_y = NA) {
  
  # LIST GATE OBJECTS ----------------------------------------------------------
  x <- unlist(x)
  
  # CALL LIST METHOD -----------------------------------------------------------
  text_xy <- .cyto_gate_center(x,
                               channels = channels,
                               text_x = text_x,
                               text_y = text_y
  )
  
  # RETURN GATE CENTER MATRIX --------------------------------------------------
  return(text_xy)
  
}

#' @noRd
.cyto_gate_center.list <- function(x,
                                   channels,
                                   text_x = NA,
                                   text_y = NA) {
  
  # LIST GATE OBJECTS ----------------------------------------------------------
  
  # WATCH OUT FOR FILTERS
  x <- unlist(x)
  
  # REPEAT ARGUMENTS -----------------------------------------------------------
  
  # GATE COUNT
  NG <- .cyto_gate_count(x, total = FALSE)
  TNG <- sum(NG)
  
  # REPEAT ARGUMENTS
  text_x <- rep_len(text_x, TNG)
  text_y <- rep_len(text_y, TNG)
  
  # SPLIT ARGUMENTS
  text_x <- split(text_x, rep(seq_len(length(x)), times = NG))
  text_y <- split(text_y, rep(seq_len(length(x)), times = NG))
  
  # GATE CENTERS ---------------------------------------------------------------
  
  # LIST GATE CENTER MATRICES
  text_xy <- mapply(function(x,
                             text_x,
                             text_y) {
    .cyto_gate_center(x,
                      channels = channels,
                      text_x = text_x,
                      text_y = text_y
    )
  }, x,
  text_x,
  text_y,
  SIMPLIFY = FALSE
  )
  
  # # GATE CENTER MATRIX
  # text_xy <- do.call("rbind", text_xy)
  
  # RETURN GATE CENTER MATRIX
  return(text_xy)
}

## .CYTO_GATE_COUNT ------------------------------------------------------------

#' Compute number of gated populations
#'
#' @param gate list of gate objects including negated filters.
#' @param drop logical indicating whether non-gate objects should be dropped
#'   from the count, set to FALSE by default.
#' @param total return the total number of gated populations
#'
#' @importFrom methods is
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_gate_count <- function(gate = NA,
                             drop = FALSE,
                             total = TRUE){
  
  # NO GATE - SINGLE POPULATION
  if(.all_na(gate)){
    return(1)
  }
  
  # POPULATION COUNT
  P <- LAPPLY(seq_len(length(gate)), function(z){
    # GATE OBJECT
    if(grepl("gate", cyto_class(gate[[z]]), ignore.case = TRUE)) {
      if(cyto_class(gate[[z]], "quadGate")) {
        return(4)
      } else {
        return(1)
      }
    # FILTER  
    } else {
      if(drop) {
        return(NULL)
      } else {
        return(1)
      }
    }
  })
  
  # RETURN TOTAL POPULATIONS
  if(total == TRUE){
    return(sum(P))
  }else{
    return(P)
  }
  
}

## .CYTO_GATE_QUAD_CONVERT -----------------------------------------------------


#' Convert between quadGate to rectangleGates 
#' @return list of rectangleGates or a quadGate.
#' @importFrom flowCore rectangleGate quadGate parameters
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @noRd
.cyto_gate_quad_convert <- function(gate,
                                    channels) {
  # RECTANGLEGATES SUPPLIED
  if(cyto_class(gate, "list", TRUE)) {
    # LIST OF RECTANGLEGATES
    if(all(LAPPLY(gate, function(z){
      cyto_class(z, "rectangleGate")
    }))){
      # ORDER QUADRANTS - GATINGSET EXTRACTED GATES
      if(all(LAPPLY(gate, function(z){
        any(grepl("quad", names(attributes(z))))
      }))){
        # QUADRANTS
        quads <- names(attributes(gate[[1]])[["quadrants"]])
      }else{
        quads <- LAPPLY(gate, function(z){
          z@filterId
        })
      }
      # CHANNELS
      chans <- as.character(parameters(gate[[1]]))
      # CO-ORDINATES
      coords <- .cyto_gate_coords(gate, channels = chans)
      coords <- LAPPLY(chans, function(z){
        unique(coords[, z][is.finite(coords[, z])])
      })
      names(coords) <- chans
      # NORMAL GATES
      if(length(coords) > 2){
        return(gate)
      }
      # QUADGATE
      qg <- quadGate(filterId = paste(quads, collapse = "|"),
                     .gate = coords)
      return(qg)
    }
  }else if(cyto_class(gate, "quadGate")){
    # INPUT CO-ORDINATE
    xcoord <- gate@boundary[1]
    ycoord <- gate@boundary[2]
    # ALIAS
    alias <- unlist(strsplit(gate@filterId, "\\|"))
    if(length(alias) != 4){
      alias <- c("Q1", "Q2", "Q3","Q4")
    }
    # TOP LEFT
    coords <- list(c(-Inf, xcoord), c(ycoord, Inf))
    names(coords) <- channels
    Q1 <- rectangleGate(coords, filterId = alias[1])
    # TOP RIGHT
    coords <- list(c(xcoord, Inf), c(ycoord, Inf))
    names(coords) <- channels
    Q2 <- rectangleGate(coords, filterId = alias[2])
    # BOTTOM RIGHT
    coords <- list(c(xcoord, Inf), c(-Inf, ycoord))
    names(coords) <- channels
    Q3 <- rectangleGate(coords, filterId = alias[3])
    # BOTTOM LEFT
    coords <- list(c(-Inf, xcoord), c(-Inf, ycoord))
    names(coords) <- channels
    Q4 <- rectangleGate(coords, filterId = alias[4])
    Q <- list(Q1, Q2, Q3, Q4)
    names(Q) <- alias
    return(Q)
  }
  return(gate)
}

## .CYTO_GATE_COORDS -----------------------------------------------------------

#' Extract gate co-ordinates from a list of gate objects
#' 
#' @param x list of gate objects
#' @param channels vector of channel names used to construct the plot.
#' 
#' @importFrom flowCore parameters
#' @importFrom methods as is
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.cyto_gate_coords <- function(x, channels){
  
  # LIST OF GATES
  coords <- lapply(channels, function(z){
    LAPPLY(x, function(y){
      # RECTANGLEGATE
      if(cyto_class(y, "rectangleGate")){
        if(length(y@min) == 2){
          coords <- as.numeric(c(y@min[z], y@max[z]))
        }else{
          coords <- as.numeric(c(y@min, y@max))
        }
      # POLYGONGATE
      }else if(cyto_class(y, "polygonGate")){
        coords <- as.numeric(y@boundaries[, z])
      # ELLIPSOIDGATE
      }else if(cyto_class(y, "ellipsoidGate")){
        y <- as(y, "polygonGate")
        coords <- as.numeric(y@boundaries[, z])
      # QUADGATE
      }else if(cyto_class(y, "quadGate")){
        coords <- as.numeric(y@boundary[z])
      # NOT SUPPORTED
      } else {
        coords <- NULL
      }
      return(coords)
    })
  })
  coords[sapply(coords, "is.null")] <- NULL
  
  # COORD MATRIX
  if(length(channels) == 1){
    coords <- matrix(coords[[1]], ncol = 1)
    colnames(coords) <- channels
  }else if(length(channels) == 2){
    coords <- as.matrix(do.call("cbind", coords))
    colnames(coords) <- channels
  }
  
  # GATE COORDS
  return(coords)
  
}
