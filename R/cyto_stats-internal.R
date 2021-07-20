# CYTO_STATS FUNCTIONS ---------------------------------------------------------

# All cyto_stat() functions below accept a pre-processed matrix prepared by
# cyto_apply() and round the computed statistics to 2 decimal places.

## DISPATCH --------------------------------------------------------------------

#' Prepare FUN to dispatch to cyto_stat function
#' @noRd
.cyto_stat_dispatch <- function(FUN) {
  
  if(is.character(FUN)) {
    if(any(LAPPLY(c("count",
                    "freq",
                    "mean",
                    "geomean",
                    "median",
                    "mode",
                    "sd",
                    "rsd",
                    "cv",
                    "rcv",
                    "quantile",
                    "auc",
                    "range",
                    "bin"), function(z){
                      grepl(paste0("^", z, "$"), FUN, ignore.case = TRUE)
                    }))) {
      FUN <- paste0("cyto_stat_", tolower(FUN))
    } else if(grepl("^percent$", FUN, ignore.case = TRUE)) {
      FUN <- "cyto_stat_freq"
    } else if(grepl("^quant$", FUN, ignore.case = TRUE)) {
      FUN <- "cyto_stat_quantile"
    }
  }
  return(FUN)
  
}

## COUNT -----------------------------------------------------------------------

#' Count
#' @param x a matrix
#' @noRd
cyto_stat_count <- function(x, 
                            ...) {
  return(structure(nrow(x), names = "count"))
}

## MEAN ------------------------------------------------------------------------

#' Mean
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_mean <- function(x, 
                           round = 2,
                           ...) {
  return(round(colMeans(x, na.rm = TRUE, ...), round))
}

## GEOMETRIC MEAN --------------------------------------------------------------

#' Geometric Mean
#' @param x a linear matrix
#' @param round numeric
#' @noRd
cyto_stat_geomean <- function(x,
                              round = 2,
                              ...) {
  
  res <- suppressWarnings(round(exp(colMeans(log(x))), round))
  # if(any(is.nan(res))){
  #   stop(
  #     "Geometric mean only works for data greater than zero."
  #   )
  # }
  return(res)
  
}

## MEDIAN ----------------------------------------------------------------------

#' Median
#' @param x a matrix
#' @param round numeric
#' @importFrom robustbase colMedians
#' @noRd
cyto_stat_median <- function(x,
                             round = 2,
                             ...) {
  return(round(colMedians(x, na.rm = TRUE, ...), round))
}

## MODE ------------------------------------------------------------------------

#' Mode - current scale
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_mode <- function(x,
                           smooth = 1,
                           bandwidth = NA,
                           bins = 256,
                           round = 2,
                           ...) {
  
  # DENSITY MATRIX
  res <- cyto_stat_density(x,
                           smooth = smooth,
                           bandwidth = bandwidth,
                           bins = bins, 
                           ...)
  
  # MODE
  res <- LAPPLY(res, function(z){
    round(z$x[z$y == max(z$y)], round)
  })
  names(res) <- colnames(x)
  return(res)
  
}

## STANDARD DEVIATION ----------------------------------------------------------

#' Standard Deviation
#' @param x a matrix
#' @param round numeric
#' @importFrom stats sd
#' @noRd
cyto_stat_sd <- function(x,
                         round = 2,
                         ...) {
  return(round(apply(x, 2, sd, na.rm = TRUE, ...), round))
}

## ROBUST STANDARD DEVIATION ---------------------------------------------------

#' Robust Standard Deviation
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_rsd <- function(x, 
                          round = 2,
                          ...) {
  
  # MEDIANS
  md <- cyto_stat_median(x)
  
  # ROBUST STANDARD DEVIATIONS
  res <- LAPPLY(names(md), function(z){
    round(median(abs(x[, z] - md[z])) * 1.4826, round)
  })
  names(res) <- names(md)
  return(res)
  
}

## COEFFICIENT OF VARIATION ----------------------------------------------------

#' Coefficient of Varaition
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_cv <- function(x, 
                         round = 2,
                         ...){
  
  # MEDIANS
  md <- cyto_stat_median(x)
  
  # STANDARD DEVIATIONS
  sd <- cyto_stat_sd(x)
  
  # COEFFICIENT OF VARIATION
  return(round(sd/md * 100, round))
  
} 

## ROBUST COEFFICIENT OF VARIATION ---------------------------------------------

#' Robust Coefficient of Varaition
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_rcv <- function(x, 
                          round = 2,
                          ...){
  
  # MEDIANS
  md <- cyto_stat_median(x)
  
  # ROBUST STANDARD DEVIATIONS - PREVENT DUPLICATE MEDIAN COMPUTATION
  rSD <- LAPPLY(names(md), function(z){
    median(abs(x[, z] - md[z])) * 1.4826
  })
  names(rSD) <- names(md)
  
  # ROBUST COEFFICIENT OF VARIATION
  res <- LAPPLY(names(md), function(z){
    round(rSD[z]/md[z] * 100, round)
  })
  names(res) <- names(md)
  return(res)
  
} 

## QUANTILES -------------------------------------------------------------------

#' Quantiles
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_quantile <- function(x,
                               round = 2,
                               ...) {
  
  return(round(apply(x, 2, quantile, na.rm = TRUE, ...), round))
  
}

## AREA UNDER CURVE ------------------------------------------------------------

#' Area Under Curve
#' @param x a matrix
#' @param round numeric
#' @importFrom stats integrate splinefun
#' @noRd
cyto_stat_auc <- function(x,
                          round = 2,
                          smooth = 1,
                          bandwidth = NA,
                          bins = 256,
                          limits = NA,
                          method = "natural",
                          min = NULL,
                          max = NULL,
                          ...) {
  
  # DENSITY - LIST
  d <- suppressWarnings(
    cyto_stat_density(x,
                      smooth = smooth,
                      stat = "count", # cyto_stat_compute conflict
                      bandwidth = bandwidth,
                      bins = bins,
                      limits = limits)
  )
  
  # AREA UNDER CURVE
  res <- LAPPLY(d, function(z){
    if(.all_na(z)) {
      return(NA)
    }
    round(
      integrate(
        splinefun(z$x, 
                  z$y, 
                  method = method),
        lower = switch(as.character(is.null(min)), 
                       "TRUE" = min(z$x), 
                       "FALSE" = min),
        upper = switch(as.character(is.null(max)),
                       "TRUE" = max(z$x),
                       "FALSE" = max),
        subdivisions = 2000,
        ...
      )$value, round)
  })
  names(res) <- names(d)
  return(res)
  
}

## RANGE -----------------------------------------------------------------------

#' Range
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_range <- function(x,
                            round = 2,
                            ...) {
  res <- suppressWarnings(round(apply(x, 2, range, ...), round))
  rownames(res) <- c("min", "max")
  return(res)
}

## DENSITY ---------------------------------------------------------------------

# BANDWIDTH COMPUTED AT CYTOFRAME/CYTOSET LEVEL USING INSTRUMENT RANGE.

#' Density
#' @param x a matrix
#' @param stat "percent", "density" or "count"
#' @param smooth numeric set to 1
#' @param bins 256
#' @param bandwidth width of bins
#' @param limits matrix with rows min and max for each column of x
#' @importFrom stats density
#' @noRd
cyto_stat_density <- function(x,
                              stat = "density",
                              smooth = 1,
                              bins = 256,
                              limits = NA,
                              bandwidth = NA,
                              ...) {
  
  # TOO FEW EVENTS
  if(nrow(x) <= 2){
    warning("Insufficient events to compute kernel density.")
    res <- rep(list(NA), ncol(x))
    names(res) <- colnames(x)
    return(res)
  }
  
  # CONVERT BINS TO BANDWIDTH
  if(.all_na(bandwidth)) {
    # BINS
    bins <- rep(bins, length.out = ncol(x))
    names(bins) <- colnames(x)
    # BANDWIDTH
    cnt <- 0
    bandwidth <- apply(x, 2, function(z){
      # COUNTER
      cnt <<- cnt + 1
      # BANDWIDTH
      if(length(z) == 0) {
        return(0)
      } else {
        if(.all_na(limits)) {
          rng <- range(z) # use data range
        } else {
          rng <- limits[, colnames(x)[cnt]]
        }
        return((rng[2] - rng[1])/bins)
      }
    })
  } else {
    if(length(bandwidth) != ncol(x)) {
      stop("A bandwidth is required for each channel!")
    }
  }
  
  # SMOOTH
  if(smooth < 1) {
    # SMOOTH < 1 RESULTS IN INACCURATE COUNTS
    stop("'smooth' must be greater than or equal to 1!")
  }
  smooth <- rep(smooth, ncol(x))
  names(smooth) <- colnames(x)
  
  # DENSITY - LIST
  cnt <- 0
  res <- apply(x, 2, function(z,
                              st = stat,
                              sm = smooth,
                              bw = bandwidth,
                              ...){
    # COUNTER
    cnt <<- cnt + 1
    # RESTRICT DATA TO LIMITS - DATA OUTSIDE PLOT LIMITS MESSES UP BANDWIDTH
    if(!.all_na(limits)) {
      z <- z[z > limits[, colnames(x)[cnt]][1] &
               z < limits[, colnames(x)[cnt]][2]]
    }
    # KERNEL DENSITY
    kd <- stats::density(z,
                         adjust = sm[cnt], # smoothing per channel
                         bw = bw[cnt], # bandwidth per channel
                         ...)
    
    # PERCENT
    if(grepl("percent", st, ignore.case = TRUE)){
      # MODAL
      kd$y <- (kd$y / max(kd$y)) * 100
      # ATTACH RANGE
      kd$range <- c(0, 100)
      # COUNT
    }else if(grepl("count", st, ignore.case = TRUE)){
      # COUNT
      kd$y <- kd$y * kd$n * kd$bw
      # ATTACH RANGE
      kd$range <- c(0, max(kd$y))
      # DENSITY  
    }else{
      # ATTACH RANGE
      kd$range <- c(0, max(kd$y))
    }
    return(kd)
  })
  
  # RETURN LIST OF DENSITY OBJECTS
  return(res)
  
}

## CYTO_STAT_BIN ---------------------------------------------------------------

#' Bin cytometry data
#'
#' @param x a matrix of values to bin.
#' @param bins the number of bins to use, set to 400 by default.
#' @param type whether to use \code{"count"} or \code{"freq"}.
#' @param limits matrix or named list containing the minimum and maximum values
#'   for each channel on the current scale, defaults to the data range if not
#'   supplied.
#' @param ... not in use.
#'
#' @return matrix of binned data per channel with either counts of frequencies.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
cyto_stat_bin <- function(x, 
                          bins = 400,
                          type = "count",
                          limits = NULL,
                          ...) {
  
  # BIN DATA
  cnt <- 0
  apply(
    x,
    2,
    function(y){
      # COUNTER
      cnt <<- cnt + 1
      # CHANNEL
      x_chan <- colnames(x)[cnt]
      # PREPARE LIMITS
      if(is.null(dim(limits))) {
        limits <- do.call("cbind", limits)
      }
      # MINIMUM
      x_min <- ifelse(
        is.numeric(limits[,x_chan][1]),
        limits[,x_chan][1],
        min(y)
      )
      # MAXIMUM
      x_max <- ifelse(
        is.numeric(limits[,x_chan][2]),
        limits[,x_chan][2],
        max(y)
      )
      # BREAKS
      breaks <- seq(
        x_min,
        x_max,
        (x_max - x_min) / bins,
      )
      # USE CUT() TO CREATE BINS
      x_bin <- cut(
        y, 
        breaks = breaks,
        labels = seq_len(bins),
        include.lowest = TRUE
      )
      # COUNTS PER BIN
      x_bin <- table(
        x_bin
      )
      # COUNTS -> FREQUENCY
      if(type != "count") {
        x_bin <- x_bin/length(x)
      }
      return(x_bin)
    }
  )
  
}

## CYTO_STAT_SCALE -------------------------------------------------------------

#' Channel-wise re-scaling for cytometry data
#'
#' @param x object of class \code{matrix}.
#' @param type indicates the type of re-scaling to perform, options include
#'   \code{"range"}, \code{"mean"} or \code{"zscore"}.
#' @param ... not in use.
#'
#' @return matrix with data re-scaled per channel.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
cyto_stat_scale <- function(x,
                            type = "range",
                            ...) {
  
  # SCALING
  apply(
    x,
    2, 
    function(z) {
      # RANGE
      if(grepl("^r", type, ignore.case = TRUE)) {
        # DATA LIMITS
        zmin <- min(z, na.rm = TRUE)
        zmax <- max(z, na.rm = TRUE)
        # RANGE SCALING
        return((z - zmin)/(zmax - zmin))
      # MEAN
      } else if(grepl("^m", type, ignore.case = TRUE)) {
        # DATA LIMITS
        zmin <- min(z, na.rm = TRUE)
        zmax <- max(z, na.rm = TRUE)
        # MEAN
        zmean <- mean(z, na.rm = TRUE)
        # MEAN SCALING
        return((z - zmean)/(zmax-zmin))
      # Z SCORE
      } else if(grepl("^z", type, ignore.case = TRUE)) {
        # MEAN 
        zmean <- mean(z, na.rm = TRUE)
        # SD
        zsd <- sd(z, na.rm = TRUE)
        # Z-SCORE SCALING
        return((z - zmean)/zsd)
      # UNSUPPORTED SCALING METHOD
      } else {
        stop(
          "'type' must be either 'range', 'mean' or 'zscore'!"
        )
      }
    }
  )
  
}
