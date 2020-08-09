# CYTO_STATS FUNCTIONS ---------------------------------------------------------

# All cyto_stat() functions below accept a pre-processed matrix prepared by
# cyto_apply() and round the computed statistics to 2 decimal places.

## DISPATCH --------------------------------------------------------------------

#' Prepare FUN to dispatch to cyto_stat function
#' @noRd
cyto_stat_dispatch <- function(FUN) {
  
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
                    "range"), function(z){
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
  if(any(is.nan(res))){
    stop(
      "Geometric mean only works for data greater than zero."
    )
  }
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
                           smooth = 0.6,
                           bandwidth = NULL,
                           round = 2,
                           ...) {
  
  # DENSITY MATRIX
  res <- cyto_apply(x,
                    "cyto_stat_density",
                    input = "matrix",
                    channels = channels,
                    smooth = smooth,
                    bandwidth = bandwidth,
                    simplify = TRUE,
                    ...)
  # MODE
  res[1, ] <- LAPPLY(seq_len(ncol(res)), function(z){
    round(z[[1]]$x[z[[1]]$y == max(z[[1]]$y)], round)
  })
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
                          smooth = 0.6,
                          bandwidth = NULL,
                          method = "natural",
                          ...) {
  
  # DENSITY - LIST
  d <- cyto_stat_density(x,
                         smooth = smooth,
                         stat = "count", # cyto_stat_compute conflict
                         bandwidth = bandwidth)
  
  # AREA UNDER CURVE
  res <- LAPPLY(d, function(z){
    round(
      integrate(
        splinefun(z$x, 
                  z$y, 
                  method = method),
        lower = min(z$x),
        upper = max(z$x),
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

## BANDWIDTH -------------------------------------------------------------------

#' Bandwidth
#' @param x a matrix
#' @param round numeric
#' @importFrom stats bw.nrd0
#' @noRd
cyto_stat_bandwidth <- function(x,
                                round = 2,
                                ...) {
  return(round(apply(x, 2, "bw.nrd0", ...), round))
}

## DENSITY ---------------------------------------------------------------------

#' Density
#' @param x a matrix
#' @param stat "percent", "density" or "count"
#' @param smooth numeric set to 0.6
#' @param bandwidth NULL
#' @importFrom stats density
#' @noRd
cyto_stat_density <- function(x,
                              stat = "density",
                              smooth = 0.6,
                              bandwidth = NULL,
                              ...) {
  
  # TOO FEW EVENTS
  if(nrow(x) <= 2){
    warning("Insufficient events to compute kernel density.")
    res <- rep(list(NA), ncol(x))
    names(res) <- colnames(x)
    return(res)
  }
  
  # BANDWIDTH
  if(is.null(bandwidth)){
    bandwidth <- rep("nrd0", ncol(x))
    names(bandwidth) <- colnames(x)
  }else{
    bandwidth <- rep(bandwidth, ncol(x))
    names(bandwidth) <- colnames(x)
  }
  
  # SMOOTH
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
    cnt <- get("cnt", envir = parent.frame(2)) + 1
    assign("cnt", cnt, envir = parent.frame(2))
    
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
      attributes(kd)$range <- c(0, 100)
      # COUNT
    }else if(grepl("count", st, ignore.case = TRUE)){
      # COUNT
      kd$y <- kd$y * kd$n
      # ATTACH RANGE
      attributes(kd)$range <- c(0, max(kd$y))
      # DENSITY  
    }else{
      # ATTACH RANGE
      attributes(kd)$range <- c(0, max(kd$y))
    }
    return(kd)
  })
  
  # RETURN LIST OF DENSITY OBJECTS
  return(res)
  
}
