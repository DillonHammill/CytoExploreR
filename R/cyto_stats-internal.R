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
