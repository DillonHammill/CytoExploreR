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
                    "skewness",
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
#' @param x a vector or matrix
#' @noRd
cyto_stat_count <- function(x, 
                            ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c("count" = length(x))
    )
  # MATRIX
  } else {
    return(
      c("count" = nrow(x))
    )
  }
  
}

## MEAN ------------------------------------------------------------------------

#' Mean
#' @param x a vector or matrix
#' @param round numeric
#' @noRd
cyto_stat_mean <- function(x, 
                           round = 2,
                           ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c("mean" = round(
        mean(x, na.rm = TRUE, ...),
        round)
      )
    )
  # MATRIX - COLMEANS FOR SPEED
  } else {
    return(
      round(
        colMeans(
          x, 
          na.rm = TRUE,
          ...
        ), 
        round
      )
    )
  }
  
}

## GEOMETRIC MEAN --------------------------------------------------------------

#' Geometric Mean
#' @param x a linear vector or matrix
#' @param round numeric
#' @noRd
cyto_stat_geomean <- function(x,
                              round = 2,
                              ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c("geomean" = suppressWarnings(round(exp(mean(log(x))), round)))
    )
  # MATRIX
  } else {
    return(
      suppressWarnings(round(exp(colMeans(log(x))), round))
    )
  }
  
}

## MEDIAN ----------------------------------------------------------------------

#' Median
#' @param x a vector or matrix
#' @param round numeric
#' @importFrom robustbase colMedians
#' @noRd
cyto_stat_median <- function(x,
                             round = 2,
                             ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      "median" = round(median(x, na.rm = TRUE, ...), round)
    )
  # MATRIX
  } else {
    return(
      round(colMedians(x, na.rm = TRUE, ...), round)
    )
  }
  
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
                           limits = c(NA, NA),
                           round = 2,
                           ...) {
  
  # PREPARE LIMITS
  if(is.null(dim(limits))) {
    # LIST
    if(cyto_class(limits, "list")) {
      limits <- do.call(
        "cbind",
        limits
      )
      # VECTOR
    } else {
      limits <- matrix(
        limits,
        ncol = ifelse(
          is.null(dim(x)),
          1,
          ncol(x)
        ),
        nrow = 2,
        dimnames = list(
          NULL,
          if(is.null(dim(x))) {
            NULL
          } else{
            colnames(x)
          }
        )
      )
    }
  }
  
  # VECTOR
  if(is.null(dim(x))) {
    kde <- cyto_stat_density(
      x,
      smooth = smooth,
      bandwidth = bandwidth,
      bins = bins,
      limits = limits[, 1],
      ...
    )
    return(
      if(.all_na(kde)) {
        c("mode" = NA)
      } else {
        c(
          "mode" = round(
            kde$x[kde$y == max(kde$y)],
            round
          )
        )
      }
    )
  # MATRIX
  } else {
    # SORT LIMITS
    if(all(colnames(limits) %in% colnames(x))) {
      limits <- limits[, colnames(x), drop = FALSE]
    }
    # REPEAT ARGUMENTS
    bins <- rep(bins, length.out = ncol(x))
    bandwidth <- rep(bandwidth, length.out = ncol(x))
    smooth <- rep(smooth, length.out = ncol(x))
    round <- rep(round, length.out = ncol(x))
    # APPLY CYTO_STAT_MODE()
    cnt <- 0
    apply(
      x,
      2, 
      function(z){
        cnt <<- cnt + 1
        cyto_stat_mode(
          z,
          smooth = smooth[cnt],
          bandwidth = bandwidth[cnt],
          bins = bins[cnt],
          limits = limits[, cnt, drop = FALSE],
          round = round,
          ...
        )
      }
    )
  }
  
}

## STANDARD DEVIATION ----------------------------------------------------------

#' Standard Deviation
#' @param x a vector or matrix
#' @param round numeric
#' @importFrom stats sd
#' @noRd
cyto_stat_sd <- function(x,
                         round = 2,
                         ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c(
        "sd" = round(
          sd(x, na.rm = TRUE, ...),
          round
        )
      )
    )
  # MATRIX
  } else {
    apply(
      x,
      2,
      "cyto_stat_sd",
      round = round
    )
  }

}

## ROBUST STANDARD DEVIATION ---------------------------------------------------

#' Robust Standard Deviation
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_rsd <- function(x, 
                          round = 2,
                          ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c(
        "rsd" = round(
          median(
            abs(x - median(x))
          ) * 1.4826,
          round
        )
      )
    )
  # MATRIX
  } else {
    # MEDIANS - COLMEDIANS FOR SPEED
    md <- cyto_stat_median(x)
    # ROBUST STANDARD DEVIATIONS
    cnt <- c(0)
    apply(
      x,
      2,
      function(z){
        cnt <<- cnt + 1
        round(
          median(
            abs(z - md[cnt])
          ) * 1.4826,
          round
        )
      }
    )
  }

}

## COEFFICIENT OF VARIATION ----------------------------------------------------

#' Coefficient of Varaition
#' @param x a vector or matrix
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
  cv <- round(sd/md * 100, round)
  
  # NAMES
  if(is.null(dim(x))) {
    names(cv) <- "cv"
  } 
  return(cv)
  
} 

## ROBUST COEFFICIENT OF VARIATION ---------------------------------------------

#' Robust Coefficient of Varaition
#' @param x a vector or matrix
#' @param round numeric
#' @noRd
cyto_stat_rcv <- function(x, 
                          round = 2,
                          ...){
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c(
        "rcv" = round(
          (median(abs(x - median(x))) * 1.4826)/median(x) * 100
          , round
        )
      )
    )
  # MATRIX
  } else {
    # MEDIANS
    md <- cyto_stat_median(x)
    cnt <- c(0)
    apply(
      x,
      2,
      function(z){
        round(
          (median(abs(z - md[cnt])) * 1.4826)/md[cnt] * 100
          , round
        )
      }
    )
  }
  
} 

## QUANTILES -------------------------------------------------------------------

#' Quantiles
#' @param x a vector or matrix
#' @param round numeric
#' @noRd
cyto_stat_quantile <- function(x,
                               round = 2,
                               ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      round(
        quantile(
          x,
          na.rm = TRUE,
          ...
        ), 
        round
      )
    )
  # MATRIX
  } else {
    return(
      apply(
        x,
        2,
        "cyto_stat_quantile",
        round = round,
        ...
      )
    )
  }
  
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
                          limits = c(NA, NA),
                          method = "natural",
                          min = NULL,
                          max = NULL,
                          ...) {
  
  # PREPARE LIMITS
  if(is.null(dim(limits))) {
    # LIST
    if(cyto_class(limits, "list")) {
      limits <- do.call(
        "cbind",
        limits
      )
      # VECTOR
    } else {
      limits <- matrix(
        limits,
        ncol = ifelse(
          is.null(dim(x)),
          1,
          ncol(x)
        ),
        nrow = 2,
        dimnames = list(
          NULL,
          if(is.null(dim(x))) {
            NULL
          } else{
            colnames(x)
          }
        )
      )
    }
  }
  
  # VECTOR
  if(is.null(dim(x))) {
    # KERNEL DENSITY ESTIMATE
    kde <- cyto_stat_density(
      x,
      smooth = smooth,
      stat = "count",
      bandwidth = bandwidth,
      bins = bins,
      limits = limits[, 1]
    )
    # AREA UNDER CURVE
    if(.all_na(kde)) {
      return(
        c("auc" = NA)
      )
    }
    return(
      c(
        "auc" = round(
          integrate(
            splinefun(x$x, 
                      x$y, 
                      method = method),
            lower = if(.all_na(min)) {
              min(x$x, na.rm = TRUE)
            } else {
              min
            },
            upper = if(.all_na(max)) {
              max(x$x, na.rm = TRUE)
            } else {
              max
            },
            subdivisions = 2000,
            ...
          )$value, round)
      )
    )
  # MATRIX
  } else {
    # SORT LIMITS
    if(all(colnames(limits) %in% colnames(x))) {
      limits <- limits[, colnames(x), drop = FALSE]
    }
    # REPEAT ARGUMENTS
    smooth <- rep(smooth, length.out = ncol(x))
    bandwidth <- rep(bandwidth, length.out = ncol(x))
    bins <- rep(bins, length.out = ncol(x))
    method <- rep(method, length.out = ncol(x))
    min <- rep(min, length.out = ncol(x))
    max <- rep(max, length.out = ncol(x))
    # APPLY CYTO_STAT_AUC()
    cnt <- 0
    apply(
      x,
      2,
      function(z){
        cnt <<- cnt + 1
        cyto_stat_auc(
          z,
          round = round,
          smooth = smooth[cnt],
          bandwidth = bandwidth[cnt],
          bins = bins[cnt],
          limits = limits[, cnt, drop = FALSE],
          method = method[cnt],
          min = min[cnt],
          max = max[cnt],
          ...
        )
      }
    )
  }
  
}

## RANGE -----------------------------------------------------------------------

#' Range
#' @param x a matrix
#' @param round numeric
#' @noRd
cyto_stat_range <- function(x,
                            round = 2,
                            ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      suppressWarnings(
        round(
          structure(
            range(
              x,
              na.rm = TRUE,
              ...
            ),
            names = c("min", "max")
          ),
          round
        )
      )
    )
  # MATRIX
  } else {
    return(
      apply(
        x,
        2,
        "cyto_stat_range",
        round = round,
        ...
      )
    )
  }
}

## DENSITY ---------------------------------------------------------------------

# BANDWIDTH COMPUTED AT CYTOFRAME/CYTOSET LEVEL USING INSTRUMENT RANGE.

#' Density
#' @param x a vector or matrix
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
                              limits = c(NA, NA),
                              bandwidth = NA,
                              ...) {
  
  # PREPARE LIMITS
  if(is.null(dim(limits))) {
    # LIST
    if(cyto_class(limits, "list")) {
      limits <- do.call(
        "cbind",
        limits
      )
      # VECTOR
    } else {
      limits <- matrix(
        limits,
        ncol = ifelse(
          is.null(dim(x)),
          1,
          ncol(x)
        ),
        nrow = 2,
        dimnames = list(
          NULL,
          if(is.null(dim(x))) {
            NULL
          } else{
            colnames(x)
          }
        )
      )
    }
  }

  # VECTOR
  if(is.null(dim(x))) {
    # X - NUMERIC
    x <- as.numeric(x)
    # TOO FEW EVENTS FOR KDE
    if(length(x) < 2) {
      warning(
        "Insufficient events to compute kernel density!"
      )
      return(
        NA
      )
    }
    # CONVERT BINS TO BANDWIDTH
    if(.all_na(bandwidth)) {
      # BINS
      if(length(x) == 0) {
        bandwidth <- 0
      } else {
        if(.all_na(limits[, 1])) {
          rng <- range(x, na.rm = TRUE)
        } else {
          rng <- limits[, 1]
        }
        bandwidth <- diff(rng) / bins
      }
    }
    # SMOOTH - INACCURATE COUNTS OTHERWISE
    if(smooth < 1) {
      stop("'smooth' must be greater than or equal to 1!")
    }
    # RESTRICT DATA TO LIMITS - DATA OUTSIDE PLOT LIMITS MESSES UP BANDWIDTH
    if(!.all_na(limits[, 1])) {
      x <- x[x > min(limits[, 1], na.rm = TRUE) &
               x < max(limits[, 1], na.rm = TRUE)]
    }
    # KERNEL DENSITY
    kd <- stats::density(
      x[!is.na(x)],
      adjust = smooth, 
      bw = bandwidth, 
      ...
    )
    # PERCENT
    if(grepl("percent", stat, ignore.case = TRUE)){
      # MODAL
      kd$y <- (kd$y / max(kd$y)) * 100
      # ATTACH RANGE
      kd$range <- c(0, 100)
      # COUNT
    }else if(grepl("count", stat, ignore.case = TRUE)){
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
  # MATRIX
  } else {
    # SORT LIMITS
    if(all(colnames(limits) %in% colnames(x))) {
      limits <- limits[, colnames(x), drop = FALSE]
    }
    # REPEAT ARGUMENTS
    bins <- rep(bins, length.out = ncol(x))
    bandwidth <- rep(bandwidth, length.out = ncol(x))
    smooth <- rep(smooth, length.out = ncol(x))
    stat <- rep(stat, length.out = ncol(x))
    # APPLY CYTO_STAT_DENSITY()
    cnt <- 0
    apply(
      x,
      2, 
      function(z) {
        cnt <<- cnt + 1
        cyto_stat_density(
          z,
          stat = stat[cnt],
          smooth = smooth[cnt],
          bandwidth = bandwidth[cnt],
          limits = limits[, cnt, drop = FALSE],
          ...
        )
      }
    )
  }
  
}

## CYTO_STAT_BIN ---------------------------------------------------------------

#' Bin cytometry data
#'
#' @param x a vector or matrix of values to bin.
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
                          limits = c(NA, NA),
                          ...) {
  
  # PREPARE LIMITS
  if(is.null(dim(limits))) {
    # LIST
    if(cyto_class(limits, "list")) {
      limits <- do.call(
        "cbind",
        limits
      )
    # VECTOR
    } else {
      limits <- matrix(
        limits,
        ncol = ifelse(
          is.null(dim(x)),
          1,
          ncol(x)
        ),
        nrow = 2,
        dimnames = list(
          NULL,
          if(is.null(dim(x))) {
            NULL
          } else{
            colnames(x)
          }
        )
      )
    }
  }
  
  # VECTOR
  if(is.null(dim(x))) {
    # X - NUMERIC
    x <- as.numeric(x)
    # MINIMUM
    xmin <- min(limits[, 1])
    if(.all_na(xmin)) {
      xmin <- min(x, na.rm = TRUE)
    }
    # MAXIMUM
    xmax <- max(limits[, 1])
    if(.all_na(xmax)) {
      xmax <- max(x, na.rm = TRUE)
    }
    # BREAKS
    breaks <- seq(
      xmin,
      xmax,
      (xmax - xmin) / bins,
    )
    # USE CUT() TO CREATE BINS
    xbin <- cut(
      x, 
      breaks = breaks,
      labels = seq_len(bins),
      include.lowest = TRUE
    )
    # COUNTS PER BIN
    xbin <- table(
      xbin
    )
    # COUNTS -> FREQUENCY
    if(type != "count") {
      xbin <- xbin/length(x)
    }
    return(xbin)
  # MATRIX
  } else {
    # SORT LIMITS
    if(all(colnames(limits) %in% colnames(x))) {
      limits <- limits[, colnames(x), drop = FALSE]
    }
    # APPLY CYTO_STAT_BIN OVER COLUMNS
    cnt <- c(0)
    apply(
      x,
      2,
      function(z){
        cnt <<- cnt + 1
        cyto_stat_bin(
          z,
          bins = bins,
          type = type,
          limits = limits[, cnt, drop = FALSE]
        )
      }
    )
  }
  
}

## CYTO_STAT_SCALE -------------------------------------------------------------

#' Channel-wise re-scaling for cytometry data
#'
#' @param x object of class \code{matrix}.
#' @param type indicates the type of re-scaling to perform, options include
#'   \code{"range"}, \code{"mean"}, \code{"median"} or \code{"zscore"}.
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
  
  # VECTOR
  if(is.null(dim(x))) {
    # RANGE
    if(grepl("^r", type, ignore.case = TRUE)) {
      # DATA LIMITS
      xmin <- min(x, na.rm = TRUE)
      xmax <- max(x, na.rm = TRUE)
      # RANGE SCALING
      return((x - xmin)/(xmax - xmin))
      # MEAN
    } else if(grepl("^mea", type, ignore.case = TRUE)) {
      # DATA LIMITS
      xmin <- min(x, na.rm = TRUE)
      xmax <- max(x, na.rm = TRUE)
      # MEAN
      xmean <- mean(x, na.rm = TRUE)
      # MEAN SCALING
      return((x - xmean)/(xmax-xmin))
      # MEDIAN
    } else if(grepl("^med", type, ignore.case = TRUE)) {
      # DATA LIMITS
      xmin <- min(x, na.rm = TRUE)
      xmax <- max(x, na.rm = TRUE)
      # MEAN
      xmed <- median(x, na.rm = TRUE)
      # MEAN SCALING
      return((x - xmed)/(xmax-xmin))
      # Z SCORE
    } else if(grepl("^z", type, ignore.case = TRUE)) {
      # MEAN 
      xmean <- mean(x, na.rm = TRUE)
      # SD
      xsd <- sd(x, na.rm = TRUE)
      # Z-SCORE SCALING
      return((x - xmean)/xsd)
      # UNSUPPORTED SCALING METHOD
    } else {
      stop(
        "'type' must be either 'range', 'mean' or 'zscore'!"
      )
    }
  } else {
    apply(
      x,
      2,
      "cyto_stat_scale"
    )
  }
  
}

## CYTO_STAT_SKEWNESS ----------------------------------------------------------

#' Compute skewness of distributions
#' 
#' @param x a matrix.
#' @param ... not in use.
#' 
#' @return matrix with skewness values per channel.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
cyto_stat_skewness <- function(x,
                               ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    x <- x[!is.na(x)]
    n <- length(x)
    return(
      (sum((x-mean(x))^3)/n) / 
        (sum((x-mean(x))^2)/n) ^ (3/2)
    )
  # MATRIX
  } else {
    apply(
      x,
      2,
      "cyto_stat_skewness"
    )
  }
  
}

## CYTO_STAT_BKDE --------------------------------------------------------------

#' Compute binned 2D kernel density estimate
#'
#' Modified version of \code{KernSmooth::bkde()} that can return 2D binned
#' counts as well as smooth kernel density estimate.
#'
#' @param x a 2D matrix.
#' @param bins number of bins to use for x and y values, set to 250 by default.
#' @param bandwidth a vector of length 2 containing the bandwidth to use for the
#'   kernel density estimate for each column in \code{x}, uses
#'   \code{KernSmooth::dpik()} to compute bandwidths if not manually supplied.
#' @param limits list containing the ranges of each column in x to truncate the
#'   grid.
#' @param smooth logical indicating whether binned counts should be smoothed
#'   using kernel density estimates, set to TRUE by default. If only binned
#'   counts are required set this argument to FALSE to prevent running KDE code.
#' @param ... not in use.
#'
#' @return list with slots counts, bkde and bins.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @import KernSmooth
#' @importFrom stats fft dnorm
#'
#' @noRd
cyto_stat_bkde2d <- function(x,
                             bins = 250,
                             bandwidth = c(NA, NA),
                             limits = list(NA, NA),
                             smooth = TRUE,
                             ...) {
  
  # SETUP PARAMETERS
  n <- nrow(x)
  M <- rep(bins, length.out = 2)
  h <- bandwidth
  tau <- 3.4             # bivariate normal kernel
  
  # LIMITS - MATRIX
  if(!is.null(dim(limits))) {
    limits <- lapply(1:ncol(limits), function(z){limits[, z, drop = TRUE]})
  }
  
  # DEFAULT LIMITS - MATCH KERNSMOOTH
  limits <- lapply(seq_along(limits), function(z){
    if(any(is.na((limits[[z]])))) {
      return(
        c(
          min(x[, z]) - 1.5 * h[z],
          max(x[, z]) + 1.5 * h[z]
        )
      )
    }
    return(limits[[z]])
  })
  
  # COMPUTE BANDWIDTH USING PLUGIN METHOD
  bandwidth <- LAPPLY(seq_along(bandwidth), function(z){
    if(.all_na(bandwidth[z])) {
      if(length(x[, z]) < 2) {
        return(NA)
      } else {
        return(
          suppressWarnings(
            dpik(
              x[, z],
              gridsize = M[z],
              range.x = limits[[z]],
              truncate = TRUE
            )
          )
        )
      }
    } else {
      # BANDWIDTH > 0
      if(bandwidth[z] <= 0) {
        stop("'bandwidth' must be strictly positive!")
      }
      return(
        bandwidth[z]
      )
    }
  })
  h <- bandwidth
  
  # COMPUTE GRID POINTS
  a <- LAPPLY(limits, "min")
  b <- LAPPLY(limits, "max")
  xpts <- seq(
    a[1],
    b[1],
    length = M[1]
  )
  ypts <- seq(
    a[2],
    b[2],
    length = M[2]
  )
  
  # LINEAR BINNING - INTERNAL LINBIN2D()
  cnts <- cyto_func_call(
    "KernSmooth:::linbin2D",
    list(
      X = x,
      gpoints1 = xpts,
      gpoints2 = ypts
    )
  )
  
  # KERNEL DENSITY ESIMATE SMOOTHING
  if(smooth & n >= 2) {
    # COMPUTE KERNEL WEIGHTS
    L <- c(0, 0)
    kapid <- list(0, 0)
    lapply(seq_len(2), function(z){
      L[z] <<- min(floor(tau*h[z]*(M[z]-1)/(b[z]-a[z])), M[z] - 1L)
      lvecid <- seq(0, L[z])
      facid <- (b[z] - a[z])/(h[z]*(M[z]-1L))
      w <- matrix(dnorm(lvecid*facid)/h[z])
      tot <- sum(c(w, rev(w[-1L]))) * facid * h[z]
      kapid[[z]] <<- w/tot
    })
    kapp <- kapid[[1L]] %*% (t(kapid[[2L]]))/n
    
    # # GRIDSIZE TOO SMALL
    # if(min(L) == 0) {
    #   warning(
    #     "Binning grid too coarse for current bandwidth: increase 'bins'."
    #   )
    # }
    
    # COMBINE WEIGHTS & COUNTS TO GET ESTIMATE (FFT)
    P <- 2^(ceiling(log(M+L)/log(2)))   # smallest powers of 2 >= M+L
    L1 <- L[1L] ; L2 <- L[2L]
    M1 <- M[1L] ; M2 <- M[2L]
    P1 <- P[1L] ; P2 <- P[2L]
    
    # WRAP AROUND VERSION OF KAPP
    rp <- matrix(0, P1, P2)
    rp[1L:(L1+1), 1L:(L2+1)] <- kapp
    if (L1) rp[(P1-L1+1):P1, 1L:(L2+1)] <- kapp[(L1+1):2, 1L:(L2+1)]
    if (L2) rp[, (P2-L2+1):P2] <- rp[, (L2+1):2]
    
    # ZERO PADDED COUNTS
    sp <- matrix(0, P1, P2)
    sp[1L:M1, 1L:M2] <- cnts
    
    # INVERT ELEMET-WISE PRODUCT FFTs - TRUNCATE & NORMALISE
    rp <- fft(rp)                       
    sp <- fft(sp)
    rp <- Re(fft(rp*sp, inverse = TRUE)/(P1*P2))[1L:M1, 1L:M2]
    
    # NON-NEGATIVE
    rp <- rp * matrix(as.numeric(rp>0), nrow(rp), ncol(rp))
    
    # RETURN
    return(
      list(
        counts = cnts,
        bkde = rp,
        bins = list("x" = xpts,
                    "y" = ypts)
      )
    )
  # COUNTS ONLY
  } else {
    return(
      list(
        counts = cnts,
        bkde = NA,
        bins = list("x" = xpts,
                    "y" = ypts)
      )
    )
  }
  
}

## CYTO_STAT_RESCALE -----------------------------------------------------------

#' Rescale values within a newly defined range
#'
#' @param x values to be rescaled in the form of a vector.
#' @param scale min max values of new scale.
#' @param limits desired range for data on the current scale, values outside
#'   this range will be set to these limits prior to rescaling.
#'
#' @return vector or matrix of rescaled values within range [0,1].
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
cyto_stat_rescale <- function(x,
                              scale = c(0,1),
                              limits = c(NA,NA)) {
  
  # MATRIX
  if(!is.null(dim(x))) {
    # PREPARE LIMITS VECTOR -> MATRIX
    if(is.null(dim(limits))) {
      limits <- suppressWarnings(
        matrix(
          limits,
          ncol = ncol(x),
          nrow = 2,
          dimnames = list(NULL, colnames(x))
        )
      )
    }
    # CHECK LIMITS MATRIX
    if(!setequal(dim(limits), c(2, ncol(x)))) {
      stop(
        "'limits' must a matrix of min/max values for each column in x!"
      )
    }
    # PREPARE SCALE VECTOR -> MATRIX
    if(is.null(dim(scale))) {
      scale <- suppressWarnings(
        matrix(
          scale,
          ncol = ncol(x),
          nrow = 2,
          dimnames = list(
            NULL,
            colnames(x)
          )
        )
      )
    }
    # CHECK SCALE MARIX
    if(!setequal(dim(scale), c(2, ncol(x)))) {
      stop(
        "'scale' must a matrix of min/max values for each column in x!"
      )
    }
    # RESCALE EACH COLUMN
    cnt <- 0
    return(
      apply(
        x,
        2,
        function(z){
          cnt <<- cnt + 1
          cyto_stat_rescale(
            z,
            scale = if(!is.null(colnames(x)) & !is.null(colnames(scale))) {
              scale[, colnames(x)[cnt]]
            } else {
              scale[, cnt]
            },
            limits = if(!is.null(colnames(x)) & !is.null(colnames(limits))) {
              limits[, colnames(x)[cnt]]
            } else {
              limits[, cnt]
            }
          )
        }
      )
    )
  # VECTOR
  } else {
    # COMPUTE LIMITS
    if(any(is.na(limits))) {
      limits[is.na(limits)] <- range(x)[is.na(limits)]
    }
    # RESTRICT
    x[x < min(limits)] <- min(limits)
    x[x > max(limits)] <- max(limits)
    # RESCALE
    x <- min(scale) + ((x-min(limits))/diff(limits))*diff(scale)
    # RETURN RESCALED DATA
    return(x)
  }
  
}
