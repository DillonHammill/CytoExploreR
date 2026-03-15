# CYTO_STATS FUNCTIONS ---------------------------------------------------------

# All cyto_stat() functions below accept a pre-processed matrix prepared by
# cyto_apply() and round the computed statistics to 2 decimal places.

# NOTE: for speed we don't handle NA values users will need to do that before
# passing their data to these functions.

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
                    "geomedian",
                    "mode",
                    "sd",
                    "rsd",
                    "cv",
                    "rcv",
                    "skewness",
                    "quantile",
                    "auc",
                    "range",
                    "bin",
                    "hex"), function(z){
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
                            parent = "root",
                            ...) {
  
  # FLOWSET - NROW NOT WORK FOR FLOWSETS
  if(cyto_class(x, "flowSet", TRUE)) {
    return(
      structure(
        sapply(
          seq_along(x),
          function(id) {
            nrow(x[[id]])
          }
        ),
        names = cyto_names(x)
      )
    )
  # CYTOSET OR GATINGSET
  } else if(cyto_class(x, c("cytoset", "GatingSet"), FALSE)) {
    if(cyto_class(x, "GatingSet")) {
      x <- gs_pop_get_data(x, parent)
    }
    return(
      unlist(
        nrow(x)
      )
    )
  # FLOWFRAME
  } else if(cyto_class(x, "flowFrame", FALSE)) {
    return(
      c("count" = nrow(x))
    )
  # VECTOR
  } else if(is.null(dim(x))) {
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
          mean_cpp(x),
          round
        )
      )
    )
  # MATRIX - COLMEANS FOR SPEED
  } else {
    return(
      round(
        col_mean_cpp(
          x
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
      c(
        "geomean" = suppressWarnings(
          round(
            geomean_cpp(x),
            round
          )
        )
      )
    )
  # MATRIX
  } else {
    return(
      suppressWarnings(
        round(
          col_geomean_cpp(x), 
          round
        )
      )
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
      c(
        "median" = round(
          median_cpp(x),
          round
        )
      )
    )
  # MATRIX
  } else {
    return(
      round(
        col_median_cpp(x),
        round
      )
    )
  }
  
}

## GEOMEDIAN -------------------------------------------------------------------

#' Geometric Median
#' @param x a vector or matrix
#' @param round numeric
#' @noRd
cyto_stat_geomedian <- function(x,
                                round = 2,
                                ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c(
        "geomedian" = round(
          unname(geometric_median_cpp(matrix(x, ncol = 1))),
          round
        )
      )
    )
    # MATRIX
  } else {
    return(
      round(
        geometric_median_cpp(x),
        round
      )
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
    future_apply(
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
          sd_cpp(x),
          round
        )
      )
    )
  # MATRIX
  } else {
    return(
      round(
        col_sd_cpp(x),
        round
      )
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
          rsd_cpp(x),
          round
        )
      )
    )
  # MATRIX
  } else {
    return(
      round(
        col_rsd_cpp(x),
        round
      )
    )
  }

}

## COEFFICIENT OF VARIATION ----------------------------------------------------

#' Coefficient of Variation
#' @param x a vector or matrix
#' @param round numeric
#' @noRd
cyto_stat_cv <- function(x, 
                         round = 2,
                         ...){
  
  # NOTE: CVs are decimals not percentages
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      c(
        "cv" = round(
          cv_cpp(x),
          round
        )
      )
    )
    # MATRIX
  } else {
    return(
      round(
        col_cv_cpp(x),
        round
      )
    )
  }
  
} 

## ROBUST COEFFICIENT OF VARIATION ---------------------------------------------

#' Robust Coefficient of Variation
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
          rcv_cpp(x),
          round
        )
      )
    )
  # MATRIX
  } else {
    return(
      round(
        col_rcv_cpp(x),
        round
      )
    )
  }
  
} 

## QUANTILES -------------------------------------------------------------------

#' Quantiles
#' @param x a vector or matrix
#' @param probs quantiles to compute
#' @param round numeric
#' @noRd
cyto_stat_quantile <- function(x,
                               probs = 0.5,
                               round = 2,
                               ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      round(
        quantile_cpp(
          x,
          probs
        ), 
        round
      )
    )
  # MATRIX
  } else {
    return(
      round(
        col_quantile_cpp(
          x,
          probs
        ),
        round
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
            splinefun(
              x$x, 
              x$y, 
              method = method
            ),
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
          )$value, 
          round
        )
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
    future_apply(
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
          range_cpp(x),
          round
        )
      )
    )
  # MATRIX
  } else {
    return(
      round(
        col_range_cpp(x),
        round
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
          rng <- range_cpp(x)
        } else {
          rng <- limits[, 1]
        }
        bandwidth <- diff(rng) / bins
      }
    }
    # SMOOTH - INACCURATE COUNTS OTHERWISE
    if(smooth < 1) {
      warning("'smooth' should be greater than or equal to 1!")
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
    future_apply(
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
    future_apply(
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
                            probs = c(0.01, 0.99),
                            ...) {
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      scale_cpp(
        x,
        type = type,
        probs = probs
      )
    )
    # MATRIX
  } else {
    return(
      col_scale_cpp(
        x,
        type = type,
        probs = probs
      )
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
  
  # JITTER SAMPLE-ID BARCODES
  ind <- grep("^Sample\\-ID$", colnames(x))
  if(length(ind) == 1) {
    # SET SEED FOR REPRODUCIBLE SAMPLING - REQUIRED
    set.seed(42)
    x[, ind] <- LAPPLY(
      unique(x[, ind]),
      function(w) {
        rnorm(
          n = length(
            x[x[, ind] == w, ind]
          ),
          mean = w,
          sd = 0.1
        )
      }
    )
    # RESET SEED
    rm(list=".Random.seed", envir=globalenv())
  }
  
  # LIMITS - MATRIX
  if(!is.null(dim(limits))) {
    limits <- lapply(
      1:ncol(limits), 
      function(z){
        limits[, z, drop = TRUE]
      }
    )
  }
  
  # DEFAULT LIMITS - MATCH KERNSMOOTH
  limits <- lapply(
    seq_along(limits),
    function(z){
      if(any(is.na((limits[[z]])))) {
        return(
          c(
            min(x[, z]) - 1.5 * h[z],
            max(x[, z]) + 1.5 * h[z]
          )
        )
      }
      return(limits[[z]])
    }
  )
  
  # COMPUTE BANDWIDTH USING PLUGIN METHOD
  bandwidth <- LAPPLY(
    seq_along(bandwidth), 
    function(z){
      if(.all_na(bandwidth[z])) {
        if(length(x[, z]) < 2) {
          return(NA)
        } else {
          bw <- tryCatch(
            suppressWarnings(
              dpik(
                x[, z],
                gridsize = M[z],
                range.x = limits[[z]],
                truncate = TRUE
              )
            ),
            error = function(e) {
              return(diff(limits[[z]])/M[z])
            }
          )
          # REDUCE BANDWIDTH FOR SAMPLE-ID
          if(grepl("^Sample\\-ID$", colnames(x)[z])) {
            bw <- bw/(8*length(unique(x[, z])))
          }
          return(
            bw
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
    }
  )
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
    lapply(
      seq_len(2), 
      function(z){
        L[z] <<- min(floor(tau*h[z]*(M[z]-1)/(b[z]-a[z])), M[z] - 1L)
        lvecid <- seq(0, L[z])
        facid <- (b[z] - a[z])/(h[z]*(M[z]-1L))
        w <- matrix(dnorm(lvecid*facid)/h[z])
        tot <- sum(c(w, rev(w[-1L]))) * facid * h[z]
        kapid[[z]] <<- w/tot
      }
    )
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
        bins = list(
          "x" = xpts,
          "y" = ypts
        ),
        data = x
      )
    )
  # COUNTS ONLY
  } else {
    return(
      list(
        counts = cnts,
        bkde = NA,
        bins = list(
          "x" = xpts,
          "y" = ypts
        ),
        data = x
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
  
  # VECTOR
  if(is.null(dim(x))) {
    return(
      rescale_cpp(
        x,
        scale = scale,
        limits = limits
      )
    )
    # MATRIX
  } else {
    return(
      col_rescale_cpp(
        x,
        scale = scale,
        limits = limits
      )
    )
  }
  
}

#' Compute hex bins for cyto_plot
#'
#' @param x xy matrix to compute hexbins.
#' @param limits truncate data prior to binning.
#' @param bins number of bins to use for hex binning, set to 256 by default.
#' @param smooth logical indicating whether the hexbins counts should be
#'   smoothed, set to TRUE by default.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
cyto_stat_hex <- function(x,
                          limits = list(NA, NA),
                          bins = 256,
                          smooth = TRUE) {
  
  # TODO: support hexbin count smoothing
  
  # handle zero event data
  if(nrow(x) == 0) {
    return(
      list(
        x = NA,
        y = NA,
        cell = NA,
        count = NA,
        width = NA,
        height = NA,
        density = NA,
        coords = NA
      )
    )
  }
  
  # LIMITS - MATRIX
  if(!is.null(dim(limits))) {
    limits <- lapply(
      1:ncol(limits), 
      function(z){
        limits[, z, drop = TRUE]
      }
    )
  }
  
  # XLIM
  if(.all_na(limits[[1]])) {
    limits[[1]] <- range(x[, 1])
  }
  
  # YLIM
  if(.all_na(limits[[2]])) {
    limits[[2]] <- range(x[, 2])
  }
  
  # TRIM DATA TO LIMITS
  x <- x[
    (x[, 1] >= min(limits[[1]]) &
      x[, 1] <= max(limits[[1]])) &
    (x[, 2] >= min(limits[[2]]) &
       x[, 2] <= max(limits[[2]])),
  ]
  
  # JITTER SAMPLE-ID BARCODES
  ind <- grep("^Sample\\-ID$", colnames(x))
  # POINT FOR HEXBINS ALREADY JITTERED ABOVE
  if(length(ind) == 1) {
    # SET SEED FOR REPRODUCIBLE SAMPLING - REQUIRED
    set.seed(42)
    x[, ind] <- LAPPLY(
      unique(x[, ind]),
      function(w) {
        rnorm(
          n = length(
            x[x[, ind] == w, ind]
          ),
          mean = w,
          sd = 0.1
        )
      }
    )
    # RESET SEED
    rm(list=".Random.seed", envir=globalenv())
  }
  
  # require hexbin package
  cyto_require(
    "hexbin"
  )
  
  # rounding
  round_any <- function(y, accuracy, f = round) {
    f(y/accuracy) * accuracy
  }
  
  # compute binwidth
  bw <- c(
    diff(limits[[1]])/bins,
    diff(limits[[2]])/bins
  )
  
  # x hexbin bounds
  xb <- c(
    round_any(min(limits[[1]]), bw[1], floor) - 1e-6,
    round_any(max(limits[[1]]), bw[1], ceiling) + 1e-6
  )
  xbins <- diff(xb)/bw[1]
  
  # y hexbin bounds
  yb <- c(
    round_any(min(limits[[2]]), bw[2], floor) - 1e-6,
    round_any(max(limits[[2]]), bw[2], ceiling) + 1e-6
  )
  ybins <- diff(yb)/bw[2]
  
  # call hexbin
  hex <- cyto_func_call(
    "hexbin::hexbin",
    list(
      x[, 1],
      xbnds = xb,
      xbins = bins,
      x[, 2],
      ybnds = yb,
      shape = 1,
      IDs = TRUE
    )
  )
  
  # NOTE: HEXBIN PACKAGE SMOOTHING IS SLOWER AND REDUCES RESOLUTION
  # smooth hexbin counts
  if(smooth) {
    hex <- .hexbin_smooth(
      hex,
      wts = c(1,5)
    )
    # hex <- cyto_func_call(
    #   "hexbin::smooth.hexbin",
    #   list(
    #     hex,
    #     wts = c(20, 5, 1)
    #   )
    # )
  }
  
  # compute hexagon widths
  sx <- hex@xbins / diff(hex@xbnds)
  sy <- (hex@xbins * hex@shape)/ diff(hex@ybnds)
  dx <- 1/(2 * sx)
  dy <- 1/(2* sqrt(3) * sy)
  
  # get hexbin stats
  hb <- cyto_func_call(
    "hexbin::hcell2xy",
    list(hex)
  )
  hb$cell <- hex@cell
  hb$id <- hex@cID
  hb$counts <- hex@count
  hb$width <- bw[1]
  hb$height <- bw[2]
  hb$density <- hb$count/sum(hb$count)
  hb$coords = cyto_func_call(
    "hexbin::hexcoords",
    list(
      dx,
      dy,
      1,
      sep = NA
    )
  )
  
  # return computed hexbin
  return(hb)
  
}

#' Fast Hexbin Smoothing (Resolution Preserving)
#'
#' This function smoothes hexbin counts using a weighted average of neighbors.
#' Unlike the standard smooth.hexbin(), this function:
#' 1. Preserves the exact grid resolution and plot dimensions.
#' 2. Is optimized for speed (pure R, faster than the original Fortran).
#'
#' @param bin A hexbin object.
#' @param wts Numeric vector of length 2: c(center_weight, neighbor_weight).
#' @return A smoothed hexbin object.
#' Fast & Resolution-Preserving Hexbin Smoothing
#'
#' @param bin A hexbin object.
#' @param wts Numeric vector of length 2: c(center_weight, neighbor_weight).
#' @param normalize Logical. If TRUE, scales counts back to original magnitude to prevent "fat" hexagons.
.hexbin_smooth <- function(bin, wts = c(48, 4), normalize = TRUE) {
  if (!inherits(bin, "hexbin")) stop("Input must be a hexbin object")
  
  # 1. SETUP (Transposed Layout for Speed)
  dims <- bin@dimen
  nrow <- dims[1]
  ncol <- dims[2]
  
  # Padded Matrix
  mat <- matrix(0L, nrow = ncol + 2L, ncol = nrow + 2L)
  
  # 2. CALCULATE INDICES (Store these for retrieval later)
  # We map the original 1-based cells to our padded, transposed matrix
  cell <- bin@cell - 1L
  c_idx <- cell %% ncol + 2L 
  r_idx <- cell %/% ncol + 2L
  
  # Fill the matrix at the original locations
  fill_locs <- cbind(c_idx, r_idx)
  mat[fill_locs] <- bin@count
  
  # 3. ACCUMULATE (The Smoothing Step)
  w1 <- as.integer(wts[1])
  w2 <- as.integer(wts[2])
  
  ix_c <- 2L:(ncol + 1L)
  ix_r <- 2L:(nrow + 1L)
  
  accum <- mat[ix_c, ix_r] * w1
  
  if (w2 > 0L) {
    mat_w2 <- mat * w2
    
    # Horizontal Neighbors
    accum <- accum + mat_w2[1L:ncol, ix_r] + mat_w2[3L:(ncol + 2L), ix_r]
    
    # Vertical Neighbors
    src_top <- mat_w2[, 1L:nrow]
    src_bot <- mat_w2[, 3L:(nrow + 2L)]
    is_odd_col <- (seq_len(nrow) %% 2L == 1L)
    is_even_col <- !is_odd_col
    
    # Top Neighbors
    accum[, is_odd_col] <- accum[, is_odd_col] + 
      src_top[1L:ncol, is_odd_col] + src_top[2L:(ncol+1L), is_odd_col]
    accum[, is_even_col] <- accum[, is_even_col] + 
      src_top[2L:(ncol+1L), is_even_col] + src_top[3L:(ncol+2L), is_even_col]
    
    # Bottom Neighbors
    accum[, is_odd_col] <- accum[, is_odd_col] + 
      src_bot[1L:ncol, is_odd_col] + src_bot[2L:(ncol+1L), is_odd_col]
    accum[, is_even_col] <- accum[, is_even_col] + 
      src_bot[2L:(ncol+1L), is_even_col] + src_bot[3L:(ncol+2L), is_even_col]
  }
  
  # 4. EXTRACT (The Fix)
  # Instead of finding NEW non-zero cells, we sample 'accum' ONLY at the 
  # original fill locations.
  # Note: 'accum' is unpadded (ncol x nrow), so we shift indices back by -1L
  extract_locs <- cbind(c_idx - 1L, r_idx - 1L)
  new_counts <- accum[extract_locs]
  
  # 5. NORMALIZE & UPDATE
  if (normalize) {
    total_weight <- w1 + (6 * w2)
    new_counts <- as.integer((new_counts + total_weight/2) / total_weight)
    # Ensure original events remain visible (don't round down to 0)
    new_counts[new_counts < 1L] <- 1L
  }
  
  # Update ONLY the counts. 
  # bin@cell, bin@cID, and bin@xbins remain completely strictly untouched.
  bin@count <- new_counts
  
  return(bin)
}