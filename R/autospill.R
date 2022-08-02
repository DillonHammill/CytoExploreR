## AUTOSPILL -------------------------------------------------------------------

#' Simplified version of autospill for use within CytoExploreR
#'
#' All functions below prefixed with \code{.cyto_asp_} are modified versions of
#' functions that exist in the autospill (\code{asp}) package.
#'
#' @param x object of class cytoset or GatingSet containing compensation single
#'   colour controls.
#' @param channels names of the channels for which the spillover coefficients
#'   should be computed.
#'
#' @references Roca et al (2021), AutoSpill is a principled framework that
#'   simplifies the analysis of multichromatic flow cytometry data, Nature
#'   Communications 12:2890
#'
#' @noRd
.cyto_asp_spill <- function(x,
                            channels = NULL,
                            trans = NA,
                            iter = 100,
                            trim = 0.001) {
  
  # AUTOSPILL - COMPUTE SPILLOVER COEFFICIENTS ---------------------------------
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    # DEFAULT - ALL AREA PARAMETERS
    channels <- channels[grepl("-A$", channels, ignore.case = TRUE)]
  } else {
    channels <- unique(cyto_channels_extract(x, channels))
  }
  
  # INITIAL SPILLOVER COEFFICIENTS - LINEAR DATA REQUIRED
  spill <- .cyto_asp_spill_init(
    x,
    channels = channels,
    trim = trim
  )
  
  # REFINE SPILLOVER COEFFICIENTS
  spill <- .cyto_asp_spill_refine(
    x,
    spill,
    channels = channels,
    trans = trans,
    iter = iter,
    trim = trim
  )
  
  # RETURN SPILLOVER MATRIX
  return(spill)
  
}

#' Fit robust linear model
#'
#' The following code a modified version of the \code{fit_robust_linear_model}
#' function in the autospill package.
#'
#' @param x object of class cytoset containing a single cytoframe. The channel
#'   associated with this cytoframe must be annotated in experiment details.
#' @param x_chan channel for x axis, defaults to channel stored in experiment
#'   details.
#' @param y_chan channel for y axis, must be supplied manually.
#' @param trim proportion of events to exclude from the top and bottom prior to
#'   fitting robust linear models, set to \code{0.01} by default.
#'
#' @return matrix with coefficients and p values as separate columns, rows are
#'   values for intercept and slope respectively.
#'
#' @importFrom MASS rlm
#' @importFrom stats pt lm
#'
#' @noRd
.cyto_asp_rlm <- function(x,
                          x_chan = NULL,
                          y_chan = NULL,
                          trim = 0.001) {
  
  # X CHANNEL
  if(is.null(x_chan)) {
    x_chan <- cyto_details(x)$channel
  }
  # EXTRACT DATA
  exprs <- cyto_exprs(
    x[[1]],
    channels = c(
      x_chan, 
      y_chan
    )
  )
  # COMPUTE QUANTILES & REMOVE EXTREME VALUES
  if(trim != 0) {
    exprs <- exprs[
      do.call(
        "intersect",
        lapply(
          seq_len(ncol(exprs)), 
          function(p) {
            v <- exprs[, p]
            q <- quantile(
              v,
              probs = c(trim, 1 - trim)
            )
            return(
              which(v > q[1] & v < q[2])
            )
          }
        )
      )
      ,]
  }
  # FIT ROBUST LINEAR MODEL
  rlm <- suppressWarnings(
    rlm(
      exprs[, 2] ~ exprs[, 1]
    )
  )
  # RLM CONVERGENCE
  if(rlm$converged) {
    rlm_coef <- rlm$coefficients
    rlm_t <- summary(rlm)$coefficients[, 3]
    rlm_df <- summary(rlm)$df[2]
    rlm_pval <- 2*(pt(abs(rlm_t), rlm_df, lower.tail = FALSE))
  # RLM NOT CONVERGED - USE LM INSTEAD
  } else {
    rlm <- lm(exprs[, 2] ~ exprs[, 1])
    rlm_coef <- rlm$coefficients
    rlm_pval <- summary(rlm)$coefficients[, 4]
  }
  # # PLOT DATA
  # plot(
  #   exprs[, 1:2],
  #   xlab = x_chan,
  #   ylab = y_chan,
  #   pch = 16,
  #   cex = 0.8
  #   # xlim = c(0, 262144),
  #   # ylim = c(0, 262144)
  # )
  # # PLOT RLM MODEL
  # abline(
  #   rlm,
  #   col = "red",
  #   lwd = 2
  # )
  # RETURN RLM PARAMETERS
  return(
    cbind(
      rlm_coef,
      rlm_pval
    )
  )

}

#' Compute initial spillover coefficients
#'
#' The following code is a modified version of the \code{get_marker_spillover}
#' function in the autospill package.
#'
#' @param x object of class cytoset containing single colour compensation
#'   controls with any unstained controls removed.
#' @param channels names of the channels for which the spillover coefficients
#'   should be computed.
#'   
#' @return list of matrices with regression intercepts and coefficients.  
#'
#' @noRd
.cyto_asp_spill_init <- function(x,
                                 channels,
                                 trim = 0.001) {
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # NUMBER OF CHANNELS
  n <- length(channels)
  
  # INITIALISE COEFFICIENTS
  spill_zero <- rep(0, n)
  names(spill_zero) <- channels
  
  # # CREATE PLOT LAYOUT
  # par(mfrow = c(3,4),
  #     mar = c(4.1,4.1,2.1,2.1))
  
  # FIT RLM & EXTRACT INTERCEPTS & COEFFICIENTS
  spill_coef <- do.call(
    "rbind",
    structure(
      lapply(
        pd$channel, 
        function(z) {
          # RESET INTERCEPTS & COEFFICIENTS
          spill_int <- spill_zero
          spill_coef <- spill_zero
          # PREPARE DATA
          cs <- cyto_select(
            x,
            channel = z
          )
          # FIT RLM MODELS & UPDATE INTERCEPTS & COEFFICIENTS
          lapply(
            channels, 
            function(w) {
              # DIAGONAL PARAMETERS COMBINATIONS
              if(z == w) {
                spill_coef[w] <<- 1
                # NON-DIAGONAL PARAMETER COMBINATIONS
              } else {
                # DATA PRIMARY & SECONDARY CHANNEL
                rlm_params <- .cyto_asp_rlm(
                  cs,
                  x_chan = z,
                  y_chan = w,
                  trim = trim
                )
                # STORE COEFFICIENT & INTERCEPT
                spill_int[w] <<- rlm_params[1, 1]
                spill_coef[w] <<- rlm_params[2, 1]
              }
            }
          )
          return(
            c(spill_int, spill_coef)
          )
        }
      ),
      names = pd$channel
    )
  )
  
  # # RESET PLOT LAYOUT
  # mtext("INITIAL SPILLOVER", outer = TRUE)
  # par(mfrow = c(1,1),
  #     oma = c(0,0,0,0),
  #     mar = c(4.1,4.1,2.1,2.1))
  
  # LIST OF INTERCEPTS & COEFFICIENTS 
  return(
    list(
      int = spill_coef[, 1:n, drop = FALSE],
      coef = spill_coef[, (n+1):(2*n), drop = FALSE]
    )
  )

}

#' Refine spillover matrix coefficients
#'
#' The following code is a modified version of the \code{refine_spillover}
#' function in the autospill package.
#'
#' @param x object of class cytoset containing untransformed uncompensated
#'   single colour controls.
#' @param spill object returned by .cyto_asp_spill_init().
#' @param trans transformerList.
#' @param channels names of the channels for which spillover coefficients should
#'   be computed.
#'
#' @return spillover matrix computed and refined using the autospill approach.
#'
#' @noRd
.cyto_asp_spill_refine <- function(x,
                                   spill = NULL,
                                   trans = NA,
                                   channels,
                                   iter = 100,
                                   trim = 0.001) {
  
  # INITIALISE PARAMETERS
  rs_convergence <- FALSE
  rs_exit <- FALSE
  rs_iter <- 0
  rs_iter_last <- FALSE
  rs_iter_width <- floor(log10(100)) + 1 # asp$rs.iter.max
  rs_lambda <- 1 # aps$rs.lambda.coarse
  rs_delta <- -1
  rs_delta_threshold <- 0.01 # asp$rs.delta.threshold.untr
  rs_delta_history <- rep(-1, 10) # asp$rs.delta.history.n
  rs_scale_untransformed <- TRUE
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # NUMBER OF SAMPLES
  n <- length(channels)
  
  # INITIAL SQUARE SPILLOVER MATRIX
  spill_curr <- matrix(
    0, 
    ncol = length(channels),
    nrow = length(channels),
    dimnames = list(
      channels,
      channels
    )
  )
  diag(spill_curr) <- 1
  
  spill_update <- spill_curr
  
  # FILL UPDATE MATRIX
  cnt <- 0
  apply(
    spill$coef,
    1,
    function(z) {
      cnt <<- cnt + 1
      spill_update[
        match(rownames(spill$coef)[cnt], rownames(spill_update)),
      ] <<- z[match(colnames(spill$coef), colnames(spill_update))]
    }
  )
  diag(spill_update) <- 0
  
  # REFINE SPILLOVER MATRIX
  while(!rs_exit) {
    # UPDATE SPILLOVER MATRIX
    spill_curr <- spill_curr + spill_update
    spill_curr <- sweep(
      spill_curr,
      1, 
      diag(spill_curr),
      "/"
    )
    # COMPENSATION ERROR
    spill_error <- .cyto_asp_spill_error(
      x,
      spill = spill_curr,
      trans = trans,
      transform = !rs_scale_untransformed,
      channels = channels,
      trim = trim
    )
    # SLOPE ERROR -> SQUARE
    spill_slope_error <- matrix(
      0,
      ncol = length(channels),
      nrow = length(channels),
      dimnames = list(
        channels,
        channels
      )
    )
    diag(spill_slope_error) <- 1
    # FILL SQUARE SLOPE ERROR MATRIX
    cnt <- 0
    apply(
      spill_error$slope,
      1,
      function(z) {
        cnt <<- cnt + 1
        spill_slope_error[
          match(rownames(spill_error$slope)[cnt], rownames(spill_slope_error)),
        ] <<- z[match(colnames(spill_error$slope), colnames(spill_slope_error))]
      }
    )
    spill_slope_error <- spill_slope_error - diag(ncol(spill_slope_error))
    # UPDATE DELTA PARAMETERS
    rs_delta_prev <- rs_delta
    # DROP EMPTY ROWS HERE
    rs_delta <- sd(
      spill_slope_error
    )
    rs_delta_max <- max(
      abs(
        spill_slope_error
      )
    )
    # UPDATE DELTA HISTORY
    if(rs_delta_prev >= 0) {
      # asp$rs.delta.history.n = 10
      rs_delta_history[rs_iter %% 10 + 1] <- rs_delta - rs_delta_prev
    } else {
      rs_delta_history[rs_iter %% 10 + 1] <- -1
    }
    # CHANGE IN DELTA
    rs_delta_change <- mean(rs_delta_history)
    # SWITCH TO BIEXPONENTIAL SCALE
    if(rs_scale_untransformed && rs_delta_max < rs_delta_threshold) {
      # RESET LAMBDA & DELTA HISTORY
      rs_scale_untransformed <- FALSE
      rs_delta_threshold <- 1e-04 # asp$rs.delta.threshold.tran
      rs_lambda <- 1
      rs_delta <- -1
      rs_delta_history <- rep(-1, 10)
      rs_delta_change <- -1
    }
    # REDUCE LAMBDA & RESET DELTA HISTORY
    if(rs_delta_change > -1e-06 && rs_lambda == 1) {
      rs_lambda <- 0.1
      rs_delta <- -1
      rs_delta_history <- rep(-1, 10) # asp$rs.delta.history.n
      rs_delta_change <- -1
    }
    # CONVERGENCE
    rs_convergence <- !rs_scale_untransformed && 
      (rs_delta_max < rs_delta_threshold || rs_delta_change > -1e-06)
    # EXIT
    rs_exit <- (rs_convergence && rs_iter_last) ||
      (rs_delta_change > -1e-06 && rs_scale_untransformed) ||
      (!rs_convergence && rs_iter == iter) || # max iterations = 100
      rs_iter > iter
    rs_iter_last <- rs_convergence
    rs_iter <- rs_iter + 1
    # UPDATE SPILLOVER MATRIX
    spill_update <- rs_lambda * (spill_slope_error %*% spill_curr)
  }
  
  # # CHECK CONVERGENCE
  # if(!rs_convergence) {
  #   stop(
  #     "Autospill failed to converge when computing spillover matrix."
  #   )
  # }
  
  # RETURN REFINED SPILLOVER MATRIX
  return(spill_curr)
  
}

#' Compute compensation error
#'
#' The following code is a modified version of the \code{get_compensation_error}
#' function in the autospill package.
#'
#' @param x object of class cytoset containing uncompensated single colour
#'   controls.
#' @param spill matrix containing the spillover coefficient estimates to be
#'   applied to the data.
#' @param trans object of class transformerList containing the transformation
#'   definitions for fluorescent channels to apply to the data when
#'   \code{transform = TRUE}.
#' @param transform logical indicating whether biexponential transformations
#'   should be applied prior to computing spillover error, set to FALSE by
#'   default,
#' @param channels the names of the channels for which spillover coefficients
#'   are to be computed.
#'
#' @return list of matrices describing compensation error, with intercepts,
#'   coefficients, slopes, and skewness.
#'
#' @noRd
.cyto_asp_spill_error <- function(x,
                                  spill = NULL,
                                  trans = NA,
                                  transform = FALSE,
                                  channels,
                                  trim = 0.001) {
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # NUMBER OF CONTROLS
  n <- length(channels)
  
  # COMPENSATE DATA
  x <- cyto_compensate(
    x,
    spillover = spill,
    copy = TRUE
  )
  
  # TRANSFORM DATA
  if(transform) {
    # TRANSFORMERS
    if(.all_na(trans)) {
      trans <- cyto_transformers_define(
        x,
        channels = channels,
        type = "biex",
        widthBasis = -100,
        plot = FALSE,
        progress = FALSE
      )
    }
    # APPLY TRANSFORMATIONS
    x <- cyto_transform(
      x,
      trans = trans,
      plot = FALSE,
      quiet = TRUE
    )
  }
  
  # INITIALISE COEFFICIENTS
  spill_zero <- rep(0, n)
  names(spill_zero) <- channels
  
  # # CREATE PLOT LAYOUT
  # par(mfrow = c(3,4),
  #     mar = c(4.1,4.1,2.1,2.1))
  
  # FIT RLM - EXTRACT INTERCEPTS, COEFFICIENTS, SLOPES & SKEWNESS
  spill_params <- do.call(
    "rbind",
    structure(
      lapply(
        pd$channel, 
        function(z){
          # RESET PARAMETERS
          spill_int <- spill_zero
          spill_coef <- spill_zero
          spill_slope <- spill_zero
          spill_skew <- spill_zero
          # EXTRACT DATA
          cs <- cyto_select(
            x,
            channel = z
          )
          # FIT RLM MODELS & UPDATE PARAMETERS
          lapply(
            channels,
            function(w) {
              # DIAGONAL PARAMETERS COMBINATIONS
              if(z %in% w) {
                spill_coef[w] <<- 1
                spill_slope[w] <<- 1
                # NON-DIAGONAL PARAMETER COMBINATIONS
              } else {
                # DATA PRIMARY & SECONDARY CHANNEL
                rlm_params <- .cyto_asp_rlm(
                  cs,
                  x_chan = z,
                  y_chan = w,
                  trim = trim
                )
                # STORE COEFFICIENT & INTERCEPT
                spill_int[w] <<- rlm_params[1, 1]
                spill_coef[w] <<- rlm_params[2, 1]
                # SLOPE/SKEWNESS - LINEAR DATA
                if(!transform) {
                  spill_slope[w] <<- spill_coef[w]
                  # SLOPE/SKEWNESS - TRANSFORMED DATA
                } else {
                  yp <- cyto_apply(
                    cs, 
                    "quantile",
                    probs = c(0.01, 0.99),
                    input = "matrix",
                    channels = z,
                    copy = FALSE
                  )[, 1]
                  xp <- c(spill_int[w] + spill_coef[w] * yp[1],
                          spill_int[w] + spill_coef[w] * yp[2])
                  if(yp[1] == yp[2] || xp[1] == xp[2]) {
                    spill_slope[w] <<- 0
                  } else {
                    ypt <- .cyto_transform(
                      yp,
                      trans = trans,
                      channel = z,
                      inverse = TRUE
                    )
                    xpt <- .cyto_transform(
                      xp,
                      trans = trans,
                      channel = w,
                      inverse = TRUE
                    )
                    spill_slope[w] <<- spill_coef[w] * (xpt[2] - xpt[1]) * 
                      (yp[2] - yp[1]) / ((xp[2] - xp[1])*(ypt[2] - ypt[1]))
                  }
                }
              }
            }
          )
          # RETURN RLM PARAMETERS
          return(
            c(
              spill_int,
              spill_coef,
              spill_slope
            )
          )
        }
      ),
      names = pd$channel
    )
  )
    
  # # RESET PLOT LAYOUT
  # mtext("SPILLOVER REFINE", outer = TRUE)
  # par(mfrow = c(1,1),
  #     oma = c(0,0,3,0),
  #     mar = c(4.1,4.1,2.1,2.1))
  
  # SPILLOVER ERROR PARAMETERS
  return(
    list(
      int = spill_params[, 1:n, drop = FALSE],
      coef = spill_params[, (n+1):(2*n), drop = FALSE],
      slope = spill_params[, (2*n+1):(3*n), drop = FALSE]
    )
  )
  
}
