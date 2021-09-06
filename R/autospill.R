## AUTOSPILL -------------------------------------------------------------------

#' Simplified version of autospill for use within CytoExploreR
#'
#' All functions below prefixed with \code{.cyto_asp_} are modified versions of
#' functions that exist in the autospill (\code{asp}) package.
#'
#' @param x object of class cytoset or GatingSet containing compensation single
#'   colour controls.
#' @param parent name of the parent population to extract from GatingSet.
#'
#' @references Roca et al (2021), AutoSpill is a principled framework that
#'   simplifies the analysis of multichromatic flow cytometry data, Nature
#'   Communications 12:2890
#'
#' @noRd
.cyto_asp_spill <- function(x,
                            parent = "root") {
  
  # HANDLE DUPLICATE CONTROLS
  # CHANNEL MATCH REQUIRED
  
  # SAMPLE PREPARATION ---------------------------------------------------------
  
  # DATA MUST BE LINEAR
  
  # EXTRACT DATA - REMOVE UNSTAINED CONTROL
  cs <- cyto_data_extract(
    x,
    parent = parent,
    format = "cytoset",
    select = list(
      channel = "Unstained",
      exclude = TRUE
    ),
    copy = FALSE
  )[[1]]
  
  # AUTOSPILL - COMPUTE SPILLOVER COEFFICIENTS ---------------------------------
  
  # INITIAL SPILLOVER COEFFICIENTS
  spill <- .cyto_asp_spill_init(cs)
  
  # REFINE SPILLOVER COEFFICIENTS
  spill <- .cyto_asp_spill_refine(
    cs,
    spill
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
#' @param skew logical indicating whether the skewness of the data should be
#'   included in the returned matrix, set to FALSE by default.
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
                          trim = 0.01,
                          skew = FALSE) {
  
  # X CHANNEL
  if(is.null(x_chan)) {
    x_chan <- cyto_details(x)$channel
  }
  
  # EXTRACT DATA
  exprs <- cyto_exprs(
    x[[1]],
    channels = c(x_chan, y_chan)
  )
  # COMPUTE QUANTILES & REMOVE EXTREME VALUES
  exprs <- exprs[
    do.call(
      "intersect",
      lapply(seq_len(ncol(exprs)), function(p){
        v <- exprs[, p]
        q <- quantile(
          v,
          probs = c(trim, 1 - trim)
        )
        return(
          which(v > q[1] & v < q[2])
        )
      })
    )
    ,]
  # FIT ROBUST LINEAR MODEL
  rlm <- rlm(exprs[, 2] ~ exprs[, 1])
  # RLM CONVERGENCE
  if(rlm$converged) {
    rlm_coef <- rlm$coefficients
    rlm_t <- summary(rlm)$coefficients[, 3]
    rlm_df <- summary(rlm)$df[2]
    rlm_pval <- 2*(pt(abs(rlm_t), rlm_df, lower.tail = FALSE))
    # RLM NOT CONVERGED - USE LM INSTEAD
  } else {
    warning(
      paste0(
        "fitting rlm() to ",
        y_chan, " ~ ", x_chan,
        " did not converge - resorting to lm() instead."
      )
    )
    rlm <- lm(exprs[, 2] ~ exprs[, 1])
    rlm_coef <- rlm$coefficients
    rlm_pval <- summary(rlm)$coefficients[, 4]
  }
  # SKEWNESS
  rlm_skew <- c(0,0)
  if(skew) {
    rlm_skew <- apply(
      exprs,
      2,
      function(z){
        n <- length(z)
        (sum((z-mean(z))^3)/n)/(sum((z-mean(z))^2)/n)^(3/2)
      }
    )
  }
  # RETURN RLM PARAMETERS
  return(
    cbind(
      rlm_coef,
      rlm_pval,
      rlm_skew
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
#'   
#' @return list of matrices with regression intercepts and coefficients.  
#'
#' @noRd
.cyto_asp_spill_init <- function(x) {
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # NUMBER OF CONTROLS
  n <- length(x)
  
  # INITIALISE COEFFICIENTS
  spill_zero <- rep(0, n)
  names(spill_zero) <- pd$channel
  
  # FIT RLM & EXTRACT INTERCEPTS & COEFFICIENTS
  spill_coef <- do.call(
    "rbind",
    structure(
      lapply(pd$channel, function(z) {
        # RESET INTERCEPTS & COEFFICIENTS
        spill_int <- spill_zero
        spill_coef <- spill_zero
        # PREPARE DATA
        cs <- cyto_select(
          x,
          channel = z
        )
        # FIT RLM MODELS & UPDATE INTERCEPTS & COEFFICIENTS
        lapply(pd$channel, function(w){
          # DIAGONAL PARAMETERS COMBINATIONS
          if(z == w) {
            spill_coef[z] <<- 1
            # NON-DIAGONAL PARAMETER COMBINATIONS
          } else {
            # DATA PRIMARY & SECONDARY CHANNEL
            rlm_params <- .cyto_asp_rlm(
              cs,
              x_chan = z,
              y_chan = w,
              trim = 0.01
            )
            # STORE COEFFICIENT & INTERCEPT
            spill_int[w] <<- rlm_params[1, 1]
            spill_coef[w] <<- rlm_params[2, 1]
          }
        })
        return(
          c(spill_int, spill_coef)
        )
      }),
      names = pd$channel
    )
  )
  
  # LIST OF INTERCEPTS & COEFFICIENTS 
  return(
    list(
      int = spill_coef[, 1:n],
      coef = spill_coef[, tail(1:ncol(spill_coef), n)]
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
#'
#' @return spillover matrix computed and refined using the autospill approach.
#'
#' @noRd
.cyto_asp_spill_refine <- function(x,
                                   spill = NULL) {
  
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
  n <- length(x)
  
  # INITIAL SPILLOVER VALUES
  spill_curr <- diag(n)
  colnames(spill_curr) <- pd$channel
  rownames(spill_curr) <- pd$channel
  spill_update <- spill$coef - diag(n)
  
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
    # SPILLOVER ORIGINAL
    spill_curr_original <- spill_curr
    # COMPENSATION ERROR
    spill_error <- .cyto_asp_spill_error(
      x,
      spill = spill_curr_original
    )
    # SLOPE ERROR
    spill_slope_error <- spill_error$slope - diag(n)
    # UPDATE DELTA PARAMETERS
    rs_delta_prev <- rs_delta
    rs_delta <- sd(spill_slope_error)
    rs_delta_max <- max(abs(spill_slope_error))
    # UPDATE DELTA HISTORY
    if(rs_delta_prev >= 0) {
      # asp$rs.delta.history.n = 10
      rs_delta_history[rs_iter %% 10 + 1] <- rs_delta - rs_delta_prev
    } else {
      rs_delta_history[rs_iter %% 10 + 1] <- -1
    }
    # CHANGE IN DELTA
    rs_delta_change <- mean(rs_delta_history)
    # # UPDATE PARAMETERS
    # if(rs_scale_untransformed && rs_delta_max < rs_delta_threshold) {
    #   # SWITCH TO BIEXPONENTIAL SCALE & RESET DELTA PARAMETERS
    #   rs_scale_untransformed <- FALSE
    #   rs_delta_threshold <- 1e-04 # asp$rs.delta.threshold.tran
    #   rs_lambda <- 1
    #   rs_delta <- -1
    #   rs_delta_history <- rep(-1, 10) # asp$rs.delta.history.n
    #   rs_delta_change <- -1
    # }
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
      (!rs_convergence && rs_iter == 100) || # max iterations = 100
      rs_iter > 100
    rs_iter_last <- rs_convergence
    rs_iter <- rs_iter + 1
    # UPDATE SPILLOVER MATRIX
    spil_update <- rs_lambda * (spill_slope_error %*% spill_curr)
  }
  
  # CHECK CONVERGENCE
  if(!rs_convergence) {
    stop(
      "Autospill failed to converge when computing spillover matrix."
    )
  }
  
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
#'
#' @return list of matrices describing compensation error, with intercepts,
#'   coefficients, slopes, and skewness.
#'
#' @noRd
.cyto_asp_spill_error <- function(x,
                                  spill = NULL) {
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # NUMBER OF CONTROLS
  n <- length(x)
  
  # COMPENSATE DATA
  x <- cyto_compensate(
    x,
    spillover = spill,
    copy = TRUE
  )
  
  # INITIALISE COEFFICIENTS
  spill_zero <- rep(0, n)
  names(spill_zero) <- pd$channel
  
  # FIT RLM - EXTRACT INTERCEPTS, COEFFICIENTS, SLOPES & SKEWNESS
  spill_params <- do.call(
    "rbind",
    structure(
      lapply(pd$channel, function(z){
        # RESET PARAMETERS
        spill_int <- spill_zero
        spill_coef <- spill_zero
        spill_slope <- spill_zero
        spill_skew <- spill_zero
        # PREPARE DATA
        cs <- cyto_select(
          x,
          channel = z
        )
        # FIT RLM MODELS & UPDATE PARAMETERS
        lapply(pd$channel, function(w){
          # DIAGONAL PARAMETERS COMBINATIONS
          if(z == w) {
            spill_coef[z] <<- 1
            spill_slope[z] <<- 1
            # NON-DIAGONAL PARAMETER COMBINATIONS
          } else {
            # DATA PRIMARY & SECONDARY CHANNEL
            rlm_params <- .cyto_asp_rlm(
              cs,
              x_chan = z,
              y_chan = w,
              trim = 0.01,
              skew = TRUE
            )
            # STORE COEFFICIENT & INTERCEPT
            spill_int[w] <<- rlm_params[1, 1]
            spill_coef[w] <<- rlm_params[2, 1]
            spill_slope[w] <<- spill_coef[w]
            spill_skew[w] <<- rlm_params[2, 3] # skewness in secondary channel
          }
        })
        # RETURN RLM PARAMETERS
        return(
          c(
            spill_int,
            spill_coef,
            spill_slope,
            spill_skew
          )
        )
      }),
      names = pd$channel
    )
  )
    
  # SPILLOVER ERROR PARAMETERS
  return(
    list(
      int = spill_params[, 1:n],
      coef = spill_params[, (n+1):(2*n)],
      slope = spill_params[, (2*n+1):(3*n)],
      skew = spill_params[, (3*n+1):(4*n)]
    )
  )
  
}
