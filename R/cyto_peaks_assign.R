# CYTO_PEAKS_ASSIGN ------------------------------------------------------------

#' Assign peak detectors to pre-gated single colour controls
#'
#' @param x object of class \code{cytoset}, \code{GatingHierarchy} or
#'   \code{GatingSet}.
#' @param events number of events to extract from unstained and stained controls
#'   to identify peak detectors, set to 500 events by default.
#' @param overwrite logical indicating whether existing detector assignments in
#'   \code{cyto_details(x)$channel} be overwritten, set to FALSE by default.
#'   2param plot logical indicating whether the estimated spectrum for each
#'   control should be displayed for visual inspection, set to FALSE by default.
#' @param plot logical indicating whether the results should be plotted for
#'   troubleshooting, set to FALSE by default.
#'
#' @return updated \code{"channel"} column in \code{cyto_details(x)}.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_channel_match}}
#'
#' @export
cyto_peaks_assign <- function(x,
                              type = "nmf",
                              events = 500,
                              overwrite = FALSE,
                              plot = FALSE) {
  
  # RSVD PACKAGE REQUIRED
  cyto_require(
    "rsvd",
    source = "CRAN",
    ref = paste0(
      "Erichson NB, Voronin S, Brunton SL, Kutz JN (2019). Randomized ",
      "Matrix Decompositions Using R. Journal of Statistical Software, ",
      "89(11), 1–48. doi: 10.18637/jss.v089.i11."
    )
  )
  
  # EXPERIMENT DETAILS
  pd <- pData(x)
  
  # GROUP VARIABLE REQUIRED
  if(!"group" %in% colnames(pd)) {
    message(
      "'group' is missing from cyto_details(x) - grouping all samples!"
    )
    pData(x)$group <- "control"
  }
  
  # ADD CHANNEL VARIABLE
  if(!"channel" %in% colnames(pd)) {
    message(
      "Adding 'channel' variable to cyto_details(x) to store peak detectors."
    )
    pData(x)$channel <- NA
  }
  
  # PARENT VARIABLE REQUIRED
  if(!"parent" %in% colnames(pd)) {
    message(
      "Adding 'parent' variable to cyto_details(x) to select parent populations."
    )
    pData(x)$parent <- "root"
  }
  
  # SAMPLE GROUPING
  x_list <- cyto_group_by(
    x,
    group_by = "group"
  )
  
  # PROGRESS BAR
  pb <- cyto_progress(
    label = "cyto_peaks_assign()",
    total = sum(
      !(pd$channel %in% "Unstained" | 
          grepl("Unstained|NIL", cyto_names(x), ignore.case = TRUE))
    ),
    clear = FALSE
  )
  
  # LOOP THROUGH GROUPS
  peaks <- lapply(
    seq_along(x_list),
    function(id) {
      # GATINGSET
      gs <- x_list[[id]]
      # GROUP DETAILS
      pd <- pData(gs)
      # AUTOFLUORESCENCE SPECTRUM
      auto <- NULL
      # GROUP UNSTAINED CONTROL
      unst_idx <- which(
        pd$channel %in% "Unstained" | 
        grepl("Unstained|NIL", cyto_names(gs), ignore.case = TRUE)
      )
      # STAINED CONTROLS
      idx <- seq_along(gs)
      if(length(unst_idx) > 0) {
        idx <- idx[-unst_idx]
      }
      # LOCATE PEAK DETECTORS
      unlist(
        lapply(
          idx,
          function(w) {
            # CURRENT PEAK
            peak <- pData(gs)$channel[w]
            # NOTE: DONT CHECK AGAINST CHANNELS FOR UNMIXED CONTROLS
            # RETAIN CURRENT PEAK
            if(!peak %in% c(NA, "NA", "") & !overwrite) {
              peak <- structure(
                list(peak),
                names = cyto_names(gs)[w]
              )
            # LOCATE PEAK
            } else {
              # AUTOFLUORESCENCE SUBTRACTION
              auto <- NULL
              if(length(unst_idx) > 0) {
                # ESTIMATE SPECTRUM FOR UNSTAINED
                exprs <- cyto_data_extract(
                  gs[unst_idx],
                  parent = pd$parent[w],
                  channels = cyto_fluor_channels(gs),
                  format = "matrix",
                  inverse = TRUE,
                  copy = TRUE,
                  coerce = TRUE,
                  events = events
                )[[1]][[1]]
                # DROP NEGATIVE VALUES
                exprs[exprs < 0] <- 0
                # NMF
                nmf <- RcppML::nmf(
                  exprs,
                  k = 1
                )
                auto <- nmf$h
                # auto <- t(
                #   rsvd::rsvd(
                #     exprs,
                #     k = 1
                #   )$v
                # )
                # rownames(H) <- 1
                # colnames(H) <- colnames(exprs)
                # AUTOFLUORESCENCE SPECTRUM
                
                rownames(auto) <- "auto"
                colnames(auto) <- colnames(exprs)
              }
              # EXTRACT & COERCE UNSTAINED & STAINED DATA
              exprs <- cyto_data_extract(
                gs[c(unst_idx, w)],
                parent = pd$parent[w],
                channels = cyto_fluor_channels(gs),
                format = "matrix",
                inverse = TRUE,
                copy = TRUE,
                coerce = TRUE,
                events = events
              )[[1]][[1]]
              # TODO: ADD WARNING FOR TOO FEW EVENTS
              # COMPUTE SVD 
              if(nrow(exprs) > 0) {
                # rsvd <- rsvd::rsvd(
                #   exprs,
                #   k = 5
                # )
                # res <- t(abs(rsvd$v))
                # colnames(res) <- colnames(exprs)
                # rownames(res) <- seq_len(nrow(res))
                # print(res)
                # cyto_plot_spectra(
                #   list(res)
                # )
                # peak <- structure(
                #   c(colnames(exprs)[which.max(abs(rsvd$v[, 1]))]),
                #   names = cyto_names(gs)[w]
                # )
                exprs[exprs < 0] <- 0
                nmf <- RcppML::nmf(
                  exprs,
                  k = 1
                )
                H <- nmf$h
                rownames(H) <- cyto_names(gs)[w]
                colnames(H) <- colnames(exprs)
                
                # H <- t(
                #   rsvd::rsvd(
                #     exprs,
                #     k = 1
                #   )$v
                # )
                # rownames(H) <- 1
                # colnames(H) <- colnames(exprs)
                # AUTOFLUORESCENCE SUBTRACTION
                if(!is.null(auto)) {
                  H[1, ] <- residuals(
                    lm(
                      as.numeric(H[1, ]) ~ as.numeric(auto[1, ]) - 1
                    )
                  )
                  H[H < 0] <- 0
                }
                peak <-  structure(
                  colnames(exprs)[
                    which.max(H[1, ])
                  ],
                  names = cyto_names(gs)[w]
                )
                if(plot) {
                  cyto_plot_spectra(
                    list(
                      H
                    )
                  )
                }
              } else {
                peak <- structure(
                  NA,
                  names = cyto_names(gs)[w]
                )
              }
            }
            # UPDATE PROGRESS BAR
            cyto_progress(pb)
            # EXTRACT PEAK DETECTOR
            return(
              peak
            )
          }
        )
      )
    }
  )
  peaks <- do.call("c", peaks)
  # UPDATE PEAK DETECTORS IN GATINGSET
  pData(x)$channel[
    match(names(peaks), cyto_names(x))
  ] <- peaks
  # RETURN UPDATED GATINGSET
  return(x)
}