# CYTO_PANEL_DESIGN ------------------------------------------------------------

#' Compute panel design statistics from single colour controls
#'
#' @param x object of class \code{GatingSet} or a list of \code{GatingSet}
#'   objects containing a set(s) of pre-gated single colour controls.
#' @param specra optional pre-computed spillover or unmixing matrices either as
#'   a vector of CSV file names or a list of matrix objects.
#' @param type algorithm to use to compute the spectra, if the supplied
#'   GatingSets haven't been through \code{cyto_spillover_compute()} or
#'   \code{cyto_unmix_compute()} and the results supplied through either
#'   \code{spillover} or \code{unmix}.
#' @param channels names of channels to computed the statistics over, set to the
#'   column names of \code{spectra} by default.
#' @param ... additional arguments passed to either
#'   \code{cyto_spillover_compute} or \code{cyto_unmix_compute()}.
#'
#' @return a list containing spectra, similarity, complexity, brightness,
#'   interference, spread, total_spread, hotspot and purity and panel complexity
#'   per panel.
#'
#' @seealso \code{\link{cyto_spillover_compute}}
#' @seealso \code{\link{cyto_unmix_compute}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # import single color controls
#' gs <- GatingSet(Compensation)
#' gs <- cyto_gatingTemplate_apply(gs)
#' gs <- cyto_transform(gs)
#'
#' # compute spillover matrix
#' spill <- cyto_spillover_compute(gs)
#'
#' # compute panel design statistics
#' cyto_panel_design(gs, spillover = spill)
#' }
#'
#' @export
cyto_panel_design <- function(x,
                              spectra = NULL,
                              type = "CytoDecode",
                              channels = NULL,
                              details = NULL,
                              ...) {
  
  # TODO: ADD REFERENCES
  
  # NOTE: X SHOULD BE UNMIXED FOR SPECTRAL - SPREAD ON UNMIXED SCALE
  # NOTE: GATES ARE NAMED <MARKER> DYE-/+
  
  # NOTE: FOR SPECTRAL DATA IF NOT UNMIXED - CYTO_SPILLOVER_COMPUTE() CALLED
  # TO GET GATES THEN CYTO_UNMIX_COMPUTE() CALLED TO GET UNMIXED GATES
  
  # USE CASES:
  # 1. Fixed panel and we want to compute statistics
  # 2. Different panels and we want to compute statistics
  # 3. Titrated controls and we want to compute statistics
  
  # WORKFLOW:
  # import, transform and pre-gate controls
  # estimate spillover
  # apply spillover (use original transformations)
  # compute panel design statistics
  
  # NOTE:
  # If we modify transformations gates will shift are cause issues.
  # We need inverse transform to put gates on linear scale then apply new
  # transformations
  
  
  # CHECK X - GATINGSET OR LIST (GatingSetList)
  if(!cyto_class(x, c("list", "GatingSet"), TRUE)) {
    stop(
      "'x' must be either a GatingSet, list of GatingSets!"
    )
  }
  
  # CHECK LIST CONTAINS GATINGSETS ONLY
  if(cyto_class(x, "list", TRUE)) {
    if(!all(LAPPLY(x, cyto_class, "GatingSet", TRUE))) {
      stop(
        "'x' must be a list of GatingSet objects!"
      )
    }
    if(is.null(names(x))) {
      names(x) <- paste0("panel-", 1:length(x))
    }
  }
  
  # GATINGSET -> LIST
  if(cyto_class(x, "GatingSet", TRUE)) {
    x <- list("panel-1" = x)
  }
  
  # SPILLOVER MATRICES SUPPLIED?
  spillover <- FALSE
  
  # PREPARE SPECTRA
  if(!is.null(spectra)) {
    # SPECTRA CSV FILES PER BATCH
    if(is.character(spectra)) {
      # UNMIXING MATRIX OR SPILLOVER MATRIX
      spectra <- structure(
        lapply(
          spectra,
          function(z) {
            data.matrix(
              read_from_csv(
                z,
                data.table = FALSE
              )
            )
          }
        ), names = spectra
      )
    }
    # CONVERT SPECTRA TO LIST
    if(cyto_class(spectra, c("matrix", "data.frame"), TRUE)) {
      spectra <- structure(
        rep(list(data.matrix(spectra)), length(x)),
        names = names(x)
      )
    }
    # A UNIQUE SPECTRAL MATRIX REQUIRED PER BATCH
    if(length(spectra) != length(x)) {
      stop(
        "'spectra' must be supplied for each GatingSet!"
      )
    }
    # CONVERT -> DATA.MATRIX
    spill <- c()
    spectra <- structure(
      lapply(
        spectra,
        function(z) {
          z <- data.matrix(z)
          if(nrow(z)!=ncol(z)) {
            if(all(diag(z) == 1)) {
              spill <<- c(spill, TRUE)
              if(is.null(rownames(z))) {
                rownames(z) <- colnames(z)
              }
            } else {
              if(is.null(rownames(z))) {
                stop(
                  "Matrix of spectra must have row names!"
                )
              }
              spill <<- c(spill, FALSE)
            }
          }
          return(z)
        }
      ),
      names = names(spectra)
    )
    # SPILLOVER MATRICES?
    spillover <- all(spill)
  }
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(
      x,
      exclude = "\\-[HW]$"
    )
  } else {
    channels <- unique(
      cyto_channels_extract(
        x, 
        channels
      )
    )
  }
  
  # LOOP THROUGH GATINGSETS
  stats <- structure(
    lapply(
      seq_along(x),
      function(id) {
        # PROGRESS MESSAGE
        message(
          "[",
          id,
          "/",
          length(x),
          "]: ",
          names(x)[id],
          "\n"
        )
        # GATINGSET
        gs <- x[[id]]
        # NODES
        nodes <- cyto_nodes(
          gs, 
          path = "auto"
        )
        # MATCH CHANNELS
        pd <- cyto_channel_match(
          gs,
          channels = channels,
          save_as = details
        )
        pd <- pd[
          match(
            rownames(cyto_details(gs)), 
            rownames(pd)
          ),
          , 
          drop = FALSE
        ]
        # LOCATE SELECTED UNSTAINED CONTROLS
        unst_idx <- which(
          pd$channel %in% c("Unstained", "unstained", NA, "NA", "") &
            pd$select %in% "TRUE"
        )
        if(length(unst_idx) > 0) {
          names(unst_idx) <- cyto_names(gs)[unst_idx]
        }
        # LOCATE SELECTED STAINED CONTROLS
        idx <- which(
          !pd$channel %in% c("Unstained", "unstained", NA, "NA", "") &
          pd$select %in% "TRUE"
        )
        names(idx) <- cyto_names(gs)[idx]
        # CHECK SPECTRA
        rerun <- FALSE
        if(is.null(spectra)) {
          spec <- NULL
          rerun <- TRUE
        } else {
          spec <- spectra[[id]]
        }
        # CHECK WHETHER POSITIVE & NEGATIVE GATES EXIST
        dyes <- paste0(
          "<",
          pd$marker[idx],
          "> ",
          pd$label[idx]
        )
        # REPLACE ILLEGAL CHARACTERS
        dyes <- gsub(
          "[\\|\\&|\\:|\\,|\\/]",
          "-",
          dyes
        )
        # SEARCH FOR GATES
        LAPPLY(
          dyes, 
          function(dye) {
            if(!any(grepl(paste0(dye, "\\-"), nodes))) {
              rerun <<- TRUE
            }
            if(!any(grepl(paste0(dye, "\\+"), nodes))) {
              rerun <<- TRUE
            }
          }
        )
        # RECOMPUTE SPECTRA & GET GATES
        # NOTE: CHANNEL MATCHING WILL BE PERFORMED TWICE
        if(rerun) {
          if(spillover) {
            spectra <- cyto_spillover_compute(
              gs,
              type = type,
              heatmap = FALSE,
              save_as = NA,
              ...
            )
          } else {
            if(!all(gsub("<.*> ", "", rownames(spec)) %in% cyto_channels(gs))) {
              stop("GatingSets should be unmixed for spectral data!")
            }
            if(is.null(spec)) {
              stop(
                "An unmixing matrix must be supplied to `spectra`!"
              )
            }
            # RUN CYTO_SPILLOVER_COMPUTE() TO GET POS|NEG GATES ON UNMIXED SCALE
            cyto_spillover_compute(
              gs,
              type = type,
              heatmap = FALSE,
              save_as = NA,
              ...
            )
          }
        }
        # COMPUTE STAIN INDICES
        message("Computing stain indices...")
        stain_index <- do.call(
          "rbind",
          lapply(
            seq_along(gs)[idx],
            function(z) {
              # LOCATE GROUP UNSTAINED
              unst <- which(
                pd$group == pd$goup[z] &
                pd$channel %in% c("Unstained", "unstained", "NA", NA, "") &
                pd$select == "TRUE"
              )
              if(length(unst) == 0) {
                unst <- z
              }
              # DYE NAME
              dye <- paste0(
                "<",
                pd$marker[z],
                "> ",
                pd$label[z]
              )
              # REPLACE ILLEGAL CHARACTERS
              dye <- gsub(
                "[\\|\\&|\\:|\\,|\\/]",
                "-",
                dye
              )
              # COMPUTE NEGATIVE MEDFI
              neg_medfi <- cyto_apply(
                gs[unst],
                FUN = "cyto_stat_median",
                round = 3,
                parent = cyto_nodes_convert(
                  gs,
                  nodes = paste0(dye, "-"),
                  anchor = pd$parent[z]
                ),
                channels = pd$channel[z],
                input = "matrix",
                inverse = TRUE,
                copy = TRUE
              )
              # COMPUTE NEGATIVE RSD
              neg_rsd <- cyto_apply(
                gs[unst],
                FUN = "cyto_stat_rsd",
                round = 3, 
                parent = cyto_nodes_convert(
                  gs,
                  nodes = paste0(dye, "-"),
                  anchor = pd$parent[z]
                ),
                channels = pd$channel[z],
                input = "matrix",
                inverse = TRUE,
                copy = TRUE
              )
              # COMPUTE POSITIVE MEDFI
              pos_medfi <- cyto_apply(
                gs[z],
                FUN = "cyto_stat_median",
                round = 3,
                parent = cyto_nodes_convert(
                  gs,
                  nodes = paste0(dye, "+"),
                  anchor = pd$parent[z]
                ),
                channels = pd$channel[z],
                input = "matrix",
                inverse = TRUE,
                copy = TRUE
              )
              # COMPUTE STAIN INDEX
              return(
                c(
                  "neg_medFI" = neg_medfi,
                  "neg_rsd" =  neg_rsd,
                  "pos_medFI" = pos_medfi,
                  "RSI" = (pos_medfi - neg_medfi)/(2*neg_rsd)
                )
              )
            }
          )
        )
        stain_index <- cbind(
          cyto_details(gs)[seq_along(gs)[idx], , drop = FALSE],
          stain_index
        )
        # SORT BY STAIN INDEX
        stain_index <- stain_index[
          order(stain_index[, "RSI"], decreasing = TRUE),
          , 
          drop = FALSE
        ]
        # LOCATE BRIGHTEST CONTROL PER DYE
        keep_idx <- sapply(
          unique(stain_index[, "label"]),
          function(z) {
            idx <- match(z, stain_index[, "label"])
            idx[which.max(stain_index[idx, "RSI"])]
          }
        )
        # SPILLOVER MATRICES WE NEED TO DROP DETECTORS WITHOUT DYES
        chans <- channels
        if(spillover) {
          chans <- colnames(spec)[
            match(
              stain_index[, "channel"],
              colnames(spec)
            )
          ]
        }
        # STAIN INDEX REDUCTION
        message("Computing stain index reduction matrix...")
        # NOTE: REQUIRE COMPENSATED OR UNMIXED DATA WITH GATES
        # APPLY COMPENSATION - NOTE UNMIXING ALREADY APPLIED
        if(spillover){
          if(is.null(cyto_spillover_extract(gs))) {
            gs <- suppressMessages(
              cyto_compensate(
                gs,
                spillover = spec,
                copy = TRUE
              )
            )
          }
        }
        # WE ALREADY HAVE STAIN INDEX PER PEAK - WE NEED SECONDARY STAIN INDEX
        reduce_stain_index <- do.call(
          "rbind",
          lapply(
            keep_idx,
            function(z) {
              sapply(
                keep_idx,
                function(w) {
                  if(z == w) {
                    return(0)
                  }
                  # LOCATE NEGATIVE POPULATION FOR DYE Z
                  gs_pos_idx_z <- match(
                    rownames(stain_index)[z], 
                    rownames(pd)
                  )
                  # LOCATE POSITIVE CONTROL FOR DYE W
                  gs_pos_idx_w <- match(
                    rownames(stain_index)[w], 
                    rownames(pd)
                  )
                  # DYE Z NAME
                  dye_z <- paste0(
                    "<",
                    stain_index[z, "marker"],
                    "> ",
                    stain_index[z, "label"]
                  )
                  # REPLACE ILLEGAL CHARACTERS
                  dye_z <- gsub(
                    "[\\|\\&|\\:|\\,|\\/]",
                    "-",
                    dye_z
                  )
                  # DYE W NAME
                  dye_w <- paste0(
                    "<",
                    stain_index[w, "marker"],
                    "> ",
                    stain_index[w, "label"]
                  )
                  # REPLACE ILLEGAL CHARACTERS
                  dye_w <- gsub(
                    "[\\|\\&|\\:|\\,|\\/]",
                    "-",
                    dye_w
                  )
                  # LOOK UP POSITIVE MEDFI
                  pos_medfi <- stain_index[w, "pos_medFI"]
                  # COMPUTE NEGATIVE MEDFI - SECONDARY DETECTOR
                  neg_medfi <- cyto_apply(
                    gs[gs_pos_idx_z],
                    FUN = "cyto_stat_median",
                    round = 3,
                    parent = cyto_nodes_convert(
                      gs,
                      nodes = paste0(dye_z, "+"),
                      anchor = pd$parent[gs_pos_idx_z]
                    ),
                    channels = pd$channel[gs_pos_idx_w],
                    input = "matrix",
                    inverse = TRUE,
                    copy = TRUE
                  )
                  # COMPUTE NEGATIVE RSD - SECONDARY DETECTOR
                  neg_rsd <- cyto_apply(
                    gs[gs_pos_idx_z],
                    FUN = "cyto_stat_rsd",
                    round = 3,
                    parent = cyto_nodes_convert(
                      gs,
                      nodes = paste0(dye_z, "+"),
                      anchor = pd$parent[gs_pos_idx_z]
                    ),
                    channels = pd$channel[gs_pos_idx_w],
                    input = "matrix",
                    inverse = TRUE,
                    copy = TRUE
                  )
                  dp_stain_index <- (pos_medfi - neg_medfi)/(2*neg_rsd)
                  # DP_STAIN_INDEX < 0 - MAX IMPACT
                  dp_stain_index[dp_stain_index < 0] <- 0
                  dp_stain_index <- 1 - (dp_stain_index/stain_index[w, "RSI"])
                  # DP_STAIN_INDEX > REFERENCE - NO IMPACT
                  dp_stain_index[dp_stain_index < 0] <- 0
                  # PROPORTION OF ORIGINAL STAIN INDEX LOST
                  return(
                    dp_stain_index
                  )
                }
              )
            }
          )
        )
        message("Computing spread...")
        SSM <- list()
        TSSM <- list()
        lapply(
          idx,
          function(z) {
            # LOCATE GROUP UNSTAINED
            unst <- which(
              pd$group == pd$goup[z] &
              pd$channel %in% c("Unstained", "unstained", "NA", NA, "") &
              pd$select == "TRUE"
            )
            if(length(unst) == 0) {
              unst <- z
            }
            # DYE NAME
            dye <- paste0(
              "<",
              pd$marker[z],
              "> ",
              pd$label[z]
            )
            # REPLACE ILLEGAL CHARACTERS
            dye <- gsub(
              "[\\|\\&|\\:|\\,|\\/]",
              "-",
              dye
            )
            # COMPUTE NEGATIVE RSD
            neg_sd <- cyto_apply(
              gs[unst],
              FUN = "cyto_stat_quantile",
              probs = c(0.84, 0.5),
              round = 3,
              parent = cyto_nodes_convert(
                gs,
                nodes = paste0(dye, "-"),
                anchor = pd$parent[z]
              ),
              channels = unique(pd$channel[idx]),
              input = "matrix",
              inverse = TRUE,
              copy = TRUE
            )
            neg_medfi <- neg_sd[2, pd$channel[z]]
            neg_sd <- neg_sd[1, ] - neg_sd[2, ]
            # COMPUTE POSITIVE RSD
            pos_sd <- cyto_apply(
              gs[z],
              FUN = "cyto_stat_quantile",
              probs = c(0.84, 0.5),
              round = 3,
              parent = cyto_nodes_convert(
                gs,
                nodes = paste0(dye, "+"),
                anchor = pd$parent[z]
              ),
              channels = unique(pd$channel[idx]),
              input = "matrix",
              inverse = TRUE,
              copy = TRUE
            )
            pos_medfi <- pos_sd[2, pd$channel[z]]
            pos_sd <- pos_sd[1, ] - pos_sd[2, ]
            # COMPUTE SPREAD
            spread <- (pos_sd^2) - (neg_sd^2)
            spread[match(pd$channel[z], unique(pd$channel[idx]))] <- 0
            spread[spread < 0] <- 0
            spread <- sqrt(spread)
            # TOTAL SPILLOVER SPREAD
            TSSM[[cyto_names(gs)[z]]] <<- spread
            SSM[[cyto_names(gs)[z]]] <<- spread/sqrt(pos_medfi - neg_medfi)
            return(NULL)
          }
        )
        # PREPARE SSM
        SSM <- do.call("rbind", SSM)
        rownames(SSM) <- cyto_names(gs)[idx]
        # SORT SSM
        ord <- order(
          match(
            pd$channel[idx], 
            unique(pd$channel[idx])
          )
        )
        SSM <- SSM[ord, , drop = FALSE]
        # PREPARE TSSM
        TSSM <- do.call("rbind", TSSM)
        rownames(TSSM) <- cyto_names(gs)[idx]
        # SORT TSSM
        TSSM <- TSSM[ord, , drop = FALSE]
        # DROP EMPTY SPECTRA
        spec <- spec[rowSums(spec) != 1, , drop = FALSE]
        # SPECTRAL SIMILARITY
        message("Computing cosine similarities...")
        sim <- cyto_spectra_compare(
          spec,
          type = "cosine",
          heatmap = FALSE
        )
        # PANEL COMPLEXITY
        message("Computing panel complexity...")
        complexity <- kappa(spec)
        # HOTSPOT MATRIX
        message("Computing unmixing error hotspot matrix...")
        hotspot <- sqrt(abs(solve(sim)))
        # PURITY
        message("Computing spectral purity...")
        purity <- CytoExploreR:::row_purity_cpp(spec)
        purity <- purity[order(purity, decreasing = TRUE)]
        # PANEL RESULTS
        return(
          list(
            "spectra" = spec,
            "similarity" = sim, 
            "complexity" = complexity, 
            "brightness" = stain_index, 
            "interference" = reduce_stain_index, 
            "spread" = SSM, 
            "total_spread" =  TSSM, 
            "purity" = purity,  
            "hotspot" = hotspot 
          )
        )
      }
    ),
    names = names(x)
  )
  
  # RETURN COMPUTED STATISTICS
  return(stats)
  
}