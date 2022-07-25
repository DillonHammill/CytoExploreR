## CYTO_REBUILD ----------------------------------------------------------------

#' Rebuild a GatingSet onto an extracted cytoset
#'
#' @param x reference object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing all
#'   events.
#' @param y extracted object of class
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}.
#' @param inverse indicates whether the extracted data has already be inverse
#'   transformed onto a linear scale, set to FALSE by default.
#'
#' @return a rebuilt GatingSet.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom purrr transpose
#' @importFrom flowWorkspace GatingSet
#'
#'
#' @noRd
cyto_rebuild <- function(x,
                         y,
                         inverse = FALSE) {
  
  # X - REFERENCE GATINGSET
  # Y - DOWNSAMPLED GATINGSET OR CYTOSET
  
  # CHECKS ---------------------------------------------------------------------
  
  # X - CLASS
  if(!cyto_class(x, "GatingSet")) {
    stop(
      "'x' must be a reference GatingSet!"
    )
  }
  
  # Y - CLASS
  if(!cyto_class(x, c("GatingSet")) ) {
    stop(
      "'y' must be either a downsampled cytoset or GatingSet!"
    )
  }
  
  # MATCHING SAMPLES REQUIRED
  if(!all(cyto_names(y) %in% cyto_names(x))) {
    stop(
      "'y' must contain the same samples as the reference GatingSet 'x'!"
    )
  }
  
  # REFERENCE COMPENSATION
  spill <- cyto_spillover_extract(x)
  
  # RESTRICT TO SAMPLES
  if(!is.null(spill)) {
    spill <- spill[cyto_names(y)]
  }
  
  # REFERENCE TRANSFORMERS
  trans <- cyto_transformers_extract(x)
  
  # MESSAGE
  message(
    "Building a new GatingSet..."
  )
  # INVERSE TRANSFORM
  if(!.all_na(trans) & !inverse) {
    y <- cyto_transform(
      y,
      trans = trans,
      inverse = TRUE,
      plot = FALSE,
      quiet = TRUE
    )
  }
  # DECOMPENSATE
  if(!is.null(spill)) {
    y <- cyto_compensate(
      y,
      spill = spill,
      remove = TRUE,
      quiet = TRUE
    )
  }
  # NEW GATINGSET
  y <- GatingSet(y)
  # RE-APPLY COMPENSATION
  if(!is.null(spill)) {
    message(
      "Compensating for fluorescent spillover..."
    )
    y <- cyto_compensate(
      y,
      spillover = spill,
      remove = FALSE,
      quiet = TRUE
    )
  }
  # APPLY TRANSFORMATIONS
  if(!.all_na(trans)) {
    message(
      "Re-applying data transformations..."
    )
    y <- cyto_transform(
      y,
      trans = trans,
      inverse = FALSE,
      plot = FALSE,
      quiet = TRUE
    )
  }
  # GATING
  message(
    "Recomputing gates..."
  )
  # GATETEMPLATE
  gt <- transpose(
    cyto_gateTemplate(
      cyto_select(
        x,
        cyto_names(y)
      )
    )
  )
  # APPLY GATES
  lapply(
    seq_along(gt),
    function(z) {
      gates <- lapply(gt[[z]], `[[`, "gate")
      parent <- unique(LAPPLY(gt[[z]], `[[`, "parent"))
      alias <- unique(LAPPLY(gt[[z]], `[[`, "alias"))
      type <- unique(LAPPLY(gates, "cyto_class"))
      # BOOLEAN GATES - USE LOGICAL INDICES
      if(type %in% "booleanFilter") {
        gates <- structure(
          lapply(
            cyto_names(y),
            function(z) {
              # REFERENCE EVENT-IDs
              x_ids <- cyto_data_extract(
                x,
                select = z,
                parent = parent,
                channels = "Event-ID",
                format = "matrix",
                copy = FALSE
              )[[1]][[1]][, 1]
              # SUBSET EVENT IDs
              x_ids <- x_ids[
                cyto_gate_indices(
                  x,
                  select = z,
                  parent = parent,
                  nodes = alias
                )[[1]][, 1]
              ]
              # NEW EVENT-IDs
              y_ids <- cyto_data_extract(
                y,
                select = z,
                parent = parent,
                channels = "Event-ID",
                format = "matrix",
                copy = FALSE
              )[[1]][[1]][, 1]
              # LOGICAL VECTOR
              gate <- rep(FALSE, length(y_ids))
              gate[y_ids %in% x_ids] <- TRUE
              return(gate)
            }
          ),
          names = cyto_names(y)
        )
      }
      # APPLY GATES
      gs_pop_add(
        y,
        gate = gates,
        parent = parent,
        name = alias,
        validityCheck = FALSE
      )
      # RECOMPUTE
      suppressMessages(
        recompute(y)
      )
      
    }
  )
  
  # RETURN REBUILT GATINGSET
  return(y)
  
}