## CYTO_GATE_TRANSFER ----------------------------------------------------------

#' Transfer all gates from one GatingHierarchy/GatingSet to another
#'
#' @param from object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} from which the gates
#'   should be extracted.
#' @param to object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} which should inherit
#'   the transferred gates.
#' 
#' @return a GatingHierarchy or GatingSet containing the transferred gates.
#' 
#' @importFrom flowWorkspace gh_pop_get_gate gs_pop_add
#' @importFrom openCyto gh_generate_template
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_gate_transfer
NULL

#' @noRd
#' @export
cyto_gate_transfer <- function(from,
                               to) {
  UseMethod("cyto_gate_transfer")
}

#' @rdname cyto_gate_transfer
#' @export
cyto_gate_transfer.GatingHierarchy <- function(from,
                                               to) {
  
  # GATINGHIERARCHY/GATINGSET SUPPORTED IN 'TO' - NO CLASS CHECKS - INTERNAL USE
  
  # GENERATE GATINGTEMPLATE
  gt <- gh_generate_template(from)
  
  # TRANSFER GATES
  if(nrow(gt) > 0) {
    # SPLIT GATINGTEMPLATE INTO CHUNKS - PARENT & DIMS
    gt_chunks <- split(
      gt,
      do.call(
        "paste0",
        c(gt[, c("parent", "dims")]) 
      )
    )
    # EXTRACT GATES
    ind <- c()
    gates <- lapply(
      gt_chunks,
      function(z) {
        # EXTRACT GATES - LIST OF GATE OBJECTS
        gate <- cyto_gate_extract(
          from,
          parent = z$parent[1],
          alias = z$alias,
          bool = FALSE # LEAVE BOOLEANFILTERS ALONE
        )[[1]]
        # ORDER
        ind <<- c(
          ind,
          LAPPLY(
            seq_along(
              gate
            ),
            function(v) {
              # MULTIPLE POPULATIONS PER GATE
              a <- strsplit(
                names(gate)[v],
                "\\|"
              )[[1]]
              # MINIMUM ROW INDEX FOR MULTIGATES
              min(
                which(
                  gt$alias %in% a &
                    gt$parent == z$parent[1] &
                      gt$dims == z$dims[1]
                )
              )
            }
          )
        )
        # RETURN GATE OBJECTS
        return(gate)
      }
    )
    # ORDER GATES
    gates <- structure(
      unlist(
        gates,
        recursive = FALSE
      ),
      names = LAPPLY(gates, "names")
    )[order(ind)]
    ind <- ind[order(ind)]
    # TRANSFER GATES
    lapply(
      seq_along(gates),
      function(z) {
        # ADD GATES TO NEW GATINGHIERARCHY
        gs_pop_add(
          to,
          gate = gates[[z]],
          parent = gt$parent[ind[z]],
          name = strsplit(
            names(gates)[z],
            "\\|"
          )[[1]], # WATCH OUT MULTIGATES
          validityCheck = FALSE
        )
      }
    )
  }
  
  # GATING - RECOMPUTE
  suppressMessages(recompute(to))
  
  # RETURN UPDATED GATINGHIERARCHY
  return(to)

}

#' @rdname cyto_gate_transfer
#' @export
cyto_gate_transfer.GatingSet <- function(from, 
                                         to) {
  
  # TRANSFER GATES
  lapply(
    seq_along(from),
    function(z) {
      cyto_gate_transfer(
        from[[z]],
        to[[z]]
      )
    }
  )
  
  # UPDATED GATINGSET
  return(to)
  
}
