## CYTO_GATE_CLUSTER -----------------------------------------------------------

# Apply clustering algorithm to data
# Pass data through cyto_apply() and apply clustering function to each group.
# Cluster data and assign labels for GatingSet
# Prepare gatingTemplate entries
# Expect a numeric/factor/category label per cluster

#' @noRd
cyto_gate_cluster <- function(x,
                              parent = "root",
                              FUN = NULL,
                              input = "cytoframe",
                              copy = TRUE,
                              channels = channels,
                              trans = NA,
                              inverse = FALSE,
                              slot = NULL,
                              merge_by = "all",
                              ...) {
  
  # SPLIT SAMPLES INTO GROUPS
  x <- cyto_group_by(
    x,
    group_by = "merge_by"
  )
  
  # PERFORM CLUSTERING ON EACH GROUP
  res <- structure(
    lapply(
      x,
      function(z) {
        # APPLY CLUSTERING ALGORITHM
        gate <- cyto_apply(
          z,
          parent = parent,
          input = input,
          channels = channels,
          copy = copy,
          trans = trans,
          inverse = inverse,
          slot = slot,
          FUN = FUN,
          ...
        )
        # PREPARE GATE TO APPLY TO GATINGSET -> FACTOR FILTERRESULT
        
        
      }
    ),
    names = names(x)
  )

  
}