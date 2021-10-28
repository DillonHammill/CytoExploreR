## CYTO_GATE_CLUST -------------------------------------------------------------

#' Use clustering algorithm to gate populations
#' 
#' @param x object of class \code{GatingSet}.
#' 
#' @noRd
cyto_gate_clust <- function(x,
                            parent = NULL,
                            alias = "Cluster",
                            channels = NULL,
                            type = NULL, 
                            group_by = NULL,
                            gatingTemplate = NULL,
                            ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # ALIAS MISSING
  if (is.null(alias)) {
    stop("Supply the name(s) for the gated population(s) to 'alias'.")
  }
  
  # CLUSTERING FUNCTION
  if(is.null(type)){
    stop("Supply the name of the clustering function to 'type'.")
  }
  
  # ACTIVE GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    gatingTemplate <- cyto_gatingTemplate_active()
  }
  
  # MISSING GATINGTEMPLATE
  if (is.null(gatingTemplate)) {
    stop("Supply the name of the gatingTemplate csv file to save the gate(s),")
  }
  
  # EXISTING ENTRIES IN GATINGTEMPLATE
  gt <- .cyto_gatingTemplate_check(parent, alias, gatingTemplate)
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
  # CLUSTERING
  pops <- suppressWarnings(
    gs_add_gating_method(
      gs = x,
      alias = paste(alias, collapse = ","),
      parent = parent,
      pop = "*",
      dims = paste(channels, collapse = ","),
      gating_method = "cyto_gate_clust",
      gating_args = as.list(...),
      groupBy = group_by,
      collapseDataForGating = TRUE
    )
  )
  
  # COMBINE NEW ENTRIES WITH EXISTING ONES
  gt <- rbind(gt, pops)
  
  # SAVE UPDATED GATINGTEMPLATE
  write_to_csv(gt, gatingTemplate)
  
  # RETURN GATINGSET
  return(x)
  
}

#' openCyto plugin
#' @noRd
.cyto_gate_clust <- function(fr,
                             pp_res = NULL,
                             channels = NULL,
                             type = NULL,
                             ...){
  
  # MESSAGE
  message(
    paste(
      "Using", type, "to cluster populations."
    )
  )
  
  # CLUSTERING FUNCTION
  cluster <- eval(parse(text = type))
  
  # CLUSTERING RESULT
  cluster_result <- cluster(
    cyto_data_extract(
      fr,
      raw = TRUE,
      channels = channels
      ), 
    ...)
  
  # PULL OUT CLUSTERS
  if(is.vector(cluster_result)){
    clusters <- cluster_result
    # LOCATE CLUSTERING VECTOR
  }else{
    clusters <- c()
    for(i in length(cluster_result)){
      if(length(clusters) > 0){
        break()
      }
      # CLUSTER VECTOR - WHOLE NUMBERS
      if(length(cluster_result[i]) == nrow(cyto_exprs(fr)) &
         all(cluster_result[i] %% 1 == 0)){
        clusters <<- cluster_result[[i]]
      }
    }
  }
  
  # FACTORISE
  clusters <- as.factor(clusters)
  
  # FACTOR LEVELS
  levels(clusters) <- paste("cluster", seq_len(length(levels(clusters))))
  
  # RETURN FILTER RESULT
  as(clusters, "filterResult")
  
}
