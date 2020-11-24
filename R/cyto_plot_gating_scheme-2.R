## CYTO_PLOT_GATING_SCHEME -----------------------------------------------------

# Cluster channels

#' Plot cytometry gating strategies
#' @name cyto_plot_gating_scheme2
NULL

#' @export
#' @noRd
cyto_plot_gating_scheme2 <- function(x, 
                                    ...){
  UseMethod("cyto_plot_gating_scheme")
}

#' @rdname cyto_plot_gating_scheme
#' @export
cyto_plot_gating_scheme2.GatingHierarchy <- function(x,
                                                    overlay = NA,
                                                    gatingTemplate = NULL,
                                                    popup = FALSE,
                                                    legend = TRUE,
                                                    ...){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- .par(c("mfrow", "oma"))
  
  # CHECKS ---------------------------------------------------------------------
  
  
  
}

#' Split gatingTemplate into gating chunks
#' @param x either  GatingHierarchy/GatingSet
#' @param drop logical indicating whether to drop boolean gates that cannot be
#'   plotted, set to TRUE by default.
#' @noRd
.cyto_gatingTemplate_split <- function(x, 
                                       drop = TRUE){
  
  # GATINGTEMPLATE
  gt <- .cyto_gatingTemplate_generate(x)
  
  # BOOLEAN GATES - GATING_ARGS ADDED ABOVE
  bool_ind <- which(LAPPLY(gt$dims, ".empty") & 
                      !LAPPLY(gt$gating_args, "is.na"))
  if(length(bool_ind) > 0){
    gt_bool <- gt[bool_ind, ]
    gt <- gt[-bool_ind, ]
  }else{
    gt_bool <- NULL
  }
  
  # CLUSTERS - NO GATING_ARGS
  clust_ind <- which(LAPPLY(gt$dims, ".empty") & 
                       LAPPLY(gt$gating_args, "is.na"))
  if(length(clust_ind) > 0){
    gt_clust <- gt[clust_ind, ]
    gt <- gt[-clust_ind, ]
  }else{
    gt_clust <- NULL
  }
  
  # SPLIT BY PARENT/CHANNEL COMBINATION
  gt$plot <- paste0(gt$parent, "_", gt$dims)
  gt_chunks <- lapply(unique(gt$plot), function(z){
    gt[gt$plot == z, -match("plot", colnames(gt))] # remove plot column
  })
  
  # AVAILABLE POPULATIONS AFTER EACH CHUNK
  pops <- c("root")
  gt_pops <- lapply(gt_chunks, function(z){
    pops <<- c(pops, z[, "alias_auto"])
    return(pops)
  })
  
  # CHECK BOOLEAN ENTRIES AGAINST GT_POPS
  if(!is.null(gt_bool)){
    lapply(seq_len(nrow(gt_bool)), function(z){
      gt_parent <- cyto_nodes_convert(x,
                                      nodes = gt_bool[z, "parent"],
                                      path = "auto")
      gt_alias <- unlist(strsplit(gt_bool[z, "gating_args"], "\\||\\&|\\!"))
      gt_alias <- gt_alias[!LAPPLY(gt_alias, ".empty")]
      gt_ind <- min(which(LAPPLY(gt_pops, function(y){
        all(c(gt_parent, gt_alias) %in% y)
      })))
      # INCLUDE BOOLEAN GATE
      #if(){
      gt_chunks[[gt_ind]] <<- rbind(gt_chunks[[gt_ind]], gt_bool[z,])
      #}
      
    })
  }
  
  # CLUSTERS
  if(!is.null(gt_clust)){
    # SPLIT BY PARENT
    gt_clust <- lapply(unique(gt_clust$parent), function(z){
      return(gt_clust[gt_clust$parent == z, ])
    })
    # ADD TO END
    gt_chunks <- c(gt_chunks, gt_clust)
  }
  
  # RETURN SPLIT GATINGTEMPLATE
  return(gt_chunks)
  
}

#' Generate gatingTemplate and include boolean logic
#' @param x GatingHierarchy or GatingSet.
#' @importFrom openCyto gh_generate_template
#' @noRd
.cyto_gatingTemplate_generate <- function(x){
  
  # GATINGHIERARCHY/SET
  if(is(x, "GatingSet")){
    gh <- x[[1]]
  }else{
    gh <- x
  }
  
  # GATINGTEMPLATE
  gt <- gh_generate_template(gh)
  
  # NODES
  nodes_auto <- cyto_nodes(gh,
                           path = "auto")
  
  # UNIQUE ALIAS
  gt$alias_auto <- cyto_nodes_convert(gh,
                                      nodes = paste0(gt$parent,
                                                     "/",
                                                     gt$alias),
                                      path = "auto")
  
  # SORT GATINGTEMPLATE EXTRIES
  gt <- gt[match(nodes_auto[-match("root", nodes_auto)], gt$alias_auto), ]
  
  # CHECK FOR BOOLEAN/CLUSTER ENTRIES
  bool_ind <- which(LAPPLY(gt$dims, ".empty"))
  
  # BOOLEAN/CLUSTER ENTRIES EXIST
  if(length(bool_ind) != 0){
    # BOOLEAN GATING ARGUMENTS
    lapply(bool_ind, function(z){
      parent <- gt$parent[bool_ind]
      alias <- gt$alias[bool_ind]
      alias_auto <- cyto_nodes_convert(gh,
                                       nodes = alias,
                                       anchor = parent,
                                       path = "auto")
      gate <- gh_pop_get_gate(gh, 
                              alias_auto)
      # BOOLEAN FILTER
      if(is(gate, "booleanFilter")){
        gt[z, "gating_args"] <<- gate@deparse
        # MULTIDIMENSIONAL CLUSTER
      }else{
        gt[z, "gating_args"] <<- NA
      }
    })
  }
  
  # RETURN GATINGTEMPLATE DATA.FRAME
  return(gt)
  
}
