# CYTO_PLOT_GATING_TREE --------------------------------------------------------

#' Plot Gating Trees
#' 
#' @importFrom openCyto templateGen
#' @importFrom magrittr %>%
#' @importFrom visNetwork visNetwork visEdges
#' 
#' @name cyto_plot_gating_tree
NULL

#' @noRd
#' @export
cyto_plot_gating_tree <- function(x, ...){
  UseMethod("cyto_plot_gating_tree")
}

#' @rdname cyto_plot_gating_tree
#' @export
cyto_plot_gating_tree.GatingSet <- function(x, ...) {
  
}

#' @rdname cyto_plot_gating_tree
#' @export
cyto_plot_gating_tree.GatingHierarchy <- function(x,
                                                  stat = NULL,
                                                  ...) {
  
  # Extract gatingTemplate from GatingHierarchy
  gt <- templateGen(x)
  
  # Extract names of all nodes
  pops <- c("root", gt$alias)
  
  # Extract nodes
  nodes <- rbind("root", gt[, c("alias","alias"), drop = FALSE])
  colnames(nodes) <- c("id","label")
  
  # Add group column for colours
  nodes$group <- nodes$id
  
  # Extract alias and parent columns
  edges <- gt[, c("alias","parent")]
  
  # Rename columns  for visNetwork
  colnames(edges) <- c("to", "from")
  
  # Convert parent to basename
  edges[, "from"] <- basename(edges[, "from"])
  
  # Scale nodes by frequency & add labels
  if(!is.null(stat)){
    
    # Calculate counts for each node
    node_counts <- cyto_stats_compute(x,
                                      alias = nodes$id,
                                      stat = "count",
                                      format = "long")
    
    # Normalise as a percentage of "root"
    if(stat == "count"){
      
      # Extract counts
      stats <- node_counts$Count
  
    }
    
    # Normalise as a percentage of parent
    if(stat == "percent"){
      
      # Order counts based on parent names
      stats <- node_counts$Count/
        node_counts[match(c("root", edges$from), 
                          node_counts$Population), "Count"] * 100
      stats <- LAPPLY(stats, function(z){.round(z, 2)})
      print(stats)
      stats <- paste(stats, "%")
      
    }
    
    print(stats)
    
    # Add value column to adjust node sizes
    nodes$value <- node_counts$Count/node_counts$Count[1]
    
    # Add percent labels to edges
    edges$label <- stats[-1]
    
    print(edges)
    
  }

  # Call to visNetwork
  visNetwork(nodes, edges) %>%
    visEdges(arrows = "to", color = "black")
  
}

#' @rdname cyto_plot_gating_tree
#' @export
cyto_plot_gating_tree.gatingTemplate <- function(x, ...){
  
}