## CYTO_PLOT_GATING_TREE -------------------------------------------------------

#' Plot Gating Trees
#'
#' \code{cyto_plot_gating_tree} provides a simpler visualisation of the gating
#' scheme for \code{GatingHierarchy}, \code{GatingSet} and \code{gatingTemplate}
#' objects. The \code{GatingHierachy} method is also capable of displaying
#' population statistics such as frequency of parent or count.
#'
#' @param x object of class \code{GatingHierarchy}, \code{GatingSet} or
#'   \code{gatingTemplate}.
#' @param stat used in \code{GatingHierachy} method to add either "percent" or
#'   "count" statistics onto the gating tree, set to NULL by default to exclude
#'   statistics.
#'
#' @importFrom openCyto gh_generate_template
#' @importFrom magrittr %>%
#' @importFrom visNetwork visNetwork visEdges
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @examples 
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # Add samples to GatingSet
#' gs <- GatingSet(fs)
#' 
#' # Apply compensation
#' gs <- cyto_compensate(gs)
#' 
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs, select = "Stim-D", trans_type = "logicle")
#' 
#' # Gating
#' gt <- Activation_gatingTemplate
#' cyto_gatingTemplate_apply(gs, gt)
#' 
#' # Visualise gating tree using gatingTemplate
#' cyto_plot_gating_tree(gt)
#' 
#' # Visualise gating tree for GatingSet (same output as gatingTemplate)
#' cyto_plot_gating_tree(gs)
#' 
#' # Visualise gating tree for GatingHierarchy
#' cyto_plot_gating_tree(gs[[32]], stat = "percent")
#' cyto_plot_gating_tree(gs[[32]], stat = "count")
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
cyto_plot_gating_tree.GatingHierarchy <- function(x,
                                                  stat = NULL) {
  
  # Extract gatingTemplate from GatingHierarchy
  gt <- gh_generate_template(x)
  
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
      stats <- paste(stats, "%")
      
    }
    
    # Add value column to adjust node sizes
    nodes$value <- node_counts$Count/node_counts$Count[1]
    
    # Add percent labels to edges
    edges$label <- stats[-1]
    
  }

  # Call to visNetwork
  visNetwork(nodes, edges) %>%
    visEdges(arrows = "to", color = "black")
  
}

#' @rdname cyto_plot_gating_tree
#' @export
cyto_plot_gating_tree.GatingSet <- function(x) {
  
  # Generate template based on first sample
  gt <- gh_generate_template(x[[1]])
  
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
  
  # Call to visNetwork
  visNetwork(nodes, edges) %>%
    visEdges(arrows = "to", color = "black")
  
}

#' @rdname cyto_plot_gating_tree
#' @export
cyto_plot_gating_tree.gatingTemplate <- function(x){
  
  # Convert gatingTemplate to data.table
  gt <- as.data.table.gatingTemplate(x)
  
  # Preprocess gatingTemplate 
  gt <- .preprocess_csv(gt)
  
  # Convert preprocessed gt to data.frame
  gt <- as.data.frame(gt[, c("alias","parent"), with = FALSE])
  
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
  
  # Call to visNetwork
  visNetwork(nodes, edges) %>%
    visEdges(arrows = "to", color = "black")
  
}
