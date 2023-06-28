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
#' @param point_cols default colour palette from which colours should be chosen
#'   for each node in the gating tree.
#' @param point_shape vector of shapes to use for each node, can be either
#'   "circle", "ellipse", "database", "box" or "text". Set to circles by
#'   default.
#' @param point_col vectors of colours to select from if a colour is not
#'   supplied for each node through \code{point_col}.
#' @param point_col vector of colours to use for each node.
#' @param point_col_alpha fill transparency for each node, set to 0.5 by
#'   default.
#' @param ... not in use.
#'
#' @importFrom openCyto gh_generate_template CytoExploreR_.preprocess_csv
#' @importFrom rhandsontable %>%
#' @importFrom visNetwork visNetwork visEdges visGroups
#' @importFrom data.table as.data.table
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' library(CytoExploreRData)
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
                                                  stat = NULL,
                                                  point_shape = "circle",
                                                  point_cols = NA,
                                                  point_col = NA,
                                                  point_col_alpha = 0.5,
                                                  ...) {
  
  # INHERIT CYTO_PLOT_THEME ----------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)[c("x", 
                                           "stat", 
                                           "point_cols", 
                                           "point_shape")]
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # GATINGTEMPLATE -------------------------------------------------------------
  
  # EXTRACT GATINGTEMPLATE
  gt <- gh_generate_template(x)
  
  # UNIQUE ALIAS
  gt$alias_unique <- LAPPLY(seq_len(nrow(gt)), function(z){
    cyto_nodes_convert(x,
                       nodes = gt$alias[z],
                       anchor = gt$parent[z],
                       path = "auto")
  })
  
  # UNIQUE PARENT
  gt$parent_unique <- LAPPLY(seq_len(nrow(gt)), function(z){
    cyto_nodes_convert(x,
                       nodes = gt$parent[z],
                       path = "auto")
  })
  
  # NODES
  nodes <- rbind("root", gt[, c("alias_unique","alias"), drop = FALSE])
  colnames(nodes) <- c("id","label")
  
  # GROUP COLUMN FOR COLOURS 
  nodes$group <- nodes$id
  
  # EDGES
  edges <- gt[, c("alias_unique","parent_unique")]
  colnames(edges) <- c("to", "from")
  
  # SCALE NODES & ADD LABELS
  if(!is.null(stat)){
    
    # NODE COUNTS
    node_counts <- cyto_stats_compute(x,
                                      alias = nodes$id,
                                      stat = "count",
                                      format = "long")
    
    # COUNT STATISTICS
    if(stat == "count"){
      stats <- node_counts[, ncol(node_counts)]
    }
    
    # FREQUENCY STATISTICS
    if(stat %in% c("percent","freq")){
      
      # ORDER COUNTS BASED ON PARENT
      stats <- node_counts$Count/
        node_counts[match(c("root", edges$from), 
                          node_counts$Population), "Count"] * 100
      stats <- LAPPLY(stats, function(z){.round(z, 2)})
      stats <- paste(stats, "%")
      
    }
    
    # SCALE NODE SIZE TO ROOT
    nodes$value <- node_counts[, ncol(node_counts)]/
      node_counts[, ncol(node_counts)][1]
    
    # STATISTIC LABELS TO EDGES
    edges$label <- stats[-1]
    
  }
  
  # NODE BORDER COLOURS
  if(.all_na(point_col)){
    if(.all_na(point_cols)){
      point_cols <- .cyto_plot_colour_palette(type = "point_cols")
    }
    point_col_palette <- colorRampPalette(point_cols)
    point_col <- point_col_palette(nrow(nodes))
    # SOME NODE COLOURS SUPPLIED
  }else{
    # COLOURS SUPPLIED PER NODE
    if(length(point_col) == nrow(nodes)){
      
    }else{
      
    }
    
  }
  
  # NODE FILL COLOURS
  point_col_alpha <- rep(point_col_alpha, length.out = nrow(nodes))
  point_fill <- LAPPLY(seq_len(nrow(nodes)), function(z){
    adjustcolor(point_col[z], point_col_alpha[z])
  })
  
  # NODE SHAPES
  point_shape <- rep(point_shape, length.out = nrow(nodes))
  
  # NODE COLOURS
  nodes$color.background <- point_fill
  nodes$color.border <- point_col
  nodes$shape <- point_shape
  
  # VISNETWORK
  visNetwork(nodes, edges) %>%
    visEdges(arrows = "to", 
             color = "black")
  
}

#' @rdname cyto_plot_gating_tree
#' @export
cyto_plot_gating_tree.GatingSet <- function(x, ...) {
  
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
cyto_plot_gating_tree.gatingTemplate <- function(x, ...){
  
  # Convert gatingTemplate to data.table
  gt <- as.data.table(x)
  
  # Preprocess gatingTemplate 
  gt <- CytoExploreR_.preprocess_csv(gt)
  
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
