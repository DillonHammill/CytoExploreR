#' convert a gatingTemplate object to a data.table
#' 
#' It is the inverse function of gatingTemplate constructor.
#' 
#' @param x gatingTemplate object
#' @param keep.rownames not used
#' @return a data.table
#' @importFrom data.table data.table rbindlist
#' @importFrom methods extends
#' @importFrom openCyto gt_get_gate gt_get_parent
#' @noRd
as.data.table.gatingTemplate <- function(x, keep.rownames = FALSE){
  
  # gate each node 
  gt_nodes <- gt_get_nodes(x, order = "tsort")[-1]
  
  res <- lapply(gt_nodes, function(node){
    
    # get parent node to gate
    nodePath <- node@id
    parent <- gt_get_parent(x, nodePath)
    # extract gate method from one edge(since multiple edge to the same node is
    # redudant)
    this_gate <- gt_get_gate(x, parent, nodePath)
    gating_args <- parameters(this_gate)
    #collapse into string when neccessary
    split <- !extends(class(this_gate), "refGate")
    gating_args <- .argDeparser(gating_args, split)
    
    #get preprocessing method
    this_ppm <- ppMethod(x, parent, nodePath)
    
    #preprocessing
    if(class(this_ppm) == "ppMethod")
    {
      preprocessing_method <- names(this_ppm)
      preprocessing_args <- parameters(this_ppm)
      preprocessing_args <- .argDeparser(preprocessing_args)
    }else{
      preprocessing_args <- preprocessing_method <- ""
    }
    
    data.table(alias = .alias(node)
               , pop = names(node)
               , parent = parent
               , dims = this_gate@dims
               , gating_method = names(this_gate)
               , gating_args =  gating_args
               , collapseDataForGating = .isCollapse(this_gate)
               , groupBy = .groupBy(this_gate)
               , preprocessing_method = preprocessing_method
               , preprocessing_args = preprocessing_args
    )    
    
  })
  
  rbindlist(res)
  
}

#' deparse a list(named) of expression into a string inverse function of
#' argParser
#'
#' @param args gatingTemplate arguments as named list to be deparsed (e.g.
#'   list(gate = gateobj)).
#' @param split logical.
#'
#' @return deparsed arguments for storage in gatingTemplate csv file
#'
#' @noRd
.argDeparser <- function(args, split = TRUE) {
  if (split) {
    args <- LAPPLY(names(args), function(argn) {
      argv <- deparse(args[[argn]])
      argv <- gsub("\"", "'", argv) # restore dquote to squote
      argv <- paste(argv, collapse = "")
      paste(argn, argv, sep = " = ")
    })


    paste(args, collapse = ", ")
  } else {
    as.character(args[[1]])
  }
}