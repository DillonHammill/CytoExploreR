## CYTO_KNN --------------------------------------------------------------------

#' Compute k-nearest neighbours
#'
#' @param x a matrix containing the data to use to construct the kNN graph.
#' @param k maximum number of nearest neighbours to compute.
#' @param method indicates whether to use the \code{"kd"} or \code{"ball"} tree
#'   algorithm to construct the graph, set to \code{"kd"} by default.
#' @param n_jobs specifies the number of cores that should be used when using
#'   sklearn to compute the kNN graph, set to \code{-2L} by default to use all
#'   except one core.
#' @param ... additional arguments passed to sklearn or RANN kNN constructor.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_knn <- function(x,
                      k = 25L,
                      method = "kd_tree",
                      n_jobs = -2L,
                      ...) {
  
  # SKLEARN PYTHON MODULE - FASTER + PARALLELISATION
  skl <- cyto_require(
    "sklearn",
    python = TRUE,
    pip = TRUE,
    import = "sklearn.neighbors"
  )
  
  # SKLEARN AVAILABLE
  if(!is.null(skl)) {
    # KD - METHOD
    if(grepl("kd", method, ignore.case = TRUE)) {
      method <- "kd_tree"
    # BALL - METHOD
    } else if(grepl("ball", method, ignore.case = TRUE)) {
      method <- "ball_tree"
    }
    # KNN
    knn <- skl$NearestNeighbors(
      n_neighbors = k,
      algorithm = method,
      n_jobs = n_jobs,
      ...
    )
    knn <- knn$fit(
      data.matrix(x)
    )
    res <- knn$kneighbors(
      data.matrix(x)
    )
    names(res) <- c("nn.dists", "nn.idx")
  # RANN 
  } else {
    # RANN
    cyto_require(
      "RANN",
      source = "CRAN"
    )
    # KD - METHOD
    if(grepl("kd", method, ignore.case = TRUE)) {
      method <- "kd"
      # BALL - METHOD
    } else if(grepl("ball", method, ignore.case = TRUE)) {
      method <- "ball"
    }
    # KNN 
    res <- cyto_func_execute(
      "RANN::nn2",
      list(
        data = x,
        k = k,
        treetype = method,
        ...
      )
    )
  }
  
  # KNN
  return(res)
  
}
