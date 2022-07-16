## CYTO_IMPUTE -----------------------------------------------------------------

#' Impute values or labels using k nearest neighbours
#'
#' @param train object of class \code{matrix}, \code{data.frame},
#'   \code{cytoframe}, \code{cytoset}, \code{GatingHierarchy} or
#'   \code{GatingSet} containing the training data.
#' @param test object of class \code{matrix}, \code{data.frame},
#'   \code{cytoframe}, \code{cytoset}, \code{GatingHierarchy} or
#'   \code{GatingSet} containing the testing data.
#' @param labels a factor of labels or a vector of values for every cell within
#'   the training data. Alternatively, multiple labels can be supplied in the
#'   form of a matrix or a named list.
#' @param parent name of the population to extract from \code{train} and
#'   \code{test} if \code{GatingHierarchy} or \code{GatingSet} objects are
#'   supplied, set to the \code{"root"} node by default.
#' @param channels names of the markers or channels to use to constructed the
#'   kNN graph, must be supplied manually as no default has been set.
#' @param k number of nearest neighbours to compute when constructing the kNN
#'   graph, set to 5 by default.
#' @param method determines whether to use \code{"dist"} or \code{"vote"} when
#'   performing kNN classification, set to \code{"dist"} by default.
#' @param scale method to use to scale the channels of training and testing data
#'   prior to constructing the kNN graph, set to \code{"range"} by default.
#' @param ... not in use.
#'
#' @return a mtrix containing the imputed labels or values.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' \dontrun{
#'
#' # Traning and testing data
#' ind <- 1:nrow(iris)
#' train_ind <- sample(ind, ceiling(0.8 * nrow(iris)))
#' test_ind <- ind[!ind %in% train_ind]
#'
#' # labels
#' cyto_impute(
#'   iris[train_ind, ],
#'   iris[test_ind, ],
#'   factor(iris[train_ind, 5]),
#'   channels = colnames(iris)[1:4]
#' )
#'
#' # values
#' cyto_impute(
#'   iris[train_ind, ],
#'   iris[test_ind, ],
#'   iris[train_ind, 4],
#'   channels = colnames(iris)[1:3]
#' )
#'
#' }
#'
#' @export
cyto_impute <- function(train = NULL,
                        test = NULL,
                        labels = NULL,
                        parent = "root",
                        channels = NULL,
                        k = 5,
                        method = "dist",
                        scale = "range",
                        ...) {
  
  # TODO: MAKE EFFICIENT BY CONSTRUCTING KNN CLASSIFIER ONCE
  
  # CHANNELS
  if(is.null(channels)) {
    stop(
      "Supply the channels to be used to construct the kNN graph!"
    )
  } else {
    channels <- cyto_channels_extract(
      train,
      channels = channels
    )
  }

  # EXTRACT TRAINING DATA
  if(!cyto_class(train, c("matrix", "data.frame"))) {
    train <- cyto_data_extract(
      train,
      parent = parent,
      channels = channels,
      coerce = TRUE,
      events = 1,
      format = "matrix"
    )[[1]][[1]]
  }
  
  # EXTRACT TESTING DATA
  if(!cyto_class(test, c("matrix", "data.frame"))) {
    test <- cyto_data_extract(
      test,
      parent = parent,
      channels = channels,
      coerce = TRUE,
      events = 1,
      format = "matrix"
    )[[1]][[1]]
  }
  
  # NUMERIC DATA REQUIRED
  train <- data.matrix(train)
  test <- data.matrix(test)
  
  # SCALE
  if (!is.null(scale)) {
    train <- cyto_stat_scale(
      train,
      type = scale
    )
    test <- cyto_stat_scale(
      test,
      type = scale
    )
  }
  
  # PREPARE LABELS
  if(cyto_class(labels, c("factor", "numeric", "integer"))) {
    labels <- list(labels)
  } else if(cyto_class(labels, "character")) {
    labels <- list(factor(labels))
  } else if(!cyto_class(labels, "list")) {
    labels <- as.list(labels)
  }
  
  # IMPUTE LABELS
  labels <- do.call(
    "cbind",
    structure(
      lapply(
        seq_along(labels),
        function(z) {
          # LABEL
          label <- labels[[z]]
          label_name <- names(labels)[z]
          # CLASSIFICATION
          if(is.factor(label)) {
            # REQUIRE RANN 
            cyto_require(
              "RANN",
              source = "CRAN"
            )
            # NEAREST NEIGHBOURS
            knn <- cyto_func_call(
              "RANN::nn2",
              list(
                data = train,
                query = test,
                k = k,
                treetype = "kd",
                searchtype = "standard"
              )
            )
            rm(list = c("test", "train"))
            # LABEL MEMBERSHIP PROBABILITIES
            label_mat <- matrix(
              label[knn$nn.idx],
              ncol = k
            )
            knn.prob <- switch(
              method,
              ## P(y_j | x_i) = sum(1/d(nn_i) * (y(nn_i) == y_j)) / sum(1/d(nn_i))
              'dist' = {
                sapply(
                  levels(label), 
                  function(cl, d, y) {
                    rowSums(1/d * (y == cl)) / rowSums(1/d)
                  }, 
                  d = pmax(knn$nn.dists, 1e-15), 
                  y = label_mat, 
                  simplify=FALSE, 
                  USE.NAMES=TRUE
                )
              },
              ## P(y_j | x_i) = sum(y(nn_i) == y_j) / k
              'vote' = {
                sapply(
                  levels(label), 
                  function(cl, y) {
                    rowSums(y == cl) / ncol(y)
                  },
                  y = label_mat, 
                  simplify=FALSE, 
                  USE.NAMES=TRUE
                )
              }
            )
            knn.prob <- as.matrix(do.call('cbind.data.frame', knn.prob))
            knn.prob <- sweep(knn.prob, 1, rowSums(knn.prob), "/")
            rm(list = c('knn', 'label_mat'))
            gc()
            # ASSIGN label
            res <- levels(label)[max.col(knn.prob, ties.method = "first")]
            res <- factor(res, levels(label))
            # REGRESSION
          } else {
            # REQUIRE FNN
            cyto_require(
              "FNN",
              source = "CRAN"
            )
            # KNN REGRESSION
            res <- cyto_func_call(
              "FNN::knn.reg",
              list(
                train = train,
                test = test,
                y = label,
                k = k,
                algorithm = "kd_tree"
              )
            )
            rm(list = c("test", "train"))
            # NEW VALUES
            res <- res$pred
          }
          return(res)
        }
      ),
      names = names(labels)
    )
  )
  
  # RETURN NEW LABELS | VALUES
  return(labels)
  
}