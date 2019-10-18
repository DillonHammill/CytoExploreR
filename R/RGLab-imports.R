## Import Unexported Functions from RGLab Suite --------------------------------

## flowCore Imports ------------------------------------------------------------

#' flowCore .estimateLogicle
#' @importFrom utils getFromNamespace
#' @noRd
.estimateLogicle <- getFromNamespace(".estimateLogicle", "flowCore")

## openCyto Imports ------------------------------------------------------------

#' openCyto .preprocess_csv
#' @importFrom utils getFromNamespace
#' @noRd
.preprocess_csv <- getFromNamespace(".preprocess_csv", "openCyto")

#' openCyto alias
#' @importFrom utils getFromNamespace
#' @noRd
.alias <- getFromNamespace("alias", "openCyto")

#' openCyto groupBy
#' @importFrom utils getFromNamespace
#' @noRd
.groupBy <- getFromNamespace("groupBy", "openCyto")

#' openCyto isCollapse
#' @importFrom utils getFromNamespace
#' @noRd
.isCollapse <- getFromNamespace("isCollapse", "openCyto")

#' openCyto ppMethod
#' @importFrom utils getFromNamespace
#' @noRd
.ppMethod <- getFromNamespace("ppMethod", "openCyto")
