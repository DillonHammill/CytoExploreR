## Import Unexported Functions from RGLab Suite --------------------------------

## flowCore Imports ------------------------------------------------------------

#' flowCore .estimateLogicle
#' @importFrom utils getFromNamespace
#' @noRd
.estimateLogicle <- getFromNamespace(".estimateLogicle", "flowCore")

## openCyto Imports ------------------------------------------------------------

#' openCyto .preprocess.csv
#' @importFrom utils getFromNamespace
#' @noRd
.preprocess_csv <- getFromNamespace(".preprocess_csv", "openCyto")

#' openCyto .argDeparser
#' @importFrom utils getFromNamespace
#' @noRd
.argDeparser <- getFromNamespace(".argDeparser", "openCyto")