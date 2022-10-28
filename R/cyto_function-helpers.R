## CYTO_FUNC_NAME --------------------------------------------------------------

#' Get the name of a function as a charcter string
#'
#' @param FUN character string or function for which the function name should be
#'   extracted.
#'
#' @return a character string with the name of a function.
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @examples 
#' cyto_func_name("stats::density")
#'
#' @export
cyto_func_name <- function(FUN) {
  
  # FUNCTION -> CHARACTER
  if(cyto_class(FUN, "function")) {
    FUN <- as.character(
      substitute(
        FUN
      )
    )
  # PARSE CHARACTER FUNCTION
  } else {
    FUN <- unlist(strsplit(FUN, ":{2,3}"))
  }
  return(FUN[length(FUN)])
  
}

## CYTO_FUNC_MATCH -------------------------------------------------------------

#' Convert function name to valid function
#'
#' @param FUN a function or name of a function to be called. Namespaced
#'   functions are also supported (i.e. "CytoExploreR::cyto_plot").
#' @param descend logical to control whether to search past non-function
#'   objects.
#'
#' @return a valid function.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' cyto_func_match("uwot::umap")
#'
#' @export
cyto_func_match <- function(FUN, 
                            descend = TRUE){
  # NAMESPACED CHARACTER
  if(is.character(FUN)) {
    if(grepl(":{2,3}", FUN)) {
      FUN <- unlist(strsplit(FUN, ":{2,3}"))
      FUN <- get(
        FUN[2],
        envir = asNamespace(FUN[1]),
        mode = "function"
      )
    }
  }
  match.fun(FUN, descend = descend)
}

## CYTO_FUNC_ARGS --------------------------------------------------------------

#' Extract the formal arguments of a function
#'
#' @param FUN a function or name of a function to be called. Namespaced
#'   functions are also supported (i.e. "CytoExploreR::cyto_plot").
#' @param drop a vector of arguments to be removed from the returned vector of
#'   argument names.
#' @param logical indicating whether a named list of default arguments should be
#'   returned instead of a vector of argument names, set to FALSE by default.
#'
#' @return vector of argument names or a named list of argument defaults for
#'   \code{FUN} with \code{drop} arguments removed.
#'
#' @importFrom methods formalArgs
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' cyto_fun_args("cyto_plot")
#'
#' @export
cyto_func_args <- function(FUN, 
                           drop = NA,
                           defaults = FALSE) {
  
  # ARGUMENT NAMES + VALUES
  if(defaults) {
    args <- formals(
      cyto_func_match(
        FUN
      )
    )
    if(!.all_na(drop)) {
      args <- args[
        !names(args) %in% drop
      ]
    }
  # ARGUMENT NAMES
  } else {
    args <- formalArgs(
      cyto_func_match(
        FUN
      )
    )
    if(!.all_na(drop)) {
      args <- args[!args %in% drop]
    }
  }

  return(args)
}

## CYTO_FUNC_CALL --------------------------------------------------------------

#' Call a function on a set of arguments
#'
#' This function is used within CytoExploreR to call functions from external
#' packages that are not included as imports.
#'
#' @param FUN a function or name of a function to be called. Namespaced
#'   functions are also supported (i.e. "CytoExploreR::cyto_plot").
#' @param args a list of named arguments on which \code{FUN} should be called.
#'   Alternatively, these arguments can be supplied separately by names to
#'   \code{...} below.
#' @param ... named arguments that will be passed to the specified function when
#'   a named list of arguments is not supplied to \code{args}.
#'
#' @return output of the specified function when once called on the supplied
#'   arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' cyto_fun_call("cyto_channels", gs)
#'
#' @export
cyto_func_call <- function(FUN,
                           args = NULL,
                           ...) {
  
  # FUNCTION
  FUN <- cyto_func_match(FUN)
  
  # ARGUMENTS
  if(is.null(args)) {
    args <- list(...)
  }
  
  # CALL FUNCTION
  do.call(FUN, args)
  
}

## CYTO_FUNC_EXECUTE -----------------------------------------------------------

#' Call a function on its formal arguments
#'
#' \code{cyto_func_execute()} is similar to \code{cyto_fun_call()} with the
#' exception that \code{FUN} is only called on arguments returned by
#' \code{cyto_fun_args()}. This mean the function will execute successfully if
#' erroneous arguments are supplied but additional arguments won't be passed
#' through \code{...}.
#' 
#' @param FUN a function or name of a function to be called. Namespaced
#'   functions are also supported (i.e. "CytoExploreR::cyto_plot").
#' @param args a list of named arguments on which \code{FUN} should be called.
#'   Alternatively, these arguments can be supplied separately by names to
#'   \code{...} below.
#' @param drop a vector of arguments to be removed from the returned vector of
#'   argument names.
#' @param ... named arguments that will be passed to the specified function when
#'   a named list of arguments is not supplied to \code{args}.
#'
#' @return output of the specified function when once called on the supplied
#'   arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @examples 
#' #' library(CytoExploreRData)
#'
#' # Activation Gatingset
#' gs <- load_gs(system.file("extdata/Activation-GatingSet",
#'                           package = "CytoExploreRData"))
#'
#' cyto_fun_execute("cyto_channels", gs)
#'
#' @export
cyto_func_execute <- function(FUN,
                              args = NULL,
                              drop = NULL,
                              ...) {
  
  # FUNCTION
  FUN <- cyto_func_match(FUN)
  
  # ARGUMENTS
  if(is.null(args)) {
    args <- list(...)
  }
  
  # CALL FUNCTION
  do.call(
    FUN,
    args[names(args) %in% cyto_func_args(FUN, drop = drop)]
  )
  
}

## CYTO_FUNC_MAPPLY ------------------------------------------------------------

#' Apply a function iteratively over a list of arguments
#'
#' \code{cyto_func_mapply()} is a convenient wrapper around
#' \code{cyto_func_match()} and \code{mapply()} with the addition of an
#' \code{INDEX} argument to apply the function over a subset of the arguments.
#'
#' @param FUN name of the function to apply, passed to \code{cyto_func_match()}
#'   to obtain a valid function.
#' @param args a named list of arguments over which \code{FUN} is to be
#'   iteratively applied.
#' @param MoreArgs additional arguments passed to \code{FUN} that should not be
#'   iterated over.
#' @param SIMPLIFY passed to \code{mapply()} to determine the output format
#'   after each function call, set to TRUE by default.
#' @param USE.NAMES passed to \code{mapply()} to control whether names should be
#'   retained in the output of each functin call, set to TRUE by default.
#' @param INDEX integers to indicate the set of arguments in \code{args} over
#'   which \code{FUN} should be applied, set to NULL by default to iterate over
#'   all arguments.
#'
#' @return output of the specified function called iteratively over the list of
#'   arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples 
#' a <- c(1, 2, 3)
#' b <- c(1, 2, 3)
#' c <- c(1, 2, 3)
#' 
#' add_nums <- function(a, b, c){
#'   a + b + c
#' }
#'
#' cyto_func_mapply(
#'   "add_nums",
#'   list(a = a, b = b, c = c)
#' )
#' 
#' cyto_func_mapply(
#'   "add_nums",
#'   list(a = a, b = b, c = c),
#'   INDEX = 2
#' )
#'
#' @export
cyto_func_mapply <- function(FUN,
                             args,
                             MoreArgs = NULL,
                             SIMPLIFY = TRUE,
                             USE.NAMES = TRUE,
                             INDEX = NULL) {
  
  # PREPARE ARGUMENTS - INDEX
  if(!is.null(INDEX)) {
    args <- structure(
      lapply(
        args,
        function(arg) {
          arg[INDEX]
        }
      ),
      names = names(args)
    )
  }
  
  # CALL FUNCTION ON ARGUMENTS
  cyto_func_call(
    "mapply",
    c(
      "FUN" = cyto_func_match(FUN),
      args,
      "MoreArgs" = MoreArgs,
      "SIMPLIFY" = SIMPLIFY,
      "USE.NAMES" = USE.NAMES
    )
  )
  
}
