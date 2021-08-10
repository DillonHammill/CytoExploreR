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
      FUN <- get(FUN[2],
                 envir = asNamespace(FUN[1]),
                 mode = "function")
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
#'
#' @return vector of argument names for \code{FUN} with \code{drop} arguments
#'   removed.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples 
#' cyto_fun_args("cyto_plot")
#'
#' @export
cyto_func_args <- function(FUN, 
                           drop = NA) {
  # FORMAL ARGUMENTS
  args <- formalArgs(FUN)
  if(!.all_na(drop)) {
    args <- args[!args %in% drop]
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
