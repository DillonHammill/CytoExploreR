## EMPTY CHARACTER STRINGS -----------------------------------------------------

#' Check if vector contains only empty chracter strings
#'
#' @param x vector.
#' 
#' @return TRUE/FALSE
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.empty <- function(x){

  if(.all_na(x)){
    return(FALSE)
  }else if(is.character(x)){
    if(all(nchar(trimws(x)) == 0)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
  
}

## MATCH_IND -------------------------------------------------------------------
  
#' Indices for matches excluding NA
#' @noRd
match_ind <- function(x, y, ...){
  ind <- match(x, y, ...)
  ind[!is.na(ind)]
}

## ALL NA ----------------------------------------------------------------------

#' Check all elements of vector are NA
#' 
#' @param x vector.
#' 
#' @return TRUE/FALSE
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.all_na <- function(x){
  if(is.null(x)){
    return(FALSE)
  }else{
    return(all(suppressWarnings(is.na(unlist(x)))))
  }
}

## ARGUMENT LIST ---------------------------------------------------------------

#' Pull down arguments from environment into list
#' 
#' Replace empty elements with empty characters "".
#' 
#' @return alist object containing arguments of parent function environment.
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.args_list <- function(...){
  
  # Pull down ... arguments
  dot_args <- list(...)
  
  # Get arguments from parental environment
  args <- as.list(parent.frame())
  
  # Combine ... args with args
  if(length(dot_args) != 0){
    args <- c(args, dot_args)
  }
  
  # Remove duplicate args
  args <- args[which(!duplicated(names(args)))]
  
  # Replace any elements with class "name" with ""
  lapply(names(args), function(x){
    if(all(class(args[[x]]) == "name")){
      args[[x]] <<- ""
    }
  })
  
  # Convert to alist
  class(args) <- "alist"
  
  # Return argument list
  return(args)
  
}

## ARGUMENT UPDATE -------------------------------------------------------------

#' Update arguments of function using a named list of arguments
#' 
#' @param x named list of arguments to assign to function environment.
#' 
#' @return update arguments in function environment.
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.args_update <- function(x){
  
  lapply(seq(1,length(x)), function(z){
    assign(names(x)[z], 
           x[[z]], envir = parent.frame(n = 3))
  })
  
}

## FILE WD CHECK ---------------------------------------------------------------

#' Check if a file exists in the current working directory
#'
#' @param name filename including file extension to be checked.
#'
#' @return TRUE/FALSE if file exists in the current working directory.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' file_wd_check("gatingTemplate.csv")
#' 
#' @export
file_wd_check <- function(name) {
  if (length(which(list.files() == name)) != 0) {
    
    # File exists in working directory
    return(TRUE)
  } else if (length(which(list.files() == name)) == 0) {
    
    # File does not exist in working directory
    return(FALSE)
  }
}

## ROUND -----------------------------------------------------------------------

#' @noRd
.round <- function(x, k = 2){
  trimws(format(round(x, k), nsmall = k))
} 

## LAPPLY ----------------------------------------------------------------------

#' Automatically flatten lapply results
#' @noRd
LAPPLY <- function(...){
  unlist(lapply(...))
}

## PAR -------------------------------------------------------------------------

#' Extract graphical parameters as list
#' @importFrom graphics par
#' @noRd
.par <- function(x){
  if(length(x) == 1){
    pars <- list(par(x))
    names(pars) <- x
  }else{
    pars <- par(x)
  }
  return(pars)
}

# SUPPRESS WARNINGS & MESSAGES -------------------------------------------------

#' Suppress warnings and messages
#' @noRd
.suppress_all_messages <- function(...){
  suppressWarnings(suppressMessages(...))
}

## SEURAT HELPER ---------------------------------------------------------------

# Helper function from seurat to get program paths.

# Get program paths in a system-agnostic way
#
# @param progs A vector of program names
# @param error Throw an error if any programs are not found
# @param add.exe Add '.exe' extension to program names that don't have it
#
# @return A named vector of program paths; missing programs are returned as
# \code{NA} if \code{error = FALSE}
#
#' @importFrom tools file_ext
#  @noRd
SysExec <- function(
  progs,
  error = ifelse(test = length(x = progs) == 1, yes = TRUE, no = FALSE),
  add.exe = .Platform$OS.type == 'windows'
) {
  cmd <- ifelse(
    test = .Platform$OS.type == 'windows',
    yes = 'where.exe',
    no = 'which'
  )
  if (add.exe) {
    missing.exe <- file_ext(x = progs) != 'exe'
    progs[missing.exe] <- paste0(progs[missing.exe], '.exe')
  }
  paths <- LAPPLY(
    X = progs,
    FUN = function(x) {
      return(tryCatch(
        expr = system2(command = cmd, args = x, stdout = TRUE)[1],
        warning = function(...) {
          return(NA_character_)
        }
      ))
    }
  )
  if (error && any(is.na(x = paths))) {
    stop(
      "Could not find the following programs: ",
      paste(names(x = paths[is.na(x = paths)]), collapse = ', '),
      call. = FALSE
    )
  }
  return(paths)
}
