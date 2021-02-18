## CYTOEXPLORER LOGO -----------------------------------------------------------

#' Path to CytoExploreR logo 
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au) 
#' 
#' @noRd
CytoExploreR_logo <- function(){
  paste0(
    "https://raw.githubusercontent.com/DillonHammill/CytoExploreR",
    "/master/man/figures/logo.png"
  )
}

## EMPTY CHARACTER STRINGS -----------------------------------------------------

#' Check if vector contains only empty character strings
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
    if(all(!nzchar(trimws(x)))){
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
#' @param x parameters to extract, extract all if missing.
#' @param set logical indicating whether only parameters that can be set should
#'   be returned.
#' @importFrom graphics par
#' @noRd
.par <- function(x,
                 set = TRUE){
  if (missing(x)) {
    pars <- par()
    if(set) {
      pars <- pars[!names(pars) %in% c("cin",
                                       "cra",
                                       "csi",
                                       "cxy",
                                       "din",
                                       "page")]
    }
  } else if(length(x) == 1) {
    pars <- list(par(x))
    names(pars) <- x
  } else {
    pars <- par(x)
  }
  return(pars)
}

## SUPPRESS WARNINGS & MESSAGES ------------------------------------------------

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

## EXTRACT NUMERIC FROM STRING -------------------------------------------------

#' Extract numeric from string
#' @noRd
str_num <- function(x) {
  as.numeric(gsub("[^0-9.-]+", "", as.character(x)))
}
