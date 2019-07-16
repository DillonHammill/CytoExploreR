# EMPTY CHARACTER STRINGS ------------------------------------------------------

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

  if(is.character(x)){
    if(all(nchar(trimws(x)) == 0)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
  
}

# ALL NA -----------------------------------------------------------------------

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
  all(suppressWarnings(is.na(x)))
}

# ARGUMENT LIST ----------------------------------------------------------------

#' Pull down arguments from environment into list
#' 
#' Replace empty elements with empty characters "".
#' 
#' @return alist object containing arguments of parent function environment.
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.args_list <- function(){
  
  # Get arguments from parental environment
  args <- as.list(parent.frame())
  
  # Replace any elements with class "name" with ""
  lapply(names(args), function(x){
    if(all(class(args[[x]]) == "name")){
      args[[x]] <<- ""
    }
  })
  
  # Convert to alist
  class(args) <- "alist"
  
  return(args)
  
}

# ARGUMENT UPDATE --------------------------------------------------------------

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

# FILE WD CHECK ----------------------------------------------------------------

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

# ROUND ------------------------------------------------------------------------

#' @noRd
.round <- function(x, k) trimws(format(round(x, k), nsmall=k))
