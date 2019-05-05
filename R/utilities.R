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
  
  empty <- all(unlist(lapply(x, function(y){
    
    if(is.character(y) &
       nchar(trimws(y) == 0)){
      return(TRUE)
    }else{
      return(FALSE)
    }
    
  })))
  
  return(empty)
  
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
  all(is.na(x))
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
#' @export
args_list <- function(){
  
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