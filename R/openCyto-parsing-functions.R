#' deparse a list(named) of expression into a string 
#' inverse function of argParser
#' 
#' @param args gatingTemplate arguments as named list to be deparsed (e.g. list(gate = gateobj)).
#' @param split logical.
#' 
#' @return deparsed arguments for storage in gatingTemplate csv file
#' 
#' @noRd
.argDeparser <- function(args, split = TRUE){
  if(split){
    args <- sapply(names(args), function(argn){
      
      argv <- deparse(args[[argn]])
      argv <- gsub("\"", "'", argv) #restore dquote to squote
      argv <- paste(argv, collapse = "")
      paste(argn, argv, sep = " = ")
    })
    
    
    paste(args, collapse = ", ")  
  }else
    as.character(args[[1]])
  
}
