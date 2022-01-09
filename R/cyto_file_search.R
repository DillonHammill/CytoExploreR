## CYTO_FILE_SEARCH ------------------------------------------------------------

#' Search for CSV file in current directory
#'
#' @param x regular expression to use when searching for files by name in
#'   specified directory.
#' @param dir location where the search should be performed.
#' @param rownames minimal requirement for row names in the file.
#' @param colnames minimal requirement for column names in the file.
#' @param ... additional named arguments indicating values to find in particular
#'   columns of the file (e.g. name = c("Activation_1.fcs")).
#'
#' @return NULL or data read in from file matching search criteria.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{ 
#'   # search for experiment details of samples
#'   cyto_file_search( 
#'     "Experiment-Details.csv$", 
#'     rownames = c("Activation_1.fcs", "Activation_2.fcs")
#'   ) 
#' }
#'
#' @export
cyto_file_search <- function(x,
                             dir = ".",
                             colnames = NULL,
                             rownames = NULL,
                             ignore.case = TRUE,
                             ...) {
  
  # FILES IN DIRECTORY
  files <- list.files(
    dir,
    pattern = x,
    recursive = TRUE,
    include.dirs = FALSE,
    full.names = TRUE,
    ignore.case = ignore.case
  )
  
  # FILES FOUND
  if(length(files) > 0) {
    # COLUMNS TO SEARCH
    cols <- list(...)
    
    # SEARCH FILES
    files <- structure(
      lapply(
        files,
        function(z) {
          # READ FILE
          f <- read_from_csv(
            z
          )
          # CHECK FILE COLUMN NAMES
          if(!is.null(colnames)) {
            if(!all(colnames %in% colnames(f))) {
              return(NULL)
            }
          }
          # CHECK FILE ROW NAMES
          if(!is.null(rownames)) {
            if(!all(rownames %in% rownames(f))) {
              return(NULL)
            }
          }
          # CHECK COLUMNS
          if(length(cols) > 0) {
            stop <- FALSE
            lapply(
              names(cols),
              function(w) {
                if(!all(cols[[w]] %in% f[, cols[[w]]])) {
                  stop <<- TRUE
                }
              }
            )
            if(stop) {
              return(NULL)
            }
          }
          return(f)
        }
      ),
      names = files
    )
    
    # RETURN DATA
    files <- files[!LAPPLY(files, "is.null")]
  }
  
  # NO FILES FOUND
  if(length(files) == 0) {
    files <- NULL
  }
  
  # READ IN DATA IN LIST
  return(files)
  
}
