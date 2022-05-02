## CYTO_FILE_SEARCH ------------------------------------------------------------

#' Search for CSV file in current directory
#'
#' @param x regular expression to use when searching for files by name in
#'   specified directory.
#' @param dir location where the search should be performed.
#' @param rownames minimal requirement for row names in the file.
#' @param colnames minimal requirement for column names in the file.
#' @param ignore.cse logical passed to \code{grepl()} to indicate whether to
#'   ignore the case when matching the file name, set to TRUE by default.
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
                if(!all(cols[[w]] %in% f[, w])) {
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

## CYTO_FILE_NAME --------------------------------------------------------------

#' Create a unique name for a file
#'
#' @param x desired name for a file optionally including the extension, if the
#'   extension is not manually supplied to \code{ext}.
#' @param ext the file extension prefixed by \code{"."}.
#' @param ... not in use.
#' 
#' @return a unique file name.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' cyto_file_name("Experiment-Details.csv")
#'
#' @export
cyto_file_name <- function(x,
                           ext = NULL,
                           ...) {
  
  # FILE EXTENSION
  if(!nzchar(file_ext(x))) {
    if(is.null(ext)) {
      stop(
        paste0(
          x,
          " does not contain a valid file extension!"
        )
      )
    }
    x <- paste0(x, ".", file_ext(x))
  }
  
  # FILE DOES NOT EXIST
  if(!file_exists(x)) {
    return(x)
  }
  
  # EXTENSION
  ext <- file_ext(x)
  
  # CHECK SUFFIX
  id <- gsub(
    paste0(
      ".*-([0-9]+)\\.",
      ext, 
      "$"
    ), 
    "\\1", 
    x
  )
  
  # NEW SUFFIX
  if(id == x) {
    id = 0
  } else {
    id <- as.numeric(id)
  }
  
  # INCREMENT SUFFIX
  while(file_exists(x)) {
    id <- id + 1
    if(id == 1) {
      x <- paste0(
        gsub(
          paste0(
            "(.*)\\.",
            ext,
            "$"
          ),
          "\\1",
          x
        ),
        "-",
        id,
        ".",
        ext
      )
    } else {
      x <- paste0(
        gsub(
          paste0(
            "(.*)-[0-9]+\\.",
            ext,
            "$"
          ),
          "\\1",
          x
        ),
        "-",
        id,
        ".",
        ext
      )
    }

  }

  return(x)
  
}
