## CYTO_FILE_SEARCH ------------------------------------------------------------

#' Search for CSV file with required information in current directory
#'
#' \code{cyto_file_search()} performs a search in the current working directory
#' for a CSV file matching the supplied file name as well as the specified
#' column and row name requirements. If multiple files are located and
#' CytoExploreR is run in interactive mode, the user will be prompted to select
#' the file manually. Otherwise, \code{cyto_file_search()} will resort to
#' importing data from the first file name which matches all the requirements.
#' \code{cyto_file_search()} will automatically import the data from the
#' requested file and return it in the form of a \code{data.frame} or
#' \code{data.table} for downstream use.
#'
#' @param x regular expression to use when searching for files by name in
#'   specified directory.
#' @param dir location where the search should be performed.
#' @param rownames minimal requirement for row names in the file.
#' @param colnames minimal requirement for column names in the file.
#' @param ignore.cse logical passed to \code{grepl()} to indicate whether to
#'   ignore the case when matching the file name, set to TRUE by default.
#' @param data.table logical indicating whether the imported data should be
#'   returned as a \code{data.table} instead of a \code{data.frame}, set to
#'   FALSE by default.
#' @param type for internal use only to specify the type of data for which the
#'   search has been performed, simply used to provide informative messages when
#'   multiple files are located.
#' @param files optional vector of file names in the current working directory
#'   from which data should be imported. Supplying the file names manually will
#'   bypass the internal file search but still perform the same checks on the
#'   imported data.
#' @param ... additional named arguments indicating values to find in particular
#'   columns of the file (e.g. name = c("Activation_1.fcs")).
#'
#' @return NULL or a list named with the file name containing a
#'   \code{data.frame} or \code{data.table} of the data imported from the
#'   selected file.
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
                             data.table = FALSE,
                             type = NULL,
                             files = NULL,
                             ...) {
  
  # SEARCH FOR FILES
  if(is.null(files)) {
    files <- list.files(
      dir,
      pattern = x,
      recursive = TRUE,
      include.dirs = FALSE,
      full.names = TRUE,
      ignore.case = ignore.case
    )
  # FILES SUPPLIED
  } else {
    file_exists(
      file,
      error = TRUE
    )
  }
  
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
            z,
            data.table = data.table
          )
          # CHECK FILE COLUMN NAMES
          if(!is.null(colnames)) {
            if(!all(colnames %in% colnames(f))) {
              return(NULL)
            }
          }
          # CHECK FILE ROW NAMES
          if(!data.table & !is.null(rownames)) {
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
  # SINGLE FILE FOUND
  } else if(length(files) == 1) {
    message(
      paste0(
        "Importing ",
        if(is.null(type)) {
          "saved data"
        } else {
          type
        },
        " from ",
        names(files),
        "..."
      )
    )
    files <- files[1]
  # MULTIPLE FILES LOCATED
  } else {
    # INTERACTIVE FILE SELECTION
    if(interactive & cyto_option("CytoExploreR_interactive")) {
      message(
        paste0(
          "Multiple files located matching the ",
          type,
          " search criteria. ",
          "Which file would you like to import ",
          if(is.null(type)) {
            "data"
          } else {
            type
          },
          " from?"
        )
      )
      # FILE OPTIONS
      message(
        paste0(
          names(files),
          sep = "\n"
        )
      )
      # FILE SELECTION
      opt <- cyto_enquire(NULL)
      opt <- tryCatch(
        as.numeric(opt),
        warning = function(w){
          return(
            match(opt, names(pd))
          )
        }
      )
      files <- files[opt]
    # DEFAULT TO FIRST FILE
    } else {
      message(
        paste0(
          "Multiple files located matching the ",
          type,
          " search criteria. ",
          "Resorting to importing  ",
          if(is.null(type)) {
            "data"
          } else {
            type
          },
          " from ",
          names(files)[1],
          "..."
        )
      )
      files <- files[1]
    }
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
