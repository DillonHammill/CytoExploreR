## FILE HELPERS ----------------------------------------------------------------

# Below are a series of functions to assist in the manipulation, checking,
# sorting or extracting elements from file names.

## FILE WD CHECK ---------------------------------------------------------------

#' Check if a file exists
#'
#' @param x filename including file extension to be checked.
#' @param error logical indicating whether an informative error should be thrown
#'   if the file does not exist, set to FALSE by default.
#'
#' @return TRUE/FALSE if file exists.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' file_exists("gatingTemplate.csv")
#' @noRd
file_exists <- function(x,
                        error = FALSE) {
  ind <- which(!file.exists(x))
  if(length(ind) > 0) {
    if(error) {
      stop(
        paste(paste(x[ind], collapse = " & "),
              "may not exist or lack the required permissions.")
      )
    }
    return(FALSE)
  }
  return(TRUE)
}

## FILE_EXT --------------------------------------------------------------------

#' Get file extensions from file names
#' @noRd
file_ext <- function(x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

## FILE_EXT_REMOVE -------------------------------------------------------------

#' Remove file extension from file names
#' @noRd
file_ext_remove <- function(x,
                            compression = FALSE) {
  if (compression) {
    x <- sub("[.](gz|bz2|xz)$", "", x)
  }
  sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)
}

## FILE_EXT_APPEND -------------------------------------------------------------

#' Append file name with an extension
#'
#' @param x vector of filenames.
#' @param ext vector of extensions to append.
#'
#' @importFrom tools file_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
file_ext_append <- function(x,
                            ext = "csv") {
  
  # REPEAT EXTENSION
  ext <- rep(ext, length(x))
  
  # ADD EXTENSIONS TO FILE NAMES WITHOUT EXTENSIONS
  LAPPLY(seq_along(x), function(z) {
    # PREPARE EXTENSION
    if (!grepl(".", ext[z])) {
      ext[z] <- paste0(".", ext[z])
    }
    # APPEND EXTENSION
    if (.empty(file_ext(x[z]))) {
      paste0(x[z], ext)
    } else {
      x[z]
    }
  })
}

## FILE_SORT -------------------------------------------------------------------

# Modified version of mixedsort from gtools.

#' Sort file names
#' @noRd
file_sort <- function(x,
                      decreasing = FALSE,
                      na.last = TRUE,
                      blank.last = FALSE) {
  
  # - Split each each character string into an vector of strings and
  #   numbers
  # - Separately rank numbers and strings
  # - Combine orders so that strings follow numbers
  
  if (length(x) < 1) {
    return(NULL)
  } else if (length(x) == 1) {
    return(1)
  }
  
  if (!is.character(x)) {
    return(order(x, decreasing = decreasing, na.last = na.last))
  }
  
  delim <- "\\$\\@\\$"

  regex <- paste0(
    "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])",
    "(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)",
    "(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))") # uses PERL syntax
  numeric <- function(x) as.numeric(x)
  
  nonnumeric <- function(x) {
    ifelse(is.na(numeric(x)), toupper(x), NA)
  }
  
  x <- as.character(x)
  
  which.nas <- which(is.na(x))
  which.blanks <- which(x == "")
  
  ####
  # - Convert each character string into an vector containing single
  #   character and numeric values.
  ####
  
  # find and mark numbers in the form of +1.23e+45.67
  delimited <- gsub(regex,
                    paste(delim, "\\1", delim, sep = ""),
                    x,
                    perl = TRUE
  )
  
  # separate out numbers
  step1 <- strsplit(delimited, delim)
  
  # remove empty elements
  step1 <- lapply(step1, function(x) x[x > ""])
  
  # create numeric version of data
  suppressWarnings(step1.numeric <- lapply(step1, numeric))
  
  # create non-numeric version of data
  suppressWarnings(step1.character <- lapply(step1, nonnumeric))
  
  # now transpose so that 1st vector contains 1st element from each
  # original string
  maxelem <- max(sapply(step1, length))
  
  step1.numeric.t <- lapply(
    1:maxelem,
    function(i) {
      sapply(
        step1.numeric,
        function(x) x[i]
      )
    }
  )
  
  step1.character.t <- lapply(
    1:maxelem,
    function(i) {
      sapply(
        step1.character,
        function(x) x[i]
      )
    }
  )
  
  # now order them
  rank.numeric <- sapply(step1.numeric.t, rank)
  rank.character <- sapply(
    step1.character.t,
    function(x) as.numeric(factor(x))
  )
  
  # and merge
  rank.numeric[!is.na(rank.character)] <- 0 # mask off string values
  
  rank.character <- t(
    t(rank.character) +
      apply(matrix(rank.numeric), 2, max, na.rm = TRUE)
  )
  
  rank.overall <- ifelse(is.na(rank.character), rank.numeric, rank.character)
  
  order.frame <- as.data.frame(rank.overall)
  if (length(which.nas) > 0) {
    if (is.na(na.last)) {
      order.frame[which.nas, ] <- NA
    } else if (na.last) {
      order.frame[which.nas, ] <- Inf
    } else {
      order.frame[which.nas, ] <- -Inf
    }
  }
  
  if (length(which.blanks) > 0) {
    if (is.na(blank.last)) {
      order.frame[which.blanks, ] <- NA
    } else if (blank.last) {
      order.frame[which.blanks, ] <- 1e99
    } else {
      order.frame[which.blanks, ] <- -1e99
    }
  }
  
  order.frame <- as.list(order.frame)
  order.frame$decreasing <- decreasing
  order.frame$na.last <- NA
  
  retval <- do.call("order", order.frame)
  
  return(x[retval])
}

## READ_FROM_CSV ---------------------------------------------------------------

#' Use data.table to read in a csv file
#' @importFrom data.table fread
#' @noRd
read_from_csv <- function(x,
                          data.table = FALSE,
                          ...){
  x <- file_ext_append(x, ".csv")
  ind <- which(!file.exists(x))
  if(length(ind) > 0){
    stop(
      paste(
        paste(x[ind], collapse = " & "),
        "do not exist or lack the required permissions."
      )
    )
  }
  dt <- suppressMessages(
    fread(x, 
          data.table = TRUE,
          ...)
  )
  if(data.table == FALSE) {
    dt <- as.data.frame(dt)
    if(colnames(dt)[1] == "V1") {
      rownames(dt) <- dt[, 1]
      dt <- dt[, -1, drop = FALSE]
    }
  }
  return(dt)
}

## WRITE_TO_CSV ----------------------------------------------------------------

#' Use data.table to write data to a csv file
#' @importFrom data.table fwrite
#' @noRd
write_to_csv <- function(x,
                         file = NULL,
                         ...){
  
  if(is.null(file)){
    stop("Supply a name for the csv file to 'file'.")
  }
  file <- file_ext_append(file, ".csv")
  suppressMessages(
    fwrite(x,
           file = file,
           ...)
  )
  
}
