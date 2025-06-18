## CYTO_SPECTRA_COMPARE --------------------------------------------------------

#' Compute spectral similarity between spectra
#'
#' @param x object of class \code{"matrix"} or list of matrix objects with each
#'   row representing a uniquely named spectrum and each column representing a
#'   detector or channel.
#' @param select vector of row to select from each matrix by name, set to NULL
#'   by default to compare all rows.
#' @param save_as name of a CSV file to which the similarity matrix should be
#'   saved, set to NULL by default to bypass saving.
#' @param heatmap logical indicating whether the computed similarity scores
#'   should be displayed in a heatmap, set to to TRUE by default.
#' @param ... additional arguments passed to \code{HeatmapR::heat_map()} to
#'   cusomise the displayed heatmap.
#'
#' @return computed similarity matrix, heatmap and exported CSV file if
#'   requested.
#'   
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @importFrom HeatmapR heat_map
#' 
#' @seealso \code{\link{cyto_plot_spectra}}
#' @seealso \code{\link{cyto_unmix_compute}}
#' @seealso \code{\link{cyto_spillover_compute}}
#'
#' @export
cyto_spectra_compare <- function(x,
                                 select = NULL,
                                 save_as = NULL,
                                 heatmap = TRUE,
                                 ...) {
  
  # MATRIX
  if(cyto_class(x, "matrix", FALSE)) {
    # SELECT ROWS
    if(!is.null(select)) {
      idx <- LAPPLY(
        select,
        function(z) {
          id <- match(z, rownames(x))
          if(.all_na(id)) {
            id <- grep(
              z, 
              rownames(x),
              ignore.case = TRUE
            )
          }
          return(id)
        }
      )
      x <- x[idx, , drop = FALSE]
    }
    # CHECK ROWS
    if(!nrow(x) > 1) {
      stop(
        "Multiple rows required in 'x' to compute similarity scores!"
      )
    }
    # CHECK ROWNAMES
    if(is.null(rownames(x))) {
      stop(
        "Rownames must be included in 'x' to identify each fluorchrome or event!"
      )
    }
    # SIMILARITY MATRIX
    cs <- diag(
      1, 
      nrow = nrow(x),
      ncol = nrow(x)
    )
    colnames(cs) <- rownames(x)
    rownames(cs) <- colnames(cs)
    for(i in 1:nrow(x)) {
      for(j in i:nrow(x)) {
        cs[j, i] <- cs[i, j] <- .cosine(
          as.numeric(x[i, ]),
          as.numeric(x[j, ])
        )
      }
    }
  # LIST
  } else if(cyto_class(x, "list", "TRUE")) {
    # LIST OF MATRICES
    if(!all(LAPPLY(x, cyto_class, expect = "matrix", class = "FALSE"))) {
      stop("'x' must be a list of matrices to compute similarity scores!")
    }
    # CHECK NAMES
    if(is.null(names(x))) {
      names(x) <- paste0("matrix-", seq_along(x))
    }
    # SELECT ROWS
    if(!is.null(select)) {
      x <- structure(
        lapply(
          x,
          function(z) {
            idx <- LAPPLY(
              select,
              function(w) {
                id <- match(w, rownames(z))
                if(.all_na(id)) {
                  id <- grep(
                    w, 
                    rownames(z),
                    ignore.case = TRUE
                  )
                }
                return(id)
              }
            )
            # FILTER MATRIX
            z <- z[idx, , drop = FALSE]
            # CHECK ROWS
            if(!nrow(z) > 1) {
              stop(
                "Multiple rows required in 'x' to compute similarity scores!"
              )
            }
            # CHECK ROWNAMES
            if(is.null(rownames(z))) {
              stop(
                "Rownames must be included in 'x' to identify each fluorchrome or event!"
              )
            }
            return(z)
          }
        ),
        names = names(x)
      )
    }
    # SELECT ROWS THAT OCCUR IN ALL MATRICES
    cnt <- table(
      LAPPLY(
        x,
        rownames
      )
    )
    keep <- names(cnt)[cnt == length(x)]
    # SIMILARITY MATRIX - COMPUTE RELATIVE TO FIRST MATRIX
    cs <- matrix(
      1,
      ncol = length(keep),
      nrow = length(x)
    )
    colnames(cs) <- keep
    rownames(cs) <- names(x)
    for(i in keep) {
      for(j in seq_along(x)) {
        cs[j, match(i, colnames(cs))] <- .cosine(
          as.numeric(x[[1]][match(i, rownames(x[[1]])), ]),
          as.numeric(x[[j]][match(i, rownames(x[[j]])), ])
        )
      }
    }
  # UNSUPPORTED OBJECT
  } else {
    stop(
      "'x' must be either a matrix or list of matrices!"
    )
  }

  # PLOT HEATMAPS
  if(heatmap) {
    HeatmapR::heat_map(
      cs,
      ...
    )
  }
  
  # RETURN COSINE SIMILARITY MATRICES
  return(cs)

}

#' Internal function to compute cosine similarity of two vectors
#'
#' @param x n-dimensional numeric vector.
#' @param y n-dimensional numeric vector.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#'
#' @examples
#' \dontrun{
#' .cosine(
#'   c(1, 2, 3, 4),
#'   c(1, 4, 5, 7)
#' )
#' }
#'
#' @noRd
.cosine <- function(x,
                    y) {
  if (length(x) != length(y)) {
    stop(
      "Input vectors must have the same length to compute cosine similarity!"
    )
  }
  
  cos <- crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
  
  return(cos[1, 1])
}
