## CYTO_SPECTRA_COMPARE --------------------------------------------------------

#' Compute cosine similarity or unmixing error hotspot matrix
#'
#' @param x object of class \code{"matrix"} or list of matrix objects with each
#'   row representing a uniquely named spectrum and each column representing a
#'   detector or channel.
#' @param select vector of row to select from each matrix by name, set to NULL
#'   by default to compare all rows.
#' @param type indicates whether to compute the \code{"cosine"} similarity,
#'   spectral \code{"purity"}, \code{"hotspot"} matrix or panel
#'   \code{"complexity"}, set to \code{"cosine"} by default. Option
#'   \code{"hotspot"} is only available when a single matrix has been supplied
#'   to \code{x}.
#' @param save_as name of a CSV file to which the similarity matrix should be
#'   saved, set to NULL by default to bypass saving.
#' @param heatmap logical indicating whether the computed similarity scores
#'   should be displayed in a heatmap, set to to TRUE by default.
#' @param title text to display in the header of the heatmap when \code{heatmap
#'   = TRUE}, defaults to either \code{"Cosine Similarity Matrix"},
#'   \code{"Spectral Purity Matrix"} or \code{"Unmixing Error Hotspot"} matrix
#'   depending on \code{type}.
#' @param ... additional arguments passed to \code{HeatmapR::heat_map()} to
#'   customise the displayed heatmap.
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
                                 type = "cosine",
                                 save_as = NULL,
                                 heatmap = TRUE,
                                 title = NULL,
                                 ...) {
  
  # MATRIX
  if(cyto_class(x, "matrix", FALSE)) {
    # SELECT ROWS
    if(!is.null(select)) {
      idx <- ulapply(
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
    # SIMILARITY MATRIX REQUIRED
    if(grepl("^cos|^sim|^h", type, ignore.case = TRUE)) {
      # COSINE SIMILARITY MATRIX
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
      # HOTSPOT MATRIX
      if(grepl("^h", type)) {
        cs <- sqrt(abs(solve(cs)))
      }
    # SPECTRAL PURITY
    } else if(grepl("^p", type)) {
      cs <- rbind(
        "purity" = row_purity_cpp(x)
      )
    # CONDITION NUMBER
    } else if(grepl("^com|^cond", type, ignore.case = TRUE)) {
      cs <- kappa(x)
    } else {
      stop(
        "Unsupported 'type'!"
      )
    }
  # LIST
  } else if(cyto_class(x, "list", "TRUE")) {
    # LIST OF MATRICES
    if(!all(ulapply(x, cyto_class, expect = "matrix", class = "FALSE"))) {
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
            idx <- ulapply(
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
      ulapply(
        x,
        rownames
      )
    )
    keep <- names(cnt)[cnt == length(x)]
    # COSINE SIMILARITY
    if(grepl("^cos|^sim", type, ignore.case = TRUE)) {
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
    # SPECTRAL PURITY
    } else if(grepl("^p", type, ignore.case = TRUE)) {
      cs <- do.call(
        "rbind",
        lapply(
          x,
          function(z) {
            row_purity_cpp(z)
          }
        )
      )
      rownames(cs) <- names(x)
    # COMPLEXITY 
    } else if(grepl("^com|^cond", type, ignore.case = TRUE)) {
      cs <- ulapply(
        x,
        "kappa"
      )
    # UNSUPPORTED TYPE
    } else {
      stop(
        "Only cosine similarity and spectral purity are supported for lists!"
      )
    }
  # UNSUPPORTED OBJECT
  } else {
    stop(
      "'x' must be either a matrix or list of matrices!"
    )
  }

  # PLOT HEATMAPS
  if(heatmap & !grep("^com|^cond", type, ignore.case = TRUE)) {
    HeatmapR::heat_map(
      cs,
      title = if(is.null(title)) {
        if(grepl("^cos|^sim", type)) {
          "Cosine Similarity Matrix"
        } else if(grepl("^p", type)) {
          "Spectral Purity Matrix"
        } else {
          "Unmixing Error Hotspot Matrix"
        }
      } else {
        title
      },
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
