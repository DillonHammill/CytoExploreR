## CYTO_STATS_COMPUTE ----------------------------------------------------------

#' Compute, export and save statistics
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param alias names of the populations in the \code{GatingHerarchy} or
#'   \code{GatingSet} for which statistics should be computed.
#' @param parent names of the parent populations in the \code{GatingHierarchy}
#'   or \code{GatingSet} required when computing frequency statistics. When
#'   \code{stat = "freq"} the frequency of each \code{alias} population will be
#'   computed as a proportion of each \code{parent} population.
#' @param channels names of the channels or markers for which statistics should
#'   be computed, set to all channels with marker assignments by default.
#' @param stat name of a statistic to compute, options include: \itemize{
#'   \item{\code{count} - number of events} \item{\code{freq/percent} -
#'   frequency of population(s) in parent population(s)} \item{\code{mean} -
#'   arithmetic mean} \item{\code{geomean} - geometric (graphical) mean computed
#'   as the inverse of the arithmetic mean on the transformed scale}
#'   \item{\code{median} - 50th quantile} \item{\code{mode} - inverse of mode on
#'   transformed scale computed using \code{smooth} and \code{bandwidth} to
#'   control the smoothness of the density distribution} \item{\code{sd} -
#'   standard deviation} \item{\code{rsd} - robust standard deviation (less
#'   influence of outliers)} \item{\code{cv} - coefficient of variation}
#'   \item{\code{rcv} - robust coefficient of variation (less influence of
#'   outliers)} \item{\code{quantile} - quantiles with probabilities supplied to
#'   \code{probs}} \item{\code{range} - minima and maxima} \item{\code{auc} -
#'   area under curve using a combination of
#'   \code{\link[stats:density]{density()}},
#'   \code{\link[stats:splinefun]{splinefun()}} and
#'   \code{\link[stats:integrate]{integrate()}}}} \code{stat} dispatches through
#'   \code{\link{cyto_apply}} so any custom function can be named through this
#'   argument as well. See \code{\link{cyto_apply}} for details.
#' @param trans an object of class \code{transformerList} containing the
#'   transformers used to transform the channels of the supplied data. The
#'   transformerList will be automatically extracted from the
#'   \code{GatingHierarchy} or \code{GatingSet}, so this argument is only
#'   required for \code{cytoframes} or \code{cytosets}. The transformerList is
#'   passed to \code{\link{cyto_apply}} which will make sure the data is
#'   appropriately transformed based on \code{inverse} prior to passing the data
#'   to the desired statistical function.
#' @param inverse logical passed to \code{\link{cyto_apply}} to indicate whether
#'   transformations applied to the data should be reversed prior to passing the
#'   data to the desired statistical function, set to TRUE by default.
#' @param round numeric indicating the number of decimal places to round the
#'   computed statistic, set to 2 decimal places by default.
#' @param format can be either \code{"wide"} or \code{"long"} to control the
#'   format of the returned statistics, set to \code{"wide"} by default.
#' @param tibble logical indicating whether the statistics should be returned as
#'   tibble instead of a data.frame, set to FALSE by default.
#' @param input passed to \code{\link{cyto_apply}} to control how the data is
#'   formatted prior to passing it to the statistical function, options include:
#'   \itemize{ \item{1 - \code{"cytoframe"}} \item{2 - \code{"matrix"}} \item{3
#'   - \code{"column"} or \code{"channel"}} \item{4 - \code{"row"} or
#'   \code{"cell"}}} Set to "matrix" by default.
#' @param smooth numeric smoothing parameter passed to \code{stats:density} when
#'   computing mode and area under the curve statistics, set to 1 by default.
#' @param bins number of bins to use for histograms, set to 256 by default.
#' @param details logical indicating whether to include the \code{cyto_details}
#'   in the output, set to TRUE by default.
#' @param markers logical indicating whether channels should be converted to
#'   markers where possible in the output, set to TRUE by default.
#' @param save_as name of a csv file to which the output should be saved, set to
#'   NULL by default to bypass this saving step.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{cytoset} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details.
#' @param merge_by vector of \code{cyto_details} column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be merged
#'   prior to computing statistics. Set this argument to "all" or NA to merge
#'   all samples. Set to \code{"name"} by default to compute statistics on each
#'   sample separately.
#' @param ... additional arguments passed to the desired statistical function.
#'
#' @return data.frame or tibble containing the computed statistics in the
#'   desired format and optionally a csv file containing the computed
#'   statistics.
#'
#' @importFrom tidyr gather
#' @importFrom dplyr all_of
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transformations
#' gs <- cyto_transform(gs)
#'
#' # Gating
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#'
#' # Compute statistics - median
#' cyto_stats_compute(gs,
#'   alias = "T Cells",
#'   channels = c("Alexa Fluor 488-A", "PE-A"),
#'   stat = "median"
#' )
#'
#' # Name csv to save results
#' tempfile <- paste0(tempdir(),
#' .Platform$file.sep,
#' "Population-Frequencies.csv")
#'
#' # Compute frequencies and save to tempfile
#' cyto_stats_compute(gs,
#'   alias = c("CD4 T Cells", "CD8 T Cells"),
#'   parent = c("Live Cells", "T Cells"),
#'   stat = "freq",
#'   save_as = tempfile
#' )
#'
#' @export
cyto_stats_compute <- function(x,
                               alias = NULL,
                               parent = NULL,
                               channels = NULL,
                               trans = NA,
                               stat = NULL,
                               inverse = TRUE,
                               round = 2,
                               format = "wide",
                               tibble = FALSE,
                               input = "matrix",
                               smooth = 1,
                               bins = 256,
                               details = TRUE,
                               markers = TRUE,
                               save_as = NULL,
                               select = NULL,
                               merge_by = "name",
                               ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # STATISTIC
  if(is.null(stat)) {
    stop("Supply the name of the statistical function to 'stat'.")
  }
  
  # DISPATCH
  stat <- .cyto_stat_dispatch(stat)
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- names(cyto_markers(x))
  }
  
  # TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # ALIAS
  if(is.null(alias)) {
    stop(
      "Supply the names of the population(s) to 'alias'!"
    )
  }
  
  # EXTRACT DATA ---------------------------------------------------------------
  
  # UNMERGED ALIAS
  if(all(merge_by %in% "name")) {
    alias <- cyto_data_extract(
      x,
      parent = alias,
      select = select,
      channels = channels,
      copy = TRUE
    )
  # MERGED ALIAS
  } else {
    alias <- structure(
      lapply(
        alias,
        function(z) {
          # MERGED CYTOSETS
          cs_list <- cyto_merge_by(
            x,
            select = select,
            parent = z,
            merge_by = merge_by,
            format = "cytoset",
            barcode = TRUE,
            channels = channels,
            copy = TRUE
          )
          # EXTRACT EXPERIMENT DETAILS - (MERGING GROUPS)
          pd <- do.call(
            "rbind",
            lapply(
              cs_list,
              "cyto_details"
            ) 
          )
          # MERGE CYTOSETS
          cs <- cytoset(
            structure(
              lapply(
                cs_list,
                function(z) {
                  z[[1]]
                }
              ),
              names = names(cs_list)
            )
          )
          # TRANSFER EXPERIMENT DETAILS - (MERGING GROUPS)
          cyto_details(cs) <- cbind(
            cyto_details(cs), 
            pd[, merge_by]
          )
          return(cs)
        }
      ),
      names = alias
    )
  }
  
  # ALIAS NAMES
  if(is.null(names(alias))) {
    names(alias) <- paste0("alias-", 1:length(alias))
  }
  
  # EXTRACT PARENT POPULATIONS
  if(cyto_class(x, "GatingSet")) {
    if(!is.null(parent)) {
      # UNMERGED PARENT
      if(all(merge_by %in% "name")) {
        parent <- cyto_data_extract(
          x,
          parent = parent,
          select = select,
          channels = channels,
          copy = TRUE
        )
      # MERGED PARENT
      } else {
        parent <- structure(
          lapply(
            parent,
            function(z) {
              # MERGED CYTOSETS
              cs_list <- cyto_merge_by(
                x,
                select = select,
                parent = z,
                merge_by = merge_by,
                format = "cytoset",
                barcode = TRUE,
                channels = channels,
                copy = TRUE
              )
              # EXTRACT EXPERIMENT DETAILS - (MERGING GROUPS)
              pd <- do.call(
                "rbind",
                lapply(
                  cs_list,
                  "cyto_details"
                ) 
              )
              # MERGE CYTOSETS
              cs <- cytoset(
                structure(
                  lapply(
                    cs_list,
                    function(z) {
                      z[[1]]
                    }
                  ),
                  names = names(cs_list)
                )
              )
              # TRANSFER EXPERIMENT DETAILS - (MERGING GROUPS)
              cyto_details(cs) <- cbind(
                cyto_details(cs), 
                pd[, merge_by]
              )
              return(cs)
            }
          ),
          names = parent
        )
      }
      # PARENT NAMES
      if(is.null(names(parent))) {
        names(parent) <- paste0("parent-", 1:length(parent))
      }
    }
    # PARENT POPULATIONS IN LIST
  } else {
    if(!all(LAPPLY(parent, "cyto_class", "flowSet"))) {
      stop("'parent' must be list of cytosets!")
    }
  }
  
  # COMPUTE STATISTICS ---------------------------------------------------------
  
  # CYTO_STAT FUNCTION
  if (ifelse(is.character(stat),
             grepl("cyto_stat_", stat),
             FALSE)) {
    # TRANSFORM
    if (grepl("count", stat) | grepl("geomean", stat)) {
      inverse <- FALSE
    } else {
      inverse <- TRUE
    }
    # STAT STRIPPED
    stat_strip <- gsub("cyto_stat_", "", stat)
    # COMPUTE STATISTICS
    if (any(LAPPLY(c(
      "count",
      "mean",
      "median",
      "sd",
      "rsd",
      "cv",
      "rcv",
      "quantile",
      "range",
      "bin"
    ), function(z) {
      grepl(paste0("^", z, "$"), stat_strip)
    }))) {
      res <- structure(
        lapply(
          alias, 
          function(z) {
            cyto_apply(
              z,
              stat,
              input = "matrix",
              copy = FALSE,
              trans = trans,
              inverse = inverse,
              round = round,
              simplify = TRUE,
              ...
            )
          }
        ), 
        names = names(alias)
      )
      # FREQUENCY
    } else if (grepl("^freq$", stat_strip)) {
      # PARENTAL COUNTS
      parent_counts <- structure(
        lapply(
          seq_along(parent), 
          function(z){
            res <- cyto_apply(
              parent[[z]],
              "nrow",
              input = "matrix",
              copy = FALSE
            )
            colnames(res) <- names(parent)[z]
            return(res)
          }
        ), 
        names = names(parent)
      )
      parent_counts <- do.call("cbind", parent_counts)
      # ALIAS COUNTS
      alias_counts <- structure(
        lapply(
          alias, 
          function(z){
            res <- cyto_apply(
              z,
              "nrow",
              input = "matrix",
              copy = FALSE
            )
            res <- res[, rep(1, ncol(parent_counts)), drop = FALSE]
            colnames(res) <- colnames(parent_counts)
            return(res)
          }
        ), names = names(alias)
      )
      # FREQUENCY
      res <- lapply(
        alias_counts, 
        function(z) {
          round(z/parent_counts*100, round)
        }
      )
      # GEOMEMTRIC MEAN
    } else if (grepl("^geomean$", stat_strip)) {
      res <- structure(
        lapply(
          alias, 
          function(z) {
            y <- cyto_data_extract(
              z, 
              format = "matrix"
            )[[1]]
            do.call(
              "rbind", 
              lapply(
                y, 
                function(q) {
                  LAPPLY(
                    unname(colnames(q)),
                    function(r) {
                      if (.all_na(trans)) {
                        cyto_stat_geomean(q[, r, drop = FALSE], 
                                          round = round)
                      } else if (!.all_na(trans)) {
                        if (r %in% names(trans)) {
                          .cyto_transform(
                            cyto_stat_mean(
                              q[, r, drop = FALSE], 
                              round = round
                            ),
                            channel = r,
                            trans = trans,
                            inverse = TRUE
                          )
                        } else {
                          cyto_stat_geomean(
                            q[, r, drop = FALSE],
                            round = round
                          )
                        }
                      }
                    }
                  )
                }
              )
            )
          }
        ), 
        names = names(alias)
      )
      # MODE
    } else if (grepl("^mode$", stat_strip)) {
      res <- structure(
        lapply(
          alias, 
          function(z) {
          # DENSITY - CURRENT SCALE
          d <- suppressWarnings(
            cyto_apply(
              z,
              "cyto_stat_density",
              input = "matrix",
              smooth = smooth,
              bins = bins,
              limits = range(z[[1]], type = "instrument"),
              simplify = TRUE,
              copy = FALSE,
              ...
            )
          )
          # MODE
          cnt <- 0
          m <- apply(
            d, 
            2, 
            function(y) {
              # CHANNEL INDEX
              cnt <<- cnt + 1
              round(
                LAPPLY(
                  names(y), 
                  function(w){
                    if(.all_na(y[[w]])) {
                      return(NA)
                    }
                    if(inverse == TRUE & colnames(d)[cnt] %in% names(trans)) {
                      .cyto_transform(
                        y[[w]]$x[y[[w]]$y == max(y[[w]]$y)],
                        channel = colnames(d)[cnt],
                        trans = trans,
                        inverse = TRUE
                      )
                    } else {
                      y[[w]]$x[y[[w]]$y == max(y[[w]]$y)]
                    }
                  }
                ), 
                round
              )
            }
          )
          rownames(m) <- rownames(d)
          return(m)
          }
        ), 
        names = names(alias)
      )
      # AUC
    } else if (grepl("^auc$", stat_strip)) {
      res <- structure(
        lapply(
          alias, 
          function(z){
            cyto_apply(
              z,
              stat,
              input = "matrix",
              trans = trans,
              inverse = inverse,
              round = round,
              simplify = TRUE,
              smooth = smooth,
              bins = bins,
              limits = range(z[[1]], type = "instrument"),
              ...
            )
          }
        )
      )
    }
    # CUSTOM STATISTIC FUNCTION
  } else {
    res <- structure(
      lapply(
        alias, 
        function(z){
          round(
            cyto_apply(
              z,
              stat,
              input = input,
              trans = trans,
              inverse = inverse,
              simplify = TRUE,
              ...
            ),
            round
          )
        }
      ), 
      names = names(alias)
    )
  }
  
  # FORMAT ---------------------------------------------------------------------
  
  # STATISTIC NAME
  if (is.function(stat)) {
    stat <- deparse(substitute(stat))
    # UNNAMED FUNCTION
    if (length(stat) > 1) {
      stat <- "FUN"
    }
  }
  
  # POPULATION & EXPERIMENT DETAILS
  res <- do.call(
    "rbind",
    c(
      structure(
        lapply(seq_along(res), function(z){
          # MARKERS
          if(markers) {
            colnames(res[[z]]) <- LAPPLY(colnames(res[[z]]), function(w) {
              if (w %in% cyto_channels(x)) {
                return(cyto_markers_extract(x, w))
              } else {
                return(w)
              }
            })
          }
          # ABSORB ROWNAMES
          if(grepl("\\|", rownames(res[[z]])[1])) {
            new_cols <- do.call("rbind", 
                                strsplit(rownames(res[[z]]), 
                                         "\\|"))[, -1, drop = FALSE] # name
            if(ncol(new_cols) == 1) {
              colnames(new_cols) <- gsub("cyto_stat_", "", stat)
            } else {
              colnames(new_cols) <- paste0(gsub("cyto_stat_", "", stat), "-",
                                           1:ncol(new_cols))
            }
            res[[z]] <- data.frame(new_cols, 
                                   res[[z]],
                                   check.names = FALSE,
                                   stringsAsFactors = FALSE)
          }
          # EXPERIMENT DETAILS
          if(details) {
            res[[z]] <- data.frame(
              cyto_details(alias[[z]], drop = TRUE),
              "alias" = rep(names(res)[z], nrow(res[[z]])),
              res[[z]],
              check.names = FALSE,
              stringsAsFactors = FALSE)
          } else {
            res[[z]] <- data.frame(
              "name" = cyto_names(alias[[z]]),
              "alias" = rep(names(res)[z], nrow(res[[z]])),
              res[[z]],
              check.names = FALSE,
              stringsAsFactors = FALSE)
          }
          # CANNOT HAVE ROWNAMES IN MERGED DATAFRAME
          rownames(res[[z]]) <- NULL
          return(res[[z]])
        }), names = names(res)
      ), 
      list("make.row.names" = FALSE,
           stringsAsFactors = FALSE)
    )
  )

  # FORMAT
  if(format == "long") {
    # LONG FORMAT REQUIRES MULTIPLE COLUMNS
    cols <- colnames(res)[!colnames(res) %in% c("name", 
                                                stat, 
                                                gsub("cyto_stat_", "", stat),
                                                "FUN", 
                                                colnames(cyto_details(x)),
                                                "alias"), 
                          drop = FALSE]
    # KEY - CHANNEL
    if(all(cols %in% cyto_channels(x))) {
      key <- "channel"
    } else if(all(cols %in% cyto_markers(x))) {
      key <- "marker"
    } else if(all(cols %in% names(parent))) {
      key <- "parent"
    } else {
      key <- "input"
    }
    
    # LONG FORMAT
    if (length(cols) > 1) {
      res <- gather(res,
                    key = {{ key }},
                    value = "value",
                    all_of(cols)
      )
    }
  }

  # SAVE
  if(!is.null(save_as)) {
    write_to_csv(res, save_as)
  }
  
  # TIBBLE
  if(tibble) {
    class(res) <- c("tbl_df", "tbl", "data.frame")
  }
  
  # COMPUTED STATISTICS
  return(res)
  
}
