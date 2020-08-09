## CYTO_STATS_COMPUTE ----------------------------------------------------------

#' Compute, export and save statistics
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
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
#'   \item{\code{median} - 50th quantile} \item{\code{mode} - mode computed
#'   using \code{smooth} and \code{bandwidth} to control the smoothness of the
#'   density distribution} \item{\code{sd} - standard deviation}
#'   \item{\code{rsd} - robust standard deviation (less influence of outliers)}
#'   \item{\code{cv} - coefficient of variation} \item{\code{rcv} - robust
#'   coefficient of variation (less influence of outliers)}
#'   \item{\code{quantile} - quantiles with probabilities supplied to
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
#' @param gate a \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate} object to apply to each \code{cytoframe} prior to
#'   computing the desired statistic. This argument has been included for
#'   backwards compatibility only and is only valid for \code{cytoframe} or
#'   \code{cytoset} objects.
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
#'   computing mode and area under the curve statistics, set to 0.6 by default.
#' @param bandwidth numeric passed to \code{stats:density} to set the bandwidth
#'   when computing mode or area under the curve statistics, set to NULL by
#'   default. If the bandwidth is not supplied or NULL, a bandwidth will be
#'   estimated based on all samples supplied to \code{cyto_stats_compute()}.
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
#' @param ... additional arguments passed to the desired statistical function.
#'
#' @return data.frame or tibble containing the computed statistics in the
#'   desired format and optionally a csv file containing the computed
#'   statistics.
#'
#' @importFrom tidyr gather
#' @importFrom flowCore Subset
#' @importFrom methods is
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
#' @name cyto_stats_compute
NULL

#' @noRd
#' @export
cyto_stats_compute <- function(x,
                               ...) {
  UseMethod("cyto_stats_compute")
}

#' @rdname cyto_stats_compute
#' @export
cyto_stats_compute.GatingSet <- function(x,
                                         alias = NULL,
                                         parent = NULL,
                                         channels = NULL,
                                         trans = NULL,
                                         stat = NULL,
                                         gate = NULL,
                                         inverse = TRUE,
                                         round = 2,
                                         format = "wide",
                                         tibble = FALSE,
                                         input = "matrix",
                                         smooth = 0.6,
                                         bandwidth = NULL,
                                         details = TRUE,
                                         markers = TRUE,
                                         save_as = NULL,
                                         select = NULL,
                                         ...) {

  # STATISTIC MISSING
  if (is.null(stat)) {
    stop("Supply the name of a statistical function to 'stat'.")
  }
  
  # DISPATCH
  stat <- cyto_stat_dispatch(stat)
  
  # SELECT
  if (!is.null(select)) {
    x <- cyto_select(x, select)
  }

  # TRANSFORMERS
  trans <- cyto_transformer_extract(x)

  # FREQUENCY - GATINGHIERARCHY METHOD
  if (switch(is.character(stat),
             "TRUE" = grepl("cyto_stat_freq", stat),
             "FALSE" = FALSE)) {
    res <- lapply(seq_along(x), function(z) {
      cyto_stats_compute(x[[z]],
        alias = alias,
        parent = parent,
        channels = channels,
        trans = trans,
        stat = stat,
        inverse = inverse,
        round = round,
        format = format,
        tibble = FALSE,
        input = "matrix",
        smooth = smooth,
        bandwidth = bandwidth,
        details = details,
        markers = markers,
        save_as = NULL,
        ...
      )
    })
  # OTHER - CYTOSET METHOD (BANDWIDTH)
  } else {
    # FLOWSET METHOD
    res <- lapply(alias, function(z){
      st <- cyto_stats_compute(cyto_extract(x, z),
                               channels = channels,
                               trans = trans,
                               stat = stat,
                               inverse = inverse,
                               round = round,
                               format = format,
                               tibble = FALSE,
                               input = input,
                               smooth = smooth,
                               bandwidth = bandwidth,
                               details = FALSE,
                               markers = markers,
                               save_as = NULL,
                               ...)
      # POPULATION
      if(is.function(stat)) {
        cols <- c("name", "FUN")
      }else {
        cols <- c("name", stat)
      }
      st <- cbind(st[, colnames(st) %in% cols, drop = FALSE],
                  "alias" = rep(z, nrow(st)),
                  st[, !colnames(st) %in% cols, drop = FALSE])
      # DETAILS
      if (details == TRUE) {
        pd <- cyto_details(x,
                           convert = TRUE,
                           drop = TRUE)
        new_cols <- do.call("rbind", lapply(st$name, function(z) {
          pd[pd$name == z, ]
        }))
        st <- cbind(new_cols, st[, -match("name", colnames(st)), drop = FALSE])
      }
      return(st)
    })
  }
  res <- do.call("rbind", res)
  rownames(res) <- NULL

  # TIBBLE
  if (tibble) {
    class(res) <- c("tbl_df", "tbl", "data.frame")
  }

  # SAVE
  if (!is.null(save_as)) {
    write_to_csv(res, save_as)
  }

  # STATISTICS
  return(res)
}

#' @rdname cyto_stats_compute
#' @export
cyto_stats_compute.GatingHierarchy <- function(x,
                                               alias = NULL,
                                               parent = NULL,
                                               channels = NULL,
                                               trans = NULL,
                                               stat = NULL,
                                               inverse = TRUE,
                                               round = 2,
                                               format = "wide",
                                               tibble = FALSE,
                                               input = "matrix",
                                               smooth = 0.6,
                                               bandwidth = NULL,
                                               details = TRUE,
                                               markers = TRUE,
                                               save_as = NULL,
                                               ...) {

  # STATISTIC MISSING
  if (is.null(stat)) {
    stop("Supply the name of a statistical function to 'stat'.")
  }

  # DISPATCH
  stat <- cyto_stat_dispatch(stat)

  # TRANSFORMERS
  trans <- cyto_transformer_extract(x)

  # ALIAS
  if (is.null(alias)) {
    stop("Supply the name of the population to 'alias'.")
  }

  # ALIAS POPULATIONS
  alias_frames <- structure(lapply(alias, function(z) {
    cyto_extract(x, z)
  }),
  names = alias
  )

  # PARENT POPULATIONS
  if (switch(is.character(stat),
    "TRUE" = grepl("cyto_stat_freq", stat),
    "FALSE" = FALSE
  )) {
    if (is.null(parent)) {
      message(
        paste(
          "Calculating frequency of 'root' as no 'parent' population(s) were",
          "specified."
        )
      )
      parent <- "root"
    }
    parent_frames <- structure(lapply(
      parent,
      function(z) {
        cyto_extract(x, z)
      }
    ),
    names = parent
    )
  }

  # FREQUENCY
  if (switch(is.character(stat),
    "TRUE" = grepl("cyto_stat_freq", stat),
    "FALSE" = FALSE
  )) {

    # PARENT COUNTS
    parent_counts <- cyto_apply(parent_frames,
      "cyto_stat_count",
      input = 2,
      copy = FALSE
    )
    # ALIAS COUNTS
    alias_counts <- cyto_apply(alias_frames,
      "cyto_stat_count",
      input = 2,
      copy = FALSE
    )

    # FREQUENCIES
    res <- lapply(parent_counts, function(z) {
      st <- alias_counts / z * 100
      if(any(is.nan(st))) {
        st[is.nan(st)] <- 0
      }
      round(st, round)
    })
    res <- do.call("cbind", res)
    colnames(res) <- parent
    # POPULATION
    res <- cbind(
      "name" = rep(cyto_names(x), nrow(res)),
      "alias" = alias,
      res
    )
    # DATA.FRAME
    res <- as.data.frame(res)
    # FORMAT
    if (format == "long") {
      res <- pivot_longer(res,
        names_to = "parent",
        values_to = "freq"
      )
    }

    # STATISTICS
  } else {
    res <- lapply(names(alias_frames), function(z) {
      st <- cyto_stats_compute(alias_frames[[z]],
                               channels = channels,
                               trans = trans,
                               stat = stat,
                               inverse = inverse,
                               format = format,
                               input = input,
                               smooth = smooth,
                               bandwidth = bandwidth,
                               markers = markers,
                               round = round,
                               tibble = FALSE,
                               save_as = NULL,
                               ...)
      if(is.function(stat)) {
        cols <- c("name", "FUN")
      }else {
        cols <- c("name", stat)
      }
      cbind(st[, colnames(st) %in% cols, drop = FALSE],
            "alias" = rep(z, nrow(st)),
            st[, !colnames(st) %in% cols, drop = FALSE])
    })
    res <- do.call("rbind", res)
  }

  # DETAILS
  if (details == TRUE) {
    pd <- cyto_details(x,
                       convert = TRUE,
                       drop = TRUE)
    new_cols <- do.call("rbind", lapply(res$name, function(z) {
      pd[pd$name == z, ]
    }))
    res <- cbind(new_cols, res[, -match("name", colnames(res)), drop = FALSE])
  }

  # TIBBLE
  if (tibble) {
    class(res) <- c("tbl_df", "tbl", "data.frame")
  }

  # SAVE
  if (!is.null(save_as)) {
    write_to_csv(res, save_as)
  }

  # STATISTICS
  return(res)
}

#' @rdname cyto_stats_compute
#' @export
cyto_stats_compute.flowSet <- function(x,
                                       channels = NULL,
                                       trans = NA,
                                       stat = NULL,
                                       gate = NA,
                                       inverse = TRUE,
                                       round = 2,
                                       format = "wide",
                                       tibble = FALSE,
                                       input = "matrix",
                                       smooth = 0.6,
                                       bandwidth = NULL,
                                       details = TRUE,
                                       markers = TRUE,
                                       save_as = NULL,
                                       select = NULL,
                                       ...) {

  # STATISTIC MISSING
  if (is.null(stat)) {
    stop("Supply the name of a statistical function to 'stat'.")
  }

  # DISPATCH
  stat <- cyto_stat_dispatch(stat)

  # SAME BANDWIDTH
  if (switch(is.character(stat),
    "TRUE" = grepl("auc", stat) | grepl("mode", stat),
    "FALSE" = FALSE
  )) {
    if (is.null(bandwidth)) {
      bandwidth <- cyto_apply(x,
        "cyto_stat_bandwidth",
        input = "matrix",
        channels = channels,
        trans = trans,
        inverse = inverse,
        simplify = TRUE
      )
      bandwidth <- colMeans(bandwidth)
    }
  }

  # SELECT
  if (!is.null(select)) {
    x <- cyto_select(x, select)
  }

  # STATISTIC - CANNOT USE CYTO_APPLY - NEED CHANNELS ETC
  res <- lapply(seq_along(x), function(z){
    cyto_stats_compute(x[[z]],
                       channels = channels,
                       trans = trans,
                       stat = stat,
                       inverse = inverse,
                       format = format,
                       input = input,
                       smooth = smooth,
                       bandwidth = bandwidth,
                       markers = markers,
                       round = round,
                       tibble = FALSE,
                       save_as = NULL,
                       gate = gate,
                       ...)
  })
  res <- do.call("rbind", res)

  # DETAILS
  if (details == TRUE) {
    pd <- cyto_details(x,
                       convert = TRUE,
                       drop = TRUE)
    new_cols <- do.call("rbind", lapply(res$name, function(z) {
      pd[pd$name == z, ]
    }))
    res <- cbind(new_cols, res[, -match("name", colnames(res)), drop = FALSE])
  }

  # TIBBLE
  if (tibble) {
    class(res) <- c("tbl_df", "tbl", "data.frame")
  }

  # SAVE
  if (!is.null(save_as)) {
    write_to_csv(res, save_as)
  }

  # STATISTICS
  return(res)
}

#' @rdname cyto_stats_compute
#' @export
cyto_stats_compute.flowFrame <- function(x,
                                         channels = NULL,
                                         trans = NA,
                                         stat = NULL,
                                         gate = NA,
                                         inverse = TRUE,
                                         round = 2,
                                         format = "wide",
                                         input = "matrix",
                                         smooth = 0.6,
                                         bandwidth = NULL,
                                         markers = TRUE,
                                         tibble = FALSE,
                                         save_as = NULL,
                                         ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (is.null(channels)) {
    channels <- names(cyto_markers(x))
  }

  # STATISTIC MISSING
  if (is.null(stat)) {
    stop("Supply the name of a statistical function to 'stat'.")
  }

  # DISPATCH
  stat <- cyto_stat_dispatch(stat)

  # GATING ---------------------------------------------------------------------

  # GATE - EITHER GATE OBJECT OR CYTOFRAME
  if (!.all_na(gate)) {
    if (!any(c(
      is(gate, "rectangleGate"),
      is(gate, "polygonGate"),
      is(gate, "ellipsoidGate"),
      is(gate, "flowFrame"),
      is(gate, "flowSet")
    ))) {
      stop(
        paste(
          "Only rectangleGate, polygonGate and ellipsoidGate objects are",
          "supported."
        )
      )
    }
    if (!is(gate, "flowFrame") & !is(gate, "flowSet")) {
      if (switch(is.character(stat),
        "TRUE" = !grepl("cyto_stat_freq", stat),
        "FALSE" = FALSE
      )) {
        x <- Subset(x, gate)
      }
    }
  }

  # COMPUTE STATISTICS ---------------------------------------------------------
  
  # CYTO_STAT FUNCTION
  if (switch(is.character(stat),
             "TRUE" = grepl("cyto_stat_", stat),
             "FALSE" = FALSE)) {
    # COPY
    if (grepl("count", stat)) {
      copy <- FALSE
    } else {
      copy <- TRUE
    }
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
      "quant",
      "range"
    ), function(z) {
      grepl(paste0("^", z, "$"), stat_strip)
    }))) {
      res <- cyto_apply(x,
        stat,
        input = "matrix",
        channels = channels,
        copy = copy,
        trans = trans,
        inverse = inverse,
        round = round,
        simplify = TRUE,
        ...
      )
      # FREQUENCY
    } else if (grepl("^freq$", stat_strip)) {
      # GATE REQUIRED
      if (is(gate, "flowFrame")) {
        y <- gate
      } else if (is(gate, "flowSet")) {
        y <- gate[[cyto_names(x)]]
      } else {
        y <- Subset(x, gate)
      }
      res <- round(
        cyto_apply(y,
          "cyto_stat_count",
          input = "matrix",
          copy = FALSE,
          simplify = TRUE
        ) /
          cyto_apply(x,
            "cyto_stat_count",
            input = "matrix",
            copy = FALSE,
            simplify = TRUE
          ) * 100, round
      )
      if(any(is.nan(res))) {
        is.nan(res) <- 0
      }
      # GEOMEMTRIC MEAN
    } else if (grepl("^geomean$", stat_strip)) {
      # EXTRACT DATA
      raw <- cyto_extract(x,
        raw = TRUE,
        channels = channels
      )[[1]]
      res <- LAPPLY(colnames(raw), function(z) {
        if (.all_na(trans)) {
          cyto_stat_geomean(raw[, z, drop = FALSE], 
                            round = round)
        } else if (!.all_na(trans)) {
          if (z %in% names(trans)) {
            .cyto_transform(cyto_stat_mean(raw[, z, drop = FALSE], 
                                           round = round),
              channel = z,
              trans = trans,
              inverse = TRUE
            )
          } else {
            cyto_stat_geomean(raw[, z, drop = FALSE],
                              round = round)
          }
        }
      })
      res <- matrix(res,
        ncol = ncol(raw),
        dimnames = list(
          cyto_names(x),
          colnames(raw)
        )
      )
      # MODE
    } else if (grepl("^mode$", stat_strip)) {
      # DENSITY - TRANSFORMED SCALE
      res <- cyto_apply(x,
        "cyto_stat_density",
        input = "matrix",
        channels = channels,
        copy = TRUE,
        trans = trans,
        inverse = inverse,
        smooth = smooth,
        bandwidth = bandwidth,
        simplify = TRUE,
        ...
      )
      col_names <- colnames(res)
      # MODE
      res <- apply(res, 1, function(z){
        LAPPLY(names(z), function(w){
          if(.all_na(z[[w]])) {
            return(NA)
          }
          if (inverse == TRUE) {
            .cyto_transform(z[[w]]$x[z[[w]]$y == max(z[[w]]$y)],
                            channel = w,
                            trans = trans,
                            inverse = TRUE
            )
          } else {
            z[[w]]$x[z[[w]]$y == max(z[[w]]$y)]
          }
        })
      })
      res <- round(res, round)
      res <- t(res)
      colnames(res) <- col_names
      # AUC
    } else if (grepl("^auc$", stat_strip)) {
      res <- cyto_apply(x,
        stat,
        input = input,
        channels = channels,
        trans = trans,
        inverse = inverse,
        round = round,
        copy = TRUE,
        simplify = TRUE,
        smooth = smooth,
        bandwidth = bandwidth,
        ...
      )
    }
    # CUSTOM STATISTIC FUNCTION
  } else {
    res <- cyto_apply(x,
      stat,
      input = input,
      channels = channels,
      trans = trans,
      inverse = inverse,
      copy = TRUE,
      simplify = TRUE,
      ...
    )
    res <- round(res, round)
  }
  
  # FORMAT ---------------------------------------------------------------------

  # MARKERS
  if (markers == TRUE) {
    colnames(res) <- LAPPLY(colnames(res), function(z) {
      if (z %in% cyto_channels(x)) {
        return(cyto_markers_extract(x, z))
      } else {
        return(z)
      }
    })
  }

  # DATAFRAME
  res <- data.frame(res,
    check.names = FALSE
  )
  
  # STAT NAME
  if (is.function(stat)) {
    stat <- deparse(substitute(stat))
    # UNNAMED FUNCTION
    if (length(stat) > 1) {
      stat <- "FUN"
    }
  }

  # STRIP CYTO_STAT
  if (grepl("cyto_stat_", stat)) {
    stat <- gsub("cyto_stat_", "", stat)
  }

  # ABSORB ROWNAMES
  new_cols <- do.call(
    "rbind",
    strsplit(rownames(res), "\\|")
  )
  colnames(new_cols) <- c("name", stat)[seq_len(ncol(new_cols))]
  res <- cbind(new_cols, res)
  rownames(res) <- NULL
  
  # FORMAT
  if (format == "long") {
    # LONG FORMAT REQUIRES MULTIPLE COLUMNS
    cols <- colnames(res)[!colnames(res) %in% c("name", stat, "FUN"), 
                          drop = FALSE]
    if (length(cols)) {
      res <- gather(res,
        key = "input",
        value = "output",
        cols
      )
    }
  }

  # TIBBLE
  if (tibble) {
    class(res) <- c("tbl_df", "tbl", "data.frame")
  }

  # SAVE
  if (!is.null(save_as)) {
    write_to_csv(res, save_as)
  }

  # STATISTICS
  return(res)
}
