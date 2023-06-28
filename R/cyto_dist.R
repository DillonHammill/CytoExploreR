## CYTO_DIST -------------------------------------------------------------------

#' Compute distance matrices between samples or distributions
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population upon which comparisons should be made,
#'   set to the \code{"root"} node by default.
#' @param select vector or list of sample selection criteria passed to
#'   \code{cyto_select()} to extract the samples for which the specified
#'   \code{type} statistic should be computed, set to NULL by default to use all
#'   samples in the supplied \code{cytoset} or \code{GatingSet}.
#' @param alias names of the populations to use when comparing the compositions
#'   of samples using the \code{Aitchison} distance.
#' @param channels names of the channel(s) or marker(s) for which comparisons on
#'   the \code{parent} population should be made.
#' @param type indicates the type of metric to use when comparing samples. For
#'   compositional comparisons when \code{alias} is supplied, \code{type} is set
#'   to \code{"Aitchison"} to compute the Aitchison distances. If
#'   \code{channels} are specified comparison of distributions across samples
#'   can be made using either the default \code{"Wasserstein"} distance,
#'   \code{"Jensen-Shannon"} divergence or \code{"Kolmogorov-Smirnov"} test. If
#'   \code{alias} and \code{channels} are not specified, comparisons of the
#'   \code{parent} populations can be made either using the default \code{"kNN"}
#'   or alternatively \code{"SOM"} entropies. The proportion of non-overlapping
#'   events between distributions can be computed using \code{type = "overlap"}.
#' @param merge_by a vector of experiment variables by which samples should be
#'   merged prior to computing the desired statistic, set to \code{"name"} by
#'   default to compare individual samples to one another.
#' @param scale optional argument to scale each channel prior to computing
#'   distances, options include \code{"global"}, \code{"range"}, \code{"mean"},
#'   \code{"median"} or \code{"zscore"}. Set to \code{"range"} by default,
#'   scaling can be turned off by setting this argument to FALSE. Global scaling
#'   preserves the structure of the data as expected on the linear scale without
#'   having to apply inverse data transformations to the entire dataset.
#' @param trans an object of class \code{transformerList} containing the
#'   transformer definitions for transformations already applied to the supplied
#'   cytoset.
#' @param inverse logical indicating whether inverse data transformations should
#'   be applied to the data from the \code{parent} population prior to computing
#'   the desired statistic, set to FALSE by default.
#' @param events indicates the number or proportion of events to extract from
#'   the \code{parent} population in each sample prior to computing the
#'   statistic specified by \code{type}. Setting \code{events = NA} will extract
#'   the minimum number of events across groups from each group.
#' @param k number of nearest neighbours to use when constructing the kNN graph
#'   for comparison of \code{kNN} entropies, set to 30 by default.
#' @param grid dimensions to use for the SOM grid when computing SOM entropies,
#'   set to \code{c(14, 14)}.
#' @param heatmap logical indicating whether to plot a heatmap for each computed
#'   distance matrix, set to TRUE by default.
#' @param markers logical indicating whether the named list of matrices should
#'   be labelled with the marker names when available, set to TRUE by default.
#'   Otherwise, the list of matrices will be labelled with the names of the
#'   channels.
#' @param ... additional arguments passed to \code{HeatmapR::heat_map()}.
#'
#' @return a distance matrix for compositional or sample comparisons, or a list
#'   of distance matrices for each channel for distribution comparisons.
#'
#' @importFrom stats ks.test as.dist
#' @importFrom HeatmapR heat_map
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_dist <- function(x,
                      parent = "root",
                      select = NULL,
                      alias = NULL,
                      channels = NULL,
                      type = "aitchison",
                      merge_by = "name",
                      scale = "range",
                      trans = NA,
                      inverse = FALSE,
                      events = NA,
                      k = 30,
                      grid = c(14,14),
                      heatmap = TRUE,
                      markers = TRUE,
                      ...) {
  
  # TODO: ADD HELLINGER DISTANCE
  # TODO: BREGMAN DIVERGENCE
  # TODO: ADD SUPPORT FOR HEATMAPS
  
  # COMPOSITION: AITCHISON
  # DISTRIBUTIONS: WASSERSTEIN | JENSEN-SHANNON | KOLMOGOROV-SMIRNOV
  # SAMPLES: KNN | SOM | MS
  
  # CHECKS ---------------------------------------------------------------------
  
  # CYTOSET | GATINGSET
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    stop(
      "cyto_compare() only accepts objects of class cytoset or GatingSet!"
    )
  }
  
  # CHANNELS REQUIRED
  if(is.null(channels)) {
    if(is.null(alias)) {
      stop(
        paste0(
          "Supply the name(s) of the channel(s) or marker(s) over which ",
          "the distance calcultaions should be made."
        )
      )
    }
  # CHANNELS SUPPLIED
  } else {
    channels  <- cyto_channels_extract(
      x,
      channels = channels
    )
  }
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(
      x,
      select
    )
  }
  
  # SPLIT INTO GROUPS
  x_list <- cyto_group_by(
    x,
    group_by = merge_by
  )
  
  # AITCHISON
  if(!is.null(alias)) {
    type <- "Aitchison"
  }
  
  # COMPUTE DISTANCES ----------------------------------------------------------
  
  # AITCHISON DISTANCE
  if(!is.null(alias)) {
    # GATINGSET ONLY
    if(!cyto_class(x, "GatingSet")) {
      stop(
        paste0(
          "'x' must be a GatingSet object to compare sample compostitions ",
          "using Aitchison distance!"
        )
      )
    }
    # REQUIRE ROBCOMPOSITIONS - ALTERNATIVE CODA.BASE
    cyto_require(
      "robCompositions",
      source = "CRAN"
    )
    # COMPUTE COUNTS PER GROUP - GROUP AS ROW & ALIAS AS COLUMNS
    cnts <- do.call(
      "rbind",
      structure(
        lapply(
          x_list,
          function(z) {
            # LABEL COUNTS
            label_cnts <- rep(0, length(alias))
            names(label_cnts) <- alias
            # LABELS
            lapply(
              seq_along(z),
              function(w) {
                res <- table(
                  cyto_gate_indices(
                    z,
                    parent = parent,
                    select = w,
                    nodes = alias,
                    labels = TRUE
                  )[[1]]
                )
                label_cnts[names(res)] <<- label_cnts[names(res)] + res
              }
            )
            return(label_cnts)
          }
        ),
        names = names(x_list)
      )
    )
    # AITCHISON DISTNACE REQUIRES NON-ZERO DATA - ADD VALUE BELOW DETECTION
    cnts <- cnts + 0.1
    # AITCHISON DISTANCE MATRIX
    res <- matrix(
      NA,
      ncol = length(x_list),
      nrow = length(x_list),
      dimnames = list(
        names(x_list),
        names(x_list)
      )
    )
    # PROGRESS BAR
    pb <- cyto_progress(
      label = "cyto_dist()",
      total = sum(seq(length(x_list), 1))
    )
    # COMPUTE AITCHSION DISTANCES
    for(i in 1:length(x_list)) {
      for(j in i:length(x_list)) {
        res[i, j] <- res[j, i] <- cyto_func_call(
          "robCompositions::aDist",
          list(
            "x" = cnts[c(i, j), ]
          )
        )
        # INCREMENT PROGRESS BAR
        cyto_progress(pb)
      }
    }
    # CONVERT TO DISTNACE MATRIX
    res <- as.matrix(
      as.dist(
        res
      )
    )
  # DISTRIBUTIONS | PARENT
  } else {
    # EVENTS
    if(.all_na(events)) {
      events <- min(
        LAPPLY(
          x_list,
          function(z) {
            cnts <- cyto_stats_compute(
              z,
              parent = parent,
              stat = "count"
            )
            sum(
              cnts[, ncol(cnts)]
            )
          }
        )
      )
    }
    # EXTRACT DATA - LIST OF CYTOSETS - > LINEAR IF REQUIRED
    x_list <- structure(
      lapply(
        x_list,
        function(z) {
          cyto_data_extract(
            z,
            parent = parent,
            channels = channels,
            format = "cytoset",
            coerce = TRUE,
            trans = trans,
            inverse = inverse,
            events = events,
            copy = inverse,
            seed = 2022
          )[[1]]
        }
      ),
      names = names(x_list)
    )
    # COMPARE PARENT - KNN | SOM
    if(grepl("^(KNN|SOM|MS)$", type, ignore.case = TRUE)) {
      # PROGRESS BAR
      pb <- cyto_progress(
        label = "cyto_dist()",
        total = sum(seq(length(x_list), 1))
      )
      # KNN ENTROPIES
      if(grepl("^kNN", type, ignore.case = TRUE)) {
        # ENTROPY DISTANCES
        res <- matrix(
          0,
          nrow = length(x_list),
          ncol = length(x_list),
          dimnames = list(
            names(x_list),
            names(x_list)
          )
        )
        # TODO: PROGRESS BAR?
        # PAIRWISE KNN ENTROPIES _ WASSERSTEIN DISTANCE
        for(i in seq_along(x_list)) {
          for(j in i:length(x_list)) {
            if(i == j) {
              cyto_progress(pb)
              next
            }
            # COMBINED DATA TO BUILD KNN GRAPH
            x_mat <- do.call(
              "rbind",
              lapply(
                c(i, j),
                function(z) {
                  m <- cyto_data_extract(
                    x_list[[z]],
                    format = "matrix"
                  )[[1]][[1]]
                  m <- cbind(
                    m,
                    "**group**" = z
                  )
                  return(m)
                }
              )
            )
            # SCALE DATA
            if(!scale %in% FALSE) {
              x_mat <- cyto_stat_scale(
                x_mat,
                type = scale
              )
            }
            # BUILD COMBINED KNN GRAPH
            knn <- .cyto_knn(
              x_mat,
              k = k + 1
            )
            # KNN INDICES
            knn <- knn$nn.idx[, -1]
            # MAP KNN INDICES TO GROUP INDICES
            knn <- apply(
              knn,
              1,
              function(z) {
                # GROUP COUNTS
                tbl <- table(x_mat[, "**group**"][z])
                cnts <- rep(0, length(x_list))
                names(cnts) <- names(x_list)
                cnts[names(x_list)[as.numeric(names(tbl))]] <- tbl
                # COMPUTE ENTROPY
                -sum(
                  sapply(
                    cnts,
                    function(w) {
                      if(!is.finite(w)) {
                        return(0)
                      } else {
                        if(w == 0) {
                          return(0)
                        } else {
                          return(w*log2(w))
                        }
                      }
                    }
                  )
                )
              }
            )
            # TODO: SYMMETRIC?
            # COMPUTE WSD
            res[i, j] <- res[j, i] <- WSD(
              rep(0, nrow(x_mat)),
              knn,
              p = 2
            )
            # UPDATE PROGRESS BAR
            cyto_progress(pb)
          }
        }
        # ENTROPY DISTANCES TO DISTANCE MATRIX
        res <- as.matrix(
          as.dist(
            res
          )
        )
      }
    # COMPARE PARENT CHANNEL DISTRIBUTIONS - LIST OF DISTANCE MATRICES
    } else {
      # DEFAULT STATISTIC -> WASSERSTEIN
      if(!grepl("^(j|k|w|o)", type, ignore.case = TRUE)) {
        type <- "ws"
      }
      # COMPUTE CHANNEL RANGES ACROSS GROUPS - SET BANDWIDTH FOR JSD
      if(grepl("^j", type, ignore.case = TRUE)) {
        limits <- do.call(
          "rbind",
          lapply(
            x_list,
            function(w) {
              cyto_apply(
                w,
                input = "matrix",
                FUN = "cyto_stat_range",
                simplify = TRUE
              )
            }
          )
        )
        limits <- apply(
          limits,
          2,
          function(w){
            range(w, na.rm = TRUE)
          }
        )
      }
      # OVERLAPPLING PACKAGE REQUIRED
      if(grepl("^o", type, ignore.case = TRUE)) {
        cyto_require(
          "overlapping",
          source = "CRAN",
          version = "2.0.0"
        )
      }
      # PROGRESS BAR
      pb <- cyto_progress(
        label = "cyto_dist()",
        total = sum(seq(length(x_list), 1)) * length(channels)
      )
      # CHANNEL-WISE DISTANCE MATRIX
      res <- structure(
        lapply(
          channels,
          function(z) {
            # DISTANCE MATRIX
            d <- matrix(
              0,
              ncol = length(x_list),
              nrow = length(x_list),
              dimnames = list(
                names(x_list),
                names(x_list)
              )
            )
            # COMPUTE DISTANCES
            for(i in 1:length(x_list)){
              for(j in i :length(x_list)) {
                # A - VECTOR
                a <- cyto_exprs(
                  x_list[[i]][[1]],
                  channels = z,
                  drop = TRUE
                )
                # B - VECTOR
                b <- cyto_exprs(
                  x_list[[j]][[1]],
                  channels = z,
                  drop = TRUE
                )
                # TODO: 512 DENSITY ESTIMATE REQUIRED? USE LESS?
                # JENSEN-SHANNON -> NON-ZERO REQUIRED
                if(grepl("^j", type, ignore.case = TRUE)) {
                  d[i, j] <- d[j, i] <- JSD(
                    cyto_stat_density(
                      a,
                      bins = 256,
                      smooth = smooth,
                      limits = limits[, z, drop = TRUE]
                    )$y,
                    cyto_stat_density(
                      b,
                      bins = 256,
                      smooth = smooth,
                      limits = limits[, z, drop = TRUE]
                    )$y
                  )
                # KOLMOGOROV-SMIRNOV
                } else if(grepl("^k", type, ignore.case = TRUE)) {
                  d[i, j] <- d[j, i] <- suppressWarnings(
                    ks.test(
                      a,
                      b
                    )$statistic
                  )
                # WASSERSTEIN 
                } else if(grepl("^w", type, ignore.case = TRUE)) {
                  d[i, j] <- d[j, i] <- WSD(
                    a,
                    b
                  )
                # OVERLAP - NON-OVERLAPPING PORTION (DISTANCE)
                } else {
                  d[i, j] <- d[j, i] <- 1 - cyto_func_call(
                    "overlap",
                    list(
                      list(a, b),
                      type = "1"
                    )
                  )$OV
                }
                # UPDATE PROGRESS BAR
                cyto_progress(pb)
              }
            }
            # CONVERT TO DISTANCE MATRIX
            d <- as.matrix(
              as.dist(d)
            )
            # HEATMAP
            
            return(d)
          }
        ),
        names = channels
      )
      # MARKER NAMES
      if(markers) {
        markers <- cyto_markers(x)
        nms <- names(res)
        ind <- which[nms %in% names(markers)]
        if(length(ind)>0) {
          nms[ind] <- markers[nms[ind]]
          names(res) <- nms
        }
      }
    }
  }
  
  # HEATMAP --------------------------------------------------------------------
  
  # HEATMAPR - TRYCATCH FOR LEGEND
  if(heatmap) {
    if(is.list(res)) {
      lapply(
        seq_along(res),
        function(z) {
          tryCatch(
            heat_map(
              res[[z]],
              title = names(res)[z],
              ...
            ),
            error = function(e) {
              return(NULL)
            }
          )
        }
      )
    } else {
      tryCatch(
        heat_map(
          res,
          ...
        ),
        error = function(e) {
          return(NULL)
        }
      )
    }
  }
  
  # DISTANCE MATRICES ----------------------------------------------------------
  
  # RETURN DISTANCE MATRICES
  return(res)
  
}

#' Wasserstein distance
#' @importFrom graphics hist
#' @noRd
WSD <- function(a,
                b,
                p = 1,
                wa = NULL,
                wb = NULL) {
  
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  if (m == n && is.null(wa) && is.null(wb)) {
    return(mean(abs(sort(b)-sort(a))^p)^(1/p))
  }
  stopifnot(is.null(wa) || length(wa) == m)
  stopifnot(is.null(wb) || length(wb) == n)
  if (is.null(wa)) {
    wa <- rep(1,m)
  } else { # remove points with zero weight
    wha <- wa > 0
    wa <- wa[wha]
    a <- a[wha]
    m <- length(a)
  }
  if (is.null(wb)) {
    wb <- rep(1,n)
  } else { # remove points with zero weight
    whb <- wb > 0
    wb <- wb[whb]
    b <- b[whb]
    n <- length(b)
  }
  
  orda <- order(a)
  ordb <- order(b)
  a <- a[orda]
  b <- b[ordb]
  wa <- wa[orda]
  wb <- wb[ordb]
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))  
  cub <- c(cumsum(ub))  
  arep <- hist(cub, breaks = c(-Inf, cua, Inf), plot = FALSE)$counts + 1
  brep <- hist(cua, breaks = c(-Inf, cub, Inf), plot = FALSE)$counts + 1
  # we sum over rectangles with cuts on the vertical axis each time one of the 
  # two ecdfs makes a jump 
  # arep and brep tell us how many times each of the a and b data have
  # to be repeated in order to get the points on the horizontal axis
  # note that sum(arep)+sum(brep) = m+n-1 
  # (we do not count the height-zero final rectangle where both ecdfs jump to 1)
  
  aa <- rep(a, times=arep)
  bb <- rep(b, times=brep)
  
  uu <- sort(c(cua,cub))
  uu0 <- c(0,uu)
  uu1 <- c(uu,1)
  areap <- sum((uu1-uu0)*abs(bb-aa)^p)^(1/p)
  #  print(rbind(uu1-uu0, pmax(aa,bb)-pmin(aa,bb)))
  return(areap)
  
}

#' Kullback-Leibler divergence
#' @noRd
KLD <- function(A, B) {
  sum(A * log(A/B))
}

#' Jensen-Shannon divergence
#' @noRd
JSD <- function(P, Q) {
  M <- (P + Q)/2
  jsd <- 0.5 * KLD(P, M) + 0.5 * KLD(Q, M)
  return(jsd)
}
