## UMAP ------------------------------------------------------------------------

#' UMAP
#' @noRd
.cyto_map_umap <- function(x,
                           ...) {
  
  # CHECK PYTHON MODULE
  py_umap <- cyto_require(
    "umap-learn",
    python = TRUE,
    pip = TRUE,
    ref = paste0(
      "McInnes L., Healy J. & Melville J. (2018) UMAP: Uniform Manifold ",
      "Approximation and Projection for Dimension Reduction. ",
      "arXiv:1802.03426v3"
    ),
    import = "umap.umap_"
  )
  
  # PYTHON MODULE LOADED
  if(!is.null(py_umap)) {
    # MESSAGE
    message(
      "Using umap-learn python module to compute UMAP co-ordinates..."
    )
    # CONFIGURE UMAP - MODIFY PARAMETERS AS NECESSARY
    umap <- py_umap$UMAP(...)
    # APPLY UMAP TO DATASET
    res <- umap$fit_transform(data.matrix(x))
  # RESORT TO UWOT
  } else {
    # UWOT
    cyto_require(
      "uwot",
      source = "CRAN",
      ref = paste0(
        "McInnes L., Healy J. & Melville J. (2018) UMAP: Uniform Manifold ",
        "Approximation and Projection for Dimension Reduction. ",
        "arXiv:1802.03426v3"
      )
    )
    # MESSAGE
    message(
      "Using uwot::umap() to compute UMAP co-ordinates..."
    )
    # UMAP CO-ORDINATES
    res <- cyto_func_call(
      "uwot::umap",
      list(x, ...)
    )
  }
  
  # PREPARE CO-ORDINATES
  colnames(res) <- paste0(
    "UMAP-",
    1:ncol(res)
  )
  
  # UMAP CO-ORDINATES
  return(res)
  
}

## PACMAP ----------------------------------------------------------------------

#' PaCMAP
#' @noRd
.cyto_map_pacmap <- function(x,
                             ...) {
  
  # CHECK PYTHON MODULE
  py_pacmap <- cyto_require(
    "pacmap",
    python = TRUE,
    pip = TRUE,
    ref = paste0(
      "Wang Y. et al. (2021) Understanding How Dimension Reduction Tools ",
      "Work: An Empirical Approach to Deciphering t-SNE, UMAP, TriMap, and ",
      "PaCMAP for Data Visualization. JMLR 22(201)"
    ),
    import = "pacmap"
  )
  
  # PACMAP UNAVAILABLE
  if(is.null(py_pacmap)) {
    stop(
      "Unable to import required pacmap module from current python environment!"
    )
  # PACMAP
  } else {
    # MESSAGE
    message(
      "Using pacmap python module to compute PaCMAP co-ordinates..."
    )
    # CONFIGURE PACMAP - MODIFY PARAMETERS AS NECESSARY
    pacmap <- py_pacmap$PaCMAP(
      ...
    )
    # APPLY PACMAP TO DATASET
    res <- pacmap$fit_transform(
      data.matrix(x), 
      init = "pca"
    )
    colnames(res) <- paste0(
      "PaCMAP-",
      1:ncol(res)
    )
    return(res)
  }
  
  
}

## TRIMAP ----------------------------------------------------------------------

#' TriMap
#' @noRd
.cyto_map_trimap <- function(x,
                             ...) {
  
  # CHECK PYTHON MODULE
  py_trimap <- cyto_require(
    "trimap",
    python = TRUE,
    ref = paste0(
      "Amid E. & Warmuth M. (2022) TriMap: Large-scale Dimensionality ","
      Reduction Using Triplets. arXiv:1910.00204"
    ),
    import = "trimap",
    pip = TRUE
  )
  
  # TRIMAP UNAVAILABLE
  if(is.null(py_trimap)) {
    stop(
      "Unable to import required trimap module from current python environment!"
    )
  # TRIMAP
  } else {
    # MESSAGE
    message(
      "Using trimap python module to compute TriMap co-ordinates..."
    )
    # CONFIGURE TRIMAP - MODIFY PARAMETERS AS NECESSARY
    trimap <- py_trimap$TRIMAP(
      ...
    )
    # APPLY TRIMAP TO DATASET
    res <- trimap$fit_transform(
      data.matrix(x)
    )
    colnames(res) <- paste0(
      "TriMap-",
      1:ncol(res)
    )
    return(res)
  }
  
}

## ISOMAP ----------------------------------------------------------------------

#' IsoMap
#' @noRd
.cyto_map_isomap <- function(x,
                             n_jobs = -2L,
                             ...) {
  
  # CHECK PYTHON MODULE
  py_isomap <- cyto_require(
    "sklearn",
    python = TRUE,
    ref = NULL,
    import = "sklearn.manifold",
    pip = TRUE
  )
  
  # ISOMAP UNAVAILABLE
  if(is.null(py_isomap)) {
    stop(
      "Unable to import required sklearn module from current python environment!"
    )
  # ISOMAP
  } else {
    # MESSAGE
    message(
      "Using isomap sklearn module to compute isomap co-ordinates..."
    )
    # CONFIGURE ISOMAP - MODIFY PARAMETERS AS NECESSARY
    isomap <- py_isomap$Isomap(
      n_jobs = n_jobs,
      ...
    )
    # APPLY ISOMAP TO DATASET
    res <- isomap$fit_transform(
      data.matrix(x)
    )
    colnames(res) <- paste0(
      "IsoMap-",
      1:ncol(res)
    )
    return(res)
  }
  
}

## PHATE -----------------------------------------------------------------------

#' PHATE
#' @noRd
.cyto_map_phate <- function(x,
                            n_jobs = -2L,
                            ...) {
  
  # CHECK PYTHON MODULE
  py_phate <- cyto_require(
    "phate",
    python = TRUE,
    ref = paste0(
      "Moon, van Dijk, Wang, Gigante et al. Visualizing Transitions and ",
      "Structure for Biological Data Exploration. 2019. Nature Biotechnology."
    ),
    import = "phate",
    pip = TRUE
  )
  
  # PHATE UNAVAILABLE
  if(is.null(py_phate)) {
    stop(
      "Unable to import required phate module from current python environment!"
    )
  # PHATE
  } else {
    # MESSAGE
    message(
      "Using phate python module to compute PHATE co-ordinates..."
    )
    # CONFIGURE TRIMAP - MODIFY PARAMETERS AS NECESSARY
    phate <- py_phate$PHATE(
      n_jobs = n_jobs,
      ...
    )
    # APPLY TRIMAP TO DATASET
    res <- phate$fit_transform(
      data.matrix(x)
    )
    colnames(res) <- paste0(
      "PHATE-",
      1:ncol(res)
    )
    return(res)
  }
  
}

## TSNE -----------------------------------------------------------------------

#' t-SNE
#' @noRd
.cyto_map_tsne <- function(x,
                           n_jobs = -2L,
                           num_threads = 0,
                           ...) {
  
  # CHECK PYTHON MODULE
  py_tsne <- cyto_require(
    "openTSNE",
    python = TRUE,
    pip = TRUE,
    ref = paste0(
      "Maaten, Laurens van der, and Geoffrey Hinton. Visualizing data ",
      " using t-SNE. Journal of machine learning research 9.Nov ",
      "(2008): 2579-2605", "\n",
      "Linderman, George C., et al. Efficient Algorithms for t-distributed ",
      "Stochastic Neighborhood Embedding. arXiv:1712.09005 (2017)"
    ),
    import = "openTSNE"
  )

  # PYTHON MODULE LOADED
  if(!is.null(py_tsne)) {
    # MESSAGE
    message(
      "Using openTSNE python module to compute t-SNE co-ordinates..."
    )
    # CONFIGURE UMAP - MODIFY PARAMETERS AS NECESSARY
    tsne <- py_tsne$TSNE(
      n_jobs = n_jobs,
      ...
    )
    res <- tsne$fit(
      data.matrix(x)
    )
  # RESORT TO RTSNE
  } else {
    # RTSNE
    cyto_require(
      "Rtsne",
      source = "CRAN",
      ref = paste0(
        "Maaten, Laurens van der, and Geoffrey Hinton. Visualizing data ",
        " using t-SNE. Journal of machine learning research 9.Nov ",
        "(2008): 2579-2605"
      )
    )
    # MESSAGE
    message(
      "Using Rtsne::Rtsne() to compute t-SNE co-ordinates..."
    )
    # TSNE CO-ORDINATES
    res <- cyto_func_call(
      "Rtsne::Rtsne",
      list(
        x, 
        num_threads = 0,
        check_duplicates = FALSE,
        ...
      )
    )
  }
  
  # PREPARE CO-ORDINATES
  colnames(res) <- paste0(
    "t-SNE-",
    1:ncol(res)
  )
  
  # TSNE CO-ORDINATES
  return(res)
  
}

## EMBEDSOM --------------------------------------------------------------------

#' EmbedSOM
#' @noRd
.cyto_map_embedsom <- function(x,
                               parallel = TRUE,
                               ...) {
  
  # EmbedSOM
  cyto_require(
    "EmbedSOM",
    source = "CRAN",
    ref = paste0(
      "Kratchovil M., Koladiya A. & Vondrasek J. (2019) Generalised ",
      "EmbedSOM on quadtree-structured self-organising maps. F1000 ",
      "Research (8:2120)"
    )
  )
  # MESSAGE
  message(
    "Using EmbedSOM::EmbedSOM() to compute EmbedSOM co-ordinates..."
  )
  # ARUMENTS
  args <- list(
    "data" = x,
    "parallel" = parallel,
    ...
  )
  # DATA
  names(args)[1] <- "data"
  # CREATE SOM - FLOWSOM NOT SUPPLIED
  if(!any(c("fsom", "map") %in% names(args))) {
    # DEFAULT SOM GRID SIZE - X
    if(!"xdim" %in% names(args)) {
      args[["xdim"]] <- 24
    }
    # DEFAULT SOM GRID SIZE - Y
    if(!"ydim" %in% names(args)) {
      args[["ydim"]] <- 24
    }
    # SOM
    args[["map"]] <- cyto_func_execute(
      "EmbedSOM::SOM",
      args
    )
  }
  # EMBEDSOM
  res <- cyto_func_execute(
    "EmbedSOM::EmbedSOM",
    args
  )
  colnames(res) <- paste0(
    "EmbedSOM-",
    1:ncol(res)
  )
  
  # TSNE CO-ORDINATES
  return(res)
  
}

## PCA -------------------------------------------------------------------------

#' PCA
#' @noRd
.cyto_map_pca <- function(x,
                           ...) {
  
  # RTSNE
  cyto_require(
    "rsvd",
    source = "CRAN",
    ref = paste0(
      "Erichson NB, Voronin S, Brunton SL, Kutz JN (2019). Randomized ",
      "Matrix Decompositions Using R. Journal of Statistical Software, ",
      "89(11), 1–48. doi: 10.18637/jss.v089.i11."
    )
  )
  # MESSAGE
  message(
    "Using rsvd::rpca() to compute PCA co-ordinates..."
  )
  # TSNE CO-ORDINATES
  res <- cyto_func_call(
    "rsvd::rpca",
    list(x, ...)
  )$x
  
  # PREPARE CO-ORDINATES
  colnames(res) <- paste0(
    "PCA-",
    1:ncol(res)
  )
  
  # TSNE CO-ORDINATES
  return(res)
  
}

## FIt-SNE ---------------------------------------------------------------------

#' FIt-SNE
#' @noRd
.cyto_map_fitsne <- function(x,
                             ...) {
  
  # RTSNE
  message(
    "Linderman, George C., et al. Efficient Algorithms for t-distributed ",
    "Stochastic Neighborhood Embedding. arXiv:1712.09005 (2017)"
  )
  # MESSAGE
  message(
    "Using FIt-SNE::fftRtsne() to compute FIt-SNE co-ordinates..."
  )
  # FITSNE CO-ORDINATES
  res <- cyto_func_call(
    "CytoExploreR::fftRtsne",
    list(x, ...)
  )
  
  # PREPARE CO-ORDINATES
  colnames(res) <- paste0(
    "FIt-SNE-",
    1:ncol(res)
  )
  
  # FITSNE CO-ORDINATES
  return(res)
  
}
