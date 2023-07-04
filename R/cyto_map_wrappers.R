## HAISU -----------------------------------------------------------------------

#' HAISU
#' @noRd
.cyto_map_haisu <- function(x,
                            nodes = NULL,
                            factor = 0.15,
                            metric = "euclidean",
                            n_jobs = -1L,
                            ...) {
  
  # NODES REQUIRED
  if(is.null(nodes)) {
    stop(
      "'nodes' must be supplied to use Haisu for dimenionality reduction!"
    )
  }
  
  # IGRAPH REQUIRED - MST GRAPH & ADJACENCY MATRIX
  cyto_require(
    "igraph",
    source = "CRAN"
  )
  
  # NODES
  nodes <- as.character(nodes)
  
  # PARTITION DATA
  x_split <- structure(
    lapply(
      unique(nodes),
      function(z) {
        x[
          match(z, nodes), , drop = FALSE
        ]
      }
    ),
    names = unique(nodes)
  )
  
  # COMPUTE CENTROIDS
  x_split <- do.call(
    "rbind",
    lapply(
      x_split,
      function(z) {
        apply(
          z,
          2,
          "median"
        )
      }
    )
  )
  
  # COMPUTE DISTANCE MATRIX
  adj <- cyto_func_call(
    "stats::dist",
    list(
      x_split
    )
  )
  
  # BUILD ADJACENCY GRAPH
  adj <- cyto_func_call(
    "igraph::graph.adjacency",
    list(
      as.matrix(adj),
      mode = "undirected",
      weighted = TRUE
    )
  )
  
  # MST TREE
  adj <- cyto_func_call(
    "igraph::minimum.spanning.tree",
    list(
      adj,
      algorithm = "prim"
    )
  )
  
  # EXTRACT ADJACENCY MATRIX
  adj <- cyto_func_call(
    "igraph::as_adjacency_matrix",
    list(
      adj
    )
  )
  
  # INSTALL HAISU REQUIREMENTS -------------------------------------------------
  cyto_require(
    "numpy",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "sys",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "networkx",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "sklearn.metrics",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "scipy",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "shapely",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "scipy.spatial",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "multiprocessing",
    python = TRUE,
    pip = TRUE
  )
  
  cyto_require(
    "ctypes",
    python = TRUE,
    pip = TRUE
  )
  
  HAISU <- function(){}
  
  # IMPORT HAISU PYTHON SCRIPT
  reticulate::source_python(
    "https://raw.githubusercontent.com/Cobanoglu-Lab/Haisu/master/Haisu.py"
  )

  # INITIALISE HAISU
  h <- HAISU(
    unique(nodes),
    adj
  )
  
  # RUN HAISU
  h$get_pairwise_matrix(
    x,
    as.character(nodes),
    factor = factor,
    transpose = TRUE,
    normalize = FALSE,
    metric = metric,
    n_jobs = n_jobs
  )

}

## UMAP ------------------------------------------------------------------------

#' UMAP
#' @noRd
.cyto_map_umap <- function(x,
                           nodes = NULL,
                           metric = "euclidean",
                           factor = 0.15,
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
  
  # HAISU
  if(!is.null(nodes) & !is.null(py_umap)) {
    # COMPUTE PAIRWISE MATRIX
    x <- cyto_func_execute(
      ".cyto_map_haisu",
      list(
        x = x,
        nodes = nodes,
        metric = metric,
        factor = factor,
        ...
      )
    )
    # SET METRIC TO PRECOMPUTED
    metric <- "precomputed"
  }
  
  # PYTHON MODULE LOADED
  if(!is.null(py_umap)) {
    # MESSAGE
    message(
      "Using umap-learn python module to compute UMAP co-ordinates..."
    )
    # CONFIGURE UMAP - MODIFY PARAMETERS AS NECESSARY
    umap <- py_umap$UMAP(
      metric = metric,
      ...
    )
    # APPLY UMAP TO DATASET
    res <- umap$fit_transform(
      data.matrix(x)
    )
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
    # HAISU NOT SUPPORTED - METRIC != "PRECOMPUTED"
    # UMAP CO-ORDINATES
    res <- cyto_func_call(
      "uwot::umap",
      list(
        x,
        metric = metric,
        ...
      )
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

## UMATO -----------------------------------------------------------------------

#' UMATO
#' @noRd
.cyto_map_umato <- function(x,
                            hub_num = 20L,
                            ...) {
  
  # CHECK PYTHON MODULE
  py_umato <- cyto_require(
    "umato",
    python = TRUE,
    ref = NULL,
    import = "umato",
    pip = TRUE
  )
  
  # UMATO UNAVAILABLE
  if(is.null(py_umato)) {
    stop(
      "Unable to import required pacmap module from current python environment!"
    )
    # UMATO
  } else {
    # MESSAGE
    message(
      "Using umato python module to compute UMATO co-ordinates..."
    )
    # CONFIGURE PACMAP - MODIFY PARAMETERS AS NECESSARY
    umato <- py_umato$UMATO(
      hub_num = as.integer(hub_num),
      ...
    )
    # APPLY UMATO TO DATASET
    res <- umato$fit_transform(
      data.matrix(x)
    )
    colnames(res) <- paste0(
      "UMATO-",
      1:ncol(res)
    )
    return(res)
  }
  
}

## MDS -------------------------------------------------------------------------

#' MDS
#' @noRd
.cyto_map_mds <- function(x,
                          n_jobs = -2,
                          ...) {
  
  # CHECK PYTHON MODULE
  py_mds<- cyto_require(
    "sklearn",
    python = TRUE,
    ref = NULL,
    import = "sklearn.manifold",
    pip = TRUE
  )
  
  # PYTHON
  if(!is.null(py_mds)) {
    # MESSAGE
    message(
      "Using MDS sklearn module to compute MDS co-ordinates..."
    )
    # CONFIGURE MDS - MODIFY PARAMETERS AS NECESSARY
    mds <- py_mds$MDS(
      n_jobs = n_jobs,
      ...
    )
    # APPLY MDS TO DATASET
    res <- mds$fit_transform(
      data.matrix(x)
    )
  # R
  } else {
    # MESSAGE
    message(
      "Using stats::cmdscale() to compute MDS co-ordinates..."
    )
    # DISTANCE MATRIX
    x <- cyto_func_execute(
      "stats::dist",
      list(
        x = x,
        ...
      )
    )
    # MDS
    res <- cyto_func_execute(
      "stats::cmdscale",
      list(
        "d" = x,
        "list." = FALSE,
        "eig" = FALSE,
        ...
      )
    )
  }
  
  # PREPARE CO-ORDINATES
  colnames(res) <- paste0(
    "MDS-",
    1:ncol(res)
  )
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
                            nodes = NULL,
                            factor = 0.15,
                            knn_dist = "euclidean",
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
  
  # HAISU
  if(!is.null(nodes) & !is.null(py_phate)) {
    # COMPUTE PAIRWISE MATRIX
    x <- cyto_func_execute(
      ".cyto_map_haisu",
      list(
        x = x,
        nodes = nodes,
        metric = knn_dist,
        factor = factor,
        ...
      )
    )
    # SET METRIC TO PRECOMPUTED
    knn_dist <- "precomputed"
  }
  
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
    # CONFIGURE PHATE - MODIFY PARAMETERS AS NECESSARY
    phate <- py_phate$PHATE(
      knn_dist = knn_dist,
      n_jobs = n_jobs,
      ...
    )
    # APPLY PHATE TO DATASET
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
                           nodes = NULL,
                           factor = 0.15,
                           metric = "euclidean",
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
  
  # HAISU
  if(!is.null(nodes) & !is.null(py_tsne)) {
    # COMPUTE PAIRWISE MATRIX
    x <- cyto_func_execute(
      ".cyto_map_haisu",
      list(
        x = x,
        nodes = nodes,
        metric = metric,
        factor = factor,
        ...
      )
    )
    # SET METRIC TO PRECOMPUTED
    metric <- "precomputed"
  }
  
  # PYTHON MODULE LOADED
  if(!is.null(py_tsne)) {
    # MESSAGE
    message(
      "Using openTSNE python module to compute t-SNE co-ordinates..."
    )
    # CONFIGURE UMAP - MODIFY PARAMETERS AS NECESSARY
    tsne <- py_tsne$TSNE(
      metric = metric,
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
    # HAISU NOT SUPPORTED - CANNOT PASS PRECOMPUTED MATRIX
    # USE IS_DISTANCE?
    # TSNE CO-ORDINATES
    res <- cyto_func_call(
      "Rtsne::Rtsne",
      list(
        x, 
        num_threads = 0,
        check_duplicates = FALSE,
        ...
      )
    )$Y
  }
  
  # PREPARE CO-ORDINATES
  colnames(res) <- paste0(
    "t-SNE-",
    1:ncol(res)
  )
  
  # TSNE CO-ORDINATES
  return(res)
  
}

## NCVIS -----------------------------------------------------------------------

#' NCVis
#' @noRd
.cyto_map_ncvis <- function(x,
                            ...) {
  
  # CHECK PYTHON MODULE
  py_ncvis <- cyto_require(
    "ncvis",
    python = TRUE,
    ref = paste0(
      "Artemenkov, A. & Panov, M. (2020) NCVis: Noise ",
      "Contrastive Approach for Scalable Visualization. In Proceedings ",
      "of The Web Conference 2020 (WWW '20). Association for Computing ",
      "Machinery, New York, NY, USA, 2941–2947."
    ),
    import = "ncvis",
    pip = TRUE
  )
  
  # NCVIS UNAVAILABLE
  if(is.null(py_ncvis)) {
    stop(
      "Unable to import required ncvis module from current python environment!"
    )
  # NCVIS
  } else {
    # MESSAGE
    message(
      "Using ncvis module to compute NCVis co-ordinates..."
    )
    # CONFIGURE  - MODIFY PARAMETERS AS NECESSARY
    ncvis <- py_ncvis$NCVis(
      ...
    )
    # APPLY NCVIS TO DATASET
    res <- ncvis$fit_transform(
      data.matrix(x)
    )
    colnames(res) <- paste0(
      "NCVis-",
      1:ncol(res)
    )
    return(res)
  }
  
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

## KNN -------------------------------------------------------------------------

#' KNN
#' @noRd
.cyto_map_knn <- function(x,
                          k = 25L,
                          ...) {
  
  # COMPUTE KNN
  knn <- .cyto_knn(
    x,
    k = as.integer(k),
    ...
  )
  
  # FORMAT KNN EDGES
  knn <- do.call(
    "rbind",
    lapply(
      seq_len(nrow(knn$nn.idx)),
      function(i) {
        matrix(
          c(rep(knn$nn.idx[i, 1], ncol(knn$nn.idx) - 1),
            knn$nn.idx[i, -1]),
          ncol = 2,
          nrow = ncol(knn$nn.idx) - 1,
          dimnames = list(
            NULL,
            c("node", "matches")
          )
        )
      }
    )
  )
  
  # REQUIRE IGRAPH
  cyto_require(
    "igraph",
    source = "CRAN"
  )
  
  # CONSTRUCT GRAPH
  g <- cyto_func_call(
    "igraph::graph_from_edgelist",
    list(
      el = knn,
      directed = FALSE
    )
  )
  
  # LAYOUT
  set.seed(2022)
  res <- cyto_func_execute(
    "igraph::layout_with_fr",
    list(
      graph = g,
      ...
    )
  )
  colnames(res) <- paste0(
    "kNN-",
    1:ncol(res)
  )
  return(res)
  
}

## SNN -------------------------------------------------------------------------

#' SNN
#' @noRd
.cyto_map_snn <- function(x,
                          k = 25,
                          ...) {
  
  # RANN
  cyto_require(
    "HGC",
    source = "github",
    repo = "XuegongLab/HGC@HGC4oldRVersion",
    ref = paste0(
      "Ziheng Zou, Kui Hua, Xuegong Zhang (2021) HGC: fast hierarchical ",
      "clustering for large-scale single-cell data, Bioinformatics ",
      "37(21)"
    )
  )
  
  # IGRAPH
  cyto_require(
    "igraph",
    source = "CRAN"
  )
  
  # ARGUMENTS
  args <- list(
    mat = x,
    k = k,
    ...
  )
  
  # BUILD SNN
  snn <- cyto_func_execute(
    "HGC::SNN.Construction",
    args
  )
  
  # CONSTRUCT GRAPH
  g <- cyto_func_call(
    "igraph::graph_from_adjacency_matrix",
    list(
      snn,
      mode = "undirected",
      weighted = TRUE
    )
  )
  
  # LAYOUT
  set.seed(2022)
  args$graph <- g
  res <- cyto_func_execute(
    "igraph::layout_with_fr",
    args
  )
  colnames(res) <- paste0(
    "sNN-",
    1:ncol(res)
  )
  return(res)
  
}

## MST -------------------------------------------------------------------------

#' MST
#' @noRd
.cyto_map_mst <- function(x,
                          method = "euclidean",
                          ...) {
  
  # IGRAPH
  cyto_require(
    "igraph",
    source = "CRAN"
  )
  
  # DISTANCE MATRIX
  d <- as.matrix(
    dist(
      x,
      method = method
    )
  )
  
  # GRAPH
  g <- cyto_func_call(
    "igraph::graph.adjacency",
    list(
      d,
      mode = "undirected",
      weighted = TRUE
    )
  )
  
  # MST
  mst <- cyto_func_call(
    "igraph::minimum.spanning.tree",
    list(
      g
    )
  )
  
  # WEIGHTS
  w <- cyto_func_call(
    "igraph::edge.attributes",
    list(
      mst
    )
  )$weight
  w <- w/mean(w)
  
  # UPDATE WEIGHTS
  cyto_func_call(
    "igraph::edge.attributes<-",
    list(
      graph = mst,
      value = list(
        weight = w
      )
    )
  )
  
  # LAYOUT
  res <- cyto_func_call(
    "igraph::layout.kamada.kawai",
    list(
      coords = as.matrix(
        expand.grid(
          seq_len(ceiling(sqrt(nrow(x)))),
          seq_len(ceiling(sqrt(nrow(x))))
        )
      ),
      mst
    )
  )
  colnames(res) <- c("MST-1", "MST-2")
  
  return(res)
  
}
