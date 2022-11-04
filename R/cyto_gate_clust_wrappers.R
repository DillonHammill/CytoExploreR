## CYTO_CLUST_SOM --------------------------------------------------------------

#' Self organising maps
#' @noRd
.cyto_clust_SOM <- function(x,
                            xdim = 14,
                            ydim = 14,
                            ...) {
  
  # LOAD KOHONEN
  cyto_require(
    "kohonen",
    source = "CRAN",
    ref = paste0(
      "Wehrens, R., & Buydens, L. M. C. (2007). Self- and Super-organizing ",
      "Maps in R: The kohonen Package. Journal of Statistical Software, ",
      "21(5), 1–19."
    )
  )
  
  # PREPARE ARGUMENTS
  args <- .args_list(...)
  names(args)[match("x", names(args))] <- "data"
  
  # DATA REQUIRED
  args[["keep.data"]] <- TRUE
  
  # PREPARE SOM GRID
  args$grid <- cyto_func_call(
    "kohonen::somgrid",
    list(
      xdim = xdim,
      ydim = ydim
    )
  )
  
  # TRAIN SOM - RETURN LABELS
  return(
    cyto_func_execute(
      "kohonen::som",
      args = args
    )$unit.classif
  )
  
}

## CYTO_CLUST_FLOWSOM ----------------------------------------------------------

#' FlowSOM 
#' @noRd
.cyto_clust_FlowSOM <- function(x,
                                compensate = FALSE,
                                transform = FALSE,
                                silent = FALSE,
                                seed = 2022,
                                nClus = NULL,
                                k = NULL,
                                xdim = 14,
                                ydim = 14,
                                ...) {
  
  # LOAD FLOWSOM
  cyto_require(
    "FlowSOM",
    source = "BioC",
    repo = "SofieVG/FlowSOM",
    ref = paste0(
      "Van Gassen S, et al. (2015). FlowSOM: Using self-organising maps ",
      "for visualisation and interpretation of cytometry data. Cytometry ",
      "A 87(7) 10.1002/cyto.a.22625"
    )
  )
  
  # ARGUMENTS
  args <- .args_list(...)
  
  # K -> NCLUS
  ind <- grep("^k$", names(args), ignore.case = TRUE)
  if(length(ind) > 0) {
    if(!is.null(args[[ind]])) {
      args$nClus <- args[[ind]]
      args <- args[-ind]
    }
  }
  
  # MESSAGE
  message(
    "Using FlowSOM:FlowSOM() to compute cluster labels..."
  )
  
  # ARGUMENTS
  args <- c(
    list(
      input = x,
      colsToUse = cyto_channels(x)
    ),
    args
  )
  
  # READINPUT
  fsom <- cyto_func_execute(
    "FlowSOM::ReadInput",
    args
  )
  
  # UPDATE ARGUMENTS
  args$fsom <- fsom
  
  # BUILD SOM
  args$fsom <- cyto_func_execute(
    "FlowSOM::BuildSOM",
    args = args
  )
  
  # LABELS
  labels <- args$fsom$map$mapping[, 1]
  
  # DATA
  args$data <- args$fsom$map$codes
  
  # AUTOMATED NCLUS SELECTION
  if(is.null(args$nClus)) {
    args$labels <- cyto_func_execute(
      "FlowSOM::MetaClustering",
      args = args
    )
  # MANUAL NCLUS
  } else {
    args$k <- args$nClus
    args$labels <- cyto_func_execute(
      "FlowSOM::metaClustering_consensus",
      args = args
    )
  }
  
  # CLUSTER ASSIGNMENTS
  return(args$labels[labels])
  
}

## CYTO_CLUST_PYPHENOGRAPH -------------------------------------------------------

#' Phenograph
#' @noRd
.cyto_clust_pyphenograph <- function(x,
                                     k = 30L,
                                     ...) {
  
  # CHECK PYTHON MODULE
  py_phenograph <- cyto_require(
    "PhenoGraph",
    python = TRUE,
    pip = TRUE,
    ref = paste0(
      "Levine, J. et al. (2018) Data-Driven Phenotypic Dissection of AML ",
      "Reveals Progenitor-like Cells that Correlate with Prognosis. ",
      "Cell 162(1)"
    ),
    import = "phenograph"
  )
  
  # PYTHON MODULE REQUIRED
  if(is.null(py_phenograph)) {
    stop(
      paste0(
        "Unable to import required PhenoGraph module from current python ",
        "environment!"
      )
    )
  # PYTHON MODULE LOADED
  } else {
    # MESSAGE
    message(
      "Using PhenoGraph python module to compute cluster labels..."
    )
    # PHENOGRAPH
    py_phenograph$cluster(
      x,
      k = as.integer(k),
      ...
    )$communities
  }
  
}

## CYTO_CLUST_RPHENOGRAPH ------------------------------------------------------

#' RPhenograph
#' @noRd
.cyto_clust_rphenograph <- function(x,
                                    k = 30L,
                                    ...) {
  
  # LOAD RPHENOGRAPH
  cyto_require(
    "Rphenograph",
    source = "GitHUb",
    repo = "JinmiaoChenLab/Rphenograph",
    ref = paste0(
      "Levine, J. et al. (2018) Data-Driven Phenotypic Dissection of AML ",
      "Reveals Progenitor-like Cells that Correlate with Prognosis. ",
      "Cell 162(1)"
    )
  )
  
  # RPHENOGRAPH
  res <- cyto_func_call(
    "Rphenograph::Rphenograph",
    args = list(
      x,
      k = as.integer(k),
      ...
    )
  )
  
  # PREPARE GATE - ORDER & OUTLIERS -> IN CASE INCORRECT RPHENOGRAPH
  gate <- rep(0, nrow(x))
  gate[
    as.numeric(
      cyto_func_call(
        "igraph::V",
        args = list(
          res[[1]]
        )
      )$name
    )
  ] <- res[[2]]$membership
  
  return(gate)
  
}

## CYTO_CLUST_DEPECHE ----------------------------------------------------------

#' Depeche
#' @noRd
.cyto_clust_depeche <- function(x,
                                createOutput = FALSE,
                                ...) {
  
  # LOAD DEPECHER
  cyto_require(
    "DepecheR",
    source = "BioC",
    ref = paste0(
      "Theorell A, Bryceson Y, Theorell J (2019). Determination of ",
      "essential phenotypic elements of clusters in high-dimensional ",
      "entities - DEPECHE. PLoS One."
    )
  )
  
  # DEPECHER
  cyto_func_call(
    "DepecheR::depeche",
    args = list(
      x,
      createOutput = createOutput,
      ...
    )
  )$clusterVector
  
}

## CYTO_CLUST_IMMUNOCLUST ------------------------------------------------------

#' ImmunoClust
#' @noRd
.cyto_clust_immunoclust <- function(x,
                                    ...) {
  
  # LOAD IMMUNOCLUST
  cyto_require(
    "immunoClust",
    source = "BioC",
    ref = paste0(
      "Sorensen T, et al. (2015). immunoClust--An automated analysis ",
      "pipeline for the identification of immunophenotypic signatures ",
      "in high-dimensional cytometric datasets. Cytometry A 87(7) ",
      "10.1002/cyto.a.22626"
    )
  )
  
  # IMMUNOCLUST
  cyto_func_call(
    "immunoClust::cell.process",
    args = list(
      x,
      parameters = cyto_channels(x),
      ...
    )
  )@label
  
}

## CYTO_CLUST_DBSCAN -----------------------------------------------------------

#' DBSCAN
#' @noRd
.cyto_clust_dbscan <- function(x,
                               type = "dbscan",
                               ...) {
  
  # LOAD DBSCAN
  cyto_require(
    "dbscan",
    source = "CRAN",
    ref = paste0(
      "Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast Density-Based ",
      "Clustering with R. Journal of Statistical Software, 91(1), 1–30."
    )
  )
  
  # HDSCAN
  if(grepl("^hdbscan$", type, ignore.case = TRUE)) {
    type <- "dbscan::hdbscan"
  # DBSCAN
  } else if(grepl("^dbscan$", type, ignore.case = TRUE)) {
    type <- "dbscan::dbscan"
  # SNNclust
  } else if(grepl("^SNNclust", type, ignore.case = TRUE)) {
    type <- "dbscan::sNNclust"
  # JARVIS-PATRICK CLUSTERING
  } else if(grepl("^jpclust", type, ignore.case = TRUE)) {
    type <- "dbscan::jpclust"
  }
  
  # DBSCAN
  cyto_func_call(
    type,
    args = list(
      data.matrix(x),
      ...
    )
  )$cluster
  
}

## CYTO_CLUST_HGC --------------------------------------------------------------

#' HGC
#' @noRd
.cyto_clust_hgc <- function(x,
                            tree = "SNN",
                            n = NULL,
                            ...) {
  
  # LOAD HGC
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
  
  # ARGUMENTS
  args <- list(
    mat = x,
    ...
  )
  
  # SNN
  if(.grepl("^SNN$", tree, ignore.case = TRUE)) {
    tree <- cyto_func_execute(
      "HGC::SNN.Construction",
      args = args
    )
  # KNN
  } else if(.grepl("^kNN$", tree, ignore.case = TRUE)) {
    tree <- cyto_func_execute(
      "HGC::KNN.Construction",
      args = args
    )
  # RNN
  } else if(.grepl("^RNN$", tree, ignore.case = TRUE)) {
    tree <- cyto_func_execute(
      "HGC::RNN.Construction",
      args = args
    )
  # CKNN
  } else if(.grepl("^CKNN$", tree, ignore.case = TRUE)) { 
    tree <- cyto_func_execute(
      "HGC::CKNN.Construction",
      args = args
    )
  # MST
  } else if(.grepl("^MST$", tree, ignore.case = TRUE)) { 
    tree <- cyto_func_execute(
      "HGC::MST.Construction",
      args = args
    )
  # PMST
  } else if(.grepl("^PMST$", tree, ignore.case = TRUE)) { 
    tree <- cyto_func_execute(
      "HGC::PMST.Construction",
      args = args
    )
  }
  
  # COMMUNITY DETECTION
  if(is.null(n)) {
    # IGRAPH REQUIRED
    cyto_require(
      "igraph",
      source = "CRAN"
    )
    # CREATE UNDIRECTED WEIGHTED GRAPH
    g <- cyto_func_call(
      "igraph::graph.adjacency",
      args = list(
        tree,
        mode = "undirected",
        weighted = TRUE
      )
    )
    # COMMUNITY DETECTION
    cyto_func_call(
      "igraph::cluster_louvain",
      args = list(
        graph = g
      )
    )$membership
  # HIERARCHICAL CLUSTERING
  } else {
    # ARGUMENTS CREATE DENDROGRAM
    args$G <- tree
    # DENDROGRAM
    tree <- cyto_func_execute(
      "HGC::HGC.dendrogram",
      args
    )
    # K - CLUSTERS
    args$k <- args$n
    args$tree <- tree
    # CUT DENDROGRAM
    cyto_func_execute(
      "stats::cutree",
      args = args
    )
  }
}

## CYTO_CLUST_RCLUSTERPP -------------------------------------------------------

#' Rclusterpp
#' @noRd
.cyto_clust_rclusterpp <- function(x,
                                   ...) {
  
  # LOAD RCLUSTERPP
  cyto_require(
    "Rclusterpp",
    source = "CRAN"
  )
  
  # ARGUMENTS
  args <- list(x, ...)
  
  # CHECK ARGUMENTS
  if(!"k" %in% names(args)) {
    stop(
      paste0(
        "Rclusterpp requires argument 'k' to indicate the desired number ",
        "of clusters!"
      )
    )
  }
  
  # RCLUSTERPP
  tree <- cyto_func_execute(
    "Rclusterpp::Rclusterpp.hclust",
    args = args
  )
  
  # CUT TREE
  cyto_func_execute(
    "stats::cutree",
    args = c(
      args, 
      list(
        "tree" = tree
      )
    )
  )
  
}


## CYTO_CLUST_FLOWPEAKS --------------------------------------------------------

#' FlowPeaks
#' @noRd
.cyto_clust_flowpeaks <- function(x,
                                  ...) {
  
  # LOAD FLOWPEAKS
  cyto_require(
    "flowPeaks",
    source = "BioC",
    ref = paste0(
      "Ge Y. et al (2012) flowPeaks: a fast unsupervised clustering for ",
      "flow cytometry data via K-means and density peak fnding. ",
      "Bioinformatics 8(15):2052-8."
    )
  )
  
  # FLOWPEAKS
  cyto_func_call(
    "flowPeaks::flowPeaks",
    args = list(
      x,
      ...
    )
  )$peaks.cluster
  
}

## CYTO_CLUST_MSTKNN -----------------------------------------------------------

#' MST-KNN
#' @noRd
.cyto_clust_mstknn <- function(x,
                               ...) {
  
  # LOAD FLOWPEAKS
  cyto_require(
    "mstknnclust",
    source = "CRAN",
    ref = paste0(
      "Inostroza-Ponta, M. (2008) An integrated and scalable approach based ",
      "on combinatorial optimization techniques for the analysis of ",
      "microarray data, University of Newcastle, ISBN."
    )
  )
  
  # COMPUTE DISTANCE MATRIX
  d <- as.matrix(
    cyto_func_execute(
      "stats::dist",
      args = list(
        "x" = x,
        ...
      )
    )
  )
  
  # MSTKNN CLUSTERING
  cyto_func_execute(
    "mstknnclust::mst.knn",
    args = list(
      distance.matrix = d,
      ...
    )
  )$cluster
  
}

## CYTO_CLUST_SPECTRUM ---------------------------------------------------------

#' Spectrum spectral clustering
#' @noRd
.cyto_clust_spectrum <- function(x,
                                 ...) {
  
  # LOAD SPECTRUM
  cyto_require(
    "Spectrum",
    source = "CRAN"
  )
  
  # SPECTRAL CLUSTERING
  cyto_func_call(
    "Spectrum::Spectrum",
    args = list(
      t(x[1:nrow(x), ]),
      ...
    )
  )$assignments
  
}

## CYTO_CLUST_KMEANS -----------------------------------------------------------

#' Kmeans clustering
#' @noRd
.cyto_clust_kmeans <- function(x,
                               k = NULL,
                               centers = NULL,
                               ...) {
  

  # KMEANS CLUSTERING
  cyto_func_call(
    "stats::kmeans",
    args = list(
      x,
      centers = c(k, centers)[1],
      ...
    )
  )$cluster
  
}
