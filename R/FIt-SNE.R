#' FIt-SNE Based on Kluger Lab FIt-SNE
#'
#' Modified version of \url{https://github.com/KlugerLab/FIt-SNE/}
#' implementation to expose argument names and defaults within CytoExploreR.
#' This function should not be used directly, data should instead be mapped
#' using \code{\link{cyto_map}}.
#'
#' @param X a matrix containing the data to be mapped.
#' @param dims dimensionality of the embedding, set to 2 by default.
#' @param perplexity used to determine the bandwidth of the Gaussian kernel in
#'   the input space, set to 30 by default.
#' @param theta set to 0 for exact t-SNE. If non-zero, then will use either
#'   Barnes Hut or FIt-SNE based on nbody_algo. If Barnes Hut, then this
#'   determines the accuracy of BH approximation. Set to 0.5 by default.
#' @param max_iter number of iterations of t-SNE to run, set to 750 by default.
#' @param fft_not_bh if theta is nonzero, this determines whether to use FIt-SNE
#'   or Barnes Hut approximation, set to TRUE by default for FIt-SNE.
#' @param ann_not_vptree use vp-trees (as in bhtsne) or approximate nearest
#'   neighbors (default). Set to be TRUE for approximate nearest neighbours.
#' @param exaggeration_factor coefficient for early exaggeration (>1), set to 12
#'   by default.
#' @param no_momentum_during_exag set to 0 to use momentum and other
#'   optimization tricks. Can be set to 1 to do plain, vanilla gradient descent
#'   (useful for testing large exaggeration coefficients).
#' @param stop_early_exag_iter when to switch off early exaggeration, set to 250
#'   by default.
#' @param start_late_exag_iter when to start late exaggeration, set to -1 by
#'   default to not use late exaggeration.
#' @param late_exag_coeff late exaggeration coefficient, set to 1 by default to
#'   not use late exaggeration.
#' @param mom_switch_iter iteration number to switch from momentum to
#'   final_momentum, set to 250 by default.
#' @param momentum initial value of momentum, set to 0.5 by default.
#' @param final_momentum value of momentum to use later in the optimisation, set
#'   to 0.8 by default.
#' @param learning_rate set to desired learning rate or 'auto', which sets
#'   learning rate to N/exaggeration_factor where N is the sample size, or to
#'   200 if N/exaggeration_factor < 200.
#' @param n_trees when using Annoy, the number of search trees to use, set to 50
#'   by default.
#' @param search_k When using Annoy, the number of nodes to inspect during
#'   search. Default is 3*perplexity*n_trees (or K*n_trees when using fixed
#'   sigma).
#' @param rand_seed seed for random initialisation, set to -1 by default to
#'   initialise random number generator with current time.
#' @param nterms if using FIt-SNE, this is the number of interpolation points
#'   per sub-interval, set to 3 by default.
#' @param intervals_per_integer see min_num_intervals.
#' @param min_num_intervals let maxloc = ceil(max(max(X))) and minloc =
#'   floor(min(min(X))). i.e. the points are in a [minloc]^no_dims by
#'   [maxloc]^no_dims interval/square. The number of intervals in each dimension
#'   is either min_num_intervals or ceil((maxloc -
#'   minloc)/intervals_per_integer), whichever is larger. min_num_intervals must
#'   be an integer >0, and intervals_per_integer must be >0. Defaults are
#'   min_num_intervals=50 and intervals_per_integer = 1.
#' @param sigma fixed sigma value to use when perplexity==-1, set to -30 by
#'   default.
#' @param K number of nearest neighbours to get when using fixed sigma, set to
#'   -1 by default.
#' @param initialization 'pca', 'random', or N x no_dims array to intialize the
#'   solution, set to 'pca' by default.
#' @param max_step_norm maximum distance that a point is allowed to move on one
#'   iteration. Larger steps are clipped to this value. This prevents possible
#'   instabilities during gradient descent. Set to -1 to switch it off.
#'   Set to 5 by default.
#' @param load_affinities if 1, input similarities are loaded from a file and
#'   not computed. If 2, input similarities are saved into a file. If 0,
#'   affinities are neither saved nor loaded.
#' @param fast_tsne_path path to FItSNE executable.
#' @param nthreads number of threads to use, set to use all available threads by
#'   default.
#' @param perplexity_list if perplexity==0 then perplexity combination will be
#'   used with values taken from perplexity_list. Default: NULL df - Degree of
#'   freedom of t-distribution, must be greater than 0. Values smaller than 1
#'   correspond to heavier tails, which can often resolve substructure in the
#'   embedding. See Kobak et al. (2019) for details. Default is 1.0.
#' @param get_costs logical indicating whether the KL-divergence costs computed
#'   every 50 iterations should be returned, set to FALSE by default.
#' @param df positive numeric that controls the degree of freedom of
#'   t-distribution. The actual degree of freedom is 2*df-1. The standard t-SNE
#'   choice of 1 degree of freedom corresponds to df=1. Large df approximates
#'   Gaussian kernel. df<1 corresponds to heavier tails, which can often resolve
#'   substructure in the embedding. See Kobak et al. (2019) for details. Default
#'   is 1.0.
#'
#' @importFrom utils file_test
#' @importFrom stats sd
#'
#' @seealso \code{\link{cyto_map}}
#'
#' @references Linderman, G., Rachh, M., Hoskins, J., Steinerberger, S.,
#'   Kluger., Y. (2019). Fast interpolation-based t-SNE for improved
#'   visualization of single-cell RNA-seq data. Nature Methods.
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6402590/}.
#'
#' @export
fftRtsne <- function(X, 
                     dims = 2, 
                     perplexity = 30, 
                     theta = 0.5,
                     max_iter = 750,
                     fft_not_bh = TRUE,
                     ann_not_vptree = TRUE,
                     stop_early_exag_iter = 250,
                     exaggeration_factor = 12.0, 
                     no_momentum_during_exag = FALSE,
                     start_late_exag_iter = -1, 
                     late_exag_coeff = 1.0,
                     mom_switch_iter = 250, 
                     momentum = 0.5, 
                     final_momentum = 0.8, 
                     learning_rate = 'auto',
                     n_trees = 50, 
                     search_k = -1, 
                     rand_seed = -1,
                     nterms = 3, 
                     intervals_per_integer = 1, 
                     min_num_intervals = 50, 
                     K = -1, 
                     sigma = -30, 
                     initialization = 'pca',
                     max_step_norm = 5,
                     load_affinities = NULL,
                     fast_tsne_path = NULL, 
                     nthreads = 0, 
                     perplexity_list = NULL, 
                     get_costs = FALSE, df = 1.0) {
  
  version_number <- '1.2.1'
  
  if (is.null(fast_tsne_path)) {
    if(.Platform$OS.type == "windows"){
      fast_tsne_path <- SysExec(progs = ifelse(
        test = .Platform$OS.type == "windows",
        yes = "FItSNE.exe",
        no = "fast_tsne"
      ))
      if (length(fast_tsne_path) == 0) {
        stop("fast_tsne binary is not in the search path.")
      }
    }else{
      # RStudio Server look in home/rstudio/ directory (hidden)
      if(dir.exists("/home/rstudio")){
        if(".FIt-SNE" %in% list.files("/home/rstudio", all.files = TRUE)){
          fast_tsne_path <- "/home/rstudio/.FIt-SNE/bin/fast_tsne"
        }
      }else{
        fast_tsne_path <- system2('which', 'fast_tsne', stdout = TRUE)
      }
    }
  }
  fast_tsne_path <- normalizePath(fast_tsne_path)
  if (!file_test('-x', fast_tsne_path)) {
    stop(fast_tsne_path, " does not exist or is not executable.")
  }  
  
  data_path <- tempfile(pattern = 'fftRtsne_data_', fileext = '.dat')
  result_path <- tempfile(pattern = 'fftRtsne_result_', fileext = '.dat')
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  
  if (!is.numeric(theta) || (theta < 0.0) || (theta > 1.0) ) { 
    stop("Incorrect theta.")
  }
  if (nrow(X) - 1 < 3 * perplexity) { 
    stop("Perplexity is too large.")
  }
  if (!is.matrix(X)) { 
    stop("Input X is not a matrix")
  }
  if (!(max_iter > 0)) { 
    stop("Incorrect number of iterations.")
  }
  if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter < 0) {
    stop("stop_early_exag_iter should be a positive integer")
  }
  if (!is.numeric(exaggeration_factor)) {
    stop("exaggeration_factor should be numeric")
  }
  if (!is.numeric(df)) { 
    stop("df should be numeric")
  }
  if (!is.wholenumber(dims) || dims <= 0) { 
    stop("Incorrect dimensionality.")
  }
  if (search_k == -1) {
    if (perplexity > 0) {
      search_k <- n_trees * perplexity * 3
    } else if (perplexity == 0) {
      search_k <- n_trees * max(perplexity_list) * 3
    } else { 
      search_k <- n_trees * K
    }
  }
  
  if (is.character(learning_rate) && learning_rate =='auto') {
    learning_rate = max(200, nrow(X)/exaggeration_factor)
  }
  if (is.character(start_late_exag_iter) && start_late_exag_iter =='auto') {
    if (late_exag_coeff > 0) {
      start_late_exag_iter = stop_early_exag_iter
    }else {
      start_late_exag_iter = -1
    }
  }
  
  if (is.character(initialization) && initialization =='pca') {
    if (rand_seed != -1)  {
      set.seed(rand_seed)
    }
    # RSVD
    cyto_require(
      "rsvd",
      source = "CRAN",
      ref = paste0(
        "Erichson NB, Voronin S, Brunton SL, Kutz JN (2019). Randomized ",
        "Matrix Decompositions Using R. Journal of Statistical Software, ",
        "89(11), 1–48. doi: 10.18637/jss.v089.i11."
      )
    )
    message('Using rsvd() to compute the top PCs for initialization.')
    X_c <- scale(X, center=T, scale=F)
    rsvd_out <- cyto_func_call("rsvd::rsvd", list(A = X_c, k = dims))
    X_top_pcs <- rsvd_out$u %*% diag(rsvd_out$d, nrow=dims)
    initialization <- 0.0001*(X_top_pcs/sd(X_top_pcs[,1])) 
  }else if (is.character(initialization) && initialization == 'random'){
    message('Random initialization')
    initialization = NULL
  }
  
  if (fft_not_bh) {
    nbody_algo <- 2
  } else {
    nbody_algo <- 1
  }
  
  if (is.null(load_affinities)) {
    load_affinities <- 0
  } else {
    if (load_affinities == 'load') {
      load_affinities <- 1
    } else if (load_affinities == 'save') {
      load_affinities <- 2
    } else {
      load_affinities <- 0
    }
  }
  
  if (ann_not_vptree) {
    knn_algo <- 1
  } else {
    knn_algo <- 2
  }
  tX <- as.numeric(t(X))
  
  f <- file(data_path, "wb")
  n <- nrow(X)
  D <- ncol(X)
  writeBin(as.integer(n), f, size = 4)
  writeBin(as.integer(D), f, size = 4)
  writeBin(as.numeric(theta), f, size = 8) #theta
  writeBin(as.numeric(perplexity), f, size = 8)
  
  if (perplexity == 0) {
    writeBin(as.integer(length(perplexity_list)), f, size = 4)
    writeBin(perplexity_list, f) 
  }
  
  writeBin(as.integer(dims), f, size = 4)
  writeBin(as.integer(max_iter), f, size = 4)
  writeBin(as.integer(stop_early_exag_iter), f, size = 4)
  writeBin(as.integer(mom_switch_iter), f, size = 4)
  writeBin(as.numeric(momentum), f, size = 8)
  writeBin(as.numeric(final_momentum), f, size = 8)
  writeBin(as.numeric(learning_rate), f, size = 8)
  writeBin(as.numeric(max_step_norm), f, size = 8)
  writeBin(as.integer(K), f, size = 4) #K
  writeBin(as.numeric(sigma), f, size = 8) #sigma
  writeBin(as.integer(nbody_algo), f, size = 4)  #not barnes hut
  writeBin(as.integer(knn_algo), f, size = 4) 
  writeBin(as.numeric(exaggeration_factor), f, size = 8) #compexag
  writeBin(as.integer(no_momentum_during_exag), f, size = 4) 
  writeBin(as.integer(n_trees), f, size = 4) 
  writeBin(as.integer(search_k), f, size = 4) 
  writeBin(as.integer(start_late_exag_iter), f, size = 4) 
  writeBin(as.numeric(late_exag_coeff), f, size = 8) 
  
  writeBin(as.integer(nterms), f, size = 4) 
  writeBin(as.numeric(intervals_per_integer), f, size = 8) 
  writeBin(as.integer(min_num_intervals), f, size = 4) 
  writeBin(tX, f) 
  writeBin(as.integer(rand_seed), f, size = 4) 
  writeBin(as.numeric(df), f, size = 8)
  writeBin(as.integer(load_affinities), f, size = 4) 
  if (!is.null(initialization) ) { writeBin( c(t(initialization)), f) }		
  close(f) 
  
  flag <- system2(command = fast_tsne_path, 
                  args = c(version_number, data_path, result_path, nthreads))
  if (flag != 0) {
    stop('tsne call failed')
  }
  f <- file(result_path, "rb")
  n <- readBin(f, integer(), n = 1, size = 4)
  d <- readBin(f, integer(), n = 1, size = 4)
  Y <- readBin(f, numeric(), n = n * d)
  Y <- t(matrix(Y, nrow = d))
  if (get_costs) {
    readBin(f, integer(), n = 1, size = 4)
    costs <- readBin(f, numeric(), n = max_iter, size = 8)
    Yout <- list(Y = Y, costs = costs)
  } else {
    Yout <- Y
  }
  close(f)
  file.remove(data_path)
  file.remove(result_path)
  Yout
}