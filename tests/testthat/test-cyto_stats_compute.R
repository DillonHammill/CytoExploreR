context("cyto_stats_compute")

# Prepare data for testing -----------------------------------------------------------------------

# Inverse transformed data -
inv <- cyto_trans_check(trans, inverse = TRUE)
V <- transform(Va2, inv)

# cyto_stats_compute flowFrame method ------------------------------------------------------------

test_that("cyto_stats_compute flowFrame method", {
  
  # Median -
  sts <- matrix(colMedians(exprs(Va2[[1]])), nrow = 1, ncol = length(chans))
  rownames(sts) <- nms[[1]]
  colnames(sts) <- chans
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "median"), sts, tolerance = 0.01)
  
  sts <- matrix(colMedians(exprs(Va2[[1]])), nrow = 1, ncol = length(chans))
  colnames(sts) <- chans
  sts <- matrix(sts[, c("Alexa Fluor 700-A","Alexa Fluor 488-A")], nrow = 1)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- c("Alexa Fluor 700-A","Alexa Fluor 488-A")
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "median", channels = c("CD4","CD8")), sts, tolerance = 0.01)
  
  sts <- matrix(colMedians(exprs(V[[1]])), nrow = 1, ncol = length(chans))
  rownames(sts) <- nms[[1]]
  colnames(sts) <- chans
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "median", trans = trans), sts, tolerance = 0.01)
  
  # Mean -
  sts <- matrix(colMeans(exprs(Va2[[1]])), nrow = 1, ncol = length(chans))
  rownames(sts) <- nms[[1]]
  colnames(sts) <- chans
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "mean"), sts, tolerance = 0.01)
  
  sts <- matrix(colMeans(exprs(Va2[[1]])), nrow = 1, ncol = length(chans))
  colnames(sts) <- chans
  sts <- matrix(sts[, c("Alexa Fluor 700-A","Alexa Fluor 488-A")], nrow = 1)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- c("Alexa Fluor 700-A","Alexa Fluor 488-A")
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "mean", channels = c("CD4","CD8")), sts, tolerance = 0.01)
  
  sts <- matrix(colMeans(exprs(V[[1]])), nrow = 1, ncol = length(chans))
  rownames(sts) <- nms[[1]]
  colnames(sts) <- chans
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "mean", trans = trans), sts, tolerance = 0.01)
  
  # Geometric mean -
  sts <- lapply(chans, function(x){
    
    if(!x %in% names(trans)){
      
      exp(mean(log(exprs(Va2[[1]])[, x])))
      
    }else{
      
      inv@transforms[[x]]@f(mean(exprs(Va2[[1]])[,x]))
      
    }
  })
  sts <- do.call("cbind", sts)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- chans
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "geo mean", trans = trans), sts, tolerance = 0.01)
  
  sts <- suppressWarnings(lapply(chans, function(x){
    
    exp(mean(log(exprs(Va2[[1]])[, x])))
    
  }))
  sts <- do.call("cbind", sts)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- chans
  
  expect_error(cyto_stats_compute(Va2[[1]], stat = "geo mean"),"Some channels have been transformed, supply the transformList to calculate the geometric mean.")
  
  # Mode -
  sts <- lapply(chans, function(x){
    d <- density(exprs(Va2[[1]])[,x], adjust = 1.5)
    d$x[d$y == max(d$y)]
  })
  sts <- do.call("cbind", sts)
  colnames(sts) <- chans
  rownames(sts) <- nms[[1]]
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "mode"), sts, tolerance = 0.01)
  
  sts <- matrix(sts[, c("Alexa Fluor 700-A", "Alexa Fluor 488-A")], nrow = 1)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "mode", channels = c("CD4","CD8")), sts, tolerance = 0.01)
  
  sts <- lapply(chans, function(x){
    d <- density(exprs(V[[1]])[,x], adjust = 1.5)
    d$x[d$y == max(d$y)]
  })
  sts <- do.call("cbind", sts)
  colnames(sts) <- chans
  rownames(sts) <- nms[[1]]
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "mode", trans = trans), sts, tolerance = 0.01)
  
  # CV -
  sts <- lapply(chans, function(x){
    (sd(exprs(Va2[[1]])[, x])/mean(exprs(Va2[[1]])[, x]))*100
  })
  sts <- do.call("cbind", sts)
  colnames(sts) <- chans
  rownames(sts) <- nms[[1]]
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "CV"), sts, tolerance = 0.01)
  
  sts <- matrix(sts[, c("Alexa Fluor 700-A", "Alexa Fluor 488-A")], nrow = 1)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "CV", channels = c("CD4","CD8")), sts, tolerance = 0.01)
  
  sts <- lapply(chans, function(x){
    (sd(exprs(V[[1]])[, x])/mean(exprs(V[[1]])[, x]))*100
  })
  sts <- do.call("cbind", sts)
  colnames(sts) <- chans
  rownames(sts) <- nms[[1]]
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "CV", trans = trans), sts, tolerance = 0.01)
  
  # CVI -
  sts <- lapply(chans, function(x){
    1/(sd(exprs(Va2[[1]])[, x])/mean(exprs(Va2[[1]])[, x]))*100
  })
  sts <- do.call("cbind", sts)
  colnames(sts) <- chans
  rownames(sts) <- nms[[1]]
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "CVI"), sts, tolerance = 0.01)
  
  sts <- matrix(sts[, c("Alexa Fluor 700-A", "Alexa Fluor 488-A")], nrow = 1)
  rownames(sts) <- nms[[1]]
  colnames(sts) <- c("Alexa Fluor 700-A", "Alexa Fluor 488-A")
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "CVI", channels = c("CD4","CD8")), sts, tolerance = 0.01)
  
  sts <- lapply(chans, function(x){
    1/(sd(exprs(V[[1]])[, x])/mean(exprs(V[[1]])[, x]))*100
  })
  sts <- do.call("cbind", sts)
  colnames(sts) <- chans
  rownames(sts) <- nms[[1]]
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "CVI", trans = trans), sts, tolerance = 0.01)
  
  # Count -
  sts <- data.frame(nrow(Va2[[1]]))
  rownames(sts) <- nms[1]
  colnames(sts) <- "count"
  
  expect_equal(cyto_stats_compute(Va2[[1]], stat = "count"), sts, tolerance = 0.01)
  
  
})

# cyto_stats_compute flowSet method --------------------------------------------------------------

test_that("cyto_stats_compute flowSet method returns the correct statistics",{
  
  sts <- fsApply(V, nrow, use.exprs = TRUE)
  colnames(sts) <- "count"
  sts <- cbind(pData(Va2), sts)
  sts <- sts[,-1]
  
  expect_equal(cyto_stats_compute(Va2, stat = "count", trans = trans)[[1]], sts, tolerance = 0.01)
  
  sts <- fsApply(V, colMedians, use.exprs = TRUE)
  sts <- cbind(pData(Va2), sts)
  sts <- sts[,-1]
  
  expect_equal(cyto_stats_compute(Va2, stat = "median", trans = trans)[[1]], sts, tolerance = 0.01)
  
})

# cyto_stats_compute GatingSet method ------------------------------------------------------------

test_that("cyto_stats_compute GatingSet method returns the correct statistics",{
  
  sts <- fsApply(V, colMedians, use.exprs = TRUE)
  sts <- cbind(pData(Va2), sts)
  sts <- sts[,-1]
  
  expect_equal(cyto_stats_compute(gs, alias = "T Cells", stat = "median", trans = trans, save = FALSE)[[1]], sts)
  
  expect_equal(cyto_stats_compute(gs, alias = "T Cells", stat = "median", save = FALSE)[[1]], sts)
  
  pops <- c("CD4 T Cells", "CD8 T Cells", "root", "Live Cells")
  cnts <- lapply(pops, function(pop){
    
    fsApply(getData(gs, pop), nrow)
    
  })
  CD4 <- do.call("cbind", cnts[c(1,3,4)])
  CD4[,2] <- (CD4[,1]/CD4[,2])*100
  CD4[,3] <- (CD4[,1]/CD4[,3])*100
  CD4 <- cbind(pData(Va2), CD4)
  CD4 <- CD4[,-1]
  CD4 <- data.frame(CD4)
  colnames(CD4) <- c("OVAConc", "Treatment", "count","root","Live Cells")
  
  CD8 <- do.call("cbind", cnts[c(2,3,4)])
  CD8[,2] <- (CD8[,1]/CD8[,2])*100
  CD8[,3] <- (CD8[,1]/CD8[,3])*100
  CD8 <- cbind(pData(Va2), CD8)
  CD8 <- CD8[,-1]
  CD8 <- data.frame(CD8)
  colnames(CD8) <- c("OVAConc", "Treatment", "count","root","Live Cells")
  
  sts <- list(CD4, CD8)
  names(sts) <- c("CD4 T Cells","CD8 T Cells")

  expect_equal(cyto_stats_compute(gs, alias = c("CD4 T Cells", "CD8 T Cells"), parent = c("root", "Live Cells"), stat = "freq"), sts, tolerance = 0.01)
  expect_true(.file_wd_check(paste(format(Sys.Date(), "%d%m%y"),"-CD4 T Cells-freq.csv", sep = "")))
  expect_true(.file_wd_check(paste(format(Sys.Date(), "%d%m%y"),"-CD8 T Cells-freq.csv", sep = "")))
  
})

base::unlink(paste(format(Sys.Date(), "%d%m%y"),"-CD4 T Cells-freq.csv", sep = ""))
base::unlink(paste(format(Sys.Date(), "%d%m%y"),"-CD8 T Cells-freq.csv", sep = ""))
