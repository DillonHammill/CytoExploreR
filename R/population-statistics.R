#' Compute Population Level Statistics for a GatingSet.
#'
#' \code{computePopStats} function allows easy extraction of population level statistics from a \code{GatingSet} object,
#' including mean, geometirc mean and median fluorescent intensities as well as population counts and frequencies. 
#' \code{PopStats} searches \code{GatingSet} object for logicle transformation parameters and if present performs the inverse
#' logicle transformation prior to calculation of relevant statistic. Users can also supply their own transformerList to the
#' "trans" argument which will be used for the inverse transformation. By default if no channels are supplied \code{computePopStats} 
#' returns statistics for channels containing \code{markernames}, if no \code{markernames} are found statistics are performed on all 
#' channels.
#'
#' @param gs a \code{GatingSet} object defining populations for which statistics should be calulated.
#' @param parent name of the parent population to use for frequency calculations.
#' @param alias a vector containing the \code{node} names of the populations of interest.
#' @param stat specify the statistic of interest by name, can be either "mean", "median" or "count". By default we use the more robust
#' median fluorescent intensity (MedFI).
#' @param channels a vector of channel names for which population statistics should be calculated, all fluorescent channels by default.
#' @param trans object of class transformerList to be used for inverse logicle transformation.
#'
#' @return a \code{list} object containing the calculated statistics for specified populations in all samples.
#'
#' @keywords population, statistics, mean, median, geometric mean, MFI, count, frequency
#' 
#' @importFrom BiocGenerics colnames
#' 
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
computePopStats <- function(gs, parent = NULL, alias = NULL, stat = "median", channels = NULL, inv.trans = TRUE, trans = NULL){
  
  # gs is not a GatingSet
  if(class(gs) != "GatingSet"){
    
    stop("computePopStats can only be applied to an object of class GatingSet.")
    
  }
  
  # No alias supplied
  if(is.null(alias)){
    
    stop("Please supply the names of the population(s) as the alias argument.")
    
  }
  
  # No channels supplied use all fluorescent channels
  if(is.null(channels)){
    
    channels <- getChannels(gs)
    
  }
  
  # Extract relevant function for calculation
  if(stat %in% c("mean","median")){
    
    func <- match.fun(stat)
    
  }else if(stat == "count"){
    
    func <- match.fun("length")
    channels <- c("FSC-A")
    inv.trans <- FALSE
    
  }else if(stat == "mode"){
    
    func <- function(x) {density(x)$x[which.max(density(x, adjust = 1.5)$y)]}
    
  }else if(stat == "geomean"){
    
    func <- function(x) {exp(mean(log(x)))}
    
  }else if(stat == "freq"){
    
    func <- match.fun("length")
    channels <- c("FSC-A")
    inv.trans <- FALSE
    
  }else{
    
    func <- match.fun(stat)
    inv.trans <- FALSE
    
  }
  
  # For stat == "freq" get counts for parent and all pops - later divide alias by parent * 100
  if(stat == "freq"){
    
    alias <- c(parent, alias)
    
  }
  
  # Extract populations from GatingSet and store results for each population as separate elements of a list
  results <- list()  
  
  for (pop in alias) {
    popData <- getData(gs, pop)
    
    # Inverse logicle transformation if required
    if(inv.trans == TRUE){
      
      if(is.null(trans)){
        
        if(is.null(gs[[1]]@transformation)){
          
          stop("Please supply a transformerList object to 'trans' argument to perform inverse transformations.")
          
        }
        
        trans <- gs[[1]]@transformation
        inv <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
        popData <- transform(popData, inv)
        
      }else if(!is.null(trans)){
        
        inv <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
        popData <- transform(popData, inv)
        
      }
      
    }else{
      
      # No inverse transformation required
      
    }
    
    # Calculate statistic for each supplied channel by name or inidices
    MedFI <- fsApply(popData, function(fr) apply(exprs(fr[, channels]), 2, func))
    results[[pop]] <- MedFI
    
  }
  
  # if stat == "freq" divide each alias by parent *100
  if(stat == "freq"){
    
    # Repeat columns 1 for each parent
    results <- lapply(results, function(x){
      
      x <- matrix(x[,1], nrow = length(fs), ncol = length(parent))
      colnames(x) <- parent
      rownames(x) <- sampleNames(fs)
      return(x)
      
    })
    
    # Parents
    prnts <- matrix(nrow = length(fs), ncol = length(parent))
    colnames(test) <- parent
    rownames(test) <- sampleNames(fs)
    for(i in 1:length(parent)){
      
      prnts[,i] <- results[[i]][,1]
      
    }
    
    # Divide each column by parent col1 by parent1 etc.
    for(i in 1:length(results)){
      
      results[[i]] <- (results[[i]] / prnts)*100
      
    }
    
    results <- results[-c(1:length(parent))]
  }
  
  # Change colnames for count data to be pop.Count
  if(stat == "count"){
    
    results <- do.call("cbind",results)
    titles <- c()
    
    for(pop in pops){
      
      titles <- c(titles,paste(pop,"Count", sep = "."))
      
    }
    
    base::colnames(results) <- titles
    
  }
  
  return(results)
  
}