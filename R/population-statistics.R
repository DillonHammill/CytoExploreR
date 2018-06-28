#' Extract Population Level Statistics from GatingSet.
#'
#' \code{PopStats} function allows easy extraction of population level statistics from a \code{GatingSet} object,
#' including geometric mean and median fluorescent intensities. \code{PopStats} searches \code{GatingSet} object for logicle
#' transformation parameters and if present performs the inverse logicle transformation prior to calculation of relevant
#' statistic. By default if no channels are supplied \code{PopStats} returns statistics for channels containing \code{markernames},
#' if no \code{markernames} are found statistics are performed on all fluorescent channels.
#'
#' @param gs a \code{GatingSet} object defining populations for which statistics should be calulated.
#' @param pops a vector containing the \code{node} names of the populations of interest.
#' @param stat specify the statistic of interest by name, can be either "mean", "median" or "count". By default we use the more robust
#' median fluorescent intensity (MedFI).
#' @param channels a vector of channel names for which population statistics should be calculated, all channels by default.
#'
#' @return a \code{list} object containing the calculated statistics for specified populations in all samples.
#'
#' @keywords population, statistics, median, geometric mean, MFI
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
getPopMFI <- function(gs, pops, stat = "median", channels = NULL){
  
  # If no markernames present use all fluorescent channels
  if(is.null(channels)){
    if(length(markernames(gs)) == 0){
      channels <- FlowDJ::fluor_channels(gs@data[[1]])
    } else {
      channels <- c(gs@data[[1]]@parameters$name[!is.na(gs@data[[1]]@parameters$desc)])
    }
  } else {
    
  }
  
  # Check whether logicle transformation has been applied to GatingSet
  if(length(gs@transformation) == 0){
    message("No transformation parameters found in GatingSet: inverse logicle transformation has not been applied.")
    inv_trans <- FALSE
  } else {
    message("Inverse logicle transformation has been applied.")
    inv_trans <- TRUE
  }
  
  # Extract relevant statistic for calculation
  if(stat %in% c("mean","median")){
    func <- match.fun(stat)
  } else if (stat == "count"){
    func <- match.fun("length")
    channels <- c("FSC-A")
  }
  
  # Extract populations from GatingSet and store results for each population as separate elements of a list
  results <- list()
  for (pop in pops) {
    popData <- getData(gs, pop)
    
    # Inverse logicle transformation if required
    if(inv_trans == TRUE){
      inv.trans <- transformList(names(gs@transformation), lapply(gs@transformation, `[[`, "inverse"))
      popData <- transform(popData, inv.trans)
    }
    
    # Calculate statistic for each supplied channel by name or inidices
    MedFI <- fsApply(popData, function(fr) apply(exprs(fr[, channels]), 2, func))
    results[[pop]] <- MedFI
  }
  
  # Change colnames for count data to be pop.Count
  if(stat == "count"){
    results <- do.call("cbind",results)
    titles <- c()
    for(pop in pops){
      titles <- c(titles,paste(pop,"Count", sep = "."))
    }
    colnames(results) <- titles
  }
  
  return(results)
}

#' Get Population Frequency
#' 
#' @param gs an object of class \code{GatingSet}.
#' @param alias name of the population.
#' @param parent name of the parent population to use for calculation
#' 
#' @importFrom flowCore fsApply
#' @importFrom flowCore nrow
#' 
#' @export
getFreq <- function(gs, alias, parent, type = "percent"){
  
  # Extract Population Events
  pop <- getData(gs, alias)
  pop.events <- fsApply(pop, function(fr){
    nrow(fr)
  })
  
  # Extract Parent Population
  prnt <- getData(gs, parent)
  prnt.events <- fsApply(prnt, function(fr){
    nrow(fr)
  })
  
  # Calculate Percentage
  prcnt <- (pop.events/prnt.events)*100
  
  return(prcnt)
  
}
