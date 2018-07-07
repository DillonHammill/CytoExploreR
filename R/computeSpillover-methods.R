#' Compute Spillover Matrix
#' 
#' Compute spillover matrix using single stain compensation controls and an unstained control.
#' 
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @references C. B. Bagwell \& E. G. Adams (1993). Fluorescence spectral overlap compensation for any number of flow cytometry parameters. in: Annals of the New York Academy of Sciences, 677:167-184.
#' 
#' @export
setGeneric(name="computeSpillover",
           def=function(x, ...){standardGeneric("computeSpillover")}
)

#' Caluclate Spillover Matrix Using Pre-gated Compensation Controls
#' 
#' @param x object of class \code{flowSet} containing pre-gated compensation single stain controls as well as an unstained control.
#' @param alias name of the gated population to use for downstream calculations, set to the last node of the GatingSet by default.
#' @param pdfile \code{pData} csv file containing additional column \code{"channel"} indicating the fluorescent channel associated
#' with each sample. This channel should be set to \code{"Unstained"} for unstained controls.
#' @param spfile name of the output spillover csv file, set to \code{"Spillover Matrix.csv"} by default.
#' 
#' @return spillover matrix object and \code{"Spillover Matrix.csv"} file as side effect.
#' 
#' @export
setMethod(computeSpillover, signature = "flowSet", definition = function(x, pdfile = NULL, spfile = NULL, ...){
  
  fs <- x
  
  # Extract fluorescent channels
  channels <- getChannels(fs)
  
  # Select a fluorescent channel for each compensation control
  if(is.null(pdfile)){
    
    pData(fs)$channel <- selectChannels(fs)
    write.csv(pData(fs), "Compensation Control Channels.csv", row.names = FALSE)
    
  }else{
    
    pData(fs) <- read.csv(pdfile, header = TRUE, row.names = 1)
    
  }
  
  # Transform fluorescent channels
  trans <- estimateLogicle(fs[[1]], channels)
  fs <- transform(fs, trans)
  
  # Extract unstained control based on selected channels in pData(fs)
  NIL <- fs[[match("Unstained", pData(fs)$channel)]]
  fs <- fs[-match("Unstained", pData(fs)$channel)]
  
  # Assign channel to each flowFrame to description slot called "Channel"
  sapply(1:length(pData(fs)$channel), function(x){
    fs[[x]]@description$Channel <- paste(pData(fs)$channel[x])
  })
  
  # Gate positive populations
  pops <- fsApply(fs, function(fr){
    
    # Call drawGate on each flowFrame using interval gate on selected channel
    gt <- drawGate(x = fr, alias = paste(fr@description$Channel,"+"), channels = fr@description$Channel, gate_type = "interval", adjust = 1.5)
    fr <- Subset(fr, gt[[1]])
    
  }, simplify = TRUE)
  
  # Inverse logicle transformation
  inv <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
  pops <- transform(pops, inv)
  NIL <- transform(NIL, inv)
  
  # Calculate MedFI for all channels for unstained control
  neg <- each_col(NIL, median)[channels]
  
  # Calculate MedFI for all channels for all stained controls
  pos <- fsApply(pops, each_col, median)[,channels]
  
  # Subtract background fluorescence
  signal <- sweep(pos, 2, neg)
  
  # Construct spillover matrix - only include values for which there is a control
  spill <- diag(x = 1, nrow = length(channels), ncol = length(channels))  
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  # Normalise each row to stained channel
  for(i in 1:nrow(signal)){
    
    signal[i, ] <- signal[i, ]/signal[i, match(fs[[i]]@description$Channel, colnames(spill))]
  } 
  
  # Insert values into appropriate rows
  rws <- match(pData(fs)$channel, rownames(spill))
  spill[rws,] <- signal
  
  write.csv(spill, "Spillover Matrix.csv", row.names = FALSE)
  return(spill)
  
})

#' Caluclate Spillover Matrix - flowSet Method
#' 
#' @param x object of class \code{GatingSet} containing gated compensation single stain controls as well as an unstained control.
#' @param alias name of the gated population to use for downstream calculations, set to the last node of the GatingSet by default.
#' @param pdfile \code{pData} csv file containing additional column \code{"channel"} indicating the fluorescent channel associated
#' with each sample. This channel should be set to \code{"Unstained"} for unstained controls.
#' @param spfile name of the output spillover csv file, set to \code{"Spillover Matrix.csv"} by default.
#' 
#' @importFrom flowWorkspace getData
#' 
#' @return spillover matrix and \code{"Spillover Matrix.csv"} file.
#' 
#' @export
setMethod(computeSpillover, signature = "GatingSet", definition = function(x, alias = NULL, pdfile = NULL, spfile = NULL, ...){
  
  gs <- x
  
  # Extract fluorescent channels
  channels <- getChannels(gs)
  
  # Select a fluorescent channel for each compensation control
  if(is.null(pdfile)){
    
    pData(gs)$channel <- paste(selectChannels(gs))
    write.csv(pData(gs), "Compensation Controls Channels.csv", row.names = FALSE)
    
  }else{
    
    pd <- read.csv(pdfile, header = TRUE, row.names = 1)
    pData(gs)$channel <- paste(pd$channel)
    
  }
  
  # Check if fluorescent channels have been transformed
  if(length(gs@transformation) == 0){
    
    # No transformation has been applied - get transform relevant channel in each sample
    indx <- seq(1:length(pData(gs)$channel))
    if(!is.na(match("Unstained", pData(gs)$channel))){
      
      indx <- indx[-match("Unstained", pData(gs)$channel)]
      
    }
    
    translist <- lapply(indx, FUN = function(x){
      
      estimateLogicle(gs[[x]], pData(gs)$channel[x])[[pData(gs)$channel[x]]]
      
    })
    trans <- transformerList(pData(gs)$channel[indx], translist)
    gs <- transform(gs, trans)
    
  }else if(length(gs@transformation) != 0){
    
    # Check which channels have been transformed
    chans <- names(gs@transformation[[1]])
    
    if(all(channels %in% chans)){
      
      # All fluorescent channels have been transformed
      
    }else{
      
      # Not all fluorescent channels are transformed - transform remaining channels
      # Which channels need transform? Use appropriate sample to get transform (indx of samples)
      chans <- names(gs@transformation[[1]])
      chans <- channels[is.na(match(channels,chans))]
      
      indx <- na.omit(match(pData(gs)$channel, chans))
      translist <- lapply(indx, FUN = function(x){
        
        estimateLogicle(gs[[x]], pData(gs)$channel[x])[[pData(gs)$channel[x]]]
        
      })
      trans <- transformerList(pData(gs)$channel[indx], translist)
      gs <- transform(gs, trans)
      
    }
    
  }
  
  # Extract Population for Downstream Analyses
  if(!is.null(alias)){
    
    fs <- flowWorkspace::getData(gs, alias)
    
  }else if(is.null(alias)){
    
    fs <- flowWorkspace::getData(gs, getNodes(gs)[length(getNodes(gs))])
    
  }
  
  # Extract unstained control based on selected channels in pData(fs)
  NIL <- fs[[match("Unstained", pData(fs)$channel)]]
  fs <- fs[-match("Unstained", pData(fs)$channel)]
  
  # Assign channel to each flowFrame to description slot called "channel"
  sapply(1:length(pData(fs)$channel), function(x){
    fs[[x]]@description$channel <- paste(pData(fs)$channel[x])
  })
  
  # Gate positive populations
  pops <- fsApply(fs, function(fr){
    
    # Call drawGate on each flowFrame using interval gate on selected channel
    gt <- drawGate(x = fr, alias = paste(fr@description$channel,"+"), channels = fr@description$channel, gate_type = "interval", adjust = 1.5)
    fr <- Subset(fr, gt[[1]])
    
  }, simplify = TRUE)
  
  # Inverse logicle transformation
  inv <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
  pops <- transform(pops, inv)
  NIL <- transform(NIL, inv)
  
  # Calculate MedFI for all channels for unstained control
  neg <- each_col(NIL, median)[channels]
  
  # Calculate MedFI for all channels for all stained controls
  pos <- fsApply(pops, each_col, median)[,channels]
  
  # Subtract background fluorescence
  signal <- sweep(pos, 2, neg)
  
  # Construct spillover matrix - only include values for which there is a control
  spill <- diag(x = 1, nrow = length(channels), ncol = length(channels))  
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  # Normalise each row to stained channel
  for(i in 1:nrow(signal)){
    
    signal[i, ] <- signal[i, ]/signal[i, match(fs[[i]]@description$channel, colnames(spill))]
  } 
  
  # Insert values into appropriate rows
  rws <- match(pData(fs)$channel, rownames(spill))
  spill[rws,] <- signal
  
  write.csv(spill, "Spillover Matrix.csv", row.names = FALSE)
  return(spill)
  
})