#' Caluclate Spillover Matrix
#' 
#' @param fs object of class \code{flowSet} containing compensation single stain controls as well as an unstained control.
#' 
#' @return spillover matrix and \code{"Spillover Matrix.csv"} file.
#' 
#' @export
computeSpillover <- function(fs, file = NULL){
  
  # Extract fluorescent channels
  channels <- getChannels(fs)
  
  # Select a fluorescent channel for each compensation control
  if(is.null(file)){
 
  pData(fs)$channel <- selectChannels(fs)
  write.csv(pData(fs), "Compensation Control Channels.csv", row.names = 1)
  
  }else{
    
    pData(fs) <- read.csv(file, header = TRUE)
    
  }
  
  # Gate compensation controls to get single cells
  gs <- GatingSet(fs)
  
  # Transform fluorescent channels
  trans <- estimateLogicle(gs[[1]], channels)
  gs <- transform(gs, trans)
  
  gs <- drawGate(gs, alias = "Cells", channels = c("FSC-A","SSC-A"), gate_type = "polygon", subSample = 250000,file = "Compensation gatingTemplate.csv")
  gs <- drawGate(gs, alias = "Single Cells", channels = c("SSC-W","SSC-H"), gate_type = "polygon", subSample = 250000,file = "Compensation gatingTemplate.csv")
  fs <- getData(gs, "Single Cells")
  
  # Extract unstained control based on selected channels in pData(fs)
  NIL <- fs[match("Unstained", pData(fs)$channel)]
  fs <- fs[-match("Unstained", pData(fs)$channel)]
  
  # Assign channel to each flowFrame to description slot called "Channel"
  sapply(1:length(pData(fs)$channel), function(x){
    fs[[x]]@description$Channel <- pData(fs)$channel[x]
  })
  
  # Gate positive populations
  pops <- fsApply(fs, function(fr){
    
    # Call drawGate on each flowFrame using interval gate on selected channel
    gt <- drawGate(fr, alias = paste(fr@description$Channel,"+"), channels = fr@description$Channel, gate_type = "interval")
    fr <- Subset(fr, gt[[1]])
    
  }, simplify = TRUE)
  
  # Inverse logicle transformation
  inv <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
  pops <- transform(pops, inv)
  NIL <- transform(NIL, inv)
  NIL <- NIL[[1]]
  
  # Calculate MedFI for all channels for unstained control
  neg <- each_col(NIL, median)[channels]
  print(neg)
  
  # Calculate MedFI for all channels for all stained controls
  pos <- fsApply(pops, each_col, median)[,channels]
  print(pos)
  
  # Subtract background fluorescence
  signal <- sweep(pos, 2, neg)
  print(signal)
  
  # Construct spillover matrix - only include values for which there is a control
  spill <- diag(x = 1, nrow = length(channels), ncol = length(channels))  
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  # Normalise each row to stained channel
  sapply(1:nrow(signal), function(x){
    
    signal[x, ] <- signal[x, ]/signal[x, match(fs[[x]]@description$Channel, colnames(spill))]
    
  })
  print(signal)
  
  # Insert values into appropriate rows
  rws <- match(pData(fs)$channel, rownames(spill))
  spill[rws,] <- signal
  
  write.csv(spill, "Spillover Matrix.csv", row.names = FALSE)
  return(spill)
  
}
