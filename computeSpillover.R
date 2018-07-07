files <- list.files(path = "Single Stain Controls", full.names = TRUE)
fs <- read.ncdfFlowSet(files = files)

gs <- GatingSet(fs)
drawGate(gs, alias = "Cells", channels = c("FSC-A","SSC-A"), gate_type = "polygon")
drawGate(gs, parent = "Cells", alias = "Single Cells", channels = c("SSC-W","SSC-H"), gate_type = "polygon", template = "Template")


setMethod(computeSpillover, signature = "GatingSet", definition = function(x, gtfile = NULL, alias = NULL, pdfile = NULL, spfile = NULL, ...){
  
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
computeSpillover(gs, alias = "Single Cells", pdfile = "Compensation Control Channels.csv")
