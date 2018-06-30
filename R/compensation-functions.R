#' Caluclate Spillover Matrix - flowSet Method
#' 
#' @param fs object of class \code{flowSet} containing compensation single stain controls as well as an unstained control.
#' @param gtfile openCyto \code{gatingTemplate} csv file to be used to pre-gate samples prior to spillover calculation. If no
#' \code{gatingTemplate} is supplied samples will be gated FSC-A/SSC-A singlets gated in SSC-W/SSC-H using \code{drawGate}. The gating
#' strategy will be saved in gatingTemplate csv file called \code{"Compensation gatingTemplate.csv"}.
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
computeSpillover <- function(fs, gtfile = NULL, alias = NULL, pdfile = NULL, spfile = NULL){
  
  # Extract fluorescent channels
  channels <- getChannels(fs)
  
  # Select a fluorescent channel for each compensation control
  if(is.null(pdfile)){
 
  pData(fs)$channel <- selectChannels(fs)
  write.csv(pData(fs), "Compensation Control Channels.csv", row.names = FALSE)
  
  }else{
    
    pData(fs) <- read.csv(pdfile, header = TRUE, row.names = 1)
    
  }
  
  # Gate compensation controls to get single cells
  gs <- GatingSet(fs)
  
  # Transform fluorescent channels
  trans <- estimateLogicle(gs[[1]], c(channels,"SSC-A","SSC-H","SSC-W"))
  gs <- transform(gs, trans)
  
  # Gate Samples
  if(!is.null(gtfile)){
    
    gt <- gatingTemplate(gtfile)
    gating(gt,gs)
    
  }else if(is.null(gtfile)){
    
    drawGate(gs, alias = "Cells", channels = c("FSC-A","SSC-A"), gate_type = "polygon", subSample = 250000,file = "Compensation gatingTemplate.csv")
    drawGate(gs, alias = "Single Cells", channels = c("SSC-W","SSC-H"), gate_type = "polygon", subSample = 250000,file = "Compensation gatingTemplate.csv")

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
  
}