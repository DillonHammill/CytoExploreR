library(CytoRSuite)
library(CytoRSuiteData)

# Needed to bypass wd checks for files
options("CytoRSuite_wd_check" = FALSE)

# Data directory
datadir <- system.file("extdata", package = "CytoRSuiteData")

# Sample for speed -
Compensation <- read.flowSet(path = paste0(datadir,"/Compensation-Controls"), pattern = ".fcs")

Comp <- flowSet(lapply(1:length(Compensation), function(x) {
  Compensation[[x]][1:1000, ]
}))
sampleNames(Comp) <- sampleNames(Compensation)
pData(Comp)$name <- sampleNames(Compensation)

# GatingSet -
gsc <- GatingSet(Comp)

# gatingTemplate -
gtc <- gatingTemplate(paste0(datadir,"/Compensation-gatingTemplate.csv"))

# Gating -
gating(gtc, gsc)

# Spillover
spfile <- system.file("extdata", "Ref-Spillover-matrix.csv", package = "CytoRSuite")

# Channel match
cmfile <- system.file("extdata", "Compensation-channels.csv", package = "CytoRSuite")
  
spillover_edit(gsc,
               parent = "Single Cells",
               spillover = spfile,
               channel_match = cmfile)