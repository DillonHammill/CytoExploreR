# Load Packages ----------------------------------------------------------------
library(CytoExploreRData)

# Load in vdiffr for image comparison -
library(vdiffr)

# Load in grDevices for recordPlot() -
library(grDevices)

# Load in robustbase (colMedians) -
library(robustbase)

# Libraries for data manipulation
library(tibble)
library(dplyr)

# mixedsort
library(gtools)

# Data directory CytoExploreRData ----------------------------------------------
datadir <- system.file("extdata", package = "CytoExploreRData")

# Activation GatingSet ---------------------------------------------------------

# Assign Activation dataset to fs -
files <- mixedsort(list.files(path = paste0(datadir, "/Activation"), 
                              pattern = ".fcs", 
                              full.names = TRUE))
fs <- read.flowSet(files = files)

# Sample for speed -
fs <- cyto_sample(fs, display = 2000, seed = 20)

# Samplenames
nms <- sampleNames(fs)

# Extract channels
chans <- colnames(fs)

# pData information -
pData(fs)$OVAConc <- c(rep(c(0, 0, 0.005, 0.005, 0.05, 0.05, 0.5, 0.5), 4), 0)
pData(fs)$Treatment <- c(
  rep("Stim-A", 8),
  rep("Stim-B", 8),
  rep("Stim-C", 8),
  rep("Stim-D", 8),
  "NA"
)
pData(fs)$Treatment <- factor(pData(fs)$Treatment, levels = c(
  "Stim-A",
  "Stim-B",
  "Stim-C",
  "Stim-D", "NA"
))
chnls <- c(
  "Alexa Fluor 405-A",
  "Alexa Fluor 430-A",
  "APC-Cy7-A", "PE-A",
  "Alexa Fluor 488-A",
  "Alexa Fluor 700-A",
  "Alexa Fluor 647-A",
  "7-AAD-A"
)
markers <- c(
  "Hoechst-405",
  "Hoechst-430",
  "CD11c",
  "Va2",
  "CD8",
  "CD4",
  "CD44",
  "CD69"
)
names(markers) <- chnls
markernames(fs) <- markers

# GatingSet -
gs <- GatingSet(fs)

# Compensation -
gs <- compensate(gs, fs[[1]]@description$SPILL)

# Transformation -
trans <- estimateLogicle(gs[[32]], cyto_fluor_channels(gs))
gs <- transform(gs, trans)

# gatingTemplate -
gt <- gatingTemplate(paste0(datadir, "/Activation-gatingTemplate.csv"))

# gatingTemplate file -
gtf <- read.csv(system.file("extdata",
  "Activation-gatingTemplate.csv",
  package = "CytoExploreRData"
))

# Gating -
suppressWarnings(gt_gating(gt, gs))

# Extract T Cells Population -
Va2 <- gs_pop_get_data(gs, "T Cells")

# Extract a flowFrame/flowSet for testing -
fr_test <- Va2[[32]]
fs_test <- Va2[c(26, 28, 30, 32)]
gs_test <- gs[c(26, 28, 30, 32)]

# Compensation GatingSet -------------------------------------------------------
Compensation <- read.flowSet(
  path = paste0(
    datadir,
    "/Compensation"
  ),
  pattern = ".fcs"
)

# Sample for speed -
Comp <- fsApply(Compensation, function(fr) {
  cyto_sample(fr, display = 500)
})

# GatingSet -
gsc <- GatingSet(Comp)

# Transformed GatingSet -
gsct <- cyto_transform(gs_clone(gsc),
                       trans_type = "biex",
                       plot = FALSE)

# gatingTemplate -
gtc <- gatingTemplate(paste0(datadir, "/Compensation-gatingTemplate.csv"))

# gatingTemplate file
gtcf <- read.csv(system.file("extdata",
  "Compensation-gatingTemplate.csv",
  package = "CytoExploreRData"
))

# Gating -
gt_gating(gtc, gsc)
gt_gating(gtc, gsct)

# Path to channel_match file
channel_match_file <- paste0(datadir, "/Compensation-Channels.csv")

# GATE OBJECTS -----------------------------------------------------------------

# 1D RECTANGLEGATE
coords <- matrix(c(0, 50000), ncol = 1)
colnames(coords) <- "FSC-A"
rownames(coords) <- c("min","max")
rg1 <- rectangleGate(.gate = coords, filterId = "Cells")

# 2D RECTANGLEGATE
coords <- matrix(c(0, 50000, 0, 50000), ncol = 2, byrow = FALSE)
colnames(coords) <- c("FSC-A", "SSC-A")
rownames(coords) <- c("min","max")
rg2 <- rectangleGate(.gate = coords, filterId = "Cells")

# POLYGONGATE
coords <- matrix(c(0, 50000, 50000, 0, 0, 0, 50000, 50000), 
                 ncol = 2, byrow = FALSE)
colnames(coords) <- c("FSC-A", "SSC-A")
pg <- polygonGate(.gate = coords, filterId = "Cells")

# ELLIPSOIDGATE
cov <- matrix(c(900000000, 0.00000008304362, 0.00000008304362, 2256250000),
              ncol = 2,
              dimnames = list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A"))
)
mean <- c("FSC-A" = 65000, "SSC-A" = 51250)
eg <- ellipsoidGate(filterId = "Cells", .gate = cov, mean = mean)

# QUADGATE
coords <- matrix(c(50000, 50000), ncol = 2)
colnames(coords) <- c("FSC-A","SSC-A")
qg <- quadGate(.gate = coords, filterId = "A|B|C|D")
