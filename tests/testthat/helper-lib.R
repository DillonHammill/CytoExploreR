# Load Packages ----------------------------------------------------------------
library(CytoRSuiteData)

# Load in vdiffr for image comparison -
library(vdiffr)

# Load in grDevices for recordPlot() -
library(grDevices)

# Load in robustbase (colMedians) -
library(robustbase)

# Libraries for data manipulation
library(tibble)
library(dplyr)

# Data directory CytoRSuiteData ------------------------------------------------
datadir <- system.file("extdata", package = "CytoRSuiteData")

# Activation GatingSet ---------------------------------------------------------

# Assign Activation dataset to fs -
fs <- read.flowSet(path = paste0(datadir, "/Activation"), pattern = ".fcs")

# Fix file order -
fs <- fs[c(
  1, 12, 23, 28, 29, 30, 31, 32, 33, 2, 3, 4, 5, 6, 7, 8, 9, 10,
  11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27
)]

# Sample for speed -
fs <- cyto_sample(fs, display = 0.04, seed = 20)

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
  package = "CytoRSuiteData"
))

# Gating -
gating(gt, gs)

# Extract T Cells Population -
Va2 <- getData(gs, "T Cells")

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
  cyto_sample(fr, display = 0.04, seed = 20)
})

# GatingSet -
gsc <- GatingSet(Comp)

# gatingTemplate -
gtc <- gatingTemplate(paste0(datadir, "/Compensation-gatingTemplate.csv"))

# gatingTemplate file
gtcf <- read.csv(system.file("extdata",
  "Compensation-gatingTemplate.csv",
  package = "CytoRSuiteData"
))

# Gating -
gating(gtc, gsc)

# Spillover Matrix -------------------------------------------------------------

spill <- read.csv("Ref-Spillover-Matrix.csv", header = TRUE, row.names = 1)
colnames(spill) <- rownames(spill)

# Construct gate objects -------------------------------------------------------

# rectangleGate -
coords <- matrix(c(25000, 150000, 5000, 150000), ncol = 2, nrow = 2)
colnames(coords) <- c("FSC-A", "SSC-A")
rownames(coords) <- c("min", "max")
rg <- rectangleGate(filterId = "Cells", .gate = coords)

# polygonGate -
coords <- matrix(c(
  50000, 100000, 100000, 75000, 50000,
  10000, 10000, 60000, 85000, 60000
),
ncol = 2, nrow = 5
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg <- polygonGate(filterId = "Cells", .gate = coords)

# intervalGate -
# 1D x axis -
coords <- matrix(c(25000, 150000), ncol = 1, nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", "max")
igx <- rectangleGate(filterId = "Cells", .gate = coords)

# 2D x axis-
coords <- matrix(c(25000, 150000, -Inf, Inf), ncol = 2, nrow = 2)
colnames(coords) <- c("FSC-A", "SSC-A")
rownames(coords) <- c("min", "max")
ig <- rectangleGate(filterId = "Cells", .gate = coords)

# 2D y axis -
coords <- matrix(c(-Inf, Inf, 25000, 150000), ncol = 2, nrow = 2)
colnames(coords) <- c("FSC-A", "SSC-A")
rownames(coords) <- c("min", "max")
igy <- rectangleGate(filterId = "Cells", .gate = coords)

# thresholdGate -
# 2D -
coords <- matrix(c(25000, Inf, 5000, Inf), ncol = 2, nrow = 2)
colnames(coords) <- c("FSC-A", "SSC-A")
rownames(coords) <- c("min", "max")
tg <- rectangleGate(filterId = "Cells", .gate = coords)

# 1D -
coords <- matrix(c(25000, Inf), ncol = 1, nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", "max")
tg1 <- rectangleGate(filterId = "Cells", .gate = coords)

# boundaryGate -
# 2D -
coords <- matrix(c(-Inf, 200000, -Inf, 200000), ncol = 2, nrow = 2)
colnames(coords) <- c("FSC-A", "SSC-A")
rownames(coords) <- c("min", "max")
bg <- rectangleGate(filterId = "Cells", .gate = coords)

# 1D -
coords <- matrix(c(-Inf, 200000), ncol = 1, nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", "max")
bg1 <- rectangleGate(filterId = "Cells", .gate = coords)

# ellipsoidGate -
cov <- matrix(c(900000000, 0.00000008304362, 0.00000008304362, 2256250000),
  ncol = 2,
  dimnames = list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A"))
)
mean <- c("FSC-A" = 65000, "SSC-A" = 51250)
eg <- ellipsoidGate(filterId = "Cells", .gate = cov, mean = mean)

# drawQuadrants -
coords <- matrix(c(-Inf, 150000, -Inf, 150000),
  ncol = 2,
  dimnames = list(c("min", "max"), c("FSC-A", "SSC-A"))
)
rg1 <- rectangleGate(filterId = "A", .gate = coords)
coords <- matrix(c(150000, Inf, -Inf, 150000),
  ncol = 2,
  dimnames = list(c("min", "max"), c("FSC-A", "SSC-A"))
)
rg2 <- rectangleGate(filterId = "B", .gate = coords)
coords <- matrix(c(150000, Inf, 150000, Inf),
  ncol = 2,
  dimnames = list(c("min", "max"), c("FSC-A", "SSC-A"))
)
rg3 <- rectangleGate(filterId = "C", .gate = coords)
coords <- matrix(c(-Inf, 150000, 150000, Inf),
  ncol = 2,
  dimnames = list(c("min", "max"), c("FSC-A", "SSC-A"))
)
rg4 <- rectangleGate(filterId = "D", .gate = coords)
qg <- filters(list(rg1, rg2, rg3, rg4))

# drawWeb -
coords <- matrix(c(
  120627.90, 4610.20, 4610.20, 63013.78,
  147367.85, 104838.23, 729.81, 729.81
),
ncol = 2,
nrow = 4
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg1 <- polygonGate(filterId = "A", .gate = coords)
coords <- matrix(c(120627.90, 63013.78, 168882.41, 147367.85, 729.81, 729.81),
  ncol = 2,
  nrow = 3
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg2 <- polygonGate(filterId = "B", .gate = coords)
coords <- matrix(c(120627.90, 168882.41, 248647.95, 147367.85, 729.81, 729.81),
  ncol = 2,
  nrow = 3
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg3 <- polygonGate(filterId = "C", .gate = coords)
coords <- matrix(c(
  120627.90, 248647.95, 262143.00, 262143.00, 147367.85,
  729.81, 729.81, 184232.42
),
ncol = 2,
nrow = 4
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg4 <- polygonGate(filterId = "D", .gate = coords)
coords <- matrix(c(
  120627.9, 262143.0, 262143.0, 233708.6,
  147367.9, 184232.4, 262143.0, 262143.0
),
ncol = 2,
nrow = 4
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg5 <- polygonGate(filterId = "E", .gate = coords)
coords <- matrix(c(
  120627.9, 233708.6, 107041.4, 147367.9,
  262143.0, 262143.0
),
ncol = 2,
nrow = 3
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg6 <- polygonGate(filterId = "F", .gate = coords)
coords <- matrix(c(
  120627.9, 107041.4, 4610.2, 4610.2, 147367.9, 262143.0,
  262143.0, 237418.4
),
ncol = 2,
nrow = 4
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg7 <- polygonGate(filterId = "G", .gate = coords)
coords <- matrix(c(
  120627.9, 4610.2, 4610.2, 147367.9,
  237418.4, 104838.2
),
ncol = 2,
nrow = 3
)
colnames(coords) <- c("FSC-A", "SSC-A")
pg8 <- polygonGate(filterId = "H", .gate = coords)

wg <- filters(list(pg1, pg2, pg3, pg4, pg5, pg6, pg7, pg8))
