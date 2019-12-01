# Load Packages ----------------------------------------------------------------
library(CytoExploreRData)

# Load in mockery for mocking
library(mockery)

# Load in vdiffr for image comparison -
library(vdiffr)

# Load in grDevices for recordPlot() -
library(grDevices)

# Load in robustbase (colMedians) -
library(robustbase)

# Libraries for data manipulation
library(tibble)
library(dplyr)

# Data directory CytoExploreRData ----------------------------------------------
datadir <- system.file("extdata", package = "CytoExploreRData")

# Activation GatingSet ---------------------------------------------------------

# Assign Activation dataset to fs -
fs <- cyto_load(paste0(datadir, "/Activation"))

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
chans <- c(
  "Alexa Fluor 405-A",
  "Alexa Fluor 430-A",
  "APC-Cy7-A", 
  "PE-A",
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
names(markers) <- chans
markernames(fs) <- markers

# Sample for speed -
fs <- cyto_sample(fs, display = 2000, seed = 20)

# GatingSet -
gs <- GatingSet(fs)

# Compensation -
gs <- cyto_compensate(gs)

# Transformation -
trans <- cyto_transformer_logicle(gs)
gs <- transform(gs, trans, plot = FALSE)

# gatingTemplate -
gt <- gatingTemplate(paste0(datadir, "/Activation-gatingTemplate.csv"))

# gatingTemplate file -
gt_file <- read.csv(system.file("extdata",
  "Activation-gatingTemplate.csv",
  package = "CytoExploreRData"
))

# Gating -
suppressWarnings(gt_gating(gt, gs))

# Compensation GatingSet -------------------------------------------------------
fs_comp <- cyto_load(paste0(datadir,"/Compensation"))

# Sample for speed -
fs_comp <- cyto_sample(fs_comp, display = 2000, seed = 20)

# GatingSet -
gs_comp <- GatingSet(fs_comp)

# Transformed GatingSet -
gs_comp_trans <- cyto_transform(gs_clone(gs_comp),
                                type = "logicle",
                                plot = FALSE)

# gatingTemplate -
gt_comp <- gatingTemplate(paste0(datadir, "/Compensation-gatingTemplate.csv"))

# gatingTemplate file
gt_comp_file <- read.csv(system.file("extdata",
  "Compensation-gatingTemplate.csv",
  package = "CytoExploreRData"
))

# Gating -
gt_gating(gt_comp, gs_comp)

# Path to channel_match file
channel_match_file <- paste0(datadir, "/Compensation-Channels.csv")

# MODIFIED ACTIVATION GATE OBJECTS ---------------------------------------------

# CD4 T Cells
coords <- matrix(c(1.824,
                   3.190, 
                   -0.227, 
                   2.111), 
                 ncol = 2, 
                 nrow = 2)
colnames(coords) <- c("Alexa Fluor 700-A","Alexa Fluor 488-A")
rownames(coords) <- c("min","max")
CD4 <- rectangleGate(filterId = "CD4 T Cells", .gate = coords)

# CD8 T Cells
coords <- matrix(c(-0.3134,
                   1.8687, 
                   2.2221,
                   4.1452), 
                 ncol = 2, 
                 nrow = 2)
colnames(coords) <- c("Alexa Fluor 700-A","Alexa Fluor 488-A")
rownames(coords) <- c("min","max")
CD8 <- rectangleGate(filterId = "CD8 T Cells", .gate = coords)

# GATE OBJECTS -----------------------------------------------------------------

# 1D RECTANGLEGATE (INTERVAL)
coords <- matrix(c(0, 
                   50000), 
                 ncol = 1)
colnames(coords) <- "FSC-A"
rownames(coords) <- c("min",
                      "max")
rg1 <- rectangleGate(.gate = coords, 
                     filterId = "Cells")

# 2D RECTANGLEGATE (RECTANGLE)
coords <- matrix(c(0, 
                   50000,
                   0, 
                   50000), 
                 ncol = 2, 
                 byrow = FALSE)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
rownames(coords) <- c("min",
                      "max")
rg2 <- rectangleGate(.gate = coords,
                     filterId = "Cells")

# POLYGONGATE
coords <- matrix(c(0,
                   50000,
                   50000,
                   0, 
                   0, 
                   0, 
                   50000, 
                   50000), 
                 ncol = 2, 
                 byrow = FALSE)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
pg <- polygonGate(.gate = coords, 
                  filterId = "Cells")

# ELLIPSOIDGATE
cov <- matrix(c(900000000,
                0.00000008304362,
                0.00000008304362,
                2256250000),
              ncol = 2,
              dimnames = list(c("FSC-A",
                                "SSC-A"), 
                              c("FSC-A", 
                                "SSC-A"))
)
mean <- c("FSC-A" = 65000, 
          "SSC-A" = 51250)
eg <- ellipsoidGate(filterId = "Cells", 
                    .gate = cov, mean = mean)

# QUADGATE
coords <- matrix(c(50000,
                   50000), 
                 ncol = 2)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
qg <- quadGate(.gate = coords, 
               filterId = "A|B|C|D")

# GATE OBJECT TYPES ------------------------------------------------------------

# intervalGate -
# 1D x axis -
coords <- matrix(c(25000, 
                   150000), 
                 ncol = 1, 
                 nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", 
                      "max")
igx <- rectangleGate(filterId = "Cells", 
                     .gate = coords)

# 2D x axis-
coords <- matrix(c(25000,
                   150000, 
                   -Inf, 
                   Inf), 
                 ncol = 2, 
                 nrow = 2)
colnames(coords) <- c("FSC-A", 
                      "SSC-A")
rownames(coords) <- c("min",
                      "max")
ig <- rectangleGate(filterId = "Cells", 
                    .gate = coords)

# 2D y axis -
coords <- matrix(c(-Inf,
                   Inf, 
                   25000, 
                   150000),
                 ncol = 2,
                 nrow = 2)
colnames(coords) <- c("FSC-A", 
                      "SSC-A")
rownames(coords) <- c("min", 
                      "max")
igy <- rectangleGate(filterId = "Cells",
                     .gate = coords)

# thresholdGate -
# 2D -
coords <- matrix(c(25000,
                   Inf, 
                   5000, 
                   Inf), 
                 ncol = 2, 
                 nrow = 2)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
rownames(coords) <- c("min", 
                      "max")
tg <- rectangleGate(filterId = "Cells",
                    .gate = coords)

# 1D -
coords <- matrix(c(25000, 
                   Inf), 
                 ncol = 1, 
                 nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", 
                      "max")
tg1 <- rectangleGate(filterId = "Cells",
                     .gate = coords)

# boundaryGate -
# 2D -
coords <- matrix(c(-Inf,
                   200000, 
                   -Inf, 
                   200000), 
                 ncol = 2, 
                 nrow = 2)
colnames(coords) <- c("FSC-A", 
                      "SSC-A")
rownames(coords) <- c("min", 
                      "max")
bg <- rectangleGate(filterId = "Cells",
                    .gate = coords)

# 1D -
coords <- matrix(c(-Inf, 
                   200000), 
                 ncol = 1, 
                 nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", 
                      "max")
bg1 <- rectangleGate(filterId = "Cells", 
                     .gate = coords)

# quadrantGates -
coords <- matrix(c(-Inf, 
                   150000, -Inf, 
                   150000), 
                 ncol = 2, 
                 dimnames = list(c("min", 
                                   "max"), 
                                 c("FSC-A", 
                                   "SSC-A")))
q1 <- rectangleGate(filterId = "A",
                     .gate = coords)
coords <- matrix(c(150000, 
                   Inf, 
                   -Inf, 
                   150000), 
                 ncol = 2, 
                 dimnames = list(c("min",
                                   "max"), 
                                 c("FSC-A",
                                   "SSC-A")))
q2 <- rectangleGate(filterId = "B", 
                     .gate = coords)
coords <- matrix(c(150000,
                   Inf, 
                   150000, 
                   Inf),
                 ncol = 2, 
                 dimnames = list(c("min",
                                   "max"), 
                                 c("FSC-A", 
                                   "SSC-A")))
q3 <- rectangleGate(filterId = "C",
                     .gate = coords)
coords <- matrix(c(-Inf, 
                   150000, 
                   150000, 
                   Inf), 
                 ncol = 2, 
                 dimnames = list(c("min",
                                   "max"), 
                                 c("FSC-A",
                                   "SSC-A")))
q4 <- rectangleGate(filterId = "D",
                     .gate = coords)
qg2 <- filters(list(q1,
                    q2, 
                    q3, 
                    q4))

# web (EXPERIMENTAL) -
coords <- matrix(c(120627.90, 
                   4610.20, 
                   4610.20, 
                   63013.78, 
                   147367.85,
                   104838.23, 
                   729.81,
                   729.81), 
                 ncol = 2, 
                 nrow = 4)
colnames(coords) <- c("FSC-A", "SSC-A")
pg1 <- polygonGate(filterId = "A", 
                   .gate = coords)
coords <- matrix(c(120627.90, 
                   63013.78, 
                   168882.41, 
                   147367.85, 
                   729.81, 
                   729.81), 
                 ncol = 2, 
                 nrow = 3)
colnames(coords) <- c("FSC-A","SSC-A")
pg2 <- polygonGate(filterId = "B", 
                   .gate = coords)
coords <- matrix(c(120627.90,
                   168882.41,
                   248647.95, 
                   147367.85, 
                   729.81,
                   729.81), 
                 ncol = 2, 
                 nrow = 3)
colnames(coords) <- c("FSC-A", "SSC-A")
pg3 <- polygonGate(filterId = "C", 
                   .gate = coords)
coords <- matrix(c(120627.90,
                   248647.95,
                   262143.00,
                   262143.00,
                   147367.85,
                   729.81,
                   729.81,
                   184232.42), 
                 ncol = 2, 
                 nrow = 4)
colnames(coords) <- c("FSC-A", "SSC-A")
pg4 <- polygonGate(filterId = "D", 
                   .gate = coords)
coords <- matrix(c(120627.9,
                   262143.0,
                   262143.0,
                   233708.6,
                   147367.9,
                   184232.4,
                   262143.0,
                   262143.0), 
                 ncol = 2, 
                 nrow = 4)
colnames(coords) <- c("FSC-A", "SSC-A")
pg5 <- polygonGate(filterId = "E", 
                   .gate = coords)
coords <- matrix(c(120627.9,
                   233708.6,
                   107041.4,
                   147367.9, 
                   262143.0,
                   262143.0), 
                 ncol = 2, 
                 nrow = 3)
colnames(coords) <- c("FSC-A", "SSC-A")
pg6 <- polygonGate(filterId = "F", 
                   .gate = coords)
coords <- matrix(c(120627.9, 
                   107041.4, 
                   4610.2, 
                   4610.2, 
                   147367.9, 
                   262143.0,
                   262143.0, 
                   237418.4), 
                 ncol = 2, 
                 nrow = 4)
colnames(coords) <- c("FSC-A", "SSC-A")
pg7 <- polygonGate(filterId = "G", .gate = coords)
coords <- matrix(c(120627.9, 
                   4610.2, 
                   4610.2, 
                   147367.9,
                   237418.4,
                   104838.2), 
                 ncol = 2, 
                 nrow = 3)
colnames(coords) <- c("FSC-A", "SSC-A")
pg8 <- polygonGate(filterId = "H",
                   .gate = coords)

wg <- filters(list(pg1, pg2, pg3, pg4, pg5, pg6, pg7, pg8))
