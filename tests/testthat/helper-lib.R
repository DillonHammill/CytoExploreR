# Load Packages ----------------------------------------------------------------
library(CytoExploreRData)

# Load in mockery for mocking
library(mockery)

# Load in vdiffr for image comparison -
# library(vdiffr)

# Load in grDevices for recordPlot() -
library(grDevices)

# Load in robustbase (colMedians) -
library(robustbase)

# Libraries for data manipulation
library(tibble)
library(dplyr)

# Turn off interactive mode
options("CytoExploreR_interactive" = FALSE)

# Data directory CytoExploreRData ----------------------------------------------
datadir <- system.file("extdata", package = "CytoExploreRData")

# Activation GatingSet ---------------------------------------------------------

# Assign Activation dataset to fs -
fs <- cyto_load(paste0(datadir, "/Activation"))

# pData information -
pData(fs)$OVAConc <- c(rep(c(0, 0, 5, 5, 50, 50, 500, 500), 4), 0)
class(pData(fs)$OVAConc) <- "numeric"
pData(fs)$Treatment <- c(
  rep("Stim-A", 8),
  rep("Stim-B", 8),
  rep("Stim-C", 8),
  rep("Stim-D", 8),"NA"
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
trans <- cyto_transformer_logicle(gs, plot = FALSE)
gs <- transform(gs, trans)

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
gt_gating(gt_comp, gs_comp_trans)

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
                     filterId = "P")

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
                     filterId = "A")

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
                  filterId = "F")

# ELLIPSOIDGATE
cov <- matrix(c(625249999.99999988079,
                0.00000000311,
                0.00000000311,
                676000000.00000000000),
              ncol = 2,
              dimnames = list(c("FSC-A",
                                "SSC-A"), 
                              c("FSC-A", 
                                "SSC-A"))
)
mean <- c("FSC-A" = 50000, 
          "SSC-A" = 50250)
eg <- ellipsoidGate(filterId = "G", 
                    .gate = cov, mean = mean)

# QUADGATE
coords <- matrix(c(50000,
                   50000), 
                 ncol = 2)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
qg <- quadGate(.gate = coords, 
               filterId = "H|I|J|K")

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
ig <- rectangleGate(filterId = "B", 
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
igy <- rectangleGate(filterId = "C",
                     .gate = coords)

# thresholdGate -
# 2D -
coords <- matrix(c(25000,
                   Inf, 
                   50000, 
                   Inf), 
                 ncol = 2, 
                 nrow = 2)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
rownames(coords) <- c("min", 
                      "max")
tg <- rectangleGate(filterId = "D",
                    .gate = coords)

# 1D -
coords <- matrix(c(25000, 
                   Inf), 
                 ncol = 1, 
                 nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", 
                      "max")
tg1 <- rectangleGate(filterId = "Q",
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
bg <- rectangleGate(filterId = "E",
                    .gate = coords)

# 1D -
coords <- matrix(c(-Inf, 
                   200000), 
                 ncol = 1, 
                 nrow = 2)
colnames(coords) <- c("FSC-A")
rownames(coords) <- c("min", 
                      "max")
bg1 <- rectangleGate(filterId = "R", 
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
q1 <- rectangleGate(filterId = "A", .gate = coords)
coords <- matrix(c(150000, 
                   Inf, 
                   -Inf, 
                   150000), 
                 ncol = 2, 
                 dimnames = list(c("min",
                                   "max"), 
                                 c("FSC-A",
                                   "SSC-A")))
q2 <- rectangleGate(filterId = "B", .gate = coords)
coords <- matrix(c(150000,
                   Inf, 
                   150000, 
                   Inf),
                 ncol = 2, 
                 dimnames = list(c("min",
                                   "max"), 
                                 c("FSC-A", 
                                   "SSC-A")))
q3 <- rectangleGate(filterId = "C", .gate = coords)
coords <- matrix(c(-Inf, 
                   150000, 
                   150000, 
                   Inf), 
                 ncol = 2, 
                 dimnames = list(c("min",
                                   "max"), 
                                 c("FSC-A",
                                   "SSC-A")))
q4 <- rectangleGate(filterId = "D", .gate = coords)
qg2 <- filters(list(q1,
                    q2, 
                    q3, 
                    q4))

# Web gates
coords <- matrix(c(100000.00,
                   -21810.30,
                   -21810.30,
                   160905.10,
                   100000.00,
                   18793.13,
                   -21810.30,
                   21810.30),
                 ncol = 2, 
                 nrow = 4,
                 byrow = FALSE)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
wg1 <- polygonGate(.gate = coords, 
                    filterId = "L")

coords <- matrix(c(100000.00,
                   160905.10,
                   283953.30,
                   283953.30,
                   100000.00,
                   -21810.30,
                   -21810.30,
                   165697.60),
                 ncol = 2, 
                 nrow = 4,
                 byrow = FALSE)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
wg2 <- polygonGate(.gate = coords, 
                   filterId = "M")

coords <- matrix(c(100000.00,
                   283953.30,
                   283953.30,
                   38682.23,
                   100000.00,
                   165697.60,
                   283953.30,
                   283953.30),
                 ncol = 2, 
                 nrow = 4,
                 byrow = FALSE)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
wg3 <- polygonGate(.gate = coords, 
                   filterId = "N")

coords <- matrix(c(100000.00,
                   38682.23,
                   -21810.30,
                   -21810.3,
                   100000.00,
                   283953.30,
                   283953.30,
                   18793.13),
                 ncol = 2, 
                 nrow = 4,
                 byrow = FALSE)
colnames(coords) <- c("FSC-A",
                      "SSC-A")
wg4 <- polygonGate(.gate = coords, 
                   filterId = "O")
wg <- filters(list(wg1,wg2,wg3,wg4))