# Required Packages ------------------------------------------------------------
library(CytoExploreRData)
library(flowWorkspace)
library(mockery)

# DIRECTORIES ------------------------------------------------------------------

temp_dir <- paste0(tempdir(), .Platform$file.sep)

# Activation GatingSet ---------------------------------------------------------

gs <- cyto_load(
  system.file("extdata/Activation-GatingSet", 
              package = "CytoExploreRData")
)

cs <- cyto_data_extract(gs, "root")[["root"]]

gs_sub <- cyto_sample(gs, 
                      display = 2000,
                      seed = 56)

cs_sub <- cyto_data_extract(gs_sub, "root")[["root"]]

# Compensation GatingSet -------------------------------------------------------

gs_comp <- cyto_load(
  system.file("extdata/Compensation-GatingSet",
              package = "CytoExploreRData")
)

cs_comp <- cyto_data_extract(gs_comp, "root")[["root"]]
