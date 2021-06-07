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

gs <- cyto_copy(gs)
gs <- cyto_barcode(gs, "events")

cs <- cyto_data_extract(gs, 
                        parent = "root",
                        copy = TRUE)[["root"]]

gs_sub <- cyto_sample(gs, 
                      display = 2000,
                      seed = 56)

cs_sub <- cyto_data_extract(gs_sub, 
                            parent = "root",
                            copy = TRUE)[["root"]]

# Compensation GatingSet -------------------------------------------------------

gs_comp <- cyto_load(
  system.file("extdata/Compensation-GatingSet",
              package = "CytoExploreRData")
)

gs_comp <- cyto_copy(gs_comp)
gs_comp <- cyto_barcode(gs_comp, "events")

cs_comp <- cyto_data_extract(gs_comp, 
                             parent = "root",
                             copy = TRUE)[["root"]]
