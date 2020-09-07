# Required Packages ------------------------------------------------------------
library(CytoExploreRData)

# Non-interactive mode
options("CytoExploreR_interactive" = FALSE)

# Activation GatingSet ---------------------------------------------------------

gs <- cyto_load(
  system.file("extdata/Activation-GatingSet", 
              package = "CytoExploreRData")
  )

cs <- cyto_extract(gs, "root")

# Compensation GatingSet -------------------------------------------------------

gs_comp <- cyto_load(
  system.file("extdata/Compensation-GatingSet",
              package = "CytoExploreRData")
  )

cs_comp <- cyto_extract(gs_comp, "root")