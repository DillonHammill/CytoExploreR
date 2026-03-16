library(CytoExploreRData)

# Load datasets explicitly into the test environment
data("Activation",                 package = "CytoExploreRData", envir = environment())
data("Activation_gatingTemplate",  package = "CytoExploreRData", envir = environment())

# Activation GatingSet ---------------------------------------------------------

gs <- GatingSet(Activation)
gs <- cyto_transform(gs)
gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
gs <- cyto_barcode(gs, "events")

# Root-level cytoset for testing cyto_apply.flowSet directly
cs <- cyto_data_extract(gs, parent = "root", copy = TRUE)[["root"]]
