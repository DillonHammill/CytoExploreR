## CYTO_GATE_EXTRACT -----------------------------------------------------------

test_that("cyto_gate_extract", {
  # GH/GS - BOOLEANFILTER
  gates <- cyto_gate_extract(gs,
                             alias = c("Dead Cells", "Live Cells"),
                             merge_by = "Treatment",
                             select = 1:32)
  expect_setequal(names(gates), c("Stim-A", "Stim-B", "Stim-C", "Stim-D"))
  expect_setequal(names(gates[[1]]), c("Live Cells", "Dead Cells"))
  expect_setequal(LAPPLY(unlist(gates[[1]]), "cyto_class"),
                  c("rectangleGate", "complementFilter"))
  # GATINGTEMPLATE
  gates <- cyto_gate_extract(gatingTemplate =  Activation_gatingTemplate,
                             parent = "Single Cells",
                             alias = c("Dead Cells", "Live Cells"))
  expect_setequal(names(gates), c("Live Cells", "Dead Cells"))
  expect_setequal(LAPPLY(gates, "cyto_class"),
                  c("rectangleGate", "complementFilter"))
})