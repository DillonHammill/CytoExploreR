context("cyto_gate_draw-helpers")

# .CYTO_GATE_TYPE --------------------------------------------------------------

test_that(".cyto_gate_type", {
  # ABBREVIATIONS
  expect_equal(.cyto_gate_type(type = c("r","e","p","b","i","t"),
                               channels = c("FSC-A","SSC-A"),
                               alias = c("A","B","C","D","E","F"),
                               negate = FALSE),
               c("rectangle",
                 "ellipse",
                 "polygon",
                 "boundary",
                 "interval",
                 "threshold"))
  # DEFAULT TYPES
  expect_equal(.cyto_gate_type(channels = "FSC-A",
                               alias = c("A","B","C"),
                               negate = FALSE),
               rep("interval", 3))
  expect_equal(.cyto_gate_type(channels = c("FSC-A","SSC-A"),
                               alias = c("A","B"),
                               negate = FALSE),
               rep("polygon", 2))
  # NEGATE GATE TYPE LENGTH
  expect_equal(.cyto_gate_type(channels = c("FSC-A","SSC-A"),
                               alias = c("A","B", "C"),
                               negate = TRUE),
               c(rep("polygon", 2), NA))
  # UNSUPPORTED TYPES
  expect_error(.cyto_gate_type(type = "test",
                               channels = c("FSC-A","SSC-A"),
                               alias = "A"),
               "test is not a valid gate type for cyto_gate_draw.")
  expect_error(.cyto_gate_type(type = c("testA","testB"),
                               channels = c("FSC-A","SSC-A"),
                               alias = c("A","B")),
               "testA & testB are not valid gate types for cyto_gate_draw.")
  # MIXED GATE TYPES NOT SUPPORTED - SPLIT METHODS
  expect_error(.cyto_gate_type(type = c("w","r"),
                               channels = c("FSC-A","SSC-A"),
                               alias = c("A","B","C"),
                               negate = FALSE),
               "Mixed gates are not supported for web gate types.")
  # NEGATE - SPLIT METHODS EXCLUDED
  expect_error(.cyto_gate_type(type = "q",
                               channels = c("FSC-A","SSC-A"),
                               alias = c("A","B","C","D"),
                               negate = TRUE),
               "Cannot negate quadrant gate types.")
  # SUPPORTED 1D GATE TYPES
  expect_error(.cyto_gate_type(type = c("p","r","e"),
                               channels = "FSC-A",
                               alias = c("A","B","C")),
          "Supported 1D gate types include interval, boundary and threshold.")
})

# .CYTO_ALIAS ------------------------------------------------------------------

test_that(".cyto_alias", {
  # SPLIT ALIAS
  expect_equal(.cyto_alias(alias = c("A","B","C"),
                           type = c("polygon", "ellipse","rectangle")),
               list("A", "B", "C"))
  expect_equal(.cyto_alias(alias = c("A","B","C","D"),
                           type = c("quadrant")),
               list(c("A","B","C","D")))
  # NEGATE
  expect_equal(.cyto_alias(alias = c("A","B","C","D"),
                           type = c("polygon","polygon","polygon", NA)),
               list("A","B","C","D"))
  # MISSING NEGATE LABEL
  expect_error(.cyto_alias(alias = c("A","B","C"),
                           type = c("polygon","polygon","polygon", NA)),
               "Supply a name for each of the population(s) to 'alias'.",
               fixed = TRUE)
})
