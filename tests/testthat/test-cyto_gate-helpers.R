context("cyto_gate-helpers")

# Turn off graphics device - missing coords are inherited from par("usr)
dev.off()

# CYTO_GATE_CONVERT ------------------------------------------------------------

test_that("cyto_gate_convert", {
  # Error
  expect_error(cyto_gate_convert(rg),
               "Supply the channels in which the gate will be plotted")
  expect_error(cyto_gate_convert(pg),
               "Supply the channels in which the gate will be plotted")
  expect_error(cyto_gate_convert(eg),
               "Supply the channels in which the gate will be plotted")
  
  # rectangleGate --------------------------------------------------------------
  
  # 1D -> 1D
  expect_equal(cyto_gate_convert(rg["FSC-A"], channels = "FSC-A"), rg["FSC-A"])
  
  # 2D -> 2D
  expect_equal(cyto_gate_convert(rg, channels = c("FSC-A","SSC-A")), rg)
  
  # 2D -> 1D
  expect_equal(cyto_gate_convert(rg, channels = "FSC-A"), rg["FSC-A"])
  
  # 1D -> 2D
  coords <- matrix(c(0,1), ncol = 1, nrow = 2)
  colnames(coords) <- "SSC-A"
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(rg["FSC-A"], channels = c("FSC-A","SSC-A")),
               rg["FSC-A"]*rg1)
  
  # polygonGate ----------------------------------------------------------------
  
  # 2D -> 2D
  expect_equal(cyto_gate_convert(pg, channels = c("FSC-A","SSC-A")), pg)
  
  # 2D -> 1D
  coords <- matrix(c(50000, 100000), ncol = 1, nrow = 2)
  colnames(coords) <- "FSC-A"
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(pg, channels = "FSC-A"), rg1)
  
  # ellipsoidGate --------------------------------------------------------------
  
  # 2D -> 2D
  expect_equal(cyto_gate_convert(eg, channels = c("FSC-A","SSC-A")), eg)
  
  # 2d -> 1D
  coords <- matrix(c(35061.64, 95000), ncol = 1, nrow = 2)
  colnames(coords) <- "FSC-A"
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(eg, channels = "FSC-A"), rg1, tolerance = 0.01)
  
})