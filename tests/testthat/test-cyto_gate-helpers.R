context("cyto_gate-helpers")

# Turn off graphics device - missing coords are inherited from par("usr)
dev.off()

# CYTO_GATE_CONVERT ------------------------------------------------------------

test_that("cyto_gate_convert", {
  
  # Error
  expect_error(cyto_gate_convert(rg),
               "Supply the channels in which the gate will be plotted")
  expect_error(cyto_gate_convert(rg, channels = c("PE-A","Alexa Fluor 488-A")),
               "Cannot plot gates not constructed in the supplied channels.")
  expect_error(cyto_gate_convert(pg),
               "Supply the channels in which the gate will be plotted")
  expect_error(cyto_gate_convert(pg, channels = c("PE-A","Alexa Fluor 488-A")),
               "Cannot plot gates not constructed in the supplied channels.")
  expect_error(cyto_gate_convert(eg),
               "Supply the channels in which the gate will be plotted")
  expect_error(cyto_gate_convert(eg, channels = c("PE-A","Alexa Fluor 488-A")),
               "Cannot plot gates not constructed in the supplied channels.")
  
  
  
  # rectangleGate --------------------------------------------------------------
  
  # 1D -> 1D
  expect_equal(cyto_gate_convert(rg["FSC-A"], channels = "FSC-A"), rg["FSC-A"])
  
  # 2D -> 2D
  expect_equal(cyto_gate_convert(rg, channels = c("FSC-A","SSC-A")), rg)
  
  # 2D  -> 2D X Channel Match
  coords <- matrix(c(0,
                     1,
                     as.numeric(rg@min["SSC-A"]),
                     as.numeric(rg@max["SSC-A"])), 
                   ncol = 2, 
                   nrow = 2, 
                   byrow = FALSE)
  colnames(coords) <- c("PE-A","SSC-A")
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(rg, channels = c("PE-A","SSC-A")),
               rg1)
  
  # 2D  -> 2D Y Channel Match
  coords <- matrix(c(as.numeric(rg@min["SSC-A"]),
                     as.numeric(rg@max["SSC-A"]),
                     0,
                     1), 
                   ncol = 2, 
                   nrow = 2, 
                   byrow = FALSE)
  colnames(coords) <- c("SSC-A","PE-A")
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(rg, channels = c("SSC-A","PE-A")),
               rg1)
  
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
  
  # 2D -> 2D X Channel Match
  coords <- matrix(c(50000, 100000, 0, 1), ncol = 2, nrow = 2, byrow = FALSE)
  colnames(coords) <- c("FSC-A","PE-A")
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(pg, channels = c("FSC-A","PE-A")), rg1)
  
  # 2D -> 2D Y Channel Match
  coords <- matrix(c(0 ,1, 50000, 100000), ncol = 2, nrow = 2, byrow = FALSE)
  colnames(coords) <- c("PE-A","FSC-A")
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(pg, channels = c("PE-A","FSC-A")), rg1)
  
  # 2D -> 1D
  coords <- matrix(c(50000, 100000), ncol = 1, nrow = 2)
  colnames(coords) <- "FSC-A"
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(pg, channels = "FSC-A"), rg1)
  
  # ellipsoidGate --------------------------------------------------------------
  
  # 2D -> 2D
  expect_equal(cyto_gate_convert(eg, channels = c("FSC-A","SSC-A")), eg)
  
  # 2D -> 2D X Channel Match
  coords <- matrix(c(35061.64, 95000, 0 , 1), ncol = 2, nrow = 2, byrow = FALSE)
  colnames(coords) <- c("FSC-A","PE-A")
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(eg, channels = c("FSC-A","PE-A")), 
               rg1, tolerance = 0.01)
  
  # 2D -> 2D Y Channel Match
  coords <- matrix(c(0, 1, 35061.64, 95000), ncol = 2, nrow = 2, byrow = FALSE)
  colnames(coords) <- c("PE-A","FSC-A")
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(eg, channels = c("PE-A","FSC-A")), 
               rg1, tolerance = 0.01)
  
  # 2d -> 1D
  coords <- matrix(c(35061.64, 95000), ncol = 1, nrow = 2)
  colnames(coords) <- "FSC-A"
  rownames(coords) <- c("min","max")
  rg1 <- rectangleGate(filterId = "Cells", .gate = coords)
  expect_equal(cyto_gate_convert(eg, channels = "FSC-A"), rg1, tolerance = 0.01)
  
  # list and filters -----------------------------------------------------------
  
  # Complex gate objects -> list of gate objects
  gts <- filters(list(rg, pg, eg))
  gts <- list(gts, rg, pg, eg)
  expect_equal(cyto_gate_convert(gts, channels = c("FSC-A","SSC-A")), 
               list(rg,pg,eg,rg,pg,eg), tolerance = 0.01)
  
})