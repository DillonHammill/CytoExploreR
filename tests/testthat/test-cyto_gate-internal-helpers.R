context("cyto_gate-internal-helpers")

# .CYTO_GATE_CENTER ------------------------------------------------------------

test_that(".cyto_gate_center", {
  
  # 1D PLOT
  cyto_plot(Va2[[1]], 
          channels = c("FSC-A"), 
          xlim = c(0, 250000),
          density_modal = TRUE)
  
  # 1D RECTANGLEGATE - 1D PLOT
  expect_equal(.cyto_gate_center(rg1,
                                 channels = "FSC-A"),
               matrix(c(25000, 50), ncol = 2,
                      dimnames = list(NULL, c("x","y"))))
  # 2D RECTANGLEGATE - 1D PLOT
  expect_equal(.cyto_gate_center(rg2,
                                 channels = "FSC-A"),
               matrix(c(25000, 50), ncol = 2,
                      dimnames = list(NULL, c("x","y"))))
  # POLYGONGATE - 1D PLOT
  expect_equal(.cyto_gate_center(pg,
                                 channels = "FSC-A"),
               matrix(c(25000, 50), ncol = 2,
                      dimnames = list(NULL, c("x","y"))))
  # ELLIPSOIDGATE - 1D PLOT
  expect_equal(.cyto_gate_center(eg,
                                 channels = "FSC-A"),
               matrix(c(65031, 50), ncol = 2,
                      dimnames = list(NULL, c("x","y"))),
               tolerance = 0.2)
  
  # 1D PLOT
  cyto_plot(Va2[[1]], 
            channels = c("FSC-A","SSC-A"), 
            xlim = c(0, 250000),
            ylim = c(0, 250000))
  
  # 1D RECTANGLEGATE 2D PLOT
  expect_equal(.cyto_gate_center(rg1,
                                 channels = c("FSC-A","SSC-A")),
               matrix(c(25000,125000), ncol = 2, 
                      dimnames = list(NULL, c("x","y"))))
  # 2D RECTANGLEGATE 2D PLOT
  expect_equal(.cyto_gate_center(rg2,
                                 channels = c("FSC-A","SSC-A")),
               matrix(c(25000,25000), ncol = 2,
                      dimnames = list(NULL, c("x","y"))))
  # 2D POLYGONGATE 2D PLOT
  expect_equal(.cyto_gate_center(pg,
                                 channels = c("FSC-A","SSC-A")),
               matrix(c(25000,25000), ncol = 2,
                      dimnames = list(NULL, c("x","y"))))
  # 2D POLYGONGATE 2D PLOT
  expect_equal(.cyto_gate_center(eg,
                                 channels = c("FSC-A","SSC-A")),
               matrix(c(65000,51250), ncol = 2,
                      dimnames = list(NULL, c("x","y"))))
  # QUADGATE 2D PLOT
  expect_equal(.cyto_gate_center(qg,
                                 channels = c("FSC-A","SSC-A")),
               matrix(c(21298, 153702, 
                        153702, 21298,
                        153702, 153702, 
                        21298, 21298), ncol = 2,
                      dimnames = list(NULL, c("x","y"))),
               tolerance = 0.1)
  
})

# .CYTO_GATE_COUNT -------------------------------------------------------------

test_that(".cyto_gate_count", {
  
  # GATES ONLY - VECTOR
  expect_equal(.cyto_gate_count(list(rg1,
                                     pg,
                                     eg,
                                     qg),
                                total = FALSE), 
  c(1,1,1,4))
  
  # GATES ONLY - TOTAL
  expect_equal(.cyto_gate_count(list(rg1,
                                     pg,
                                     eg,
                                     qg),
                                total = TRUE), 
               7)
  
  # NEGATE - VECTOR
  expect_equal(.cyto_gate_count(list(rg1,
                                     pg,
                                     eg),
                                total = FALSE,
                                negate = TRUE), 
               c(1,1,1,1))
  
  # NEGATE - TOTAL
  expect_equal(.cyto_gate_count(list(rg1,
                                     pg,
                                     eg),
                                total = TRUE,
                                negate = TRUE), 
               4)
  
})

# .CYTO_GATE_QUAD_CONVERT ------------------------------------------------------

test_that(".cyto_gate_quad_convert", {
  
  # Q1
  coords <- list(c(-Inf, 50000), c(50000, Inf))
  names(coords) <- c("FSC-A","SSC-A")
  Q1 <- rectangleGate(coords, filterId = "A")
  
  # Q2
  coords <- list(c(50000, Inf), c(50000, Inf))
  names(coords) <- c("FSC-A","SSC-A")
  Q2 <- rectangleGate(coords, filterId = "B")
  
  # Q3
  coords <- list(c(50000, Inf), c(-Inf, 50000))
  names(coords) <- c("FSC-A","SSC-A")
  Q3 <- rectangleGate(coords, filterId = "C")
  
  # Q4
  coords <- list(c(-Inf, 50000), c(-Inf, 50000))
  names(coords) <- c("FSC-A","SSC-A")
  Q4 <- rectangleGate(coords, filterId = "D")
  
  # QUADGATE TO RECTANGLE GATES
  expect_equal(.cyto_gate_quad_convert(qg, 
                                       channels = c("FSC-A","SSC-A")),
               list(Q1,Q2,Q3,Q4))
  
  # RECTANGLEGATES TO QUADGATES
  expect_equal(.cyto_gate_quad_convert(list(Q1,Q2,Q3,Q4), 
                                       channels = c("FSC-A","SSC-A")),
               qg)
  
})

# CLOSE GARPHICS DEVICE
dev.off()
# . CYTO_GATE_COORDS -----------------------------------------------------------

test_that(".cyto_gate_coords", {
  
  # RECTANGLEGATE & POLYGONGATE & QUADGATE
  expect_equal(.cyto_gate_coords(list(rg2,pg, qg), 
                                 channels = c("FSC-A","SSC-A")),
               matrix(c(0,50000,
                        0,50000,
                        50000, 0,
                        50000, 50000,
                        0, 50000,
                        0,0,
                        50000, 50000,
                        50000,50000),
                      ncol = 2,
                      byrow = FALSE,
                      dimnames = list(NULL, c("FSC-A","SSC-A"))))
  
})