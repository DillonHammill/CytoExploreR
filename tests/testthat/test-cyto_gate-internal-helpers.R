context("cyto_gate-internal-helpers")

# .CYTO_GATE_CENTER ------------------------------------------------------------

test_that(".cyto_gate_center", {
  
  # REACTNGLE GATES
  
  # POLYGON GATES
  
  # ELLIPSOID GATES
  
  # QUADRANT GATES
  
  # FILTERS
  
  # LIST
  
})

# .CYTO_GATE_COUNT -------------------------------------------------------------

test_that(".cyto_gate_count", {
  
  expect_equal(.cyto_gate_count(list(rg1, rg2, pg, eg, qg),
                                negate = FALSE,
                                total = FALSE),
               c(1,1,1,1,4))
  
  expect_equal(.cyto_gate_count(list(rg1, rg2, pg, eg),
                                negate = TRUE,
                                total = TRUE),
               5)
  
})

# .CYTO_GATE_QUAD_CONVERT ------------------------------------------------------

test_that(".cyto_gate_quad_convert", {
  
  r1 <- rectangleGate(list("FSC-A" = c(-Inf, 50000),
                           "SSC-A" = c(50000, Inf)),
                      filterId = "H")
  r2 <- rectangleGate(list("FSC-A" = c(50000, Inf),
                           "SSC-A" = c(50000, Inf)),
                      filterId = "I")
  r3 <- rectangleGate(list("FSC-A" = c(50000, Inf),
                           "SSC-A" = c(-Inf, 50000)),
                      filterId = "J")
  r4 <- rectangleGate(list("FSC-A" = c(-Inf, 50000),
                           "SSC-A" = c(-Inf, 50000)),
                      filterId = "K")
  r <- list("H" = r1,
            "I" = r2,
            "J" = r3,
            "K" = r4)
  
  expect_equal(.cyto_gate_quad_convert(qg, 
                                       channels = c("FSC-A","SSC-A")),
               r)
  
  expect_equal(.cyto_gate_quad_convert(r,
                                       channels = c("FSC-A", "SSC-A")), qg)
  
})

# .CYTO_GATE_COORDS ------------------------------------------------------------

test_that(".cyto_gate_coords", {
  
  coords <- matrix(c(0.00, 0.00,
                     50000.00, 50000.00,
                     0.00, 0.00,
                     50000.00, 50000.00,
                     0.00, 0.00,
                     50000.00, 0.00,
                     50000.00, 50000.00,
                     0.00, 50000.00,
                     75005.00, 50250.00,
                     74799.71, 53574.81,
                     74187.21, 56845.02,
                     73177.55, 60006.94,
                     71787.32, 63008.66,
                     70039.35, 65800.87,
                     67962.33, 68337.75,
                     65590.36, 70577.62,
                     62962.41, 72483.71,
                     60121.61, 74024.73,
                     57114.61, 75175.36,
                     53990.80, 75916.73,
                     50801.45, 76236.64,
                     47598.94, 76129.86,
                     44435.86, 75598.13,
                     41364.15, 74650.18,
                     38434.23, 73301.58,
                     35694.22, 71574.48,
                     33189.12, 69497.23,
                     30960.04, 67103.94,
                     29043.61, 64433.91,
                     27471.27, 61530.98,
                     26268.86, 58442.81,
                     25456.11, 55220.12,
                     25046.38, 51915.83,
                     25046.38, 48584.17,
                     25456.11, 45279.88,
                     26268.86, 42057.19,
                     27471.27, 38969.02,
                     29043.61, 36066.09,
                     30960.04, 33396.06,
                     33189.12, 31002.77,
                     35694.22, 28925.52,
                     38434.23, 27198.42,
                     41364.15, 25849.82,
                     44435.86, 24901.87,
                     47598.94, 24370.14,
                     50801.45, 24263.36,
                     53990.80, 24583.27,
                     57114.61, 25324.64,
                     60121.61, 26475.27,
                     62962.41, 28016.29,
                     65590.36, 29922.38,
                     67962.33, 32162.25,
                     70039.35, 34699.13,
                     71787.32, 37491.34,
                     73177.55, 40493.06,
                     74187.21, 43654.98,
                     74799.71, 46925.19,
                     75005.00, 50250.00,
                     50000.00, 50000.00), 
                   ncol = 2, 
                   byrow = TRUE)
  colnames(coords) <- c("FSC-A", "SSC-A")
  
  expect_equal(.cyto_gate_coords(list(rg1,
                                      rg2,
                                      pg,
                                      eg,
                                      qg),
                                 channels = "FSC-A"),
               coords[, 1, drop = FALSE], tolerance = 0.01)
  
  expect_equal(.cyto_gate_coords(list(rg1,
                                      rg2,
                                      pg,
                                      eg,
                                      qg),
                                 channels = c("FSC-A", "SSC-A")),
               coords, tolerance = 0.01)
  
})