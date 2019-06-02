context("gating functions")

cyto_plot(Va2[[1]],
          channels = c("FSC-A","SSC-A"))

xmin <- par("usr")[1]
xmax <- par("usr")[2]
ymin <- par("usr")[3]
ymax <- par("usr")[4]

# .cyto_gate_polygon_draw ------------------------------------------------------

test_that(".cyto_gate_polygon_draw", {
  
  # Alias error
  expect_error(.cyto_gate_polygon_draw(Va2[[1]],
                                       channels = c("FSC-A","SSC-A")),
        "Supply a name for the gated population(s) to the 'alias' argument.",
        fixed = TRUE)
  
  # Too few points selected
  mock_locator <- mock(list("x" = c(50000, 100000), "y" = c(50000, 100000)))
  testthat::with_mock(
    locator = mock_locator,
    expect_error(.cyto_gate_polygon_draw(Va2[[1]],
                                         alias = "Cells",
                                         channels = c("FSC-A","SSC-A"),
                                         plot = FALSE),
             "A minimum of 3 points is required to construct a polygon gate."))
  
  # 2D -------------------------------------------------------------------------
  pg <- matrix(c(
    50000,
    75000,
    100000,
    75000,
    50000,
    50000,
    100000,
    100000
  ),
  ncol = 2,
  nrow = 4,
  byrow = FALSE
  )
  colnames(pg) <- c("FSC-A", "SSC-A")
  pg <- polygonGate(filterId = "Cells", .gate = pg)
  pg <- filters(list(pg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      75000,
      100000,
      75000
    ),
    "y" = c(
      50000,
      50000,
      100000,
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_polygon_draw(Va2[[1]],
        alias = "Cells",
        channels = c("FSC-A", "SSC-A"),
        plot = FALSE
      ),
      pg
    )
  )
  
})

# .cyto_gate_rectangle_draw ----------------------------------------------------

test_that(".cyto_gate_rectangle_draw", {
  
  # Alias error
  expect_error(.cyto_gate_rectangle_draw(Va2[[1]],
                                         channels = c("FSC-A","SSC-A")),
        "Supply a name for the gated population(s) to the 'alias' argument.",
        fixed = TRUE)

  # 2D -------------------------------------------------------------------------
  rg <- matrix(c(
    50000,
    100000,
    50000,
    100000
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A", "SSC-A")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      100000
    ),
    "y" = c(
      50000,
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_rectangle_draw(Va2[[1]],
        alias = "Cells",
        channels = c("FSC-A", "SSC-A"),
        plot = FALSE
      ),
      rg
    )
  )
  
})

# .cyto_gate_interval_draw -----------------------------------------------------

test_that(".cyto_gate_interval_draw", {
  
  # Alias error
  expect_error(.cyto_gate_interval_draw(Va2[[1]],
                                         channels = c("FSC-A","SSC-A")),
          "Supply a name for the gated population(s) to the 'alias' argument.",
          fixed = TRUE)
  
  # 1D - x axis ----------------------------------------------------------------
  rg <- matrix(c(
    50000,
    100000
  ),
  ncol = 1,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      100000
    ),
    "y" = c(
      50000,
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_interval_draw(Va2[[1]],
                               alias = "Cells",
                               channels = c("FSC-A"),
                               axis = "x",
                               plot = FALSE
      ),
      rg
    )
  )

  # 1D - y axis (error) --------------------------------------------------------
  rg <- matrix(c(
    50000,
    100000
  ),
  ncol = 1,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      100000
    ),
    "y" = c(
      50000,
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_error(
      .cyto_gate_interval_draw(Va2[[1]],
                               alias = "Cells",
                               channels = c("FSC-A"),
                               axis = "y",
                               plot = FALSE
      ),
      "Cannot gate y axis if a single channel is supplied."
    )
  )
  
  # 2D - x axis ----------------------------------------------------------------
  rg <- matrix(c(
    50000,
    100000,
    -Inf,
    Inf
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A", "SSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      100000
    ),
    "y" = c(
      50000,
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_interval_draw(Va2[[1]],
                                alias = "Cells",
                                channels = c("FSC-A", "SSC-A"),
                                axis = "x",
                                plot = FALSE
      ),
      rg
    )
  )
  
  # 2D - y axis ----------------------------------------------------------------
  rg <- matrix(c(
    -Inf,
    Inf,
    50000,
    100000
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A", "SSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      100000
    ),
    "y" = c(
      50000,
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_interval_draw(Va2[[1]],
                               alias = "Cells",
                               channels = c("FSC-A", "SSC-A"),
                               axis = "y",
                               plot = FALSE
      ),
      rg
    )
  )
  
})

# .cyto_gate_threshold_draw ----------------------------------------------------

test_that(".cyto_gate_threshold_draw", {
  
  # Multiple alias error
  expect_error(.cyto_gate_threshold_draw(Va2[[1]],
                                         alias = c("A","B","C"),
                                         channels = c("FSC-A","SSC-A"),
                                         plot = FALSE),
          "Multiple threshold gates are not supported.",
          fixed = TRUE)
  
  # Alias error
  expect_error(.cyto_gate_threshold_draw(Va2[[1]],
                                         channels = c("FSC-A","SSC-A")),
          "Supply a name for the gated population(s) to the 'alias' argument.",
          fixed = TRUE)
  
  # 1D -------------------------------------------------------------------------
  rg <- matrix(c(
    50000,
    Inf
  ),
  ncol = 1,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000
    ),
    "y" = c(
      50000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_threshold_draw(Va2[[1]],
                               alias = "Cells",
                               channels = c("FSC-A"),
                               plot = FALSE
      ),
      rg
    )
  )
  
  # 2D -------------------------------------------------------------------------
  rg <- matrix(c(
    50000,
    Inf,
    50000,
    Inf
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A", "SSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      50000
    ),
    "y" = c(
      50000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_threshold_draw(Va2[[1]],
                                alias = "Cells",
                                channels = c("FSC-A", "SSC-A"),
                                plot = FALSE
      ),
      rg
    )
  )
  
})

# .cyto_gate_boundary_draw -----------------------------------------------------

test_that(".cyto_gate_boundary_draw", {
  
  # Multiple alias error
  expect_error(.cyto_gate_boundary_draw(Va2[[1]],
                                         alias = c("A","B","C"),
                                         channels = c("FSC-A","SSC-A"),
                                         plot = FALSE),
               "Multiple boundary gates are not supported.",
               fixed = TRUE)
  
  # Alias error
  expect_error(.cyto_gate_boundary_draw(Va2[[1]],
                                         channels = c("FSC-A","SSC-A")),
          "Supply a name for the gated population(s) to the 'alias' argument.",
          fixed = TRUE)
  
  # 1D -------------------------------------------------------------------------
  rg <- matrix(c(
    -Inf,
    100000
  ),
  ncol = 1,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      100000
    ),
    "y" = c(
      50
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_boundary_draw(Va2[[1]],
                                alias = "Cells",
                                channels = c("FSC-A"),
                                plot = FALSE
      ),
      rg
    )
  )
  
  # 2D -------------------------------------------------------------------------
  rg <- matrix(c(
    -Inf,
    100000,
    -Inf,
    100000
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(rg) <- c("FSC-A", "SSC-A")
  rownames(rg) <- c("min","max")
  rg <- rectangleGate(filterId = "Cells", .gate = rg)
  rg <- filters(list(rg))
  mock_locator <- mock(list(
    "x" = c(
      100000
    ),
    "y" = c(
      100000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_boundary_draw(Va2[[1]],
                                alias = "Cells",
                                channels = c("FSC-A", "SSC-A"),
                                plot = FALSE
      ),
      rg
    )
  )
  
})

# .cyto_gate_ellipse_draw ------------------------------------------------------

test_that(".cyto_gate_ellipse_draw", {
  
  # NOTE: Horizontal Major axis not supported!
  
  # Alias error
  expect_error(.cyto_gate_ellipse_draw(Va2[[1]],
                                       channels = c("FSC-A","SSC-A")),
          "Supply a name for the gated population(s) to the 'alias' argument.",
          fixed = TRUE)
  
  # 2D -------------------------------------------------------------------------
  cov <- matrix(c(850000000,
                  400000000,
                  400000000,
                  850000000),
                ncol = 2,
                nrow = 2,
                byrow = FALSE)
  rownames(cov) <- c("FSC-A","SSC-A")
  colnames(cov) <- c("FSC-A", "SSC-A")
  cnt <- c(75000, 75000)
  eg <- ellipsoidGate(filterId = "Cells",
                      mean = cnt,
                      .gate = cov)
  eg <- filters(list(eg))
  mock_locator <- mock(list(
    "x" = c(
      50000,
      90000,
      100000,
      60000
    ),
    "y" = c(
      50000,
      60000,
      100000,
      90000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_ellipse_draw(Va2[[1]],
                              alias = "Cells",
                              channels = c("FSC-A", "SSC-A"),
                              plot = FALSE
      ),
      eg
    )
  )
  
})

# .cyto_gate_quadrant_draw -----------------------------------------------------

test_that(".cyto_gate_quadrant_draw", {
  
  # Alias error
  expect_error(.cyto_gate_quadrant_draw(Va2[[1]],
                                       channels = c("FSC-A","SSC-A")),
          "Supply a name for the gated population(s) to the 'alias' argument.",
          fixed = TRUE)
  
  # length(alias) != 4
  expect_error(.cyto_gate_quadrant_draw(Va2[[1]],
                                        alias = c("A","B"),
                                        channels = c("FSC-A","SSC-A"),
                                        plot = FALSE),
               "'alias' must contain 4 population names for quadrant gates.",
               fixed = TRUE)
  
  # 2D -------------------------------------------------------------------------
  # Q1 - bottom left
  q1 <- matrix(c(
    -Inf,
    75000,
    -Inf,
    75000
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(q1) <- c("FSC-A", "SSC-A")
  rownames(q1) <- c("min","max")
  q1 <- rectangleGate(filterId = "A", .gate = q1)
  
  # Q2
  q2 <- matrix(c(
    75000,
    Inf,
    75000,
    -Inf
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(q2) <- c("FSC-A", "SSC-A")
  rownames(q2) <- c("min","max")
  q2 <- rectangleGate(filterId = "B", .gate = q2)
  
  # Q3
  q3 <- matrix(c(
    75000,
    Inf,
    75000,
    Inf
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(q3) <- c("FSC-A", "SSC-A")
  rownames(q3) <- c("min","max")
  q3 <- rectangleGate(filterId = "C", .gate = q3)
  
  # Q4
  q4 <- matrix(c(
    75000,
    -Inf,
    75000,
    Inf
  ),
  ncol = 2,
  nrow = 2,
  byrow = FALSE
  )
  colnames(q4) <- c("FSC-A", "SSC-A")
  rownames(q4) <- c("min","max")
  q4 <- rectangleGate(filterId = "D", .gate = q4)
  
  # combined quadrants
  q <- filters(list(q1,q2,q3,q4))
  
  mock_locator <- mock(list(
    "x" = c(
      75000
    ),
    "y" = c(
      75000
    )
  ))
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_quadrant_draw(Va2[[1]],
                              alias = c("A","B","C","D"),
                              channels = c("FSC-A","SSC-A"),
                              plot = FALSE
      ),
      q
    )
  )
  
})

# .cyto_gate_web_draw ----------------------------------------------------------

test_that(".cyto_gate_web_draw", {
  
  # Alias error
  expect_error(.cyto_gate_web_draw(Va2[[1]],
                                        channels = c("FSC-A","SSC-A")),
          "Supply a name for the gated population(s) to the 'alias' argument.",
          fixed = TRUE)
  
  # W1
  w1 <- matrix(c(
    133872.17,
    xmin,
    xmin,
    167267.80,
    155121.92,
    81904.57,
    ymin,
    ymin
  ),
  ncol = 2,
  nrow = 4,
  byrow = FALSE
  )
  colnames(w1) <- c("FSC-A", "SSC-A")
  w1 <- polygonGate(filterId = "A", .gate = w1)

  # W2
  w2 <- matrix(c(
    133872.2,
    167267.8,
    xmax,
    xmax,
    155121.92,
    ymin,
    ymin,
    216622.7
  ),
  ncol = 2,
  nrow = 4,
  byrow = FALSE
  )
  colnames(w2) <- c("FSC-A", "SSC-A")
  w2 <- polygonGate(filterId = "B", .gate = w2)
  
  # W3
  w3 <- matrix(c(
    133872.17,
    xmax,
    xmax,
    77791.76,
    155121.9,
    216622.7,
    ymax,
    ymax
  ),
  ncol = 2,
  nrow = 4,
  byrow = FALSE
  )
  colnames(w3) <- c("FSC-A", "SSC-A")
  w3 <- polygonGate(filterId = "C", .gate = w3)
  
  # W4
  w4 <- matrix(c(
    133872.17,
    77791.76,
    xmin,
    xmin,
    155121.92,
    ymax,
    ymax,
    81904.57
  ),
  ncol = 2,
  nrow = 4,
  byrow = FALSE
  )
  colnames(w4) <- c("FSC-A", "SSC-A")
  w4 <- polygonGate(filterId = "D", .gate = w4)
  
  w <- filters(list(w1,w2,w3,w4))
  
  mock_locator <- mock(list("x" = 125000, 
                            "y" = 125000),
                       list("x" = xmin,
                            "y" = 80000),
                       list("x" = 170000,
                            "y" = ymin),
                       list("x" = xmax,
                            "y" = 210000),
                       list("x" = 80000,
                            "y" = ymax))
  
  testthat::with_mock(
    locator = mock_locator,
    expect_equal(
      .cyto_gate_web_draw(Va2[[1]],
                               alias = c("A","B","C","D"),
                               channels = c("FSC-A","SSC-A"),
                               plot = FALSE
      ),
      w, tolerance = 0.2
    )
  )
  
})
