context("gating functions")

# .cyto_gate_polygon_draw ------------------------------------------------------

test_that(".cyto_gate_polygon_draw", {
  
  # Alias error
  expect_error(.cyto_gate_polygon_draw(Va2[[1]],
                                       channels = c("FSC-A","SSC-A")),
            "Supply a name for the gated population to the 'alias' argument.",
            fixed = TRUE)
  
  # Too few point selected
  mock_locator <- mock(list("x" = c(50000, 100000), "y" = c(50000, 100000)))
  testthat::with_mock(
    locator = mock_locator,
    expect_error(.cyto_gate_polygon_draw(Va2[[1]],
                                         alias = "Cells",
                                         channels = c("FSC-A","SSC-A"),
                                         plot = FALSE),
             "A minimum of 3 points is required to construct a polygon gate."))
  
  # Valid gate construction
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

# .cyto_gate_interval_draw -----------------------------------------------------
