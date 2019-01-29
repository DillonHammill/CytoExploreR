context("Gating Functions")

## gate_polygon_draw -----------------------------------------------------------------------

test_that("gate_polygon_draw", {
  expect_error(gate_polygon_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")

  dp <- gate_polygon_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(dp, "filters")
  expect_s4_class(dp[[1]], "polygonGate")
  expect_equal(dp[[1]], pg)
  expect_equal(parameters(dp[[1]]), parameters(pg))
})

## gate_rectangle_draw ---------------------------------------------------------------------

test_that("gate_rectangle_draw", {
  expect_error(gate_rectangle_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")

  dr <- gate_rectangle_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(dr, "filters")
  expect_s4_class(dr[[1]], "rectangleGate")
  expect_equal(dr[[1]], rg)
  expect_equal(parameters(dr[[1]]), parameters(rg))
})

## gate_interval_draw ----------------------------------------------------------------------

test_that("gate_interval_draw", {

  # 1-D ---------------------------------------------------------------
  di <- gate_interval_draw(fs[[1]], alias = "Cells", channels = "FSC-A")

  expect_s4_class(di, "filters")
  expect_s4_class(di[[1]], "rectangleGate")
  expect_equal(di[[1]], igx)
  expect_equal(parameters(di[[1]]), parameters(igx))

  expect_error(gate_interval_draw(fs[[1]], alias = "Cells", channels = "FSC-A", axis = "y"), "Cannot gate y axis if a single channel is supplied.")
  expect_error(gate_interval_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")

  # 2-D x axis --------------------------------------------------------
  di <- gate_interval_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(di, "filters")
  expect_s4_class(di[[1]], "rectangleGate")
  expect_equal(di[[1]], ig)
  expect_equal(parameters(di[[1]]), parameters(ig))

  # 2-D y axis --------------------------------------------------------
  di <- gate_interval_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05, axis = "y")

  expect_s4_class(di, "filters")
  expect_s4_class(di[[1]], "rectangleGate")
  expect_equal(di[[1]], igy)
  expect_equal(parameters(di[[1]]), parameters(igy))
})

## gate_threshold_draw ------------------------------------------------------

test_that("gate_threshold_draw", {
  expect_error(gate_threshold_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")
  expect_error(gate_threshold_draw(fs[[1]], alias = c("A", "B"), channels = "FSC-A"), "Multiple threshold gates are not supported.")

  # 1-D ---------------------------------------------------------------
  dt <- gate_threshold_draw(fs[[1]], alias = "Cells", channels = "FSC-A")

  expect_s4_class(dt, "filters")
  expect_s4_class(dt[[1]], "rectangleGate")
  expect_equal(dt[[1]], tg1)
  expect_equal(parameters(dt[[1]]), parameters(tg1))

  # 2-D ---------------------------------------------------------------
  dt <- gate_threshold_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(dt, "filters")
  expect_s4_class(dt[[1]], "rectangleGate")
  expect_equal(dt[[1]], tg)
  expect_equal(parameters(dt[[1]]), parameters(tg))
})

## gate_boundary_draw -------------------------------------------------------

test_that("gate_boundary_draw", {
  expect_error(gate_boundary_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")
  expect_error(gate_boundary_draw(fs[[1]], alias = c("A", "B"), channels = "FSC-A"), "Multiple boundary gates are not supported.")

  # 1-D ---------------------------------------------------------------
  db <- gate_boundary_draw(fs[[1]], alias = "Cells", channels = "FSC-A")

  expect_s4_class(db, "filters")
  expect_s4_class(db[[1]], "rectangleGate")
  expect_equal(db[[1]], bg1)
  expect_equal(parameters(db[[1]]), parameters(bg1))

  # 2-D ---------------------------------------------------------------
  db <- gate_boundary_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(db, "filters")
  expect_s4_class(db[[1]], "rectangleGate")
  expect_equal(db[[1]], bg)
  expect_equal(parameters(db[[1]]), parameters(bg))
})

## gate_ellipse_draw --------------------------------------------------------

test_that("gate_ellipse_draw", {
  expect_error(gate_ellipse_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")

  de <- gate_ellipse_draw(fs[[1]], alias = "Cells", channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(de, "filters")
  expect_s4_class(de[[1]], "ellipsoidGate")
  expect_equal(de[[1]]@mean, eg@mean)
  expect_equal(round(de[[1]]@cov, 8), round(eg@cov, 8))
})

## gate_quadrant_draw -------------------------------------------------------

test_that("gate_quadrant_draw", {
  expect_error(gate_quadrant_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")
  expect_error(gate_quadrant_draw(fs[[1]], alias = "A", channels = "FSC-A"), "Supply 4 population names as the alias argument to construct a set of quadrant gates.")

  dq <- gate_quadrant_draw(fs[[1]], alias = c("A", "B", "C", "D"), channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(dq, "filters")
  expect_length(dq, 4)
  expect_equal(as.vector(sapply(dq, class)), rep("rectangleGate", 4))
  expect_equal(qg, dq)
})

## gate_web_draw -------------------------------------------------------------

test_that("gate_web_draw", {
  expect_error(gate_web_draw(fs[[1]], channels = "FSC-A"), "Please supply a name for a the gated population as the alias argument.")

  dw <- gate_web_draw(fs[[1]], alias = c("A", "B", "C", "D", "E", "F", "G", "H"), channels = c("FSC-A", "SSC-A"), display = 0.05)

  expect_s4_class(dw, "filters")
  expect_length(dw, 8)
  expect_equal(as.vector(sapply(dw, class)), rep("polygonGate", 8))
  expect_equal(wg, dw, tolerance = 0.01)
})
