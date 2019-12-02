context("cyto_gate_draw")

# CYTO_GATE_DRAW ---------------------------------------------------------------

test_that("cyto_gate_draw", {

  # 2D GATES -------------------------------------------------------------------
  
  # 2D GATE TYPES
  gate_types <- c(
    "rectangle",
    "interval", # X AXIS
    "interval", # Y AXIS
    "threshold",
    "boundary",
    "polygon",
    "ellipse",
    "quadrant",
    "web"
  )

  # 2D GATE NAMES
  gate_names <- list(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    c("H", "I", "J", "K"),
    c("L", "M", "N", "O")
  )

  # 2D GATE COORDS
  gate_coords <- list(
    list(
      "x" = c(0, 50000),
      "y" = c(0, 50000)
    ),
    list(
      "x" = c(25000, 150000),
      "y" = c(50000, 50000)
    ), # IRRELEVANT
    list(
      "x" = c(50000, 50000), # IRRELEVANT
      "y" = c(25000, 150000)
    ),
    list(
      "x" = c(25000),
      "y" = c(50000)
    ),
    list(
      "x" = c(200000),
      "y" = c(200000)
    ),
    list(
      "x" = c(0, 50000, 50000, 0),
      "y" = c(0, 0, 50000, 50000)
    ),
    list(
      "x" = c(25000, 50000, 75000, 50000),
      "y" = c(49000, 25000, 50000, 77000)
    ),
    list(
      "x" = c(50000),
      "y" = c(50000)
    ),
    list(
      list("x" = c(100000),
           "y" = c(100000)),
      list("x" = c(10000),
           "y" = c(40000)),
      list("x" = c(150000),
           "y" = c(0)),
      list("x" = c(240000),
           "y" = c(150000)),
      list("x" = c(50000),
           "y" = c(250000))
    )
  )

  # 2D GATE AXIS
  gate_axes <- c(
    "x",
    "x",
    "y",
    "x",
    "x",
    "x",
    "x",
    "x",
    "x"
  )

  # GATES
  gates <- list(
    rg2,
    ig,
    igy,
    tg,
    bg,
    pg,
    eg,
    qg,
    list(wg1,wg2,wg3,wg4)
  )

  # TESTS
  mapply(
    function(gate_name,
             gate_type,
             gate_coord,
             gate_axis,
             gate) {
      # INPUT COORDS
      if(gate_type == "web"){
        mock_locator <- mock(gate_coord[[1]],
                             gate_coord[[2]],
                             gate_coord[[3]],
                             gate_coord[[4]],
                             gate_coord[[5]])
        gate <- filters(gate)
      }else{
        mock_locator <- mock(gate_coord)
        gate <- filters(list(gate))
      }
      testthat::with_mock(locator = mock_locator,{
        expect_equivalent(
          cyto_gate_draw(fs,
            channels = c("FSC-A", "SSC-A"),
            alias = gate_name,
            type = gate_type,
            axis = gate_axis
          ),
          gate,
          tolerance = 0.0000000001)})
    },
    gate_names,
    gate_types,
    gate_coords,
    gate_axes,
    gates
  )
  
  # 1D GATES -------------------------------------------------------------------
  
  # 1D GATE TYPES
  gate_types <- c("")
  
  
})
