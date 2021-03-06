context("cyto_gate_draw")

# CYTO_GATE_DRAW ---------------------------------------------------------------

# WEB GATE NAMES DON'T MATCH

test_that("cyto_gate_draw", {

  # GATE TYPES
  gate_types <- c(
    "rectangle",
    "interval", # X AXIS
    "interval", # Y AXIS
    "threshold",
    "boundary",
    "polygon",
    "ellipse",
    "quadrant",
    "web",
    "interval",
    "threshold",
    "boundary"
  )

  # GATE NAMES
  gate_names <- list(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    c("H", "I", "J", "K"),
    c("L", "M", "N", "O"),
    "P",
    "Q",
    "R"
  )

  # GATE COORDS
  gate_coords <- list(
    list(
      list("x" = 0,
           "y" = 0),
      list("x" = 50000,
           "y" = 50000)
    ),
    list(
      list("x" = 25000,
           "y" = 50000),
      list("x" = 150000,
           "y" = 50000)
    ), 
    list(
      list("x" = 50000,
           "y" = 25000),
      list("x" = 50000,
           "y" = 150000)
    ),
    list(
      "x" = c(25000),
      "y" = c(50000)
    ),
    list("x" = c(200000),
         "y" = c(200000)
    ),
    list(
      list("x" = 0,
           "y" = 0),
      list("x" = 50000,
           "y" = 0),
      list("x" = 50000,
           "y" = 50000),
      list("x" = 0,
           "y" = 50000)
    ),
    list(
      list("x" = 25000,
           "y" = 49000),
      list("x" = 50000,
           "y" = 25000),
      list("x" = 75000,
           "y" = 50000),
      list("x" = 50000,
           "y" = 77000)
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
    ),
    list(
      list("x" = 0,
           "y" = 50),
      list("x" = 50000,
           "y" = 50)
    ),
    list(
      "x" = c(25000),
      "y" = c(50)
    ),
    list(
       "x" = c(200000),
       "y" = c(50)
    )
  )

  # GATE AXIS
  gate_axes <- c(
    "x",
    "x",
    "y",
    "x",
    "x",
    "x",
    "x",
    "x",
    "x",
    "x",
    "x",
    "x"
  )

  # GATES
  gates <- list(
    list("A" = rg2),
    list("B" = ig),
    list("C" = igy),
    list("D" = tg),
    list("E" = bg),
    list("F" = pg),
    list("G" = eg),
    list("H|I|J|K" = qg),
    list("L" = wg1,
         "M" = wg2,
         "N" = wg3,
         "O" = wg4),
    list("P" = rg1),
    list("Q" = tg1),
    list("R" = bg1)
  )

  # GATE CHANNELS
  gate_channels <- c(rep(list(c("FSC-A","SSC-A")), 9),
                     list("FSC-A"),
                     list("FSC-A"),
                     list("FSC-A"))
  
  # TESTS
  mapply(
    function(gate_channels,
             gate_name,
             gate_type,
             gate_coord,
             gate_axis,
             gate) {
      # INPUT COORDS
      if(gate_type %in% c("threshold","boundary", "quadrant")){
        mock_locator <- mock(gate_coord)
      }else if(gate_type %in% c("interval", "rectangle")){
        mock_locator <- mock(gate_coord[[1]],
                             gate_coord[[2]])
      }else if(gate_type %in% c("polygon", "ellipse")){
        mock_locator <- mock(gate_coord[[1]],
                             gate_coord[[2]],
                             gate_coord[[3]],
                             gate_coord[[4]], cycle = TRUE)
      }else if(gate_type == "web"){
        mock_locator <- mock(gate_coord[[1]],
                             gate_coord[[2]],
                             gate_coord[[3]],
                             gate_coord[[4]],
                             gate_coord[[5]])
      }
      # GATE
      gate <- list(gate)
      names(gate) <- "Combined Events"
      # print(gate)
      testthat::with_mock(locator = mock_locator,{
        expect_equivalent(
          cyto_gate_draw(fs,
            channels = gate_channels,
            alias = gate_name,
            type = gate_type,
            axis = gate_axis,
            display = 100
          ),
          gate,
          tolerance = 0.1)})
    },
    gate_channels,
    gate_names,
    gate_types,
    gate_coords,
    gate_axes,
    gates
  )
  
})
