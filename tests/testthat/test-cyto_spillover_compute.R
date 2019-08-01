context("cyto_spillover_compute")

# cyto_spillover_compute -------------------------------------------------------

test_that("cyto_stats_compute", {
  
  # flowSet Un-Transformed Error -----------------------------------------------
  expect_error(
    cyto_spillover_compute(cyto_extract(gsc, 
                                        parent = "Single Cells")),
    paste("Fluorescent channels MUST be transformed. Supply the transformerList",
          "to 'axes_trans'."),
    fixed = TRUE
  )
  
  # GatingSet Un-Transformed Error ---------------------------------------------
  expect_error(
    cyto_spillover_compute(gsc,
                           parent = "Single Cells"),
    "Transform ALL fluorescent channels using cyto_transform().",
    fixed = TRUE
  )
  
  # GatingSet Universal Unstained Compensation Control -------------------------
  mock_locator <- mock(list("x" = c(2.099,4.422),
                            "y" = c(50,50)),
                       list("x" = c(2.010,4.1561),
                            "y" = c(50,50)),
                       list("x" = c(1.8297,3.9772),
                            "y" = c(50,50)),
                       list("x" = c(2.5032,4.6266),
                            "y" = c(50,50)),
                       list("x" = c(2.5282,4.1284),
                            "y" = c(50,50)),
                       list("x" = c(2.7357,4.2116),
                            "y" = c(50,50)))
    
  mock_menu <- mock(4,10,11,9,1,2,12)
  
  spill <- read.csv("Universal-Spillover-Matrix.csv", 
                    header = TRUE,
                    row.names = 1)
  colnames(spill) <- rownames(spill)
  spill <- as.matrix(spill)
  
  testthat::with_mock(
    locator = mock_locator,
    menu = mock_menu,
    expect_equal(
      cyto_spillover_compute(gsct,
                             parent = "Single Cells",
                             spillover = "Universal-Spillover-Matrix.csv"),
      spill, tolerance = 0.01
    )  
  )
  
  expect_true(file_wd_check("Universal-Spillover-Matrix.csv"))
  
  # GatingSet Internal Unstained Reference Populations -------------------------
  mock_locator <- mock(list("x" = c(-0.5,1.2),
                            "y" = c(50,50)),
                       list("x" = c(2.099,4.422),
                            "y" = c(50,50)),
                       list("x" = c(-0.5,1.2),
                            "y" = c(50,50)),
                       list("x" = c(2.010,4.1561),
                            "y" = c(50,50)),
                       list("x" = c(-0.5,1.2),
                            "y" = c(50,50)),
                       list("x" = c(1.8297,3.9772),
                            "y" = c(50,50)),
                       list("x" = c(-0.5,1.2),
                            "y" = c(50,50)),
                       list("x" = c(2.5032,4.6266),
                            "y" = c(50,50)),
                       list("x" = c(-0.5,1.2),
                            "y" = c(50,50)),
                       list("x" = c(2.5282,4.1284),
                            "y" = c(50,50)),
                       list("x" = c(-0.5,1.2),
                            "y" = c(50,50)),
                       list("x" = c(2.7357,4.2116),
                            "y" = c(50,50)))
  
  mock_menu <- mock(4,10,11,9,1,2)

  spill <- read.csv("Internal-Spillover-Matrix.csv", 
                    header = TRUE,
                    row.names = 1)
  colnames(spill) <- rownames(spill)
  spill <- as.matrix(spill)
  
  testthat::with_mock(
    locator = mock_locator,
    menu = mock_menu,
    expect_equal(
      cyto_spillover_compute(gsct[-7],
                             parent = "Single Cells",
                             spillover = "Internal-Spillover-Matrix.csv"),
      spill, tolerance = 0.01
    )  
  )
  
  expect_true(file_wd_check("Internal-Spillover-Matrix.csv"))
  
})

# Close plot windows
graphics.off()
  
unlink("Compensation-Channels.csv")
unlink("Spillover-Matrix.csv")