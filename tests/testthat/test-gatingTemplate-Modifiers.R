context("gatingTemplate modifiers")

gs3 <- clone(gs)

# gate_remove ------------------------------------------------------------------

test_that("gate_remove", {
  
  expect_error(gate_remove(gs3), 
               "Please supply the name of the population to be removed.")
  expect_error(gate_remove(gs3, alias = "Test"), 
               "Supplied alias does not exist in the GatingSet.")
  
  expect_error(gate_remove(gs3, alias = "T Cells"), 
        "Supply the name of the gatingTemplate csv file to remove the gate.")
  expect_error(gate_remove(gs3, 
                           alias = "T Cells", 
                           gatingTemplate = "Test.csv"), 
               "Test.csv is not in this working directory.")
  
  gate_remove(gs3, 
              alias = "CD4 T Cells", 
              gatingTemplate = "Activation-gatingTemplate.csv")
  expect_equal(basename(getNodes(gs3)), 
               c("root",
                 "Cells",
                 "Single Cells",
                 "Live Cells",
                 "Dendritic Cells",
                 "T Cells",
                 "CD8 T Cells",
                 "CD69+ CD8 T Cells"))
  expect_equal(basename(gatingTemplate("Activation-gatingTemplate.csv")@nodes),
               c("root",
                 "Cells",
                 "Single Cells",
                 "Live Cells",
                 "T Cells",
                 "Dendritic Cells",
                 "CD8 T Cells",
                 "CD69+ CD8 T Cells"))

  gate_remove(gs3, 
              alias = c("Dendritic Cells", "T Cells"), 
              gatingTemplate = "Activation-gatingTemplate.csv")
  expect_equal(basename(getNodes(gs3)), 
               c("root",
                 "Cells",
                 "Single Cells",
                 "Live Cells"))
  expect_equal(basename(gatingTemplate("Activation-gatingTemplate.csv")@nodes), 
               c("root",
                 "Cells",
                 "Single Cells",
                 "Live Cells"))
  
  write.csv(gtf, "Activation-gatingTemplate.csv", row.names = FALSE)
  
})

# gate_extract -----------------------------------------------------------------

test_that("gate_extract", {
  
  expect_error(gate_extract(alias = "T Cells"),
               "Please supply the name of the parent population.")
  expect_error(gate_extract(parent = "T Cells"), 
               "Please supply the name(s) of the alias to extract.", 
               fixed = TRUE)
  expect_error(gate_extract(parent = "T Cells", 
                            alias = "CD4 T Cells"), 
      "Please supply the name of the gatingTemplate to extract gates from.")
  expect_error(gate_extract(parent = "T Cells", 
                            alias = "CD4 T Cells", 
                            gatingTemplate = "Test.csv"), 
               "Test.csv is not in this working directory.")
  
  gate <- gate_extract(parent = "T Cells", 
                       alias = c("CD4 T Cells","CD8 T Cells"), 
                       gatingTemplate = "Activation-gatingTemplate.csv")
  
  coords <- matrix(c(1.7257,3.5148, -0.1326, 1.9919), ncol = 2, nrow = 2)
  colnames(coords) <- c("Alexa Fluor 700-A","Alexa Fluor 488-A")
  rownames(coords) <- c("min","max")
  test1 <- rectangleGate(filterId = "CD4 T Cells", .gate = coords)
  
  coords <- matrix(c(-0.3480,1.8883, 2.1335, 4.0759), ncol = 2, nrow = 2)
  colnames(coords) <- c("Alexa Fluor 700-A","Alexa Fluor 488-A")
  rownames(coords) <- c("min","max")
  test2 <- rectangleGate(filterId = "CD8 T Cells", .gate = coords)
  
  expect_equal(gate, list(list(filters(list(test1))), 
                          list(filters(list(test2)))), 
               tolerance = 0.01)
  
})

# gate_edit --------------------------------------------------------------------

test_that("gate_edit", {
  
  gs3 <- clone(gs)
  
  expect_error(gate_edit(gs), 
               "Please supply the name of the parent population.")
  expect_error(gate_edit(gs, 
                         parent = "Test"), 
               "Supplied parent does not exist in the GatingSet.")
  expect_error(gate_edit(gs, 
                         parent = "T Cells"), 
               "Please supply the name(s) of the gates to edit to 'alias'.", 
               fixed =  TRUE)
  expect_error(gate_edit(gs, 
                         parent = "T Cells", 
                         alias = "Test"), 
               "Supplied alias does not exist in the GatingSet.")
  expect_error(gate_edit(gs, 
                         parent = "T Cells", 
                         alias = "CD4 T Cells"), 
        "Please supply the name of gatingTemplate to the 'gatingTemplate'.")
  expect_error(gate_edit(gs, 
                         parent = "T Cells", 
                         alias = "CD4 T Cells", 
                         gatingTemplate = "Test.csv"), 
               "Test.csv is not in this working directory.")
  
  gate_edit(gs3, 
            parent = "root", 
            alias = "Cells", 
            gatingTemplate = "Activation-gatingTemplate.csv")
  
  expect_equal(getGate(gs3, "Cells")[[1]], pg)
  expect_equal(gate_extract(parent = "root", 
                            alias = "Cells", 
                            "Activation-gatingTemplate.csv")[[1]][[1]][[1]], 
               pg)
  
  gate_edit(gs3, 
            parent = "root", 
            alias = "Cells", 
            overlay = "CD4 T Cells", 
            gatingTemplate = "Activation-gatingTemplate.csv", 
            type = "r")
  
  expect_equal(getGate(gs3, "Cells")[[1]], rg)
  expect_equal(gate_extract(parent = "root", 
                            alias = "Cells", 
                            "Activation-gatingTemplate.csv")[[1]][[1]][[1]], 
               rg)
  
  write.csv(gtf, "Activation-gatingTemplate.csv", row.names = FALSE)
  
})

# gate_type --------------------------------------------------------------------

test_that("gate_type", {
  
  expect_equal(gate_type(filters(list(rg))), "rectangle")
  expect_equal(gate_type(filters(list(pg))), "polygon")
  expect_equal(gate_type(filters(list(igx))), "interval")
  expect_equal(gate_type(filters(list(ig))), "interval")
  expect_equal(gate_type(filters(list(igy))), "interval")
  expect_equal(gate_type(filters(list(tg))), "threshold")
  expect_equal(gate_type(filters(list(tg1))), "threshold")
  expect_equal(gate_type(filters(list(bg))), "boundary")
  expect_equal(gate_type(filters(list(bg1))), "boundary")
  expect_equal(gate_type(filters(list(eg))), "ellipse")
  expect_equal(gate_type(qg), "quadrant")
  expect_equal(gate_type(wg), "web")
  expect_equal(gate_type(filters(list(ig,tg,bg))), 
               c("interval",
                 "threshold",
                 "boundary"))
  expect_equal(gate_type(filters(list(eg, eg, eg, eg))), 
               c("ellipse",
                 "ellipse",
                 "ellipse",
                 "ellipse"))
  expect_equal(gate_type(filters(list(igx, igx, igx, igx))), 
               c("interval",
                 "interval",
                 "interval",
                 "interval"))
  expect_equal(gate_type(filters(list(pg, pg, pg, pg))), 
               c("polygon",
                 "polygon",
                 "polygon",
                 "polygon"))
  expect_equal(gate_type(filters(list(rg,pg,ig,tg,bg,eg))), 
               c("rectangle",
                 "polygon",
                 "interval",
                 "threshold",
                 "boundary",
                 "ellipse"))
  expect_equal(gate_type(filters(list(rg,ig, tg,bg))), 
               c("rectangle",
                 "interval",
                 "threshold",
                 "boundary"))
  
})
