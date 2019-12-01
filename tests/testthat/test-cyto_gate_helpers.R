context("cyto_gate_helpers")

# CLONE GATINGSET
gs_clone <- gs_clone(gs)

# DUPLICATE GATINGTEMPLATE
file.copy("Reference-Activation-gatingTemplate.csv", 
          "Activation-gatingTemplate.csv")

# CYTO_GATE_EDIT ---------------------------------------------------------------

test_that("cyto_gate_edit", {
  
  expect_error(cyto_gate_edit(gs), 
               "Supply the name of the parent population.")
  expect_error(cyto_gate_edit(gs, 
                              parent = "Test"), 
               "Supplied parent does not exist in the GatingSet.")
  expect_error(cyto_gate_edit(gs, 
                              parent = "T Cells"), 
               "Supply the name(s) of the gates to edit to 'alias'.", 
               fixed =  TRUE)
  expect_error(cyto_gate_edit(gs, 
                              parent = "T Cells", 
                              alias = "Test"), 
               "Supplied alias does not exist in the GatingSet.")
  expect_error(cyto_gate_edit(gs, 
                              parent = "T Cells", 
                              alias = "CD4 T Cells"), 
               "Supply the name of the gatingTemplate to edit gate(s).",
               fixed = TRUE)
  
  # EDIT CD4 T Cells & CD8 T Cells Gates
  mock_cyto_gate_draw <- mock(list(filters(list(CD4)),
                                   filters(list(CD8))))
  testthat::with_mock(cyto_gate_draw = mock_cyto_gate_draw,
                      cyto_gate_edit(gs_clone,
                                     parent = "T Cells",
                                     alias = c("CD4 T Cells",
                                               "CD8 T Cells"),
                                     channels = c("CD4","CD8"),
                                     gatingTemplate = "Activation-gatingTemplate.csv"))
  expect_equal(gh_pop_get_gate(gs_clone[[1]], "CD4 T Cells"),
               CD4)
  expect_equal(gh_pop_get_gate(gs_clone[[1]], "CD8 T Cells"),
               CD8)
  expect_equivalent(cyto_gate_extract(parent = "T Cells",
                                 alias = c("CD4 T Cells",
                                           "CD8 T Cells"),
                                 gatingTemplate = "Activation-gatingTemplate.csv"),
                    list(filters(list(CD4)),
                         filters(list(CD8))))
  
  # CLOSE GRAPHICS DEVICE
  dev.off()
  
})

# CYTO_GATE_REMOVE -------------------------------------------------------------

test_that("cyto_gate_remove", {
  
  expect_error(cyto_gate_remove(gs_clone,
                                parent = "T Cells"), 
               "Supply the name of the population to be removed to 'alias'.")
  expect_error(cyto_gate_remove(gs_clone, 
                                alias = "Test"), 
               "Supplied alias does not exist in the GatingSet.")
  
  expect_error(cyto_gate_remove(gs_clone, alias = "T Cells"), 
               "Supply the name of the gatingTemplate to remove gate(s).",
               fixed  =TRUE)

  cyto_gate_remove(gs_clone, 
                   parent = "T Cells",
                   alias = "CD4 T Cells", 
                   gatingTemplate = "Activation-gatingTemplate.csv")
  expect_equal(cyto_nodes(gs_clone, path = "auto"), 
               c("root",
                 "Cells",
                 "Single Cells",
                 "Dead Cells",
                 "Live Cells",
                 "T Cells",
                 "CD8 T Cells",
                 "CD69+ CD8 T Cells",
                 "Dendritic Cells"))
  expect_equal(basename(gatingTemplate("Activation-gatingTemplate.csv")@nodes),
               c("root",
                 "Cells",
                 "Single Cells",
                 "Dead Cells",
                 "Live Cells",
                 "Dendritic Cells",
                 "T Cells",
                 "CD8 T Cells",
                 "CD69+ CD8 T Cells"))
  
  cyto_gate_remove(gs_clone, 
                   alias = c("Dendritic Cells", "T Cells"), 
                   gatingTemplate = "Activation-gatingTemplate.csv")
  expect_equal(cyto_nodes(gs_clone, path = "auto"), 
               c("root",
                 "Cells",
                 "Single Cells",
                 "Dead Cells",
                 "Live Cells"))
  expect_equal(basename(gatingTemplate("Activation-gatingTemplate.csv")@nodes), 
               c("root",
                 "Cells",
                 "Single Cells",
                 "Dead Cells",
                 "Live Cells"))
  
})

# CYTO_GATE_RENAME -------------------------------------------------------------

test_that("cyto_gate_rename", {
  
  expect_error(cyto_gate_rename(gs_clone,
                                alias = "T Cells"),
               "Supply the name of the gatingTemplate to rename gate(s).",
               fixed = TRUE)
  expect_error(cyto_gate_rename(gs_clone,
                                alias = "T Lymphocytes",
                                gatingTemplate = "Activation-gatingTemplate.csv"),
               "Supplied gate(s) do not exist in this GatingSet.",
               fixed = TRUE)
  
  cyto_gate_rename(gs_clone,
                   alias = c("Single Cells",
                             "Live Cells",
                             "Dead Cells"),
                   names = c("Singlets",
                             "Live",
                             "Dead"),
                   gatingTemplate = "Activation-gatingTemplate.csv")
  expect_equal(cyto_nodes(gs_clone, path = "auto"),
               c("root",
                 "Cells",
                 "Singlets",
                 "Dead",
                 "Live"))
  expect_equal(basename(gatingTemplate("Activation-gatingTemplate.csv")@nodes), 
               c("root",
                 "Cells",
                 "Singlets",
                 "Dead",
                 "Live"))
  
})

# CYTO_GATE_EXTRACT ------------------------------------------------------------

test_that("cyto_gate_extract", {
  
  expect_error(cyto_gate_extract(alias = "T Cells"),
               "Supply the name of the parent population.")
  expect_error(cyto_gate_extract(parent = "T Cells"), 
               "Supply the name(s) of the alias to extract.", 
               fixed = TRUE)
  expect_error(cyto_gate_extract(parent = "T Cells", 
                                 alias = "CD4 T Cells"), 
               "Supply the name of the gatingTemplate to extract gate(s).",
               fixed = TRUE)
  
  gate <- cyto_gate_extract(parent = "T Cells", 
                            alias = c("CD4 T Cells","CD8 T Cells"), 
                            gatingTemplate = "Reference-Activation-gatingTemplate.csv")
  expect_equivalent(gate, list(filters(list(CD4)), 
                               filters(list(CD8))), 
               tolerance = 0.01)
  
})

# CYTO_GATE_TYPE ---------------------------------------------------------------

test_that("cyto_gate_type", {
  
  expect_equal(cyto_gate_type(filters(list(rg2))), "rectangle")
  expect_equal(cyto_gate_type(filters(list(pg))), "polygon")
  expect_equal(cyto_gate_type(filters(list(igx))), "interval")
  expect_equal(cyto_gate_type(filters(list(ig))), "interval")
  expect_equal(cyto_gate_type(filters(list(igy))), "interval")
  expect_equal(cyto_gate_type(filters(list(tg))), "threshold")
  expect_equal(cyto_gate_type(filters(list(tg1))), "threshold")
  expect_equal(cyto_gate_type(filters(list(bg))), "boundary")
  expect_equal(cyto_gate_type(filters(list(bg1))), "boundary")
  expect_equal(cyto_gate_type(filters(list(eg))), "ellipse")
  expect_equal(cyto_gate_type(qg), "quadrant")
  expect_equal(cyto_gate_type(wg), "web")
  expect_equal(cyto_gate_type(filters(list(ig,tg,bg))), 
               c("interval",
                 "threshold",
                 "boundary"))
  expect_equal(cyto_gate_type(filters(list(eg, eg, eg, eg))), 
               c("ellipse",
                 "ellipse",
                 "ellipse",
                 "ellipse"))
  expect_equal(cyto_gate_type(filters(list(igx, igx, igx, igx))), 
               c("interval",
                 "interval",
                 "interval",
                 "interval"))
  expect_equal(cyto_gate_type(filters(list(pg, pg, pg, pg))), 
               c("polygon",
                 "polygon",
                 "polygon",
                 "polygon"))
  expect_equal(cyto_gate_type(filters(list(rg2,pg,ig,tg,bg,eg))), 
               c("rectangle",
                 "polygon",
                 "interval",
                 "threshold",
                 "boundary",
                 "ellipse"))
  expect_equal(cyto_gate_type(filters(list(rg2,ig, tg,bg))), 
               c("rectangle",
                 "interval",
                 "threshold",
                 "boundary"))
  
})

# REMOVE GATINGTEMPLATE
base::unlink("Activation-gatingTemplate.csv")
