context("gate_draw")

# flowFrame method -------------------------------------------------------------

test_that("gate_draw flowFrame method", {
  
  gts <- gate_draw(fs[[1]],
                   alias = c("Cells","Cells","Cells","Cells","Cells","Cells"), 
                   channels = c("FSC-A","SSC-A"), 
                   type = c("r","p","e","i","t","b"), 
                   display = 0.1)
  
  expect_s4_class(gts, "filters")
  expect_equal(lapply(gts,function(x) {class(x)[1]}), 
               list(rectangle = "rectangleGate", 
                    polygon = "polygonGate", 
                    ellipse = "ellipsoidGate", 
                    interval = "rectangleGate", 
                    threshold = "rectangleGate", 
                    boundary = "rectangleGate"))
  expect_length(gts, 6)
  expect_equal(unname(gts), filters(list(rg,pg,eg,ig,tg,bg)), 
               tolerance = 0.01)
  
  gts <- gate_draw(fs[[1]], 
                   alias = c("A","B","C","D"), 
                   channels = c("FSC-A","SSC-A"), 
                   type = "q", 
                   display = 0.1)
  
  expect_s4_class(gts, "filters")
  expect_equal(lapply(gts,function(x) {class(x)[1]}), 
               list("rectangleGate",
                    "rectangleGate",
                    "rectangleGate",
                    "rectangleGate"))
  expect_length(gts, 4)
  expect_equal(gts, qg)
  
  gts <- gate_draw(fs[[1]], 
                   alias = c("A","B","C","D","E","F","G","H"), 
                   channels = c("FSC-A","SSC-A"), 
                   type = "w", 
                   display = 0.1)
  
  expect_s4_class(gts, "filters")
  expect_equal(unlist(lapply(gts,function(x) {class(x)[1]})), 
               rep("polygonGate",8))
  expect_length(gts, 8)
  expect_equal(gts, wg, tolerance = 0.01)
  
  gts <- gate_draw(fs[[1]], alias = "Cells", channels = "FSC-A")
  
  expect_equal(unname(gts), filters(list(igx)))
  
})

# flowSet method ---------------------------------------------------------------

test_that("gate_draw flowSet method", {
    
    gts <- gate_draw(fs, 
                     select = 1, 
                     alias = c("Cells",
                               "Cells",
                               "Cells",
                               "Cells",
                               "Cells",
                               "Cells"), 
                     channels = c("FSC-A","SSC-A"), 
                     type = c("r","p","e","i","t","b"), 
                     display = 0.1)
    
    expect_s4_class(gts, "filters")
    expect_equal(lapply(gts,function(x) {class(x)[1]}), 
                 list(rectangle = "rectangleGate", 
                      polygon = "polygonGate", 
                      ellipse = "ellipsoidGate", 
                      interval = "rectangleGate", 
                      threshold = "rectangleGate", 
                      boundary = "rectangleGate"))
    expect_length(gts, 6)
    expect_equal(unname(gts), filters(list(rg,pg,eg,ig,tg,bg)), 
                 tolerance = 0.01)
    
    expect_error(gate_draw(fs, 
                           select = "A", 
                           alias = "Cells", 
                           channels = c("FSC-A","SSC-A")), 
        "'select' must contain the numeric indicies of the samples to plot.")
    
    gts <- gate_draw(fs, 
                     alias = c("A","B","C","D"), 
                     channels = c("FSC-A","SSC-A"), 
                     type = "q", 
                     display = 0.1)
    
    expect_s4_class(gts, "filters")
    expect_equal(lapply(gts,function(x) {class(x)[1]}), 
                 list("rectangleGate",
                      "rectangleGate",
                      "rectangleGate",
                      "rectangleGate"))
    expect_length(gts, 4)
    expect_equal(gts, qg)
    
    gts <- gate_draw(fs, 
                     alias = c("A","B","C","D","E","F","G","H"), 
                     channels = c("FSC-A","SSC-A"), 
                     type = "w", 
                     display = 0.1)
    
    expect_s4_class(gts, "filters")
    expect_equal(unlist(lapply(gts,function(x) {class(x)[1]})), 
                 rep("polygonGate",8))
    expect_length(gts, 8)
    expect_equal(gts, wg, tolerance = 0.01)
    
    gts <- gate_draw(fs, alias = "Cells", channels = "FSC-A")
    
    expect_equal(unname(gts), filters(list(igx)))
    
  })

# GatingSet method -------------------------------------------------------------

test_that("gate_draw GatingSet method", {
  
  gs1 <- GatingSet(fs)
    
  gate_draw(gs1, 
            parent = "root", 
            alias = c("x","y","z"), 
            channels = c("FSC-A","SSC-A"), 
            type = c("e","p","i"), 
            display = 0.1, 
            gatingTemplate = "gatingTemplate.csv")
  
  expect_equal(basename(getNodes(gs1)), c("root","x","y","z"))
  
  gate_draw(gs1, 
            parent = "root", 
            alias = c("K","L","M"), 
            channels = c("FSC-A","SSC-A"),
            type = "w", 
            display = 0.1)
  
  expect_equal(basename(getNodes(gs1)), c("root","x","y","z","K","L","M"))
  
  gby <- gate_draw(gs1, 
                   group_by = "OVAConc", 
                   parent = "root", 
                   alias = c("X"), 
                   channels = c("FSC-A","SSC-A"), 
                   type = "p", 
                   display = 0.1)

  expect_equal(basename(getNodes(gs1)), c("root","x","y","z","K","L","M","X"))
  
})
unlink("gatingTemplate.csv")