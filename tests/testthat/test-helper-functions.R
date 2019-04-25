context("Helper Functions")

# cyto_marker_extract ----------------------------------------------------------

test_that("cyto_marker_extract", {
  
  expect_equal(cyto_marker_extract(Va2[[4]], c("FSC-A","SSC-A")), 
               c("FSC-A","SSC-A"))
  expect_equal(cyto_marker_extract(Va2[[4]], c("FSC-A","CD4")), 
               c("FSC-A","CD4"))

})

# cyto_channel_extract ---------------------------------------------------------

test_that("cyto_channel_extract", {
  
  expect_error(cyto_channel_extract(Va2[[4]], c("FSC-A","SSC-A")),
               "Supplied markers are not valid.")
  expect_equal(cyto_channel_extract(Va2[[4]], c("CD8","CD4")), 
               c("Alexa Fluor 488-A","Alexa Fluor 700-A"))
  
})

# cyto_fluor_channels ----------------------------------------------------------

test_that("cyto_fluor_channels", {
  
  # flowFrame ------------------------------------------------------------------
  expect_equal(cyto_fluor_channels(fs[[1]]), 
               c("Alexa Fluor 488-A", 
                 "PE-A", "PE-Texas Red-A", 
                 "7-AAD-A", "PE-Cy7-A", 
                 "Alexa Fluor 405-A", 
                 "Alexa Fluor 430-A", 
                 "Qdot 605-A", 
                 "Alexa Fluor 647-A", 
                 "Alexa Fluor 700-A", 
                 "APC-Cy7-A"))
  
  # flowSet --------------------------------------------------------------------
  expect_equal(cyto_fluor_channels(fs), 
               c("Alexa Fluor 488-A", 
                 "PE-A", 
                 "PE-Texas Red-A",
                 "7-AAD-A", "PE-Cy7-A", 
                 "Alexa Fluor 405-A", 
                 "Alexa Fluor 430-A", 
                 "Qdot 605-A", 
                 "Alexa Fluor 647-A", 
                 "Alexa Fluor 700-A", 
                 "APC-Cy7-A"))
  
  # GatingHierarchy ------------------------------------------------------------
  expect_equal(cyto_fluor_channels(gs[[1]]), 
               c("Alexa Fluor 488-A", 
                 "PE-A", 
                 "PE-Texas Red-A",
                 "7-AAD-A", "PE-Cy7-A", 
                 "Alexa Fluor 405-A", 
                 "Alexa Fluor 430-A", 
                 "Qdot 605-A", 
                 "Alexa Fluor 647-A", 
                 "Alexa Fluor 700-A", 
                 "APC-Cy7-A"))
  
  # GatingSet ------------------------------------------------------------------
  expect_equal(cyto_fluor_channels(gs), 
               c("Alexa Fluor 488-A", 
                 "PE-A", 
                 "PE-Texas Red-A", 
                 "7-AAD-A", 
                 "PE-Cy7-A", 
                 "Alexa Fluor 405-A",
                 "Alexa Fluor 430-A", 
                 "Qdot 605-A",
                 "Alexa Fluor 647-A", 
                 "Alexa Fluor 700-A", 
                 "APC-Cy7-A"))
  
})

# cyto_channel_select ----------------------------------------------------------

test_that("cyto_channel_select", {
  
  # flowFrame ------------------------------------------------------------------
  mock_menu <- mock(5)
  testthat::with_mock(menu = mock_menu, 
            expect_equal(cyto_channel_select(Comp[[1]]), 
                         "PE-Cy7-A"))
  
  # flowSet --------------------------------------------------------------------
  mock_menu <- mock(4,7,11,12,5,2, cycle = TRUE)
  testthat::with_mock(menu = mock_menu,
                      expect_equal(cyto_channel_select(Comp),
                                   c("7-AAD-A",
                                     "Alexa Fluor 430-A",
                                     "APC-Cy7-A",
                                     "Unstained",
                                     "PE-Cy7-A",
                                     "PE-A")))
  
  # GatingHierarchy ------------------------------------------------------------
  mock_menu <- mock(5)
  testthat::with_mock(menu = mock_menu, 
                      expect_equal(cyto_channel_select(gsc[[1]]), 
                                   "PE-Cy7-A"))
  
  # GatingSet ------------------------------------------------------------------
  mock_menu <- mock(4,7,11,12,5,2, cycle = TRUE)
  testthat::with_mock(menu = mock_menu,
                      expect_equal(cyto_channel_select(gsc),
                                   c("7-AAD-A",
                                     "Alexa Fluor 430-A",
                                     "APC-Cy7-A",
                                     "Unstained",
                                     "PE-Cy7-A",
                                     "PE-A")))
  
})

# cyto_sample ------------------------------------------------------------------

test_that("cyto_sample returns subsetted flowFrame", {
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 1))), 2000)
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 0.5))), 1000)
})
