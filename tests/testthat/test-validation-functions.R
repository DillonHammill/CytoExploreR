context("Validation Functions")

## file_wd_check ---------------------------------------------------------------

test_that("file_wd_check", {
  
  expect_false(file_wd_check("Test.csv"))
  
})

## cyto_channels_extract -------------------------------------------------------

test_that("cyto_channels_extract flowFrame method", {

  # flowFrame
  expect_equal(cyto_channels_extract(fs[[1]], 
                                     c("CD4","CD8","PE-A")),
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "PE-A"))
  
  expect_error(cyto_channels_extract(fs[[1]], 
                                     c("CD4","CD8","PE-A"),
                                     plot = TRUE),
               "Invalid number of supplied channels.")
  
  expect_error(cyto_channels_extract(fs[[1]], 
                                     c("CD4","CD8","P-A")),
               "P-A is not a valid channel/marker", fixed = TRUE)
  
  # flowSet
  expect_equal(cyto_channels_extract(fs, 
                                     c("CD4","CD8","PE-A")),
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "PE-A"))
  
  expect_error(cyto_channels_extract(fs, 
                                     c("CD4","CD8","PE-A"),
                                     plot = TRUE),
               "Invalid number of supplied channels.")
  
  # GatingHierarchy
  expect_equal(cyto_channels_extract(gs[[1]], 
                                     c("CD4","CD8","PE-A")),
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "PE-A"))
  
  expect_error(cyto_channels_extract(gs[[1]], 
                                     c("CD4","CD8","PE-A"),
                                     plot = TRUE),
               "Invalid number of supplied channels.")
  
  # GatingSet
  expect_equal(cyto_channels_extract(gs, 
                                     c("CD4","CD8","PE-A")),
               c("Alexa Fluor 700-A", "Alexa Fluor 488-A", "PE-A"))
  
  expect_error(cyto_channels_extract(gs, 
                                     c("CD4","CD8","PE-A"),
                                     plot = TRUE),
               "Invalid number of supplied channels.")
  
})

# cyto_markers_extract ---------------------------------------------------------

test_that("cyto_markers_extract", {
  
  #flowFrame
  expect_equal(cyto_markers_extract(fs[[1]],
                                    c("CD4","PE-A")),
               c("CD4","Va2"))
  
  expect_error(cyto_markers_extract(fs[[1]], 
                                    c("CD4","CD8","PE-A"),
                                    plot = TRUE),
               "Invalid number of supplied channels.")
  
  expect_warning(cyto_markers_extract(fs[[1]], 
                                      c("CD","P-A"),
                                      plot = TRUE),
               "'channels' contains invalid channel or marker names.")
  
  expect_equal(cyto_markers_extract(fs[[1]], 
                                    c("CD4","CD8","PE-A","FSC-A")),
               c("CD4","CD8","Va2","FSC-A"))
  
  # flowSet
  expect_equal(cyto_markers_extract(fs,
                                    c("CD4","PE-A")),
               c("CD4","Va2"))
  
  expect_error(cyto_markers_extract(fs, 
                                    c("CD4","CD8","PE-A"),
                                    plot = TRUE),
               "Invalid number of supplied channels.")
  
  # GatingHierarchy
  expect_equal(cyto_markers_extract(gs[[1]],
                                    c("CD4","PE-A")),
               c("CD4","Va2"))
  
  expect_error(cyto_markers_extract(gs[[1]], 
                                    c("CD4","CD8","PE-A"),
                                    plot = TRUE),
               "Invalid number of supplied channels.")
  
  # GatingSet
  expect_equal(cyto_markers_extract(gs,
                                    c("CD4","PE-A")),
               c("CD4","Va2"))
  
  expect_error(cyto_markers_extract(gs, 
                                    c("CD4","CD8","PE-A"),
                                    plot = TRUE),
               "Invalid number of supplied channels.")
  
})