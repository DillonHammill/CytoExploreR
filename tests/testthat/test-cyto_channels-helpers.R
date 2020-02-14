context("cyto_channels-helpers")

# CYTO_CHANNELS ----------------------------------------------------------------

test_that("cyto_channels", {
  
  # ALL CHANNELS
  expect_equal(cyto_channels(fs[[1]]),
               colnames(fs[[1]]))
  # SELECT CHANNELS
  expect_equal(cyto_channels(fs[[1]], select = "Alexa"),
               c("Alexa Fluor 488-A",
                 "Alexa Fluor 405-A",
                 "Alexa Fluor 430-A",
                 "Alexa Fluor 647-A",
                 "Alexa Fluor 700-A"))
  # EXCLUDE CHANNELS
  expect_equal(cyto_channels(fs[[1]], exclude = c("FSC","SSC")),
               colnames(fs[[1]])[-seq_len(6)])
})

# CYTO_MARKERS -----------------------------------------------------------------

test_that("cyto_markers", {
  
  # GatingSet -> flowFrame
  markers <- c("CD8",
               "Va2",
               "CD69",
               "Hoechst-405",
               "Hoechst-430",
               "CD44",
               "CD4",
               "CD11c")
  names(markers) <- c("Alexa Fluor 488-A",
                      "PE-A",
                      "7-AAD-A",
                      "Alexa Fluor 405-A",
                      "Alexa Fluor 430-A",
                      "Alexa Fluor 647-A",
                      "Alexa Fluor 700-A",
                      "APC-Cy7-A")
  expect_equal(cyto_markers(gs),
               markers)
  
  # GatingHierachy -> flowFrame
  expect_equal(cyto_markers(gs[[1]]),
               markers)
  
  # GatingSet - > flowFrame
  expect_equal(cyto_markers(gs[1]),
               markers)
  
  # flowSet
  expect_equal(cyto_markers(fs),
               markers)
  
})

# CYTO_FLUOR_CHANNELS ----------------------------------------------------------

test_that("cyto_fluor_channels",{
  # FLUORESCENT CHANNELS
  expect_equal(cyto_fluor_channels(fs[[1]]),
               cyto_channels(fs[[1]], exclude = c("Time",
                                                  "FSC",
                                                  "SSC")))
})

# CYTO_CHANNELS_EXTRACT --------------------------------------------------------

test_that("cyto_channels_extract", {
  # MARKERS TO CHANNELS
  expect_equal(cyto_channels_extract(gs, c("CD4","FSC-A")),
               c("Alexa Fluor 700-A", "FSC-A"))
  # INVALID NUMBER OF CHANNELS
  expect_error(cyto_channels_extract(gs, 
                                     channels = c("CD4","FSC-A", "Va2"),
                                       plot = TRUE),
               "Invalid number of supplied channels.")
  # INVALID MARKER
  expect_error(cyto_channels_extract(fs, c("B220","FSC-A")),
               "B220 is not a valid channel/marker.",
               fixed = TRUE)
})

# CYTO_MARKERS_EXTRACT ---------------------------------------------------------

test_that("cyto_markers_extract", {
  # CHANNELS TO MARKERS
  expect_equal(cyto_markers_extract(gs, c("Alexa Fluor 700-A","FSC-A")),
               c("CD4", "FSC-A"))
  # INVALID NUMBER OF CHANNELS
  expect_error(cyto_markers_extract(gs, 
                                     channels = c("CD4","FSC-A", "Va2"),
                                     plot = TRUE),
               "Invalid number of supplied channels.")
  # INVALID MARKER
  expect_error(cyto_markers_extract(gs, c("B220","FSC-A")),
               "'channels' contains invalid channel or marker names.",
               fixed = TRUE)
})

# CYTO_CHANNEL_SELECT ----------------------------------------------------------

test_that("cyto_channel_select", {
  
  mock_menu <- mock(4, 10, 11, 9, 1, 2, 12, cycle = TRUE)
  testthat::with_mock(menu = mock_menu,
                      expect_equal(cyto_channel_select(gs_comp),
                                   c("7-AAD-A",
                                     "Alexa Fluor 700-A",
                                     "APC-Cy7-A",
                                     "Alexa Fluor 647-A",
                                     "Alexa Fluor 488-A",
                                     "PE-A",
                                     "Unstained")))
  
})

# CYTO_CHANNEL_RESTRICT --------------------------------------------------------

test_that("cyto_channels_restrict", {
  
  fs_restricted <- cyto_channels_restrict(fs)
  expect_equal(cyto_channels(fs_restricted), 
               c("FSC-A",
                 "FSC-H",
                 "FSC-W",
                 "SSC-A",
                 "SSC-H",
                 "SSC-W",
                 "Alexa Fluor 488-A",
                 "PE-A",
                 "7-AAD-A",
                 "Alexa Fluor 405-A",
                 "Alexa Fluor 430-A",
                 "Alexa Fluor 647-A",
                 "Alexa Fluor 700-A",
                 "APC-Cy7-A",
                 "Time"))
  
})
