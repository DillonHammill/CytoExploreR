context("cyto_channels-helpers")

# CYTO_CHANNELS ----------------------------------------------------------------

test_that("cyto_channels", {
  # ALL CHANNELS
  expect_equal(cyto_channels(fs[[1]]),
               chans)
  # EXCLUDE CHANNELS
  expect_equal(cyto_channels(fs[[1]], exclude = c("FSC","SSC")),
               chans[-seq_len(6)])
})

# CYTO_FLUOR_CHANNELS ----------------------------------------------------------

test_that("cyto_fluor_channels",{
  # FLUORESCENT CHANNELS
  expect_equal(cyto_fluor_channels(fs[[1]]),
               chans[-c(seq_len(6),length(chans))])
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
                      expect_equal(cyto_channel_select(Comp),
                                   c("7-AAD-A",
                                     "Alexa Fluor 700-A",
                                     "APC-Cy7-A",
                                     "Alexa Fluor 647-A",
                                     "Alexa Fluor 488-A",
                                     "PE-A",
                                     "Unstained")))
  
})
