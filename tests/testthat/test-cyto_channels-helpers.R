# CYTO_CHANNELS ----------------------------------------------------------------

test_that("cyto_channels", {
  # GATINGSET
  expect_snapshot_output(cyto_channels(gs))
  # GATINGHIERARCHY
  expect_snapshot_output(cyto_channels(gs[[1]]))
  # CYTOSET
  expect_snapshot_output(cyto_channels(cs))
  # CYTOFRAME
  expect_snapshot_output(cyto_channels(cs[[1]]))
})

# CYTO_MARKERS -----------------------------------------------------------------

test_that("cyto_markers", {
  # GATINGSET
  expect_snapshot_output(cyto_markers(gs))
  # GATINGHIERARCHY
  expect_snapshot_output(cyto_markers(gs[[1]]))
  # CYTOSET
  expect_snapshot_output(cyto_markers(cs))
  # CYTOFRAME
  expect_snapshot_output(cyto_markers(cs[[1]]))
})

# CYTO_FLUOR_CHANNELS ----------------------------------------------------------

test_that("cyto_fluor_channels", {
  # GATINGSET
  expect_snapshot_output(cyto_fluor_channels(gs))
})

# CYTO_CHANNELS_EXTRACT --------------------------------------------------------

test_that("cyto_channels_extract", {
  # GATINGSET
  expect_snapshot_output(cyto_channels_extract(gs, 
                                               c("cd6",     # part marker
                                                 "CD4",     # marker
                                                 "FSC",     # part channel
                                                 "SSC-A"))) # channel
})

# CYTO_MARKERS_EXTRACT ---------------------------------------------------------

test_that("cyto_markers_extract", {
  # GATINGSET
  expect_snapshot_output(cyto_markers_extract(gs,
                                              c("cd6",     # part marker
                                                "CD4",     # marker
                                                "FSC",     # part channel
                                                "SSC-A"))) # channel
})

# CYTO_CHANNEL_SELECT ----------------------------------------------------------

test_that("cyto_channel_select", {
  # MOCK DATA_EDIT
  stub(cyto_channel_select,
       "data_edit",
         data.frame("name" = cyto_names(gs_comp),
                    "channel" = c("7-AAD-A",
                                  "Alexa Fluor 700-A",
                                  "APC-Cy7-A",
                                  "Alexa Fluor 647-A",
                                  "Alexa Fluor 488-A",
                                  "PE-A",
                                  "Unstained"),
                    stringsAsFactors = FALSE))
  expect_equal(
    cyto_channel_select(gs_comp),
    c("7-AAD-A",
      "Alexa Fluor 700-A",
      "APC-Cy7-A",
      "Alexa Fluor 647-A",
      "Alexa Fluor 488-A",
      "PE-A",
      "Unstained")
  )
})

# CYTO_CHANNEL_MATCH -----------------------------------------------------------

test_that("cyto_channel_match", {
  stub(cyto_channel_match,
       "cyto_channel_select",
       c("7-AAD-A",
         "Alexa Fluor 700-A",
         "APC-Cy7-A",
         "Alexa Fluor 647-A",
         "Alexa Fluor 488-A",
         "PE-A",
         "Unstained"))
  expect_equal(cyto_channel_match(gs_comp,
                                  save_as = paste0(temp_dir, 
                                                   "Channel_Match.csv")),
               data.frame("name" = c("Compensation-7AAD.fcs",
                                     "Compensation-AF700.fcs",    
                                     "Compensation-APC-Cy7.fcs",
                                     "Compensation-APC.fcs",      
                                     "Compensation-FITC.fcs",
                                     "Compensation-PE.fcs",       
                                     "Compensation-Unstained.fcs"),
                          "channel" = c("7-AAD-A",
                                        "Alexa Fluor 700-A",
                                        "APC-Cy7-A",
                                        "Alexa Fluor 647-A",
                                        "Alexa Fluor 488-A",
                                        "PE-A",
                                        "Unstained"),
                          stringsAsFactors = FALSE))
  expect_true(file.exists(paste0(temp_dir, "Channel_Match.csv")))
  unlink(paste0(temp_dir, "Channel_Match.csv"))
})
