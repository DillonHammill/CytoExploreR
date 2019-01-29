context("Helper Functions")

## ------------------------------------------------------------------
# cyto_fluor_channels -
test_that("cyto_fluor_channels returns the correct channels for flowFrames, flowSets and GatingSets", {
  expect_equal(cyto_fluor_channels(fs[[1]]), c("Alexa Fluor 488-A", "PE-A", "PE-Texas Red-A", "7-AAD-A", "PE-Cy7-A", "Alexa Fluor 405-A", "Alexa Fluor 430-A", "Qdot 605-A", "Alexa Fluor 647-A", "Alexa Fluor 700-A", "APC-Cy7-A"))
  expect_equal(cyto_fluor_channels(fs), c("Alexa Fluor 488-A", "PE-A", "PE-Texas Red-A", "7-AAD-A", "PE-Cy7-A", "Alexa Fluor 405-A", "Alexa Fluor 430-A", "Qdot 605-A", "Alexa Fluor 647-A", "Alexa Fluor 700-A", "APC-Cy7-A"))
  expect_equal(cyto_fluor_channels(gs), c("Alexa Fluor 488-A", "PE-A", "PE-Texas Red-A", "7-AAD-A", "PE-Cy7-A", "Alexa Fluor 405-A", "Alexa Fluor 430-A", "Qdot 605-A", "Alexa Fluor 647-A", "Alexa Fluor 700-A", "APC-Cy7-A"))
})

## ------------------------------------------------------------------
# cyto_channel_select -

test_that("cyto_channel_select", {
  print(pData(Comp))
  expect_equal(cyto_channel_select(Comp[[1]]), "PE-Cy7-A")
  expect_equal(cyto_channel_select(Comp), c("7-AAD-A", "Alexa Fluor 430-A", "APC-Cy7-A", "Unstained", "PE-Cy7-A", "PE-A"))
  expect_equal(cyto_channel_select(gsc), c("7-AAD-A", "Alexa Fluor 430-A", "APC-Cy7-A", "Unstained", "PE-Cy7-A", "PE-A"))
})

## ------------------------------------------------------------------
# cyto_sample -
test_that("cyto_sample returns subsetted flowFrame", {
  expect_equal(nrow(exprs(cyto_sample(fs[[1]], 0.5))), 1000)
})
