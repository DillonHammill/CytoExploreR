## CYTO_MATCH ------------------------------------------------------------------

test_that("cyto_match - NULL returns all indices", {
  expect_equal(cyto_match(gs, NULL), seq_along(gs))
})

test_that("cyto_match - numeric indices pass through", {
  expect_equal(cyto_match(gs, 1:5), 1:5)
})

test_that("cyto_match - exact rowname match", {
  nm <- rownames(cyto_details(gs))[1]
  expect_equal(cyto_match(gs, nm), 1L)
})

test_that("cyto_match - partial name match", {
  expect_equal(cyto_match(gs, "Activation_1"), 1L)
})

test_that("cyto_match - variable match", {
  pd       <- cyto_details(gs)
  expected <- which(pd[, "Treatment"] == "Stim-A")
  expect_equal(cyto_match(gs, Treatment = "Stim-A"), expected)
})

test_that("cyto_match - exclude returns negative indices", {
  pd       <- cyto_details(gs)
  expected <- -which(pd[, "Treatment"] == "Stim-A")
  expect_equal(cyto_match(gs, Treatment = "Stim-A", exclude = TRUE), expected)
})

test_that("cyto_match - invalid variable errors", {
  expect_error(cyto_match(gs, FakeVar = "foo"))
})

test_that("cyto_match - invalid level errors", {
  expect_error(cyto_match(gs, Treatment = "NonExistent"))
})

## CYTO_NODES ------------------------------------------------------------------

test_that("cyto_nodes - returns all 11 nodes including root", {
  nodes <- cyto_nodes(gs)
  expect_equal(length(nodes), 11L)
  expect_true("root" %in% nodes)
})

test_that("cyto_nodes - terminal returns only leaf nodes", {
  nodes <- cyto_nodes(gs, terminal = TRUE)
  expect_equal(length(nodes), 4L)
  expect_false("root" %in% nodes)
  expect_setequal(nodes, c(
    "/Cells/Single Cells/Dead Cells",
    "/Cells/Single Cells/Live Cells/T Cells/CD8 T Cells/CD69+ CD8 T Cells",
    "/Cells/Single Cells/Live Cells/T Cells/CD4 T Cells/CD69+ CD4 T Cells",
    "/Cells/Single Cells/Live Cells/Dendritic Cells"
  ))
})

test_that("cyto_nodes - select filters to matching nodes", {
  nodes <- cyto_nodes(gs, select = "T Cells")
  expect_equal(length(nodes), 5L)
  expect_true(all(grepl("T Cells", nodes, fixed = TRUE)))
})

test_that("cyto_nodes - exclude removes matching nodes", {
  nodes <- cyto_nodes(gs, exclude = "root")
  expect_false("root" %in% nodes)
  expect_equal(length(nodes), 10L)
})

test_that("cyto_nodes - select by numeric index", {
  all_nodes <- cyto_nodes(gs)
  expect_equal(cyto_nodes(gs, select = 1L), all_nodes[1])
})

## CYTO_GROUPS -----------------------------------------------------------------

test_that("cyto_groups - group_by all returns Combined Events", {
  groups <- cyto_groups(gs, group_by = "all")
  expect_equal(unname(groups), "Combined Events")
  expect_equal(as.integer(names(groups)), length(gs))
})

test_that("cyto_groups - group_by variable returns one group per level", {
  groups <- cyto_groups(gs, group_by = "Treatment")
  expect_setequal(unname(groups), c("NA", "Stim-A", "Stim-B", "Stim-C", "Stim-D"))
  expect_equal(sum(as.integer(names(groups))), length(gs))
})

test_that("cyto_groups - list controls level order and appends missing levels", {
  levels_specified <- c("Stim-A", "Stim-B", "Stim-C", "Stim-D")
  groups <- cyto_groups(gs, group_by = list(Treatment = levels_specified))
  expect_equal(unname(groups[seq_along(levels_specified)]), levels_specified)
  expect_equal(unname(groups[length(groups)]), "NA")
})

test_that("cyto_groups - details returns split pData as list", {
  pd_split <- cyto_groups(gs, group_by = "Treatment", details = TRUE)
  expect_true(is.list(pd_split))
  expect_setequal(names(pd_split), c("NA", "Stim-A", "Stim-B", "Stim-C", "Stim-D"))
  expect_true(all(sapply(pd_split, is.data.frame)))
})

test_that("cyto_groups - invalid variable errors", {
  expect_error(cyto_groups(gs, group_by = "FakeVar"))
})

test_that("cyto_groups - invalid factor level errors", {
  expect_error(cyto_groups(gs, group_by = list(Treatment = c("Stim-A", "NonExistent"))))
})

## CYTO_APPLY ------------------------------------------------------------------

test_that("cyto_apply - result names match cyto_names", {
  res <- cyto_apply(cs, nrow, simplify = FALSE)
  expect_equal(names(res), cyto_names(cs))
})

test_that("cyto_apply - cytoframe input applies FUN to each frame", {
  res <- cyto_apply(cs, cyto_stat_count)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), length(cs))
  expect_equal(rownames(res), cyto_names(cs))
})

test_that("cyto_apply - matrix input applies FUN to expression matrix", {
  res <- cyto_apply(cs, nrow, input = "matrix", simplify = FALSE)
  expect_equal(names(res), cyto_names(cs))
  expect_true(all(sapply(res, is.numeric)))
})

test_that("cyto_apply - numeric input aliases resolve correctly", {
  res_str <- cyto_apply(cs, cyto_stat_count, input = "cytoframe")
  res_num <- cyto_apply(cs, cyto_stat_count, input = 2)
  expect_identical(res_str, res_num)
})
