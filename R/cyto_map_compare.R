# CYTO_MAP_COMPARE -------------------------------------------------------------
# 
# cyto_map_compare <- function(x,
#                              select = NULL,
#                              parent = "root",
#                              merge_by = NA,
#                              channels = NULL,
#                              base_test = "ks",
#                              plot = TRUE,
#                              ...) {
#   
#   # SELECT
#   if(!is.null(select)) {
#     x <- cyto_select(
#       x,
#       select
#     )
#   }
#   
#   # CROSS ENTROPY
#   if(is.null(channels)) {
#     channels <- cyto_channels(
#       x,
#       select = "entropy"
#     )
#     # MISSING CROSS ENTROPY
#     if(length(channels) == 0) {
#       stop(
#         paste0(
#           "Cannot locate channel of cross-entropy values - either supply the ",
#           "channel manually to 'channels' or call cyto_map() with ",
#           "cross_entropy = TRUE."
#         )
#       )
#     # SELECT CROSS ENTROPY
#     } else if(length(channels) > 1) {
#       # REQUEST CROSS ENTROPY
#       if(interactive() & cyto_option("CytoExploreR_interactive")) {
#         message(
#           paste0(
#             "Multiple channels containing cross-entropy values located. ",
#             "Which channel do you want to use?\n",
#             paste0(
#               1:length(channels),
#               ". ",
#               channels,
#               sep = "\n"
#             )
#           )
#         )
#         # SELECTION
#         opt <- cyto_enquire(
#           NULL
#         )
#         if(is.numeric(opt)) {
#           channels <- channels[opt]
#         }
#       } else {
#         warning(
#           paste0(
#             "Multiple channels containing cross-entropy values located! ",
#             "Resorting to using ",
#             channels[1], "."
#           )
#         )
#         channels <- channels[1]
#       }
#     }
#   # PREPARE CHANNELS
#   } else {
#     channels <- cyto_channels_extract(
#       x,
#       channels = channels
#     )
#   }
#   
#   # SPLIT INTO GROUPS
#   x <- cyto_group_by(
#     x,
#     group_by = merge_by
#   )
#   
#   # SINGLE GROUP
#   if(length(x) == 1) {
#     stop(
#       "Cannot test for cross-entropy differences on a single group!"
#     )
#   }
#   
#   # COMPARISONS PER POPULATION
#   res <- structure(
#     lapply(
#       parent,
#       function(z) {
#         # CREATE MATRIX TO STORE P VALUES
#         pval <- matrix(
#           1,
#           ncol = length(x),
#           nrow = length(x),
#           dimnames = list(names(x), names(x))
#         )
#         # CREATE MATRIX TO STORE DISTANCES
#         dst <- pval
#         diag(dst) <- 0
#         # COMPARE CROSS ENTROPY VALUES - KS TEST
#         if(base_test == "ks") {
#           ks_tests <- list()
#           ks_pval <- c()
#           for(i in 1:(length(x) - 1)) {
#             for(j in (i + 1):length(x)) {
#               # COMPUTE KOLMOGOROV-SMIRNOV TEST
#               ks_test <- suppressWarnings(
#                 ks.test(
#                   cyto_data_extract(
#                     x[[i]],
#                     parent = z,
#                     channels = channels,
#                     format = "matrix",
#                     coerce = TRUE,
#                     copy = FALSE
#                   )[[1]][[1]][, 1],
#                   cyto_data_extract(
#                     x[[j]],
#                     parent = z,
#                     channels = channels,
#                     format = "matrix",
#                     coerce = TRUE,
#                     copy = FALSE
#                   )[[1]][[1]][, 1]
#                 )
#               )
#               # STORE P VALUES FOR CORRECTION
#               ks_pval <- c(ks_pval, ks_test$p.value)
#               # UPDATE SLOTS
#               ks_test$p.adj.value <- ks_test$p.value
#               # STORE CO-ORDINATE PAIRS
#               ks_test$coords <- c(i, j)
#               # STORE KOLMOGOROV-SMIRNOV TEST RESULTS
#               ks_tests <- c(
#                 ks_tests,
#                 structure(
#                   list(
#                     ks_test
#                   ),
#                   names = paste0(
#                     names(x)[c(i,j)],
#                     collapse = " x "
#                   )
#                 )
#                 
#               )
#             }
#           }
#           # CORRECT FOR MULTIPLE COMPARISONS
#           if(length(x) > 2) {
#             # HOLMS ADJUSTMENT
#             ks_adj_pval <- p.adjust(ks_pval, "holm")
#             # STORE ADJUSTED P VALUES
#             ks_tests <- structure(
#               lapply(
#                 seq_along(ks_tests),
#                 function(v) {
#                   ks_tests[[v]]$p.adj.value <- ks_adj_pval[v]
#                   return(ks_tests[[v]])
#                 }
#               ),
#               names = names(ks_tests)
#             )
#           }
#           # STORE P VALUES & DISTANCES
#           for(i in 1:length(ks_tests)) {
#             # P VALUES
#             pval[ks_tests[[i]]$coords[1], ks_tests[[i]]$coords[2]] <-
#               ks_tests[[i]]$p.adj.value
#             pval[ks_tests[[i]]$coords[2], ks_tests[[i]]$coords[1]] <-
#               ks_tests[[i]]$p.adj.value
#             # DISTANCES
#             dst[ks_tests[[i]]$coords[1], ks_tests[[i]]$coords[2]] <-
#               ks_tests[[i]]$statistic
#             dst[ks_tests[[i]]$coords[2], ks_tests[[i]]$coords[1]] <-
#               ks_tests[[i]]$statistic
#           }
#         # COMPARE CROSS ENTROPY VALUES  - RANK TEST
#         } else if(base_test == "rank") {
#           # PAIRWISE COMPARISON - WILCOX TEST
#           if(length(x) == 2) {
#            # WILCOX TEST
#             rank_test <- wilcox.test(
#               cyto_data_extract(
#                 x[[1]],
#                 parent = z,
#                 channels = channels,
#                 format = "matrix",
#                 coerce = TRUE,
#                 copy = FALSE
#               )[[1]][[1]][, 1],
#               cyto_data_extract(
#                 x[[2]],
#                 parent = z,
#                 channels = channels,
#                 format = "matrix",
#                 coerce = TRUE,
#                 copy = FALSE
#               )[[1]][[1]][, 1]
#             )
#             # STORE P VALUES
#             pval[1, 2] <- rank_test$p.value
#             pval[2, 1] <- rank_test$p.value
#             # COMPUTE DISTANCE
#             d <- abs(
#               median(
#                 cyto_data_extract(
#                   x[[1]],
#                   parent = z,
#                   channels = channels,
#                   format = "matrix",
#                   coerce = TRUE,
#                   copy = FALSE
#                 )[[1]][[1]][, 1]
#               ) - 
#               median(
#                 cyto_data_extract(
#                   x[[2]],
#                   parent = z,
#                   channels = channels,
#                   format = "matrix",
#                   coerce = TRUE,
#                   copy = FALSE
#                 )[[1]][[1]][, 1]
#               )
#             )
#             # STORE DISTANCES
#             dst[1, 2] <- d
#             dst[2, 1] <- d
#           # MULTIPLE COMPARISONS - KRUSKAL-WALLIS + DUNN CORRECTION
#           } else {
#             # DUNN.TEST PACKAGE REQUIRED
#             cyto_require(
#               "dunn.test",
#               source = "CRAN"
#             )
#             # DUNN MULTIPLE COMPARISONS
#             rank_test <- cyto_func_call(
#               "dunn.test::dunn.test",
#               args = list(
#                 lapply(
#                   x,
#                   function(v) {
#                     cyto_data_extract(
#                       v,
#                       parent = z,
#                       channels = channels,
#                       format = "matrix",
#                       coerce = TRUE,
#                       copy = FALSE
#                     )[[1]][[1]][, 1]
#                   }
#                 ),
#                 method = "holm",
#                 alpha = 0.05,
#                 altp = TRUE,
#                 kw = TRUE,
#                 table = FALSE,
#                 list = TRUE
#               )
#             )
#             # STORE P VALUES
#             for(i in 1:(length(x) - 1)) {
#               for(j in (i + 1):length(x)) {
#                 # STORE P VALUE
#                 pval[i, j] <- rank_test[[grep("adjusted", names(rank_test))]][
#                   match(paste(c(i, j), collapse = " - "), rank_test$comparisons)
#                 ]
#                 pval[j, i] <- rank_test[[grep("adjusted", names(rank_test))]][
#                   match(paste(c(i, j), collapse = " - "), rank_test$comparisons)
#                 ]
#                 # COMPUTE DISTANCE
#                 d <- abs(
#                   median(
#                     cyto_data_extract(
#                       x[[i]],
#                       parent = z,
#                       channels = channels,
#                       format = "matrix",
#                       coerce = TRUE,
#                       copy = FALSE
#                     )[[1]][[1]][, 1]
#                   ) - 
#                     median(
#                       cyto_data_extract(
#                         x[[j]],
#                         parent = z,
#                         channels = channels,
#                         format = "matrix",
#                         coerce = TRUE,
#                         copy = FALSE
#                       )[[1]][[1]][, 1]
#                     )
#                 )
#                 # STORE DISTANCE
#                 dst[i, j] <- d
#                 dst[j, i] <- d
#               }
#             }
#           }
#         # UNSUPPORTED TEST FOR COMPARING CROSS ENTROPY VALUES
#         } else {
#           stop(
#             "'base_test' must be either 'ks' or 'rank'!"
#           )
#         }
#         # CDF
#         # HIERARCHICAL CLUSTERING
#         hc <- hclust(
#           as.dist(
#             dst
#           )
#         )
#         # DENDROGRAM
#         if(plot == TRUE) {
#           plot(
#             hc,
#             e = "rectangle",
#             hang = -1,
#             main = paste0(
#               z,
#               " - Dendrogram"
#             ),
#             xlab = "",
#             sub = ""
#           )
#         }
#         # OUTPUT RESULTS
#         list(
#           "p" = pval,
#           "dist" = dst,
#           "hclust" = hc
#         )
#       }
#     ),
#     names = parent
#   )
#   
#   # RETURN CROSS ENTROPY TEST RESULTS
#   return(res)
#   
# }
