## CYTO_PLOT_GATING_SCHEME -----------------------------------------------------

#' Visualise cytometry gating strategies
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{cytoset} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param gatingTemplate name of the gatingTemplate csv file used to gate
#'   \code{x}. If not supplied the gating scheme will be obtained directly from
#'   the GatingHierarchy or GatingSet.
#' @param merge_by a vector of experiment variables to sort and merge samples
#'   into groups prior to plotting, set to "name" by default to prevent merging.
#'   To merge all samples set this argument to \code{TRUE} or \code{"all"}.
#' @param overlay name(s) of the populations in the \code{GatingHierarchy} or
#'   \code{GatingSet} to overlay on every plot in the gating scheme. For
#'   back-gating, \code{overlay} can also be set to \code{"children"},
#'   \code{"descendants"}, \code{"terminal children"} or \code{"terminal
#'   descendants"} in which case the nodes to overlay are determined separately
#'   for every plot based on the \code{parent} population.
#' @param order can be either \code{"groups"} or \code{"nodes"} to control the
#'   order in which the data is plotted. Setting \code{order = "nodes"} will
#'   plot each group on a single page with a plot for each set of nodes. On the
#'   other hand, setting \code{order = "groups"} will plot each set of nodes on
#'   a separate page with plots for each group.
#' @param alias name(s) of the populations for which the gating scheme is to be
#'   plotted, the order of the populations supplied here dictates the plot order
#'   within \code{cyto_plot_gating_scheme()}.
#' @param hidden logical to indicate whether nodes hidden in the
#'   \code{GatingHierarchy} or \code{GatingSet} should be included in the gating
#'   scheme, set to FALSE by default.
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param title a vector of text to label each of the plots in the gating
#'   scheme, set to the name of the parent population in each plot by default.
#' @param title_text_col colour to use in each plot for the title text, set to
#'   \code{"black"} by default. Setting \code{title_test_col = NA} will
#'   automatically match the title colour with parent colour.
#' @param header a vector of text to label each page of plots, set to the name
#'   of each group specified by \code{merge_by} by default.
#' @param header_text_font numeric to control the font of the header text, set
#'   to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param header_text_size numeric to control the size of the header text, set
#'   to 1 by default.
#' @param header_text_col colour to use for the header text, set to "black" by
#'   default.
#' @param legend logical indicating whether a should be included when
#'   populations are overlaid onto the plots, set to TRUE by default.
#'   Alternatively, for gating schemes with many nodes, users can supply a
#'   vector of integers to split the legend across multiple panels, set to 12
#'   nodes by default.
#' @param legend_text vector of text to include in the legend to override the
#'   default labels.
#' @param legend_text_font numeric to control the font of legend text, set to 1
#'   for plain font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param legend_text_size numeric to control the size of text in the legend,
#'   set to 1.2 by default.
#' @param point_col vector of colours to use for each population defined in the
#'   gating scheme, colours are automatically selected from \code{point_cols} if
#'   not supplied.
#' @param point_cols vector of colours from which a colour for each node should
#'   be selected when colours are not supplied manually through
#'   \code{point_col}, resorts to using the default \code{cyto_plot()} colour
#'   palette if not supplied.
#' @param border_line_col colour to use for the border of each plot, set to
#'   "black" by default. Setting \code{border_line_col = NA} will automatically
#'   match the border colour with parent colour.
#' @param gate_line_col colour(s) to use for gates, set to \code{"red"} by
#'   default. Setting \code{gate_line_col = NA} will automatically match the
#'   gate colour with node colour.
#' @param ... additional arguments passed to \code{cyto_plot()}.
#'
#' @importFrom flowWorkspace gh_pop_is_hidden
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \link{cyto_plot}
#' @seealso \link{cyto_plot_gating_tree}
#' @seealso \link{cyto_plot_profile}
#' @seealso \link{cyto_plot_explore}
#' @seealso \link{cyto_plot_map}
#'
#' @export
cyto_plot_gating_scheme <- function(x,
                                    select = NULL,
                                    gatingTemplate = NULL,
                                    merge_by = "name",
                                    overlay = NA,
                                    order = "nodes",
                                    alias = NULL,
                                    hidden = FALSE,
                                    layout,
                                    title,
                                    header,
                                    legend = TRUE,
                                    legend_text = NA,
                                    point_col = NA,
                                    point_cols = NA,
                                    gate_line_col = "red",
                                    border_line_col = "black",
                                    title_text_col = "black",
                                    header_text_font = 2,
                                    header_text_size = 1,
                                    header_text_col = "black",
                                    legend_text_font = 1,
                                    legend_text_size = 1.2,
                                    ...) {
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "scheme")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # PULL DOWN THEME - SET COLOURS
  cyto_plot_theme <- cyto_option("cyto_plot_theme")
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGHIERRACHY | GATINGSET ONLY
  if(!cyto_class(x, "GatingSet")) {
    stop(
      paste0(
        "'cyto_plot_gating_scheme()' only supports objects of class ",
        "GatingHierarchy or GatingSet!"
      )
    )
  }
  
  # RESTRICT DATA
  x <- cyto_select(
    x,
    select
  )
  
  # GATINGTEMPLATE -------------------------------------------------------------
  
  # APPLY GATINGTEMPLATE
  if(!is.null(gatingTemplate)) {
    cyto_gatingTemplate_apply(x, gatingTemplate)
  }
  
  # PREPARE GATINGTEMPLATE -----------------------------------------------------
  
  # EXTRACT GATINGTEMPLATE
  gt <- cyto_gatingTemplate_extract(
    x,
    bool = TRUE,
    drop = TRUE
  )
  
  # UPDATE GATINGTEMPLATE
  gt <- data.frame(
    "parent" = gt$parent,
    "parent_auto" = cyto_nodes_convert(
      x,
      nodes = gt$parent,
      path = "auto"
    ),
    "alias" = gt$alias,
    "alias_auto" = LAPPLY(
      seq_len(nrow(gt)),
      function(z) {
        cyto_nodes_convert(
          x,
          nodes = gt$alias[z],
          anchor = gt$parent[z],
          path = "auto"
        )
      }
    ),
    dims = gt$dims,
    stringsAsFactors = FALSE
  )
  
  # EXCLUDE HIDDEN NODES
  if(!hidden) {
    gt <- gt[
      LAPPLY(
        gt[, "alias_auto"],
        function(z) {
          # ONLY EXCLUDE ALIAS BUT KEEP PARENT FOR DESCENDANTS
          if(gh_pop_is_hidden(x[[1]], z)) {
            message(
              paste0(
                z,
                " is a hidden node in this ",
                cyto_class(x),
                " and will be excluded from the gating scheme."
              )
            )
            return(FALSE)
          } else {
            return(TRUE)
          }
        }
      ), , drop = FALSE
    ]
  }
  
  # SUBSET BY ALIAS
  if(!is.null(alias)) {
    # VALID AUTO NODES
    alias <- cyto_nodes_convert(
      x,
      nodes = alias,
      hidden = hidden,
      path = "auto"
    )
    # SUBSET GATINGTEMPLATE
    gt <- gt[gt[, "alias_auto"] %in% alias, , drop = FALSE]
  }

  # SPLIT BY PARENT & DIMS
  gt_split <- structure(
    lapply(
      unique(gt$parent),
      function(z) {
        gt_chunk <- gt[gt$parent %in% z, , drop = FALSE]
        structure(
          lapply(
            unique(gt_chunk$dims),
            function(w) {
              gt_chunk[gt_chunk$dims %in% w, , drop = FALSE]
            }
          ),
          names = unique(gt_chunk$dims)
        )
      }
    ),
    names = unique(gt$parent)
  )
  gt_split <- unlist(gt_split, recursive = FALSE, use.names = FALSE)
  
  # ORDER NODES BY ALIAS
  if(!is.null(alias)) {
    ind <- c()
    lapply(
      alias,
      function(z) {
        ind <<- c(
          ind, 
          LAPPLY(
            seq_along(gt_split),
            function(v) {
              if(z %in% gt_split[[v]][, "alias_auto"]) {
                return(v)
              } else {
                return(NULL)
              }
            }
          )
        )
      }
    )
    gt_split <- gt_split[unique(ind)]
  }
  
  # NODE PROPERTIES ------------------------------------------------------------
  
  # NODES
  nodes <- unique(
    LAPPLY(
      seq_len(nrow(gt)),
      function(z) {
        c(gt$parent_auto[z], gt$alias_auto[z])
      }
    )
  )
  
  # NODE COUNT
  node_count <- length(nodes)
  
  # OVERLAY
  if(!.all_na(overlay)) {
    if(!is.character(overlay)) {
      stop(
        paste0(
          "'cyto_plot_gating_scheme()' only supports the names of populations ",
          "through 'overlay'!"
        )
      )
    }
    # OVERLAY - CHECK POPULATION NAMES
    if(!all(overlay %in% c("descendants",
                           "children", 
                           "terminal descendants", 
                           "terminal children"))) {
      overlay <- cyto_nodes_convert(
        x,
        nodes = overlay,
        path = "auto"
      )
    }
  }

  # NO OVERLAY 
  if(.all_na(overlay)) {
    # PSEUDOCOLOUR DEFAULT
    node_cols <- rep(point_col, length.out = node_count)
    # NO LEGEND REQUIRED
    legend <- FALSE
  # OVERLAY
  } else {
    # DEFAULT POINT COLOURS
    if(.all_na(point_col)) {
      # CHECK THEME
      if("point_cols" %in% names(cyto_plot_theme)) {
        point_cols <- cyto_plot_theme[["point_cols"]]
      } else {
        point_cols <- .cyto_plot_colour_palette("point_cols")
      }
      # COLOUR PALETTE
      cols <- colorRampPalette(point_cols)
      # NODE COLOURS
      node_cols <- cols(node_count)
    # CUSTOM POINT COLOURS
    } else {
      node_cols <- rep(point_col, length.out = node_count)
    }
  }
  
  # NODE PROPERTIES
  node_props <- data.frame(
    "node" = nodes,
    "point_col" = node_cols,
    "gate_line_col" = if(.all_na(gate_line_col)) {
      node_cols
    } else {
      rep(gate_line_col, length.out = length(nodes))
    },
    "hist_fill" = node_cols,
    "border_line_col" = if(.all_na(border_line_col)) {
      node_cols
    } else {
      rep(border_line_col, length.out = length(nodes))
    },
    "title_text_col" = if(.all_na(title_text_col)) {
      node_cols
    } else {
      rep(title_text_col, length.out = length(nodes))
    },
    stringsAsFactors = FALSE
  )
  
  # DEFAULT LEGEND_TEXT
  if(.all_na(legend_text)) {
    legend_text <- nodes
    legend_text[legend_text %in% "root"] <- "All Events"
  # CUSTOM LEGEND_TEXT
  } else {
    legend_text <- rep(
      legend_text, 
      length.out = length(nodes)
    )
  }
  names(legend_text) <- nodes
  
  # LEGEND
  if(!all(legend %in% FALSE)) {
    # DEFAULT SPLIT - 12 NODES
    if(all(legend %in% TRUE)) {
      ind <- c(
        rep(12, floor(length(legend_text)/12)), 
        length(legend_text) %% 12
      )
    # CUSTOM SPLIT
    } else {
      # SINGLE SPLIT
      if(length(legend) == 1) {
        ind <- c(
          rep(legend, floor(length(legend_text)/legend)), 
          length(legend_text) %% legend
        )
      # VARIED SPLITS
      } else {
        # INCORRECT LEGEND LENGTH
        if(sum(legend) != length(legend_text)) {
          stop(
            paste0(
              "'legend' must add up to a total of ",
              length(legend_text), "!"
            )
          )
        }
        ind <- legend
      }
    }
    # SPLIT LEGEND
    ind <- rep(1:length(ind), times = ind)
    legend <- split(legend_text, ind)
    # FORMAT LEGEND - DOUBLE PANELS
    cnt <- 0
    ind <- c()
    LAPPLY(
      names(legend), 
      function(z){
        if(as.numeric(z) %% 2 > 0) {
          cnt <<- cnt + 1
          ind <<- c(ind, cnt)
        } else {
          ind <<- c(ind, cnt)
        }
      }
    )
    legend <- split(legend, ind)
  # NO LEGEND
  } else {
    legend = NULL
  }
  
  # PREPARE PLOT ARGUMENTS -----------------------------------------------------
  
  # PREPARE SAMPLES - LIST OF GATINGSETS
  x_list <- cyto_group_by(
    x,
    merge_by
  )
  
  # LAYOUT 
  if(missing(layout)) {
    # SWITCH LAYOUTS - SINGLE GROUP
    if(length(x_list) == 1) {
      order  <- "nodes"
    # SWITCH LAYOUTS - SINGLE NODE
    } else if(length(nodes) == 1) {
      order <- "groups"
    }
    # LAYOUT - NODE ORDER
    if(.grepl("^n", order)) {
      layout <- .cyto_plot_layout(
        rep(NA, length(gt_split) + length(legend))
      )
    # LAYOUT - GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(
        rep(NA, length(x_list) + length(legend))
      )
    }
  }
  
  # TOTAL PLOTS PER PAGE - BASED ON LAYOUT
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES - NODES / LAYOUT
  if(.grepl("^n", order)) {
    # NUMBER OF GROUPS TO PLOT
    n <- length(x_list)
    # PAGES PER GROUP
    pg <- ceiling(length(gt_split)/np)
    # TOTAL PAGES
    tpg <- n * pg
    # TOTAL PLOTS
    tp <- n * length(gt_split)
  # PAGES - GROUPS / LAYOUT
  } else {
    # NUMBER OF GROUPS TO PLOT
    n <- length(gt_split)
    # PAGES PER GROUP
    pg <- ceiling(length(x_list)/np)
    # TOTAL PAGES
    tpg <- n * pg
    # TOTAL PLOTS
    tp <- n * length(x_list)
  }
  
  # DEFAULT HEADERS 
  if(missing(header)) {
    # SEPARATE GATING SCHEMES
    if(.grepl("^n", order)) {
      header <- rep(
        names(x_list), 
        each = pg
      )
    # GATING SCHEME STEPS - LABEL WITH PARENT & CHANNELS
    } else {
      header <- rep(
        LAPPLY(
          gt_split,
          function(z) {
            paste0(
              cyto_nodes_convert(
                x,
                nodes = unique(z$parent),
                path = "auto"
              ),
              " ",
              paste(
                cyto_markers_extract(
                  x,
                  channels = unlist(
                    strsplit(
                      unique(z$dims),
                      ","
                    )
                  ),
                  append = TRUE
                ),
                collapse = " | "
              )
            )
          }
        ),
        each = pg
      )
    }
  # CUSTOM HEADERS
  } else {
    # HEADER PER GROUP
    if(length(header) == n) {
      header <- rep(header, each = pg)
    # REPEAT HEADERS
    } else {
      header <- rep(header, length.out = tpg)
    }
  }
  
  # REPEAT HEADER ARGUMENTS
  header_text_font <- rep(header_text_font, length.out = length(header))
  header_text_size <- rep(header_text_size, length.out = length(header))
  header_text_col <- rep(header_text_col, length.out = length(header))
  
  # DEFAULT TITLES
  if(missing(title)) {
    # SEPARATE GATING SCHEMES
    if(.grepl("^n", order)) {
      title <- rep(
        cyto_nodes_convert(
          x,
          nodes = LAPPLY(
            gt_split,
            function(z) {
              unique(z$parent)
            }
          ),
          path = "auto"
        ),
        times = length(x_list)
      )
    # GATING SCHEME STEPS
    } else {
      title <- LAPPLY(
        names(x_list),
        function(z) {
          paste0(
            z,
            "\n",
            cyto_nodes_convert(
              x,
              nodes = LAPPLY(
                gt_split,
                function(w) {
                  unique(w$parent)
                }
              ),
              path = "auto"
            )
          )
        }
      )
    }
  # CUSTOM TITLES
  } else {
    # TITLES PER GROUP
    if(length(title) == n) {
      title <- rep(title, each = pg)
    # REPEAT TITLES
    } else {
      title <- rep(title, length.out = tp)
    }
  }
  
  # PLOT CONSTRUCTION ----------------------------------------------------------
  
  # PROGRESS BAR - INCREMENTED BY CYTO_PLOT()
  pb <- cyto_progress(
    label = "cyto_plot_gating_scheme()",
    total = tp
  )
  
  # CALL CYTO_PLOT() - PER GROUP
  if(grepl("^n", order, ignore.case = TRUE)) {
    # LOOP THROUGH GROUPS
    plots <- structure(
      lapply(
        seq_along(x_list),
        function(z) {
          # GATINGSET
          gs <- x_list[[z]]
          # HEADER COUNTER
          header_cnt <- (z - 1) * pg
          # TITLE COUNTER
          title_cnt <- 0
          # CONSTRUCT PLOTS
          p <- structure(
            lapply(
              seq_along(gt_split),
              function(w) {
                # TITLE COUNTER
                title_cnt <<- title_cnt + 1
                # PARENT
                parent <- unique(gt_split[[w]]$parent_auto)
                # ALIAS
                alias <- unique(gt_split[[w]]$alias_auto)
                # CHANNELS
                channels <- unlist(strsplit(unique(gt_split[[w]]$dims), ","))
                # OVERLAY
                pops <- NULL
                if(!.all_na(overlay)) {
                  # DESCENDANTS | CHILDREN 
                  if(all(overlay %in% c("descendants", 
                                        "children",
                                        "terminal descendants",
                                        "terminal children"))) {
                    overlay <- cyto_nodes_kin(
                      gs,
                      nodes = parent,
                      type = overlay,
                      path = "auto",
                      hidden = hidden
                    )
                  }
                  pops <- overlay
                }
                # CALL CYTO_PLOT
                cyto_plot(
                  gs,
                  parent = parent,
                  alias = alias,
                  overlay = overlay,
                  channels = channels,
                  merge_by = "all", # TODO: FIX TITLE
                  layout = layout,
                  title = title[title_cnt],
                  header = if(!.all_na(header)) {
                    ""
                  } else {
                    NA
                  }, # HEADER ADDED MANUALLY
                  page = FALSE, # PAGE MANUALLY
                  legend = FALSE,
                  point_col = node_props[
                    node_props$node %in% c(parent, pops), "point_col"
                  ],
                  hist_fill = node_props[
                    node_props$node %in% c(parent, pops), "hist_fill"
                  ],
                  gate_line_col = node_props[
                    node_props$node %in% alias, "gate_line_col"
                  ],
                  border_line_col = node_props[
                    node_props$node %in% parent, "border_line_col"
                  ],
                  title_text_col = node_props[
                    node_props$node %in% parent, "title_text_col"
                  ],
                  ...
                )
                # LEGEND
                if(w == length(gt_split) & cyto_class(legend, "list", TRUE)) {
                  #CONSTRUCT LEGEND
                  lapply(
                    seq_along(legend),
                    function(v) {
                      # OPEN NEW DEVICE
                      if(.par("page")[[1]]) {
                        cyto_plot_new()
                      }
                      # CREATE NEW PLOT
                      plot.new()
                      # SINGLE CENTRAL LEGEND
                      if(length(legend) == 1 & length(legend[[v]]) == 1) {
                        legend(
                          "center",
                          legend = legend[[v]][[1]],
                          fill = node_props$point_col[
                            match(names(legend[[v]][[1]]), node_props$node)   
                          ],
                          bty = "n",
                          cex = legend_text_size,
                          x.intersp = 0.5
                        )
                      # MULTIPLE LEGENDS
                      } else {
                        # COMPUTE PLOT BOUNDARIES
                        usr <- .par("usr")[[1]]
                        rng_usr <- mapply(".line2user", .par("mar")[[1]], 1:4)
                        # CONSTRUCT LEGENDS
                        lapply(
                          seq_along(legend[[v]]),
                          function(g) {
                            # CREATE LEFT LEGEND
                            if(g == 1) {
                              legend(
                                "left",
                                inset = c(
                                  -0.88 * abs(
                                    diff(c(rng_usr[1], usr[1]))
                                  )/abs(diff(usr[1:2])),
                                  0
                                ),
                                legend = legend[[v]][[g]],
                                fill = node_props$point_col[
                                  match(
                                    names(legend[[v]][[g]]), node_props$node
                                  )
                                ],
                                bty = "n",
                                cex = legend_text_size,
                                x.intersp = 0.5,
                                xpd = TRUE
                              )
                            # CREATE RIGHT LEGEND 
                            } else {
                              legend(
                                "right",
                                inset = c(
                                  -0.88 * abs(
                                    diff(c(rng_usr[2], usr[2]))
                                  )/abs(diff(usr[1:2])) + usr[2],
                                  0
                                ),
                                legend = legend[[v]][[g]],
                                fill = node_props$point_col[
                                  match(
                                    names(legend[[v]][[g]]), node_props$node
                                  )
                                ],
                                bty = "n",
                                cex = legend_text_size,
                                x.intersp = 0.5,
                                xpd = TRUE
                              )
                            }
                          }
                        )
                      }
                    }
                  )
                }
                # PAGE | HEADER | RECORD
                rec <- NULL
                if(.par("page")[[1]] | w == length(gt_split)) {
                  # HEADER INDEX
                  ind <- header_cnt + ceiling(w/np)
                  # HEADER
                  if(!.all_na(header[ind])) {
                    .cyto_plot_header(
                      header[ind],
                      header_text_font = header_text_font[ind],
                      header_text_size = header_text_size[ind],
                      header_text_col = header_text_col[ind]
                    )
                  }
                  # NEW PAGE
                  cyto_plot_new_page()
                  # RECORD
                  rec <- cyto_plot_record()
                }
                return(rec)
              }
            ),
            names = NULL
          )
          # PREPARE PLOTS
          p[LAPPLY(p, "is.null")] <- NULL
          return(p)
        }
      ),
      names = names(x_list)
    )
  # CALL CYTO_PLOT() - PER GATING STEP
  } else {
    # LOOP THROUGH GATING STEPS
    plots <- structure(
      lapply(
        seq_along(gt_split),
        function(z) {
          # HEADER COUNTER
          header_cnt <- (z - 1) * pg
          # TTLE COUNTER
          title_cnt <- 0
          # GATINGTEMPLATE CHUNK
          gt_chunk <- gt_split[[z]]
          # CONSTRUCT PLOTS PER GROUP
          p <- structure(
            lapply(
              seq_along(x_list),
              function(v) {
                # TITLE COUNTER
                title_cnt <<- title_cnt + 1
                # GATINGHIERARCHY | GATINGSET
                gs <- x_list[[v]]
                # PARENT 
                parent <- unique(gt_chunk[, "parent_auto"])
                # ALIAS
                alias <- unique(gt_chunk[, "alias_auto"])
                # CHANNELS
                channels <- unlist(strsplit(unique(gt_chunk[, "dims"]), ","))
                # POPULATIONS TO OVERLAY 
                pops <- NULL
                # OVERLAY
                if(!.all_na(overlay)) {
                  # DESCENDANTS | CHILDREN
                  if(all(overlay %in% c("descendants",
                                        "children",
                                        "terminal descendants",
                                        "terminal children"))) {
                    overlay <- cyto_nodes_kin(
                      gs,
                      nodes = parent,
                      type = overlay,
                      path = "auto",
                      hidden = hidden
                    )
                  }
                  pops <- overlay
                }
                # CALL CYTO_PLOT()
                cyto_plot(
                  gs,
                  parent = parent,
                  alias = alias,
                  overlay = overlay,
                  channels = channels,
                  merge_by = "all", # TODO: FIX TITLE
                  layout = layout,
                  title = title[title_cnt],
                  header = if(!.all_na(header)) {
                    ""
                  } else {
                    NA
                  }, # HEADER ADDED MANUALLY
                  page = FALSE, # PAGE MANUALLY
                  legend = FALSE,
                  point_col = node_props[
                    node_props$node %in% c(parent, pops), "point_col"
                  ],
                  hist_fill = node_props[
                    node_props$node %in% c(parent, pops), "hist_fill"
                  ],
                  gate_line_col = node_props[
                    node_props$node %in% alias, "gate_line_col"
                  ],
                  border_line_col = node_props[
                    node_props$node %in% parent, "border_line_col"
                  ],
                  title_text_col = node_props[
                    node_props$node %in% parent, "title_text_col"
                  ],
                  ...
                )
                # LEGEND - FULL LEGEND ON EACH PAGE
                if(v == length(x_list) & cyto_class(legend, "list", TRUE)) {
                  # CONSTRUCT LEGEND
                  lapply(
                    seq_along(legend),
                    function(w) {
                      # OPEN NEW DEVICE
                      if(.par("page")[[1]]) {
                        cyto_plot_new()
                      }
                      # CREATE NEW PLOT
                      plot.new()
                      # SINGLE CENTRAL LEGEND
                      if(length(legend) == 1 & length(legend[[w]]) == 1) {
                        legend(
                          "center",
                          legend = legend[[w]][[1]],
                          fill = node_props$point_col[
                            match(names(legend[[w]][[1]]), node_props$node)   
                          ],
                          bty = "n",
                          cex = legend_text_size,
                          x.intersp = 0.5
                        )
                      # MULTIPLE LEGENDS
                      } else {
                        # COMPUTE PLOT BOUNDARIES
                        usr <- .par("usr")[[1]]
                        rng_usr <- mapply(".line2user", .par("mar")[[1]], 1:4)
                        # CONSTRUCT LEGENDS
                        lapply(
                          seq_along(legend[[w]]),
                          function(g) {
                            # CREATE LEFT LEGEND
                            if(g == 1) {
                              legend(
                                "left",
                                inset = c(
                                  -0.88 * abs(
                                    diff(c(rng_usr[1], usr[1]))
                                  )/abs(diff(usr[1:2])),
                                  0
                                ),
                                legend = legend[[w]][[g]],
                                fill = node_props$point_col[
                                  match(
                                    names(legend[[w]][[g]]), node_props$node
                                  )
                                ],
                                bty = "n",
                                cex = legend_text_size,
                                x.intersp = 0.5,
                                xpd = TRUE
                              )
                              # CREATE RIGHT LEGEND 
                            } else {
                              legend(
                                "right",
                                inset = c(
                                  -0.88 * abs(
                                    diff(c(rng_usr[2], usr[2]))
                                  )/abs(diff(usr[1:2])) + usr[2],
                                  0
                                ),
                                legend = legend[[w]][[g]],
                                fill = node_props$point_col[
                                  match(
                                    names(legend[[w]][[g]]), node_props$node
                                  )
                                ],
                                bty = "n",
                                cex = legend_text_size,
                                x.intersp = 0.5,
                                xpd = TRUE
                              )
                            }
                          }
                        )
                      }
                    }
                  )
                }
                # PAGE | HEADER | RECORD
                rec <- NULL
                if(.par("page")[[1]] | v == length(x_list)) {
                  # HEADER INDEX
                  ind <- header_cnt + ceiling(v/np)
                  # HEADER
                  if(!.all_na(header[ind])) {
                    .cyto_plot_header(
                      header[ind],
                      header_text_font = header_text_font[ind],
                      header_text_size = header_text_size[ind],
                      header_text_col = header_text_col[ind]
                    )
                  }
                  # NEW PAGE
                  cyto_plot_new_page()
                  # RECORD
                  rec <- cyto_plot_record()
                }
                return(rec)
              }
            ),
            names = names(x_list)
          )
          # PREPARE PLOTS
          p[LAPPLY(p, "is.null")] <- NULL
          return(p)
        }
      ),
      names = names(gt_split)
    )
    
  }

  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}