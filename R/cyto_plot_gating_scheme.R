## CYTO_PLOT_GATING_SCHEME -----------------------------------------------------

#' Plot Cytometry Gating Strategies
#' 
#' \code{cyto_plot_gating_scheme} automatically plots the entire gating scheme
#' and has full support for gate tracking and back-gating through
#' \code{gate_track} and \code{back_gate}.
#' 
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param gatingTemplate name of the gatingTemplate csv file used to gate
#'   \code{x}. If not supplied the gating scheme will be obtained directly from
#'   the GatingSet.
#' @param back_gate names of the population(s) to back-gate, set to \code{FALSE}
#'   by default to turn off back-gating. To back-gate all populations set this
#'   argument to \code{"all"}.
#' @param gate_track logical indicating whether gate colour should be tracked
#'   throughout gating scheme, set to TRUE by default.
#' @param show_all logical indicating whether every population should be
#'   included in every plot in the gating scheme, set to \code{FALSE} by
#'   default.
#' @param header character string to use as the header for the plot layout, set
#'   to "Gating Scheme" by default.
#' @param title vector of titles to use above each plot.
#' @param popup logical indicating whether the gating scheme should be plotted
#'   in a pop-up window, set to FALSE by default.
#' @param layout a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param point_col colour of points in 2D plots set to NA to use default
#'   red-blue colour scale. Control the colour of overlays by supplying multiple
#'   colours to this argument (e.g. c("blue","red")).
#' @param density_stack numeric [0,1] indicating the degree of offset for 1-D
#'   density distributions with overlay, set to 0.5 by default.
#' @param density_fill fill colour for 1D density distributions. Control the
#'   colour of overlays by supplying multiple colours to this argument (e.g.
#'   c(NA,"red")).
#' @param gate_line_col vector of colours to use for gates. Individual gate
#'   colours can only be controlled when \code{gate_track} is set to TRUE.
#' @param border_line_col line colour for plot border, set to "black" by
#'   default.
#' @param border_line_width line width for plot border, set to 3 when gate_track
#'   is TRUE.
#' @param legend logical indicating whether a legend should be included when an
#'   overlay is supplied.
#' @param legend_text vector of character strings to use for legend when an
#'   overlay is supplied.
#' @param legend_text_size character expansion for legend text, set to 1.2 by
#'   default.
#' @param title_text_col colour for plot title.
#' @param label_text_size numeric to control the size of text in the plot
#'   labels, set to 0.8 by default.
#' @param ... extra arguments passed to \code{\link{cyto_plot}}.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @examples
#' library(CytoRSuiteData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#'
#' # Gate using gate_draw
#' gt <- Activation_gatingTemplate
#' gating(gt, gs)
#'
#' # Gating scheme
#' cyto_plot_gating_scheme(gs[[4]])
#'
#' # Back-gating
#' cyto_plot_gating_scheme(gs[[32]],
#'   back_gate = TRUE
#' )
#'
#' # Gate-tracking
#' cyto_plot_gating_scheme(gs[1:4],
#'   gate_track = TRUE
#' )
#' 
#' @name cyto_plot_gating_scheme
NULL

#' @noRd
#' @export
cyto_plot_gating_scheme <- function(x, ...) {
  UseMethod("cyto_plot_gating_scheme")
}

#' @rdname cyto_plot_gating_scheme
#' @export
cyto_plot_gating_scheme.GatingHierarchy <- function(x,
                                                    gatingTemplate = NULL,
                                                    back_gate = FALSE,
                                                    gate_track = FALSE,
                                                    show_all = FALSE,
                                                    header,
                                                    title,
                                                    popup = FALSE,
                                                    layout,
                                                    point_col,
                                                    density_stack = 0.5,
                                                    density_fill,
                                                    gate_line_col,
                                                    border_line_col = NA,
                                                    border_line_width,
                                                    legend = TRUE,
                                                    legend_text,
                                                    legend_text_size = 1.2,
                                                    title_text_col = NA,
                                                    label_text_size = 0.8, ...) {


  # Set plot method
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "Gating/GatingHierarchy")
  }

  # Assign x to gh
  gh <- x

  # Gating template supplied - apply to GatingHierarchy
  if (!is.null(gatingTemplate)) {
    if (getOption("CytoExploreR_wd_check") == TRUE) {
      if (.file_wd_check(gatingTemplate)) {
        gt <- gatingTemplate(gatingTemplate)
        gating(gt, gh)
      } else {
        stop(paste(
          gatingTemplate,
          "is not in this working directory."
        ))
      }
    } else {
      gt <- gatingTemplate(gatingTemplate)
      gating(gt, gh)
    }
  }

  # Back-gating
  if (back_gate[1] == TRUE) {
    back_gate <- "all"
  }

  # Populations
  pops <- cyto_nodes(gh, path = "auto")

  # Number of Populations
  npop <- length(pops)

  # Back gating
  if (back_gate[1] != FALSE) {
    if (back_gate[1] == "all") {
      overlay <- pops[-1]
    } else {
      overlay <- back_gate
    }
  } else {
    overlay <- NULL
  }

  # Number of overlays
  ovn <- length(overlay)

  # Border thickness
  if (!gate_track) {
    if (missing(border_line_width)) {
      border_line_width <- 1
    }
  } else {
    if (missing(border_line_width)) {
      border_line_width <- 3
    }
  }

  # Colours
  cols <- c(
    "cyan",
    "deepskyblue",
    "navyblue",
    "turquoise4",
    "springgreen3",
    "green",
    "darkgreen",
    "goldenrod4",
    "orange",
    "firebrick2",
    "red",
    "darkred",
    "deeppink2",
    "darkmagenta",
    "purple4",
    "magenta"
  )
  cols <- colorRampPalette(cols)

  # Get colour for each population
  if (gate_track & back_gate != FALSE) {

    # Point col
    if (missing(point_col)) {
      point_col <- c("grey32", cols(npop - 1))

      if (back_gate[1] != "all") {
        point_col[!pops %in% back_gate] <- point_col[1]
      }
    } else {
      point_col <- c(point_col, cols(ovn))[seq_len(npop)]
    }

    # Density fill colour
    if (missing(density_fill)) {
      density_fill <- point_col
    } else {
      density_fill <- c(density_fill, cols(ovn))[seq_len(npop)]
    }

    # Gate line colour
    if (missing(gate_line_col)) {
      gate_line_col <- c("grey32", cols(npop - 1))
    } else {
      gate_line_col <- c(gate_line_col, cols(ovn))[seq_len(npop)]
    }
  } else if (gate_track) {

    # Point colour
    if (missing(point_col)) {
      point_col <- "grey32"

      if (back_gate[1] != "all") {
        point_col[!pops %in% back_gate] <- point_col[1]
      }
    } else {
      point_col <- rep(point_col[1], npop)
    }

    # Density fill colour
    if (missing(density_fill)) {
      density_fill <- point_col
    } else {
      density_fill <- rep(density_fill[1], npop)
    }

    # Gate line colour
    if (missing(gate_line_col)) {
      gate_line_col <- c("grey32", cols(npop - 1)) # include root
    } else {
      gate_line_col <- c(gate_line_col, cols(ovn))[seq_len(npop)]
    }
  } else if (back_gate != FALSE) {

    # Point colour
    if (missing(point_col)) {
      point_col <- c("grey32", cols(npop - 1))

      if (back_gate[1] != "all") {
        point_col[!pops %in% back_gate] <- point_col[1]
      }
    } else {
      point_col <- c(point_col, cols(ovn))[seq_len(npop)]
    }

    # Density fill colour
    if (missing(density_fill)) {
      density_fill <- point_col
    } else {
      density_fill <- c(density_fill, cols(ovn))[seq_len(npop)]
    }

    # Gate line colour
    if (missing(gate_line_col)) {
      gate_line_col <- "red"
    } else {
      gate_line_col <- rep(gate_line_col[1], npop)
    }
  } else {

    # Point colour
    if (missing(point_col)) {
      point_col <- NA
    } else {
      point_col <- point_col[1]
    }

    # Density fill colour
    if (missing(density_fill)) {
      density_fill <- point_col
    } else {
      density_fill <- density_fill[1]
    }

    # Gate line colour
    if (missing(gate_line_col)) {
      gate_line_col <- "red"
    } else {
      gate_line_col <- gate_line_col[1]
    }
  }
  popcols <- data.frame(
    "node" = pops,
    "ptcol" = rep(point_col,
      length.out = npop
    ),
    "gtcol" = rep(gate_line_col,
      length.out = npop
    ),
    "denscol" = rep(density_fill,
      length.out = npop
    ),
    stringsAsFactors = FALSE
  )

  # Header - add spaces to center
  if (missing(header)) {
    header <- "       Gating Scheme"
  }

  # Use GatingHierarchy directly to get gating scheme
  gt <- templateGen(gh)
  gts <- data.frame(
    parent = basename(gt$parent),
    alias = gt$alias,
    xchannel = do.call(
      "rbind",
      strsplit(gt$dims,
        ",",
        fixed = TRUE
      )
    )[, 1],
    ychannel = do.call(
      "rbind",
      strsplit(gt$dims,
        ",",
        fixed = TRUE
      )
    )[, 2],
    stringsAsFactors = FALSE
  )

  # Extract unique parents for plotting
  parents <- unique(gts$parent)

  # Pop-up window?
  cyto_plot_new(popup = popup)

  # Number of plots
  np <- nrow(unique(gts[, c("parent", "xchannel", "ychannel")]))

  # Calculate layout parameters based on number of parents
  if (missing(layout)) {
    if (legend == FALSE) {
      layout <- c(
        n2mfrow(np)[2],
        n2mfrow(np)[1]
      )
    } else {
      if (back_gate == FALSE) {
        layout <- c(
          n2mfrow(np)[2],
          n2mfrow(np)[1]
        )
      } else {
        layout <- c(
          n2mfrow(np + 1)[2],
          n2mfrow(np + 1)[1]
        )
      }
    }
    par(mfrow = layout)
  } else {
    if (layout[1] != FALSE) {
      par(mfrow = layout)
    }
  }

  if (!is.null(header)) {
    par(oma = c(0, 0, 3, 0))
  }

  # Titles
  if (missing(title)) {
    prnts <- parents
    if ("root" %in% prnts) {
      prnts[prnts %in% "root"] <- "All Events"
    }
    title <- prnts
  }

  # New plot per parent
  mapply(function(parents, title) {
    gt <- gts[gts$parent == parents, ]

    # Parent may have gates in different channels
    # construct a plot for each channel set
    lapply(
      seq(1, nrow(unique(gt[, c("xchannel", "ychannel")]))),
      function(x) {
        parent <- as.character(parents)
        xchannel <- as.character(unique(gt[, "xchannel"])[x])
        ychannel <- as.character(unique(gt[, "ychannel"])[x])
        channels <- c(xchannel, ychannel)
        alias <- as.vector(gt[gt$parent == parents &
          gt$xchannel == xchannel &
          gt$ychannel == ychannel, "alias"])

        # Histograms - stacked & gated
        if (channels[1] == channels[2]) {
          channels <- channels[1]
        }

        # Back-gating
        if (back_gate[1] != FALSE) {
          if (back_gate[1] == "all") {
            if (!show_all) {
              overlay <- c(alias, LAPPLY(
                seq_along(alias),
                function(x) {
                  getDescendants(gh, alias[x], path = "auto")
                }
              ))
            }
          } else {
            if (!show_all) {
              if (any(back_gate %in%
                c(alias, LAPPLY(seq_along(alias), function(x) {
                  getDescendants(gh, alias[x], path = "auto")
                })))) {
                ind <- c(
                  alias,
                  LAPPLY(
                    seq_along(alias),
                    function(x) {
                      getDescendants(gh, alias[x], path = "auto")
                    }
                  )
                )
                overlay <- back_gate[back_gate %in% ind]
              } else {
                overlay <- NULL
              }
            }
          }
        }

        # Point colour
        point_col <- popcols[, "ptcol"][match(c(parent, overlay), pops)]

        # Gate line colour
        gate_line_col <- popcols[, "gtcol"][match(alias, pops)]

        # Density fill colour
        density_fill <- popcols[, "denscol"][match(
          c(parent, overlay),
          pops
        )]

        # Border line col
        if (gate_track) {
          if (is.na(border_line_col)) {
            border_line_col <- popcols[, "gtcol"][match(parent, pops)]
          }
        } else {
          if (is.na(border_line_col)) {
            border_line_col <- "black"
          }
        }

        # Title text colour
        if (is.na(title_text_col)) {
          title_text_col <- border_line_col
        }

        # Skip boolean gates
        if (any(
          LAPPLY(alias, function(x) {
            flowWorkspace:::isNegated(gh, x)
          })
        )) {
          message("Skipping boolean gates.")
          alias <- alias[!unlist(
            lapply(
              alias,
              function(x) {
                flowWorkspace:::isNegated(gh, x)
              }
            )
          )]
        }

        # Call to cyto_plot
        cyto_plot(gh,
          parent = parent,
          alias = alias,
          overlay = overlay,
          channels = channels,
          legend = FALSE,
          legend_text = NA,
          title = title,
          point_col = point_col,
          density_stack = density_stack,
          density_fill = density_fill,
          gate_line_col = gate_line_col,
          border_line_col = border_line_col,
          border_line_width = border_line_width,
          title_text_col = title_text_col,
          label_text_size = label_text_size, ...
        )
      }
    )
  }, parents, title)

  # header
  if (!is.null(header)) {
    mtext(header, outer = TRUE, cex = 1, font = 2)
  }

  # Legend
  if (!is.null(overlay)) {
    if (legend == TRUE) {

      # Legend Text
      if (missing(legend_text)) {
        if (class(overlay) == "character") {
          legend_text <- overlay
        } else {
          stop("Please supply vector of names to use in the legend.")
        }
      }

      # Add dummy plot
      plot.new()

      # Add legend
      legend("center",
        legend = legend_text,
        fill = popcols[, "ptcol"][match(overlay, pops)],
        bty = "n",
        cex = legend_text_size,
        x.intersp = 0.5
      )
    }
  }

  # Return default plot layout
  par(mfrow = c(1, 1))
  par(oma = c(0, 0, 0, 0))

  # Turn off graphics device for saving
  if (getOption("cyto_plot_save")) {
    if (inherits(
      x,
      basename(getOption("cyto_plot_method"))
    )) {

      # Close graphics device
      dev.off()

      # Reset cyto_plot_save
      options("cyto_plot_save" = FALSE)

      # Reset cyto_plot_method
      options("cyto_plot_method" = NULL)
    }
  }
}

#' @rdname cyto_plot_gating_scheme
#' @export
cyto_plot_gating_scheme.GatingSet <- function(x,
                                              gatingTemplate = NULL,
                                              group_by,
                                              back_gate = FALSE,
                                              gate_track = FALSE,
                                              show_all = FALSE,
                                              display = NULL,
                                              header = NA,
                                              title = NULL,
                                              popup = FALSE,
                                              layout = NULL,
                                              point_col = NULL,
                                              density_stack = 0.5,
                                              density_fill = NULL,
                                              gate_line_col = NULL,
                                              border_line_col = NULL,
                                              border_line_width = NULL,
                                              legend = TRUE,
                                              legend_text = NULL,
                                              legend_text_size = 1.2,
                                              title_text_col = NULL,
                                              label_text_size = 0.8, ...){
  
  # Set plot method
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "Gating/GatingSet")
  }
  
  # gatingTemplate supplied - apply to GatingSet
  if (!is.null(gatingTemplate)) {
    if (getOption("CytoExploreR_wd_check") == TRUE) {
      if (.file_wd_check(gatingTemplate)) {
        gt <- gatingTemplate(gatingTemplate)
        gating(gt, x)
      } else {
        stop(paste(gatingTemplate, "is not in this working directory"))
      }
    } else {
      gt <- gatingTemplate(gatingTemplate)
      gating(gt, x)
    }
  }
  
  # Back-gating
  if (back_gate[1] == TRUE) {
    back_gate <- "all"
  }
  
  # Extract pData for group_by
  pd <- cyto_details(x)
  
  # group_by all samples by default
  if (missing(group_by)) {
    group_by <- "all"
  }
  
  # Sort pd by group_by colnames
  if (group_by[1] != "all") {
    pd <- pd[do.call("order", pd[group_by]), ]
  }
  
  # Convert GatingSet into list of GatingSets split by group_by
  if (group_by == "all") {
    
    # Add group_by column to pd
    pd$group_by <- rep("all", length(x))
    
    # All GatingSet in same group
    gs.lst <- list(x)
    names(gs.lst) <- "all"
  } else if (length(group_by) == 1) {
    
    # Check group_by is pData variable
    if (!group_by %in% colnames(pd)) {
      stop("group_by should contain the name(s) of pData variable(s).")
    }
    
    # Add group_by column to pd
    pd$group_by <- pd[, group_by]
    
    # Split GatingSet by group_by into list of GatingSets
    gs.lst <- lapply(unique(pd$group_by), function(y) {
      x[pd$name[pd$group_by == y]]
    })
  } else if (length(group_by) > 1) {
    
    # Check group_by is pData variable
    if (!all(group_by %in% colnames(pd))) {
      stop("group_by should contain the name(s) of pData variable(s).")
    }
    
    # Add group_by column to pd
    pd$group_by <- do.call("paste", pd[, group_by])
    
    # Split GatingSet by group_by into list of GatingSets
    gs.lst <- lapply(unique(pd$group_by), function(y) {
      x[pd$name[pd$group_by == y]]
    })
  }
  
  # Plot gating scheme for each gatingSet in gs.lst
  lapply(gs.lst, function(gs) {
    
    # Populations
    pops <- cyto_nodes(gs[[1]], path = "auto")
    
    # Number of Populations
    npop <- length(pops)
    
    # Back gating
    if (back_gate[1] != FALSE) {
      if (back_gate[1] == "all") {
        overlay <- pops[-1]
      } else {
        overlay <- back_gate
      }
    } else {
      overlay <- NULL
    }
    
    # Number of overlays
    ovn <- length(overlay)
    
    # Border thickness
    if (!gate_track) {
      if (is.null(border_line_width)[1]) {
        border_line_width <- 1
      }
    } else {
      if (is.null(border_line_width)[1]) {
        border_line_width <- 3
      }
    }
    
    # Colours
    cols <- c(
      "cyan",
      "deepskyblue",
      "navyblue",
      "turquoise4",
      "springgreen3",
      "green",
      "darkgreen",
      "goldenrod4",
      "orange",
      "firebrick2",
      "red",
      "darkred",
      "deeppink2",
      "darkmagenta",
      "purple4",
      "magenta"
    )
    cols <- colorRampPalette(cols)
    
    # Get colour for each population
    if (gate_track & back_gate != FALSE) {
      
      # Point col
      if (is.null(point_col)) {
        point_col <- c("grey32", cols(npop - 1))
        
        if (back_gate[1] != "all") {
          point_col[!pops %in% back_gate] <- point_col[1]
        }
      } else {
        point_col <- c(point_col, cols(ovn))[seq_len(npop)]
      }
      
      # Density fill colour
      if (is.null(density_fill)) {
        density_fill <- point_col
      } else {
        density_fill <- c(density_fill, cols(ovn))[seq_len(npop)]
      }
      
      # Gate line colour
      if (is.null(gate_line_col)) {
        gate_line_col <- c("grey32", cols(npop - 1))
      } else {
        gate_line_col <- c(gate_line_col, cols(ovn))[seq_len(npop)]
      }
    } else if (gate_track) {
      
      # Point colour
      if (is.null(point_col)) {
        point_col <- "grey32"
        
        if (back_gate[1] != "all") {
          point_col[!pops %in% back_gate] <- point_col[1]
        }
      } else {
        point_col <- rep(point_col[1], npop)
      }
      
      # Density fill colour
      if (is.null(density_fill)) {
        density_fill <- point_col
      } else {
        density_fill <- rep(density_fill[1], npop)
      }
      
      # Gate line colour
      if (is.null(gate_line_col)) {
        gate_line_col <- c("grey32", cols(npop - 1)) # include root
      } else {
        gate_line_col <- c(gate_line_col, cols(ovn))[seq_len(npop)]
      }
    } else if (back_gate != FALSE) {
      
      # Point colour
      if (is.null(point_col)) {
        point_col <- c("grey32", cols(npop - 1))
        
        if (back_gate[1] != "all") {
          point_col[!pops %in% back_gate] <- point_col[1]
        }
      } else {
        point_col <- c(point_col, cols(ovn))[seq_len(npop)]
      }
      
      # Density fill colour
      if (is.null(density_fill)) {
        density_fill <- point_col
      } else {
        density_fill <- c(density_fill, cols(ovn))[seq_len(npop)]
      }
      
      # Gate line colour
      if (is.null(gate_line_col)) {
        gate_line_col <- "red"
      } else {
        gate_line_col <- rep(gate_line_col[1], npop)
      }
    } else {
      
      # Point colour
      if (is.null(point_col)) {
        point_col <- NA
      } else {
        point_col <- point_col[1]
      }
      
      # Density fill colour
      if (is.null(density_fill)) {
        density_fill <- point_col
      } else {
        density_fill <- density_fill[1]
      }
      
      # Gate line colour
      if (is.null(gate_line_col)) {
        gate_line_col <- "red"
      } else {
        gate_line_col <- gate_line_col[1]
      }
    }
    popcols <- data.frame(
      "node" = pops,
      "ptcol" = rep(point_col, length.out = npop),
      "gtcol" = rep(gate_line_col, length.out = npop),
      "denscol" = rep(density_fill, length.out = npop),
      stringsAsFactors = FALSE
    )
    
    # Header - add spaces to center
    if (is.na(header)) {
      if (group_by == "all") {
        header <- paste("       ", "Combined Gating Scheme")
      } else {
        header <- paste(
          "       ",
          pd$group_by[pd$name %in% sampleNames(gs)][1]
        )
      }
    } else {
      header <- paste("       ", header)
    }
    
    # Use GatingSet directly to get gating scheme - template from first member
    gt <- templateGen(gs[[1]])
    gts <- data.frame(
      parent = basename(gt$parent),
      alias = gt$alias,
      xchannel = do.call(
        rbind,
        strsplit(gt$dims, ",",
                 fixed = TRUE
        )
      )[, 1],
      ychannel = do.call(
        rbind,
        strsplit(gt$dims, ",",
                 fixed = TRUE
        )
      )[, 2],
      stringsAsFactors = FALSE
    )
    
    # Extract unique parents for plotting
    parents <- unique(gts$parent)
    
    # Pop-up window?
    if (popup == TRUE) {
      .cyto_plot_window()
    }
    
    # Number of plots
    np <- nrow(unique(gts[, c("parent", "xchannel", "ychannel")]))
    
    # Calculate layout parameters based on number of parents
    if (is.null(layout)) {
      if (legend == FALSE) {
        layout <- c(
          n2mfrow(np)[2],
          n2mfrow(np)[1]
        )
      } else {
        if (back_gate == FALSE) {
          layout <- c(
            n2mfrow(np)[2],
            n2mfrow(np)[1]
          )
        } else {
          layout <- c(
            n2mfrow(np + 1)[2],
            n2mfrow(np + 1)[1]
          )
        }
      }
      par(mfrow = layout)
    } else {
      if (layout[1] != FALSE) {
        par(mfrow = layout)
      }
    }
    
    if (!is.null(header)) {
      par(oma = c(0, 0, 3, 0))
    }
    
    # Titles
    if (is.null(title)) {
      prnts <- parents
      if ("root" %in% prnts) {
        prnts[prnts %in% "root"] <- "All Events"
      }
      title <- prnts
    }
    
    mapply(function(parents, title) {
      gt <- gts[gts$parent == parents, ]
      
      # Parent may have gates in different channels
      # construct a plot for each channel set
      lapply(seq(1, nrow(
        unique(gt[, c("xchannel", "ychannel")])
      )), function(x) {
        parent <- as.character(parents)
        xchannel <- as.character(unique(gt[, "xchannel"])[x])
        ychannel <- as.character(unique(gt[, "ychannel"])[x])
        channels <- c(xchannel, ychannel)
        alias <- as.vector(gt[gt$parent == parents &
                                gt$xchannel == xchannel &
                                gt$ychannel == ychannel, "alias"])
        
        if (channels[1] == channels[2]) {
          channels <- channels[1]
        }
        
        # Back-gating
        if (back_gate[1] != FALSE) {
          if (back_gate[1] == "all") {
            if (!show_all) {
              overlay <- c(alias, LAPPLY(
                seq_along(alias),
                function(x) {
                  getDescendants(gs[[1]], alias[x], path = "auto")
                }
              ))
            }
          } else {
            if (!show_all) {
              if (any(back_gate %in%
                      c(alias, LAPPLY(seq_along(alias), function(x) {
                        getDescendants(gs[[1]], alias[x], path = "auto")
                      })))) {
                overlay <- back_gate[back_gate %in%
                                       c(
                                         alias,
                                         LAPPLY(
                                           seq_along(alias),
                                           function(x) {
                                             getDescendants(
                                               gs[[1]],
                                               alias[x],
                                               path = "auto"
                                             )
                                           }
                                         )
                                       )]
              } else {
                overlay <- NULL
              }
            }
          }
        }
        
        # Point colour
        point_col <- popcols[, "ptcol"][match(
          c(parent, overlay),
          pops
        )]
        
        # Gate line colour
        gate_line_col <- popcols[, "gtcol"][match(alias, pops)]
        
        # Density fill colour
        density_fill <- popcols[, "denscol"][match(
          c(parent, overlay),
          pops
        )]
        
        # Border line col
        if (gate_track) {
          if (is.null(border_line_col)) {
            border_line_col <- popcols[, "gtcol"][match(parent, pops)]
          }
        } else {
          if (is.null(border_line_col)) {
            border_line_col <- "black"
          }
        }
        
        # Title text colour
        if (is.null(title_text_col)) {
          title_text_col <- border_line_col
        }
        
        # Number of events to display
        if (is.null(display)) {
          display <- 1 / length(gs)
        }
        
        # Skip boolean gates
        if (any(LAPPLY(
          alias,
          function(x) {
            flowWorkspace:::isNegated(
              gs[[1]],
              x
            )
          }
        ))) {
          message("skipping boolean gates.")
          alias <- alias[!unlist(
            lapply(alias, function(x) {
              flowWorkspace:::isNegated(gs[[1]], x)
            })
          )]
        }
        
        # Call to cyto_plot
        cyto_plot(gs,
                  parent = parent,
                  group_by = "all",
                  alias = alias,
                  overlay = overlay,
                  channels = channels,
                  legend = FALSE,
                  legend_text = NA,
                  title = title,
                  point_col = point_col,
                  density_stack = density_stack,
                  density_fill = density_fill,
                  gate_line_col = gate_line_col,
                  border_line_col = border_line_col,
                  border_line_width = border_line_width,
                  title_text_col = title_text_col,
                  label_text_size = label_text_size,
                  layout = FALSE,
                  display = display, ...
        )
      })
    }, parents, title)
    
    # header
    if (!is.na(header)) {
      mtext(header, outer = TRUE, cex = 1, font = 2)
    }
    
    # Legend
    if (!is.null(overlay)) {
      if (legend == TRUE) {
        
        # Legend Text
        if (is.null(legend_text)) {
          if (class(overlay) == "character") {
            legend_text <- overlay
          } else {
            stop("Please supply vector of names to use in the legend.")
          }
        }
        
        # Add dummy plot
        plot.new()
        
        # Add legend
        legend("center",
               legend = legend_text,
               fill = popcols[, "ptcol"][match(overlay, pops)],
               bty = "n",
               cex = legend_text_size,
               x.intersp = 0.5
        )
      }
    }
  })
  
  # Return default plot layout
  par(mfrow = c(1, 1))
  par(oma = c(0, 0, 0, 0))
  
  # Turn off graphics device for saving
  if (getOption("cyto_plot_save")) {
    if (inherits(
      x,
      basename(getOption("cyto_plot_method"))
    )) {
      
      # Close graphics device
      dev.off()
      
      # Reset cyto_plot_save
      options("cyto_plot_save" = FALSE)
      
      # Reset cyto_plot_method
      options("cyto_plot_method" = NULL)
    }
  }
}
