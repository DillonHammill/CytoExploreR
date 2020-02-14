## CYTO_PLOT_LABEL HELPERS -----------------------------------------------------

# A collection of functions to calculate statistics and offset locations for
# labels.

## .CYTO_LABEL_POPS ------------------------------------------------------------

#' Get a list of gated populations to label
#'
#' @param x object of class flowFrame
#' @param gate list of gate objects to apply to x
#' @param negate logical indicating whether negated population should be
#'   included.
#'
#' @return list of flowFrames
#'
#' @importFrom flowCore Subset split quadGate
#' @importFrom methods is
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_label_pops <- function(x,
                             gate,
                             negate = FALSE) {

  # NO GATES -------------------------------------------------------------------

  # RETURN X
  if (.all_na(gate)) {
    return(list(x))
  }

  # PREPARE GATES --------------------------------------------------------------

  # LIST OF GATE OBJECTS
  if (is(gate)[1] == "list") {
    if (all(LAPPLY(gate, "is") %in% c(
      "rectangleGate",
      "polygonGate",
      "ellipsoidGate",
      "quadGate",
      "filters"
    ))) {
      gate <- unlist(gate)
    }
  } else if (is(gate)[1] == "filters") {
    gate <- unlist(gate)
  } else if (is(gate)[1] %in% c(
    "rectangleGate",
    "polygonGate",
    "ellipsoidGate",
    "quadGate",
    "filters"
  )) {
    gate <- list(gate)
  }

  # List of RectangleGates to QuadGate (MUST BE 4 RECTANGLES)
  if(length(gate) == 4 & all(LAPPLY(gate, function(z){
    is(z, "rectangleGate") & any(grepl("quad", names(attributes(z))))
    }))){
      # CHANNELS
      chans <- as.character(parameters(gate[[1]]))
      quad_order <- LAPPLY(gate, function(z){z@filterId})
      gate <- list(.cyto_gate_quad_convert(gate, channels = chans))
  }
  
  # GATES ----------------------------------------------------------------------

  # NEGATED GATE - QUADGATES EXCLUDED
  if (negate == TRUE & !any(LAPPLY(gate, function(z) {
    is(z, "quadGate")
  }))) {
    if (length(gate) > 1) {
      gate <- c(gate, list(do.call("|", gate)))
    } else {
      gate <- c(gate, gate)
    }
  }
  
  # POPULATIONS ----------------------------------------------------------------

  # CAREFUL CANNOT NEGATE INDIVIDUAL QUADRANTS EITHER
  
  # ARGUMENTS
  args <- .args_list()
  
  # GATING
  pops <- LAPPLY(seq_len(length(gate)), function(z) {
    # NEGATED POPULATION
    if (negate == TRUE & z == length(gate)) {
      split(x, gate[[z]])[[2]]
      # GATED POPULATIONS
    } else {
      # QUADGATES RETURN MULTIPLE POPULATIONS
      if (is(gate[[z]], "quadGate")) {
        if("quad_order" %in% names(args)){
          quads <- unlist(strsplit(gate[[z]]@filterId, "|"))
          quads <- quads[quads != "|"]
          split(x, gate[[z]])[c(2, 1, 3, 4)][match(quad_order,
                                                   quads)]# FIX ORDER
        }else{
          split(x, gate[[z]])[c(2, 1, 3, 4)]# FIX ORDER
        }
        # SINGLE POPULATIONS
      } else {
        # RECTANGLE BELONGS TO QUADGATE
        if(is(gate[[z]], "rectangleGate") &
           any(grepl("quad", names(attributes(gate[[z]]))))){
          q <- names(attributes(gate[[z]])[["quadrants"]])
          coords <- .cyto_gate_coords(gate[z], 
                                      channels = as.character(parameters(gate[[z]])))
          chans <- colnames(coords)
          coords <- lapply(colnames(coords), function(y){
            unique(coords[, y][is.finite(coords[, y])])
          })
          names(coords) <- chans
          qg <- quadGate(filterId = paste(q, collapse = "|"), 
                         .gate = coords)
          p <- split(x, qg)[c(2, 1, 3, 4)] # FIX ORDER
          names(p) <- q
          p[[match(gate[[z]]@filterId, names(p))]]
        }else{
          Subset(x, gate[[z]])
        }
      }
    }
  })

  # RETURN LIST OF GATED POPULATIONS
  return(pops)
}


## .CYTO_LABEL_STAT ------------------------------------------------------------

#' Compute and prepare statistics for labels
#'
#' @param x list of parental flowFrame objects.
#' @param pops list of populations to label.
#' @param channels vector of channels used to construct the plot.
#' @param label_stat names of statistics to include in labels. Supplied per
#'   layer.
#'
#' @return computed statistics to include in labels.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_label_stat <- function(x,
                             pops,
                             channels,
                             axes_trans = NA,
                             label_stat,
                             density_smooth) {

  # CHECKS ---------------------------------------------------------------------

  # NO LABEL_STAT
  if (.all_na(label_stat)) {
    return(label_stat)
  }

  # VALID LABEL_STAT
  label_stat <- LAPPLY(label_stat, ".cyto_stat_check")

  # SUPPORTED LABEL_STAT
  if (length(channels) == 2) {
    if (!any(label_stat %in% c("freq", "count"))) {
      stop("Only count and frequency statistics are supported in 2D plots.")
    }
  }

  # GENERAL --------------------------------------------------------------------

  # SAMPLES
  SMP <- length(x)

  # POPULATIONS PER LAYER
  NP <- length(pops) / SMP

  # TOTAL POPULATIONS - SPLIT INDICES
  TNP <- seq_len(NP * SMP)

  # SPLIT TNP
  TNP <- split(TNP, rep(seq_len(SMP), each = NP))

  # COMPUTE LABEL_STAT ---------------------------------------------------------

  # STATISTICS
  LABEL_STAT <- LAPPLY(seq_len(SMP), function(z) {
    # LABEL_STAT
    ST <- lapply(TNP[[z]], function(y) {
      # STATISTIC SUPPLIED
      if (!.all_na(label_stat[y])) {
        # FREQUENCY STATISTIC
        if (grepl("freq", label_stat[y], ignore.case = TRUE)) {
          # PERCENT OF BASE LAYER (WITHOUT GATES)
          if(length(pops) == length(x)){
            st <- .cyto_count(pops[[y]]) / .cyto_count(x[[1]]) * 100
          # PERCENT GATED PER LAYER
          }else{
            st <- .cyto_count(pops[[y]]) / .cyto_count(x[[z]]) * 100
          }
          st <- paste(.round(st, 2), "%")
          # CV STATISTIC
        } else if (grepl("CV", label_stat[y], ignore.case = TRUE)) {
          st <- suppressMessages(
            .cyto_CV(pops[[y]],
              channels = channels,
              trans = axes_trans
            )
          )
          st <- paste(.round(st, 2), "%")
          # OTHER STATISTIC
        } else {
          st <- suppressMessages(
            cyto_stats_compute(pops[[y]],
              channels = channels,
              trans = axes_trans,
              stat = label_stat[y],
              format = "long",
              density_smooth = density_smooth
            )
          )
          st <- as.numeric(st[, ncol(st)])
          if (!grepl("count", label_stat[y], ignore.case = TRUE)) {
            st <- .round(st, 2)
          }
        }
        # NO STATISTIC
      } else {
        st <- NA
      }
      return(st)
    })
    return(ST)
  })

  # RETURN COMPUTED STATISTICS -------------------------------------------------
  return(LABEL_STAT)
}

## .CYTO_LABEL_TEXT ------------------------------------------------------------

#' Prepare text for labels to include statistics
#'
#' @param label_text text to include in the labels.
#' @param label_stat vector of computed statistics to include in the labels.
#'
#' @return vector of finalised labels incorporating both label_text and
#'   label_stat.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_label_text <- function(label_text,
                             label_stat) {

  # MERGE LABEL_TEXT & LABEL_STAT
  LAPPLY(seq_len(length(label_text)), function(z) {
    # NO LABEL_TEXT & NO LABEL_STAT
    if (.all_na(label_text[z]) & .all_na(label_stat[z])) {
      label_text[z] <<- NA
      # NO LABEL_TEXT & LABEL_STAT
    } else if (.all_na(label_text[z]) & !.all_na(label_stat[z])) {
      label_text[z] <<- label_stat[z]
      # LABEL_TEXT & NO LABEL_STAT
    } else if (!.all_na(label_text[z]) & .all_na(label_stat[z])) {
      label_text[z] <<- label_text[z]
      # LABEL_TEXT & LABEL_STAT
    } else if (!.all_na(label_text[z]) & !.all_na(label_stat[z])) {
      label_text[z] <<- paste(label_text[z], label_stat[z], sep = "\n")
    }
  })

  # RETURN LABEL_TEXT
  return(label_text)
}

## .CYTO_LABEL_COORDS ----------------------------------------------------------

#' Compute offset label co-ordinates
#'
#' Used internally within cyto_plot to compute offset label co-ordinates. Only
#' called if label_position and label are set to TRUE.
#'
#' @param args list of named cyto_plot arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @noRd
.cyto_label_coords <- function(args) {

  # GRAPHICAL PARAMETERS -------------------------------------------------------

  # PLOT LIMITS
  lims <- par("usr")

  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04) / 2 # 2% BUFFER
  xmin <- xmin + xpad
  xmax <- xmax - xpad
  xrng <- xmax - xmin

  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04) / 2 # 2% BUFFER
  ymin <- ymin + ypad
  ymax <- ymax - ypad
  yrng <- ymax - ymin

  # DENSITY_STACK --------------------------------------------------------------

  # 1D PLOT
  if (length(args[["channels"]]) == 1) {
    # NO STACKING
    if (args[["density_stack"]] == 0) {
      # DEFAULT Y COORD - 50% Y RANGE
      y_coords <- rep(0.5 * yrng, length(args[["label_text"]]))
      # STACKING
    } else if (args[["density_stack"]] != 0) {
      # Y_MAX PER LAYER
      y_max <- as.numeric(unlist(
        strsplit(names(args[["fr_dens_list"]])[1], "-")
      )[2])
      # STACK LEVELS
      stk <- LAPPLY(seq(0, args[["SMP"]], 1), function(z) {
        args[["density_stack"]] * y_max * z
      })
      # DEFAULT LABEL Y CO-ORDINATES
      y_coords <- LAPPLY(seq_len(args[["SMP"]]), function(z) {
        # HALFWAY BETWEEN HORIZONTAL LINES
        stk[z] + 0.5 * stk[2]
      })
      # REPEAT (Y_COORDS/LAYER)
      y_coords <- rep(y_coords, each = args[["NP"]])
    }
  }

  # COMPUTE LABEL CO-ORDINATES -------------------------------------------------

  # SPLIT TNP LABELS
  label_ind <- split(
    seq_len(args[["TNP"]]),
    rep(seq_len(args[["SMP"]]), each = args[["NP"]])
  )

  # GATE CENTERS - (IGNORE SUPPLIED LABEL COORDS)
  if (!.all_na(args[["gate"]])) {
    gate_centers <- .cyto_gate_center(args[["gate"]],
      channels = args[["channels"]],
      text_x = rep(NA, .cyto_gate_count(args[["gate"]], total = TRUE)),
      text_y = rep(NA, .cyto_gate_count(args[["gate"]], total = TRUE))
    )
  }
  
  # LABEL_TEXT_XY - MATRIX
  label_text_xy <- lapply(seq_len(args[["SMP"]]), function(z) {

    # TEMPORARY STORAGE VECTORS
    text_x <- c()
    text_y <- c()
    # COMPUTE LABEL CO-ORDINATES
    lapply(label_ind[[z]], function(y) {
      # LABEL INDEX (GATE & TEXT_X/Y)
      ind <- match(y, label_ind[[z]])
      # LABEL
      if (!.all_na(args[["label_text"]][y])) {
        # ID PLOT - CENTER OF RANGE
        if (length(args[["channels"]]) == 1) {
          # GATE
          if (!.all_na(args[["gate"]])) {
            # GATED POPULATION
            if (ind <= nrow(gate_centers)) {
              # X COORD - GATE CENTER
              if(.all_na(args[["label_text_x"]][y])){
                text_x[ind] <<- gate_centers[ind, "x"]
              }else{
                text_x[ind] <<- args[["label_text_x"]][y]
              }
              # Y COORD - STACKS/LIMITS
              if (.all_na(args[["label_text_y"]][y])) {
                text_y[ind] <<- y_coords[y]
              }
              # NEGATED POPULATION
            } else if (ind > nrow(gate_centers)) {
              # X COORD - RANGE CENTER
              if (.all_na(args[["label_text_x"]][y])) {
                # NO EVENTS
                if (.cyto_count(args[["pops"]][[y]]) == 0) {
                  # RANGE CENTER - PLOT LIMITS
                  text_x[ind] <<- mean(c(xmin, xmax))
                } else {
                  # RANGE
                  rng <- suppressMessages(
                    .cyto_range(args[["pops"]][[y]],
                      channels = args[["channels"]],
                      axes_limits = args[["axes_limits"]],
                      buffer = 0,
                      anchor = FALSE
                    )[, args[["channels"]]]
                  )
                  # UNIMODAL - 50% RANGE
                  if (abs(diff(rng)) <= 0.6 * xrng) {
                    text_x[ind] <<- quantile(rng, 0.5)
                    # UMULTIMODAL - 56% RANGE
                  } else {
                    text_x[ind] <<- quantile(rng, 0.56)
                  }
                }
                # X COORD MANUALLY SUPPLIED
              } else if (!.all_na(args[["label_text_x"]][y])) {
                text_x[ind] <<- args[["label_text_x"]][y]
              }
              # Y COORD - STACKS/LIMITS
              if (.all_na(args[["label_text_y"]][y])) {
                text_y[ind] <<- y_coords[y]
                # Y COORD MANUALLY SUPPLIED
              } else if (!.all_na(args[["label_text_y"]][y])) {
                text_y[ind] <<- args[["label_text_y"]][y]
              }
            }
            # NO GATE
          } else if (.all_na(args[["gate"]])) {
            # X COORD - RANGE CENTER
            if (.all_na(args[["label_text_x"]][y])) {
              # NO EVENTS
              if (.cyto_count(args[["pops"]][[y]]) == 0) {
                # RANGE CENTER - PLOT LIMITS
                text_x[ind] <<- mean(c(xmin, xmax))
              } else {
                # RANGE
                rng <- suppressMessages(
                  .cyto_range(args[["pops"]][[y]],
                    channels = args[["channels"]],
                    axes_limits = args[["axes_limits"]],
                    buffer = 0,
                    anchor = FALSE
                  )[, args[["channels"]]]
                )
                # UNIMODAL - 50% RANGE
                if (abs(diff(rng)) <= 0.6 * xrng) {
                  text_x[ind] <<- quantile(rng, 0.5)
                  # MULTIMODAL - 56% RANGE
                } else {
                  text_x[ind] <<- quantile(rng, 0.56)
                }
              }
              # X COORD MANUALLY SUPPLIED
            } else if (!.all_na(args[["label_text_x"]][y])) {
              text_x[ind] <<- args[["label_text_x"]][y]
            }
            # Y COORD - STACKS/LIMITS
            if (.all_na(args[["label_text_y"]][y])) {
              text_y[ind] <<- y_coords[y]
              # Y COORD MANUALLY SUPPLIED
            } else if (!.all_na(args[["label_text_y"]][y])) {
              text_y[ind] <<- args[["label_text_y"]][y]
            }
          }
          # 2D PLOT - MODE
        } else if (length(args[["channels"]]) == 2) {
          # GATE
          if (!.all_na(args[["gate"]])) {
            # GATED POPULATION
            if (ind <= nrow(gate_centers)) {
              # X COORD - GATE CENTER
              if(.all_na(args[["label_text_x"]][y])){
                text_x[ind] <<- gate_centers[ind, "x"]
              }else{
                text_x[ind] <<- args[["label_text_x"]][y]
              }
              # Y COORD - GATE CENTER
              if(.all_na(args[["label_text_y"]][y])){
                text_y[ind] <<- gate_centers[ind, "y"]
              }else{
                text_y[ind] <<- args[["label_text_y"]][y]
              }
              # NEGATED POPULATION
            } else if (ind > nrow(gate_centers)) {
              # X COORD - MODE/RANGE CENTER
              if (.all_na(args[["label_text_x"]][y])) {
                # NO EVENTS
                if (.cyto_count(args[["pops"]][[y]]) < 2) {
                  # RANGE CENTER
                  text_x[ind] <<- mean(c(xmin, xmax))
                } else {
                  # MODE
                  text_x[ind] <<- suppressMessages(
                    .cyto_mode(args[["pops"]][[y]],
                      channels = args[["channels"]][1],
                      density_smooth = args[["density_smooth"]]
                    )
                  )
                }
                # X COORD MANUALLY SUPPLIED
              } else if (!.all_na(args[["label_text_x"]][y])) {
                text_x[ind] <<- args[["label_text_x"]][y]
              }
              # Y COORD - MODE/RANGE CENTER
              if (.all_na(args[["label_text_y"]][y])) {
                # NO EVENTS
                if (.cyto_count(args[["pops"]][[y]]) == 0) {
                  # RANGE CENTER
                  text_y[ind] <<- mean(c(ymin, ymax))
                } else {
                  # MODE
                  text_y[ind] <<- suppressMessages(
                    .cyto_mode(args[["pops"]][[y]],
                      channels = args[["channels"]][2],
                      density_smooth = args[["density_smooth"]]
                    )
                  )
                }
                # Y COORD MANUALLY SUPPLIED
              } else if (!.all_na(args[["label_text_y"]][y])) {
                text_y[ind] <<- args[["label_text_y"]][y]
              }
            }
            # NO GATE
          } else if (.all_na(args[["gate"]])) {
            # X COORD - MODE/RANGE CENTER
            if (.all_na(args[["label_text_x"]][y])) {
              # NO EVENTS
              if (.cyto_count(args[["pops"]][[y]]) < 2) {
                # RANGE CENTER
                text_x[ind] <<- mean(c(xmin, xmax))
              } else {
                # MODE
                text_x[ind] <<- suppressMessages(
                  .cyto_mode(args[["pops"]][[y]],
                    channels = args[["channels"]][1],
                    density_smooth = args[["density_smooth"]]
                  )
                )
              }
              # X COORD SUPPLIED MANUALLY
            } else if (!.all_na(args[["label_text_x"]][y])) {
              text_x[ind] <<- args[["label_text_x"]][y]
            }
            # Y COORD - MODE
            if (.all_na(args[["label_text_y"]][y])) {
              # NO EVENTS
              if (.cyto_count(args[["pops"]][[y]]) < 2) {
                text_y[ind] <<- mean(c(ymin, ymax))
              } else {
                # MODE
                text_y[ind] <<- suppressMessages(
                  .cyto_mode(args[["pops"]][[y]],
                    channels = args[["channels"]][2],
                    density_smooth = args[["density_smooth"]]
                  )
                )
              }
              # Y COORD MANUALLY SUPPLIED
            } else if (!.all_na(args[["label_text_y"]][y])) {
              text_y[ind] <<- args[["label_text_y"]][y]
            }
          }
        }
        # NO LABEL
      } else if (.all_na(args[["label_text"]][y])) {
        text_x[ind] <<- NA
        text_y[ind] <<- NA
      }
    })
    # MATRIX
    text_xy <- matrix(c(text_x, text_y),
      ncol = 2,
      byrow = FALSE
    )
    colnames(text_xy) <- c("x", "y")
    return(text_xy)
  })
  label_text_xy <- do.call("rbind", label_text_xy)

  # UPDATE LABEL_TEXT_X & LABEL_TEXT_Y
  args[["label_text_x"]] <- unlist(label_text_xy[, "x"])
  args[["label_text_y"]] <- unlist(label_text_xy[, "y"])

  # OFFSET LABEL CO-ORDINATES --------------------------------------------------

  # LABEL DIMENSIONS
  label_dims <- lapply(seq_len(length(args[["label_text"]])), function(z) {
    # COMPUTE LABEL DIMENSIONS
    if (!.all_na(args[["label_text"]][z])) {
      .cyto_label_dims(
        label_text = args[["label_text"]][z],
        label_text_x = args[["label_text_x"]][z],
        label_text_y = args[["label_text_y"]][z],
        label_text_size = args[["label_text_size"]][z]
      )
    } else {
      matrix(rep(NA, 4),
        ncol = 2,
        dimnames = list(NULL, c("x", "y"))
      )
    }
  })
  
  # OFFSET BY LAYER
  if (length(args[["channels"]]) == 1 & args[["density_stack"]] != 0) {
    # SPLIT LABEL_DIMS BY LAYER
    label_dims <- split(
      label_dims,
      rep(seq_len(args[["SMP"]]), each = args[["NP"]])
    )
    # OFFSET PER LAYER
    lapply(seq_len(args[["SMP"]]), function(z) {
      # LABEL OVERLAP
      if (.cyto_label_overlap(label_dims[[z]])) {
        # LABEL HEIGHT - OFFSETTING
        label_height <- max(LAPPLY(label_dims[[z]], function(y) {
          max(y[, "y"]) - min(y[, "y"])
        }), na.rm = TRUE)
        # LABEL HEIGHT BUFFERING
        label_height <- 1.18 * label_height
        # Y COORDS TO OFFSET
        text_y <- args[["label_text_y"]][label_ind[[z]]]
        # OFFSET Y CO-ORDINATES (EXCLUDE NA)
        args[["label_text_y"]][label_ind[[z]]][!is.na(text_y)] <<- tryCatch({
          .suppress_all_messages(
            .spread.labels(text_y[!is.na(text_y)],
              mindiff = label_height,
              min = ymin,
              max = ymax
            )
          )
        }, error = function(e) {
          return(text_y[!is.na(text_y)])
        })
      }
    })
    # OFFSET ALL LABELS
  } else {
    # LABEL OVERLAP
    if (.cyto_label_overlap(label_dims)) {
      # LABEL HEIGHT - OFFSETTING
      label_height <- max(LAPPLY(label_dims, function(y) {
        max(y[, "y"]) - min(y[, "y"])
      }), na.rm = TRUE)
      # LABEL HEIGHT BUFFERING
      label_height <- 1.18 * label_height
      # OFFSET Y CO-ORDINATES (EXCLUDE NA)
      args[["label_text_y"]][!is.na(args[["label_text_y"]])] <-
        .suppress_all_messages(
          .spread.labels(args[["label_text_y"]][!is.na(args[["label_text_y"]])],
            mindiff = label_height,
            min = ymin,
            max = ymax
          )
        )
    }
  }

  # RETURN COMPUTED LABEL CO-ORDINATES -----------------------------------------

  # LABEL CO-ORDINATE MATRIX
  label_text_xy <- matrix(c(args[["label_text_x"]], args[["label_text_y"]]),
    ncol = 2,
    byrow = FALSE
  )
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
}

## .CYTO_LABEL_DIMS ------------------------------------------------------------

#' Compute label dimensions
#'
#' Labels should already contain statistic as well. Co-ordiante for each label
#' must already be computed.
#'
#' @importFrom graphics strwidth strheight
#'
#' @return upper left and bottom right x and y co-ordinates of labels
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_label_dims <- function(label_text,
                             label_text_x,
                             label_text_y = NA,
                             label_text_size = 1,
                             xpad = 1.2,
                             ypad = 1.2,
                             adj = 0.5) {

  # GRAPHICAL PARAMETERS -------------------------------------------------------

  # RESET CEX & XPD
  old_pars <- par(c("cex", "xpd"))
  on.exit(par(old_pars))

  # SET CEX & XPD
  par(cex = label_text_size)
  par(xpd = TRUE)

  # LABEL DIMENSIONS -----------------------------------------------------------

  # BOX ADJUSTMENT
  box_adj <- adj + (xpad - 1) * label_text_size * (0.5 - adj)

  # BOX DIMENSIONS
  lwidths <- strwidth(label_text)
  rwidths <- lwidths * (1 - box_adj)
  lwidths <- lwidths * box_adj
  bheights <- theights <- strheight(label_text) * 0.5

  # BOX X COORDS
  xr <- label_text_x - lwidths * xpad
  xl <- label_text_x + lwidths * xpad

  # BOX Y COORDS
  yb <- label_text_y - bheights * ypad
  yt <- label_text_y + theights * ypad

  # LABEL DIMENSIONS MATRIX ----------------------------------------------------

  # MATRIX - TOP LEFT THEN BOTTOM RIGHT
  coords <- matrix(c(
    min(c(xl, xr)),
    max(c(yb, yt)),
    max(c(xl, xr)),
    min(c(yb, yt))
  ),
  ncol = 2,
  byrow = TRUE
  )
  colnames(coords) <- c("x", "y")

  # RETURN LABEL DIMENSIONS ----------------------------------------------------
  return(coords)
}

## .CYTO_LABEL_OVERLAP ---------------------------------------------------------

#' Check if any labels are overlapping.
#'
#' @param x list of label dimensions.
#'
#' @return TRUE or FALSE based on whether any overlapping labels are found.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_label_overlap <- function(x) {

  # For each rectangle in x
  overlaps <- LAPPLY(seq_len(length(x)), function(y) {

    # Check if other rectangles overlap
    LAPPLY(seq_len(length(x))[-y], function(z) {

      # Co-ordinates of reference label
      x1 <- x[[y]][, "x"]
      y1 <- x[[y]][, "y"]

      # Co-ordinates of comparison label
      x2 <- x[[z]][, "x"]
      y2 <- x[[z]][, "y"]

      # MISSING LABELS - NO OVERLAP
      if (any(is.na(c(x1, x2, y1, y2)))) {
        return(FALSE)
      }

      # X co-ordinates are overlapping
      if (min(x2) >= min(x1) & min(x2) <= max(x1) |
        max(x2) >= min(x1) & max(x2) <= max(x1)) {
        # Y co-ordinates are also overlapping
        if (min(y2) >= min(y1) & min(y2) <= max(y1) |
          max(y2) >= min(y1) & max(y2) <= max(y1)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }

      # Non-overlapping x and y co-ordinates
      return(FALSE)
    })
  })

  # RETURN TRUE OR FALSE
  return(any(overlaps))
}
