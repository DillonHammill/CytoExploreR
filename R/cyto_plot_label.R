# CYTO_PLOT_LABEL --------------------------------------------------------------

#' Add labels to cyto_plot
#'
#' Interactively label existing plots with boxed labels containing text and
#' statistics.
#'
#' @param x object of class \code{"flowFrame"}.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param gate object of class
#'   \code{\link[flowCore:rectangleGate-class]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate-class]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}},
#'   \code{\link[flowCore:quadGate-class]{quadGate}}, \code{"list"} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param label_text character string to include in the label above the
#'   statistic (e.g. population name(s)).
#' @param label_stat indicates the type of statistic to include in the label,
#'   can be either code{"count"}, \code{"median"}, \code{"mean"}, \code{"mode"},
#'   \code{"geo mean"} or \code{"CV"}. Only count and percent statistics are
#'   supported for 2D plots.
#' @param label_text_x vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param label_text_y vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param trans object of class
#'   \code{\link[flowWorkspace]{transformerList}} which was used to
#'   transform the fluorescent channels of the supplied flowFrame.
#' @param negate logical indicating whether a label should be included for the
#'   neagted population (i.e. events outside the gate). Set to FALSE by default
#'   to only calculate statistics for events within the gate.
#' @param label_text_font integer [1,2,3,4] passed to \code{text} to alter the
#'   font, set to \code{2} by default for a bold font.
#' @param label_text_size numeric character expansion used to control the size
#'   of the text in the labels, set to \code{0.8} by default. See \code{?text}
#'   for details.
#' @param label_text_col specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param label_fill fill colour to use for labels, set to "white" by default.
#' @param label_fill_alpha numeric [0,1] controls the transparency of the fill
#'   colour, set to \code{0.6} by default.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param density_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust kernel density for mode
#'   calculation.
#'
#' @return add a boxed text label to an existing cyto_plot.
#'
#' @importFrom flowCore Subset
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
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(fs))
#' gs <- transform(gs, trans)
#'
#' # Gate using gate_draw
#' gating(Activation_gatingTemplate, gs)
#'
#' # Plot
#' cyto_plot(gs[[32]],
#'   parent = "CD4 T Cells",
#'   channels = c("CD69")
#' )
#'
#' # Label - median fluorescent intensity
#' cyto_plot_label(getData(gs, "CD4 T Cells")[[32]],
#'   channels = "CD69",
#'   label_stat = "median",
#'   label_text = "MedFI",
#'   trans = trans,
#'   label_text_x = 3,
#'   label_text_y = 50
#' )
#'
#' # Label - no statistic
#' cyto_plot_label(
#'   label_text = "Activated \n Cells",
#'   label_text_x = 3,
#'   label_text_y = 25,
#'   label_fill = "red"
#' )
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_plot_label <- function(x,
                            channels = NULL,
                            trans = NA,
                            display = 1,
                            gate = NA,
                            negate = FALSE,
                            label_text = NA,
                            label_stat = NA,
                            label_text_x = NA,
                            label_text_y = NA,
                            label_text_font = 2,
                            label_text_size = 0.8,
                            label_text_col = "black",
                            label_fill = "white",
                            label_fill_alpha = 0.6,
                            density_smooth = 0.6){
  
  # CHECKS ---------------------------------------------------------------------
  
  # FLOWFRAME
  if(!inherits(x, "flowFrame")){
    stop("'x' must be a flowFrame object.")
  }
  
  # SAMPLING
  if(display != 1){
    x <- cyto_sample(x, display = display, seed = 56)
  }
  
  # STATISTIC- CHECK & SUPPORT
  if(!.all_na(label_stat)){
    # CONVERT TO VALID STATISTIC
    label_stat <- LAPPLY(label_stat, function(z){
      .cyto_stat_check(z)
    })
    # STATISTICS 2D PLOT
    if(length(channels == 2)){
      # COUNT ONLY - 2D PLOT WITHOUT GATE
      if(.all_na(gate)){
        if(!.all(label_stat == "count")){
          stop("Only 'count' is supported for 2D plots without gates.")
        }
      # COUNT & FREQ - 2D PLOT
      }else{
        if(!all(label_stat %in% c("count","freq"))){
          stop("Only 'count' and 'freq' are supported in 2D plots.")
        }
      }
    }
    # FREQ - MUST HAVE GATE
    if(any(label_stat == "freq") & .all_na(gate)){
      stop("Supply gate objects to 'gate' to calculate frequency.")
    }
  }

  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # PREPARE GATES --------------------------------------------------------------
  
  # LIST OF GATE OBJECTS
  if(!.all_na(gate)){
    if (!inherits(gate, "list")) {
      if (inherits(gate, "rectangleGate") |
          inherits(gate, "polygonGate") |
          inherits(gate, "ellipsoidGate") |
          inherits(gate, "quadGate")) {
        gate <- list(gate)
      } else if (inherits(gate, "filters")) {
        gate <- unlist(gate)
      }
    } else {
      gate <- unlist(gate)
    }
  }
  
  # CONVERT GATE - DIMENSIONS
  if(!.all_na(gate)){
    gate <- cyto_gate_convert(gate, channels = channels)
  }
  
  # POPULATIONS TO LABEL -------------------------------------------------------
  
  # APPLY GATE
  if(!.all_na(label_stat)){
    # GATE SUPPLIED
    if(!.all_na(gate)){
      # QUADGATE - NEGATE NOT SUPPORTED
      if(all(LAPPLY(gate, "class") == "quadGate")){
        # GATED POPULATIONS X4
        pops <- split(x, gate[[1]])
      # LIST OF GATE OBJECTS - NEGATE SUPPORTED
      }else{
        # GATED POPULATIONS
        pops <- lapply(gate, function(z){
          Subset(x, z)
        })
        # NEGATED POPULATION
        if(negate == TRUE){
          gate_filter <- do.call("|", gate)
          pops <- c(pops,
                    list(split(x, gate_filter)[[2]]))
        }
      }
    # NO GATE
    }else{
      pops <- list(x)
    }
  }
  
  # POPULATIONS
  PNS <- length(pops)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # REPEAT ARGUMENTS - GATE
  if(!.all_na(gate)){
   
    # REPEAT LABEL ARGS PER POPULATION
    args[["label_text"]] <- rep(c(label_text, rep(NA, PNS)), 
                              length.out = PNS)
    label_args <- formalArgs(cyto_plot_label)
    label_args <- label_args[grepl("label_", label_args)]
    label_args <- label_args[-match("label_text", label_args)]
    lapply(label_args, function(z){
      args[[z]] <<- rep(args[[z]], length.out = PNS)
    })
    
  # REPEAT ARGUMENTS - WITHOUT GATE
  }else{
    
    # REPEAT LABEL ARGS PER LABEL_TEXT
    label_args <- formalArgs(cyto_plot_label)
    label_args <- label_args[grepl("label_", label_args)]
    label_args <- label_args[-match("label_text", label_args)]
    lapply(label_args, function(z){
      args[[z]] <<- rep(args[[z]], length.out = length(label_text))
    })
    
  }

  # UPDATE ARGUMENTS
  .args_update(args)
  
  # STATISTICS PER POPULATION --------------------------------------------------
  
  # CALCULATE STATISTICS
  if (!.all_na(label_stat)) {
    st <- LAPPLY(pops, function(z) {
      if (label_stat == "freq") {
        sts <- .cyto_count(z) / .cyto_count(x) * 100
        sts <- .round(sts, 2)
        sts <- paste(sts, "%")
      }else if(label_stat == "CV"){
        sts <- .cyto_CV(z,
                        channels = channels,
                        trans = trans)
        sts <- .round(sts, 2)
        sts <- paste(sts, "%")
      } else {
        sts <- cyto_stats_compute(z,
                                  channels = channels,
                                  trans = trans,
                                  stat = label_stat,
                                  format = "long",
                                  density_smooth = density_smooth
        )
        sts <- LAPPLY(sts[, ncol(sts)], function(x){.round(x, 2)})
      }
      return(sts)
    })
  } else {
    st <- NA
  }
  
  # LABEL CONSTRUCTION ---------------------------------------------------------
  
  # MANUAL CO-ORDINATES - SUPPLIED OR SELECTED
  label_text_xy <- mapply(
    function(label_text,
             st,
             label_text_x,
             label_text_y,
             label_text_font,
             label_text_size,
             label_text_col,
             label_fill,
             label_fill_alpha) {
      
      # Co-ordinates can be interactively selected
      if (any(.all_na(c(label_text_x, label_text_y)))) {
        message("Select a location on the plot the position the label.")
        label_text_xy <- locator(n = 1)
        label_text_x <- label_text_xy[[1]]
        label_text_y <- label_text_xy[[2]]
      }
      
      # Add labels to plot
      if (!.all_na(label_text) & !.all_na(label_stat)) {
        .boxed.labels(
          x = label_text_x,
          y = label_text_y,
          labels = paste(label_text, st, sep = "\n"),
          border = FALSE,
          font = label_text_font,
          cex = label_text_size,
          col = label_text_col,
          bg = label_fill,
          alpha.bg = label_fill_alpha
        )
      } else if (!.all_na(label_text) & .all_na(label_stat)) {
        .boxed.labels(
          x = label_text_x,
          y = label_text_y,
          labels = label_text,
          border = FALSE,
          font = label_text_font,
          cex = label_text_size,
          col = label_text_col,
          bg = label_fill,
          alpha.bg = label_fill_alpha
        )
      } else if (.all_na(label_text) & !.all_na(label_stat)) {
        .boxed.labels(
          x = label_text_x,
          y = label_text_y,
          labels = st,
          border = FALSE,
          font = label_text_font,
          cex = label_text_size,
          col = label_text_col,
          bg = box_fill,
          alpha.bg = box_alpha
        )
      }
      
      # RETURN LABEL CO-ORDINATES
      return(c(label_text_X, label_text_y))
      
    },
    label_text,
    st,
    label_text_x,
    label_text_y,
    label_text_font,
    label_text_size,
    label_text_col,
    label_fill,
    label_fill_alpha,
    SIMPLIFY = FALSE
  )
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL CO-ORDINATES IN MATRIX
  label_text_xy <- do.call("rbind", label_text_xy)
  colnames(label_text_xy) <- c("x","y")
  label_text_xy <- as.matrix(label_text_xy)
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)
  
}

# .CYTO_PLOT_LABEL -------------------------------------------------------------

# This internal function handles all the labelling within cyto_plot. It is
# different to cyto_plot_label which labels by layer and expects user input to
# position labels. The exported cyto_plot_label function is to be used to label
# existing plots, and since it is impossible to know the y_max per layer for 1D
# plots - automtic label locations cannot be computed. Here automatic label
# positions can be calculated given internal access to all of cyto_plot
# arguments.

#' Internal cyto_plot_label function - used by cyto_plot
#' @param x named list of cyto_plot arguments
#' @importFrom graphics par
#' @noRd
.cyto_plot_label <- function(x){
  
  # INHERIT ARGUMENTS ----------------------------------------------------------

  .args_update(x)
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng/1.04)/2   # DEFAULT 104% * ylim
  ymin <- ymin + 0.5 * ypad # 1% BUFFER
  ymax <- ymax - 0.5 * ypad # 1% BUFFER
  yrng <- ymax - ymin
  
  # GENERAL --------------------------------------------------------------------
  
  # SAMPLES
  SMP  <- length(fr_list)
  
  # GATES
  if(!.all_na(gate)){
    GTS <- length(gate)
  }
  
  # NEGATE FILTER
  if(negate == TRUE){
    negate_filter <- do.call("|", gate)
  }
  
  # LABEL_STAT
  label_stat <- LAPPLY(label_stat, function(z){
    if(.all_na(label_stat[z])){
      return(NA)
    }else{
      .cyto_stat_check(stat[z])
    }
  })

  # POPULATIONS TO LABEL -------------------------------------------------------
  
  # CYTO_PLOT ALWAYS WANTS TO LABEL EACH LAYER!
  # 1D PLOT - UNSTACKED - UPPER LAYERS SWITCHED OFF (LABEL_TEXT & LABEL_STAT)
  # 1D PLOT - STACKED - EACH LAYER IS LABELLED BY DEFAULT
  # 2D PLOT - UPPER LAYERS ARE SWITCHED OFF (LABEL_TEXT = NA & LABEL_STAT = NA)
  
  # NUMBER OF LABELS PER LAYER
  if(!.all_na(gate)){
    L <- length(gate)
    if(negate == TRUE){
      L <- L + 1
    }
  }else if(.all_na(gate)){
    L <- 1
  }
  
  # TOTAL NUMBER OF LABELS
  TL <- length(fr_list) * L
  
  # DEFAULT LABEL CO-ORDINATES - DENSITY ---------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # Y_MAX PER LAYER
    y_max <- strsplit(names(y)[1], "-")[2]
    # NO STACKING
    if(density_stack == 0){
      y_coords <- rep(0.5 * yrng, SMP)
    # CALCULATE Y STACKING
    }else if(density_stack != 0){
      # STACK LEVELS
      stk <- lapply(seq(0, SMP, 1), function(z){
        density_stack * y_max * z
      })
      # DEFAULT LABEL Y CO-ORDINATES
      y_coords <- lapply(seq_len(SMP), function(z){
        # HALFWAY BETWEEN HORIZONTAL LINES
        0.5 * stk[(z + 1)]
      })
      # REPEAT LBLS TIMES (Y_COORDS/LAYER)
      y_coords <- rep(y_coords, each = LBLS)
    }
  }  
  
  # SUPPORTED STATISTICS -------------------------------------------------------
  
  # 2D PLOTS - COUNT & FREQ ONLY
  if(!.all_na(label_stat)){
    if(length(channels) == 2 & !all(label_stat %in% c("count","freq"))){
      message("Only frequency and count statistics are supported in 2D plots")
      # REPLACE WITH NA
      label_stat[!label_stat %in% c("count", "freq")] <- NA
    }
  }
  
  # LABEL_STAT & LABEL COORDS --------------------------------------------------
  
  # LABEL_STAT & LABEL_COORDS
  lapply(seq_len(SMP), function(z){
    # PARENTAL FLOWFRAME
    fr <- fr_list[[z]]
    # PER LABEL
    lapply(seq_len(L), function(y){
      # INDEX IN VECTOR
      ind <- z * y
      # GATING REQUIRED FOR LABEL_STAT
      if(!.all_na(label_stat[ind])){
        # GATED POPULATION
        if(!.all_na(gate)){
          # GATED POPULATION
          if(L <= GTS){
            pop <- Subset(fr, gate[[L]])
          }else{
            pop <- split(fr, negate_filter)[[2]]
          }
        # NO GATE
        }else{
          pop <-fr
        }
        # FREQUENCY STATISTIC
        if(grepl("freq", label_stat[ind], ignore.case = TRUE)){
          label_stat[ind] <<- .cyto_count(pop) / .cyto_count(fr) * 100
          label_stat[ind] <<- paste(.round(label_stat[ind], 2), "%")
        # CV STATISTIC
        }else if(grepl("CV", label_stat[ind], ignore.case = TRUE)){
          label_stat[ind] <<- .cyto_CV(pop,
                                       channels = channels,
                                       trans = axes_trans)
          label_stat[ind] <<- paste(.round(label_stat[ind], 2), "%")
        # OTHER STATISTIC
        }else{
          label_stat[ind] <<- cyto_stats_compute(pop,
                                                 channels = channels,
                                                 trans = axes_trans,
                                                 stat = label_stat[ind],
                                                 format = "long",
                                                 density_smooth = density_smooth)
        }
      }
      # LABEL CO-ORDINATES
      if(label_position == "auto"){
        # 1D PLOT
        if(length(channels) == 1){
          # GATE
          if(!.all_na(gate)){
            # GATE CENTER - GATED POPULATION
            if(L <= GTS){
              # X COORD - GATE CENTER
              label_coords <- .cyto_plot_label_center(
                gate[[L]],
                channels = channels,
                text_X = label_text_x[ind],
                text_y = label_text_y[ind]
              )
              label_text_x[ind] <<- label_coords[, "x"]
              # Y COORD
              if(.all_na(label_text_y[ind])){
                label_text_y[ind] <<- y_coords[ind]
              }
            # MODE - NEGATED POPULATION
            }else{
              # X MODE
              if(.all_na(label_text_x[ind])){
                label_text_x[ind] <<- suppressMessages(
                  .cyto_mode(pop,
                             channels = channels,
                             density_smooth = density_smooth)
                )
              }
              # Y COORD
              if(.all_na(label_text_y[ind])){
                label_text_y[ind] <<- y_coords[ind]
              }
            }
          # NO GATE
          }else if(.all_na(gate)){
            # X COORD - MODE
            if(.all_na(label_text_x[ind])){
              label_text_x[ind] <<- suppressMessages(
                .cyto_mode(pop,
                           channels = channels,
                           density_smooth = density_smooth)
              )
            }
            # Y COORD
            if(.all_na(label_text_y[ind])){
              label_text_y[ind] <<- y_coords[ind]
            }
          }
        # 2D PLOT 
        }else if(length(channels) == 2){
          # GATE
          if(!.all_na(gate)){
            # GATED POPULATION - GATE CENTER
            if(L <= GTS){
              label_coords <- .cyto_plot_label_center(
                gate[[L]],
                channels = channels,
                text_X = label_text_x[ind],
                text_y = label_text_y[ind]
              )
              label_text_x[ind] <<- label_coords[, "x"]
              label_text_y[ind] <<- label_coords[, "y"]
            # NEGATED POPULATION - MODE
            }else{
              # X COORD - MODE
              if(.all_na(label_text_x[ind])){
                label_text_x[ind] <<- suppressMessages(
                  .cyto_mode(pop,
                             channels = channels[1],
                             density_smooth = density_smooth)
                )
              }
              # Y COORD - MODE
              if(.all_na(label_text_y[ind])){
                label_text_y[ind] <<- suppressMessages(
                  .cyto_mode(pop,
                             channels = channels[2],
                             density_smooth = density_smooth)
                )
              }
            }
          # NO GATE
          }else if(.all_na(gate)){
            # X COORD - MODE
            if(.all_na(label_text_x[ind])){
              label_text_x[ind] <<- suppressMessages(
                .cyto_mode(pop,
                           channels = channels[1],
                           density_smooth = density_smooth)
              )
            }
            # Y COORD - MODE
            if(.all_na(label_text_y[ind])){
              label_text_y[ind] <<- suppressMessages(
                .cyto_mode(pop,
                           channels = channels[2],
                           density_smooth = density_smooth)
              )
            }
          }
        }
      }
    })
    # UPDATE LABEL_STAT
    label_stat[seq_len(L) * z] <<- label_stat[seq_len(L) * z]
    # UPDATE LABEL CO-ORDINATES
    label_text_x[seq_len(L) * z] <<- label_text_x[seq_len(L) * z]
    label_text_y[seq_len(L) * z] <<- label_text_y[seq_len(L) * z]
  })
  
  # PREPARE LABEL_TEXT ---------------------------------------------------------
  
  # MERGE LABEL_TEXT & LABEL_STAT
  lapply(seq_len(label_text), function(z){
    # NO LABEL
    if(.all_na(label_text[z]) & .all_na(label_stat[z])){
      label_text[z] <<- NA
    # NO LABEL_TEXT 
    }else if(.all_na(label_text[z]) & !.all_na(label_stat[z])){
      label_text[z] <<- paste(label_text[z], label_stat[z], sep = "\n")
    # NO LABEL_STAT
    }else if(!.all_na(label_text[z]) & .all_na(label_stat[z])){
      label_text[z] <<- label_stat[z]
    }
  })
    
  # OFFSET LABEL CO-ORDINATES --------------------------------------------------
  
  # OFFSET LABELS IN EACH LAYER
  # 2D OFFSET ALL
  # 1D UNSTACKED OFFSET ALL
  
  # OFFSET - LABEL_POSITION == "auto" ONLY
  if(label_position == "auto"){
    # OFFSET EACH LAYER - STACKED 1D PLOTS
    if(length(channels) == 1 & density_stack != 0){
      # UPDATE COORDS PER LAYER
      lapply(seq_len(SMP), function(z){
        # CALCULATE OFFSET LABEL LOCATIONS- MULTIPLE LABELS ONLY
        if(L > 1){
          # INDICES
          ind <- seq_len(L) * z
          # OFFSET LABEL COORDS
          label_text_xy <- .cyto_plot_label_offset(text = label_text[ind],
                                                   text_x = label_text_x[ind],
                                                   text_y = label_text_y[ind],
                                                   text_size = label_text_size[ind])
          # UPDATE LABEL CO-ORDINATES
          label_text_x[ind] <<- label_text_xy[, "x"]
          label_text_y[ind] <<- label_text_xy[, "y"]
        }
      })
    # OFFSET ALL LABELS - UNSTACKED 1D & 2D PLOTS  
    }else{
      # CALCULATE OFFSET LABEL LOCATIONS - MULTIPLE LABELS ONLY
      if(TL > 1){
      label_text_xy <- .cyto_plot_label_offset(text = label_text,
                                               text_x = label_text_x,
                                               text_y = label_text_y,
                                               text_size = label_text_size)
      }
    }
  }
  
  # CYTO_PLOT_LABEL ------------------------------------------------------------
  
  # CALL TO CYTO_PLOT_LABEL - ADDS LABEL ONLY - NO CALCULATIONS
  label_text_xy <- mapply(function(label_text,
                                   label_text_x,
                                   label_text_y,
                                   label_text_font,
                                   label_text_size,
                                   label_text_col,
                                   label_fill,
                                   label_fill_alpha){
    
    # SELECT CO-ORDINATES FOR NA 
    if(!.all_na(label_text)){
      cyto_plot_label(fr_list[[1]], # IRRELEVANT - NO CALCULATIONS
                      channels = channels,
                      trans = axes_trans,
                      gate = NA,
                      negate = FALSE,
                      label_text = label_text,
                      label_stat = NA,
                      label_text_x = label_text_x,
                      label_text_y = label_text_y,
                      label_text_font = label_text_font,
                      label_text_size = label_text_size,
                      label_text_col = label_text_col,
                      label_fill = label_fill,
                      label_fill_alpha = label_fill_alpha,
                      density_smooth = density_smooth)
    }

    
  },
  label_text,
  label_text_x,
  label_text_y,
  label_text_font,
  label_text_size,
  label_text_col,
  label_fill,
  label_fill_alpha,
  SIMPLIFY = FALSE)
  
  # LABEL CO-ORDINATES MATRIX
  label_text_xy <- do.call("rbind", label_text_xy)
  colnames(label_text_xy) <- c("x","y")
  label_text_xy <- as.matrix(label_text_xy)
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  return(label_text_xy)
  
}
