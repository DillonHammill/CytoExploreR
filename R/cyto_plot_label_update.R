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
#' cyto_plot_label2(getData(gs, "CD4 T Cells")[[32]],
#'   channels = "CD69",
#'   label_stat = "median",
#'   label_text = "MedFI",
#'   trans = trans,
#'   label_text_x = 3,
#'   label_text_y = 50
#' )
#'
#' # Label - no statistic
#' cyto_plot_label2(
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
  xmin <- lims[1]
  xmax <- lims[2]
  ymin <- lims[3]
  ymax <- lims[4]
  ypad <- ((ymax - ymin) - (ymax - ymin)/1.04)/2   # DEFAULT 104% * ylim
  yrng <- ymax - ypad
  
  # GENERAL --------------------------------------------------------------------
  
  # SAMPLES
  SMP  <- length(fr_list)
  
  # GATES
  if(!.all_na(gate)){
    GTS <- length(gate)
  }

  # POPULATIONS TO LABEL -------------------------------------------------------
  
  # GATE
  if(!.all_na(gate)){
    LBLS <- length(gate)
    if(negate == TRUE){
      LBLS <- LBLS + 1
    }
  # NO GATE  
  }else{
    # NUMBER OF UNIQUE LABELS
    LBLS <- length(label_text)
  }
  
  # PREPARE LIST OF POPULATION LISTS -------------------------------------------

  # 1D PLOT
  if(length(channels) == 1){
    # NO STACKING - LIST OF POPULATIONS
    if(density_stack == 0){
      # GATE
      if(!.all_na(gate)){
        # GATED POPULATIONS
        pops <- lapply(gate, function(z){
            Subset(fr_list[[1]], z)
        })
        # NEGATED POPULATION
        if(negate == TRUE){
          gate_filter <- do.call("|", gate)
          pops <- c(pops, split(fr_list[[1]], gate_filter)[2])
        }
        # LIST OF POPULATION LISTS
        pops <- list(pops)
      # NO GATE  
      }else{
        # SAME POPULATION FOR EACH LABEL
        pops <- rep(fr_list[1], LBLS)
        # LIST OF POPULATION LISTS
        pops <- list(pops)
      }
    # STACKING - LIST OF POPULATION LISTS  
    }else if(density_stack != 0){
      # GATE 
      if(!.all_na(gate)){
        # GATE EACH LAYER
        pops <- lapply(fr_list, function(z){
          # GATED POPULATIONS
          pops_layer <- lapply(gate, function(y){
            Subset(z, y)
          })
          # NEGATED POPULATION
          if(negate == TRUE){
            gate_filter <- do.call("|", gate)
            pops_layer <- c(pops_layer,
                            split(z, gate_filter)[2])
          }
          return(pops_layer)
        })
      # NO GATE
      }else if(.all_na(gate)){
        # SAME POPULATION FOR EACH LABEL
        pops <- rep(fr_list[1], LBLS)
        # LIST OF POPULATION LISTS
        pops <- list(pops)
      }
    }
  # 2D PLOT - BASE LAYER ONLY
  }else if(length(channels) == 2){
    # GATE
    if(!.all_na(gate)){
      # GATED POPULATIONS
      pops <- lapply(gate, function(z){
        if(is(z, "quadGate")){
          split(fr_list[[1]], z)
        }else{
          Subset(fr_list[[1]], z)
        }
      })
      # NEGATED POPULATION
      if(negate == TRUE & !any(LAPPLY(gate, "is") == "quadGate")){
        gate_filter <- do.call("|", gate)
        pops <- c(pops, split(fr_list[[1]], gate_filter)[2])
      }
      # LIST OF POPULATION LISTS
      pops <- list(pops)
      # NO GATE  
    }else{
      # SAME POPULATION FOR EACH LABEL
      pops <- rep(fr_list[1], LBLS)
      # LIST OF POPULATION LISTS
      pops <- list(pops)
    }
  }
  
  # DEFAULT LABEL CO-ORDINATES - DENSITY ---------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # Y_MAX PER LAYER
    y_max <- strsplit(names(y)[1], "-")[2]
    # NO STACKING
    if(density_stack == 0){
      y_coords <- 0.5 * yrng
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
  
  # AUTOMATIC LABEL POSITIONING ------------------------------------------------
  
  # GATED POPULATIONS - GATE CENTERS
  # NEGATED POPULATIONS - MODE
  # POPULATIONS - MODE
  
  # CALCULATE CO-ORDINATES FOR LABELS
  if(label_position == "auto"){
    # 1D PLOT
    if(length(channels) == 1){
      # NO STACKING - SINGLE LABEL/GATE + NEGATED LABEL (BASE LAYER ONLY)
      if(density_stack == 0){
        # GATE SUPPLIED - USE GATE CENTERS
        if(!.all_na(gate)){
          # GATE CENTERS
          label_coords <- .cyto_plot_label_center(
            gate,
            channels = channels,
            text_X = label_text_x[seq_len(GTS)],
            text_y = label_text_y[seq_len(GTS)]
            )
          # NEGATE
          if(negate == TRUE){
            # X COORD - MODE
            if(.all_na(label_text_x[GTS + 1])){
              label_text_x[GTS + 1] <- suppressMessages(
                .cyto_mode(pops[[1]][[length(pops)]],
                           channels = channels[1],
                           density_smooth = density_smooth)
              )
            }
            # Y COORD - 50%
            if(.all_na(label_text_y[GTS + 1])){
              label_text_y[GTS + 1] <- y_coords
            }
          }
          # ADD NEGATE CO-ORDINATES TO LABEL_COORDS
          label_coords <- rbind(label_coords, c(label_text_x[GTS + 1],
                                                label_text_y[GTS + 1]))
        # NO GATE - MODE
        }else if(.all_na(gate)){
          # X COORD - MODE (OFFSET LATER)
          lapply(seq_len(label_text_X), function(z){
            # X COORD - MODE
            if(.all_na(label_text_x[z])){
              label_text_x[z] <<- suppressMessages(
                .cyto_mode(fr_list[[1]],
                           channels = channels[1],
                           density_smooth = density_smooth)
              )
            }
          })
          # Y COORD - 50%
          lapply(seq_len(label_text_y), function(z){
            # Y COORD - 50%
            if(.all_na(label_text_y[z])){
              label_text_y[z] <<- y_coords
            }
          })
          # LABEL_COORDS
          label_coords <- cbind(label_text_x, label_text_y)
          colnames(label_coords) <- c("x","y")
        }
      # STACKING - LABEL/GATE + NEGATED LABEL PER LAYER (ALL LAYERS)
      }else if(density_stack != 0){
        # GATE SUPPLIED
        if(!.all_na(gate)){
          # X COORDS - GATE CENTERS
          label_coords <- .cyto_plot_label_center(
            gate,
            channels = channels,
            text_X = label_text_x[seq_len(GTS)],
            text_y = label_text_y[seq_len(GTS)]
          )
          # LABEL_TEXT_X
          label_text_x[seq_len(GTS)] <- label_coords[,"x"]
          # NEGATE
          if(negate == TRUE){
            # X COORD - CALCULATE MODE FOR MERGED NEGATED POPULATIONS
            if(.all_na(label_text_x[GTS + 1])){
              # NEGATED POPULATIONS
              NP <- lapply(seq_len(length(pops)), function(z){
                pops[[z]][[length(pops[[z]])]]
              })
              NP <- as(flowSet(NP), "flowFrame")
              # NEGATE X COORD
              label_text_x[GTS + 1] <- suppressMessages(
                .cyto_mode(NP,
                           channels = channels,
                           density_smooth = density_smooth)
              )
            }
          }
          # Y COORDS
          label_text_y <- y_coords
          # REPEAT LABEL_TEXT_X PER LAYER
          label_text_x <- rep(label_text_x, times = SMP)
          # LABEL_COORDS MATRIX
          label_coords <- cbind(label_text_x, label_text_y)
          colnames(label_coords) <- c("x","y")
        # NO GATE SUPPLIED
        }else if(.all_na(gate)){
          # X COORDS - MODE
          lapply(seq_len(label_text_x), function(z){
            if(.all_na(label_text_x[z])){
              label_text_x[z] <<- .cyto_mode(fr_negated,
                                             channels = channels[1],
                                             density_smooth = density_smooth)
            }
          })
          # Y COORDS - STACKS
          label_text_y <- y_coords
          # LABEL_COORDS
          label_coords <- cbind(label_text_x, label_text_y)
          colnames(label_coords) <- c("x","y")
        }
      }
    # 2D PLOT - BASE LAYER ONLY 
    }else if(length(channels) == 2){
      # GATE SUPPLIED - USE GATE CENTERS
      if(!.all_na(gate)){
        # GATE CENTERS - CONSIDER MANUALLY SUPPLIED CO-ORDINATES
        label_coords <- .cyto_plot_label_center(gate,
                                              channels = channels,
                                              text_X = label_text_x[seq_len(GTS)],
                                              text_y = label_text_y[seq_len(GTS)])
        # NEGATE LABEL - MODE IN 2D (CANNOT NEGATE QUADGATES)
        if(negate == TRUE & !any(LAPPLY(gate, "is") == "quadGate")){
          # X COORD - MODE
          if(.all_na(label_text_x[SMP + 1])){
            label_text_x[GTS + 1] <- suppressMessages(
              .cyto_mode(fr_negated,
                         channels = channels[1],
                         density_smooth = density_smooth)
            )
          }
          # Y COORD - MODE
          if(.all_na(label_text_y[GTS + 1])){
            label_text_y[GTS + 1] <- suppressMessages(
              .cyto_mode(
                fr_negated,
                channels = channels[2],
                density_smooth = density_smooth)
            )
          }
        }
        # ADD NEGATE CO-ORDINATES TO LABEL_COORDS
        label_coords <- rbind(label_coords, c(label_text_x[GTS + 1],
                                              label_text_y[GTS + 1]))
      # NO GATE - RESORT TO MODE  
      }else if(.all_na(gate)){
        # ALL LABELS SAME COORD - OFFSET LATER
        lapply(seq_len(label_text_X), function(z){
          # X COORD - MODE
          if(.all_na(label_text_x[z])){
            label_text_x[z] <<- suppressMessages(
              .cyto_mode(fr_list[[1]],
                         channels = channels[1],
                         density_smooth = density_smooth)
            )
          }
          # Y COORD - MODE
          if(.all_na(label_text_x[z])){
            label_text_y[z] <<- suppressMessages(
              .cyto_mode(fr_list[[1]],
                         channels = channels[2],
                         density_smooth = density_smooth)
            )
          }
        })
        # LABEL_COORDS
        label_coords <- cbind(label_text_x, label_text_y)
        colnames(label_coords) <- c("x","y")
      }
    }
  }
  
  # PREPARE LABEL_TEXT WITH STATISTICS -----------------------------------------
  
  # 1D PLOT
  if(length(channels) == 1){
    # NO STACKING - BASE LAYER ONLY
    if(density_stack == 0){
      # GATE - ALL STATISTICS SUPPORTED
      if(!.all_na(gate)){
        # STATISTICS PER POPULATION
        lapply(seq_len(pops[[1]]), function(z){
          # CALCULATE STATISTIC
          if(!.all_na(label_stat[z])){
            # FREQUENCY
            if(label_stat[z] == "freq"){
              label_stat[z] <<- .cyto_count(pops[[1]][[z]]) / 
                .cyto_count(fr_list[[1]]) * 100
              label_stat[z] <<- paste(.round(label_stat[z], 2),"%")
            # CV
            }else if(label_stat[z] == "CV"){
              label_stat[z] <<- .cyto_CV(pops[[1]][[z]],
                                        channels = channels,
                                        trans = axes_trans)
              label_stat[z] <<- paste(.round(label_stat[z], 2), "%")
            # OTHER STATISTICS
            }else{
              label_stat[z] <<- cyto_stats_compute(pops[[1]][[z]],
                                                   channels = channels,
                                                   trans = axes_trans,
                                                   stat = label_stat[z],
                                                   format = "long",
                                                   density_smooth = density_smooth)
              
            }
          }
        })
      # NO GATE
      }else if(.all_na(gate)){
        # STATISTICS PER LABEL
        lapply(seq_len(LBLS), function(z){
          # CALCULATE STATISTIC
          if(!.all_na(label_stat[z])){
            # FREQUENCY
            if(label_stat[z] == "freq"){
              label_stat[z] <<- .cyto_count(fr_list[[1]]) / 
                .cyto_count(fr_list[[1]]) * 100
              label_stat[z] <<- paste(.round(label_stat[z], 2),"%")
              # CV
            }else if(label_stat[z] == "CV"){
              label_stat[z] <<- .cyto_CV(fr_list[[1]],
                                         channels = channels,
                                         trans = axes_trans)
              label_stat[z] <<- paste(.round(label_stat[z], 2), "%")
              # OTHER STATISTICS
            }else{
              label_stat[z] <<- cyto_stats_compute(fr_list[[1]],
                                                   channels = channels,
                                                   trans = axes_trans,
                                                   stat = label_stat[z],
                                                   format = "long",
                                                   density_smooth = density_smooth)
              
            }
          }
        })
        
      }
    # STACKING
    }else if(density_stack != 0){
      # GATE
      if(!.all_na(gate)){
        # STATISTICS PER POPULATION LIST
        lapply(seq_len(length(pops)), function(z){
          # STATISTICS PER POPULATION
          lapply(seq_len(length(pops[[z]])), function(y){
            # INDEX IN VECTOR
            w <- z * y
            # CALCULATE STATISTICS
            if(!.all_na(label_stat[w])){
              # FREQUENCY
              if(label_stat[w] == "freq"){
                label_stat[w] <<- .cyto_count(pops[[z]][[y]]) / 
                  .cyto_count(fr_list[[z]]) * 100
                label_stat[w] <<- paste(.round(label_stat[w], 2),"%")
              # CV
              }else if(label_stat[w] == "CV"){
                label_stat[w] <<- .cyto_CV(pops[[z]][[y]],
                                           channels = channels,
                                           trans = axes_trans)
                label_stat[w] <<- paste(.round(label_stat[w], 2), "%")
              # OTHER STATISTICS
              }else{
                label_stat[w] <<- cyto_stats_compute(pops[[z]][[y]],
                                                     channels = channels,
                                                     trans = axes_trans,
                                                     stat = label_stat[w],
                                                     format = "long",
                                                     density_smooth = density_smooth)
              }
            }
          })
        })
        # NO GATE
      }else if(.all_na(gate)){
        
      }
    }
  # 2D PLOT - FLUORESCENT INTENSITY STATISTICS NOT SUPPORTED (2 DIMENSIONS)
  }else if(length(channels) == 2){
    # GATE - FREQ OR COUNT ONLY
    if(!.all_na(gate)){
      # STATISTICS PER POPULATION
      lapply(seq_len(pops[[1]]), function(z){
        # CALCULATE STATISTIC
        if(!.all_na(label_stat[z])){
          # UNSUPPORTED STATISTIC
          if(!label_stat[z] %in% c("freq","count")){
            stop("Only 'freq' and 'count' are supported for gated 2D plots.")
          }
          # FREQUENCY
          if(label_stat[z] == "freq"){
            label_stat[z] <<- .cyto_count(pops[[1]][[z]]) / 
              .cyto_count(fr_list[[1]]) * 100
            label_stat[z] <<- paste(.round(label_stat[z], 2),"%")
            # COUNT
          }else if(label_stat[z] == "count"){
            label_stat[z] <<- .cyto_count(pops[[1]][[z]])
          }
        }
      })
      # NO GATE - COUNT ONLY
    }else if(.all_na(gate)){
      # STATISTICS PER LABEL
      lapply(seq_len(LBLS), function(z){
        # CALCULATE STATISTIC
        if(!.all_na(label_stat[z])){
          # UNSUPPORTED STATISTIC
          if(!label_stat[z] == "count"){
            stop("Only 'count' is supported for un-gated 2D plots.")
          }
          # FREQUENCY
          if(label_stat[z] == "count"){
            label_stat[z] <<- .cyto_count(fr_list[[1]])
          }
        }
      })
    }
  }
  
  # PREPARE LABEL_TEXT ---------------------------------------------------------
  
  # COMBINE LABEL_TEXT & LABEL_STAT
  lapply(seq_len(label_text), function(z){
    if(!.all_na(label_text[z]) & !.all_na(label_stat[z])){
      label_text[z] <<- paste(label_text_x, label_stat[z], sep = "\n")
    }else if(.all_na(label_text[z]) & !.all_na(label_stat[z])){
      label_text[z] <<- label_stat[z]
    }
  })
  
  # LABEL OFFSET ---------------------------------------------------------------
  
  # OFFSET OVERLAPPING LABELS
  if(label_position == "auto"){
    
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
    cyto_plot_label(fr_list[[1]], # IRRELEVANT - NO CALCULATIONS
                    channels = channels,
                    trans = trans,
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
