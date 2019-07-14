# Internal Function Definition for cyto_plot -----------------------------------

# Handles construction of 1D or 2D plots based on length of channels.

# .cyto_plot acceepts a pre-processed list of flowFrames.

# ARGUMENTS: .cyto_plot inherits all arguments from cyto_plot. All data
# manipualtion and theme inheritance occurs at the cyto_plot level. The only
# data manipulation reserved for .cyto_plot is the sampling through the display
# argument (it is easier this way). All checks should have been performed at the
# upper cyto_plot level.

# ARGUMENT CONVENTIONS: Arguments that can be replaced internally are set to
# missing by default and replaced if not assigned NA. Arguments that can be
# either supplied or not are set to NA by default. Stick to using "channels"
# for arguments accepting channel name(s) this make it easier to pass arguments
# through do.call.

# MISSING ARGUMENTS:
# - title
# - xlab
# - ylab
# - label
# - label_stat

# ADDING NEW FEATURES TO CYTO_PLOT:
# 1. Add the argument to cyto_plot with an appropriate default.
# 2. Add correct implementation in .cyto_plot_args_split to repeat arguments as
# required for downstream mapply call.
# 3. Add argument to downstream mapply calls.
# 3. Modify the code in the called .cyto_plot function (below) to add
# the new feature.

#' @importFrom flowCore exprs parameters identifier
#' @importFrom flowWorkspace pData
#' @importFrom graphics plot axis title abline polygon legend par box
#' @importFrom grDevices adjustcolor
#' @importFrom stats density
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @noRd
.cyto_plot <- function(x,
                       channels,
                       gate = NA,
                       axes_trans = NA,
                       limits = "machine",
                       display = 1,
                       popup = FALSE,
                       xlim = NA,
                       ylim = NA,
                       xlab,
                       ylab,
                       title,
                       negate = FALSE,
                       title_text_font = 2,
                       title_text_size = 1.1,
                       title_text_col = "black",
                       density_modal = TRUE,
                       density_smooth = 0.6,
                       density_stack = 0.5,
                       density_cols = NA,
                       density_fill = NA,
                       density_fill_alpha = 1,
                       density_line_type = 1,
                       density_line_width = 1,
                       density_line_col = "black",
                       point_shape = ".",
                       point_size = 2,
                       point_col_scale = NA,
                       point_cols = NA,
                       point_col = NA,
                       point_col_alpha = 1,
                       contour_lines = 0,
                       contour_line_type = 1,
                       contour_line_width = 1,
                       contour_line_col = "black",
                       axes_text = list(TRUE, TRUE),
                       axes_text_font = 1,
                       axes_text_size = 1,
                       axes_text_col = "black",
                       axes_label_text_font = 1,
                       axes_label_text_size = 1.1,
                       axes_label_text_col = "black",
                       legend = FALSE,
                       legend_text = NA,
                       legend_text_font = 1,
                       legend_text_size = 1,
                       legend_text_col = "black",
                       legend_line_type = NA,
                       legend_line_width = NA,
                       legend_line_col = NA,
                       legend_box_fill = NA,
                       legend_point_col = NA,
                       gate_line_type = 1,
                       gate_line_width = 2.5,
                       gate_line_col = "red",
                       gate_fill = "white",
                       gate_fill_alpha = 0,
                       label,
                       label_text = NA,
                       label_stat,
                       label_position = "auto",
                       label_text_x = NA,
                       label_text_y = NA,
                       label_text_font = 2,
                       label_text_size = 1,
                       label_text_col = "black",                       
                       label_fill = "white",
                       label_fill_alpha = 0.6,
                       border_line_type = 1,
                       border_line_width = 1,
                       border_line_col = "black",
                       border_fill = "white",
                       border_fill_alpha = 1, ...) {

  # Get current graphics parameters and reset on exit
  pars <- par("mar")
  on.exit(par(pars))

  # ARGUMENTS ------------------------------------------------------------------

  # Pull down arguments to named list - convert missing to ""
  args <- .args_list()

  # Update arguments
  .args_update(args)

  # SAMPLING -------------------------------------------------------------------

  # Sample to display percentage/# events - set seed for consistency
  if (display != 1) {
    x <- cyto_sample(x, display = display, seed = 56)
  }

  # CONSERVED ARGUMENTS --------------------------------------------------------

  # TITLE
  title <- .cyto_plot_title(x[[1]],
    channels = channels,
    overlay = x[2:length(x)],
    title = title
  )
  # AXES LABELS
  labs <- .cyto_plot_axes_label(x[[1]],
    channels = channels,
    xlab = xlab,
    ylab = ylab,
    density_modal = density_modal
  )

  # XLAB
  xlab <- labs[[1]]

  # YLAB
  ylab <- labs[[2]]

  # LEGEND TEXT - required for setting plot margins
  if (.all_na(legend_text)) {
    legend_text <- names(x)
  }

  # 1D DENSITY DISTRIBUTIONS ---------------------------------------------------
  if (length(channels) == 1) {

    # Get kernel density for each list element
    fr_dens <- lapply(x, function(z) {
      suppressWarnings(.cyto_density(z,
        channel = channels,
        smooth = density_smooth,
        modal = density_modal
      ))
    })

    # Number of overlays
    ovn <- length(fr_dens) - 1

    # fr_dens does contain some valid density objects
    if (!.all_na(fr_dens)) {
      # Calculate the mean maximum y value for kernel densities
      if (density_modal) {
        y_max <- 100
      } else {
        y_max <- mean(unlist(lapply(fr_dens, function(d) {
          if (!.all_na(d)) {
            max(d$y)
          } else {
            NA
          }
        })), na.rm = TRUE)
      }

      # Stacked distributions require shifting of y values
      shft <- seq(
        0,
        ovn * density_stack * y_max,
        density_stack * y_max
      )

      # Adjust y values if stacking is been applied
      if (density_stack > 0 & length(x) > 1) {
        # Shift distributions for stacking
        lapply(seq_len(length(fr_dens)), function(z) {
          if (!.all_na(fr_dens[[z]])) {
            fr_dens[[z]]$y <<- fr_dens[[z]]$y + shft[z]
          }
        })
      }

      # fr_dens is composed of NA - no events in any fr_list
    } else if (.all_na(fr_dens)) {

      # Turn off y axis text
      axes_text[[2]] <- FALSE

      # Set y_max to 100
      y_max <- 100

      # Set y axis limits to 0-100
      ylim <- c(0, y_max + ovn * density_stack * y_max)

      # Stacked distributions require shifting of y values
      shft <- seq(
        0,
        ovn * density_stack * y_max,
        density_stack * y_max
      )
    }

    # EMPTY PLOT - handles margins and axes limits internally
    cyto_plot_empty(x,
      channels = channels,
      axes_trans = axes_trans,
      xlim = xlim,
      ylim = ylim,
      limits = limits,
      title = title,
      xlab = xlab,
      ylab = ylab,
      density_modal = density_modal,
      density_smooth = density_smooth,
      density_stack = density_stack,
      axes_text = axes_text,
      axes_text_font = axes_text_font,
      axes_text_size = axes_text_size,
      axes_text_col = axes_text_col,
      axes_label_text_font = axes_label_text_font,
      axes_label_text_size = axes_label_text_size,
      axes_label_text_col = axes_label_text_col,
      title_text_font = title_text_font,
      title_text_size = title_text_size,
      title_text_col = title_text_col,
      border_line_type = border_line_type,
      border_line_width = border_line_width,
      border_line_col = border_line_col,
      border_fill = border_fill,
      border_fill_alpha = border_fill_alpha,
      legend = legend,
      legend_text = legend_text,
      legend_text_size = legend_text_size
    )

    # DENSITY FILL - inherits theme internally
    density_fill <- .cyto_plot_density_fill(fr_dens,
      density_fill = density_fill,
      density_cols = density_cols,
      density_fill_alpha = density_fill_alpha
    )

    # DENSITY - no alpha adjustment here - happens above
    cyto_plot_density(fr_dens,
      density_modal = density_modal,
      density_stack = density_stack,
      density_cols = density_cols,
      density_fill = density_fill,
      density_fill_alpha = 1,
      density_line_type = density_line_type,
      density_line_width = density_line_width,
      density_line_col = density_line_col
    )
    
    # LEGEND
    if (legend != FALSE) {
      .cyto_plot_legend(
        channels = channels,
        legend = legend,
        legend_text = legend_text,
        legend_text_font = legend_text_font,
        legend_text_size = legend_text_size,
        legend_text_col = legend_text_col,
        legend_line_type = legend_line_type,
        legend_line_width = legend_line_width,
        legend_line_col = legend_line_col,
        legend_box_fill = legend_box_fill,
        legend_point_col = legend_point_col,
        density_fill = density_fill,
        density_fill_alpha = density_fill_alpha,
        density_line_type = density_line_type,
        density_line_width = density_line_width,
        density_line_col = density_line_col,
        point_shape = point_shape,
        point_size = point_size,
        point_col = point_col,
        point_col_alpha = point_col_alpha
      )
    }

    # Turn off gates and labels if no density in fr_dens
    if (.all_na(fr_dens)) {
      label <- FALSE
      gate <- NA
    }

    # GATES & LABEL - without overlay
    if (length(x) == 1) {

      # LABEL?
      if (!.all_na(gate) &
        .empty(label)) {
        label <- TRUE # turn on labels if gate
      } else if (.all_na(gate) &
        .empty(label)) {
        label <- FALSE # turn off labels if no gate
      }
      
      # STAT
      if (!.all_na(gate) &
        .empty(label_stat)) {
        label_stat <- "freq"
      } else if (.all_na(gate) &
        .empty(label_stat)) {
        label_stat <- NA
      }
      
      # GATES
      if (!.all_na(gate)) {
        gate <- cyto_plot_gate2(x[[1]],
                               channels = channels,
                               gate = gate,
                               label = label,
                               label_text = label_text,
                               label_stat = label_stat,
                               label_text_x = label_text_x,
                               label_text_y = label_text_y,
                               label_position = label_position,
                               trans = axes_trans,
                               negate = negate,
                               display = display,
                               gate_line_type = gate_line_type,
                               gate_line_width = gate_line_width,
                               gate_line_col = gate_line_col,
                               gate_fill = gate_fill,
                               gate_fill_alpha = gate_fill_alpha,
                               label_text_font = label_text_font,
                               label_text_size = label_text_size,
                               label_text_col = label_text_col,
                               label_fill = label_fill,
                               label_fill_alpha = label_fill_alpha,
                               density_smooth = density_smooth
        )
      }

      # GATE & LABEL - with overlay
    } else if(length(x) != 1) {

      # Calculate label y positions for stacked overlays
      if (density_stack != 0) {

        # Need to compute y label positions
        if (.all_na(label_text_y) & label_position == "auto") {
          # Default is half way beteen horizontal lines (shft)
          shft <- c(shft, shft[length(shft)] + shft[2])
          label_text_y <- unlist(lapply(seq(1, length(shft) - 1), function(z) {
            (shft[z] + shft[z + 1]) / 2
          }))

          # Repeat gate times if gate supplied - one label per layer if no gates
          if (!.all_na(gate)) {
            label_text_y <- rep(label_text_y, each = length(gate))
          }
        }

        # LABEL?
        if (!.all_na(gate) &
          .empty(label)) {
          label <- TRUE # turn on labels if gate
        } else if (.all_na(gate) &
          .empty(label)) {
          label <- FALSE # turn off labels if no gate
        }
        
        # STAT
        if (!.all_na(gate) &
          .empty(label_stat)) {
          label_stat <- "freq"
        } else if (.all_na(gate) &
          .empty(label_stat)) {
          label_stat <- NA
        }

        # Gate density overlays
        .cyto_plot_overlay_gate(
          x = x[[1]],
          channel = channels,
          trans = axes_trans,
          overlay = x[2:length(x)],
          gate = gate,
          density_stack = density_stack,
          density_modal = density_modal,
          label = label,
          label_text = label_text,
          label_stat = label_stat,
          label_text_size = label_text_size,
          label_text_font = label_text_font,
          label_text_col = label_text_col,
          label_text_x = label_text_x,
          label_text_y = label_text_y,
          label_fill_alpha = label_fill_alpha,
          gate_line_col = gate_line_col,
          gate_line_width = gate_line_width,
          gate_line_type = gate_line_type
        )
        
      } 
    }

    # 2D SCATTER PLOTS -----------------------------------------------------------
  } else if (length(channels) == 2) {

    # REMOVE NEGATIVE FSC/SSC EVENTS - DENSITY BUG (other linear channels?)
    lapply(seq_len(2), function(z){
      if(grepl("FSC", channels[z], ignore.case = TRUE) |
         grepl("SSC", channels[z], ignore.case = TRUE)){
        
        x <<- lapply(x, function(y){
          if(min(range(y, type = "data")[, channels[z]]) < 0){
            # Gate out negative events
            coords <- matrix(c(0, Inf), ncol = 1, nrow = 2)
            rownames(coords) <- c("min","max")
            colnames(coords) <- channels[z]
            nonDebris <- rectangleGate(.gate = coords)
            y <- Subset(y, nonDebris)
            # Fix plot limits
            if(z == 1){
              if(xlim[1] < 0){
                xlim[1] <<- 0
              }
            }else if(z == 2){
              if(ylim[1] < 0){
                ylim[1] <<- 0
              }
            }
          }
          return(y)
        })
        
        xlim <<- xlim
        ylim <<- ylim
        
      }
    })
    
    # EMPTY PLOT - handles margins and axes limits internally
    cyto_plot_empty(x,
      channels = channels,
      axes_trans = axes_trans,
      xlim = xlim,
      ylim = ylim,
      limits = limits,
      title = title,
      xlab = xlab,
      ylab = ylab,
      density_modal = density_modal,
      density_smooth = density_smooth,
      density_stack = density_stack,
      axes_text = axes_text,
      axes_text_font = axes_text_font,
      axes_text_size = axes_text_size,
      axes_text_col = axes_text_col,
      axes_label_text_font = axes_label_text_font,
      axes_label_text_size = axes_label_text_size,
      axes_label_text_col = axes_label_text_col,
      title_text_font = title_text_font,
      title_text_size = title_text_size,
      title_text_col = title_text_col,
      border_line_type = border_line_type,
      border_line_width = border_line_width,
      border_line_col = border_line_col,
      border_fill = border_fill,
      border_fill_alpha = border_fill_alpha,
      legend = legend,
      legend_text = legend_text,
      legend_text_size = legend_text_size
    )
    
    # POINT COL - list
    if (.all_na(point_col)) {
      point_col <- .cyto_plot_point_col(x,
        channels = channels,
        point_col_scale = point_col_scale,
        point_cols = point_cols,
        point_col = point_col,
        point_col_alpha = point_col_alpha
      )
    }
    
    # POINTS - base layer
    cyto_plot_point(x[1],
      channels = channels,
      point_shape = point_shape[1],
      point_size = point_size[1],
      point_col_scale = point_col_scale,
      point_cols = point_cols,
      point_col = point_col[1],
      point_col_alpha = point_col_alpha[1]
    )

    # CONTOUR LINES -  only supported for base layer
    if (contour_lines != 0) {
      cyto_plot_contour(
        x = x[[1]],
        channels = channels,
        contour_lines = contour_lines[1],
        contour_line_type = contour_line_type[1],
        contour_line_width = contour_line_width[1],
        contour_line_col = contour_line_col[1]
      )
    }

    # POINTS - overlay layers
    if (length(x) > 1) {
      cyto_plot_point(x[-1],
        channels = channels,
        point_shape = point_shape[-1],
        point_size = point_size[-1],
        point_col_scale = point_col_scale,
        point_cols = point_cols,
        point_col = point_col[-1],
        point_col_alpha = point_col_alpha[-1]
      )
    }

    # POINT DENSITY COLOUR SCALE
    point_cols <- .cyto_plot_point_cols(point_cols)

    # POINT_COL LEGEND - vector (replace density colours with first point_cols)
    point_col <- unlist(lapply(point_col, function(z) {

      # colours defined for each point
      if (length(z) > 1) {
        return(point_cols[1])
      } else {
        return(z)
      }
    }))
    
    # LEGEND
    if (legend != FALSE) {
      .cyto_plot_legend(
        channels = channels,
        legend = legend,
        legend_text = legend_text,
        legend_text_font = legend_text_font,
        legend_text_size = legend_text_size,
        legend_text_col = legend_text_col,
        legend_line_type = legend_line_type,
        legend_line_width = legend_line_width,
        legend_line_col = legend_line_col,
        legend_box_fill = legend_box_fill,
        legend_point_col = legend_point_col,
        density_fill = density_fill,
        density_fill_alpha = density_fill_alpha,
        density_line_type = density_line_type,
        density_line_width = density_line_width,
        density_line_col = density_line_col,
        point_shape = point_shape,
        point_size = point_size,
        point_col = point_col,
        point_col_alpha = point_col_alpha
      )
    }

    # If no events turn off gates and labels
    if (all(unlist(lapply(x, "nrow")) == 0)) {
      gate <- NA
      label <- FALSE
    }

    # LABEL
    if (.empty(label)) {
      label <- TRUE
    }
    
    # LABEL STAT
    if (!.all_na(gate) &
      .empty(label_stat)) {
      label_stat <- "freq"
    } else if (.all_na(gate) &
      .empty(label_stat)) {
      label_stat <- NA
    }
    
    # GATES
    if (!.all_na(gate)) {
      gate <- cyto_plot_gate2(x = x[[1]],
                             channels = channels,
                             gate = gate,
                             label = label,
                             label_text = label_text,
                             label_stat = label_stat,
                             label_text_x = label_text_x,
                             label_text_y = label_text_y,
                             label_position = label_position,
                             trans = axes_trans,
                             negate = negate,
                             display = display,
                             gate_line_type = gate_line_type,
                             gate_line_width = gate_line_width,
                             gate_line_col = gate_line_col,
                             gate_fill = gate_fill,
                             gate_fill_alpha = gate_fill_alpha,
                             label_text_font = label_text_font,
                             label_text_size = label_text_size,
                             label_text_col = label_text_col,
                             label_fill = label_fill,
                             label_fill_alpha = label_fill_alpha,
                             density_smooth = density_smooth
      )

    # LABELS WITHOUT GATES
    } else if (.all_na(gate) &
      !.all_na(label_text) &
      label == TRUE) {

      # label - limited to # layers - arg_split
      text_xy <- mapply(
        function(label_text,
                         label_stat,
                         label_text_x,
                         label_text_y,
                         label_text_font,
                         label_text_size,
                         label_text_col,
                         label_fill,
                         label_fill_alpha) {
          if (!.all_na(label_stat) & label_stat == "percent") {
            label_stat <- NA
          }
          suppressMessages(cyto_plot_label2(
            x = x[[1]],
            channels = channels,
            gate = gate,
            trans = axes_trans,
            label_text = label_text,
            label_stat = label_stat,
            label_text_x = label_text_x,
            label_text_y = label_text_y,
            trans = axes_trans,
            negate = negate,
            display = display,
            label_text_font = label_text_font,
            label_text_size = label_text_size,
            label_text_col = label_text_col,
            label_fill = label_fill,
            label_fill_alpha = label_fill_alpha,
            density_smooth = density_smooth
          ))
        }, label_text,
        label_stat,
        label_text_x,
        label_text_y,
        label_text_font,
        label_text_size,
        label_text_col,        
        label_fill,
        label_fill_alpha,
        SIMPLIFY = FALSE
      )

    }
  }
}
