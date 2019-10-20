## CYTO_PLOT CONSTRUCTORS ------------------------------------------------------

# Collection of internal wrappers used within cyto_plot to construct aspects of
# the plot. This includes wrappers for cyto_plot_empty, cyto_plot_density,
# cyto_plot_point, cyto_plot_gate and cyto_plot_label. These wrappers bypass
# some of the calculations in the exported equivalents as these calculations are
# already performed internally within cyto_plot. All of these functions inherit
# a named list of arguments from cyto_plot. Arguments must be extracted directly
# from the list to prevent R CMD CHECK global binding NOTEs.

## .CYTO_PLOT_EMPTY ------------------------------------------------------------

#' Construct empty cyto_plot
#'
#' @param x named list of cyto_plot arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_empty <- function(x) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REMOVE X ARGUMENT
  x <- x[-match("x", names(x))]

  # RENAME FR_LIST ARGUMENT TO X (DO.CALL() ON CYTO_PLOT_EMPTY)
  names(x)[match("fr_list", names(x))] <- "x"

  # UPDATE ARGUMENTS - missing converted to ""
  .args_update(x)

  # PLOT CONSTRUCTION ----------------------------------------------------------

  # Call flowFrame method - no density objects supplied
  if (.all_na(fr_dens_list)) {

    # PULL DOWN UPDATED ARGUMENTS
    args <- .args_list()

    # CYTO_PLOT_EMPTY ARGUMENTS
    ARGS <- formalArgs("cyto_plot_empty.list")

    # CYTO_PLOT_EMPTY
    do.call("cyto_plot_empty.list", args[names(args) %in% ARGS])

    # DENSITY SUPPLIED
  } else if (length(channels) == 1 & !.all_na(fr_dens_list)) {

    # GRAPHICAL PARAMETERS -----------------------------------------------------

    # Prevent scientific notation on axes - reset on exit
    scipen <- getOption("scipen")
    options(scipen = 100000000)
    on.exit(options(scipen = scipen))

    # Extract current graphics parameters
    pars <- par("mar")

    # Reset graphics parameters on exit
    on.exit(par(pars))

    # GENERAL ------------------------------------------------------------------

    # SAMPLES
    SMP <- length(x)

    # CHANNELS -----------------------------------------------------------------

    # EXTRACT CHANNELS
    channels <- cyto_channels_extract(
      x[[1]],
      channels,
      TRUE
    )

    # AXES LIMITS --------------------------------------------------------------

    # GATE COORDS
    if(!.all_na(gate)){
      gate_coords <- .cyto_gate_coords(gate, channels)
    }
    # XLIM
    if (.all_na(xlim)) {
      xlim <- .cyto_range(x,
        channels = channels[1],
        limits = limits,
        plot = TRUE
      )[, 1]
    # GATE coords
    if(!.all_na(gate)){
      # MIN & MAX GATE COORDS
      gate_xcoords <- gate_coords[, channels[1]]
      gate_xcoords <- c(min(gate_xcoords), max(gate_xcoords))
      # GATE COORDS BELOW XMIN
      if(is.finite(gate_xccords[1]) & gate_xcoords[1] < xlim[1]){
        xlim[1] <- gate_xcoords[1]
      }
      # GATE COORDS ABOVE XMAX
      if(is.finite(gate_xcoords[2]) & gate_xcoords[2] > xlim[2]){
        xlim[2] <- gate_xcoords[2]
      }
    }      
      # XLIM SUPPLIED MANUALLY
    } else {
      if (getOption("cyto_plot_method") == "flowFrame") {
        if (!.all_na(axes_trans)) {
          if (channels[1] %in% names(axes_trans)) {
            xlim <- axes_trans[[channels[1]]]$transform(xlim)
          }
        }
      }
    }

    # YLIM
    if (.all_na(ylim)) {
      # EXTRACT YLIM - NAMES
      ymin <- as.numeric(unlist(strsplit(names(fr_dens_list)[1], "-"))[1])
      ymax <- as.numeric(unlist(strsplit(names(fr_dens_list)[SMP], "-"))[2])
      ylim <- c(ymin, ymax)
    }

    # AXES TEXT ----------------------------------------------------------------

    # Convert axes_text to list - allows inheritance from cyto_plot
    if (!inherits(axes_text, "list")) {
      axes_text <- list(axes_text[1], axes_text[2])
    }

    # X axis breaks and labels -  can be inherited from cyto_plot
    if (!inherits(axes_text[[1]], "list")) {
      if (.all_na(axes_text[[1]])) {
        # NA == TRUE returns NA not T/F
      } else if (axes_text[[1]] == TRUE) {
        lims <- list(xlim)
        names(lims) <- channels[1]
        axes_text[[1]] <- .cyto_plot_axes_text(fr_list,
                                               channels = channels[1],
                                               axes_trans = axes_trans,
                                               axes_range = lims,
                                               limits = limits
        )[[1]]
      }
    }
    
    # Y axis breaks and labels - can be inherited from cyto_plot
    if (!inherits(axes_text[[2]], "list")) {
      if (.all_na(axes_text[[2]])) {
        # NA == TRUE returns NA not T/F
      } else if (axes_text[[2]] == TRUE) {
        if (length(channels) == 2) {
          lims <- list(ylim)
          names(lims) <- channels[2]
          axes_text[[2]] <- .cyto_plot_axes_text(fr_list,
                                                 channels = channels[2],
                                                 axes_trans = axes_trans,
                                                 axes_range = lims,
                                                 limits = limits
          )[[1]]
        } else {
          axes_text[[2]] <- NA
        }
      }
    }

    # Turn off y axis labels for stacked overlays
    if (!.all_na(overlay) &
      density_stack != 0 &
      length(channels) == 1) {
      axes_text <- list(axes_text[[1]], FALSE)
    }

    # AXES LABELS --------------------------------------------------------------

    # AXES LABELS - missing replaced - NA removed
    axes_labels <- .cyto_plot_axes_label(x[[1]],
      channels = channels,
      xlab = xlab,
      ylab = ylab,
      density_modal = density_modal
    )
    xlab <- axes_labels[[1]]
    ylab <- axes_labels[[2]]

    # TITLE --------------------------------------------------------------------

    # TITLE - missing replaced - NA removed
    title <- .cyto_plot_title(x[[1]],
      channels = channels,
      overlay = overlay,
      title = title
    )

    # MARGINS ------------------------------------------------------------------

    # Set plot margins - set par("mar")
    .cyto_plot_margins(x,
      legend = legend,
      legend_text = legend_text,
      legend_text_size = legend_text_size,
      title = title,
      axes_text = axes_text
    )

    # PLOT CONSTRUCTION --------------------------------------------------------

    # Plot
    graphics::plot(1,
      type = "n",
      axes = FALSE,
      xlim = xlim,
      ylim = ylim,
      xlab = "",
      ylab = "",
      bty = "n"
    )

    # X AXIS - TRANSFORMED
    if (inherits(axes_text[[1]], "list")) {
      # MINOR TICKS
      mnr_ind <- which(as.character(axes_text[[1]]$label) == "")
      axis(1,
        at = axes_text[[1]]$at[mnr_ind],
        labels = axes_text[[1]]$label[mnr_ind],
        tck = -0.015
      )

      # MAJOR TICKS - MUST BE >2% XRANGE FROM ZERO
      mjr_ind <- which(as.character(axes_text[[1]]$label) != "")
      mjr <- list(
        "at" = axes_text[[1]]$at[mjr_ind],
        "label" = axes_text[[1]]$label[mjr_ind]
      )
      # Zero included on plot
      if (any(as.character(mjr$label) == "0")) {
        zero <- which(as.character(mjr$label) == "0")
        zero_break <- mjr$at[zero]
        zero_buffer <- c(
          zero_break - 0.02 * (xlim[2] - xlim[1]),
          zero_break + 0.02 * (xlim[2] - xlim[1])
        )
        mjr_ind <- c(
          zero,
          which(mjr$at < zero_buffer[1] |
            mjr$at > zero_buffer[2])
        )
      } else {
        mjr_ind <- seq_len(length(mjr$label))
      }

      axis(1,
        at = mjr$at[mjr_ind],
        labels = mjr$label[mjr_ind],
        font.axis = axes_text_font,
        col.axis = axes_text_col,
        cex.axis = axes_text_size,
        tck = -0.03
      )
      # X AXIS - UNTRANSFORMED
    } else if (.all_na(axes_text[[1]])) {
      axis(1,
        font.axis = axes_text_font,
        col.axis = axes_text_col,
        cex.axis = axes_text_size,
        tck = -0.03
      )
    }

    # Y AXIS - TRANSFORMED
    if (inherits(axes_text[[2]], "list")) {
      # MINOR TICKS
      mnr_ind <- which(as.character(axes_text[[2]]$label) == "")
      axis(2,
        at = axes_text[[2]]$at[mnr_ind],
        labels = axes_text[[2]]$label[mnr_ind],
        tck = -0.015
      )
      # MAJOR TICKS - MUST BE >2% yrange FROM ZERO
      mjr_ind <- which(as.character(axes_text[[2]]$label) != "")
      mjr <- list(
        "at" = axes_text[[2]]$at[mjr_ind],
        "label" = axes_text[[2]]$label[mjr_ind]
      )
      # Zero included on plot
      if (any(as.character(mjr$label) == "0")) {
        zero <- which(as.character(mjr$label) == "0")
        zero_break <- mjr$at[zero]
        zero_buffer <- c(
          zero_break - 0.02 * (ylim[2] - ylim[1]),
          zero_break + 0.02 * (ylim[2] - ylim[1])
        )
        mjr_ind <- c(
          zero,
          which(mjr$at < zero_buffer[1] |
            mjr$at > zero_buffer[2])
        )
      } else {
        mjr_ind <- seq_len(length(mjr$label))
      }

      axis(2,
        at = mjr$at[mjr_ind],
        labels = mjr$label[mjr_ind],
        font.axis = axes_text_font,
        col.axis = axes_text_col,
        cex.axis = axes_text_size,
        tck = -0.03
      )
      # Y AXIS - LINEAR
    } else if (.all_na(axes_text[[2]])) {
      axis(2,
        font.axis = axes_text_font,
        col.axis = axes_text_col,
        cex.axis = axes_text_size,
        tck = -0.03
      )
    }

    # BORDER
    box(
      which = "plot",
      lty = border_line_type,
      lwd = border_line_width,
      col = border_line_col
    )


    # BORDER_FILL
    if (border_fill != "white") {
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
        col = adjustcolor(border_fill, border_fill_alpha),
        border = NA
      )
    }

    # TITLE
    if (!.all_na(title)) {
      title(
        main = title,
        cex.main = title_text_size,
        col.main = title_text_col,
        font.main = title_text_font
      )
    }

    # XLAB - position labels closer if axes text is missing
    if (!.all_na(xlab)) {
      if (inherits(axes_text[[1]], "list")) {
        title(
          xlab = xlab,
          font.lab = axes_label_text_font,
          col.lab = axes_label_text_col,
          cex.lab = axes_label_text_size
        )
      } else if (.all_na(axes_text[[1]])) {
        title(
          xlab = xlab,
          font.lab = axes_label_text_font,
          col.lab = axes_label_text_col,
          cex.lab = axes_label_text_size
        )
      } else if (axes_text[[1]] == FALSE) {
        title(
          xlab = xlab,
          font.lab = axes_label_text_font,
          col.lab = axes_label_text_col,
          cex.lab = axes_label_text_size,
          mgp = c(2, 0, 0)
        )
      }
    }

    # YLAB - position labels closer if axes text is missing
    if (!.all_na(ylab)) {
      if (inherits(axes_text[[2]], "list")) {
        title(
          ylab = ylab,
          font.lab = axes_label_text_font,
          col.lab = axes_label_text_col,
          cex.lab = axes_label_text_size
        )
      } else if (.all_na(axes_text[[2]])) {
        title(
          ylab = ylab,
          font.lab = axes_label_text_font,
          col.lab = axes_label_text_col,
          cex.lab = axes_label_text_size
        )
      } else if (axes_text[[2]] == FALSE) {
        title(
          ylab = ylab,
          font.lab = axes_label_text_font,
          col.lab = axes_label_text_col,
          cex.lab = axes_label_text_size,
          mgp = c(2, 0, 0)
        )
      }
    }

    # LEGEND
    if (legend != FALSE) {
      .cyto_plot_legend(x,
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
        density_cols = density_cols,
        density_fill = density_fill,
        density_fill_alpha = density_fill_alpha,
        density_line_type = density_line_type,
        density_line_width = density_line_width,
        density_line_col = density_line_col,
        point_shape = point_shape,
        point_size = point_size,
        point_cols = point_cols,
        point_col = point_col,
        point_col_alpha = point_col_alpha
      )
    }
  }
}

## .CYTO_PLOT_DENSITY ----------------------------------------------------------

#' Add density distributions to an existing cyto_plot
#'
#' @param x named list of cyto_plot_arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_density <- function(x) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REMOVE X ARGUMENT
  x <- x[-match("x", names(x))]

  # RENAME FR_LIST ARGUMENT TO X (DO.CALL() ON CYTO_PLOT_EMPTY)
  names(x)[match("fr_dens_list", names(x))] <- "x"

  # CYTO_PLOT_DENSITY ARGUMENTS
  ARGS <- formalArgs("cyto_plot_density.list")

  # RESTRICT SUPPLIED ARGUMENTS
  args <- x[names(x) %in% ARGS]

  # CALL CYTO_PLOT_DENSITY -----------------------------------------------------
  do.call("cyto_plot_density.list", args)
}

## .CYTO_PLOT_POINT ------------------------------------------------------------

#' Add points to an existing cyto_plot
#'
#' @param x named list of cyto_plot arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_point <- function(x) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REMOVE X ARGUMENT
  x <- x[-match("x", names(x))]

  # RENAME FR_LIST ARGUMENT TO X (DO.CALL() ON CYTO_PLOT_EMPTY)
  names(x)[match("fr_list", names(x))] <- "x"

  # REST DISPLAY ARGUMENT
  x["display"] <- 1

  # CYTO_PLOT_DENSITY ARGUMENTS
  ARGS <- formalArgs("cyto_plot_point.list")

  # RESTRICT SUPPLIED ARGUMENTS
  args <- x[names(x) %in% ARGS]

  # CALL CYTO_PLOT_DENSITY -----------------------------------------------------
  do.call("cyto_plot_point.list", args)
}

## .CYTO_PLOT_GATE -------------------------------------------------------------

#' Add gates and labels to an existing cyto_plot
#'
#' @param x named list of cyto_plot arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_gate <- function(x) {

  # INHERIT ARGUMENTS ----------------------------------------------------------

  # UPDATE ARGUMENTS
  .args_update(x)

  # GENERAL --------------------------------------------------------------------

  # GATES
  NG <- length(gate)

  # POPULATIONS PER GATE
  P <- .cyto_gate_count(gate, negate = FALSE, total = FALSE)

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # PULL DOWN ARGUMENTS
  args <- .args_list()

  # CYTO_PLOT_GATE & CYTO_PLOT_LABELLER ARGUMENTS
  gate_args <- formalArgs("cyto_plot_gate.list")
  label_args <- formalArgs("cyto_plot_labeller")

  # RESTRICT ARGUMENTS
  args <- args[names(args) %in% c(gate_args, label_args)]

  # SPLIT GATE_FILL ARGUMENTS BY POPULATIONS
  lapply(gate_args[which(grepl("gate_fill", gate_args))], function(z) {
    args[[z]] <<- split(args[[z]], rep(seq_len(NG), times = P))
  })

  # SPLIT LABEL ARGUMENTS BY POPULATIONS PER LAYER
  lapply(label_args, function(z) {
    args[[z]] <<- rep(args[[z]], length.out = TNP)
  })
  lapply(label_args, function(z) {
    args[[z]] <<- split(args[[z]], rep(seq_len(SMP), each = NP))
  })

  # RE-ARRANGE LABEL COORDS PER GATE
  lapply(label_args, function(z) {
    args[[z]] <<- lapply(seq_len(NP), function(y) {
      LAPPLY(args[[z]], `[[`, y)
    })
  })

  # GATE ARGUMENTS
  gate_args <- args[names(args) %in% gate_args]

  # LABEL ARGUMENTS
  label_args <- args[names(args) %in% label_args]

  # GATE & ASSOCIATED LABELS ---------------------------------------------------

  # GATES & LABELS
  label_text_xy <- lapply(seq_len(NP), function(z) {
    # GATED POPULATION(S) - GATE & LABEL(S)
    if (z <= NG) {
      # PLOT GATE
      do.call(
        "cyto_plot_gate",
        c(
          list("channels" = NULL),
          lapply(
            gate_args[!grepl("channels", names(gate_args))],
            `[[`, z
          )
        )
      )
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (!.all_na(label_text[[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(label_text[z])),
          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
      # NEGATED POPULATION(S) - LABEL(s) ONLY
    } else {
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (!.all_na(label_text[[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(label_text[[z]])),
          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
    }
    return(text_xy)
  })
  label_text_xy <- do.call("rbind", label_text_xy)

  # RE-ARRANGE LABEL ARGUMENTS -------------------------------------------------

  # UPDATE LABEL_TEXT_X & LABEL_TEXT_Y
  label_text_x <- label_text_xy[, "x"]
  label_text_y <- label_text_xy[, "y"]

  # REVERT LABEL_TEXT_X & LABEL_TEXT_Y TO ORIGINAL FORMAT
  if (SMP > 1) {
    label_text_x <- LAPPLY(seq_len(SMP), function(z) {
      label_text_x[names(label_text_x) == z]
    })
    label_text_y <- LAPPLY(seq_len(SMP), function(z) {
      label_text_y[names(label_text_y) == z]
    })
  }

  # RETURN LABEL CO-ORDINATES --------------------------------------------------

  # LABEL_COORDS MATRIX
  label_text_xy <- matrix(c(label_text_x, label_text_y),
    ncol = 2,
    byrow = FALSE
  )
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
}

## .CYTO_PLOT_LABEL ------------------------------------------------------------

#' Add labels to an existing cyto_plot
#'
#' @param x named list of cyto_plot arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label <- function(x) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # CYTO_PLOT_LABEL ARGUMENTS
  ARGS <- formalArgs("cyto_plot_labeller")

  # LABEL CONSTRUCTION ---------------------------------------------------------

  # CALL CYTO_PLOT_LABELLER
  label_text_xy <- do.call("cyto_plot_labeller", x[ARGS])

  # RETURN LABEL CO-ORDINATES
  return(label_text_xy)
}
