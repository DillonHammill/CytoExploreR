#' Plot gates
#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2 <- function(x,
                            channels,
                            gate = NA, ...) {
  UseMethod("cyto_plot_gate2", gate)
}

#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2.default <- function(x,
                                    channels,
                                    gate = NA,
                                    label = TRUE,
                                    label_text,
                                    label_stat = "freq",
                                    label_text_x,
                                    label_text_y,
                                    label_position = "auto",
                                    trans = NA,
                                    negate = FALSE,
                                    gate_line_type = 1,
                                    gate_line_width = 2.5,
                                    gate_line_col = "red",
                                    gate_fill = "white",
                                    gate_fill_alpha = 0,
                                    label_text_font,
                                    label_text_size,
                                    label_text_col,
                                    label_fill = "white",
                                    label_fill_alpha = 0.6,
                                    density_smooth = 0.6) {

  # No gate supplied
  if (.all_na(gate)) {
    stop("No gate object(s) supplied for plotting.")
  }
}

#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2.rectangleGate <- function(x,
                                          channels,
                                          gate = NA,
                                          label = TRUE,
                                          label_text,
                                          label_stat = "freq",
                                          label_text_x = NA,
                                          label_text_y = NA,
                                          label_position = "auto",
                                          trans = NA,
                                          negate = FALSE,
                                          gate_line_type = 1,
                                          gate_line_width = 2.5,
                                          gate_line_col = "red",
                                          gate_fill = "white",
                                          gate_fill_alpha = 0,
                                          label_text_font = 2,
                                          label_text_size = 0.8,
                                          label_text_col = "black",
                                          label_fill = "white",
                                          label_fill_alpha = 0.6,
                                          density_smooth = 0.6) {

  # x must be a flowFrame
  if (!inherits(x, "flowFrame")) {
    stop("'x' must be a flowFrame object.")
  }

  # Allow 1D gate plotted in 2D
  if (!missing(channels)) {
    if (length(channels) == 2 & length(parameters(gate)) == 1) {
      rg <- matrix(c(
        as.numeric(gate@min),
        as.numeric(gate@max),
        -Inf, Inf
      ),
      ncol = 2,
      nrow = 2
      )
      colnames(rg) <- c(
        as.vector(parameters(gate)),
        channels[!channels == as.vector(parameters(gate))]
      )
      rownames(rg) <- c("min", "max")
      gate <- rectangleGate(.gate = rg)
    }
  } else if (missing(channels)) {
    channels <- as.vector(parameters(gate))
  }

  # Allow 2D gate plotted in 1D
  if (length(channels) == 1 &
    length(parameters(gate)) == 2) {
    gate <- gate[channels]
  }

  # Gate parameters must match channels for plotting
  if (!all(as.vector(parameters(gate)) %in% channels)) {
    stop(paste(
      "Channels argument should contain",
      paste(parameters(gate), collapse = " & "),
      "to plot supplied gate(s)."
    ))
  }

  # 2D rectangle Gate
  if (length(gate@min) == 2) {

    # Replace -Inf x values for plotting
    if (is.infinite(gate@min[channels[1]]) |
      gate@min[channels[1]] < (par("usr")[1] - 0.13 * par("usr")[1])) {
      gate@min[channels[1]] <- (par("usr")[1] - 0.13 * par("usr")[1])
    }

    # Replace Inf x values for plotting
    if (is.infinite(gate@max[channels[1]]) |
      gate@max[channels[1]] > 0.98 * par("usr")[2]) {
      gate@max[channels[1]] <- 0.98 * par("usr")[2]
    }

    # Replace -Inf y values for plotting
    if (is.infinite(gate@min[channels[2]]) |
      gate@min[channels[2]] < (par("usr")[3] - 0.13 * par("usr")[3])) {
      gate@min[channels[2]] <- (par("usr")[3] - 0.13 * par("usr")[3])
    }

    # Replace Inf y values for plotting
    if (is.infinite(gate@max[channels[2]]) |
      gate@max[channels[2]] > 0.98 * par("usr")[4]) {
      gate@max[channels[2]] <- 0.98 * par("usr")[4]
    }

    rect(
      xleft = gate@min[channels[1]],
      ybottom = gate@min[channels[2]],
      xright = gate@max[channels[1]],
      ytop = gate@max[channels[2]],
      border = gate_line_col,
      lwd = gate_line_width,
      lty = gate_line_type
    )
  } else if (length(gate@min) == 1) {

    # Replace -Inf values for plotting
    if (is.infinite(gate@min[1]) |
      gate@min[1] < (par("usr")[1] - 0.13 * par("usr")[1])) {
      gate@min[1] <- (par("usr")[1] - 0.13 * par("usr")[1])
    }

    # Replace Inf values for plotting
    if (is.infinite(gate@max[1]) |
      gate@max[1] > 0.98 * par("usr")[2]) {
      gate@max[1] <- 0.98 * par("usr")[2]
    }

    # Add rectangle to plot
    rect(
      xleft = gate@min,
      ybottom = 0.6 * par("usr")[3],
      xright = gate@max,
      ytop = 0.985 * par("usr")[4],
      border = gate_line_col,
      lwd = gate_line_width,
      lty = gate_line_type,
      col = adjustcolor(gate_fill, gate_fill_alpha)
    )
  }

  # Add label to gate?
  if (label == TRUE) {

    # Gate flowFrame
    if (negate == TRUE) {
      fr_list <- split(x, gate)
      fr <- fr_list[[1]]
      fr_negated <- fr_list[[2]]
    } else {
      fr <- Subset(x, gate)
    }

    # Need to add extra label if negate == TRUE
    if (negate == TRUE) {
      label_text <- rep(label_text, length.out = 2)
      label_stat <- rep(label_stat, length.out = 2)
      label_text_x <- rep(label_text_x, length.out = 2)
      label_text_y <- rep(label_text_y, length.out = 2)
      label_text_font <- rep(label_text_font, length.out = 2)
      label_text_size <- rep(label_text_size, length.out = 2)
      label_text_col <- rep(label_text_col, length.out = 2)
      label_fill <- rep(label_fill, length.out = 2)
      label_fill_alpha <- rep(label_fill_alpha, length.out = 2)
    }

    # Calculate statistics for gated population
    if (!.all_na(label_stat[1])) {
      # Unsupported statistics in 2D
      if (length(channels) == 2) {
        if (!label_stat[1] %in% c("count", "freq")) {
          stop(
            "Only count and frequency statistics are supported in 2D plots."
          )
        }
        if (label_stat[1] == "freq" & .all_na(gate)) {
          stop("Supply a 'gate' object to compute frequency.")
        }
      }
      # Statistic for gated population
      fr_stat <- cyto_stats_compute(x,
        channels = channels,
        trans = trans,
        stat = label_stat[1],
        gate = gate,
        format = "long",
        density_smooth = density_smooth
      )
      fr_stat <- as.numeric(fr_stat[1, ncol(fr_stat)])
      # Combine label_text with statistic
      if (!.all_na(label_text[1])) {
        if (label_stat[1] == "freq") {
          label_text[1] <- paste(label_text[1], "\n", paste(fr_stat, "%"))
        } else {
          label_text[1] <- paste(label_text[1], "\n", fr_stat)
        }
      } else {
        if (label_stat[1] == "freq") {
          label_text[1] <- paste(fr_stat, "%")
        } else {
          label_text[1] <- fr_stat
        }
      }
    }

    # Sort out co-ordinate for gated label
    if (label_position == "auto") {

      # Calculate x co-ordinate for label
      if (.all_na(label_text_x[1])) {
        label_text_x[1] <- sum(gate@min[channels[1]], gate@max[channels[1]]) / 2
      }

      # Calculate y co-ordinate for label
      if (.all_na(label_text_y[1])) {
        label_text_y[1] <- sum(gate@min[channels[2]], gate@max[channels[2]]) / 2
      }

      # Interactively select co-ordinatefor label
    } else if (label_position == "manual") {

      # Interactively select
      if (any(.all_na(c(label_text_x[1], label_text_y[1])))) {
        message(
          paste("Select a label location for the", gate@filterId, "gate:")
        )
        label_text_xy <- locator(n = 1)
        label_text_x[1] <- label_text_xy[[1]]
        label_text_y[1] <- label_text_xy[[2]]
      }
    }

    # Add label to plot
    cyto_plot_label2(
      label_text = label_text[1],
      label_text_x = label_text_x[1],
      label_text_y = label_text_y[1],
      label_text_font = label_text_font[1],
      label_text_size = label_text_size[1],
      label_text_col = label_text_col[1],
      label_fill = label_fill[1],
      label_fill_alpha = label_fill_alpha[1]
    )

    # Calculate statistics for negated population
    if (negate == TRUE) {
      # Calculate statistics for gated population
      if (!.all_na(label_stat[2])) {
        # Unsupported statistics in 2D
        if (length(channels) == 2) {
          if (!label_stat[2] %in% c("count", "freq")) {
            stop(
              "Only count and frequency statistics are supported in 2D plots."
            )
          }
          if (label_stat[2] == "freq" & .all_na(gate)) {
            stop("Supply a 'gate' object to compute frequency.")
          }
        }
        # Statistic for gated population
        if (label_stat[2] == "freq") {
          fr_negated_stat <- 100 - as.numeric(fr_stat)
        } else {
          fr_negated_stat <- cyto_stats_compute(fr_negated,
            channels = channels,
            trans = trans,
            stat = label_stat[2],
            gate = NA,
            format = "long",
            density_smooth = density_smooth
          )
          fr_negated_stat <- as.numeric(fr_negated_stat[1, ncol(fr_negated_stat)])
        }

        # Add % sign to freq statistic
        if (label_stat[2] == "freq") {
          fr_stat <- paste(fr_stat, "%")
        }
        # Combine label_text with statistic
        if (!.all_na(label_text[2])) {
          if (label_stat[2] == "freq") {
            label_text[2] <- paste(
              label_text[2], "\n",
              paste(fr_negated_stat, "%")
            )
          } else {
            label_text[2] <- paste(label_text[2], "\n", fr_negated_stat)
          }
          # No label text - statistic only
        } else {
          if (label_stat[2] == "freq") {
            label_text[2] <- paste(fr_negated_stat, "%")
          } else {
            label_text[2] <- fr_negated_stat
          }
        }
      }
      # Add label to plot
      cyto_plot_label2(
        label_text = label_text[2],
        label_text_x = label_text_x[2],
        label_text_y = label_text_y[2],
        label_text_font = label_text_font[2],
        label_text_size = label_text_size[2],
        label_text_col = label_text_col[2],
        label_fill = label_fill[2],
        label_fill_alpha = label_fill_alpha[2]
      )
    }
  }

  # Return modified gate
  invisible(gate)
}

#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2.polygonGate <- function(x,
                                        channels,
                                        gate = NA,
                                        label = TRUE,
                                        label_text,
                                        label_stat = "freq",
                                        label_text_x = NA,
                                        label_text_y = NA,
                                        label_position = "auto",
                                        trans = NA,
                                        negate = FALSE,
                                        gate_line_type = 1,
                                        gate_line_width = 2.5,
                                        gate_line_col = "red",
                                        gate_fill = "white",
                                        gate_fill_alpha = 0,
                                        label_text_font = 2,
                                        label_text_size = 0.8,
                                        label_text_col = "black",
                                        label_fill = "white",
                                        label_fill_alpha = 0.6,
                                        density_smooth = 0.6) {

  # x must be a flowFrame
  if (!inherits(x, "flowFrame")) {
    stop("'x' must be a flowFrame object.")
  }

  # Get channels from gate
  if (missing(channels)) {
    channels <- as.vector(parameters(gate))
  }

  # Allow plotting ploygonGate in 1D plot - use min/max coords in channel
  if (length(channels) == 1) {

    # Gate not constructed using supplied channel
    if (!any(as.vector(parameters(gate)) %in% channels)) {
      stop(
        paste("Supplied gate does not have co-ordinates in", channels, ".")
      )
    } else {
      coords <- range(gate@boundaries[, parameters(gate) %in% channels])
      coords <- matrix(coords,
        dimnames = list(c("min", "max"), channels),
        ncol = 1,
        nrow = 2
      )
      gate <- rectangleGate(filterId = gate@filterId, .gate = coords)

      # Call to rectangleGate method
      cyto_plot_gate2(x,
        channels = channels,
        gate = gate,
        label = label,
        label_text = label_text,
        label_stat = label_stat,
        label_text_x = label_text_x,
        label_text_y = label_text_y,
        label_position = label_position,
        trans = trans,
        negate = negate,
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
  } else {

    # Parameters of gate must match channels of plot
    if (!all(as.vector(parameters(gate)) %in% channels)) {
      stop(paste(
        "Channels argument should contain",
        paste(parameters(gate), collapse = " & "),
        "to plot supplied gate(s)."
      ))
    }

    # Replace Inf values with plot limits
    if (!all(is.finite(gate@boundaries))) {
      cnt <- 0
      lapply(seq_along(channels), function(x) {
        cnt <<- cnt + 1

        if (any(!is.finite(gate@boundaries[, channels[x]]) &
          any(gate@boundaries[, channels[x]] < 0))) {
          if (cnt == 1) {
            ind <- which(gate@boundaries[, channels[x]] < 0)
            gate@boundaries[, channels[x]][ind] <<- par("usr")[1]
          } else if (cnt == 2) {
            ind <- which(gate@boundaries[, channels[x]] < 0)
            gate@boundaries[, channels[x]][ind] <<- par("usr")[3]
          }
        }

        if (any(!is.finite(gate@boundaries[, channels[x]]) &
          any(!gate@boundaries[, channels[x]] < 0))) {
          if (cnt == 1) {
            ind <- which(!is.finite(gate@boundaries[, channels[x]]) &
              !gate@boundaries[, channels[x]] < 0)
            gate@boundaries[, channels[x]][ind] <<- par("usr")[2]
          } else if (cnt == 2) {
            ind <- which(!is.finite(gate@boundaries[, channels[x]]) &
              !gate@boundaries[, channels[x]] < 0)
            gate@boundaries[, channels[x]][ind] <<- par("usr")[4]
          }
        }
      })
    }

    polygon(gate@boundaries[, channels[1]],
      gate@boundaries[, channels[2]],
      border = gate_line_col,
      lwd = gate_line_width,
      lty = gate_line_type
    )
  }

  # Add label to gate?
  if (label == TRUE) {

    # Gate flowFrame
    if (negate == TRUE) {
      fr_list <- split(x, gate)
      fr <- fr_list[[1]]
      fr_negated <- fr_list[[2]]
    } else {
      fr <- Subset(x, gate)
    }

    # Need to add extra label if negate == TRUE
    if (negate == TRUE) {
      label_text <- rep(label_text, length.out = 2)
      label_stat <- rep(label_stat, length.out = 2)
      label_text_x <- rep(label_text_x, length.out = 2)
      label_text_y <- rep(label_text_y, length.out = 2)
      label_text_font <- rep(label_text_font, length.out = 2)
      label_text_size <- rep(label_text_size, length.out = 2)
      label_text_col <- rep(label_text_col, length.out = 2)
      label_fill <- rep(label_fill, length.out = 2)
      label_fill_alpha <- rep(label_fill_alpha, length.out = 2)
    }
    # Calculate statistics for gated population
    if (!.all_na(label_stat[1])) {
      # Unsupported statistics in 2D
      if (length(channels) == 2) {
        if (!label_stat[1] %in% c("count", "freq")) {
          stop(
            "Only count and frequency statistics are supported in 2D plots."
          )
        }
        if (label_stat[1] == "freq" & .all_na(gate)) {
          stop("Supply a 'gate' object to compute frequency.")
        }
      }
      # Statistic for gated population
      fr_stat <- cyto_stats_compute(x,
        channels = channels,
        trans = trans,
        stat = label_stat[1],
        gate = gate,
        format = "long",
        density_smooth = density_smooth
      )
      fr_stat <- as.numeric(fr_stat[1, ncol(fr_stat)])
      # Combine label_text with statistic
      if (!.all_na(label_text[1])) {
        if (label_stat[1] == "freq") {
          label_text[1] <- paste(label_text[1], "\n", paste(fr_stat, "%"))
        } else {
          label_text[1] <- paste(label_text[1], "\n", fr_stat)
        }
      } else {
        if (label_stat[1] == "freq") {
          label_text[1] <- paste(fr_stat, "%")
        } else {
          label_text[1] <- fr_stat
        }
      }
    }

    # Sort out co-ordinate for gated label
    if (label_position == "auto") {

      # Calculate x co-ordinate for label
      if (.all_na(label_text_x[1])) {
        label_text_x[1] <- sum(gate@boundaries[, channels[1]]) /
          nrow(gate@boundaries)
      }

      # Calculate y co-ordinate for label
      if (.all_na(label_text_y[1])) {
        label_text_y[1] <- sum(gate@boundaries[, channels[2]]) /
          nrow(gate@boundaries)
      }

      # Interactively select co-ordinatefor label
    } else if (label_position == "manual") {

      # Interactively select
      if (any(.all_na(c(label_text_x[1], label_text_y[1])))) {
        message(
          paste("Select a label location for the", gate@filterId, "gate:")
        )
        label_text_xy <- locator(n = 1)
        label_text_x[1] <- label_text_xy[[1]]
        label_text_y[1] <- label_text_xy[[2]]
      }
    }

    # Add label to plot
    cyto_plot_label2(
      label_text = label_text[1],
      label_text_x = label_text_x[1],
      label_text_y = label_text_y[1],
      label_text_font = label_text_font[1],
      label_text_size = label_text_size[1],
      label_text_col = label_text_col[1],
      label_fill = label_fill[1],
      label_fill_alpha = label_fill_alpha[1]
    )

    # Calculate statistics for negated population
    if (negate == TRUE) {
      # Calculate statistics for gated population
      if (!.all_na(label_stat[2])) {
        # Unsupported statistics in 2D
        if (length(channels) == 2) {
          if (!label_stat[2] %in% c("count", "freq")) {
            stop(
              "Only count and frequency statistics are supported in 2D plots."
            )
          }
          if (label_stat[2] == "freq" & .all_na(gate)) {
            stop("Supply a 'gate' object to compute frequency.")
          }
        }
        # Statistic for gated population
        if (label_stat[2] == "freq") {
          fr_negated_stat <- 100 - as.numeric(fr_stat)
        } else {
          fr_negated_stat <- cyto_stats_compute(fr_negated,
            channels = channels,
            trans = trans,
            stat = label_stat[2],
            gate = NA,
            format = "long",
            density_smooth = density_smooth
          )
          fr_negated_stat <- as.numeric(fr_negated_stat[1, ncol(fr_negated_stat)])
        }

        # Add % sign to freq statistic
        if (label_stat[2] == "freq") {
          fr_stat <- paste(fr_stat, "%")
        }
        # Combine label_text with statistic
        if (!.all_na(label_text[2])) {
          if (label_stat[2] == "freq") {
            label_text[2] <- paste(
              label_text[2], "\n",
              paste(fr_negated_stat, "%")
            )
          } else {
            label_text[2] <- paste(label_text[2], "\n", fr_negated_stat)
          }
          # No label text - statistic only
        } else {
          if (label_stat[2] == "freq") {
            label_text[2] <- paste(fr_negated_stat, "%")
          } else {
            label_text[2] <- fr_negated_stat
          }
        }
      }
      # Add label to plot
      cyto_plot_label2(
        label_text = label_text[2],
        label_text_x = label_text_x[2],
        label_text_y = label_text_y[2],
        label_text_font = label_text_font[2],
        label_text_size = label_text_size[2],
        label_text_col = label_text_col[2],
        label_fill = label_fill[2],
        label_fill_alpha = label_fill_alpha[2]
      )
    }
  }

  # Return modified gate
  invisible(gate)
}

#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2.ellipsoidGate <- function(x,
                                          channels,
                                          gate = NA,
                                          label = TRUE,
                                          label_text,
                                          label_stat = "freq",
                                          label_text_x = NA,
                                          label_text_y = NA,
                                          label_position = "auto",
                                          trans = NA,
                                          negate = FALSE,
                                          gate_line_type = 1,
                                          gate_line_width = 2.5,
                                          gate_line_col = "red",
                                          gate_fill = "white",
                                          gate_fill_alpha = 0,
                                          label_text_font = 2,
                                          label_text_size = 0.8,
                                          label_text_col = "black",
                                          label_fill = "white",
                                          label_fill_alpha = 0.6,
                                          density_smooth = 0.6) {

  # Check Channels
  if (missing(channels)) {
    channels <- as.vector(parameters(gate))
  }

  # Allow plotting ploygonGate in 1D plot - use min/max coords in channel
  if (length(channels) == 1) {

    # Gate not constructed using supplied channel
    if (!any(as.vector(parameters(gate)) %in% channels)) {
      stop(
        paste("Supplied gate does not have co-ordinates in", channels, ".")
      )
    } else {
      # convert ellipsoidGate into polygonGate
      gate <- as(gate, "polygonGate")
      coords <- range(gate@boundaries[, parameters(gate) %in% channels])
      coords <- matrix(coords,
        dimnames = list(c("min", "max"), channels),
        ncol = 1,
        nrow = 2
      )
      gate <- rectangleGate(filterId = gate@filterId, .gate = coords)

      # Call to rectangleGate method
      cyto_plot_gate2(x,
        channels = channels,
        gate = gate,
        label = label,
        label_text = label_text,
        label_stat = label_stat,
        label_text_x = label_text_x,
        label_text_y = label_text_y,
        label_position = label_position,
        trans = trans,
        negate = negate,
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
  } else {
    if (!all(as.vector(parameters(gate)) %in% channels)) {
      stop(paste(
        "Channels argument should contain",
        paste(parameters(gate), collapse = " & "),
        "to plot supplied gate(s)."
      ))
    }

    # Coerce to polygonGate
    gate <- as(gate, "polygonGate")

    # Call to polygonGate method
    cyto_plot_gate2(x,
      channels = channels,
      gate = gate,
      label = label,
      label_text = label_text,
      label_stat = label_stat,
      label_text_x = label_text_x,
      label_text_y = label_text_y,
      label_position = label_position,
      trans = trans,
      negate = negate,
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

  # Return modified gate
  invisible(gate)
}

#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2.list <- function(x,
                                 channels,
                                 gate = NA,
                                 label = TRUE,
                                 label_text,
                                 label_stat = "freq",
                                 label_text_x = NA,
                                 label_text_y = NA,
                                 label_position = "auto",
                                 trans = NA,
                                 negate = FALSE,
                                 gate_line_type = 1,
                                 gate_line_width = 2.5,
                                 gate_line_col = "red",
                                 gate_fill = "white",
                                 gate_fill_alpha = 0,
                                 label_text_font = 2,
                                 label_text_size = 0.8,
                                 label_text_col = "black",
                                 label_fill = "white",
                                 label_fill_alpha = 0.6,
                                 density_smooth = 0.6) {


  # Extract gates from any filters objects
  gate <- unlist(gate)

  # Number of gates
  n <- length(gate)

  # Number of expected arguments
  if (negate == TRUE) {
    N <- length(gate) + 1
  } else {
    N <- length(gate)
  }

  # Pull down arguments
  args <- .args_list()

  # Arguments to repeat
  args_to_repeat <- names(args)[!names(args) %in% c(
    x,
    channels,
    gate,
    trans,
    negate,
    density_smooth
  )]

  # Repeat arguments
  lapply(args_to_repeat, function(z) {
    args[[z]] <<- rep(args[[z]], length.out = N)
  })

  # Update arguments
  .args_update(args)

  # Make calls to rectangleGate, polygonGate or ellipsoidGate method
  gate <- mapply(
    function(label,
                 label_text,
                 label_stat,
                 label_text_x,
                 label_text_y,
                 label_position,
                 gate_line_type,
                 gate_line_width,
                 gate_line_col,
                 gate_fill,
                 gate_fill_alpha,
                 label_text_font,
                 label_text_size,
                 label_text_col,
                 label_fill,
                 label_fill_alpha) {

      # Add gates and labels for gates - set negate to FALSE
      cyto_plot_gate2(x,
        channels = channels,
        gate = gate,
        label = label,
        label_text = label_text,
        label_stat = label_stat,
        label_text_x = label_text_x,
        label_text_y = label_text_y,
        label_position = label_position,
        trans = trans,
        negate = FALSE,
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
    }, label[seq_len(n)],
    label_text[seq_len(n)],
    label_stat[seq_len(n)],
    label_text_x[seq_len(n)],
    label_text_y[seq_len(n)],
    label_position[seq_len(n)],
    gate_line_type[seq_len(n)],
    gate_line_width[seq_len(n)],
    gate_line_col[seq_len(n)],
    gate_fill[seq_len(n)],
    gate_fill_alpha[seq_len(n)],
    label_text_font[seq_len(n)],
    label_text_size[seq_len(n)],
    label_text_col[seq_len(n)],
    label_fill[seq_len(n)],
    label_fill_alpha[seq_len(n)]
  )

  # Label for negated population
  if (label == TRUE & negate == TRUE) {

    # Get negated population
    gate_filter <- do.call("|", gate)
    fr_list <- split(x, gate_filter)[[2]]
    fr <- fr_list[[1]]
    fr_negated <- fr_list[[2]]

    # Calculate statistics for gated population
    if (!.all_na(label_stat[N])) {
      # Unsupported statistics in 2D
      if (length(channels) == 2) {
        if (!label_stat[N] %in% c("count", "freq")) {
          stop(
            "Only count and frequencystatistics are supported in 2D plots."
          )
        }
        if (label_stat[N] == "freq" & .all_na(gate)) {
          stop("Supply a 'gate' object to compute frequency.")
        }
      }
      # Statistic for gated population
      if (label_stat[N] == "freq") {
        fr_negated_stat <- .cyto_count(fr_negated) / .cyto_count(fr) * 100
        fr_negated_stat <- round(fr_negated_stat, 2)
      } else {
        fr_negated_stat <- cyto_stats_compute(fr_negated,
          channels = channels,
          trans = trans,
          stat = label_stat[N],
          gate = NA,
          format = "long",
          density_smooth = density_smooth
        )
        fr_negated_stat <- as.numeric(fr_negated_stat[1, ncol(fr_negated_stat)])
      }

      # Add % sign to freq statistic
      if (label_stat[N] == "freq") {
        fr_stat <- paste(fr_stat, "%")
      }
      # Combine label_text with statistic
      if (!.all_na(label_text[N])) {
        if (label_stat[N] == "freq") {
          label_text[N] <- paste(
            label_text[N], "\n",
            paste(fr_negated_stat, "%")
          )
        } else {
          label_text[N] <- paste(label_text[N], "\n", fr_negated_stat)
        }
        # No label text - statistic only
      } else {
        if (label_stat[N] == "freq") {
          label_text[N] <- paste(fr_negated_stat, "%")
        } else {
          label_text[N] <- fr_negated_stat
        }
      }
    }
    # Add label to plot
    cyto_plot_label2(
      label_text = label_text[N],
      label_text_x = label_text_x[N],
      label_text_y = label_text_y[N],
      label_text_font = label_text_font[N],
      label_text_size = label_text_size[N],
      label_text_col = label_text_col[N],
      label_fill = label_fill[N],
      label_fill_alpha = label_fill_alpha[N]
    )
  }

  # Return modified gates
  invisible(gate)
}

#' @rdname cyto_plot_gate2
#' @export
cyto_plot_gate2.filters <- function(x,
                                    channels,
                                    gate = NA,
                                    label = TRUE,
                                    label_text,
                                    label_stat = "freq",
                                    label_text_x = NA,
                                    label_text_y = NA,
                                    label_position = "auto",
                                    trans = NA,
                                    negate = FALSE,
                                    gate_line_type = 1,
                                    gate_line_width = 2.5,
                                    gate_line_col = "red",
                                    gate_fill = "white",
                                    gate_fill_alpha = 0,
                                    label_text_font = 2,
                                    label_text_size = 0.8,
                                    label_text_col = "black",
                                    label_fill = "white",
                                    label_fill_alpha = 0.6,
                                    density_smooth = 0.6) {

  # Convert filters object to list and call list method
  gate <- unlist(gate)

  # Call to cyto_plot_gate2 list method
  gate <- cyto_plot_gate2(x,
    channels = channels,
    gate = gate,
    label = label,
    label_text = label_text,
    label_stat = label_stat,
    label_text_x = label_text_x,
    label_text_y = label_text_y,
    label_position = label_position,
    trans = trans,
    negate = negate,
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

  # Return list of modified gates
  invisible(gate)
}
