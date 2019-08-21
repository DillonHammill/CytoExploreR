# CYTO_PLOT_GATE ---------------------------------------------------------------

# NOTE: Gates for stacked distributions must be modified prior to cyto_plot_gate
# to include "y" co-ordinates.

#' cyto_plot_gate2
#' @return add gates and labels to \code{cyto_plot} and invisibly return label
#'   co-ordinates.
#' @export
cyto_plot_gate2 <- function(x,
                            channels,
                            overlay = NA,
                            gate = NA,
                            label = TRUE,
                            label_text = NA,
                            label_stat = "freq",
                            label_text_x = NA,
                            label_text_y = NA,
                            label_position = "auto",
                            trans = NA,
                            negate = FALSE,
                            display = 1,
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
                            density_smooth = 0.6,
                            density_stack = 0) {

  # CHECKS ---------------------------------------------------------------------

  # x must be a flowFrame
  if (!inherits(x, "flowFrame")) {
    stop("'x' must be a flowFrame object.")
  }
  
  # PREPARE DATA ---------------------------------------------------------------

  # OVERLAY - LIST OF FLOWFRAMES
  if (.all_na(overlay)) {
    overlay <- cyto_convert(overlay, "list of flowFrames")
  }

  # LIST OF FLOWFRAMES
  if (!.all_na(overlay)) {
    fr_list <- list(list(x), overlay)
  } else {
    fr_list <- list(x)
  }

  # SAMPLES
  SMP <- length(fr_list)

  # SAMPLING -------------------------------------------------------------------

  # Sample x
  if (display != 1) {
    fr_list <- cyto_sample(fr_list, display = display, seed = 56)
  }

  # PREPARE GATES --------------------------------------------------------------

  # 1D GATE - RECTANGLE ONLY
  if (length(channels) == 1 & .all_na(gate)) {

    # PREPARE GATES - 1D NO STACK
    if (density_stack == 0) {
      # LIST OF GATES
      if (class(gate) %in% c(
        "rectangleGate",
        "polygonGate",
        "ellipsoidGate"
      )) {
        gate <- list(gate)
      } else if (class(gate) == "filters") {
        gate <- unlist(gate)
      } else if (class(gate) == "list" &
        all(LAPPLY(gate, "class") %in% c(
          "rectangleGate",
          "polygonGate",
          "ellipsoidGate",
          "filters"
        ))) {
        gate <- unlist(gate)
      } else {
        stop("'gate' must be a list of gate objects.")
      }

      # DIMENSIONS - RECTANGLE GATES
      gate <- cyto_gate_convert(gate, channels = channels)

      # PREPARE GATES - 1D STACK
    } else if (density_stack != 0) {
      # LIST OF GATE LISTS
      if (class(gate) %in% c(
        "rectangleGate",
        "polygonGate",
        "ellipsoidGate"
      )) {
        gate <- rep(list(gate), SMP)
      } else if (class(gate) == "filters") {
        gate <- rep(unlist(gate), SMP)
      } else if (class(gate) == "list" &
        all(LAPPLY(gate, "class") %in% c(
          "rectangleGate",
          "polygonGate",
          "ellipsoidGate",
          "filters"
        ))) {
        gate <- rep(unlist(gate), SMP)
      } else {
        stop("'gate' must be a list of gate lists.")
      }

      # DIMENSIONS - RECTANGLE GATES
      gate <- lapply(gate, function(z) {
        cyto_gate_convert(z, chnanels - channels)
      })
    }

    # 2D GATE - RECTANGLE, POLYGON OR ELLIPSE
  } else if (length(channels) == 2 & .all_na(gate)) {

    # GATES - LIST OF GATES
    if (class(gate) %in% c(
      "rectangleGate",
      "polygonGate",
      "ellipsoidGate"
    )) {
      gate <- list(gate)
    } else if (class(gate) == "filters") {
      gate <- unlist(gate)
    } else if (class(gate) == "list" &
      all(LAPPLY(gate, "class") %in% c(
        "rectangleGate",
        "polygonGate",
        "ellipsoidGate",
        "filters"
      ))) {
      gate <- unlist(gate)
    } else {
      stop("'gate' must be a list of gate objects.")
    }

    # DIMENSIONS - 2D GATES - RECTANGLE/POLYGON/ELLIPSE
    gate <- cyto_gate_convert(gate, channels = channels)
  }

  # ARGUMENTS ------------------------------------------------------------------

  # ARGUMENTS
  args <- .args_list()

  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)

  # GATES #
  if (length(channels) == 1 & density_stack != 0) {
    gates <- length(gate)
  } else {
    gates <- length(gate[[1]])
  }

  # SPLIT ARGUMENTS - GATE & LABEL ARGUMENTS
  .cyto_plot_args_split(args[c(
    "negate",
    "gate_line_type",
    "gate_line_width",
    "gate_line_col",
    "gate_fill",
    "gate_fill_alpha",
    "label_text",
    "label_stat",
    "label_text_x",
    "label_text_y",
    "label_text_font",
    "label_text_size",
    "label_text_col",
    "label_fill",
    "label_fill_alpha"
  )],
  channels = channels,
  n = length(fr_list),
  plots = 1,
  layers = length(fr_list),
  gates = gates
  )

  # UPDATE ARGUMENTS
  .args_update(args)

  # PLOT GATES -----------------------------------------------------------------

  # ADD GATES TO PLOT
  if (plot == TRUE) {

    # 1D GATES
    if (length(channels) == 1 & !.all_na(gate)) {

      # UNSTACKED
      if (density_stack == 0) {

        # RECTANGLE
        mapply(
          function(gate,
                             gate_line_type,
                             gate_line_width,
                             gate_line_col,
                             gate_fill,
                             gate_fill_alpha) {
            rect(
              xleft = gate@min[channels[1]],
              ybottom = 0.6 * par("usr")[3],
              xright = gate@max[channels[1]],
              ytop = 0.985 * par("usr")[4],
              border = gate_line_col,
              lwd = gate_line_width,
              lty = gate_line_type,
              col = adjustcolor(gate_fill, gate_fill_alpha)
            )
          }, gate,
          gate_line_type,
          gate_line_width,
          gate_line_col,
          gate_fill,
          gate_fill_alpha
        )

        # STACKED
      } else if (density_stack != 0) {

        # RECTANGLE
        mapply(
          function(gate,
                             gate_line_type,
                             gate_line_width,
                             gate_line_col,
                             gate_fill,
                             gate_fill_alpha) {
            rect(
              xleft = gate@min[channels[1]],
              ybottom = 0.6 * par("usr")[3],
              xright = gate@max[channels[1]],
              ytop = 0.985 * par("usr")[4],
              border = gate_line_col,
              lwd = gate_line_width,
              lty = gate_line_type,
              col = adjustcolor(gate_fill, gate_fill_alpha)
            )
          }, gate[[1]],
          gate_line_type,
          gate_line_width,
          gate_line_col,
          gate_fill,
          gate_fill_alpha
        )
      }
      # 2D GATES
    } else if (length(channels) == 2 & !.all_na(gate)) {

      # RECTANGLE/POLYGON/ELLIPSE
      mapply(
        function(gate,
                         gate_line_type,
                         gate_line_width,
                         gate_line_col,
                         gate_fill,
                         gate_fill_alpha) {

          # ELLIPSE -> POLYGON
          if (class(gate) == "ellipsoidGate") {
            gate <- as(gate, "polygonGate")
          }

          # rectangleGate
          if (class(gate) == "rectangleGate") {
            rect(
              xleft = gate@min[channels[1]],
              ybottom = gate@min[channels[2]],
              xright = gate@max[channels[1]],
              ytop = gate@max[channels[2]],
              border = gate_line_col,
              lwd = gate_line_width,
              lty = gate_line_type,
              col = adjustcolor(gate_fill, gate_fill_alpha)
            )
            # polygonGate
          } else if (class(gate) == "polygonGate") {
            polygon(gate@boundaries[, channels[1]],
              gate@boundaries[, channels[2]],
              border = gate_line_col,
              lwd = gate_line_width,
              lty = gate_line_type,
              col = adjustcolor(gate_fill, gate_fill_alpha)
            )
          }
        }, gate[[1]],
        gate_line_type,
        gate_line_width,
        gate_line_col,
        gate_fill,
        gate_fill_alpha
      )
    }
  }

  # LABELS ---------------------------------------------------------------------
  
  # 1D STACKED - ALL LAYERS
  if(length(channels) == 1 & density_stack != 0 & label == TRUE){
    # NUMBER OF GATES - SEQ
    n <- seq_len(length(gate[[1]]))
    # LABLE PER LAYER
    lapply(fr_list, function(y){
      # STATISTICS - GATED
      fr_stats <- LAPPLY(n, function(z){
        # GATE
        fr <- Subset(y, gate[[z]])
        # STATISTIC
        if(label_stat[z] == "freq"){
          fr_stat <- .cyto_count(fr) / .cyto_count(y) * 100
          fr_stat <- round(fr_stat, 2)
        }else{
          fr_stat <- cyto_stats_compute(fr,
                                        channels = channels,
                                        trans = trans,
                                        stat = label_stat[z],
                                        gate = NA,
                                        format = "long",
                                        density_smooth = density_smooth)
          fr_stat <- as.numeric(fr_stat[1, ncol(fr_stat)])
          fr_stat <- LAPPLY(fr_stat, function(z) {.round(z, 2)})
        }
      })
      
      # STATISTICS - NEGATED
      if(negate == TRUE){
        # NEGATED POPULATION
        gate_filter <- do.call("|", gate)
        fr_gated_list <- flowCore::split(y, gate_filter)
        fr <- fr_gated_list[[1]]
        fr_negated <- fr_gated_list[[2]]
        # STATISTIC
        if(label_stat[length(label_stat)] == "freq"){
          fr_negated_stats <- 100 - sum(fr_stats)
        }else{
          fr_negated_stat <- cyto_stats_compute(fr_negated,
                                                channels = channels,
                                                trans = trans,
                                                stat = label_stat[length(label_stat)],
                                                gate = NA,
                                                format = "long",
                                                density_smooth = density_smooth)
          fr_negated_stat <- as.numeric(fr_negated_stat[1, ncol(fr_negated_stat)])
          fr_negated_stat <- .round(fr_negated_stat, 2)
        }
        fr_stats <- c(fr_stats, fr_negated_stats)
      }
      
      # LABEL_TEXT
      lapply(seq_len(length(label_text)), function(z){
        # LABEL TEXT + STAT
        if(!.all_na(label_text[z])){
          if(label_stat[z] == "freq"){
            label_text[z] <<- paste(label_text[z], "\n",
                                    paste(fr_stats[z], "%"))
          }else{
            label_text[z] <<- paste(label_text[z], "\n", fr_stats[z])
          }
          # LABEL STAT ONLY
        }else{
          if(label_stat[z] == "freq"){
            label_text[z] <<- paste(fr_stats[z], "%")
          }else{
            label_text[z] <<- fr_stats[z]
          }
        }
      })
      
      # LABEL CO-ORDINATES - AUTO/MANUAL (NOT BOTH)
      if(length(gate) > 1 & label_position == "auto"){
        # LABEL CO-ORDINATES -> NA
        label_text_x[n] <- rep(NA, n)
        label_text_y[n] <- rep(NA, n)
        # LABEL OFFSET
        label_text_xy <- .cyto_plot_label_offset(x = gate,
                                                 channels = channels,
                                                 text = label_text[n],
                                                 stat = label_stat[n],
                                                 text_x = label_text_x[n],
                                                 text_y = label_text_y[n],
                                                 text_size = label_text_size[n])
        # UPDATE LABEL CO-ORDINATES
        label_text_x[n] <- label_text_xy[1,]
        label_text_y[n] <- label_text_xy[2,]
      }
    })
    
  # 1D STACKED & 2D - BASE LAYER
  }else if(label == TRUE){
    # NUMBER OF GATES - SEQ
    n <- seq_len(length(gate))
    # STATISTICS - GATED
    fr_stats <- LAPPLY(seq_len(length(gate)), function(z){
      # GATE
      fr <- Subset(fr_list[[1]], gate[[z]])
      # STATISTIC
      if(label_stat[z] == "freq"){
        fr_stat <- .cyto_count(fr) / .cyto_count(fr_list[[1]]) * 100
        fr_stat <- round(fr_stat, 2)
      }else{
        fr_stat <- cyto_stats_compute(fr,
                                      channels = channels,
                                      trans = trans,
                                      stat = label_stat[z],
                                      gate = NA,
                                      format = "long",
                                      density_smooth = density_smooth)
        fr_stat <- as.numeric(fr_stat[1, ncol(fr_stat)])
        fr_stat <- LAPPLY(fr_stat, function(z) {.round(z, 2)})
      }
    })

    # STATISTICS - NEGATED
    if(negate == TRUE){
      # NEGATED POPULATION
      gate_filter <- do.call("|", gate)
      fr_gated_list <- flowCore::split(fs_list[[1]], gate_filter)
      fr <- fr_gated_list[[1]]
      fr_negated <- fr_gated_list[[2]]
      # STATISTIC
      if(label_stat[length(label_stat)] == "freq"){
        fr_negated_stats <- 100 - sum(fr_stats)
      }else{
        fr_negated_stat <- cyto_stats_compute(fr_negated,
                                              channels = channels,
                                              trans = trans,
                                              stat = label_stat[length(label_stat)],
                                              gate = NA,
                                              format = "long",
                                              density_smooth = density_smooth)
        fr_negated_stat <- as.numeric(fr_negated_stat[1, ncol(fr_negated_stat)])
        fr_negated_stat <- .round(fr_negated_stat, 2)
      }
      fr_stats <- c(fr_stats, fr_negated_stats)
    }
    
    # LABEL_TEXT
    lapply(seq_len(length(label_text)), function(z){
      # LABEL TEXT + STAT
      if(!.all_na(label_text[z])){
        if(label_stat[z] == "freq"){
          label_text[z] <<- paste(label_text[z], "\n",
                                  paste(fr_stats[z], "%"))
        }else{
          label_text[z] <<- paste(label_text[z], "\n", fr_stats[z])
        }
      # LABEL STAT ONLY
      }else{
        if(label_stat[z] == "freq"){
          label_text[z] <<- paste(fr_stats[z], "%")
        }else{
          label_text[z] <<- fr_stats[z]
        }
      }
    })
    
    # LABEL CO-ORDINATES - AUTO
    if(length(gate) > 1 & label_position == "auto"){
      # LABEL CO-ORDINATES -> NA
      label_text_x[n] <- rep(NA, n)
      label_text_y[n] <- rep(NA, n)
      # LABEL OFFSET
      n <- seq_len(length(gate))
      label_text_xy <- .cyto_plot_label_offset(x = gate,
                                               channels = channels,
                                               text = label_text[n],
                                               stat = label_stat[n],
                                               text_x = label_text_x[n],
                                               text_y = label_text_y[n],
                                               text_size = label_text_size[n])
      # UPDATE LABEL CO-ORDINATES
      label_text_x[seq_len(length(gate))] <- label_text_xy[1,]
      label_text_y[seq_len(length(gate))] <- label_text_xy[2,]
      
    }
    
  }

  # LABEL PLOT -----------------------------------------------------------------
  
  label_text_xy <- lapply(seq_len(label_text), function(z){
    # MESSAGE NEGATED LABEL
    if(z == length(label_text) & 
       negate == TRUE & 
       .all_na(c(label_text_x[z], label_text_y[z]))){
      message("Select a label location for the negated population.")
    }
    suppressMessages(cyto_plot_label2(label_text = label_text[z],
                                      label_text_x = label_text_x[z],
                                      label_text_y = label_text_y[z],
                                      label_text_font = label_text_font[z],
                                      label_text_size = label_text_size[z],
                                      label_text_col = label_text_col[z],
                                      label_fill = label_fill[z],
                                      label_fill_alpha = label_fill_alpha[z]))
  })
  
  # RETURN ---------------------------------------------------------------------
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)

}
