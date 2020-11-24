## CYTO_PLOT_MAP ---------------------------------------------------------------

#' Plot maps coloured by parameter
#'
#' @param page can be either "channel" or "sample" indicating whether each page
#'   should contain for each sample or each channel.
#'
#' @importFrom grDevices n2mfrow rgb colorRamp adjustcolor
#' @importFrom graphics mtext
#' @name cyto_plot_map
NULL

#' @export
#' @noRd
cyto_plot_map <- function(x, ...){
  UseMethod("cyto_plot_map")
}

#' @rdname cyto_plot_map
#' @export
cyto_plot_map.flowSet <- function(x,
                                  map = NULL,
                                  channels = NULL,
                                  select = NULL,
                                  group_by = "name",
                                  layout = NULL,
                                  header = NULL,
                                  header_text_font = 2,
                                  header_text_size = 1,
                                  header_text_col = "black",
                                  ...){
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # SELECT
  if(!is.null(select)){
    x <- cyto_select(x, select = select)
  }
  
  # GROUP_BY
  fr_list <- cyto_group_by(x,
                           group_by = group_by)
  
  # CALIBRATION
  message(
    paste0("Remember to calibrate the colour scales using cyto_calibrate() ",
           "before calling cyto_plot_map().")
  )
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # HEADER
  if(is.null(header)){
    header <- names(fr_list)
    if(all(header == "all")){
      header <- "Combined Events"
    }
  }else{
    header <- rep(header, length.out = length(fr_list))
  }
  
  # REPEAT HEADER ARGUMENTS
  header_text_font <- rep(header_text_font, length.out = length(fr_list))
  header_text_size <- rep(header_text_size, length.out = length(fr_list))
  header_text_col <- rep(header_text_col, length.out = length(fr_list))
  
  # CALL FLOWFRAME METHOD ------------------------------------------------------
  
  print(fr_list)
  print(seq_along(fr_list))
  
  # CYTO_PLOT_MAP
  plots <- lapply(seq_along(fr_list), function(z){
    suppressMessages(
      cyto_plot_map(fr_list[[z]],
                    map = map,
                    channels = channels,
                    layout = layout,
                    header = header[z],
                    header_text_font = header_text_font[z],
                    header_text_size = header_text_size[z],
                    header_text_col = header_text_col[z],
                    point_col_scale = point_col_scale,
                    point_col_alpha = point_col_alpha,
                    ...)
    )
  })
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE
  if (getOption("cyto_plot_save")) {
    # PLOT METHOD
    if (is(x, basename(getOption("cyto_plot_method")))) {
      # CLOSE GRAPHICS DEVICE
      cyto_plot_complete()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN LIST OF RECORDED PLOTS
  invisible(plots)
  
}

#' @rdname cyto_plot_map
#' @export
cyto_plot_map.flowFrame <- function(x,
                                    map = NULL,
                                    channels = NULL,
                                    layout = NULL,
                                    header = NULL,
                                    header_text_font = 2,
                                    header_text_size = 1,
                                    header_text_col = "black",
                                    point_col_scale = NA,
                                    point_col_alpha = 1,
                                    ...){
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- .par(c("oma", "mfrow"))
  
  # THEME ----------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CHECKS ---------------------------------------------------------------------
  
  # METHOD & RESET
  if (is.null(getOption("cyto_plot_method"))) {
    # SET PLOT METHOD
    options("cyto_plot_method" = "map/flowFrame")
    # RESET PLOT METHOD & GRAPHICAL PARAMETERS ON EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_method" = NULL)
    })
  } else {
    # RESET GRAPHICAL PARAMETERS ON EXIT
    on.exit(par(old_pars))
  }
  
  # CALIBRATION 
  cyto_cal <- .cyto_calibrate_recall()
  if(is.null(cyto_cal)){
    message(
      "Please use cyto_calibrate() to properly calibrate the colour scales."
    )
  }
  
  # CHANNELS -------------------------------------------------------------------
  
  # MAP
  if(is.null(map)){
    map <- cyto_channels(x)
    map <- map[c(length(map) - 1, length(map))]
  }
  
  # DEFAULT CHANNELS WITH MARKERS
  if(is.null(channels)){
    channels <- names(cyto_markers(x))
  }
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # HEADER
  if(is.null(header)){
    header <- cyto_names(x)
  }
  
  # HEADER SPACE
  if (!is.na(header)) {
    par(oma = c(0, 0, 3, 0))
  }
  
  # LEGEND
  legend_x <- NULL
  legend_y <- NULL
  
  # PLOT LAYOUT
  if(is.null(layout)){
    dims <- n2mfrow(length(channels))
    layout <- c(dims[2], dims[1])
  }
  
  # SET LAYOUT
  par("mfrow" = layout)
  
  # PLOTS
  plots <- lapply(seq_along(channels), function(z){
    
    # CYTO_PLOT ----------------------------------------------------------------
    cyto_plot(x,
              channels = map,
              overlay = NA,
              title = channels[z],
              point_col = channels[z],
              legend = FALSE,
              margins = c(5.1, 5.1, 4.1, 2.5),
              ...)
    
    # LEGEND -------------------------------------------------------------------
    
    # LEGEND X COORDS
    if(is.null(legend_x)){
      legend_x <- c(par("usr")[2] + 0.01 * (par("usr")[2] - par("usr")[1]),
                    par("usr")[2] + 0.04 * (par("usr")[2] - par("usr")[1]))
    }
    
    # LEGEND Y COORDS
    if(is.null(legend_y)){
      legend_y <- c(par("usr")[3], par("usr")[4])
    }
    
    # LEGEND COLOUR SCALE
    legend_col_scale <- .cyto_plot_point_col_scale(point_col_scale)
    legend_col_ramp <- colorRamp(legend_col_scale)
    legend_cols <- seq(0, 1, 1/50) # 50 boxes
    legend_cols <- legend_col_ramp(legend_cols)
    legend_cols <- rgb(legend_cols[, 1],
                       legend_cols[, 2],
                       legend_cols[, 3],
                       maxColorValue = 255)
    legend_cols <- adjustcolor(legend_cols, point_col_alpha)
    
    # LEGEND BORDER
    rect(legend_x[1],
         legend_y[1],
         legend_x[2],
         legend_y[2],
         xpd = TRUE)
    
    # LEGEND BOXES
    legend_box_x <- legend_x
    legend_box_y <- seq(par("usr")[3], 
                        par("usr")[4], 
                        (par("usr")[4] - par("usr")[3])/50)
    lapply(seq_len(50), function(z){
      rect(legend_box_x[1],
           legend_box_y[z],
           legend_box_x[2],
           legend_box_y[z + 1],
           xpd = TRUE,
           col = legend_cols[z],
           border = NA)
    })
    
    # LEGEND LABELS
    
    # HEADER & RECORD ----------------------------------------------------------
    
    # ADD HEADER & RECORD PLOT(S)
    if(z %in% c(prod(layout), length(channels))){
      # HEADER
      if (!is.null(header)) {
        mtext(header,
              outer = TRUE,
              font = header_text_font,
              cex = header_text_size,
              col = header_text_col
        )
      }
      return(cyto_plot_record())
    }else{
      return(NULL)
    }
    
  })
  
  # REMOVE EMPTY PLOTS
  plots <- plots[!LAPPLY(plots, "is.null")]
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE
  if (getOption("cyto_plot_save")) {
    # PLOT METHOD
    if (is(x, basename(getOption("cyto_plot_method")))) {
      # CLOSE GRAPHICS DEVICE
      cyto_plot_complete()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN LIST OF RECORDED PLOTS
  invisible(plots)
  
}
