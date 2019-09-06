# CYTO_PLOT_DENSITY ------------------------------------------------------------

#' Add density distributions to empty plots
#'
#' @param x object of class flowFrame or a list of density objects.
#' @param channel name of the channels to be used to construct the plot.
#' @param overlay list of flowFRame objects to overlay.
#' @param display controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param density_modal logical indicating whether density should be normalised
#'   to mode and presented as a percentage. Set to \code{TRUE} by default.
#' @param density_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust kernel density.
#' @param density_stack numeric [0,1] indicating the degree of offset for
#'   overlaid populations, set to 0.5 by default.
#' @param density_cols vector of colls to use to select density_fill colours
#' @param density_fill colour(s) used to fill polygons.
#' @param density_fill_alpha numeric [0,1] used to control fill transparency,
#'   set to 1 by default to remove transparency.
#' @param density_line_type line type(s) to use for border(s), set to solid
#'   lines by default.
#' @param density_line_width line width for border.
#' @param density_line_col colour(s) for border line, set to "black" by default.
#'
#' @importFrom methods formalArgs
#' @importFrom graphics abline polygon
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_plot_density
#'
#' @export
cyto_plot_density <- function(x, ...) {
  UseMethod("cyto_plot_density")
}

#' @rdname cyto_plot_density
#' @export
cyto_plot_density.flowFrame <- function(x,
                                        channel,
                                        overlay = NA,
                                        display = 1,
                                        density_modal = TRUE,
                                        density_smooth = 1.5,
                                        density_stack = 0.5,
                                        density_cols = NA,
                                        density_fill = NA,
                                        density_fill_alpha = 1,
                                        density_line_type = 1,
                                        density_line_width = 1,
                                        density_line_col = "black") {

  # PREPARE DATA ---------------------------------------------------------------
  
  # LIST OF FLOWFRAMES
  if (!.all_na(overlay)) {
    fr_list <- c(list(x), cyto_convert(overlay, "list of flowFrames"))
  } else {
    fr_list <- list(x)
  }

  # SAMPLE DATA ----------------------------------------------------------------
  
  # DISPLAY
  if(display != 1){
    fr_list <- cyto_sample(fr_list,
                           display = display,
                           seed = 56)
  }
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()

  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)

  # UPDATE ARGUMENTS
  .args_update(args)
  
  # DENSITY --------------------------------------------------------------------
  
  # KERNEL DENSITY
  fr_dens <- .cyto_density(fr_list,
                           channel = channel,
                           smooth = density_smooth,
                           modal = density_modal,
                           stack = density_stack)

  # DENSITY_FILL ---------------------------------------------------------------
  
  # DENSITY_FILL COLOURS
  density_fill <- .cyto_plot_density_fill(fr_dens,
                                          density_fill = density_fill,
                                          density_cols = density_cols,
                                          density_fill_alpha = density_fill_alpha)

  # REPEAT ARGUMENTS ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()[c("density_line_type",
                         "density_line_width",
                         "density_line_col")]
  
  # REPEAT ARGUMENTS
  args <- lapply(args, function(arg){
    rep(arg, length.out = length(fr_list))
  })
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # HORIZONTAL LINES -----------------------------------------------------------
  
  # YMIN PER LAYER
  ylim <- strsplit(names(fr_dens), "-")
  ymin <- as.numeric(lapply(ylim, `[[`, 1))
  
  # LINES UNDER DENSITY
  abline(
    h = ymin,
    col = density_line_col,
    lwd = density_line_width,
    lty = density_line_type
  )
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # Add density distributions - reverse plot order and colours
  if (!.all_na(overlay) & 
      density_stack == 0) {
    
    mapply(
      function(fr_dens,
                     density_fill,
                     density_line_col,
                     density_line_width,
                     density_line_type) {
        if(!.all_na(fr_dens)){
          polygon(fr_dens,
            col = density_fill,
            border = density_line_col,
            lwd = density_line_width,
            lty = density_line_type)
        }
      }, fr_dens,
      density_fill,
      density_line_col,
      density_line_width,
      density_line_type)
    
  } else {
    
    mapply(
      function(fr_dens,
                     density_fill,
                     density_line_col,
                     density_line_width,
                     density_line_type) {
        if(!.all_na(fr_dens)){
          polygon(fr_dens,
                  col = density_fill,
                  border = density_line_col,
                  lwd = density_line_width,
                  lty = density_line_type)
        }

      }, rev(fr_dens),
      rev(density_fill),
      rev(density_line_col),
      rev(density_line_width),
      rev(density_line_type))
    
  }
}


#' @param x list of density objects for plotting
#' @rdname cyto_plot_density
#' @export
cyto_plot_density.list <- function(x,
                                   density_cols = NA,
                                   density_fill = NA,
                                   density_fill_alpha = 1,
                                   density_line_type = 1,
                                   density_line_width = 1,
                                   density_line_col = "black", ...){
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # DENSITY_FILL ---------------------------------------------------------------
  
  # DENSITY_FILL COLOURS
  density_fill <- .cyto_plot_density_fill(x,
                                          density_fill = density_fill,
                                          density_cols = density_cols,
                                          density_fill_alpha = density_fill_alpha)
  
  # REPEAT ARGUMENTS -----------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()[c("density_line_type",
                         "density_line_width",
                         "density_line_col")]
  
  # REPEAT ARGUMENTS
  args <- lapply(args, function(arg){
    rep(arg, length.out = length(fr_list))
  })

  # UPDATE ARGUMENTS
  .args_update(args)
  
  # HORIZONTAL LINES -----------------------------------------------------------
  
  # YMIN PER LAYER
  ylim <- strsplit(names(fr_dens), "-")
  ymin <- as.numeric(lapply(ylim, `[[`, 1))
  
  # LINES UNDER DENSITY
  abline(
    h = ymin,
    col = density_line_col,
    lwd = density_line_width,
    lty = density_line_type,
    xpd = FALSE
  )
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # Add density distributions - reverse plot order and colours
  if (length(x) > 1 & 
      all(floor(ymin) == 0)) {
    
    mapply(
      function(x,
               density_fill,
               density_line_col,
               density_line_width,
               density_line_type) {
        if(!.all_na(x)){
          polygon(x,
                  col = density_fill,
                  border = density_line_col,
                  lwd = density_line_width,
                  lty = density_line_type)
        }
      }, x,
      density_fill,
      density_line_col,
      density_line_width,
      density_line_type)
    
  } else {
    
    mapply(
      function(x,
               density_fill,
               density_line_col,
               density_line_width,
               density_line_type) {
        if(!.all_na(x)){
          polygon(x,
                  col = density_fill,
                  border = density_line_col,
                  lwd = density_line_width,
                  lty = density_line_type)
        }
      }, rev(x),
      rev(density_fill),
      rev(density_line_col),
      rev(density_line_width),
      rev(density_line_type))
    
  }
  
}