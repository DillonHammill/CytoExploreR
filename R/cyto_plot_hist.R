## CYTO_PLOT_HIST --------------------------------------------------------------

#' Add histograms to empty plots
#'
#' @param x object of class flowFrame or a list of density objects.
#' @param channel name of the channels to be used to construct the plot.
#' @param overlay list of flowFRame objects to overlay.
#' @param display controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param hist_stat can be either \code{"count"}, \code{"percent"} or
#'   \code{"density"} to indicate the statistic to display on histograms, set to
#'   \code{"percent"} by default. The \code{"percent"} option applies modal
#'   normalisation and expresses the result as a percentage.
#' @param hist_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust the smoothness kernel
#'   density, set to \code{0.6} by default.
#' @param hist_stack numeric [0,1] indicating the degree of stacking for
#'   histograms, set to \code{0.5} by default.
#' @param hist_cols vector colours to draw from when selecting histogram fill
#'   colours if none are supplied to \code{hist_fill}.
#' @param hist_fill fill colour(s) for histograms, select from \code{hist_cols}
#'   if not supplied.
#' @param hist_fill_alpha numeric [0,1] used to control histogram fill colour
#'   transparency, set to \code{1} by default for solid colours.
#' @param hist_line_type line type(s) to use for histogram borders, set to 1 by
#'   default to use solid lines. See \code{\link[graphics:par]{lty}} for
#'   alternatives.
#' @param hist_line_width numeric to control line width(s) for histogram borders
#'   lines, set to 1 by default.
#' @param hist_line_col colour(s) for histogram borders, set to \code{"black"}
#'   by default.
#' @param ... not in use.
#'
#' @importFrom methods formalArgs
#' @importFrom graphics abline polygon
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_plot_hist
#'
#' @export
cyto_plot_hist <- function(x, ...) {
  UseMethod("cyto_plot_hist")
}

#' @rdname cyto_plot_hist
#' @export
cyto_plot_hist.flowFrame <- function(x,
                                     channel,
                                     overlay = NA,
                                     display = 1,
                                     hist_stat = "percent",
                                     hist_smooth = 0.6,
                                     hist_stack = 0.5,
                                     hist_cols = NA,
                                     hist_fill = NA,
                                     hist_fill_alpha = 1,
                                     hist_line_type = 1,
                                     hist_line_width = 1,
                                     hist_line_col = "black",
                                     ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
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
                           smooth = hist_smooth,
                           stat = hist_stat,
                           stack = hist_stack)
  
  # DENSITY_FILL ---------------------------------------------------------------
  
  # DENSITY_FILL COLOURS
  hist_fill <- .cyto_plot_density_fill(fr_dens,
                                       hist_fill = hist_fill,
                                       hist_cols = hist_cols,
                                       hist_fill_alpha = hist_fill_alpha)
  
  # REPEAT ARGUMENTS ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()[c("hist_line_type",
                         "hist_line_width",
                         "hist_line_col")]
  
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
    col = hist_line_col,
    lwd = hist_line_width,
    lty = hist_line_type
  )
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # Add density distributions - reverse plot order and colours
  if (!.all_na(overlay) & 
      hist_stack == 0) {
    
    mapply(
      function(fr_dens,
               hist_fill,
               hist_line_col,
               hist_line_width,
               hist_line_type) {
        if(!.all_na(fr_dens)){
          polygon(fr_dens,
                  col = hist_fill,
                  border = hist_line_col,
                  lwd = hist_line_width,
                  lty = hist_line_type)
        }
      }, fr_dens,
      hist_fill,
      hist_line_col,
      hist_line_width,
      hist_line_type)
    
  } else {
    
    mapply(
      function(fr_dens,
               hist_fill,
               hist_line_col,
               hist_line_width,
               hist_line_type) {
        if(!.all_na(fr_dens)){
          polygon(fr_dens,
                  col = hist_fill,
                  border = hist_line_col,
                  lwd = hist_line_width,
                  lty = hist_line_type)
        }
        
      }, rev(fr_dens),
      rev(hist_fill),
      rev(hist_line_col),
      rev(hist_line_width),
      rev(hist_line_type))
    
  }
}


#' @param x list of density objects for plotting
#' @rdname cyto_plot_hist
#' @export
cyto_plot_hist.list <- function(x,
                                hist_cols = NA,
                                hist_fill = NA,
                                hist_fill_alpha = 1,
                                hist_line_type = 1,
                                hist_line_width = 1,
                                hist_line_col = "black", 
                                ...){
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # DENSITY_FILL ---------------------------------------------------------------
  
  # DENSITY_FILL COLOURS
  hist_fill <- .cyto_plot_hist_fill(x,
                                    hist_fill = hist_fill,
                                    hist_cols = hist_cols,
                                    hist_fill_alpha = hist_fill_alpha)
  
  # REPEAT ARGUMENTS -----------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()[c("hist_line_type",
                         "hist_line_width",
                         "hist_line_col")]
  
  # REPEAT ARGUMENTS
  args <- lapply(args, function(arg){
    rep(arg, length.out = length(x))
  })
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # HORIZONTAL LINES -----------------------------------------------------------
  
  # YMIN PER LAYER
  ylim <- strsplit(names(x), "-")
  ymin <- as.numeric(lapply(ylim, `[[`, 1))
  
  # LINES UNDER DENSITY
  abline(
    h = ymin,
    col = hist_line_col,
    lwd = hist_line_width,
    lty = hist_line_type,
    xpd = FALSE
  )
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # Add density distributions - reverse plot order and colours
  if (length(x) > 1 & 
      all(floor(ymin) == 0)) {
    
    mapply(
      function(x,
               hist_fill,
               hist_line_col,
               hist_line_width,
               hist_line_type) {
        if(!.all_na(x)){
          polygon(x,
                  col = hist_fill,
                  border = hist_line_col,
                  lwd = hist_line_width,
                  lty = hist_line_type)
        }
      }, x,
      hist_fill,
      hist_line_col,
      hist_line_width,
      hist_line_type)
    
  } else {
    
    mapply(
      function(x,
               hist_fill,
               hist_line_col,
               hist_line_width,
               hist_line_type) {
        if(!.all_na(x)){
          polygon(x,
                  col = hist_fill,
                  border = hist_line_col,
                  lwd = hist_line_width,
                  lty = hist_line_type)
        }
      }, rev(x),
      rev(hist_fill),
      rev(hist_line_col),
      rev(hist_line_width),
      rev(hist_line_type))
    
  }
  
}