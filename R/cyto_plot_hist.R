## CYTO_PLOT_HIST --------------------------------------------------------------

#' Add histograms to empty plots
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}}, list of
#'   \code{\link[flowWorkspace:cytoframe]{cytoframes}} or a list of
#'   \code{\link[stats:density]{density}} objects.
#' @param parent name of the parent population to extract from GatingHierarchy
#'   or GatingSet objects.
#' @param channel name of the channels to be used to construct the plot.
#' @param overlay list of cytoframe objects to overlay.
#' @param events controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param hist_stat can be either \code{"count"}, \code{"percent"} or
#'   \code{"density"} to indicate the statistic to display on histograms, set to
#'   \code{"count"} by default. The \code{"percent"} option applies modal
#'   normalisation and expresses the result as a percentage.
#' @param hist_bins number of bins to use for histograms, set to 256 by default.
#' @param hist_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust the smoothness of the kernel
#'   density for histograms, set to \code{1} by default. Only values greater or
#'   equal to 1 are supported.
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
#' @param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'   sampling is applied with each \code{\link{cyto_plot_contour}} call, set to
#'   an arbitrary numeric by default. This behaviour can be turned off by
#'   setting this argument to NULL.
#' @param ... not in use.
#'
#' @importFrom methods is
#' @importFrom graphics abline polygon
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_plot_hist <- function(x,
                           parent = "root",
                           channel,
                           overlay = NA,
                           events = 1,
                           hist_stat = "count",
                           hist_bins = 256,
                           hist_smooth = 1,
                           hist_stack = 0.5,
                           hist_cols = NA,
                           hist_fill = NA,
                           hist_fill_alpha = 1,
                           hist_line_type = 1,
                           hist_line_width = 1,
                           hist_line_col = "black",
                           seed = 42,
                           ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # DENSITY LIST
  d <- NULL
  
  # X - CYTO_PLOT ARGUMENTS
  if(cyto_class(x, "cyto_plot")) {
    .args_update(x)
  # X - LIST OF DENSITY OBJECTS  
  } else if(cyto_class(x, "list") & all(LAPPLY(x, "cyto_class", "density"))) {
    d <- x
  # X - CYTOFRAME/CYTOSET/GATINGHIERARCHY/GATINGSET/LIST
  } else {
    # X -> LIST OF CYTOSETS
    if(cyto_class(x, c("flowFrame", "flowSet", "GatingSet"))) {
      x <- cyto_data_extract(x,
                             parent = parent,
                             format = "cytoset",
                             copy = FALSE)
      # OVERLAY
      if(!.all_na(overlay)) {
        if(!cyto_class(overlay, "list")) {
          overlay <- cyto_list(overlay)
        }
        x <- c(x, overlay)
      }
    }
    # SAMPLE
    if(events != 1){
      x <- cyto_sample(x,
                       events = events,
                       seed = seed)
    }
    # HISTOGRAMS
    args <- .args_list(...)
    args$channels <- cyto_channels_extract(x[[1]], channel)
    d <- do.call(".cyto_plot_hist", args)
  }
  
  # STACKING
  if (length(d) == 1) {
    hist_stack <- 0
  }
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)[-match("args", names(args))]
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # HIST_FILL ------------------------------------------------------------------
  
  # HIST_FILL COLOURS
  hist_fill <- .cyto_plot_hist_fill(x,
                                    hist_fill = hist_fill,
                                    hist_cols = hist_cols,
                                    hist_fill_alpha = hist_fill_alpha
  )
  
  # REPEAT ARGUMENTS -----------------------------------------------------------
  
  # ARGUMENTS TO REPEAT
  args <- .args_list()[c(
    "hist_line_type",
    "hist_line_width",
    "hist_line_col"
  )]
  
  # REPEAT ARGUMENTS
  args <- lapply(args, function(arg) {
    rep(arg, length.out = length(x))
  })
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # HORIZONTAL LINES -----------------------------------------------------------
  
  # YMIN PER LAYER
  ylim <- lapply(d, function(D){
    D$range
  })
  ymin <- lapply(ylim, `[[`, 1)
  
  # LINES UNDER DENSITY
  abline(
    h = ymin,
    col = hist_line_col,
    lwd = hist_line_width,
    lty = hist_line_type
  )
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # REVERSE PLOTTING ORDER - STACKS
  if (length(d) > 1 & hist_stack > 0) {
    d <- rev(d)
    hist_fill <- rev(hist_fill)
    hist_line_col <- rev(hist_line_col)
    hist_line_width <- rev(hist_line_width)
    hist_line_type <- rev(hist_line_type)
  }
  
  # HISTOGRAMS
  cnt <- 0
  mapply(
    function(D,
             hist_fill,
             hist_line_col,
             hist_line_width,
             hist_line_type) {
      # COUNTER
      cnt <<- cnt + 1
      # BYPASS NA
      if(!.all_na(D)) {
        # RANGES
        xrng <- c(D$x[1], D$x[length(D$x)]) # x values are sorted
        yrng <- D$range
        # PLOT - ANCHOR POLYGON
        if (!.all_na(d)) {
          polygon(
            x = c(
              xrng[1] - 0.0001 * (xrng[2] - xrng[1]),
              D$x,
              xrng[2] + 0.0001 * (xrng[2] - xrng[1])
            ),
            y = c(
              yrng[1] - 0.0001 * (yrng[2] - yrng[1]),
              D$y,
              yrng[1] - 0.0001 * (yrng[2] - yrng[1])
            ),
            col = hist_fill,
            border = hist_line_col,
            lwd = hist_line_width,
            lty = hist_line_type
          )
        }
      }
    },
    d,
    hist_fill,
    hist_line_col,
    hist_line_width,
    hist_line_type
  )
}
