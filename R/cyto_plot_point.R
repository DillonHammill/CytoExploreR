## CYTO_PLOT_POINT -------------------------------------------------------------

#' Add points and contour lines to empty cyto_plot
#'
#' @param x either an object of class
#'   \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} or a list containing
#'   cytoframe or cytoset objects.
#' @param parent name of the parent population to extract from GatingHierarchy
#'   or GatingSet objects.
#' @param channels names of the channels used to construct the plot.
#' @param overlay optional argument if x is a cytoframe to overlay a list of
#'   cytoframe objects.
#' @param display controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param point_shape shape(s) to use for points in 2-D scatterplots, set to
#'   \code{"."} by default to maximise plotting speed.  See
#'   \code{\link[graphics:par]{pch}} for alternatives.
#' @param point_size numeric to control the size of points in 2-D scatter plots
#'   set to 2 by default.
#' @param point_col_scale vector of ordered colours to use for the density
#'   colour gradient of points.
#' @param point_cols vector colours to draw from when selecting colours for
#'   points if none are supplied to point_col.
#' @param point_col colour(s) to use for points in 2-D scatter plots, set to NA
#'   by default to use a blue-red density colour scale.
#' @param point_col_alpha numeric [0,1] to control point colour transparency in
#'   2-D scatter plots, set to 1 by default to use solid colours.
#' @param point_fast logical indicating whether points should be plotted using
#'   the \code{scattermore} package, set to FALSE by default. This is an
#'   optional feature so if you intend to use it, make sure that you have the
#'   scattermore package installed \code{install.packages("scattermore")}.
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 0 by default to exclude contour lines.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 2 by default.
#' @param contour_line_col colour to use for contour lines, set to "black" by
#'   default.
#' @param contour_line_alpha numeric [0,1] to control transparency of contour
#'   lines, set to 1 by default to remove transparency.
#' @param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'   sampling is applied with each \code{\link{cyto_plot_contour}} call, set to
#'   an arbitrary numeric by default. This behaviour can be turned off by
#'   setting this argument to NULL.
#' @param ... not in use.
#'
#' @importFrom flowCore exprs polygonGate
#' @importFrom graphics par rasterImage points
#' @importFrom grDevices col2rgb dev.size
#' @importFrom stats rnorm
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_plot_point <- function(x,
                            parent = "root",
                            channels,
                            overlay = NA,
                            display = 1,
                            point_shape = ".",
                            point_size = 2,
                            point_col_scale = NA,
                            point_cols = NA,
                            point_col = NA,
                            point_col_alpha = 1,
                            point_fast = FALSE,
                            contour_lines = 0,
                            contour_line_type = 1,
                            contour_line_width = 1,
                            contour_line_col = "black",
                            contour_line_alpha = 1, 
                            seed = 42, 
                            ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # X - CYTOFRAME/CYTOSET/GATINGHIERARCHY/GATINGSET
  if(cyto_class(x, c("flowFrame", "flowSet", "GatingSet"))) {
    x <- cyto_data_extract(x, 
                           parent = parent,
                           format = "cytoset",
                           copy = FALSE)
    # OVERLAY
    if(!.all_na(overlay)){
      if(!cyto_class(overlay, "list")) {
        overlay <- cyto_list(overlay)
      }
      x <- c(x, overlay)
    }
  # CYTO_PLOT ARGUMENTS
  } else if(cyto_class(x, "cyto_plot")) {
    .args_update(x)
  # CHECK LISTS
  } else if(!all(LAPPLY(x, cyto_class, c("flowFrame", "flowSet")))) {
    stop("'x' must be a list of cytoframes or cytosets!")
  }
  
  # CHANNELS
  channels <- cyto_channels_extract(x, channels)
  
  # SAMPLE DATA ----------------------------------------------------------------
  
  # DISPLAY
  if (display != 1) {
    x <- cyto_sample(x,
                     display = display,
                     seed = seed
    )
  }
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE
  .args_update(args)
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # COLOUR BY 3RD PARAMETER
  x <- lapply(seq_along(x), function(z) {
    # COLOUR BY 3RD PARAMETER
    if (point_col[z] %in% c(cyto_channels(x[[z]]), cyto_markers(x[[z]]))) {
      chan <- cyto_channels_extract(x[[z]], point_col[z])
      # SORT EVENTS IN CYTOSET
      return(
        cytoset(
          structure(
            list(x[[z]][[1]][order(cyto_data_extract(x[[z]],
                                                channels = chan,
                                                format = "matrix",
                                                copy = FALSE)[[1]][[1]]), ]),
           names = cyto_names(x[[z]])
          )
        )
      )
    }
    return(x[[z]])
  })
  
  # POINT_COL ------------------------------------------------------------------
  point_col <- .cyto_plot_point_col(x,
                                    channels = channels,
                                    point_col_scale = point_col_scale,
                                    point_cols = point_cols,
                                    point_col = point_col,
                                    point_col_alpha = point_col_alpha
  )
  
  # REPEAT ARGUMENTS -----------------------------------------------------------
  
  # ARGUMENTS TO REPEAT
  args <- .args_list()[c(
    "point_shape",
    "point_size",
    "contour_lines",
    "contour_line_type",
    "contour_line_width",
    "contour_line_col",
    "contour_line_alpha"
  )]
  
  # REPEAT ARGUMENTS
  args <- lapply(args, function(arg) {
    rep(arg, length.out = length(x))
  })
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  usr <- par("usr")
  
  # FAST PLOTTING --------------------------------------------------------------
  
  # GLOBAL OPTION
  if (point_fast) {
    if (requireNamespace("scattermore")) {
      # DEVICE SIZE
      size <- as.integer(dev.size("px") / dev.size("in") * par("pin"))
      point_fast <- TRUE
    } else {
      message("Fast plotting requires the 'scattermore' package.")
      message("Please install it using install.packages('scattermore').")
      message("Resorting to conventional plotting...")
      point_fast <- FALSE
    }
  } else {
    point_fast <- FALSE
  }
  
  # POINT & CONTOUR LINES LAYERS -----------------------------------------------
  
  # LAYERS
  lapply(seq_len(length(x)), function(z) {
    # EXTRACT MATRIX - CANNOT RESTRICT MATRIX
    cf_exprs <- cyto_data_extract(
      x[[z]],
      format = "matrix",
      copy = FALSE
    )[[1]][[1]]
    # POINTS - SKIP NO EVENTS
    if (!is.null(nrow(cf_exprs))) {
      # POINTS
      if (nrow(cf_exprs) != 0) {
        # SAMPLE-ID?
        ind <- which(
          LAPPLY(channels, function(v){
            grepl("^Sample.*ID.*$", v) # flowJo compatibility SampleID
          })
        )
        # SAMPLE ID - JITTER BARCODES
        if(length(ind) > 0) {
          cf_exprs[, ind] <- LAPPLY(unique(cf_exprs[, ind]), function(z) {
            rnorm(
              n = length(cf_exprs[, ind][cf_exprs[, ind] == z]),
              mean = z,
              sd = 0.1
            )
          })
        }
        # PLOT DEFAULT POINTS
        if (!point_fast) {
          # CONVENTIONAL PLOTTING
          points(
            x = cf_exprs[, channels[1]],
            y = cf_exprs[, channels[2]],
            pch = point_shape[z],
            cex = point_size[z],
            col = point_col[[z]]
          )
          # SCATTERMORE POINTS - LACK PCH CONTROL
        } else {
          # RASTER
          rasterImage(
            scattermore::scattermore(
              xy = cbind(cf_exprs[, channels[1]], cf_exprs[, channels[2]]),
              size = size,
              xlim = usr[1:2],
              ylim = usr[3:4],
              cex = 0.5 * point_size[[z]],
              rgba = col2rgb(point_col[[z]], alpha = TRUE),
              output.raster = TRUE
            ),
            xleft = usr[1],
            xright = usr[2],
            ybottom = usr[3],
            ytop = usr[4]
          )
        }
      }
    }
    # CONTOUR_LINES
    if (contour_lines[z] != 0) {
      cyto_plot_contour(x[[z]],
                        channels = channels,
                        display = 1,
                        contour_lines = contour_lines[z],
                        contour_line_type = contour_line_type[z],
                        contour_line_width = contour_line_width[z],
                        contour_line_col = contour_line_col[z],
                        contour_line_alpha = contour_line_alpha[z]
      )
    }
  })
  
  # INVISIBLE NULL RETURN
  invisible(NULL)
}
