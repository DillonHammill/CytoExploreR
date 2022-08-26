## CYTO_PLOT_POINT -------------------------------------------------------------

#' Add points and contour lines to empty cyto_plot
#'
#' @param x either an object of class
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to extract from GatingHierarchy
#'   or GatingSet objects.
#' @param channels names of the channels used to construct the plot.
#' @param overlay name(s) of the populations to overlay or a \code{cytoset} or
#'   \code{list of cytosets} containing populations to be overlaid onto the
#'   plot(s). This argument can be set to "children" or "descendants" when a
#'   \code{GatingSet} or \code{GatingHierarchy} to overlay all respective nodes.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param events controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param point_shape shape(s) to use for points in 2-D scatterplots, set to
#'   \code{"."} by default to maximise plotting speed.  See
#'   \code{\link[graphics:par]{pch}} for alternatives.
#' @param point_size numeric to control the size of points in 2-D scatter plots
#'   set to 2 by default.
#' @param point_col_scale vector of ordered colours to use for the density
#'   colour gradient of points.
#' @param point_col_smooth logical indicating whether the 2D binned counts
#'   should be smoothed using kernel density estimates prior to selecting
#'   colours from \code{point_col_scale}, set to TRUE by default. Setting
#'   \code{point_col_smooth} to FALSE will significantly improve plotting speed
#'   on less powerful machines but produce more granular plots.
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
#' @importFrom flowCore exprs
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
                            select = NULL,
                            events = 1,
                            point_shape = ".",
                            point_size = 2,
                            point_col_scale = NA,
                            point_col_smooth = TRUE,
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
  if(cyto_class(x, c("flowSet", "GatingSet"))) {
    # PREPARE DATA - .CYTO_PLOT_DATA()
    # CYTO_PLOT NOT CALLED AFTER POINTS - DON'T SET CYTO_PLOT_DATA OPTION
    x <- .cyto_plot_data(
      x,
      parent = parent,
      overlay = overlay,
      merge_by = "all",
      select = select,
      events = events,
      barcode = FALSE,
      seed = seed
    )[[1]]
    # PULL DOWN ARGUMENTS
    args <- .args_list(...)
  # CYTO_PLOT ARGUMENTS
  } else if(cyto_class(x, "cyto_plot")) {
    args <- x
  # CHECK LISTS
  } else if(!all(LAPPLY(x, cyto_class, "flowSet"))) {
    stop("'x' must be a list of cytosets!")
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # REPEAT ARGUMENTS PER LAYER
  layer_args <- c(
    "point_shape",
    "point_size",
    "point_col",
    "point_col_alpha",
    names(args)[grepl("contour_", names(args))]
  )
  lapply(
    layer_args,
    function(z) {
      args[[z]] <<- rep(args[[z]], length.out = length(args$x))
    }
  )
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # BASE LAYER - BKDE2D
  if(.all_na(args$point_col[1])) {
    # BKDE
    if(!"bkde2d" %in% names(args)) {
      args$bkde2d <- list(
        counts = NA,
        bins = NA,
        bkde = NA
      )
    }
    # RECOMPUTE BKDE
    if(.all_na(args$bkde2d$bkde)) {
      # RECOMPUTE BKDE2D - DEFAULT BINS
      args$bkde2d <- cyto_apply(
        args$x[[1]],
        "cyto_stat_bkde2d",
        input = "matrix",
        channels = args$channels,
        limits = list(.par("usr")[[1]][1:2],
                      .par("usr")[[1]][3:4]),
        smooth = args$point_col_smooth,
        copy = FALSE,
        simplify = FALSE
      )[[1]]
    }
  } else {
    args$bkde2d <- list(
      counts = NA,
      bins = NA,
      bkde = NA
    )
  }
  
  # SORT EVENTS BY COLOUR AND|OR SIZE
  sort_by <- c()
  if(args$point_col[1] %in% c(cyto_channels(args$x[[1]]),
                              cyto_markers(args$x[[1]]))) {
    sort_by <- c(sort_by, "col")
  }
  if(cyto_class(args$point_size, "list", TRUE)) {
    sort_by <- c(sort_by, "size")
  }
  
  # TODO: SORTING ONLY PERFORMED ON BASE LAYER
  
  # SORTING REQUIRED
  if(length(sort_by) > 0) {
    # TWO WAY SORT - SIZE & COLOUR
    if(length(sort_by) == 2) {
      ind <- order(
        -args$point_size[[1]],
        cyto_exprs(
          args$x[[1]],
          channels = args$point_col[1],
          drop = TRUE
        )[[1]]
      )
    # SINGLE WAY SORT
    } else {
      # SIZE
      if(sort_by %in% "size") {
        ind <- order(
          -args$point_size[[1]]
        )
      # COLOUR
      } else {
        ind <- order(
          cyto_exprs(
            args$x[[1]],
            channels = args$point_col[1],
            drop = TRUE
          )[[1]]
        )
      }
    }
    # SORT EVENTS
    cyto_exprs(args$x[[1]]) <- list(
      cyto_exprs(args$x[[1]])[[1]][
        ind,
      ]
    )
    # SORT POINT_SIZE - POINT_COL HANDLED LATER
    if(cyto_class(args$point_size, "list", TRUE)) {
      args$point_size[[1]] <- args$point_size[[1]][ind]
    }
    rm(ind)
  }
  
  # POINT_COL ------------------------------------------------------------------
  
  # GET POINT COLOURS - PASS KEY_SCALE & BKDE2D
  args$point_col <- cyto_func_call(
    ".cyto_plot_point_col",
    args
  )
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  usr <- .par("usr")[[1]]
  
  # FAST PLOTTING --------------------------------------------------------------
  
  # SCATTERMORE FAST PLOTIING
  if (args$point_fast) {
    # SCATTERMORE
    cyto_require(
      "scattermore",
      source = "CRAN"
    )
    # DEVICE SIZE
    dev_size <- as.integer(dev.size("px") / dev.size("in") * .par("pin")[[1]])
  }
  
  # POINT & CONTOUR LINES LAYERS -----------------------------------------------
  
  # ADD POINTS & CONTOURS
  lapply(
    seq_along(args$x),
    function(z) {
      # EXTRACT MATRIX
      exprs <- cyto_data_extract(
        args$x[[z]],
        format = "matrix",
        copy = FALSE
      )[[1]][[1]]
      # POINTS - SKIP NO EVENTS
      if(!is.null(nrow(exprs))) {
        # POINTS - BYPASS EMPTY CYTOFRAME
        if(nrow(exprs) != 0) {
          # SAMPLE-ID
          ind <- which(
            LAPPLY(
              args$channels,
              function(v) {
                grepl("^Sample-ID$", v) # NOT FLOWJO SAMPLE IDS
              }
            )
          )
          # JITTER BARCODES FOR SAMPLE-ID
          if(length(ind) > 0) {
            exprs[, ind] <- LAPPLY(
              unique(exprs[, ind]),
              function(w){
                rnorm(
                  n = length(
                    exprs[, ind][exprs[, ind] == w]
                  ),
                  mean = w,
                  sd = 0.1
                )
              }
            )
          }
          # PLOT DEFAULT POINTS
          if (!args$point_fast) {
            # CONVENTIONAL PLOTTING
            points(
              x = exprs[, args$channels[1]],
              y = exprs[, args$channels[2]],
              pch = args$point_shape[z],
              cex = if(cyto_class(args$point_size, "list", TRUE)) {
                args$point_size[[z]]
              } else {
                args$point_size[z]
              },
              col = if(args$point_shape[z] %in% c(21:25)) {
                "black"
              } else {
                args$point_col[[z]]
              },
              bg = if(args$point_shape[z] %in% c(21:25)) {
                args$point_col[[z]]
              } else {
                "black"
              }
            )
          # SCATTERMORE POINTS - LACK PCH CONTROL
          } else {
            # RASTER
            rasterImage(
              cyto_func_call(
                "scattermore::scattermore",
                list(
                  xy = cbind(
                    exprs[, args$channels[1]], 
                    exprs[, args$channels[2]]
                  ),
                  size = dev_size,
                  xlim = usr[1:2],
                  ylim = usr[3:4],
                  cex = if(cyto_class(args$point_size, "list", TRUE)){
                    0.5 * args$point_size[[z]][1]
                  }else {
                    0.5 * args$point_size[z]
                  },
                  rgba = col2rgb(
                    args$point_col[[z]], 
                    alpha = TRUE
                  ),
                  output.raster = TRUE
                )
              ),
              xleft = usr[1],
              xright = usr[2],
              ybottom = usr[3],
              ytop = usr[4]
            )
          }
          # CONTOUR_LINES
          if (args$contour_lines[z] != 0) {
            cyto_plot_contour(
              args$x[[z]],
              channels = args$channels,
              events = 1,
              contour_lines = args$contour_lines[z],
              contour_line_type = args$contour_line_type[z],
              contour_line_width = args$contour_line_width[z],
              contour_line_col = args$contour_line_col[z],
              contour_line_alpha = args$contour_line_alpha[z],
              bkde2d = args$bkde2d
            )
          }
        }
      }
    }
  )
  
  # INVISIBLE NULL RETURN
  invisible(NULL)
}
