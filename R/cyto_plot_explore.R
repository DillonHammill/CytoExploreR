## CYTO_PLOT_EXPLORE -----------------------------------------------------------

#' Explore cytometry data in a series of bivariate scatterplots
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied, set to the \code{"root"} node by
#'   default.
#' @param channels_x vector of channels or markers to plot on the x axis of the
#'   plots, by default we use channels returned by
#'   \code{\link{cyto_fluor_channels}}. If \code{channel} metadata variable is
#'   present each sample with be plotted with the specified channel on the x
#'   axis.
#' @param channels_y vector of channels or markers to plot on the y axis of the
#'   plots, by default we use all channels as returned by
#'   \code{\link{cyto_channels}}. \code{channels_y} overlapping with
#'   \code{channels_x} will be plotted as histograms.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param merge_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to NA by default to prevent merging. To merge
#'   all samples set this argument to \code{TRUE} or \code{"all"}.
#' @param overlay name(s) of the populations to overlay or a \code{cytoset},
#'   \code{list of cytosets} or \code{list of cytoset lists} containing
#'   populations to be overlaid onto the plot(s). This argument can be set to
#'   "children" or "descendants" when a \code{GatingSet} or
#'   \code{GatingHierarchy} to overlay all respective nodes.
#' @param order can be either \code{"channels"} or \code{"groups"} to control
#'   the order in which the data is plotted. Setting \code{order = "channels"}
#'   will plot each group on a single page with a plot for each channel in
#'   \code{channels_y}. On the other hand, setting \code{order = "groups"} will
#'   plot each channel on a separate page with a plot for each group in
#'   \code{merge_by}.
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied data. \code{cyto_plot} does not
#'   support in-line transformations and as such the transformations should be
#'   applied to the data prior to plotting. The transformerList is used
#'   internally to ensure that the axes on the constructed plots are
#'   appropriately labelled.
#' @param header can be set to NA to remove headers from plots. Alternatively,
#'   users can supply a vector of headers for each page, this can be tricky to
#'   get right so it is recommended for users to use the default headers where
#'   possible.
#' @param title text to display above each plot, set to the name of the group or
#'   channels displayed in the plot by default.
#' @param page_fill colour to use to fill the page prior to adding plot panels,
#'   set to \code{"white"} by default.
#' @param page_fill_alpha numeric [0,1] to control the transparency of the page
#'   colour, set to 1 by default to remove transparency.
#' @param ylim lower and upper limits of y axes (e.g. c(0,250000)).
#' @param ... additional arguments passed to \code{\link{cyto_plot}}.
#'
#' @return list of recorded plots.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot}}
#' @seealso \code{\link{cyto_plot_profile}}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- cyto_load(
#'   system.file(
#'     "extdata/Activation-GatingSet",
#'     package = "CytoExploreRData"
#'   )
#' )
#'
#' # Plot a sample in all relevant channels
#' cyto_plot_explore(
#'   gs,
#'   parent = "Live Cells",
#'   select = 32,
#'   order = "channels",
#'   channels_x = c("Va2", "CD4", "CD8", "CD44", "CD69", "CD11c"),
#'   channels_y = c("Va2", "CD4", "CD8", "CD44", "CD69", "CD11c"),
#' )
#'
#' # Compare differences between samples in all relevant channels
#' cyto_plot_explore(
#'   gs,
#'   parent = "Live Cells",
#'   select = 32,
#'   order = "groups",
#'   channels_x = c("Va2", "CD4", "CD8", "CD44", "CD69", "CD11c"),
#'   channels_y = c("Va2", "CD4", "CD8", "CD44", "CD69", "CD11c"),
#' )
#' }
#'
#' @export
cyto_plot_explore <- function(x,
                              parent = "root",
                              channels_x = NULL,
                              channels_y = NULL,
                              select = NULL,
                              merge_by = "name",
                              overlay = NA,
                              order = "channels",
                              layout,
                              axes_trans = NA,
                              header,
                              title,
                              page_fill = "white",
                              page_fill_alpha = 1,
                              ylim = c(NA, NA), 
                              ...) {
  
  # TODO: ADD POINT_STACK ARGUMENT TO OVERLAY MERGED DATA ONTO SINGLE PANEL IF NO OVERLAY
  
  # TODO: REMOVE EXCESS ARGUMENTS - PASS THROUGH ... TO .CYTO_PLOT_DATA()
  
  # TODO: SUPPORT LEGEND PANEL FOR OVERLAYS
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "explore")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  } else {
    on.exit(
      cyto_option("cyto_plot_data", FALSE)
    )
  }
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # METADATA
  pd <- cyto_details(x)
  
  # CONTROLS PROVIDED
  controls <- FALSE
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # X CHANNELS - DEFAULT TO FLUORESCENT CHANNELS
  if(is.null(channels_x)) {
    # PEAK DETECTOR X CHANNEL FOR CONTROLS
    if("channel" %in% colnames(pd)) {
      # DROP UNSTAINED CONTROLS
      rm_idx <- which(
        pd$channel %in% c("", "NA", NA, "Unstained", "unstained")
      )
      if(length(rm_idx) > 0) {
        x <- x[-rm_idx]
      }
      # UPDATE METADATA TO EXCLUDE UNSTAINED
      pd <- cyto_details(x)
      # CHECK CHANNELS
      channels_x <- cyto_channels_extract(
        x, 
        channels = pd$channel
      )
      # CONTROLS PROVIDED
      controls <- TRUE
    } else {
      channels_x <- cyto_fluor_channels(x)
    }
  } else {
    channels_x <- cyto_channels_extract(x, channels_x)
  }
  
  # Y CHANNELS - DEFAULT ALL CHANNELS
  if(is.null(channels_y)) {
    channels_y <- cyto_channels(x)
  } else {
    channels_y <- cyto_channels_extract(x, channels_y)
  }
  
  # CALL .CYTO_PLOT_DATA() - PASS TWO CHANNELS FOR FORMATTING
  args <- .args_list(...)
  args$channels <- c(
    channels_x[1], 
    channels_y[channels_y %in% channels_x[1]][1]
  )
  
  # PREPARE DATA
  x <- cyto_func_execute(
    ".cyto_plot_data",
    args
  )
  rm(args)
  
  # SIGNAL CYTO_PLOT_DATA CALL
  cyto_option("cyto_plot_data", TRUE)
  
  # PREPARE PLOT PARAMETERS ----------------------------------------------------
  
  # TODO: HANDLE CUSTOM LAYOUTS
  
  # PREPARE LAYOUT
  if(missing(layout)) {
    # # SWITCH LAYOUTS FOR SINGLE PLOTS
    # if(length(channels_x) == 1) {
    #   order <- "channels"
    # } else if(length(x) == 1) {
    #   order <- "groups"
    # }
    # CONTROLS USE CHHANEL ORDER
    if(controls) {
      order <- "channels"
    }
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      # EXCLUDE DUPLICATE CHANNELS
      layout <- .cyto_plot_layout(
        channels_y
      )
      # GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(x)
    }
  }
  
  # TOTAL PLOTS PER PAGE - BASED ON LAYOUT
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES - X CHANNEL AGAINST Y CHANNELS * GROUPS
  if(grepl("^c", order, ignore.case = TRUE)) {
    # NUMBER OF GROUPS TO PLOT
    n <- length(x)
    # NUMBER OF X CHANNELS PER GROUP
    if(controls) {
      nx <- 1
    } else {
      nx <- length(channels_x)
    }
    # PAGES PER GROUP - X CHANNEL SEPARATE PAGE BY UNIQUE Y CHANNELS
    pg <- ceiling(length(channels_y)/np) * nx
    # TOTAL PAGES
    tpg <- n * pg
    # TOTAL PLOTS
    tp <- n * length(channels_y) * nx
  # PAGES - GROUPS FOR each X/Y CHANNEL COMBINATION
  } else {
    # NOTE: GROUPS ORDER NOT SUPPORTED FOR CONTROLS
    # NUMBER OF X CHANNELS PER GROUP
    nx <- length(channels_x)
    # NUMBER OF GROUPS TO PLOT
    n <- length(channels_x) * length(channels_y)
    # PAGES PER GROUP - LENGTH(X)
    pg <- ceiling(length(x)/np)
    # TOTAL PAGES
    tpg <- n * pg
    # TOTAL PLOTS
    tp <- n * length(x)
  }
  
  # PREPARE PAGE_FILL ARGUMENTS
  page_fill <- rep(page_fill, length.out = tpg)
  page_fill_alpha <- rep(page_fill_alpha, length.out = tpg)
  
  # PREPARE HEADERS
  if(missing(header)) {
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      header <- rep(
        names(x),
        each = pg
      )
    # GROUP ORDER
    } else {
      header <- rep(
        unname(
          ulapply(
            channels_x, 
            function(z){
              ulapply(
                channels_y,
                function(v) {
                  if(v == z) {
                    cyto_markers_extract(
                      x[[1]],
                      channels = z,
                      append = TRUE
                    )
                  } else {
                    paste0(
                      cyto_markers_extract(
                        x[[1]],
                        channels = z,
                        append = TRUE
                      ),
                      " / ",
                      cyto_markers_extract(
                        x[[1]],
                        channels = v,
                        append = TRUE
                      )
                    )
                  }
                }
              )
            }
          )
        ),
        each = pg
      )
    }
    # SUPPLIED HEADERS
  } else {
    # NO HEADERS
    if(.all_na(header)) {
      header <- rep(NA, length.out = tpg)
    } else {
      # HEADER PROVIDED PER PAGE
      if(length(header) == tpg) {
        # HEADER PREPARED ALREADY
      # HEADER PROVIDED PER GROUP
      } else if(length(header) == n) {
        header <- rep(header, each = pg)
      # HEADERS SHOULD BE PROVIDED PER GROUP or ALL PAGES
      } else {
        stop(
          paste0(
            "Supply a header for each ",
            ifelse(grepl("^c", order, ignore.case = TRUE),
                   "group",
                   "channel"),
            " or page!"
          )
        )
      }
    }
  }
  
  # PREPARE TITLES
  if(missing(title)) {
   # CHANNEL ORDER - DEFAULT TITLES
    if(grepl("^c", order, ignore.case = TRUE)) {
      title <- rep(
        "",
        length.out = length(x) * nx * length(channels_y)
      )
    # GROUP ORDER
    } else {
      title <- rep(
        names(x),
        times = nx * length(channels_y)
      )
    }
  # TITLES SUPPLIED
  } else {
    title <- rep(
      title,
      length.out = length(x) * nx * length(channels_y)
    )
  }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # PROGRESS BAR - INCREMENTED BY CYTO_PLOT()
  pb <- cyto_progress(
    label = "cyto_plot_explore()",
    total = tp
  )
  
  # CALL CYTO_PLOT - CHANNEL ORDER
  if(grepl("^c", order, ignore.case = TRUE)) {
    # CONSTRUCT & RECORD PLOTS
    plots <- structure(
      lapply(
        seq_along(x), 
        function(z){
          # HEADER COUNTER
          cnt <- (z - 1) * pg
          # RECORDED PLOTS PER GROUP - LOOP THROUGH CHANNELS_X
          structure(
            lapply(
              if(controls) {
                1
              } else {
                seq_along(channels_x)
              }, 
              function(w) {
                # X CHANNEL
                if(controls) {
                  x_chan <- pd$channel[z]
                }else {
                  x_chan <- channels_x[w]
                }
                # RESTRICT Y CHANNELS
                y_chans <- channels_y
                # LOOP THROUGH CHANNELS_Y
                p <- lapply(
                  seq_along(y_chans), 
                  function(v) {
                    # Y CHANNEL
                    y_chan <- y_chans[v]
                    if(y_chan == x_chan) {
                      y_chan <- NULL
                      ylim <- c(NA, NA) # FORCE DEFAULT LIMITS FOR HISTOGRAM Y AXIS
                    }
                    # CONSTRUCT PLOT
                    cyto_plot(
                      x[[z]], # USE LIST METHOD CYTO_PLOT_DATA CALLED 
                      parent = parent,
                      channels = c(x_chan, y_chan),
                      axes_trans = axes_trans,
                      hist_layers = length(x[[z]]), # SINGLE PANEL ONLY
                      layout = layout,
                      title = title[
                        (z - 1) * nx * length(channels_y) + 
                          (w - 1) * length(channels_y) + v
                      ],
                      header = header[cnt + ceiling(w/np)],
                      page = if(v == length(y_chans)) {
                        TRUE
                      } else {
                        FALSE
                      },
                      page_fill = page_fill[cnt + ceiling(w/np)],
                      page_fill_alpha = page_fill_alpha[cnt + ceiling(w/np)],
                      ylim = ylim,
                      ...
                    )
                  }
                )
                # PREPARE PLOTS
                p[ulapply(p, "is.null")] <- NULL
                p <- lapply(p, `[[`, 1)
                return(p)
              }
            ),
            names = if(controls) {
              "controls"
            } else {
              channels_x
            }
          )
        }
      ),
      names = names(x)
    )
  # CALL CYTO_PLOT - GROUP ORDER - NOT SUPPORTED FOR CONTROLS
  } else {
    # CONSTRUCT & RECORD PLOTS
    plots <- structure(
      lapply(
        seq_along(channels_x), 
        function(z){
          # X CHANNEL
          x_chan <- channels_x[z]
          # HEADER COUNTER
          cnt <- (z - 1) * pg * length(channels_y)
          # LOOP THROUGH CHANNELS_Y
          structure(
            lapply(
              seq_along(channels_y),
              function(w) {
                # INCREMENT HEADER COUNTER
                cnt <- cnt + (w - 1) * pg
                # Y CHANNEL
                y_chan <- channels_y[w]
                if(y_chan == x_chan) {
                  y_chan <- NULL
                  ylim <- c(NA, NA) # FORCE DEFAULT Y AXIS LIMITS FOR HISTOGRAMS
                }
                # LOOP THROUGH GROUPS
                p <- structure(
                  lapply(
                    seq_along(x),
                    function(v) {
                      # CONSTRUCT PLOT
                      cyto_plot(
                        x[[v]], # USE LIST METHOD CYTO_PLOT_DATA CALLED
                        channels = c(x_chan, y_chan),
                        axes_trans = axes_trans,
                        hist_layers = length(x[[v]]), # SINGLE PANEL ONLY
                        layout = layout,
                        title = title[
                          (z - 1) * length(channels_y) * length(x) + 
                            (w - 1) * length(x) + v
                        ],
                        header = header[cnt + ceiling(v/np)],
                        page = if(v == length(x)) {
                          TRUE
                        } else {
                          FALSE
                        },
                        page_fill = page_fill[cnt + ceiling(v/np)],
                        page_fill_alpha = page_fill_alpha[cnt + ceiling(v/np)],
                        ylim = ylim,
                        ...
                      )
                    }
                  ),
                  names = names(x)
                )
                # PREPARE PLOTS
                p[ulapply(p, "is.null")] <- NULL
                p <- lapply(p, `[[`, 1)
                return(p)
              }
            ),
            names = channels_y
          )
        }
      ),
      names = channels_x
    )
  }
  
  # FORMAT RECORDED PLOTS
  plots[ulapply(plots, "is.null")] <- NULL
  if(length(plots) == 0) {
    plots <- NULL
  }
  
  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}