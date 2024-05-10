## CYTO_PLOT_MAP ---------------------------------------------------------------

#' Plot expression profiles for dimension reduced maps
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied, set to the \code{"root"} node by
#'   default.
#' @param channels names of the channels or markers to use on the x and y axes
#'   of the plots, normally dimension reduced parameters such as
#'   \code{c("FIt-SNE-1", "FIt-SNE-2")}.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param merge_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to NA by default to prevent merging. To
#'   merge all samples set this argument to \code{TRUE} or \code{"all"}.
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
#' @param alias name of the gated population(s) to gated in the plot when a
#'   GatingHierarchy or GatingSet object is supplied. Setting alias to "" will
#'   automatically plot any gates constructed in the supplied channels. alias is
#'   equivalent to the gate argument for cytoset objects.
#' @param gate gate objects to be plotted, can be either objects of class
#'   rectangleGate, polygonGate, ellipsoidGate, quadGate or filters. Lists of
#'   these supported gate objects are also supported.
#' @param negate logical indicating whether a label should be included for the
#'   negated population when gate objects are supplied, set to FALSE by default.
#'   Setting negate = TRUE will result in the creation of a boolean filter that
#'   contains all the events outside the gates supplied to alias or gate. If
#'   such a boolean gate exists in the GatingSet/GatingHierarchy it will
#'   automatically be extracted when alias = "". In order to prevent plotting of
#'   these boolean gates, users will need to explicitly pass the names of the
#'   gates they want to display to alias.
#' @param point_col names of the channels or markers to use for the colour scale
#'   of each plot, set to all fluorescent channels by default.
#' @param point_col_scale vector of colours to use for the colour scale within
#'   each plot. \code{point_col_scale} can also be supplied as a list of length
#'   \code{point_col} to use a different colour scale for each marker/channel.
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param header can be set to NA to remove headers from plots. Alternatively,
#'   users can supply a vector of headers for each page, this can be tricky to
#'   get right so it is recommended for users to use the default headers where
#'   possible.
#' @param title text to display above each plot, set to the name of the group or
#'   channel displayed in the plot by default.
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied data. \code{cyto_plot} does not
#'   support in-line transformations and as such the transformations should be
#'   applied to the data prior to plotting. The transformerList is used
#'   internally to ensure that the axes on the constructed plots are
#'   appropriately labelled.
#' @param ... additional arguments passed to \code{\link{cyto_plot}}.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot_calibrate}}
#' @seealso \code{\link{cyto_map}}
#' @seealso \code{\link{cyto_plot}}
#' @seealso \code{\link{cyto_plot_profile}}
#' @seealso \code{\link{cyto_plot_explore}}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transformations
#' gs <- cyto_transform(gs)
#'
#' # Gates
#' cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#'
#' # Dimensionality reduction
#' gs <- cyto_map(
#'   gs[1:32],
#'   parent = "Live Cells",
#'   channels = c(
#'     "Va2",
#'     "CD4",
#'     "CD8",
#'     "CD44",
#'     "CD69",
#'     "CD11c",
#'     "FSC-A",
#'     "SSC-A"
#'   ),
#'   type = "UMAP"
#' )
#'
#' # Calibrate cyto_plot colour scales
#' cyto_plot_calibrate(
#'   gs,
#'   parent = "Live Cells"
#' )
#'
#' # Visualisation
#' cyto_plot_map(
#'   gs,
#'   parent = "Live Cells",
#'   merge_by = "all",  # global embedding
#'   channels = c("UMAP-1", "UMAP-2"),
#'   point_col = c(
#'     "Va2",
#'     "CD4",
#'     "CD8",
#'     "CD44",
#'     "CD69",
#'     "CD11c",
#'     "FSC-A",
#'     "SSC-A"
#'   )
#' )
#' }
#'
#' @export
cyto_plot_map <- function(x,
                          parent = "root",
                          channels = NULL,
                          select = NULL,
                          merge_by = "name",
                          overlay = NA,
                          order = "channels",
                          alias = NA,
                          gate = NA,
                          negate = FALSE,
                          point_col = NA,
                          point_col_scale = NA,
                          layout,
                          header,
                          title,
                          axes_trans = NA,
                          ...) {
  
  # TODO: ADD SUPPORT FOR ALIAS
  
  # CYTO_PLOT_COMPLETE ----------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "map")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  } else {
    on.exit(
      cyto_option("cyto_plot_data", FALSE)
    )
  }

  # CYTO_PLOT_CALIBRATE --------------------------------------------------------
  
  # CALIBRATION RECOMMENDED
  if(is.null(.cyto_plot_calibrate_recall())) {
    message(
      paste0(
        "Please call 'cyto_plot_calibrate()' before calling ",
        "'cyto_plot_map()' to properly calibrate colour scales."
      )
    )
  }
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # CHANNELS
  if(is.null(channels)) {
    stop(
      "Supply the names of the channels/markers to plot to 'channels'."
    )
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # POINT_COL
  if(.all_na(point_col)) {
    # DEFAULT - ALL FLUORESCENT  CHANNELS
    point_col <- cyto_fluor_channels(x)
    # POINT_COL SUPPLIED
  } else {
    point_col <- cyto_channels_extract(x, point_col)
  }
  
  # GATES
  gate <- .cyto_plot_gates(
    x,
    parent = parent,
    channels = channels,
    alias = alias, 
    gate = gate,
    merge_by = merge_by,
    select = select,
    negate = negate
  )
  
  # CALL .CYTO_PLOT_DATA() - PASS TWO CHANNELS FOR FORMATTING
  args <- .args_list(...)
  
  # PROGRESS BAR
  args$progress <- "cyto_plot_map()"
  
  # PREPARE DATA
  x <- cyto_func_execute(
    ".cyto_plot_data",
    args
  )
  rm(args)
  
  # SIGNAL CYTO_PLOT_DATA CALL
  cyto_option("cyto_plot_data", TRUE)
  
  # PREPARE PLOT PARAMETERS ----------------------------------------------------

  # POINT_COL_SCALE - LIST - DIFFERENT COLOUR SCALE PER CHANNEL
  if(cyto_class(point_col_scale, "list", TRUE)) {
    point_col_scale <- structure(
      rep(
        point_col_scale, 
        length.out = length(point_col)
      ), 
      names = point_col
    )
  # POINT_COL_SCALE - VECTOR/FUNCTION - SAME COLOUR SCALE PER CHANNEL
  } else {
    point_col_scale <- structure(
      lapply(
        seq_along(point_col), 
        function(z) {
          return(
            point_col_scale
          )
        }
      ),
      names = point_col
    )
  }
  
  # LAYOUT
  if(missing(layout)) {
    if(length(point_col) == 1) {
      order <- "groups"
    } else if(length(x) == 1) {
      order <- "channels"
    }
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      layout <- .cyto_plot_layout(point_col)
    # GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(x)
    }
  }
  
  # TOTAL PLOTS PER PAGE
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES - CHANNELS ON SAME PAGE
  if(grepl("^c", order, ignore.case = TRUE)) {
    # NUMBER OF GROUPS TO PLOT
    n <- length(x)
    # PAGES PER GROUP
    pg <- ceiling(length(point_col)/np)
    # TOTAL PAGES
    tpg <- n * pg
    # TOTAL PLOTS
    tp <- n * length(point_col)
  # PAGES - GROUPS ON SAME PAGE
  } else {
    # NUMBER OF GROUPS TO PLOT
    n <- length(point_col)
    # PAGES PER GROUP
    pg <- ceiling(length(x)/np)
    # TOTAL PAGES
    tpg <- n * pg
    # TOTAL PLOTS
    tp <- n * length(x)
  }
  
  # DEFAULT HEADERS
  if(missing(header)) {
    # CHANNEL ORDER - GROUPS AS HEADERS
    if(grepl("^c", order, ignore.case = TRUE)) {
      header <- rep(
        names(x),
        each = pg
      )
    # GROUP ORDER - CHANNELS AS HEADERS
    } else {
      header <- rep(
        cyto_markers_extract(
          x,
          channels = point_col,
          append = TRUE
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
      # REPEAT HEADERS - MULTIPLE PAGES
      if(length(header) == tpg) {
        header <- rep(header, each = pg)
      } else {
        stop(
          paste0(
            "Supply a header for each ",
            ifelse(grepl("^c", order, ignore.case = TRUE),
                   "group!",
                   "channel!")
          )
        )
      }
    }
  }
  
  # TITLES - DEFAULT
  if(missing(title)) {
    # CHANNEL ORDER - POINT_COL HEADERS
    if(grepl("^c", order, ignore.case = TRUE)) {
      title <- rep(
        cyto_markers_extract(
          x[[1]],
          channels = point_col,
          append = TRUE
        ), 
        times = length(x)
      )
    # GROUP ORDER
    } else {
      title <- rep(
        names(x), 
        times = length(point_col)
      )
    }
  # TITLES SUPPLIED
  } else {
    title <- rep(
      title, 
      length.out = length(x) * length(point_col)
    )
  }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # CREATE PROGRESS BAR - INCREMENTED BY CYTO_PLOT()
  pb <- cyto_progress(
    label = "cyto_plot_map()",
    total = tp
  )
  
  # CALL CYTO_PLOT - CHANNEL ORDER
  if(grepl("^c", order, ignore.case = TRUE)) {
    # CONSTRUCT PLOTS
    plots <- structure(
      lapply(
        seq_along(x), 
        function(z){
          # HEADER COUNTER
          cnt <- (z - 1) * pg
          # RECORD PLOTS PER GROUP
          p <- structure(
             lapply(
               seq_along(point_col),
               function(w){
                 # CONSTRUCT PLOT
                 cyto_plot(
                   x[[z]], # USE LIST METHOD - CYTO_PLOT DATA CALLED
                   channels = channels,
                   gate = gate[[z]],
                   axes_trans = axes_trans,
                   layout = layout,
                   header = header[cnt + ceiling(w/np)],
                   title = title[(z - 1) * length(point_col) + w],
                   page = if(w == length(point_col)) {
                     TRUE
                   } else {
                     FALSE
                   },
                   point_col = point_col[w],
                   point_col_scale = point_col_scale[[w]],
                   ...
                 )
               }
            ),
            names = point_col
          )
          # PREPARE RECORDED PLOTS
          p[LAPPLY(p, "is.null")] <- NULL
          p <- lapply(p, `[[`, 1)
          return(p)
        }
      ),
      names = names(x)
    )
  # CALL CYTO_PLOT - GROUP ORDER
  } else {
    # CONSTRUCT PLOTS
    plots <- structure(
      lapply(
        seq_along(point_col), 
        function(z){
          # HEADER COUNTER
          cnt <- (z - 1) * pg
          # RECORD PLOTS PER CHANNEL
          p <- structure(
            lapply(
              seq_along(x),
              function(w) {
                # CONSTRUCT PLOT
                cyto_plot(
                  x[[w]],
                  channels = channels,
                  gate = gate[[w]],
                  axes_trans = axes_trans,
                  layout = layout,
                  header = header[cnt + ceiling(w/np)],
                  title = title[(z - 1) * length(x) + w],
                  page = if(w == length(x)) {
                    TRUE
                  } else {
                    FALSE
                  },
                  point_col = point_col[z],
                  point_col_scale = point_col_scale[[z]],
                  ...
                )
              }
            ),
            names = names(x)
          )
          # PREPARE RECORDED PLOTS
          p[LAPPLY(p, "is.null")] <- NULL
          return(p)
        }
      ),
      names = point_col
    )
  }
  
  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}
