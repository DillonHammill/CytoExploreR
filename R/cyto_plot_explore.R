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
#'   \code{\link{cyto_fluor_channels}}.
#' @param channels_y vector of channels or markers to plot on the y axis of the
#'   plots, by default we use all channels as returned by
#'   \code{\link{cyto_channels}}.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param merge_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to "name" by default to prevent merging. To
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
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param header can be set to NA to remove headers from plots. Alternatively,
#'   users can supply a vector of headers for each page, this can be tricky to
#'   get right so it is recommended for users to use the default headers where
#'   possible.
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
                              header,
                              ...) {
  
  # TODO: HEADERS INCORRECT FOR LAYOUTS > DEFAULT
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "explore")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # PREPARE PLOT PARAMETERS ----------------------------------------------------
  
  # X CHANNELS - DEFAULT TO FLUORESCENT CHANNELS
  if(is.null(channels_x)) {
    channels_x <- cyto_fluor_channels(x)
  } else {
    channels_x <- cyto_channels_extract(x, channels_x)
  }
  
  # Y CHANNELS - DEFAULT ALL CHANNELS
  if(is.null(channels_y)) {
    channels_y <- cyto_channels(x)
  } else {
    channels_y <- cyto_channels_extract(x, channels_y)
  }
  
  # EXTRACT NAMES OF VARIABLES FOR MERGING
  if(cyto_class(merge_by, "list", TRUE)) {
    vars <- names(merge_by)
  } else {
    vars <- merge_by
  }
  
  # GROUPS
  grps <- cyto_groups(
    x,
    select = select,
    group_by = merge_by,
    details = TRUE
  )
  
  # LAYOUT 
  if(missing(layout)) {
    # SWITCH LAYOUTS FOR SINGLE PLOTS
    if(length(channels_x) == 1) {
      order <- "groups"
    } else if(length(grps) == 1) {
      order <- "channels"
    }
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      # LAYOUT ACCOUNTING FOR REMOVAL OF DUPLICATE CHANNELS
      layout <- .cyto_plot_layout(
        seq_len(
          max(
            LAPPLY(channels_x, function(z){
              return(
                length(channels_y[!channels_y %in% z])
              )
            })
          )
        )
      )
    # GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(grps)
    }
  }
  
  # TODO: GROUP ORDER ACCOUNT FOR DUPLICATE CHANNELS
  
  # TOTAL PLOTS PER PAGE - BASED ON LAYOUT
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES
  if(grepl("^c", order, ignore.case = TRUE)) {
    n <- length(grps)
    pg <- ceiling(length(channels_y)/np)
    tpg <- n * pg
  } else {
    n <- length(channels_y)
    pg <- ceiling(length(grps)/np)
    tpg <- n * pg
  }
  
  # PREPARE HEADERS
  if(missing(header)) {
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      header <- rep(
        names(grps),
        each = pg
      )
    # GROUP ORDER
    } else {
      header <- rep(
        cyto_markers_extract(
          x,
          channels = channels_x,
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
  
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # CALL CYTO_PLOT - CHANNEL ORDER
  if(grepl("^c", order, ignore.case = TRUE)) {
    # CONSTRUCT & RECORD PLOTS
    plots <- lapply(seq_along(grps), function(z){
      # HEADER COUNTER
      cnt <- (z - 1) * pg
      # RECORDED PLOTS PER GROUP - LOOP THROUGH CHANNELS_X
      lapply(seq_along(channels_x), function(w){
        # X CHANNEL
        x_chan <- channels_x[w]
        # RESTRICT Y CHANNELS
        y_chans <- channels_y[!channels_y %in% x_chan]
        # LOOP THROUGH CHANNELS_Y
        lapply(seq_along(y_chans), function(v){
          # Y CHANNEL
          y_chan <- y_chans[v]
          # CONSTRUCT PLOT
          p <- cyto_plot(
            x,
            parent = parent,
            select = grps[[z]][, "name"],
            merge_by = vars,
            channels = c(x_chan, y_chan),
            overlay = overlay,
            layout = layout,
            header = header[cnt + ceiling(w/np)],
            page = if(v == length(y_chans)) {
              TRUE
            } else {
              FALSE
            },
            ...
          )
          return(p)
        })
      })
    })
  # CALL CYTO_PLOT - GROUP ORDER
  } else {
    # CONSTRUCT & RECORD PLOTS
    plots <- lapply(seq_along(channels_x), function(z){
      # HEADER COUNTER
      cnt <- (z - 1) * pg
      # X CHANNEL
      x_chan <- channels_x[z]
      # LOOP THROUGH CHANNELS_Y
      plots <- lapply(seq_along(channels_y), function(w){
        # Y CHANNEL
        y_chan <- channels_y[w]
        # CONTRUCT PLOT
        p <- cyto_plot(
          x,
          parent = parent,
          select = select,
          channels = c(x_chan, y_chan),
          overlay = overlay,
          layout = layout,
          header = header[cnt + ceiling(w/np)],
          page = if(w == length(channels_y)) {
            TRUE
          } else {
            FALSE
          },
          ...
        )
        # PREPARE RECORDED PLOTS
        p[LAPPLY(p, "is.null")] <- NULL
        return(p)
      })
    })
  }
  
  # FORMAT RECORDED PLOTS
  plots[LAPPLY(plots, "is.null")] <- NULL
  if(length(plots) == 0) {
    plots <- NULL
  }
  
  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  return(plots)
  
}