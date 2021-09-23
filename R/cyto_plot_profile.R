## CYTO_PLOT_PROFILE -----------------------------------------------------------

#' Plot expression profile in multiple channels
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied, set to the \code{"root"} node by
#'   default.
#' @param channels names of the channels or markers in which the data should be
#'   plotted, set to all channels except \code{"Time"} and \code{"Event-ID"} by
#'   default.
#' @param overlay name(s) of the populations to overlay or a \code{cytoset},
#'   \code{list of cytosets} or \code{list of cytoset lists} containing
#'   populations to be overlaid onto the plot(s). This argument can be set to
#'   "children" or "descendants" when a \code{GatingSet} or
#'   \code{GatingHierarchy} to overlay all respective nodes.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param merge_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to "name" by default to prevent merging. To
#'   merge all samples set this argument to \code{TRUE} or \code{"all"}.
#' @param order can be either \code{"channels"} or \code{"groups"} to control
#'   the order in which the data is plotted. Setting \code{order = "channels"}
#'   will plot each group on a single page in all the supplied channels. On the
#'   other hand, setting \code{order = "groups"} will plot each channel on a
#'   separate page for all groups.
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param hist_layers numeric indicating the number of histograms to stack in
#'   each plot when there are no overlays, set to 1 histogram per plot by
#'   default. Each plot must contain the same number of histograms.
#' @param hist_stack numeric [0,1] indicating the degree of stacking for
#'   histograms, set to \code{0} by default.
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
#' @seealso \code{\link{cyto_plot_explore}}
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
#' # Expression profile in all channels per sample
#' cyto_plot_profile(
#'   gs,
#'   parent = "Live Cells",
#'   order = "channels"
#' )
#'
#' # Compare expression between samples in each channel
#' cyto_plot_profile(
#'   gs,
#'   parent = "Live Cells",
#'   order = "groups"
#' )
#' }
#'
#' @export
cyto_plot_profile <- function(x,
                              parent = "root",
                              channels = NULL,
                              overlay = NA,
                              select = NULL,
                              merge_by = "name",
                              order = "channels",
                              layout,
                              hist_layers = NA,
                              hist_stack = 0,
                              header,
                              ...) {
  
  # TODO: IMPROVE TITLES
  # TODO: HEADERS INCORRECT FOR LAYOUTS > DEFAULT
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "profile")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # PREPARE PLOT PARAMETERS ----------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x, 
      exclude = c("Time", "Event-ID")
    )
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
  
  # PREPARE GROUPS PER PLOT - OVERLAY
  if(!.all_na(overlay)) {
    # EACH GROUP IN SEPARATE PLOT
    hist_layers <- NA
    # PREPARE GROUP NAMES (GROUPS PER PLOT) - GRPS UNCHANGED
    grps_names <- lapply(
      seq_along(grps), 
      function(z) {
        names(grps)[z]
      }
    )
  # PREPARE GROUPS PER PLOT - NO OVERLAY
  } else {
    # SET HIST_LAYERS BASED ON HIST_STACK
    if(.all_na(hist_layers)) {
      if(all(hist_stack == 0)) {
        hist_layers <- 1
      } else {
        hist_layers <- length(grps)
      }
    # CHECK SUPPLIED HIST_LAYERS
    } else {
      if(length(grps) %% hist_layers != 0) {
        stop(
          "Each plot must have the same number of layers!"
        )
      }
    }
    # SPLIT GROUPS BY HIST_LAYERS
    if(hist_layers != length(grps)) {
      # ALL GROUPS - SAME PLOT
      if(hist_layers == length(grps)) {
        # PREPARE GROUP NAMES PER PLOT
        grps_names <- list(names(grps))
        # PREPARE GROUPS PER PLOT
        grps <- list(do.call("rbind", grps))
      # MULTIPLE GROUPS PER PLOT
      } else {
        # SPLIT GROUPS INTO PLOTS
        grps_ind <- split(
          seq_along(grps),
          rep(
            seq(
              length(grps) / hist_layers
            ),
            each = hist_layers
          )
        )
        # SPLIT GROUP NAMES
        grps_names <- lapply(
          grps_ind,
          function(z){
            names(grps)[z]
          }
        )
        # SPLIT GROUPS
        grps <- lapply(
          grps_ind,
          function(z){
            do.call(
              "rbind",
              grps[z]
            )
          }
        )
      }
    # ALL GROUPS PER PLOT
    } else {
      # GROUP NAMES ONE PER PLOT
      grps_names <- list(names(grps))
      # COMBINE GROUPS
      grps <- list(do.call("rbind", grps))
    }
  }
  
  # PREPARE LAYOUT - IS THIS CORRECT?
  if(missing(layout)) {
    # SWITCH LAYOUTS FOR SINGLE PLOTS
    if(length(channels) == 1) {
      order <- "groups"
    } else if(length(grps) == 1) {
      order <- "channels"
    }
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
       layout <- .cyto_plot_layout(channels)
    # GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(grps)
    }
  }
    
  # TOTAL PLOTS PER PAGE - BASED ON LAYOUT
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES
  if(grepl("^c", order, ignore.case = TRUE)) {
    n <- length(grps)
    pg <- ceiling(length(channels)/np)
    tpg <- n * pg
  } else {
    n <- length(channels)
    pg <- ceiling(length(grps)/np)
    tpg <- n * pg
  }
  
  # PREPARE HEADERS
  if(missing(header)) {
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      # PREPARE HEADER BASED ON GROUPS PER PLOT
      if(!.all_na(hist_layers)) {
        # SINGLE GROUP PER PLOT
        if(hist_layers == 1) {
          header <- rep(unlist(grps_names), each = pg)
        # MULTIPLE GROUPS PER PLOT
        } else {
          header <- paste0(
            "Expression Profile - ", 
            rep(
              seq_len(n),
              each = pg
            )
          )
        }
      # SINGLE GROUP PER PLOT
      } else {
        header <- rep(unlist(grps_names), each = pg)
      }
    # GROUP ORDER
    } else {
      header <- rep(
        cyto_markers_extract(
          x,
          channels = channels,
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
    # CUSTOM HEADERS
    } else {
      # REPEAT HEADERS - MULTIPLE PAGES
      if(length(header) == tp) {
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
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # CALL CYTO_PLOT IN CHANNEL ORDER
  if(grepl("^c", order, ignore.case = TRUE)) {
    # CONSTRUCT & RECORD PLOTS
    plots <- lapply(seq_along(grps), function(z){
      # HEADER COUNTER
      cnt <- (z - 1) * pg
      # RECORDED PLOTS PER GROUP
      p <- lapply(seq_along(channels), function(w){
        # CONSTRUCT PLOT
        cyto_plot(
          x,
          parent = parent,
          select = grps[[z]][, "name"],
          channels = channels[w],
          overlay = overlay,
          merge_by = vars,
          layout = layout,
          hist_layers = hist_layers,
          hist_stack = hist_stack, 
          header = header[cnt + ceiling(w/np)],
          # title = NA,
          page = if(w == length(channels)) {
            TRUE
          } else {
            FALSE
          },
          ...
        )
      })
      # PREPARE RECORDED PLOTS
      p[LAPPLY(p, "is.null")] <- NULL
      return(p)
    })
  # CALL CYTO_PLOT IN GROUP ORDER 
  } else {
    # CONSTRUCT & RECORD PLOTS
    plots <- lapply(seq_along(channels), function(z){
      # HEADER COUNTER
      cnt <- (z - 1) * pg
      # RECORDED PLOTS PER GROUP
      p <- lapply(seq_along(grps), function(w){
        # CONSTRUCT PLOT
        cyto_plot(
          x,
          parent = parent,
          select = grps[[w]][, "name"],
          channels = channels[z],
          overlay = overlay,
          merge_by = vars,
          layout = layout,
          hist_layers = hist_layers,
          hist_stack = hist_stack, 
          header = header[cnt + ceiling(w/np)],
          # title = NA,
          page = if(w == length(grps)) {
            TRUE
          } else {
            FALSE
          },
          ...
        )
      })
      # PREPARE RECORDED PLOTS
      p[LAPPLY(p, "is.null")] <- NULL
      return(p)
    })
  }
  
  # FORMAT RECORDED PLOTS
  plots[LAPPLY(plots, "is.null")] <- NULL
  if(length(plots) == 0) {
    plots <- NULL
  }

  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}