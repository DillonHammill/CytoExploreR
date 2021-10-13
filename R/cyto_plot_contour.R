## CYTO_PLOT_CONTOUR -----------------------------------------------------------

#' Add contour lines to cyto_plot
#'
#' @param x either an object of class
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to extract from GatingHierarchy
#'   or GatingSt objects, set to the "root" node by default.
#' @param overlay name(s) of the populations to overlay or a \code{cytoset} or
#'   \code{list of cytosets} containing populations to be overlaid onto the
#'   plot(s). This argument can be set to "children" or "descendants" when a
#'   \code{GatingSet} or \code{GatingHierarchy} to overlay all respective nodes.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param channels channels used to construct the existing cyto_plot.
#' @param events controls the number or percentage of events to display, set to
#'   1 by default to display all events.
#' @param contour_lines numeric indicating the number of levels to use for
#'   contour lines, set to 15 by default.
#' @param contour_line_type type of line to use for contour lines, set to 1 by
#'   default.
#' @param contour_line_width line width for contour lines, set to 1 by default.
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
#' @importFrom graphics contour par
#' @importFrom grDevices adjustcolor
#'
#' @examples
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
#' # Visualise T Cells
#' cyto_plot(gs[[1]],
#' parent = "T Cells",
#' channels = c("CD4", "CD8"))
#'
#' # Add contour lines
#' cyto_plot_contour(gs[[1]],
#' parent = "T Cells",
#' channels = c("CD4", "CD8"),
#' contour_lines = 15)
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
cyto_plot_contour <- function(x,
                              parent = "root",
                              channels,
                              overlay = NA,
                              select = NULL,
                              events = 1,
                              contour_lines = 15,
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
    # CYTO_PLOT NOT CALLED AFTER CONTOURS - DON'T SET CYTO_PLOT_DATA OPTION
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
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # CONTOUR ARGUMENTS
  contour_args <- names(args)[grepl("contour_", names(args))]
  
  # REPEAT CONTOUR ARGUMENTS
  lapply(
    contour_args,
    function(z) {
      args[[z]] <<- rep(args[[z]], length.out = length(args$x))
    }
  )
  
  # CONTOUR LINES --------------------------------------------------------------
  
  # ADD CONTOUR LINES
  invisible(
    lapply(
      seq_along(args$x),
      function(z) {
        # CONTOUR LINES REQUIRED
        if(args$contour_lines != 0) {
          # BKDE NOT PRE-COMPUTED
          if(!"bkde2d" %in% names(args)) {
            args$bkde2d <- list(
              counts = NA,
              bins = NA,
              bkde = NA
            )
          }
          # BKDE IS NOT SMOOTHED (TOO FEW EVENTS OR COUNTS ONLY)
          if(.all_na(args$bkde2d$bkde)) {
            # RECOMPUTE BKDE2D - POSSIBLE DUPLICATE BINNING
            args$bkde2d <- cyto_apply(
              args$x[[z]],
              "cyto_stat_bkde2d",
              input = "matrix",
              channels = args$channels,
              bins = c(250, 250),
              limits = list(.par("usr")[[1]][1:2],
                            .par("usr")[[1]][3:4]),
              copy = FALSE,
              simplify = FALSE
            )[[1]]
          }
          # ADD CONTOUR LINES- BKDE REQUIRED (BYPASS SAMPLES TOO FEW EVENTS)
          if(!.all_na(args$bkde2d$bkde)) {
            graphics::contour(
              z = args$bkde2d$bkde,
              x = args$bkde2d$bins$x,
              y = args$bkde2d$bins$y,
              add = TRUE,
              drawlabels = FALSE,
              nlevels = args$contour_lines[z],
              col = adjustcolor(
                args$contour_line_col[z], 
                args$contour_line_alpha[z]
              ),
              lwd = args$contour_line_width[z],
              lty = args$contour_line_type[z]
            )
          }
        }
      }
    )
  )
  
}
