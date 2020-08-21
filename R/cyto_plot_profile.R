## CYTO_PLOT_PROFILE -----------------------------------------------------------

#' Plot expression profile in all channels
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels a vector channels to use to construct the plots, set to all
#'   channels by default.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied.
#' @param layout a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param hist_stack numeric [0,1] indicating the degree of offset for 1-D
#'   density distributions with overlay, set to 0.5 by default.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to complete axes ranges.
#' @param popup logical indicating whether a pop-up window should be opened
#'   prior to plotting, set to \code{FALSE} by default.
#' @param header character string to include in the plot header, set to the
#'   sample names by default. The header can be removed by setting this argument
#'   to NA.
#' @param header_text_font numeric indicating the font to use for the header,
#'   set to 2 by default for bold font.
#' @param header_text_size numeric to control the text size of the header, set
#'   to 1 by default.
#' @param header_text_col colour to use for the header, set to "black" by
#'   default.
#' @param ... additional arguments passed to \code{\link{cyto_plot}}.
#'
#' @importFrom grDevices n2mfrow recordPlot
#' @importFrom graphics par mtext
#' @importFrom methods is
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot}}
#'
#' @examples
#'
#' # Load in CytoExploreR to access data
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- cyto_compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs, trans_type = "logicle")
#'
#' # Apply gates
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # Plot expression profile of T Cells in all channels
#' cyto_plot_profile(gs[1:9],
#'   parent = "T Cells"
#' )
#'
#' # Group samples by Treatment & select highest OVAConc
#'  cyto_plot_profile(gs[1:9],
#'   parent = "CD4 T Cells",
#'   merge_by = "Treatment",
#'   select = list("OVAConc" = 500)
#' )
#'
#' @export
cyto_plot_profile <- function(x,
                              parent = "root",
                              channels = NULL,
                              layout = NULL,
                              hist_stack = 0.5,
                              axes_limits = "machine",
                              popup = FALSE,
                              header = NULL,
                              header_text_font = 2,
                              header_text_size = 1,
                              header_text_col = "black",
                              ...) {

  # CHECKS ---------------------------------------------------------------------

  # PLOT METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "cyto_plot_profile")
  }

  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # RESET OLD PARAMETERS
  old_pars <- .par(c("mfrow","oma"))
  on.exit(par(old_pars))
  
  # PREPARE PLOTTING SAPCE -----------------------------------------------------
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(x)
  }
  
  # POPUP
  cyto_plot_new(popup)
  
  # LAYOUT DIMENSIONS
  layout <- .cyto_plot_layout(channels, layout)
  if(all(layout == FALSE)) {
    layout <- par("mfrow")
  }else {
    par("mfrow" = layout)
  }
  
  # HEADER SPACE
  if (!.all_na(header)) {
    par(oma = c(0, 0, 3, 0))
  }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # HEADER
  if (is.null(header)) {
    header <- "Expression Profile"
  }
  
  # CONSTRUCT PLOTS
  p <- lapply(seq_along(channels), function(z) {
    # PLOT
    cyto_plot(x,
              parent = parent,
              channels = channels[z],
              title = NA,
              hist_layers = NA,
              hist_stack = hist_stack,
              axes_limits = axes_limits,
              popup = FALSE,
              ...
    )
    # HEADER
    if(z %% prod(layout) == 0 |
       z == length(channels)){
      # ADD HEADER - LAYOUT ONLY
      if (!.all_na(header) & all(layout != FALSE)) {
        mtext(header, 
              outer = TRUE, 
              font = header_text_font,
              cex = header_text_size,
              col = header_text_col)
      }
      # RECORD PLOT
      cyto_plot_record()
    }
  })

  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save")) {
    # CLOSE GRAPHICS DEVICE
    dev.off()
    # RESET CYTO_PLOT_SAVE
    options("cyto_plot_save" = FALSE)
    # RESET CYTO_PLOT_METHOD
    options("cyto_plot_method" = NULL)
  }
  
  # RETURN RECORDED PLOTS
  invisible(p)
  
}
