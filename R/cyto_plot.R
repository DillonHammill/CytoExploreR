## CYTO_PLOT -------------------------------------------------------------------

# GatingHierarchy/GatingSet methods are almost identical and could probably be
# combined into a single method for simplicity. select and merge_by arguments
# are required in GatingHiewrarchy method as it calls cytoset method now.

#' cyto_plot
#'
#' Explore and visualise cytometry data.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoframe]{cytoframe}},
#'   \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied.
#' @param alias name of the gated population(s) to gated in the plot when a
#'   \code{GatingHierarchy} or \code{GatingSet} object is supplied. Setting
#'   \code{alias} to "" will automatically plot any gates constructed in the
#'   supplied channels. \code{alias} is equivalent to the \code{gate} argument
#'   for \code{cytoset} objects.
#' @param channels name of the channel(s) or marker(s) to be used to construct
#'   the plot. The length of channels determines the type of plot to be
#'   constructed, either a histogram for a single channel or a 2-D scatterplot
#'   with density colour scale for two channels.
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied flowFrame. \code{cyto_plot} does
#'   not support in-line transformations and as such the transformations should
#'   be applied to the data prior to plotting. The transformerList is used
#'   internally to ensure that the axes on the constructed plots are
#'   appropriately labelled.
#' @param merge_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to "name" by default to prevent merging. To
#'   merge all samples set this argument to \code{TRUE} or \code{"all"}.
#' @param overlay name(s) of the populations to overlay or a \code{cytoset},
#'   \code{list of cytosets} or \code{list of cytoset lists} containing
#'   populations to be overlaid onto the plot(s). This argument can be set to
#'   "children" or "descendants" when a \code{GatingSet} or
#'   \code{GatingHierarchy} to overlay all respective nodes.
#' @param gate gate objects to be plotted, can be either objects of class
#'   \code{rectangleGate}, \code{polygonGate}, \code{ellipsoidGate},
#'   \code{quadGate} or \code{filters}. Lists of these supported gate objects
#'   are also supported.
#' @param display numeric to control the number or percentage of events to
#'   display. Values [0,1] indicate the percentage of events to display (i.e.
#'   value of 1 will display all events), whilst values larger than 1 indicate
#'   the number of events to display. The default value for \code{display} is
#'   set to 25000 to display 25000 events only.
#' @param layout a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param margins a vector of length 4 to control the margins around the bottom,
#'   left, top and right of the plot, set to c(NA, NA, NA, NA) by default to let
#'   \code{cyto_plot} compute optimal margins.
#' @param popup logical indicating whether the plot should be constructed in a
#'   pop-up window, set to TRUE by default. \code{popup} will open OS-specific
#'   graphic device prior to plotting. Mac users will need to install
#'   \href{https://www.xquartz.org/}{XQuartz} for this functionality.
#' @param popup_size a vector of length 2 to control the height and width of
#'   pop-up graphics device in inches, set to \code{c(7,7)} by default.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details. Sample selection occurs prior to grouping with \code{merge_by}.
#' @param xlim lower and upper limits of x axis (e.g. c(0,250000)).
#' @param ylim lower and upper limits of y axis (e.g. c(0,250000)).
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param title title to use for the plot, set to the name of the sample by
#'   default. Title can be removed by setting this argument to \code{NA}.
#' @param negate logical indicating whether a label should be included for the
#'   negated population when gate objects are supplied, set to FALSE by default.
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
#'   histograms, set to \code{0} by default.
#' @param hist_layers numeric indicating the number of histograms to stack in
#'   each plot when there are no overlays, set to 1 histogram per plot by
#'   default. Each plot must contain the same number of histograms.
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
#'   contour lines in 2-D scatter plots, set to 0 by default to turn off contour
#'   lines.
#' @param contour_line_type integer [0,6] to control the line type of contour
#'   lines in 2-D scatter plots, set to \code{1} to draw solid lines by default.
#'   See \code{\link[graphics:par]{lty}} for alternatives.
#' @param contour_line_width numeric to control line width(s) for contour lines
#'   in 2-D scatter plots, set to 2 by default.
#' @param contour_line_col colour(s) to use for contour lines in 2-D scatter
#'   plots, set to \code{"black"} by default.
#' @param contour_line_alpha numeric [0,1] to control the transparency of
#'   contour lines, set to 1 by default to remove transparency.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"auto"} by default to use optimised axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param axes_limits_buffer decimal indicating the percentage of buffering to
#'   add to either end of the axes limits, set to 0.03 by default.
#' @param axes_text logical vector of length 2 indicating whether axis text
#'   should be included for the x and y axes respectively, set to
#'   \code{c(TRUE,TRUE)} by default to display axes text on both axes.
#' @param axes_text_font numeric to control the font of axes text, set to 1 for
#'   plain font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param axes_text_size numeric to control the size of axes text, set to 1 by
#'   default.
#' @param axes_text_col colour to use for axes text, set to \code{"black"} by
#'   default.
#' @param axes_label_text_font numeric to control the font axes labels, set to 1
#'   for plain font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param axes_label_text_size numeric to control the text size of axes labels,
#'   set to 1.1 by default.
#' @param axes_label_text_col colour to use for axes labels text, set to
#'   \code{"black"} by default.
#' @param title_text_font numeric to control the font of title text, set to 2
#'   for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param title_text_size numeric to control the text size of the plot title,
#'   set to 1.1 by default.
#' @param title_text_col colour to use for plot title text, set to
#'   \code{"black"} by default.
#' @param legend can be either \code{"line"} or \code{"fill"} to indicate
#'   whether a legend should be constructed based on the histogram \code{"line"}
#'   or \code{"fill"}, set to FALSE by default to remove the legend.
#' @param legend_text vector of labels to use in the legend.
#' @param legend_text_font numeric to control the font of legend text, set to 1
#'   for plain font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param legend_text_size numeric to control the size of text in the legend,
#'   set to 1 by default.
#' @param legend_text_col colour(s) to use for text in legend, set to
#'   \code{"black"} by default.
#' @param legend_line_type numeric to control the line type for line legends,
#'   set to 1 by default. Refer to \code{lty} in \code{\link[graphics:par]{par}}
#'   for alternatives.
#' @param legend_line_width numeric to control the line width in line legend,
#'   set to 1 by default. Refer to \code{lwd} in \code{\link[graphics:par]{par}}
#'   for alternatives.
#' @param legend_line_col colour(s) to use for the lines in 1-D plot legends
#'   when legend is set to \code{"line"}.
#' @param legend_box_fill fill colour(s) to use for the boxes in 1-D plot
#'   legends when legend is set to \code{"fill"}.
#' @param legend_point_col colour(s) to use for points in 2-D scatter plot
#'   legend.
#' @param gate_line_type integer [0,6] to control the line type of gates, set to
#'   \code{1} to draw solid lines by default. See
#'   \code{\link[graphics:par]{lty}} for alternatives.
#' @param gate_line_width numeric to control the line width(s) of gates, set to
#'   \code{2.5} by default.
#' @param gate_line_col colour(s) to use for gates, set to \code{"red"} by
#'   default.
#' @param gate_fill fill colour(s) to use for gates, set to "white by default.
#' @param gate_fill_alpha numeric to control the fill transparency of gates, set
#'   to 0 by default to remove fill colour(s).
#' @param label logical indicating whether gated populations should be labelled.
#'   To include the names of the populations in these labels, supply the
#'   population names to the \code{label_text} argument. The default statistic
#'   is \code{"percent"} for gated data. This argument must be set to TRUE in
#'   order to add labels with gates.
#' @param label_text vector of population names to use in the labels.The exclude
#'   the population names set this argument to NA.
#' @param label_stat indicates the type of statistic to include in the plot
#'   labels, can be \code{"percent"}, \code{"count"}, \code{"mean"},
#'   \code{"median"}, \code{"mode"} or \code{"geo mean"}, set to
#'   \code{"percent"} for gated data or \code{NA} to exclude statistics for
#'   un-gated data. Currently, only \code{"percent"} and \code{"count"} are
#'   supported for 2-D scatter plots.
#' @param label_position either "auto" or "manual". The "auto" option (default)
#'   positions labels will be placed in the center of gates and offset if
#'   necessary. The "manual" option will allow label positioning by mouse click.
#'   Label positions are set on a per gate basis, all samples in the same group
#'   will have the same label positions. To individually label plots users must
#'   manually supply the co-ordinates to label_text_x and label_text_y.
#' @param label_text_x vector of x co-ordinate(s) to manually adjust the
#'   position plot label(s) on the plot. To interactively position labels set
#'   either \code{label_text_x} or \code{label_text_y} to "select".
#' @param label_text_y vector of y co-ordinate(s) to manually adjust the
#'   position plot label(s) on the plot. To interactively position labels set
#'   either \code{label_text_x} or \code{label_text_y} to "select".
#' @param label_text_font numeric to control the font of text in plot labels,
#'   set to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param label_text_size numeric to control the size of text in the plot
#'   labels, set to 1 by default.
#' @param label_text_col colour(s) to use for text in plot labels, set to
#'   \code{"black"} by default.
#' @param label_text_col_alpha numeric [0, 1] to control the transparency of the
#'   text colour, set to 1 by default to remove transparency.
#' @param label_fill fill colour(s) to use for labels, set to "white" by
#'   default.
#' @param label_fill_alpha numeric to control background fill transparency of
#'   label, set to 0.6 by default to introduce some transparency.
#' @param border_line_type integer [0,6] to control the line type of plot
#'   border, set to \code{1} by default for a solid border. See
#'   \code{\link[graphics:par]{lty}} for alternatives.
#' @param border_line_width numeric to control line width for the plot border,
#'   set to 1 by default.
#' @param border_line_col colour to use for the plot border, set to "black" by
#'   default.
#' @param border_fill border_fill fill colour to use inside the plot border
#'   (i.e. background colour), set to "white" by default.
#' @param border_fill_alpha transparency to use for border_fill colour, set to 1
#'   by default for no transparency.
#' @param grid logical indicating whether to include grid lines in the plot
#'   background, set to TRUE by default. Alternatively, users can supply a
#'   integer to indicate the number of equally spaced quantiles to used for the
#'   grid lines.
#' @param grid_line_type integer [0,6] to control the line type of grid lines,
#'   set to \code{1} to draw solid lines by default. See
#'   \code{\link[graphics:par]{lty}} for alternatives.
#' @param grid_line_width numeric to control the line width(s) of grid lines,
#'   set to \code{1} by default.
#' @param grid_line_col colour to use for grid lines, set to \code{"grey95"} by
#'   default.
#' @param grid_line_alpha numeric [0,1] to control the transparency of grid
#'   lines, set to 1 by default to remove transparency.
#' @param header text to include above the plots when plotting a \code{cytoset}
#'   or \code{GatingSet}. \code{cyto_plot} expects a single header which will be
#'   added to each new page within the \code{cyto_plot} call.
#' @param header_text_font numeric to control the font of the header text, set
#'   to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param header_text_size numeric to control the size of the header text, set
#'   to 1 by default.
#' @param header_text_col colour to use for the header text, set to "black" by
#'   default.
#' @param seed numeric passed to \code{\link{set.seed}} to ensure that the same
#'   sampling is applied with each \code{\link{cyto_plot}} call, set to an
#'   arbitrary numeric by default. This behaviour can be turned off by setting
#'   this argument to NULL.
#' @param ... not currently in use.
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
#' # 2-D scatter plot with overlay & Gates
#' cyto_plot(gs[1:9],
#'   parent = "CD4 T Cells",
#'   alias = "CD69+ CD4 T Cells",
#'   channels = c("Alexa Fluor 647-A", "7-AAD-A"),
#'   overlay = "CD8 T Cells"
#' )
#'
#' # 2-D Scatter Plots with Back-Gating & Gates
#' cyto_plot(gs[1:9],
#'   parent = "T Cells",
#'   alias = c("CD4 T Cells", "CD8 T Cells"),
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   overlay = c("CD69+ CD4 T Cells", "CD69+ CD8 T Cells")
#' )
#' @importFrom graphics par
#' @importFrom grDevices recordPlot
#' @importFrom openCyto gh_generate_template
#' @importFrom flowWorkspace gh_pop_is_negated gs_pop_get_children
#'   gh_pop_get_descendants gh_pop_get_children gh_pop_get_gate
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_plot
NULL

#' @noRd
#' @export
cyto_plot <- function(x, ...) {
  UseMethod("cyto_plot")
}

#' @rdname cyto_plot
#' @export
cyto_plot.GatingSet <- function(x,
                                parent= "root",
                                alias = NA,
                                channels,
                                axes_trans = NA,
                                merge_by = "name",
                                overlay = NA,
                                gate = NA,
                                display = 25000,
                                layout,
                                margins = c(NA, NA, NA, NA),
                                popup = TRUE,
                                popup_size = c(7,7),
                                select = NULL,
                                xlim = c(NA, NA),
                                ylim = c(NA, NA),
                                xlab,
                                ylab,
                                title,
                                negate,
                                hist_stat = "percent",
                                hist_bins = 256,
                                hist_smooth = 1,
                                hist_stack = 0,
                                hist_layers = 1,
                                hist_cols = NA,
                                hist_fill = NA,
                                hist_fill_alpha = 1,
                                hist_line_type = 1,
                                hist_line_width = 1,
                                hist_line_col = "black",
                                point_shape = ".",
                                point_size = 2,
                                point_col_scale = NA,
                                point_cols = NA,
                                point_col = NA,
                                point_col_alpha = 1,
                                contour_lines = 0,
                                contour_line_type = 1,
                                contour_line_width = 1,
                                contour_line_col = "black",
                                contour_line_alpha = 1,
                                axes_limits = "auto",
                                axes_limits_buffer = 0.03,
                                axes_text = c(TRUE, TRUE),
                                axes_text_font = 1,
                                axes_text_size = 1,
                                axes_text_col = "black",
                                axes_label_text_font = 1,
                                axes_label_text_size = 1.1,
                                axes_label_text_col = "black",
                                title_text_font = 2,
                                title_text_size = 1.1,
                                title_text_col = "black",
                                legend = FALSE,
                                legend_text = NA,
                                legend_text_font = 1,
                                legend_text_size = 1,
                                legend_text_col = "black",
                                legend_line_type = NA,
                                legend_line_width = NA,
                                legend_line_col = NA,
                                legend_box_fill = NA,
                                legend_point_col = NA,
                                gate_line_type = 1,
                                gate_line_width = 2.5,
                                gate_line_col = "red",
                                gate_fill = "white",
                                gate_fill_alpha = 0,
                                label,
                                label_text = NA,
                                label_stat = "",
                                label_position = "auto",
                                label_text_x = NA,
                                label_text_y = NA,
                                label_text_font = 2,
                                label_text_size = 0.8,
                                label_text_col = "black",
                                label_text_col_alpha = 1,
                                label_fill = "white",
                                label_fill_alpha = 0.6,
                                border_line_type = 1,
                                border_line_width = 1,
                                border_line_col = "black",
                                border_fill = "white",
                                border_fill_alpha = 1,
                                grid = TRUE,
                                grid_line_type = 1,
                                grid_line_width = 1,
                                grid_line_col = "grey95",
                                grid_line_alpha = 1,
                                header = NA,
                                header_text_font = 2,
                                header_text_size = 1,
                                header_text_col = "black",
                                seed = 42, 
                                ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (missing(channels)) {
    stop("Supply the channel/marker(s) to construct the plot.")
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # GATINGSET - gs (x available for flowSet method call)
  gs <- x
  
  # SELECT
  if(!is.null(select)){
    gs <- cyto_select(gs, select)
  }
  
  # GATINGHIERARCHY
  gh <- gs[[1]]
  
  # TRANSFORMATIONS 
  axes_trans <- cyto_transformers_extract(gs)
  
  # PREPARE DATA & ARGUMENTS ---------------------------------------------------
  
  # EXTRACT PARENT POPULATIONS
  x <- cyto_data_extract(gs, 
                         parent = parent,
                         copy = FALSE)[[1]]
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(gs)
  
  # NO GATES IN GATINGSET - ALIAS NA
  if(!.all_na(alias) & length(cyto_nodes(gs)) == 1){
    message(
      "GatingSet contains no gates - setting 'alias' to NA."
    )
    alias <- NA
  }
  
  # PREPARE GATINGTEMPLATE (PARENT ENTRIES ONLY)
  if(!.all_na(alias)){
    gt <- gh_generate_template(gh)
    gt <- gt[basename(gt$parent) == parent, ]
  }
  
  # EMPTY ALIAS - BOOLEAN FILTERS NOT SUPPORTED (LACK CHANNELS)
  if (.empty(alias)) {
    # 2D PLOT - BOTH CHANNELS MATCH
    if (length(channels) == 2) {
      alias <- gt$alias[gt$dims == paste(channels, collapse = ",") |
                          gt$dims == paste(rev(channels), collapse = ",")]
      # 1D PLOT - ONE CHANNEL MATCH
    } else if (length(channels) == 1) {
      ind <- lapply(gt$dims, function(z) {
        grep(channels, z)
      })
      ind <- LAPPLY(ind, "length") != 0
      alias <- gt$alias[ind]
    }
    # NO ALIAS IN SUPPLIED CHANNELS
    if (length(alias) == 0) {
      alias <- NA
    }
    # BOOL GATE CHECK - BOOL GATES HAVE NO CHANNELS IN TEMPLATE
    if (!.all_na(alias)) {
      # CHECK FOR BOOL GATES
      if (any(LAPPLY(gt[!gt$alias %in% alias, "dims"], ".empty"))) {
        # EMPTY CHANNELS ALIAS INDEX
        ind <- which(LAPPLY(gt[!gt$alias %in% alias, "dims"], ".empty"))
        # PULL OUT ALIAS
        empty_alias <- gt$alias[!gt$alias %in% alias][ind]
        # ADD VALID BOOL GATES TO ALIAS
        valid_bool_gate <- LAPPLY(empty_alias, function(z) {
          # EXTRACT GATE
          g <- gh_pop_get_gate(gh, 
                               cyto_nodes_convert(gh, 
                                                  nodes = z, 
                                                  anchor = parent))
          # BOOL GATE
          if (cyto_class(g, "booleanFilter")) {
            # BOOLEAN LOGIC
            bool <- g@deparse
            # ONLY NOT AND BOOL GATES SUPPORTED
            if(!grepl("!", bool)){
              message("Only NOT boolean gates are supported.")
              return(FALSE)
            }else if(grepl("|", bool, fixed = TRUE)){
              message("Only NOT AND boolean gates are supported.")
              return(FALSE)
              # NOT AND GATE - CORRECT ALIAS
            } else {
              # STRIP &
              bool_alias <- strsplit(bool, "&!")[[1]]
              # STRIP !
              bool_alias <- unlist(strsplit(bool_alias, "!"))
              # REMOVE EMPTY ALIAS
              bool_alias <- bool_alias[!LAPPLY(bool_alias, ".empty")]
              # BOOL ALIAS MUST BE IN ALIAS
              if(all(bool_alias %in% alias)){
                return(TRUE)
              }else{
                return(FALSE)
              }
            }
            # NOT A BOOL GATE
          } else {
            return(FALSE)
          }
        })
        if (any(valid_bool_gate)) {
          # UPDATE ALIAS
          alias <- c(alias, empty_alias[which(valid_bool_gate)])
          # TURN ON NEGATE
          negate <- TRUE
        }
      }
    }
    # ALIAS MANUALLY SUPPLIED - CONTAINS BOOLEAN FILTER (MUST HAVE ALL ALIAS)
  }else if(!.all_na(alias)){
    # ALIAS MAY BE BOOL GATE
    if(any(LAPPLY(gt[gt$alias %in% alias, "dims"], ".empty"))){
      # EMPTY CHANNELS ALIAS INDEX
      ind <- which(LAPPLY(alias, function(z){
        .empty(gt[gt$alias == z, "dims"])}))
      empty_alias <- alias[ind]
      # VALID BOOL GATE - ALIAS CORRECT
      lapply(empty_alias, function(z){
        # EXTRACT GATE
        g <- gh_pop_get_gate(gh, 
                             cyto_nodes_convert(gh, 
                                                nodes = z, 
                                                anchor = parent))
        # BOOL GATE
        if(cyto_class(g, "booleanFilter")){
          # BOOLEAN LOGIC
          bool <- g@deparse
          # ONLY NOT AND BOOL GATES SUPPORTED
          if(!grepl("!", bool)){
            message("Only NOT boolean gates are supported.")
            # REMOVE FROM ALIAS
            alias <<- alias[-match(z, alias)]
            # NEGATE
            negate <<- FALSE
          }else if(grepl("|", bool, fixed = TRUE)){
            message("Only NOT AND boolean gates are supported.")
            # REMOVE FROM ALIAS
            alias <<- alias[-match(z, alias)]
            # NEGATE
            negate <<- FALSE
            # NOT AND GATE - CORRECT ALIAS
          } else {
            # STRIP &
            bool_alias <- unlist(strsplit(bool, "&!"))
            # STRIP !
            bool_alias <- unlist(strsplit(bool_alias, "!"))
            # REMOVE EMPTY ALIAS
            bool_alias <- bool_alias[!LAPPLY(bool_alias, ".empty")]
            # BOOL ALIAS MUST BE IN ALIAS
            if(!all(bool_alias %in% alias)){
              alias <<- unique(c(alias, bool_alias))
            }
            # NEGATE 
            negate <<- TRUE
          }
          # BOOLEAN ALIAS MUST BE LAST
          alias <<- c(alias[-match(z, alias)], z)
          # NOT BOOL GATE 
        }else{
          # REMOVE FROM ALIAS
          message(paste0("Cannot plot gate ", z,"."))
          alias <<- alias[-match(z, alias)]
        }
      })
    }
  }
  
  # EXTRACT GATE OBJECTS - BYPASS BOOLEAN FILTERS - ALIAS OVERRIDES GATE
  if (!.all_na(alias)) {
    # REMOVE DUPLICATE ALIAS
    alias <- unique(alias)
    # GATES FOR EACH GATINGHIERARCHY
    gate <- lapply(seq_along(gs), function(y) {
      gt <- lapply(alias, function(w) {
        gh_pop_get_gate(gs[[y]],
                        cyto_nodes_convert(gs[[y]],
                                           nodes = w,
                                           anchor = parent))
      })
      names(gt) <- alias
      # REMOVE BOOLEAN GATES
      gt[LAPPLY(gt, function(z){cyto_class(z, "booleanFilter")})] <- NULL
      return(gt)
    })
    names(gate) <- cyto_names(gs)
    # NEGATED GATES - SINGLE NEGATED GATE
    if(all(LAPPLY(alias, function(z){
      gh_pop_is_negated(gh, cyto_nodes_convert(gh, 
                                               nodes = z, 
                                               anchor = parent))
    }))){
      # LABEL_TEXT (NA FOR GATE - ALIAS FOR LABEL)
      alias <- c(NA, alias)
    }
  }
  
  # CAPTURE OVERLAY POPULATION NAMES
  nms <- NA
  
  # OVERLAY - POPULATION NAMES
  if (!.all_na(overlay)) {
    # POPULATION NAMES TO OVERLAY
    if (is.character(overlay)) {
      # OVERLAY DESCENDANTS
      if(any(grepl("descendants", overlay))){
        overlay <- tryCatch(gh_pop_get_descendants(gh, 
                                                   parent,
                                                   path = "auto"), 
                            error = function(e){NA}) 
        # OVERLAY CHILDREN  
      }else if(any(grepl("children", overlay))){
        overlay <- tryCatch(gh_pop_get_children(gh, 
                                                parent,
                                                path = "auto"),
                            error = function(e){NA})
      }
      # LABEL_STAT WITHOUT GATES
      if (.all_na(gate) & !.empty(label_stat)) {
        if (missing(label_text)) {
          label_text <- rep(NA, length(overlay) + 1)
          ind <- which(!is.na(rep(c(
            label_stat,
            rep(NA, length(overlay) + 1)
          ),
          length.out = length(overlay) + 1
          )))
          label_text[ind] <- c(parent, overlay)[ind]
        }
      }
      # CHECK OVERLAY - MAY BE NA ABOVE
      if(!.all_na(overlay)){
        # EXTRACT POPULATIONS
        nms <- overlay
        overlay <- lapply(overlay, function(z) {
          cyto_data_extract(gs, 
                            parent = cyto_nodes_convert(gs, 
                                                        nodes = z, 
                                                        anchor = parent),
                            copy = FALSE)[[1]]
        })
        names(overlay) <- nms
      }
    }
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # NEGATE
  if (missing(negate)) {
    negate <- FALSE
  }
  
  # LAYERS
  if(.all_na(overlay)){
    # 1D PLOTS - STACKING/LAYERS
    if(length(channels) == 1){
      # HISTOGRAM LAYERS
      if(.all_na(hist_layers)){
        # MERGE_BY
        if(!.all_na(merge_by)){
          L <- length(cyto_groups(x, merge_by))
        }else{
          L <- length(x)
        }
      }else{
        L <- hist_layers
      }
    }else{
      L <- 1
    }
  }else{
    # 2 LAYERS
    if(cyto_class(overlay, "flowSet")){
      L <- 2
    }else{
      # OVERLAY LIST OF CYTOSETS
      L <- length(overlay) + 1
    }
  }
  
  # GATES
  if(.all_na(gate)){
    G <- 0
  }else{
    G <- length(gate[[1]])
    if(negate == TRUE){
      G <- G + 1
    }
  }
  
  # LABEL_TEXT - ALIAS (PER PLOT)
  if (missing(label_text)) {
    if (!.all_na(alias)) {
      # BASE LAYER ONLY
      if (length(channels) == 1 & hist_stack == 0) {
        label_text <- rep(
          c(alias, rep(NA, abs(length(alias) - G))), 
          length.out = G
        )
        # EACH LAYER
      } else if (length(channels) == 1 & hist_stack != 0) {
        label_text <- rep(
          c(alias, rep(NA, abs(length(alias) - G))),
          length.out = G * L
        )
        # BASE LAYER ONLY
      } else {
        label_text <- rep(
          c(alias, rep(NA, abs(length(alias) - G))), 
          length.out = G
        )
      }
    } else {
      label_text <- NA
    }
  }
  
  # LEGEND_TEXT
  if (.all_na(legend_text)) {
    # PARENT ONLY
    if (.all_na(overlay)) {
      legend_text <- parent
      # PARENT & OVERLAY
    } else {
      if (!.all_na(nms)) {
        legend_text <- c(parent, nms)
      } else {
        legend_text <- parent
      }
    }
  }
  
  # TITLE
  if (missing(title)) {
    title <- cyto_groups(gs, 
                         group_by = merge_by, 
                         details = FALSE)
    if(all(title == "all")){
      title <- "Combined Events"
    }
    # PARENT
    title <- LAPPLY(title, function(z) {
      if (parent == "root") {
        pt <- "All Events"
      } else {
        pt <- parent
      }
      # 1D PLOT - STACKED NO OVERLAY - LACK SAMPLENAMES
      if (length(channels) == 1 &
          .all_na(overlay)) {
        if(ifelse(.all_na(hist_layers), FALSE, hist_layers == 1)) {
          paste(z, pt, sep = "\n")
        } else {
          pt
        }
        # 1D PLOT - STACKED OVERLAY - SAMPLENAMES ONLY
      } else if (length(channels) == 1 &
                 !.all_na(overlay)) {
        z
        # PASTE SAMPLNAME & PARENT
      } else {
        paste(z, pt, sep = "\n")
      }
    })
  }
  
  # CALL CYTO_PLOT FLOWSET METHOD ----------------------------------------------
  
  # PULL DOWN ARGUMENTS
  args <- .args_list(...)
  
  # CALL CYTOSET METHOD
  .execute("cyto_plot.flowSet", args)
  
}

#' @rdname cyto_plot
#' @export
cyto_plot.GatingHierarchy <- function(x,
                                      parent = "root",
                                      alias = NA,
                                      channels,
                                      axes_trans = NA,
                                      merge_by = "name",
                                      overlay = NA,
                                      gate = NA,
                                      select = NULL,
                                      axes_limits = "auto",
                                      display = 25000,
                                      margins = c(NA, NA, NA, NA),
                                      popup = TRUE,
                                      popup_size = c(7,7),
                                      xlim = c(NA, NA),
                                      ylim = c(NA, NA),
                                      xlab,
                                      ylab,
                                      title,
                                      negate,
                                      hist_stat = "percent",
                                      hist_bins = 256,
                                      hist_smooth = 1,
                                      hist_stack = 0,
                                      hist_layers = 1,
                                      hist_cols = NA,
                                      hist_fill = NA,
                                      hist_fill_alpha = 1,
                                      hist_line_type = 1,
                                      hist_line_width = 1,
                                      hist_line_col = "black",
                                      point_shape = ".",
                                      point_size = 2,
                                      point_col_scale = NA,
                                      point_cols = NA,
                                      point_col = NA,
                                      point_col_alpha = 1,
                                      contour_lines = 0,
                                      contour_line_type = 1,
                                      contour_line_width = 1,
                                      contour_line_col = "black",
                                      contour_line_alpha = 1,
                                      axes_limits_buffer = 0.03,
                                      axes_text = c(TRUE, TRUE),
                                      axes_text_font = 1,
                                      axes_text_size = 1,
                                      axes_text_col = "black",
                                      axes_label_text_font = 1,
                                      axes_label_text_size = 1.1,
                                      axes_label_text_col = "black",
                                      title_text_font = 2,
                                      title_text_size = 1.1,
                                      title_text_col = "black",
                                      legend = FALSE,
                                      legend_text = NA,
                                      legend_text_font = 1,
                                      legend_text_size = 1,
                                      legend_text_col = "black",
                                      legend_line_type = NA,
                                      legend_line_width = NA,
                                      legend_line_col = NA,
                                      legend_box_fill = NA,
                                      legend_point_col = NA,
                                      gate_line_type = 1,
                                      gate_line_width = 2.5,
                                      gate_line_col = "red",
                                      gate_fill = "white",
                                      gate_fill_alpha = 0,
                                      label,
                                      label_text,
                                      label_stat = "",
                                      label_position = "auto",
                                      label_text_x = NA,
                                      label_text_y = NA,
                                      label_text_font = 2,
                                      label_text_size = 1,
                                      label_text_col = "black",
                                      label_text_col_alpha = 1,
                                      label_fill = "white",
                                      label_fill_alpha = 0.6,
                                      border_line_type = 1,
                                      border_line_width = 1,
                                      border_line_col = "black",
                                      border_fill = "white",
                                      border_fill_alpha = 1,
                                      grid = TRUE,
                                      grid_line_type = 1,
                                      grid_line_width = 1,
                                      grid_line_col = "grey95",
                                      grid_line_alpha = 1,
                                      header = NA,
                                      header_text_font = 2,
                                      header_text_size = 1,
                                      header_text_col = "black",
                                      seed = 42, 
                                      ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (missing(channels)) {
    stop("Supply the channel/marker(s) to construct the plot.")
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # GATINGHIERACHY - gh (x available for flowFrame method call)
  gh <- x
  
  # TRANSFORMATIONS
  axes_trans <- cyto_transformers_extract(gh)
  
  # PREPARE DATA & ARGUMENTS ---------------------------------------------------
  
  # EXTRACT PARENT POPULATION
  x <- cyto_data_extract(gh, 
                         parent = parent,
                         copy = FALSE)[[1]]
  
  # NO GATES IN GATINHIERARCHY - ALIAS NA
  if(!.all_na(alias) & length(cyto_nodes(gh)) == 1){
    message(
      "GatingHierarchy contains no gates - setting 'alias' to NA."
    )
    alias <- NA
  }
  
  # PREPARE GATINGTEMPLATE (PARENT ENTRIES ONLY)
  if(!.all_na(alias)){
    gt <- gh_generate_template(gh)
    gt <- gt[basename(gt$parent) == parent, ]
  }
  
  # EMPTY ALIAS - BOOLEAN FILTERS NOT SUPPORTED (LACK CHANNELS)
  if (.empty(alias)) {
    # 2D PLOT - BOTH CHANNELS MATCH
    if (length(channels) == 2) {
      alias <- gt$alias[gt$dims == paste(channels, collapse = ",") |
                          gt$dims == paste(rev(channels), collapse = ",")]
      # 1D PLOT - ONE CHANNEL MATCH
    } else if (length(channels) == 1) {
      ind <- lapply(gt$dims, function(z) {
        grep(channels, z)
      })
      ind <- LAPPLY(ind, "length") != 0
      alias <- gt$alias[ind]
    }
    # NO ALIAS IN SUPPLIED CHANNELS
    if (length(alias) == 0) {
      alias <- NA
    }
    # BOOL GATE CHECK - BOOL GATES HAVE NO CHANNELS IN TEMPLATE
    if (!.all_na(alias)) {
      # CHECK FOR BOOL GATES
      if (any(LAPPLY(gt[!gt$alias %in% alias, "dims"], ".empty"))) {
        # EMPTY CHANNELS ALIAS INDEX
        ind <- which(LAPPLY(gt[!gt$alias %in% alias, "dims"], ".empty"))
        # PULL OUT ALIAS
        empty_alias <- gt$alias[!gt$alias %in% alias][ind]
        # ADD VALID BOOL GATES TO ALIAS
        valid_bool_gate <- LAPPLY(empty_alias, function(z) {
          # EXTRACT GATE
          g <- gh_pop_get_gate(gh, 
                               cyto_nodes_convert(gh, 
                                                  nodes = z, 
                                                  anchor = parent))
          # BOOL GATE
          if (cyto_class(g, "booleanFilter")) {
            # BOOLEAN LOGIC
            bool <- g@deparse
            # ONLY NOT AND BOOL GATES SUPPORTED
            if(!grepl("!", bool)){
              message("Only NOT boolean gates are supported.")
              return(FALSE)
            }else if(grepl("|", bool, fixed = TRUE)){
              message("Only NOT AND boolean gates are supported.")
              return(FALSE)
              # NOT AND GATE - CORRECT ALIAS
            } else {
              # STRIP &
              bool_alias <- strsplit(bool, "&!")[[1]]
              # STRIP !
              bool_alias <- unlist(strsplit(bool_alias, "!"))
              # REMOVE EMPTY ALIAS
              bool_alias <- bool_alias[!LAPPLY(bool_alias, ".empty")]
              # BOOL ALIAS MUST BE IN ALIAS
              if(all(bool_alias %in% alias)){
                return(TRUE)
              }else{
                return(FALSE)
              }
            }
            # NOT A BOOL GATE
          } else {
            return(FALSE)
          }
        })
        if (any(valid_bool_gate)) {
          # UPDATE ALIAS
          alias <- c(alias, empty_alias[which(valid_bool_gate)])
          # TURN ON NEGATE
          if(missing(negate)){
            negate <- TRUE
          }
        }
      }
    }
    # ALIAS MANUALLY SUPPLIED - CONTAINS BOOLEAN FILTER (MUST HAVE ALL ALIAS)
  }else if(!.all_na(alias)){
    # ALIAS MAY BE BOOL GATE
    if(any(LAPPLY(gt[gt$alias %in% alias, "dims"], ".empty"))){
      # EMPTY CHANNELS ALIAS INDEX
      ind <- which(LAPPLY(alias, function(z){
        .empty(gt[gt$alias == z, "dims"])}))
      empty_alias <- alias[ind]
      # VALID BOOL GATE - ALIAS CORRECT
      lapply(empty_alias, function(z){
        # EXTRACT GATE
        g <- gh_pop_get_gate(gh, 
                             cyto_nodes_convert(gh, 
                                                nodes = z, 
                                                anchor = parent))
        # BOOL GATE
        if(cyto_class(g, "booleanFilter")){
          # BOOLEAN LOGIC
          bool <- g@deparse
          # ONLY NOT AND BOOL GATES SUPPORTED
          if(!grepl("!", bool)){
            message("Only NOT boolean gates are supported.")
            # REMOVE FROM ALIAS
            alias <<- alias[-match(z, alias)]
            # NEGATE
            negate <<- FALSE
          }else if(grepl("|", bool, fixed = TRUE)){
            message("Only NOT AND boolean gates are supported.")
            # REMOVE FROM ALIAS
            alias <<- alias[-match(z, alias)]
            # NEGATE
            negate <<- FALSE
            # NOT AND GATE - CORRECT ALIAS
          } else {
            # STRIP &
            bool_alias <- unlist(strsplit(bool, "&!"))
            # STRIP !
            bool_alias <- unlist(strsplit(bool_alias, "!"))
            # REMOVE EMPTY ALIAS
            bool_alias <- bool_alias[!LAPPLY(bool_alias, ".empty")]
            # BOOL ALIAS MUST BE IN ALIAS
            if(!all(bool_alias %in% alias)){
              alias <<- unique(c(alias, bool_alias))
            }
            # NEGATE 
            negate <<- TRUE
          }
          # BOOLEAN ALIAS MUST BE LAST
          alias <<- c(alias[-match(z, alias)], z)
          # NOT BOOL GATE 
        }else{
          # REMOVE FROM ALIAS
          message(paste0("Cannot plot gate ", z,"."))
          alias <<- alias[-match(z, alias)]
        }
      })
    }
  }
  
  # EXTRACT GATE OBJECTS - BYPASS BOOLEAN FILTERS
  if (!.all_na(alias)) {
    # REMOVE DUPLICATE ALIAS
    alias <- unique(alias)
    # GATES
    gate <- lapply(alias, function(y) {
      gh_pop_get_gate(gh, 
                      cyto_nodes_convert(gh, 
                                         nodes = y, 
                                         anchor = parent))
    })
    names(gate) <- alias
    # REMOVE BOOLEAN GATES
    ind <- which(LAPPLY(gate, function(z){cyto_class(z, "booleanFilter")}))
    gate[ind] <- NULL
    # NEGATED GATES - SINGLE NEGATED GATE
    if(all(LAPPLY(alias, function(z){
      gh_pop_is_negated(gh, 
                        cyto_nodes_convert(gh, 
                                           nodes = z, 
                                           anchor = parent))
    }))){
      # LABEL_TEXT (NA FOR GATE - ALIAS FOR LABEL)
      alias <- c(NA, alias)
    }
  }
  
  # CAPTURE OVERLAY POPULATION NAMES
  nms <- NA
  
  # OVERLAY - POPULATION NAMES
  if (!.all_na(overlay)) {
    # POPULATION NAMES TO OVERLAY
    if (is.character(overlay)) {
      # OVERLAY DESCENDANTS
      if(any(grepl("descendants", overlay))){
        overlay <- tryCatch(gh_pop_get_descendants(gh, 
                                                   parent,
                                                   path = "auto"), 
                            error = function(e){NA}) 
        # OVERLAY CHILDREN  
      }else if(any(grepl("children", overlay))){
        overlay <- tryCatch(gh_pop_get_children(gh, 
                                                parent,
                                                path = "auto"),
                            error = function(e){NA})
      }
      # LABEL_STAT WITHOUT GATES
      if (.all_na(gate) & !.empty(label_stat)) {
        if (missing(label_text)) {
          label_text <- rep(NA, length(overlay) + 1)
          ind <- which(!is.na(rep(c(
            label_stat,
            rep(NA, length(overlay) + 1)
          ),
          length.out = length(overlay) + 1
          )))
          label_text[ind] <- c(parent, overlay)[ind]
        }
      }
      # CHECK OVERLAY - MAY BE NA ABOVE
      if(!.all_na(overlay)){
        # EXTRACT POPULATIONS
        nms <- overlay
        overlay <- lapply(overlay, function(z) {
          cyto_data_extract(gh,
                            parent = cyto_nodes_convert(gh, 
                                                        nodes = z, 
                                                        anchor = parent),
                            copy = FALSE)[[1]]
        })
        names(overlay) <- nms
      }
    }
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # NEGATE
  if (missing(negate)) {
    negate <- FALSE
  }
  
  # LAYERS
  if(.all_na(overlay)){
    L <- 1
  }else{
    if(cyto_class(overlay, "flowSet")){
      L <- 2
    }else{
      # OVERLAY LIST OF CYTOSETS
      L <- length(overlay) + 1
    }
  }
  
  # GATES 
  if(.all_na(gate)){
    G <- 0
  }else{
    G <- length(gate)
    if(negate == TRUE){
      G <- G + 1
    }
  }
  
  # LABEL_TEXT - ALIAS (PER PLOT)
  if (missing(label_text)) {
    if (!.all_na(alias)) {
      # BASE LAYER ONLY
      if (length(channels) == 1 & hist_stack == 0) {
        label_text <- rep(
          c(alias, rep(NA, length(alias) - G)), 
          length.out = G
        )
        # EACH LAYER
      } else if (length(channels) == 1 & hist_stack != 0) {
        label_text <- rep(
          c(alias, rep(NA, length(alias) - G)),
          length.out = G * L
        )
        # BASE LAYER ONLY
      } else {
        label_text <- rep(
          c(alias, rep(NA, length(alias) - G)), 
          length.out = G
        )
      }
    } else {
      label_text <- NA
    }
  }
  
  # LEGEND_TEXT
  if (.all_na(legend_text)) {
    # PARENT ONLY
    if (.all_na(overlay)) {
      legend_text <- parent
      # PARENT & OVERLAY
    } else {
      if (!.all_na(nms)) {
        legend_text <- c(parent, nms)
      } else {
        legend_text <- parent
      }
    }
  }
  
  # TITLE
  if (missing(title)) {
    # SAMPLENAME
    title <- cyto_names(x)
    # PARENT
    title <- LAPPLY(title, function(z) {
      if (parent == "root") {
        pt <- "All Events"
      } else {
        pt <- parent
      }
      # 1D PLOT - STACKED NO OVERLAY - LACK SAMPLENAMES
      if (length(channels) == 1 &
          .all_na(overlay)) {
        paste(z, pt, sep = "\n")
        # 1D PLOT - STACKED OVERLAY - SAMPLENAMES ONLY
      } else if (length(channels) == 1 &
                 !.all_na(overlay)) {
        z
        # PASTE SAMPLENAME & PARENT
      } else {
        paste(z, pt, sep = "\n")
      }
    })
  }
  
  # CALL CYTO_PLOT FLOWFRAME METHOD --------------------------------------------
  
  # PULL DOWN ARGUMENTS
  args <- .args_list(...)
  
  # CALL CYTOSET METHOD
  .execute("cyto_plot.flowSet", args)
}

#' @rdname cyto_plot
#' @export
cyto_plot.flowSet <- function(x,
                              channels,
                              axes_trans = NA,
                              merge_by = "name",
                              overlay = NA,
                              gate = NA,
                              axes_limits = "auto",
                              display = 50000,
                              layout,
                              margins = c(NA, NA, NA, NA),
                              popup = TRUE,
                              popup_size = c(7,7),
                              select = NULL,
                              xlim = c(NA, NA),
                              ylim = c(NA, NA),
                              xlab,
                              ylab,
                              title,
                              negate = FALSE,
                              hist_stat = "percent",
                              hist_bins = 256,
                              hist_smooth = 1,
                              hist_stack = 0.5,
                              hist_layers = 1,
                              hist_cols = NA,
                              hist_fill = NA,
                              hist_fill_alpha = 1,
                              hist_line_type = 1,
                              hist_line_width = 1,
                              hist_line_col = "black",
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
                              axes_limits_buffer = 0.03,
                              axes_text = c(TRUE, TRUE),
                              axes_text_font = 1,
                              axes_text_size = 1,
                              axes_text_col = "black",
                              axes_label_text_font = 1,
                              axes_label_text_size = 1.1,
                              axes_label_text_col = "black",
                              title_text_font = 2,
                              title_text_size = 1.1,
                              title_text_col = "black",
                              legend = FALSE,
                              legend_text = NA,
                              legend_text_font = 1,
                              legend_text_size = 1,
                              legend_text_col = "black",
                              legend_line_type = NA,
                              legend_line_width = NA,
                              legend_line_col = NA,
                              legend_box_fill = NA,
                              legend_point_col = NA,
                              gate_line_type = 1,
                              gate_line_width = 2.5,
                              gate_line_col = "red",
                              gate_fill = "white",
                              gate_fill_alpha = 0,
                              label,
                              label_text = NA,
                              label_stat,
                              label_position = "auto",
                              label_text_x = NA,
                              label_text_y = NA,
                              label_text_font = 2,
                              label_text_size = 0.8,
                              label_text_col = "black",
                              label_text_col_alpha = 1,
                              label_fill = "white",
                              label_fill_alpha = 0.6,
                              border_line_type = 1,
                              border_line_width = 1,
                              border_line_col = "black",
                              border_fill = "white",
                              border_fill_alpha = 1,
                              grid = TRUE,
                              grid_line_type = 1,
                              grid_line_width = 1,
                              grid_line_col = "grey95",
                              grid_line_alpha = 1,
                              header = NA,
                              header_text_font = 2,
                              header_text_size = 1,
                              header_text_col = "black",
                              seed = 42,
                              ...) {
  
  # CYTO_PLOT_EXIT -------------------------------------------------------------
  
  # SIGNAL CALL TO CYTO_PLOT
  if(is.null(cyto_option("cyto_plot_method"))) {
    cyto_option("cyto_plot_method", "cytoset")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_par(reset = TRUE)
      cyto_option("cyto_plot_method", NULL)
    })
  }
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)
  
  # REMOVE OLD PARAMETERS & LAYOUT
  args <- args[!names(args) %in% c("old_pars")]
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # CONVERT MISSING TO EMPTY
  .args_update(args)
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(.empty(args$channels)) {
    stop("Supply channel/marker(s) to 'channels' to construct the plot.")
  } else {
    args$channels <- cyto_channels_extract(args$x,
                                           channels = args$channels,
                                           plot = TRUE)
  }
  
  # LAYOUT - CUSTOM PLOT
  if(cyto_option("cyto_plot_method") == "custom") {
    args$layout <- FALSE
  }
  
  # LAYOUT - TURNED OFF
  if(!all(.empty(args$layout))) {
    if(all(args$layout == FALSE)){
      cyto_option("cyto_plot_method", "custom")
    }
  }
  
  # AXES_TRANS
  if(!.all_na(args$axes_trans)) {
    if(cyto_class(args$axes_trans, "transformList")) {
      message("'axes_trans' must be a transformerList!")
      args$axes_trans <- NA
    }
  }
  
  # DATA PREPARATION -----------------------------------------------------------
  
  # LIST OF CYTOSET LISTS
  
  # BASE & OVERLAY - NO BARCODING
  args$x <- .execute(".cyto_plot_data", args)
  
  # REMOVE DATA ARGUMENTS
  args <- args[!names(args) %in% c("overlay",
                                   "display", 
                                   "select", 
                                   "seed")]
  
  # HISTOGRAM LAYERS -----------------------------------------------------------
  
  # TODO: ALLOW HIST_LAYERS TO ACCEPT GROUPING VARIABLES?
  
  # INDIVIDUAL LAYERS (NO OVERLAY)
  if(length(args$channels) == 1 & all(LAPPLY(args$x, length) == 1)) {
    # LEGEND TEXT
    if(.all_na(args$legend_text)){
      args$legend_text <- names(args$x)
    }
    # UNPACK X
    args$x <- unlist(args$x)
    # HIST_LAYERS
    if(sum(args$hist_layers) != length(x)) {
      # SAME # LAYERS PER PLOT
      if(length(x) %% args$hist_layers != 0) {
        stop(
          stop("Each plot must have the same number of layers!")
        )
      }
      # REPEAT HIST_LAYERS
      args$hist_layers <- rep(args$hist_layers,
                              length.out = length(x)/args$hist_layers)
    }
    # LAYER INDICES
    ind <- split(seq_along(x),
                 LAPPLY(seq_along(args$hist_layers), function(z){
                   rep(z, args$hist_layers[z])
                 }))
    # SPLIT LAYERS
    args$x <- lapply(unique(ind), function(z){
      args$x[z]
    })
  }
  
  # GATES ----------------------------------------------------------------------
  
  # EXPECT LIST OF GATE OBJECT LISTS - ONE PER PLOT
  if (!.all_na(args$gate)) {
    # GATE OBJECTS
    if(is(args$gate)[1] != "list") {
      # EXTRACT FILTERS
      if(is(args$gate)[1] == "filters") {
        args$gate <- rep(list(unlist(args$gate)), 
                         length.out = length(args$x))
        # LIST GATE OBJECTS
      } else {
        args$gate <- rep(list(list(unlist(args$gate))), 
                         length.out = length(args$x))
      }
      # LIST OF GATE OBJECTS/ LIST OF GATE OBJECT LISTS
    } else if(is(args$gate)[1] == "list") {
      # LIST OF GATE OBJECT LISTS
      if(is(args$gate[[1]])[1] == "list") {
        # EXTRACT FILTERS
        args$gate <- structure(lapply(args$gate, "unlist"), 
                               names = names(args$gate))
        # GATE PER SUPPLIED CYTOFRAME
        if(length(args$gate) != length(x)) {
          args$gate <- structure(rep(args$gate, length.out = length(x)),
                                 names = cyto_names(x))
        }
        # MERGE_BY
        if(args$merge_by != "name") {
          groups <- cyto_groups(x, 
                                group_by = args$merge_by, 
                                details = TRUE)
          if(is.null(names(args$gate))) {
            names(args$gate) <- cyto_names(x)
          }
          args$gate <- lapply(groups, function(z){
            args$gate[[z$name[1]]]
          })
        }
        # HIST_LAYERS - FIRST SET OF GATES
        if(length(args$channels) == 1) {
          if(.all_na(hist_layers)) {
            args$gate <- args$gate[1]
          } else {
            args$gate <- args$gate[seq(1, length(x), by = args$hist_layers)]
          }
        }
        # LIST OF GATE OBJECTS
      } else {
        # EXTRACT FILTERS & REPEAT
        args$gate <- rep(list(unlist(args$gate)), length.out = length(args$x))
      }
    }
  }
  args <- args[!names(args) %in% c("merge_by",
                                   "hist_layers")]
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # X AXIS LIMITS - SUPPLIED ON LINEAR SCALE
  args$xlim <- .cyto_transform(args$xlim,
                               trans = args$axes_trans,
                               channel = args$channels[1],
                               inverse = FALSE)
  
  # COMPUTE MISSING X AXIS LIMITS
  if (any(is.na(args$xlim))) {
    args$xlim[is.na(args$xlim)] <- 
      .cyto_plot_axes_limits(args$x,
                             channels = args$channels[1],
                             axes_limits = args$axes_limits,
                             buffer = args$axes_limits_buffer
      )[, args$channels[1]][is.na(args$xlim)]
  }
  
  # Y AXIS LIMITS - SUPPLIED ON LINEAR SCALE
  if(length(args$channels) == 2){
    args$ylim <- .cyto_transform(args$ylim,
                                 trans = args$axes_trans,
                                 channel = args$channels[2])
  }
  
  # YLIM - 1D CALCULATED LATER
  if (any(is.na(args$ylim)) & length(args$channels) == 2) {
    args$ylim[is.na(args$ylim)] <- 
      .cyto_plot_axes_limits(args$x,
                             channels = args$channels[2],
                             axes_limits = args$axes_limits,
                             buffer = args$axes_limits_buffer
      )[, args$channels[2]][is.na(args$ylim)]
  }
  
  # X AXIS BREAKS & LABELS
  if (args$axes_text[1] == TRUE) {
    axes_text_x <- 
      .cyto_plot_axes_text(args$x[[1]],
                           channels = args$channels[1],
                           axes_trans = args$axes_trans,
                           axes_range = structure(list(args$xlim, 
                                                       args$ylim),
                                                  names = rep(c(args$channels, 
                                                                NA),
                                                              length.out = 2)),
                           axes_limits = args$axes_limits
      )[[1]]
  } else {
    axes_text_x <- FALSE
  }
  
  # Y AXIS BREAKS & LABELS
  if (args$axes_text[2] == TRUE) {
    if (length(args$channels) == 2) {
      axes_text_y <- 
        .cyto_plot_axes_text(args$x[[1]],
                             channels = args$channels[2],
                             axes_trans = args$axes_trans,
                             axes_range = structure(list(args$xlim, args$ylim),
                                                    names = rep(c(args$channels, 
                                                                  NA),
                                                                length.out = 2)),
                             axes_limits = args$axes_limits
        )[[1]]
    } else if (length(args$channels) == 1) {
      axes_text_y <- NA
    }
  } else {
    axes_text_y <- FALSE
  }
  
  # AXES_TEXT
  args$axes_text <- list(axes_text_x, axes_text_y)
  args$axes_text <- rep(list(args$axes_text), length.out = length(args$x))
  
  # LABEL_TEXT_X CO-ORDINATES - LINEAR -> TRANSFORMED SCALE
  args$label_text_x <- .cyto_transform(args$label_text_x,
                                       trans = args$axes_trans,
                                       channel = args$channels[1],
                                       inverse = FALSE)
  
  # LABEL_TEXT_Y CO-ORDINATES - LINEAR -> TRANSFORMED SCALE
  if(length(channels) == 2) {
    args$label_text_y <- .cyto_transform(args$label_text_y,
                                         trans = args$axes_trans,
                                         channel = args$channels[1],
                                         inverse = FALSE)
  }
  
  # LAYOUT MISSING - SET FOR MULTIPLE PLOTS ONLY
  if (all(.empty(args$layout))) {
    # LAYOUT DIMENSIONS
    args$layout <- .cyto_plot_layout(args$x,
                                     layout = args$layout)
    # LAYOUT TURNED OFF
  } else if (all(args$layout == FALSE) | .all_na(args$layout)) {
    # USE CURRENT DIMENSIONS
    args$layout <- par("mfrow")
  }
  
  # HEADER ARGUMENTS - REPEAT PER PAGE
  header_args <- args[grepl("^header", names(args))]
  for(i in names(header_args)) {
    header_args[[i]] <- rep(header_args[[i]], 
                            length.out = ceiling(length(args$x)/
                                                   prod(args$layout)))
  }
  
  # UPDATE HEADER ARGUMENTS
  .args_update(header_args)
  
  # REMOVE HEADER ARGUMENTS
  args <- args[!names(args) %in% names(header_args)]
  
  # HEADER SPACE
  if(!.all_na(header)) {
    oma <- c(0,0,3,0)
  } else {
    oma <- par("oma")
  }
  
  # PREPARE GRAPHICS DEVICE ----------------------------------------------------
  
  # GRAPHICS DEVICE
  if(cyto_option("cyto_plot_method") == "cytoset") {
    cyto_plot_new(args$popup,
                  popup_size = args$popup_size,
                  layout = args$layout,
                  oma = oma)
  }
  
  # REMOVE LAYOUT FROM ARGUMENTS - CANNOT SPLIT BELOW
  args <- args[!names(args) %in% c("layout")]
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # REPEAT & SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args)
  
  # TRANSPOSE ARGUMENTS
  args <- lapply(seq_along(args$x), function(z) {
    return(lapply(args, `[[`, z))
  })
  
  # RESET MEMORY
  if(cyto_option("cyto_plot_method") == "cytoset" & 
     !cyto_option("cyto_plot_save")) {
    .cyto_plot_args_remove()
  }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # ADD PLOTS 
  cnt <- 0
  memory <- list()
  plots <- lapply(seq_along(args), function(z){
    
    #  PLOT COUNTER
    cnt <<- cnt + 1
    
    # PLOT ARGUMENTS
    ARGS <- args[[z]]
    class(ARGS) <- "cyto_plot"
    
    # HISTOGRAMS ---------------------------------------------------------------
    
    # HISTOGRAMS - STAT/SMOOTH/STACK
    if(length(ARGS$channels) == 1) {
      ARGS$d <- do.call(".cyto_plot_hist", ARGS)
    } else {
      ARGS$d <- NA
    }
    
    # GATES --------------------------------------------------------------------
    
    # GATE DIMENSIONS
    if(!.all_na(ARGS$gate)) {
      ARGS$gate <- cyto_gate_prepare(ARGS$gate,
                                     ARGS$channels)
    }
    
    # LABEL TEXT & STAT ARGUMENTS ----------------------------------------------
    
    # LABEL ARGUMENTS
    label_args <- do.call(".cyto_plot_label_args", ARGS)
    ARGS[c("label",
           "label_text",
           "label_stat")] <- label_args[c("label",
                                          "label_text",
                                          "label_stat")]
    
    # LEGEND -------------------------------------------------------------------
    
    # LEGEND_TEXT
    if(.all_na(ARGS$legend_text) |
       length(unique(ARGS$legend_text)) == 1) {
      ARGS$legend_text <- LAPPLY(ARGS$x, "cyto_names") 
    }
    
    # LABEL STATISTICS ---------------------------------------------------------
    
    # STATISTICS
    if(ARGS$label == TRUE) {
      # POPULATIONS
      ARGS$pops <- do.call(".cyto_plot_label_pops", ARGS)
      # STATISTICS
      ARGS$label_stat <- do.call(".cyto_plot_label_stat", ARGS)
      # COMBINE LABEL_TEXT & LABEL_STAT
      ARGS$label_text <- do.call(".cyto_plot_label_text", ARGS)
    }
    
    # INHERIT LABEL CO-ORDINATES FROM MEMORY BETWEEN PLOTS
    if(!cyto_option("cyto_plot_save") &
       z > 1 &
       any(
         is.na(
           c(ARGS$label_text_x,
             ARGS$label_text_y)
         )
       )) {
      label_args <- c("label_text_x", "label_text_y")
      lapply(label_args, function(arg){
        arg_ind <- which(is.na(ARGS[[arg]]))
        ARGS[[arg]][arg_ind] <<- memory[[1]][[arg]][arg_ind]
      })
    }
    
    # INHERIT LABEL CO-ORDINATES WHEN SAVING
    if(cyto_option("cyto_plot_save")) {
      ARGS <- .cyto_plot_args_inherit(ARGS)
    }
    
    # PLOT CONSTRUCTION --------------------------------------------------------
    
    # BUILD PLOT
    ARGS <- .cyto_plot_build(ARGS)
    
    # PLOT MEMORY --------------------------------------------------------------
    
    # RECORD LABEL CO-ORDINATES
    memory[[z]] <<- ARGS[c("label_text_x", "label_text_y")]
    
    # GARPHICS DEVICE ----------------------------------------------------------
    
    # RECORD FULL PAGE & OPEN NEW DEVICE
    if(.par("page")[[1]] | cnt == length(ARGS)) {
      # PAGE
      pg <- ceiling(z/prod(.par("mfrow")[[1]]))
      # HEADER
      if(!.all_na(header[pg])) {
        .cyto_plot_header(header[pg],
                          header_text_font = header_text_font[pg],
                          header_text_size = header_text_size[pg],
                          header_text_col = header_text_col[pg])
      }
      # RECORD
      p <- cyto_plot_record()
      # OPEN NEW DEVICE
      if(length(args) > cnt) {
        cyto_plot_new() # USE GLOBAL SETTINGS
      }
    } else {
      p <- NULL
    }
    
    # RETURN RECORDED PLOT
    return(p)
    
  })
  
  # SAVE MEMORY ----------------------------------------------------------------
  
  # COMBINE MEMORY
  memory <- structure(
    lapply(names(memory[[1]]), function(z){
      m <- LAPPLY(memory, `[[`, z)
      names(m) <- rep(NA, length(m))
      return(m)
    }), names = names(memory[[1]])
  )
  
  # UPDATE MEMORY
  if(!cyto_option("cyto_plot_save")) {
    cyto_plot_memory <- .cyto_plot_args_recall()
    # NO MEMORY
    if(is.null(cyto_plot_memory)) {
      .cyto_plot_args_save(memory)
    } else {
      .cyto_plot_args_save(
        structure(
          lapply(names(cyto_plot_memory), function(q){
            c(cyto_plot_memory[[q]],
              memory[[q]])
          }), names = names(cyto_plot_memory)
        )
      )
    }
  }
  
  # RESET SAVE & RECORD --------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if(cyto_option("cyto_plot_save") == TRUE) {
    if(cyto_option("cyto_plot_method") == "cytoset") {
      # CLOSE GRAPHICS DEVICE
      dev.off()
      # RESET CYTO_PLOT_SAVE
      cyto_option("cyto_plot_save", FALSE)
    }
  }
  
  # RECORD
  plots[LAPPLY(plots, "is.null")] <- NULL
  if(length(plots) == 0){
    plots <- NULL
  }
  invisible(plots)
  
}