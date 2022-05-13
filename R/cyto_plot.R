## CYTO_PLOT -------------------------------------------------------------------

#' cyto_plot
#'
#' Explore and visualise cytometry data.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied, set to the \code{"root"} node by
#'   default.
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
#'   to transform the channels of the supplied data. \code{cyto_plot} does not
#'   support in-line transformations and as such the transformations should be
#'   applied to the data prior to plotting. The transformerList is used
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
#' @param events numeric to control the number or percentage of events to
#'   display. Values [0,1] indicate the percentage of events to display (i.e.
#'   value of 1 will display all events), whilst values larger than 1 indicate
#'   the number of events to display. The default value for \code{events} is set
#'   to 50000 to display 50000 events only.
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param margins a vector of length 4 to control the margins around the bottom,
#'   left, top and right of the plot, set to c(NA, NA, NA, NA) by default to let
#'   \code{cyto_plot} compute optimal margins.
#' @param popup logical indicating whether the plot should be constructed in a
#'   pop-up window, set to TRUE by default. \code{popup} will open OS-specific
#'   graphic device prior to plotting. Mac users will need to install
#'   \href{https://www.xquartz.org/}{XQuartz} for this functionality.
#' @param popup_size a vector of length 2 to control the height and width of
#'   pop-up graphics device in inches, set to \code{c(10,10)} by default.
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
#'   Setting \code{negate = TRUE} will result in the creation of a boolean
#'   filter that contains all the events outside the gates supplied to
#'   \code{alias} or \code{gate}. If such a boolean gate exists in the
#'   GatingSet/GatingHierarchy it will automatically be extracted when
#'   \code{alias = ""}. In order to prevent plotting of these boolean gates,
#'   users will need to explicitly pass the names of the gates they want to
#'   display to \code{alias}.
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
#' @param axes_limits options include \code{"auto"}, \code{"trim"},
#'   \code{"data"} or \code{"machine"} to use optimised, lower trimmed (remove
#'   lowest 1% of events), data or machine limits respectively. Set to
#'   \code{"auto"} by default to use optimised axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param axes_limits_buffer decimal indicating the percentage of buffering to
#'   add to either end of the axes limits, set to 0.03 by default.
#' @param axes_text logical vector of length 2 indicating whether axis text
#'   should be included for the x and y axes respectively, set to
#'   \code{c(TRUE,TRUE)} by default to display axes text on both axes.
#'   \code{axes_text} can also be set to NA to remove the text only but keep the
#'   axes ticks.
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
#' @param key options include \code{"scale"} to display only the colour scale,
#'   \code{"both"} to display the colour scale and text to annotate
#'   \code{counts} or \code{fluorescence units}, or \code{"none"} to remove the
#'   key from the plot(s). Set to \code{"both"} by default to display colour
#'   scale and associated text in the key.
#' @param key_scale can be either \code{"fixed"} to use the same colour scale
#'   for points between plots, or \code{"free"} to set a unique colour scale for
#'   each plot. \code{key_scale} is set to \code{"fixed"} by default create a
#'   shared colour scale which allows comparison of relative counts between
#'   plots. \code{key_scale} also accepts a vector containing the range of
#'   values for the scale, or a list of ranges for the scale in each plot.
#' @param key_text_font font to use for text in the key, set to 1 by default for
#'   plain font.
#' @param key_text_size numeric to control the size of text in the plot key, set
#'   to 0.9 by default.
#' @param key_text_col colour to use for text in the plot key, set to
#'   \code{"black"} by default.
#' @param key_text_col_alpha numeric [0, 1] to control the transparency of the
#'   text colour in the plot key, set to 1 by default to remove transparency.
#' @param key_title title to place above the key set to \code{"count"} when
#'   \code{point_col = NA} or the name of the channel/marker supplied to
#'   \code{point_col}.
#' @param key_title_text_font font to use for the key title, set to 1 by default
#'   for plain font.
#' @param key_title_text_size numeric to control the size of title text in teh
#'   key, set to 1 by default.
#' @param key_title_text_col colour to use for the key title, set to
#'   \code{"black"} by default.
#' @param key_title_text_col_alpha numeric [0, 1] to control the transparency of
#'   the key title text, set to 1 by default to remove transparency.
#' @param key_hist_line_type line type(s) to use for histogram borders in the
#'   key, set to 1 by default to use solid lines. See
#'   \code{\link[graphics:par]{lty}} for alternatives.
#' @param key_hist_line_width numeric to control line width(s) for histogram
#'   borders lines, set to 1 by default.
#' @param key_hist_line_col colour to use for the histogram border in the key,
#'   set to \code{"black"} by default.
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
#' @param label_text vector of population names to use in the labels.To exclude
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
#'   added to each full page within the \code{cyto_plot} call. For custom
#'   plotting functions, if the final page has empty panels, we need to set
#'   \code{page = TRUE} in custom plotting functions to ensure that the header
#'   is appropriately added.
#' @param header_text_font numeric to control the font of the header text, set
#'   to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param header_text_size numeric to control the size of the header text, set
#'   to 1 by default.
#' @param header_text_col colour to use for the header text, set to "black" by
#'   default.
#' @param page used internally within cyto_plot() wrapper functions to signal if
#'   a new page is required after constructing the plot(s), set to FALSE by
#'   default. This argument is only required for custom layouts which result in
#'   pages with empty panels that require a header, because cyto_plot() only
#'   adds headers to complete pages.
#' @param memory logical indicating whether \code{cyto_plot()} should remember
#'   label co-ordinates when \code{label_position = "manual"} and use those
#'   co-ordinates when \code{cyto_plot_save()} is called, set to TRUE by
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
#' @return a list containing the recorded plots.
#'
#' @export
cyto_plot <- function(x,
                      parent = "root",
                      channels,
                      alias = NA,
                      axes_trans = NA,
                      merge_by = "name",
                      overlay = NA,
                      gate = NA,
                      axes_limits = "auto",
                      events = 50000,
                      layout,
                      margins = c(NA, NA, NA, NA),
                      popup = TRUE,
                      popup_size = c(NA, NA),
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
                      hist_stack = 0,
                      hist_layers = NA,
                      hist_cols = NA,
                      hist_fill = NA,
                      hist_fill_alpha = 1,
                      hist_line_type = 1,
                      hist_line_width = 1,
                      hist_line_col = "black",
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
                      key = "both",
                      key_scale = "fixed", 
                      key_text_font = 1,
                      key_text_size = 0.9,
                      key_text_col = "black",
                      key_text_col_alpha = 1,
                      key_title = "",
                      key_title_text_font = 1,
                      key_title_text_size = 1,
                      key_title_text_col = "black",
                      key_title_text_col_alpha = 1,
                      key_line_type = 1,
                      key_line_width = 1.5,
                      key_line_col = "black",
                      gate_line_type = 1,
                      gate_line_width = 2.5,
                      gate_line_col = "red",
                      gate_fill = "white",
                      gate_fill_alpha = 0,
                      label,
                      label_text,
                      label_stat,
                      label_position = "auto",
                      label_text_x = NA,
                      label_text_y = NA,
                      label_text_font = 2,
                      label_text_size = 0.8,
                      label_text_col = "black",
                      label_text_col_alpha = 1,
                      label_fill = "white",
                      label_fill_alpha = 0.75,
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
                      page = FALSE,
                      memory = TRUE,
                      seed = 42,
                      ...) {
  
  # CYTO_PLOT_EXIT -------------------------------------------------------------
  
  # SIGNAL CALL TO CYTO_PLOT & RESET
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "cytoset")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)
  
  # REMOVE PAGE ARGUMENT
  args <- args[!names(args) %in% c("page")]
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # CONVERT MISSING TO EMPTY
  .args_update(args)
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if(.empty(args$channels)) {
    stop(
      "Supply channel/marker(s) to 'channels' to construct the plot."
    )
  } else {
    args$channels <- cyto_channels_extract(
      args$x,
      channels = args$channels,
      plot = TRUE
    )
  }
  
  # LAYOUT - CUSTOM PLOT
  if(cyto_option("cyto_plot_method") == "custom") {
    args$layout <- FALSE
  }
  
  # # LAYOUT - TURNED OFF
  # if(!all(.empty(args$layout))) {
  #   if(all(args$layout == FALSE)){
  #     cyto_option("cyto_plot_method", "custom")
  #   }
  # }
  
  # AXES_TRANS - SUPPLIED
  if(!.all_na(args$axes_trans)) {
    if(cyto_class(args$axes_trans, "transformList")) {
      message("'axes_trans' must be a transformerList!")
      args$axes_trans <- NA
    }
  # AXES_TRANS - EXTRACT
  } else {
    if(cyto_class(args$x, "GatingSet")) {
      args$axes_trans <- cyto_transformers_extract(args$x)
    }
  }
  
  # # HEADER - HEADER ADDED TO COMPLETE LAYOUT
  # if(cyto_option("cyto_plot_method") == "custom") {
  #   # CANNOT SET HEADER FOR CUSTOM LAYOUTS
  #   args$header <- NA
  # }
  
  # GATES PREPARATION ----------------------------------------------------------
  
  # PASS GATES MANUALLY IF DATA ALREADY PREPARED BY .CYTO_PLOT_DATA()
  # CANNOT EXTRACT GATES FROM CYTOSETS
  
  # GATE - LIST OF GATE OBJECT LISTS
  args$gate <- cyto_func_execute(".cyto_plot_gates", args)
  
  # DATA PREPARATION -----------------------------------------------------------
  
  # X - LIST OF CYTOSET LISTS - NO BARCODING
  args$x <- cyto_func_execute(".cyto_plot_data", args)
  
  # REMOVE DATA ARGUMENTS
  args <- args[!names(args) %in% c("parent",
                                   "alias",
                                   "overlay",
                                   "select",
                                   "merge_by",
                                   "events", 
                                   "seed",
                                   "negate")]
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # HIST_LAYERS - DATA PREPARED IN .CYTO_PLOT_DATA()
  args$hist_layers <- LAPPLY(args$x, "length")
  
  # FORMAT GATES FOR HISTOGRAMS IN SAME PLOT
  if(length(args$channels) == 1 & !all(args$hist_layers == 1)) {
    # USE FIRST SET OF GATES PER GROUP
    if(!.all_na(args$gate)) {
      cnt <- 0
      args$gate <- structure(
        lapply(
          seq_along(args$x),
          function(z) {
            ind <- seq(cnt, cnt + length(z))
            cnt <<- cnt + length(z)
            args$gate[ind][[1]]
          }
        ),
        names = names(args$x)
      )
    # NO GATES
    } else {
      args$gate <- structure(
        rep(
          list(NA),
          length.out = length(args$x)
        ),
        names = names(args$x)
      )
    }
  }
  
  # TITLE - INHERIT NAMES
  if(all(LAPPLY(args$title, ".empty"))) {
    # REPEAT
    args$title <- rep("", length(args$x))
    # HISTOGRAMS
    if(length(args$channels) == 1) {
      args$title <- LAPPLY(
        seq_along(args$title), 
        function(z){
          paste0(args$title[z], names(args$x)[z])
        }
      )
      # GROUPS
      if(!is.null(names(args$x))) {
        args$title <- names(args$x)
      }
      # POPULATIONS - SHARED
      args$title <- LAPPLY(
        seq_along(args$x), 
        function(z){
          if(length(unique(names(args$x[[z]]))) == 1) {
            return(
             paste(
               args$title[z],
               unique(names(args$x[[z]])),
               sep = "\n")
            )
        } else {
          return(args$title[z])
        }
      })
    # POINTS  
    } else {
      # GROUPS
      if(!is.null(names(args$x))) {
        args$title <- names(args$x)
      }
      # POPULATIONS/SAMPLES
      if(!is.null(LAPPLY(args$x, "names"))) {
        args$title <- paste(
          args$title,
          LAPPLY(
            args$x, 
            function(z) {
              names(z)[1] # BASE LAYER
            }
          ),
          sep = "\n"
        )
      }
    }
  }
  
  # LEGEND_TEXT
  if(.all_na(args$legend_text)) {
    # POPULATIONS
    if(!is.null(LAPPLY(args$x, "names"))) {
      args$legend_text <- LAPPLY(args$x, function(z){
        nm <- names(z)
        nm[nm == "root"] <- "All Events"
        return(nm)
      })
    }
  }
  
  # X AXIS LIMITS - SUPPLIED ON LINEAR SCALE
  args$xlim <- .cyto_transform(
    args$xlim,
    trans = args$axes_trans,
    channel = args$channels[1],
    inverse = FALSE
  )
  
  # COMPUTE MISSING X AXIS LIMITS
  if (any(is.na(args$xlim))) {
    args$xlim[is.na(args$xlim)] <- 
      .cyto_plot_axes_limits(
        args$x,
        channels = args$channels[1],
        gate = args$gate[[1]],
        axes_limits = args$axes_limits,
        buffer = args$axes_limits_buffer
      )[, args$channels[1]][is.na(args$xlim)]
  }
  
  # Y AXIS LIMITS - SUPPLIED ON LINEAR SCALE
  if(length(args$channels) == 2){
    args$ylim <- .cyto_transform(
      args$ylim,
      trans = args$axes_trans,
      channel = args$channels[2]
    )
    if(any(is.na(args$ylim))) {
      args$ylim[is.na(args$ylim)] <- 
        .cyto_plot_axes_limits(
          args$x,
          channels = args$channels[2],
          gate = args$gate[[1]],
          axes_limits = args$axes_limits,
          buffer = args$axes_limits_buffer
        )[, args$channels[2]][is.na(args$ylim)]
    }
  }
  
  # X AXIS BREAKS & LABELS
  axes_text_x <- 
    .cyto_plot_axes_text(
      args$x[[1]],
      channels = args$channels[1],
      axes_text = args$axes_text[1],
      axes_trans = args$axes_trans,
      axes_range = list(args$xlim),
      axes_limits = args$axes_limits
    )[[1]]

  
  # Y AXIS BREAKS & LABELS
  axes_text_y <- 
    .cyto_plot_axes_text(
      args$x[[1]],
      channels = args$channels[2],
      axes_text = args$axes_text[2],
      axes_trans = args$axes_trans,
      axes_range = list(args$ylim),
      axes_limits = args$axes_limits
    )[[1]]
  
  # AXES_TEXT
  args$axes_text <- list(axes_text_x, axes_text_y)
  args$axes_text <- rep(list(args$axes_text), length.out = length(args$x))
  
  # LABEL_TEXT_X CO-ORDINATES - LINEAR -> TRANSFORMED SCALE
  args$label_text_x <- .cyto_transform(
    args$label_text_x,
    trans = args$axes_trans,
    channel = args$channels[1],
    inverse = FALSE
  )
  
  # LABEL_TEXT_Y CO-ORDINATES - LINEAR -> TRANSFORMED SCALE
  if(length(channels) == 2) {
    args$label_text_y <- .cyto_transform(
      args$label_text_y,
      trans = args$axes_trans,
      channel = args$channels[2],
      inverse = FALSE
    )
  }
  
  # LAYOUT MISSING - SET FOR MULTIPLE PLOTS ONLY
  if (all(.empty(args$layout))) {
    # LAYOUT DIMENSIONS
    args$layout <- .cyto_plot_layout(args$x, layout = args$layout)
  # LAYOUT TURNED OFF
  } else if (all(args$layout == FALSE) | .all_na(args$layout)) {
    # CHECK GLOBAL PARAMETERS
    if(any(c("layout", "mfrow", "mfcol") %in% 
           names(cyto_option("cyto_plot_par")))) {
      args$layout <- cyto_option("cyto_plot_par")[c("layout",
                                                    "mfrow",
                                                    "mfcol")][[1]]
    } else {
      # USE CURRENT DIMENSIONS
      args$layout <- par("mfrow")
    }
  }
  
  # HEADER ARGUMENTS - REPEAT PER PAGE
  header_args <- args[grepl("^header", names(args))]
  for(i in names(header_args)) {
    header_args[[i]] <- rep(
      header_args[[i]], 
      length.out = if(!is.null(dim(args$layout))) {
        max(unique(unlist(args$layout)))
      }else {
        ceiling(length(args$x)/
        prod(args$layout))
      }
    )
  }
  
  # UPDATE HEADER ARGUMENTS
  .args_update(header_args)
  
  # REMOVE HEADER ARGUMENTS - CANNOT SPLIT LATER
  args <- args[!names(args) %in% names(header_args)]
  
  # KEY_SCALE 
  if(length(channels) == 2) {
    args$key_scale <- .cyto_plot_key_scale(
      args$x,
      channels = args$channels,
      xlim = args$xlim,
      ylim = args$ylim,
      point_col = args$point_col,
      key_scale = args$key_scale,
      axes_trans = args$axes_trans
    )
  # KEY_SCALE 1D - REQUIRED FOR ARGUMENT SPLITTER
  } else {
    args$key_scale <- split(
      rep("fixed", length.out = length(args$x)),
      1:length(args$x)
    )
  }
  
  # PREPARE GRAPHICS DEVICE ----------------------------------------------------
  
  # POPUP_SIZE - INHERIT CUSTOM DEVICE SIZE
  if(.all_na(args$popup_size)) {
    if("popup_size" %in% names(cyto_option("cyto_plot_par"))) {
      args$popup_size <- cyto_option("cyto_plot_par")[["popup_size"]]
    } else {
      args$popup_size <- c(10, 10)
    }
  }
  
  # HEADER SPACE
  if(.all_na(header)) {
    oma <- .par("oma")[[1]]
  } else {
    oma <- c(0, 0, 3, 0)
  }
  
  # GRAPHICS DEVICE
  cyto_plot_new(
    args$popup,
    popup_size = args$popup_size,
    layout = args$layout,
    oma = oma
  )
  
  # RESET LABEL MEMORY
  if(cyto_option("cyto_plot_method") == "cytoset" & 
     !cyto_option("cyto_plot_save") & !memory) {
    cyto_plot_memory_reset()
  }
  args <- args[!names(args) %in% "memory"]
  
  # REMOVE LAYOUT FROM ARGUMENTS - CANNOT SPLIT BELOW
  args <- args[!names(args) %in% c("layout")]
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # REPEAT & SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args)
  
  # TRANSPOSE ARGUMENTS
  args <- lapply(seq_along(args$x), function(z) {
    # lapply(seq_along(args), function(w){
    #   args[[w]][z]
    # })
    return(lapply(args, `[[`, z))
  })
  
  # if(cyto_option("cyto_plot_method") == "cytoset" & 
  #    !cyto_option("cyto_plot_save")) {
  #   .cyto_plot_args_remove()
  # }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # PLOTTING PROGRESS BAR
  pb <- cyto_progress(
    label = "cyto_plot()",
    total = length(args)
  )
  
  # ADD PLOTS 
  cnt <- 0
  m <- list()
  plots <- lapply(
    seq_along(args),
    function(z){
      
      # PLOT COUNTER
      cnt <<- cnt + 1
      
      # PLOT ARGUMENTS
      ARGS <- args[[z]]
      class(ARGS) <- "cyto_plot"
      
      # RESORT TO AUTO LABEL POSITIONS
      if(cyto_option("cyto_plot_save") & !memory) {
        ARGS$label_position <- "auto"
      }
      
      # HISTOGRAMS -------------------------------------------------------------
      
      # HISTOGRAMS - STAT/SMOOTH/STACK
      if(length(ARGS$channels) == 1) {
        ARGS$d <- do.call(".cyto_plot_hist", ARGS)
      } else {
        ARGS$d <- NA
      }
      
      # GATES ------------------------------------------------------------------
      
      # GATE DIMENSIONS
      if(!.all_na(ARGS$gate)) {
        ARGS$gate <- cyto_gate_prepare(ARGS$gate, ARGS$channels)
      }
      
      # LABEL TEXT & STAT ARGUMENTS --------------------------------------------
      
      # LABEL ARGUMENTS
      label_args <- do.call(".cyto_plot_label_args", ARGS)
      ARGS[c("label",
             "label_text",
             "label_stat")] <- label_args[c("label",
                                            "label_text",
                                            "label_stat")]
      
      # LEGEND -----------------------------------------------------------------
      
      # LEGEND_TEXT
      if(.all_na(ARGS$legend_text)) {
        ARGS$legend_text <- LAPPLY(ARGS$x, "cyto_names") 
      }
      
      # LABEL STATISTICS -------------------------------------------------------
      
      # POPULATIONS TO LABEL
      ARGS$pops <- do.call(".cyto_plot_label_pops", ARGS)
      
      # STATISTICS
      if(ARGS$label == TRUE) {
        # STATISTICS
        ARGS$label_stat <- do.call(".cyto_plot_label_stat", ARGS)
        # COMBINE LABEL_TEXT & LABEL_STAT
        ARGS$label_text <- do.call(".cyto_plot_label_text", ARGS)
      }
      
      # INHERIT LABEL CO-ORDINATES FROM MEMORY BETWEEN PLOTS
      if(!cyto_option("cyto_plot_save") &
         ARGS$label &
         z > 1 &
         any(
           is.na(
             c(ARGS$label_text_x,
               ARGS$label_text_y)
           )
         )) {
        label_args <- c("label_text_x", "label_text_y")
        lapply(
          label_args, 
          function(arg) {
            arg_ind <- which(is.na(ARGS[[arg]]))
            ARGS[[arg]][arg_ind] <<- m[[1]][[arg]][arg_ind]
          }
        )
      }
      
      # INHERIT LABEL CO-ORDINATES WHEN SAVING
      if(cyto_option("cyto_plot_save") & 
         memory & 
         ARGS$label_position == "manual") {
        ARGS <- .cyto_plot_args_inherit(ARGS)
      }
      
      # PLOT CONSTRUCTION ------------------------------------------------------
      
      # BUILD PLOT
      ARGS <- .cyto_plot_build(ARGS)
      
      # PLOT MEMORY ------------------------------------------------------------
      
      # RECORD LABEL CO-ORDINATES
      m[[z]] <<- ARGS[c("label_text_x", "label_text_y")]
      
      # GRAPHICS DEVICE --------------------------------------------------------
      
      # FILL GRAPHICS DEVICE - NEW PAGE
      if(z == length(args) & 
         (cyto_option("cyto_plot_method") == "cytoset" | page)) {
        cyto_plot_new_page()
      }
      
      # RECORD FULL GRAPHICS DEVICE
      if(.par("page") == TRUE) {
        # HEADER - EXCLUDED FROM ARGS
        if(!.all_na(header)) {
          if(!.all_na(header[1])) {
            .cyto_plot_header(
              header[1],
              header_text_font = header_text_font[1],
              header_text_size = header_text_size[1],
              header_text_col = header_text_col[1]
            )
            header <- header[-1]
            header_text_font <- header_text_font[-1]
            header_text_size <- header_text_size[-1]
            header_text_col <- header_text_col[-1]
          }
        }
        # RECORD
        p <- cyto_plot_record()
        # NEW DEVICE
        if(z < length(args)) {
          cyto_plot_new() # USE GLOBAL SETTINGS
        }
      } else {
        p <- NULL
      }
      
      # UPDATE PROGRESS BAR
      pb <- cyto_option("CytoExploreR_progress")
      if(!is.null(pb)) {
        if(.grepl("^cyto_plot", names(pb))) {
          cyto_progress(pb[[1]])
        }
      }
      
      # RETURN RECORDED PLOT
      return(p)
      
    }
  )
  
  # SAVE MEMORY ----------------------------------------------------------------
  
  # REMEMBER LABEL CO-ORDINATES
  if(memory) {
    # COMBINE MEMORY
    m <- structure(
      lapply(
        names(m[[1]]), 
        function(z){
          m <- LAPPLY(m, `[[`, z)
          names(m) <- rep(NA, length(m))
          return(m)
        }
      ), 
      names = names(m[[1]])
    )
    # UPDATE MEMORY
    if(!cyto_option("cyto_plot_save")) {
      cyto_plot_memory <- .cyto_plot_args_recall()
      # NO MEMORY
      if(is.null(cyto_plot_memory)) {
        .cyto_plot_args_save(m)
      } else {
        .cyto_plot_args_save(
          structure(
            lapply(
              names(cyto_plot_memory), 
              function(q) {
                c(cyto_plot_memory[[q]],
                  m[[q]])
              }
            ), 
            names = names(cyto_plot_memory)
          )
        )
      }
    }
  }
  
  # RETURN RECORDED PLOTS ------------------------------------------------------
  
  # RECORDED PLOTS
  plots[LAPPLY(plots, "is.null")] <- NULL
  if(length(plots) == 0){
    plots <- NULL
  }
  invisible(plots)
  
}
