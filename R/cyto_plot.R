## CYTO_PLOT -------------------------------------------------------------------

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
#'   for \code{flowFrame} and \code{flowSet} objects.
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
#' @param overlay name(s) of the populations to overlay or a \code{flowFrame},
#'   \code{flowSet}, \code{list of flowFrames}, \code{list of flowSets} or
#'   \code{list of flowFrame lists} containing populations to be overlaid onto
#'   the plot(s). This argument can be set to "children" or "descendants" when a
#'   \code{GatingSet} or \code{GatingHierarchy} to overlay all respective nodes.
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
#'   pop-up window, set to FALSE by default. \code{popup} will open OS-specific
#'   graphic device prior to plotting. Mac users will need to install
#'   \href{https://www.xquartz.org/}{XQuartz} for this functionality.
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
#' @param hist_bandwidth calculated internally based on data range and number of
#'   bins passed to \code{hist_bins}, set to NA by default and should not be
#'   changed.
#' @param hist_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust the smoothness of the kernel
#'   density for histograms, set to \code{1} by default. Only values greater or
#'   equal to 1 are supported.
#' @param hist_stack numeric [0,1] indicating the degree of stacking for
#'   histograms, set to \code{0} by default.
#' @param hist_layers numeric indicating the number of histograms to stack in
#'   each plot, set to all samples by default. Each plot must contain the same
#'   number of histograms.
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
#' @param ... not currently in use.
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load samples into GatingSet
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' gs <- transform(gs, trans)
#'
#' # Apply gatingTemplate
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
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
#' @importFrom magrittr %>%
#' @importFrom purrr transpose
#' @importFrom openCyto gh_generate_template
#' @importFrom methods formalArgs is
#' @importFrom flowWorkspace gh_pop_is_negated gs_pop_get_children
#'   gh_pop_get_descendants gh_pop_get_children gh_pop_get_gate
#' @importFrom methods is
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
                                parent,
                                alias = NA,
                                channels,
                                axes_trans = NA,
                                merge_by = "name",
                                overlay = NA,
                                gate = NA,
                                display = 25000,
                                layout,
                                margins = c(NA, NA, NA, NA),
                                popup = FALSE,
                                select = NULL,
                                xlim = c(NA, NA),
                                ylim = c(NA, NA),
                                xlab,
                                ylab,
                                title,
                                negate,
                                hist_stat = "percent",
                                hist_bins = 256,
                                hist_bandwidth = NA,
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
                                ...) {
  
  # GATINGSET METHOD - CALLS FLOWSET METHOD
  
  # CHECKS ---------------------------------------------------------------------
  
  # PARENT
  if (missing(parent)) {
    stop("Supply the name of the 'parent' population to plot.")
  }
  
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
  axes_trans <- cyto_transformer_extract(gs)
  
  # PREPARE DATA & ARGUMENTS ---------------------------------------------------
  
  # EXTRACT PARENT POPULATIONS
  x <- cyto_extract(gs, parent)
  
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
          if (is(g, "booleanFilter")) {
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
        if(is(g, "booleanFilter")){
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
      gt[LAPPLY(gt, function(z){is(z, "booleanFilter")})] <- NULL
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
          cyto_extract(gs, cyto_nodes_convert(gs, 
                                              nodes = z, 
                                              anchor = parent))
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
          L <- length(cyto_group_by(x, merge_by))
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
    if(is(overlay, "flowFrame") | is(overlay, "flowSet")){
      L <- 2
    }else{
      # OVERLAY LIST OF FLOWFRAMES OR FLOWSETS
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
  args <- .args_list()
  
  # CYTO_PLOT FLOWSET ARGUMENTS
  ARGS <- formalArgs("cyto_plot.flowSet")
  
  # RESTRICT ARGUMENTS
  args <- args[names(args) %in% ARGS]
  
  # CALL FLOWSET METHOD
  do.call("cyto_plot.flowSet", args)
}

#' @rdname cyto_plot
#' @export
cyto_plot.GatingHierarchy <- function(x,
                                      parent,
                                      alias = NA,
                                      channels,
                                      axes_trans = NA,
                                      overlay = NA,
                                      gate = NA,
                                      axes_limits = "auto",
                                      display = 25000,
                                      margins = c(NA, NA, NA, NA),
                                      popup = FALSE,
                                      xlim = c(NA, NA),
                                      ylim = c(NA, NA),
                                      xlab,
                                      ylab,
                                      title,
                                      negate,
                                      hist_stat = "percent",
                                      hist_bins = 256,
                                      hist_bandwidth = NA,
                                      hist_smooth = 1,
                                      hist_stack = 0,
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
                                      ...) {
  
  # GATINGHIERARCHY METHOD - CALLS FLOWFRAME METHOD
  
  # CHECKS ---------------------------------------------------------------------
  
  # PARENT
  if (missing(parent)) {
    stop("Supply the name of the 'parent' population to plot.")
  }
  
  # CHANNELS
  if (missing(channels)) {
    stop("Supply the channel/marker(s) to construct the plot.")
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # GATINGHIERACHY - gh (x available for flowFrame method call)
  gh <- x
  
  # TRANSFORMATIONS
  axes_trans <- cyto_transformer_extract(gh)
  
  # PREPARE DATA & ARGUMENTS ---------------------------------------------------
  
  # EXTRACT PARENT POPULATION
  x <- cyto_extract(gh, parent)
  
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
          if (is(g, "booleanFilter")) {
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
        if(is(g, "booleanFilter")){
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
    ind <- which(LAPPLY(gate, function(z){is(z, "booleanFilter")}))
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
          cyto_extract(gh,
                       cyto_nodes_convert(gh, 
                                          nodes = z, 
                                          anchor = parent))
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
    if(is(overlay, "flowFrame") | is(overlay, "flowSet")){
      L <- 2
    }else{
      # OVERLAY LIST OF FLOWFRAMES OR FLOWSETS
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
        # PASTE SAMPLNAME & PARENT
      } else {
        paste(z, pt, sep = "\n")
      }
    })
  }
  
  # CALL CYTO_PLOT FLOWFRAME METHOD --------------------------------------------
  
  # PULL DOWN ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT FLOWFRAME ARGUMENTS
  ARGS <- formalArgs("cyto_plot.flowFrame")
  
  # RESTRICT ARGUMENTS
  args <- args[names(args) %in% ARGS]
  
  # CALL FLOWFRAME METHOD
  do.call("cyto_plot.flowFrame", args)
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
                              display = 25000,
                              layout,
                              margins = c(NA, NA, NA, NA),
                              popup = FALSE,
                              select = NULL,
                              xlim = c(NA, NA),
                              ylim = c(NA, NA),
                              xlab,
                              ylab,
                              title,
                              negate = FALSE,
                              hist_stat = "percent",
                              hist_bins = 256,
                              hist_bandwidth = NA,
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
                              ...) {
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)

  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # CONVERT MISSING TO EMPTY
  .args_update(args)
  
  # CHECKS ---------------------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- par(c("mar", "mfrow", "mfcol"))
  
  # METHOD & RESET
  if (is.null(getOption("cyto_plot_method"))) {
    # SET PLOT METHOD
    options("cyto_plot_method" = "flowSet")
    # RESET PLOT METHOD & GRAPHICAL PARAMETERS ON EXIT
    on.exit({
      if (!getOption("cyto_plot_custom")) {
        par(old_pars)
      }
      options("cyto_plot_method" = NULL)
    })
  } else {
    # RESET GRAPHICAL PARAMETERS ON EXIT
    if (getOption("cyto_plot_method") == "flowSet") {
      on.exit({
        if (!getOption("cyto_plot_custom")) {
          par(old_pars)
        }
      })
    }
  }
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  if (.empty(args$channels)) {
    stop("Supply channel/marker(s) to construct the plot.")
  } else {
    args$channels <- cyto_channels_extract(args$x, 
                                           channels = args$channels, 
                                           plot = TRUE)
  }
  
  # POPUP
  if (getOption("cyto_plot_save") == TRUE) {
    args$popup <- FALSE
  }
  
  # REMOVE LAYOUT FROM ARGS
  args <- args[-match("layout", names(args))]
  
  # CUSTOM PLOT SIGNALLED
  if (getOption("cyto_plot_custom")) {
    layout <- FALSE
  }
  
  # LAYOUT TURNED OFF
  if (!all(.empty(layout))) {
    if (all(layout == FALSE)) {
      options("cyto_plot_custom" = TRUE)
    }
  }
  
  # SAMPLE PREPARATION - BASE LAYERS & OVERLAYS --------------------------------
  
  # SELECT
  if (!is.null(args$select)) {
    args$x <- cyto_select(args$x, args$select)
  }
  args <- args[-match("select", names(args))]
  
  # DATA TO LIST & GROUP - CONVERT GROUPS TO CYTOFRAMES
  args$x <- cyto_merge_by(x, args$merge_by)
  
  # OVERLAY LIST & GROUP
  if (!.all_na(args$overlay)) {
    # REPEAT FLOWFRAME PER PLOT
    if (is(args$overlay, "flowFrame")) {
      args$overlay <- rep(list(list(args$overlay)), length(arsg$x))
      # FLOWSET TO LIST OF CYTOFRAMES
    } else if (is(args$overlay, "flowSet")) {
      # GROUPING
      args$overlay <- cyto_merge_by(args$overlay, 
                                    merge_by = args$merge_by)
      # LIST OF FLOWFRAME LISTS
      args$overlay <- lapply(args$overlay, function(z) {
        list(z)
      })
      # LIST OF FLOWSETS TO LIST OF FLOWFRAME LISTS
    } else if (is(args$overlay, "list")) {
      # ALLOW LIST OF CYTOFRAMES OF LENGTH FR_LIST
      if (all(LAPPLY(unlist(args$overlay), function(z) {
        is(z, "flowFrame")
      }))) {
        # SAME LENGTH AS FR_LIST
        if (length(args$overlay) != length(args$x)) {
          stop(
            paste(
              "'overlay' must be a list of flowFrame lists -",
              "one flowFrame list per plot."
            )
          )
        }
        # LIST OF FLOWSETS
      } else if (all(LAPPLY(args$overlay, function(z) {
        is(z, "flowSet")
      }))) {
        # GROUPING
        args$overlay <- lapply(args$overlay, function(z) {
          cyto_merge_by(z, 
                        merge_by = args$merge_by)
        })
        args$overlay <- args$overlay %>% transpose()
        # OVERLAY NOT SUPPORTED
      } else {
        stop(
          paste(
            "'overlay' should be either a flowFrame, flowSet, list of flowFrame",
            "lists or list of flowSet lists."
          )
        )
      }
    }
    # COMBINE X & OVERLAY
    NMS <- names(args$x)
    args$x <- lapply(seq_len(length(args$x)), function(z) {
      c(args$x[z], args$overlay[[z]])
    })
    names(args$x) <- NMS
    # NO OVERLAY
  } else if (.all_na(args$overlay)) {
    # SAMPLENAMES
    NMS <- names(args$x)
    # LIST OF FLOWFRAME LISTS
    args$x <- lapply(seq_len(length(args$x)), function(z) {
      lst <- list(args$x[[z]])
      names(lst) <- NMS[z]
      return(lst)
    })
  }
  args <- args[-match("overlay", names(args))]
  
  # CLEANING -------------------------------------------------------------------
  
  # REMOVAL NEGATIVE FSC/SSC EVENTS - POINT_COL SCALE ISSUE
  args$x <- lapply(args$x, function(cf_list) {
    cyto_apply(cf_list, function(cf){
      lapply(seq_along(args$channels), function(z){
        if(grepl("FSC", args$channels[z], ignore.case = TRUE) |
           grepl("SSC", args$channels[z], ignore.case = TRUE)){
          if (BiocGenerics::nrow(cf) != 0) {
            if (min(range(cf, type = "data")[, args$channels[z]]) < 0) {
              coords <- matrix(c(0, Inf), 
                               nrow = 2,
                               dimnames = list(c("min", "max"),
                                               args$channels[z]))
              cf <<- Subset(cf, rectangleGate(.gate = coords))
            }
          }
        }
      })
      return(cf)
    }, simplify = FALSE, input = "cytoframe")
  })
  
  # HIST_LAYERS -------------------------------------------------------------
  # DEFAULT ALL SAMPLES SAME PLOT -> LAYERS -> STACK
  
  # GROUP SAMPLES
  if (length(args$channels) == 1 & length(args$x[[1]]) == 1) {
    # LIST(LIST OF CYTOFRAMES)
    args$x <- list(unlist(args$x))
    # LEGEND_TEXT - SAMPLENAMES
    if (.all_na(args$legend_text)) {
      args$legend_text <- NMS
    }
    # HIST_LAYERS
    if(!.all_na(hist_layers)) {
      # SAME # LAYERS PER PLOT
      if (length(args$x[[1]]) %% args$hist_layers != 0) {
        stop("Each plot must have the same number of layers!")
      }
      # INDICES
      ind <- rep(seq_len(length(args$x[[1]])),
                 each = args$hist_layers,
                 length.out = length(args$x[[1]])
      )
      # SPLITTING BY HIST_LAYERS
      args$x <- lapply(unique(ind), function(z) {
        args$x[[1]][ind == z]
      })
    }
  }
  
  # GATES PREPARATION ----------------------------------------------------------
  
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
  args <- args[-match("merge_by", names(args))]
  args <- args[-match("hist_layers", names(args))]
  
  # ARGUMENT PREPARATION -------------------------------------------------------
  
  # X AXIS LIMITS
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
  
  # Y AXIS LIMITS
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
  
  # LAYOUT MISSING - SET FOR MULTIPLE PLOTS ONLY
  if (length(args$x) > 1) {
    if (all(.empty(layout))) {
      # LAYOUT DIMENSIONS
      layout <- .cyto_plot_layout(args$x,
                                  layout = layout)
      # LAYOUT TURNED OFF
    } else if (all(layout == FALSE) | .all_na(layout)) {
      # USE CURRENT DIMENSIONS
      layout <- par("mfrow")
    }
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REPEAT & SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args)
  
  # TRANSPOSE ARGUMENTS
  args <- lapply(seq_along(args$x), function(z) {
    return(lapply(args, `[[`, z))
  })
  
  # MEMORY
  cyto_plot_memory <- .cyto_plot_args_recall()
  
  # INHERIT MEMORY - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save") == TRUE) {
    # NO MEMORY
    if (is.null(cyto_plot_memory)) {
      # RESET CYTO_PLOT_SAVE
      cyto_plot_save_reset()
      stop(
        paste(
          "cyto_plot must be called without cyto_plot_save to",
          "select label co-ordinates, once co-ordinates have been",
          "selected call cyto_plot_save and repeat the cyto_plot call."
        )
      )
    }
    # CYTO_PLOT_MEMORY INDEX
    if (is.null(names(cyto_plot_memory))) {
      memory_ind <- 1
    } else {
      memory_ind <- match(NA, names(cyto_plot_memory))
    }
    # UPDATE CYTO_PLOT_MEMORY (USED ARGS)
    names(cyto_plot_memory)[memory_ind] <- "DONE"
    # RESET CYTO_PLOT_MEMORY INDICATOR
    if (!any(is.na(names(cyto_plot_memory)))) {
      names(cyto_plot_memory) <- NULL
    }
    # INHERIT MEMORY PER PLOT
    args <- lapply(seq_along(args), function(z) {
      # ARGUMENTS
      args_to_update <- args[[z]]
      # INHERIT MEMORY
      label_args <- c("label_text_x", "label_text_y")
      lapply(label_args, function(arg) {
        label_arg_length <- length(cyto_plot_memory[[memory_ind]][[z]][[arg]])
        label_arg_memory <- cyto_plot_memory[[memory_ind]][[z]][[arg]]
        label_arg <- rep(c(args_to_update[[arg]], rep(NA, label_arg_length)),
                         length.out = label_arg_length
        )
        arg_ind <- which(is.na(label_arg))
        label_arg[arg_ind] <- label_arg_memory[arg_ind]
        args_to_update[[arg]] <<- label_arg
      })
      return(args_to_update)
    })
    # RESET MEMORY - CYTO_PLOT_SAVE
  } else if (getOption("cyto_plot_save") == FALSE) {
    .cyto_plot_args_remove()
    cyto_plot_memory <- NULL
  }
  
  # CALL CYTO_PLOT FLOWFRAME METHOD --------------------------------------------
  
  # GRAPHICS DEVICE - ONLY OPEN IF LAYOUT IS BEING SET
  if(!all(layout == FALSE)) {
    cyto_plot_new(popup, layout)
  }
  
  # PASS ARGUMENTS TO CYTO_PLOT FLOWFRAME METHOD
  cnt <- 0
  plots <- lapply(seq_along(args), function(z) {
    # COUNTER
    assign("cnt", cnt + 1, envir = parent.frame(2))
    # ARGUMENTS
    cyto_plot_args <- args[[z]]
    # PREPARE X AND OVERLAY
    if(length(cyto_plot_args$x) > 1) {
      cyto_plot_args$overlay <- cyto_plot_args$x[-1]
    } else {
      cyto_plot_args$overlay <- NA
    }
    cyto_plot_args$x <- cyto_plot_args$x[[1]]
    # TURN OFF POPUP
    cyto_plot_args$popup <- FALSE
    # LABEL CO-ORDINATE INHERITANCE BETWEEN PLOTS
    if (!getOption("cyto_plot_save") &
        any(is.na(c(
          cyto_plot_args[["label_text_x"]],
          cyto_plot_args[["label_text_y"]]
        ))) &
        cnt > 1) {
      label_args <- c("label_text_x", "label_text_y")
      lapply(label_args, function(arg) {
        arg_ind <- which(is.na(cyto_plot_args[[arg]]))
        cyto_plot_args[[arg]][arg_ind] <<- cyto_plot_memory[[arg]][arg_ind]
      })
      # CYTO_PLOT WILL TRANSFORM LABEL CO-ORDINATES WHEN NOT SAVING
      # INVERSE TRANSFORM LABEL_TEXT_X
      cyto_plot_args[["label_text_x"]] <- 
        .cyto_transform(cyto_plot_args[["label_text_x"]],
                        trans = cyto_plot_args[["axes_trans"]],
                        channel = cyto_plot_args[["channels"]][1],
                        inverse = TRUE)
      # INVERSE TRANSFORM LABEL_TEXT_Y
      if(length(channels) == 2){
        cyto_plot_args[["label_text_y"]] <- 
          .cyto_transform(cyto_plot_args[["label_text_y"]],
                          trans = cyto_plot_args[["axes_trans"]],
                          channel = cyto_plot_args[["channels"]][2],
                          inverse = TRUE)
      }
    }
    # CYTO_PLOT
    do.call("cyto_plot", cyto_plot_args)
    # UPDATE CYTO_PLOT_MEMORY
    if (cnt == 1) {
      cyto_plot_memory <<- .cyto_plot_args_recall()
      cyto_plot_memory <<- cyto_plot_memory[[length(cyto_plot_memory)]]
    }
    # RECORD FULL PAGE & OPEN NEW DEVICE
    if(par("page")) {
      if(getOption("cyto_plot_method") == "flowSet") {
        p <- cyto_plot_record()
      } else {
        p <- NULL
      }
      # OPEN NEW DEVICE 
      if(length(args) > cnt) {
        cyto_plot_new(popup) # use global layout
      }
    } else {
      p <- NULL
    }
    # RETURN RECORDED PLOT
    return(p)
  })
  
  # FORMAT CYTO_PLOT_MEMORY ----------------------------------------------------
  
  # FORMAT CYTO_PLOT_MEMORY - INITIAL SAVE
  if (!getOption("cyto_plot_save")) {
    
    # CYTO_PLOT_MEMORY
    cyto_plot_memory <- .cyto_plot_args_recall()
    
    # FORMAT
    if (length(cyto_plot_memory) == length(args)) {
      cyto_plot_memory <- list(cyto_plot_memory)
    } else {
      start <- length(cyto_plot_memory) - length(args)
      cyto_plot_memory <- c(
        cyto_plot_memory[seq_len(start)],
        list(cyto_plot_memory[seq(
          start + 1,
          start + 1 + length(args),
          1
        )])
      )
    }
    
    # UPDATE CYTO_PLOT_MEMORY
    .cyto_plot_args_save(cyto_plot_memory)
  }
  
  # CYTO_PLOT_SAVE -------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save") == TRUE) {
    if (is(x, getOption("cyto_plot_method"))) {
      if (!getOption("cyto_plot_custom")) {
        # CLOSE GRAPHICS DEVICE
        dev.off()
      }
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RESET CYTO_PLOT_LAYOUT
  if(!getOption("cyto_plot_custom")) {
    options("cyto_plot_layout" = NULL)
  }

  # RETURN RECORDED PLOT
  plots[LAPPLY(plots, "is.null")] <- NULL
  if (length(plots) == 0) {
    plots <- NULL
  }
  invisible(plots)
}

#' @rdname cyto_plot
#' @export
cyto_plot.flowFrame <- function(x,
                                channels,
                                axes_trans = NA,
                                overlay = NA,
                                gate = NA,
                                axes_limits = "auto",
                                display = 25000,
                                margins = c(NA, NA, NA, NA),
                                popup = FALSE,
                                xlim = c(NA, NA),
                                ylim = c(NA, NA),
                                xlab,
                                ylab,
                                title,
                                negate = FALSE,
                                hist_stat = "percent",
                                hist_bins = 256,
                                hist_bandwidth = NA,
                                hist_smooth = 1,
                                hist_stack = 0,
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
                                label_text,
                                label_stat,
                                label_position = "auto",
                                label_text_x = NA,
                                label_text_y = NA,
                                label_text_font = 2,
                                label_text_size = 1,
                                label_text_col = "black",
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
                                ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    # SET PLOT METHOD
    options("cyto_plot_method" = "flowFrame")
    # RESET PLOT METHOD ON EXIT
    on.exit(options("cyto_plot_method" = NULL))
  }
  
  # POPUP - CYTO_PLOT_SAVE 
  if (getOption("cyto_plot_save") == TRUE) {
    popup <- FALSE
  }
  
  # CYTO_PLOT OBJECT -----------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list(...)
  
  # THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # CHANNELS/TRANSFORMERS ------------------------------------------------------
  
  # CHANNELS
  if(.empty(args$channels)) {
    stop("Supply channel/marker(s) to construct the plot.")
  } else {
    args$channels <- cyto_channels_extract(x, 
                                           channels, 
                                           plot = TRUE)
  }
  
  # AXES_TRANS
  if(!.all_na(args$axes_trans)){
    if(is(args$axes_trans, "transformList")) {
      message("'axes_trans' must be a transformerList.")
      args$axes_trans <- NA
    }
  }
  
  # CYTO_PLOT_MEMORY -----------------------------------------------------------
  
  # MEMORY
  cyto_plot_memory <- .cyto_plot_args_recall()
  
  # FLOWFRAME METHOD
  if (getOption("cyto_plot_method") == "flowFrame") {
    # MEMORY RESET - CYTO_PLOT_SAVE
    if (getOption("cyto_plot_save") == FALSE) {
      # RESET MEMORY
      .cyto_plot_args_remove()
      cyto_plot_memory <- NULL
      # INHERIT MEMORY - CYTO_PLOT_SAVE
    } else if (getOption("cyto_plot_save") == TRUE) {
      # REPLACE NA WITH MEMORY
      if (any(is.na(c(args[["label_text_x"]], args[["label_text_y"]])))) {
        # NO MEMORY
        if (is.null(cyto_plot_memory)) {
          # RESET CYTO_PLOT_SAVE
          cyto_plot_save_reset()
          stop(
            paste(
              "cyto_plot must be called without cyto_plot_save to",
              "select label co-ordinates, once co-ordinates have been",
              "selected call cyto_plot_save and repeat the cyto_plot call."
            )
          )
        }
        # MEMORY INDEX
        if (is.null(names(cyto_plot_memory))) {
          memory_ind <- 1
        } else {
          memory_ind <- match(NA, names(cyto_plot_memory))
        }
        # UPDATE CYTO_PLOT_MEMORY (USED ARGS)
        names(cyto_plot_memory)[memory_ind] <- "DONE"
        # RESET CYTO_PLOT_MEMORY INDICATOR
        if (!any(is.na(names(cyto_plot_memory)))) {
          names(cyto_plot_memory) <- NULL
        }
        # INHERIT MEMORY
        label_args <- c("label_text_x", "label_text_y")
        lapply(label_args, function(arg) {
          label_arg_length <- length(cyto_plot_memory[[memory_ind]][[arg]])
          label_arg_memory <- cyto_plot_memory[[memory_ind]][[arg]]
          label_arg <- rep(c(args[[arg]], rep(NA, label_arg_length)),
                           length.out = label_arg_length
          )
          arg_ind <- which(is.na(label_arg))
          label_arg[arg_ind] <- label_arg_memory[arg_ind]
          args[[arg]] <<- label_arg
        })
      }
    }
  }
  
  # COMBINE X & OVERLAY - LIST OF CYTOFRAMES -----------------------------------
  
  # X - LIST
  args$x <- structure(list(x), names = cyto_names(x))
  
  # OVERLAY - LIST
  if (!.all_na(args$overlay)) {
    # overlay must be list of CYTOFRAMES
    # flowFrame overlay added to list
    if (is(args$overlay, "flowFrame")) {
      args$overlay <- list(args$overlay)
      # flowSet overlay convert to list of CYTOFRAMES
    } else if (is(args$overlay, "flowSet")) {
      args$overlay <- cyto_list(args$overlay)
      # flowFrame list overlay as is - flowSet list overlay use overlay[[1]]
    } else if (is(args$overlay, "list")) {
      # overlay should be list of CYTOFRAMES
      if (all(LAPPLY(args$overlay, function(z) {
        is(z, "flowFrame")
      }))) {
        args$overlay <- args$overlay
        # overlay list of flowSets - use first fs convert to list of CYTOFRAMES
      } else if (all(LAPPLY(args$overlay, function(z) {
        is(z, "flowSet")
      }))) {
        args$overlay <- cyto_list(args$overlay[[1]])
        # overlay not supported
      } else {
        stop(paste(
          "'overlay' should be either the names of the populations to",
          "overlay, a cytoframe, a cytoset or a list of cytoframes!"
        ))
      }
    }
  }
  
  # MERGE X & OVERLAY - REMOVE OVERLAY ARGUMENT
  if(!.all_na(overlay)) {
    args$x <- c(args$x, args$overlay)
    args <- args[-match("overlay", names(args))]
  }
  
  # SAMPLING & CLEANING --------------------------------------------------------
  
  # SAMPLING - SET SEED
  if (args$display != 1) {
    args$x <- cyto_sample(args$x,
                          display = args$display,
                          seed = 56,
                          plot = TRUE)
    args$display <- 1 # downstream APIs should not re-sample
  }
  
  # REMOVAL NEGATIVE FSC/SSC EVENTS - POINT_COL SCALE ISSUE
  args$x <- cyto_apply(args$x, function(cf){
    lapply(seq_along(args$channels), function(z){
      if(grepl("FSC", args$channels[z], ignore.case = TRUE) |
         grepl("SSC", args$channels[z], ignore.case = TRUE)){
        if (BiocGenerics::nrow(cf) != 0) {
          if (min(range(cf, type = "data")[, args$channels[z]]) < 0) {
            coords <- matrix(c(0, Inf), 
                             nrow = 2,
                             dimnames = list(c("min", "max"),
                                             args$channels[z]))
            cf <<- Subset(cf, rectangleGate(.gate = coords))
          }
        }
      }
    })
    return(cf)
  }, simplify = FALSE, input = "cytoframe")
  
  # HISTOGRAMS -----------------------------------------------------------------
  
  # HISTOGRAMS - STAT/SMOOTH/STACK
  if (length(args$channels) == 1) {
    args$d <- do.call(".cyto_plot_hist", args)
  } else {
    args$d <- NA
  }

  # GATES ----------------------------------------------------------------------
  
  # LIST OF GATE OBJECTS
  if (!.all_na(args$gate)) {
    args$gate <- cyto_gate_prepare(args$gate, 
                                   args$channels)
  }
  
  # LABELS ---------------------------------------------------------------------
  
  # LABEL ARGUMENTS
  label_args <- do.call(".cyto_plot_label_args", args)
  args[c("label", 
         "label_text", 
         "label_stat")] <- label_args[c("label", 
                                        "label_text", 
                                        "label_stat")]
  
  # TRANSFORM SUPPLIED LABEL CO-ORDINATES - NOT DURING SAVING
  if (!getOption("cyto_plot_save")) {
    # TRANSFORM LABEL_TEXT_X - SKIP CYTO_PLOT_SAVE
    args$label_text_x <- .cyto_transform(args$label_text_x,
                                         trans = args$axes_trans,
                                         channel = args$channels[1],
                                         inverse = FALSE)
    
    # TRANSFORM LABEL_TEXT_Y - SKIP CYTO_PLOT_SAVE
    if(length(args$channels) == 2){
      args$label_text_y <- .cyto_transform(args$label_text_y,
                                           trans = args$axes_trans,
                                           channel = args$channels[2],
                                           inverse = FALSE)
    }
  }
  
  # AXES LIMITS ----------------------------------------------------------------
  
  # TRANSFORM AXES LIMITS - FLOWFRAME METHOD ONLY
  if (getOption("cyto_plot_method") == "flowFrame") {
    # X AXES LIMITS
    args$xlim <- .cyto_transform(args$xlim,
                                 trans = args$axes_trans,
                                 channel = args$channels[1],
                                 inverse = FALSE)
    # Y AXES LIMITS
    if(length(args$channels) == 2){
      args$ylim <- .cyto_transform(args$ylim,
                                   trans = args$axes_trans,
                                   channel = args$channels[2],
                                   inverse = FALSE)
    }
  }
  
  # LEGEND ---------------------------------------------------------------------
  
  # LEGEND_TEXT
  if (.all_na(args$legend_text) | 
      length(unique(args$legend_text)) == 1) {
    args$legend_text <- cyto_names(args$x)
  }  
  
  # COMPUTE LABEL STATISTICS ---------------------------------------------------
  
  # STATISTICS
  if (args$label == TRUE) {
    # POPULATIONS
    args$pops <- do.call(".cyto_plot_label_pops", args)
    # STATISTICS
    args$label_stat <- do.call(".cyto_plot_label_stat", args)
    # COMBINE LABEL_TEXT & LABEL_STAT
    args$label_text <- do.call(".cyto_plot_label_text", args)
  }
  
  # ARGUMENT SPLITTING ---------------------------------------------------------
  
  # REPEAT ARGUMENTS
  args <- .cyto_plot_args_split(args)
  
  # PLOT CONSTRUCTION ----------------------------------------------------------
  
  # CYTO_PLOT_BUILD
  args <- .cyto_plot_build(args)
  
  # UPDATE MEMORY --------------------------------------------------------------
  
  # MEMORY UPDATE - CYTO_PLOT_SAVE FALSE
  if (getOption("cyto_plot_save") == FALSE) {
    if (is.null(cyto_plot_memory)) {
      .cyto_plot_args_save(list(list(
        "label_text_x" = args$label_text_x,
        "label_text_y" = args$label_text_y
      )))
    } else {
      .cyto_plot_args_save(c(
        cyto_plot_memory,
        list(list(
          "label_text_x" = args$label_text_x,
          "label_text_y" = args$label_text_y
        ))
      ))
    }
    # MEMORY UPDATE - CYTO_PLOT_SAVE
  } else if (getOption("cyto_plot_save") == TRUE) {
    .cyto_plot_args_save(cyto_plot_memory)
  }
  
  # RECORD PLOT ----------------------------------------------------------------
  
  # RECORD FLOWFRAME METHOD ONLY
  if (getOption("cyto_plot_method") == "flowFrame") {
    p <- cyto_plot_record()
  } else {
    p <- NULL
  }
  
  # CYTO_PLOT_SAVE -------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save") == TRUE) {
    if (is(x, getOption("cyto_plot_method"))) {
      if (!getOption("cyto_plot_custom")) {
        # CLOSE GRAPHICS DEVICE
        dev.off()
      }
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN RECORDED PLOT
  invisible(p)
}
