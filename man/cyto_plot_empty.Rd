% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyto_plot-helpers.R
\name{cyto_plot_empty}
\alias{cyto_plot_empty}
\alias{cyto_plot_empty.flowFrame}
\title{Create an empty cyto_plot}
\usage{
cyto_plot_empty(x, ...)

\method{cyto_plot_empty}{flowFrame}(
  x,
  channels,
  axes_trans = NA,
  overlay = NA,
  gate = NA,
  xlim = NA,
  ylim = NA,
  axes_limits = "auto",
  axes_limits_buffer = 0.03,
  title,
  xlab,
  ylab,
  margins = NULL,
  density_modal = TRUE,
  density_smooth = 1.5,
  density_stack = 0.5,
  density_cols = NA,
  density_fill = NA,
  density_fill_alpha = 1,
  density_line_type = 1,
  density_line_width = 1,
  density_line_col = "black",
  point_shape = ".",
  point_size = 2,
  point_col_scale = NA,
  point_cols = NA,
  point_col = NA,
  point_col_alpha = 1,
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
  border_line_type = 1,
  border_line_width = 1,
  border_line_col = "black",
  border_fill = "white",
  border_fill_alpha = 1,
  legend = FALSE,
  legend_text,
  legend_text_font = 1,
  legend_text_size = 1,
  legend_text_col = "black",
  legend_line_type = NA,
  legend_line_width = NA,
  legend_line_col = NA,
  legend_box_fill = NA,
  legend_point_col = NA,
  ...
)
}
\arguments{
\item{x}{object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.}

\item{...}{not in use.}

\item{channels}{name of the channel(s) or marker(s) to be used to construct
the plot. The length of channels determines the type of plot to be
constructed, either a 1-D density distribution for a single channel or a
2-D scatterplot with blue-red colour scale for two channels.}

\item{axes_trans}{object of class
\code{\link[flowCore:transformList-class]{transformList}} or
\code{\link[flowWorkspace]{transformerList}} generated by
\code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
transform the fluorescent channels of the supplied flowFrame. This
transformation object will be used internally to ensure that the axes
labels of the plot are appropriately transformed. The transformation object
will NOT be applied to the flowFrame internally and should be applied to
the flowFrame prior to plotting.}

\item{overlay}{a list of flowFrames to overlay onto the plot.}

\item{gate}{list of gate objects to be plotted, used internlaly to ensure
gate co-ordinates are taken into account when computing axes limits.}

\item{xlim}{lower and upper limits of x axis (e.g. c(0,5)).}

\item{ylim}{lower and upper limits of y axis (e.g. c(0,5)).}

\item{axes_limits}{options include \code{"auto"}, \code{"data"} or
\code{"machine"} to use optimised, data or machine limits respectively. Set
to \code{"auto"} by default to use optimised axes ranges. Fine control over
axes limits can be obtained by altering the \code{xlim} and \code{ylim}
arguments.}

\item{axes_limits_buffer}{decimal indicating the percentage of buffering to
add to either end of the axes limits, set to 0.03 by default.}

\item{title}{title to use for the plot, set to the name of the sample by
default. Title can be removed by setting this argument to \code{NA}.}

\item{xlab}{x axis label.}

\item{ylab}{y axis label.}

\item{margins}{a vector of length 4 to control the margins around the bottom,
left, top and right of the plot, set to NULL by default to let `cyto_plot`
compute optimal margins.}

\item{density_modal}{logical indicating whether density should be normalised
to mode and presented as a percentage. Set to \code{TRUE} by default.}

\item{density_smooth}{smoothing parameter passed to
\code{\link[stats:density]{density}} to adjust kernel density.}

\item{density_stack}{numeric [0,1] indicating the degree of offset for
overlaid populations, set to 0.5 by default. #' @param density_cols vector
colours to draw from when selecting density fill colours if none are
supplied to density_fill.}

\item{density_cols}{vector colours to draw from when selecting density fill
colours if none are supplied to density_fill.}

\item{density_fill}{colour(s) used to fill polygons.}

\item{density_fill_alpha}{numeric [0,1] used to control fill transparency,
set to 1 by default to remove transparency.}

\item{density_line_type}{line type(s) to use for border(s), set to solid
lines by default.}

\item{density_line_width}{line width for border.}

\item{density_line_col}{colour(s) for border line, set to "black" by default.}

\item{point_shape}{shape(s) to use for points in 2-D scatterplots, set to
\code{"."} by default to maximise plotting speed.  See
\code{\link[graphics:par]{pch}} for alternatives.}

\item{point_size}{numeric to control the size of points in 2-D scatter plots
set to 2 by default.}

\item{point_col_scale}{vector of colours to use for density gradient.}

\item{point_cols}{vector colours to draw from when selecting colours for
points if none are supplied to point_col.}

\item{point_col}{colour(s) to use for points in 2-D scatter plots, set to NA
by default to use a blue-red density colour scale.}

\item{point_col_alpha}{numeric [0,1] to control point colour transparency in
2-D scatter plots, set to 1 by default to use solid colours.}

\item{axes_text}{logical vector of length 2 indicating whether axis text
should be included for the x and y axes respectively, set to
\code{c(TRUE,TRUE)} by default to display axes text on both axes.}

\item{axes_text_font}{numeric indicating the font to use for axes, set to 1
for plain font by default. See \code{\link[graphics:par]{?par}} font for
details.}

\item{axes_text_size}{character expansion for axis text.}

\item{axes_text_col}{colour of axis text.}

\item{axes_label_text_font}{numeric indicating the font to use for title, set
to 1 for plain font by default. See \code{\link[graphics:par]{?par}} font
for details.}

\item{axes_label_text_size}{character expansion for axis labels.}

\item{axes_label_text_col}{colour of axis labels.}

\item{title_text_font}{numeric indicating the font to use for title, set to 2
for bold font by default. See \code{\link[graphics:par]{?par}} font for
details.}

\item{title_text_size}{character expansion for plot title.}

\item{title_text_col}{colour for plot title.}

\item{border_line_type}{line type to use for plot border, set to 1 by default
for a sold border.}

\item{border_line_width}{line width for plot border, set to 1 by default.}

\item{border_line_col}{line colour for plot border, set to "black" by
default.}

\item{border_fill}{colour to use for the plot background, set to "white" by
default.}

\item{border_fill_alpha}{transparency to use for border_fill colour, set to 1
by default to add no transparency.}

\item{legend}{can be either \code{"line"} or \code{"fill"} to indicate
whether a legend should be constructed based on the density \code{"line"}
or \code{"fill"}, set to FALSE by default to remove the legend.}

\item{legend_text}{vector of labels to use in the legend.}

\item{legend_text_font}{numeric to control the font of legend text, set to 1
for plain font by default. See \code{\link[graphics:par]{font}} for
alternatives.}

\item{legend_text_size}{numeric to control the size of text in the legend,
set to 1 by default.}

\item{legend_text_col}{colour(s) to use for text in legend, set to
\code{"black"} by default.}

\item{legend_line_type}{numeric to control the line type for line legends,
set to 1 by default. Refer to \code{lty} in \code{\link[graphics:par]{par}}
for alternatives.}

\item{legend_line_width}{numeric to control the line width in line legend,
set to 1 by default. Refer to \code{lwd} in \code{\link[graphics:par]{par}}
for alternatives.}

\item{legend_line_col}{colour(s) to use for the lines in 1-D plot legends
when legend is set to \code{"line"}.}

\item{legend_box_fill}{fill colour(s) to use for the boxes in 1-D plot
legends when legend is set to \code{"fill"}.}

\item{legend_point_col}{colour(s) to use for points in 2-D scatter plot
legend.}
}
\description{
\code{cyto_plot_empty} generates to base for cyto_plot by creating an empty
plot with border, axes, axes_text and titles. Data is subsequently added to
this base layer with \code{cyto_plot_point} or \code{cyto_plot_density}.
}
\examples{
library(CytoExploreRData)

# Construct an empty 2D plot with black background
cyto_plot_empty(Activation[[32]],
  channels = c("FSC-A", "SSC-A"),
  border_fill = "black"
)

# Construct an empty 1D plot
cyto_plot_empty(Activation[[32]],
  channels = c("FSC-A"),
  overlay = Activation[1:2]
)
}
\author{
Dillon Hammill (Dillon.Hammill@anu.edu.au)
}
