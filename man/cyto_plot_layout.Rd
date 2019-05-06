% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyto_plot-helpers.R
\name{cyto_plot_layout}
\alias{cyto_plot_layout}
\title{Set Panel Layout for cyto_plot}
\usage{
cyto_plot_layout(nrow, ncol, by = "row")
}
\arguments{
\item{nrow}{number of rows to include in the panel layout.}

\item{ncol}{number of columns to include in the panel layout.}

\item{by}{chracter string indicating whether the plot should be filled by
\code{"row"} or \code{"column"}, set to \code{"row"} by default.}
}
\description{
\code{cyto_plot_layout()} sets the panel layout dimensions for combining
different types of cyto_plot plots. Make a call to \code{cyto_plot_layout()}
prior to making multiple calls to \code{cyto_plot()}.
}
\examples{
library(CytoRSuiteData)

# Load samples into GatingSet
fs <- Activation
gs <- GatingSet(fs)

# Apply compensation
gs <- compensate(gs, fs[[1]]@description$SPILL)

# Transform fluorescent channels
trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
gs <- transform(gs, trans)

# Apply gatingTemplate
gt <- Activation_gatingTemplate
gating(gt, gs)

# Set out plot layout
cyto_plot_layout(1, 2)

# Add 2D plot
cyto_plot(gs[[4]],
  parent = "CD4 T Cells",
  alias = "",
  channels = c("Alexa Fluor 647-A", "7-AAD-A"),
  layout = FALSE
)

# Add 1D plot
cyto_plot(gs,
  parent = "CD4 T Cells",
  alias = "",
  channels = "7-AAD-A",
  density_stack = 0.6,
  layout = FALSE
)
}
\author{
Dillon Hammill (Dillon.Hammill@anu.edu.au)
}