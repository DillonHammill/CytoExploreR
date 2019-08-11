# CYTO_PLOT_EXPLORE ------------------------------------------------------------

#' Explore Cytometry Data in Bivariate Plots
#'
#' \code{cyto_plot_explore} is an extremely useful tool to explore all aspectes
#' of your cytometry data. \code{cyto_plot_explore} constructs a faceted plot
#' for each of the supplied \code{channels_x} by plotting them against all the
#' supplied \code{channels_y}. Thus providing a rapid means tro explore the data
#' in all available channels.
#'
#' @param x an object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierachy} or \code{GatingSet}. \code{flowSet} and
#'   \code{GatingSet} objects will be coerced to a single \code{flowFrame} prior
#'   to plotting.
#' @param parent name of the parent population to plot when a
#'   \code{GatingHierarchy} or \code{GatingSet} object is supplied, set to the
#'   last gated population by default.
#' @param channels_x vector of channels to explore, set to all fluorescent
#'   channels by default. Each channel in \code{channels_x} will be plotted
#'   against all \code{channels_y}, with each \code{channels_x} on a separate
#'   page.
#' @param channels_y vector of channels to plot each channels_x against, set to
#'   all channels by default.
#' @param axes_trans object of class transformerList passed to \code{cyto_plot}
#'   to appropriately display transformed data in plots.
#' @param layout vector of grid dimensions \code{c(#rows,#columns)} for each
#'   plot.
#' @param popup logical indicating whether plots should be constructed in a
#'   pop-up window.
#' @param title text to include above each plot, set to NA by default to remove
#'   titles.
#' @param header title to use for the plots, set to the name of the sample by
#'   default. Turn off the header by setting this argument to NA.
#' @param header_text_font font to use for header text, set to 2 by default.
#' @param header_text_size text size for header, set to 1 by default.
#' @param header_text_col colour for header text, set to "black" by default.
#' @param ... additional arguments passed to \code{cyto_plot}.
#'
#' @importFrom grDevices n2mfrow
#' @importFrom graphics mtext
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @name cyto_plot_explore
NULL

#' @noRd
#' @export
cyto_plot_explore <- function(x, ...){
  UseMethod("cyto_plot_explore")
}

#' @rdname cyto_plot_explore
#' @export
cyto_plot_explore.GatingSet <- function(x,
                                        parent = NULL,
                                        channels_x = NULL,
                                        channels_y = NULL,
                                        axes_trans = NA,
                                        layout,
                                        popup = FALSE,
                                        title = NA,
                                        header,
                                        header_text_font = 2,
                                        header_text_size = 1,
                                        header_text_col = "black",
                                        ...){
  
  # Set plot method
  if (is.null(getOption("CytoRSuite_cyto_plot_method"))) {
    options("CytoRSuite_cyto_plot_method" = "Explore/GatingSet")
  }
  
  # Default parent population
  if(is.null(parent)){
    parent <- cyto_nodes(x)[length(cyto_nodes(x))]
  }
  
  # Extract parent population
  fs <- cyto_extract(x)
  
  # Convert fs to flowFrame
  fr <- cyto_convert(fs)
  
  # Call to flowFrame method
  cyto_plot_explore(x = fr,
                    channels_x = channels_x,
                    channels_y = channels_y,
                    axes_trans = axes_trans,
                    layout = layout,
                    popup = popup,
                    title = title,
                    header = header,
                    header_text_font = header_text_font,
                    header_text_size = header_text_size,
                    header_text_col = header_text_col,
                    ...)
  
  # Turn off graphics device for saving
  if (getOption("CytoRSuite_cyto_plot_save")) {
    if (inherits(x, basename(getOption("CytoRSuite_cyto_plot_method")))) {
      
      # Close graphics device
      dev.off()
      
      # Reset CytoRSuite_cyto_plot_save
      options("CytoRSuite_cyto_plot_save" = FALSE)
      
      # Reset CytoRSuite_cyto_plot_method
      options("cytoRSuite_cyto_plot_method" = NULL)
    }
  }
  
}

#' @rdname cyto_plot_explore
#' @export
cyto_plot_explore.GatingHierarchy <- function(x,
                                              parent = NULL,
                                              channels_x = NULL,
                                              channels_y = NULL,
                                              axes_trans = NA,
                                              layout,
                                              popup = FALSE,
                                              title = NA,
                                              header,
                                              header_text_font = 2,
                                              header_text_size = 1,
                                              header_text_col = "black",
                                              ...){
  
  # Set plot method
  if (is.null(getOption("CytoRSuite_cyto_plot_method"))) {
    options("CytoRSuite_cyto_plot_method" = "Explore/GatingHierarchy")
  }
  
  # Default parent population
  if(is.null(parent)){
    parent <- cyto_nodes(x)[length(cyto_nodes(x))]
  }
  
  # Extract parent population
  fr <- cyto_extract(x)
  
  # Call to flowFrame method
  cyto_plot_explore(x = fr,
                    channels_x = channels_x,
                    channels_y = channels_y,
                    axes_trans = axes_trans,
                    layout = layout,
                    popup = popup,
                    title = title,
                    header = header,
                    header_text_font = header_text_font,
                    header_text_size = header_text_size,
                    header_text_col = header_text_col,
                    ...)
  
  # Turn off graphics device for saving
  if (getOption("CytoRSuite_cyto_plot_save")) {
    if (inherits(x, basename(getOption("CytoRSuite_cyto_plot_method")))) {
      
      # Close graphics device
      dev.off()
      
      # Reset CytoRSuite_cyto_plot_save
      options("CytoRSuite_cyto_plot_save" = FALSE)
      
      # Reset CytoRSuite_cyto_plot_method
      options("cytoRSuite_cyto_plot_method" = NULL)
    }
  }
  
}

#' @rdname cyto_plot_explore
#' @export
cyto_plot_explore.flowSet <- function(x,
                                      channels_x = NULL,
                                      channels_y = NULL,
                                      axes_trans = NA,
                                      layout,
                                      popup = FALSE,
                                      title = NA,
                                      header,
                                      header_text_font = 2,
                                      header_text_size = 1,
                                      header_text_col = "black",
                                      ...){
  
  # Set plot method
  if (is.null(getOption("CytoRSuite_cyto_plot_method"))) {
    options("CytoRSuite_cyto_plot_method" = "Explore/flowSet")
  }
  
  # Merge flowSet to flowFrame
  fr <- cyto_convert(x, "flowFrame")
  
  # Call to flowFrame method
  cyto_plot_explore(x = fr,
                    channels_x = channels_x,
                    channels_y = channels_y,
                    axes_trans = axes_trans,
                    layout = layout,
                    popup = popup,
                    title = title,
                    header = header,
                    header_text_font = header_text_font,
                    header_text_size = header_text_size,
                    header_text_col = header_text_col,
                    ...)
  
  # Turn off graphics device for saving
  if (getOption("CytoRSuite_cyto_plot_save")) {
    if (inherits(x, basename(getOption("CytoRSuite_cyto_plot_method")))) {
      
      # Close graphics device
      dev.off()
      
      # Reset CytoRSuite_cyto_plot_save
      options("CytoRSuite_cyto_plot_save" = FALSE)
      
      # Reset CytoRSuite_cyto_plot_method
      options("cytoRSuite_cyto_plot_method" = NULL)
    }
  }
  
}

#' @rdname cyto_plot_explore
#' @export
cyto_plot_explore.flowFrame <- function(x,
                                        channels_x = NULL,
                                        channels_y = NULL,
                                        axes_trans = NA,
                                        layout,
                                        popup = FALSE,
                                        title = NA,
                                        header,
                                        header_text_font = 2,
                                        header_text_size = 1,
                                        header_text_col = "black",
                                        ...){
  
  # Set plot method
  if (is.null(getOption("CytoRSuite_cyto_plot_method"))) {
    options("CytoRSuite_cyto_plot_method" = "Explore/flowFrame")
  }
  
  # Graphics parameters
  pars <- par(c("mfrow","oma"))
  on.exit(par(pars))
  
  # Extract channels
  chans <- colnames(x)
  fluor_chans <- cyto_fluor_channels(x)
  
  # Use all fluorescent channels as default for channels_x
  if(is.null(channels_x)){
    channels_x <- fluor_chans
  }
  
  # Use all channels as default for channels_y
  if(is.null(channels_y)){
    channels_y <- chans
  }
  
  # Pop up
  if(popup == TRUE){
    cyto_plot_new(popup)
  }
  
  # Layout
  if(missing(layout)){
    layout <- c(
      n2mfrow(length(channels_y))[2],
      n2mfrow(length(channels_y))[1]
    )
    par(mfrow = layout)
  }else{
    if(layout[1] == FALSE){
      # Do nothing
    }else{
      par(mfrow = layout)
    }
  }
  
  # Header
  if(missing(header)){
    header <- cyto_names(x)
  }
  
  # Space for header
  if(!.all_na(header)){
    par(oma = c(0, 0, 3, 0))
  }
  
  # Plots - loop through channels_x
  lapply(seq_len(length(channels_x)), function(z){
    lapply(seq_len(length(channels_y)), function(y){
      # Construct plot for each channels_y
      if(channels_x[z] == channels_y[y]){
        cyto_plot(x,
                  channels = channels_x[z],
                  axes_trans = axes_trans,
                  title = title,
                  legend = FALSE, ...)
      }else{
        cyto_plot(x,
                  channels = c(channels_x[z], channels_y[y]),
                  axes_trans = axes_trans,
                  title = title,
                  legend = FALSE, ...)
      }

      # Add header & call new plot
      if(channels_y[y] == channels_y[length(channels_y)]){
        # Add header
        if(!.all_na(header)){
          mtext(header,
                outer = TRUE,
                font = header_text_font,
                cex = header_text_size,
                col = header_text_col)
        }
        # Call new plot
        if(popup == TRUE){
          cyto_plot_new(popup)
          par(mfrow = layout)
          par(oma = c(0, 0, 3, 0))
        }else{
          plot.new()
          par(mfrow = layout)
          par(oma = c(0, 0, 3, 0))
        }
      }
      
      # Custom layouts require headers
      if(channels_y[y] == prod(layout)){
        if (!.all_na(header)) {
          mtext(header,
                outer = TRUE,
                cex = header_text_size,
                font = header_text_font,
                col = header_text_col
          )
        }
      }
      
    })
  })
  
  # Turn off graphics device for saving
  if (getOption("CytoRSuite_cyto_plot_save")) {
    if (inherits(x, basename(getOption("CytoRSuite_cyto_plot_method")))) {
      
      # Close graphics device
      dev.off()
      
      # Reset CytoRSuite_cyto_plot_save
      options("CytoRSuite_cyto_plot_save" = FALSE)
      
      # Reset CytoRSuite_cyto_plot_method
      options("cytoRSuite_cyto_plot_method" = NULL)
    }
  }
  
}
