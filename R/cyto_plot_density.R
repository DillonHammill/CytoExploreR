# CYTO_PLOT_DENSITY ------------------------------------------------------------

#' Add density distributions to empty plots
#'
#' @param x object of class flowFrame or a list of density objects.
#' @param channel name of the channels to be used to construct the plot.
#' @param overlay list of flowFRame objects to overlay.
#' @param density_modal logical indicating whether density should be normalised
#'   to mode and presented as a percentage. Set to \code{TRUE} by default.
#' @param density_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust kernel density.
#' @param density_stack numeric [0,1] indicating the degree of offset for
#'   overlaid populations, set to 0.5 by default.
#' @param density_cols vector of colls to use to select density_fill colours
#' @param density_fill colour(s) used to fill polygons.
#' @param density_fill_alpha numeric [0,1] used to control fill transparency,
#'   set to 1 by default to remove transparency.
#' @param density_line_type line type(s) to use for border(s), set to solid
#'   lines by default.
#' @param density_line_width line width for border.
#' @param density_line_col colour(s) for border line, set to "black" by default.
#'
#' @importFrom methods formalArgs
#' @importFrom graphics abline polygon
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @rdname cyto_plot_density
#'
#' @export
cyto_plot_density <- function(x, ...) {
  UseMethod("cyto_plot_density")
}

#' @rdname cyto_plot_density
#' @export
cyto_plot_density.flowFrame <- function(x,
                                        channel,
                                        overlay = NA,
                                        density_modal = TRUE,
                                        density_smooth = 1.5,
                                        density_stack = 0.5,
                                        density_cols = NA,
                                        density_fill = NA,
                                        density_fill_alpha = 1,
                                        density_line_type = 1,
                                        density_line_width = 1,
                                        density_line_col = "black") {

  # PREPARE DATA ---------------------------------------------------------------
  
  # LIST OF FLOWFRAMES
  if (!.all_na(overlay)) {
    fr_list <- c(list(x), cyto_convert(overlay, "flowFrame list"))
  } else {
    fr_list <- list(x)
  }

  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()

  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)

  # UPDATE ARGUMENTS
  .args_update(args)
  
  # DENSITY --------------------------------------------------------------------
  
  # OVERLAYS
  ovn <- length(fr_list) - 1
  
  # KERNEL DENSITY
  fr_dens <- lapply(fr_list, function(fr){
    suppressWarnings(.cyto_density(fr,
                                   channel,
                                   density_smooth,
                                   density_modal))
  })
  
  # fr_dens contains some density objects
  if(!.all_na(fr_dens)){
    # Get height of y axis for each density distribution
    if(density_modal){
      y_max <- 100
    }else{
      y_max <- mean(LAPPLY(fr_dens, function(z){
    
        if(!.all_na(z)){
          max(z$y)
        }else{
          NA
        }
    
      }), na.rm = TRUE)
    
    }
  
    # Get vector of density_stack values
    ofst <- seq(
      0,
      ovn * density_stack * y_max,
      density_stack * y_max
    )

    # Get a list of offset kernel densities
    fr_dens <- mapply(function(fr_dens, ofst) {

      # Adjust values for stacking
      if (ofst != 0 & !.all_na(fr_dens)) {
        fr_dens$y <- fr_dens$y + ofst
      }

      return(fr_dens)
    
    }, fr_dens, ofst, SIMPLIFY = FALSE)
    
  # fr_dens does not contain any density objects
  } else if(.all_na(fr_dens)){
    
    # Set y_max to 100
    y_max <- 100
    
    # Get positions for horizontal lines
    ofst <- seq(
      0,
      ovn * density_stack * y_max,
      density_stack * y_max
    )
    
  }

  # DENSITY_FILL
  density_fill <- .cyto_plot_density_fill(fr_dens,
                                          density_fill = density_fill,
                                          density_cols = density_cols,
                                          density_fill_alpha = density_fill_alpha)

  # Add horizontal lines
  abline(
    h = ofst,
    col = density_line_col,
    lwd = density_line_width,
    lty = density_line_type
  )
  
  # SPLIT ARGUMENTS ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args[c("density_fill",
                                       "density_fill_alpha",
                                       "density_line_type",
                                       "density_line_width",
                                       "density_line_col")],
                                channels = channels,
                                n = length(fr_list),
                                plots = 1,
                                layers = length(fr_list),
                                gates = 0)
  
  # UNLIST ARGUMENTS - SINGLE PLOT
  args <- lapply(args, function(z){z[[1]]})
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # Add density distributions - reverse plot order and colours
  if (!.all_na(overlay) & 
      density_stack == 0) {
    
    mapply(
      function(fr_dens,
                     density_fill,
                     density_line_col,
                     density_line_width,
                     density_line_type) {
        if(!.all_na(fr_dens)){
          polygon(fr_dens,
            col = density_fill,
            border = density_line_col,
            lwd = density_line_width,
            lty = density_line_type)
        }
      }, frs_dens,
      density_fill,
      density_line_col,
      density_line_width,
      density_line_type)
    
  } else {
    
    mapply(
      function(fr_dens,
                     density_fill,
                     density_line_col,
                     density_line_width,
                     density_line_type) {
        if(!.all_na(fr_dens)){
          polygon(fr_dens,
                  col = density_fill,
                  border = density_line_col,
                  lwd = density_line_width,
                  lty = density_line_type)
        }

      }, rev(frs_dens),
      rev(density_fill),
      rev(density_line_col),
      rev(density_line_width),
      rev(density_line_type))
    
  }
}

#' @rdname cyto_plot_density
#' @export
cyto_plot_density.list <- function(x,
                                   density_modal = TRUE,
                                   density_stack = 0.5,
                                   density_cols = NA,
                                   density_fill = NA,
                                   density_fill_alpha = 1,
                                   density_line_type = 1,
                                   density_line_width = 1,
                                   density_line_col = "black", ...){
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # Get density_fill colours (inherits from theme internally)
  density_fill <- .cyto_plot_density_fill(x,
                                          density_fill = density_fill,
                                          density_cols = density_cols,
                                          density_fill_alpha = density_fill_alpha)

  # Number of overlays
  ovn <- length(x) - 1
  
  # Calculate the mean maximum y value for kernel densities
  if(density_modal){
    y_max <- 100
  }else{
    y_max<- mean(LAPPLY(x, function(d){
      if(!.all_na(d)){
        max(d$y) - min(d$y)
      }else{
        NA
      }
    }), na.rm = TRUE)
  }
  
  # Empty list - all NA without density objects
  if(is.nan(y_max)){
    y_max <- 100
  }
  
  ofst <- seq(0,
              ovn * density_stack * y_max,
              density_stack * y_max)
  
  abline(
    h = ofst,
    col = density_line_col,
    lwd = density_line_width,
    lty = density_line_type,
    xpd = FALSE
  )
  
  # Minimum for each distribution
  mn <- LAPPLY(x, function(z){
    if(!.all_na(z)){
      min(z$y)
    }else{
      0
    }
  })
  
  # SPLIT ARGUMENTS ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args,
                                channels = channels,
                                n = length(x),
                                plots = 1,
                                layers = length(x),
                                gates = 0)
  
  # UNLIST ARGUMENTS
  args <- lapply(args, function(z){z[[1]]})
  
  # PLOT DENSITY ---------------------------------------------------------------
  
  # Add density distributions - reverse plot order and colours
  if (length(x) > 1 & 
      all(floor(mn) == 0)) {
    
    mapply(
      function(x,
               density_fill,
               density_line_col,
               density_line_width,
               density_line_type) {
        if(!.all_na(x)){
          polygon(x,
                  col = density_fill,
                  border = density_line_col,
                  lwd = density_line_width,
                  lty = density_line_type)
        }
      }, x,
      density_fill,
      density_line_col,
      density_line_width,
      density_line_type)
    
  } else {
    
    mapply(
      function(x,
               density_fill,
               density_line_col,
               density_line_width,
               density_line_type) {
        if(!.all_na(x)){
          polygon(x,
                  col = density_fill,
                  border = density_line_col,
                  lwd = density_line_width,
                  lty = density_line_type)
        }
      }, rev(x),
      rev(density_fill),
      rev(density_line_col),
      rev(density_line_width),
      rev(density_line_type))
    
  }
  
}