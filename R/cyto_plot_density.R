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

  # Combine x and overlay into a list - fr_list
  if (!.all_na(overlay)) {
    fr_list <- c(list(x), .cyto_convert(overlay, "flowFrame list"))
  } else {
    fr_list <- list(x)
  }

  # Pull down arguments to named list
  args <- as.list(environment())

  # Inherit cyto_plot_theme arguments
  args <- .cyto_plot_theme_inherit(args)

  # Number of overlays
  ovn <- length(fr_list) - 1
  
  # Get a list of kernel densities
  fr_dens <- lapply(fr_list, function(fr){
    
    # Calculate kernel density
    suppressWarnings(.cyto_density(fr,
                                   channel,
                                   args[["density_smooth"]],
                                   args[["density_modal"]]))
  })
  
  # Get height of y axis for each density distribution
  if(args[["density_modal"]]){
    y_max <- 100
  }else{
    y_max <- mean(unlist(lapply(fr_dens, function(z){
    
      if(!.all_na(z)){
        max(z$y)
      }else{
        NA
      }
    
    })), na.rm = TRUE)
    
  }
  
  # Get vector of density_stack values
  ofst <- seq(
    0,
    ovn * args[["density_stack"]] * y_max,
    args[["density_stack"]] * y_max
  )

  # Get a list of kernel densities
  fr_dens <- mapply(function(fr_dens, ofst) {

    # Adjust values for stacking
    if (ofst != 0 & !.all_na(fr_dens)) {
      fr_dens$y <- fr_dens$y + ofst
    }

    return(fr_dens)
    
  }, fr_dens, ofst, SIMPLIFY = FALSE)

  # Get density_fill
  .args <- formalArgs(".cyto_plot_density")
  args[["density_fill"]] <- do.call(".cyto_plot_density", 
                                    c(fr_dens, args[.args]))

  # Add horizontal lines
  abline(
    h = ofst,
    col = args[["density_line_col"]],
    lwd = args[["density_line_width"]],
    lty = args[["density_line_type"]]
  )
  
  # Add density distributions - reverse plot order and colours
  if (!.all_na(args[["overlay"]]) & 
      args[["density_stack"]] == 0) {
    
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
      args[["density_fill"]],
      args[["density_line_col"]],
      args[["density_line_width"]],
      args[["density_line_type"]])
    
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
      rev(args[["density_fill"]]),
      rev(args[["density_line_col"]]),
      rev(args[["density_line_width"]]),
      rev(args[["density_line_type"]]))
    
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
  
  # Check x contains density objects (may contain NA as well)
  if(!any(unlist(lapply(x, "class")) == "density")){
    stop("'x' should be a list of density objects.")
  }
  
  # Pull down arguments to named list
  args <- args_list()
  
  # Inherit theme arguments
  args <- .cyto_plot_theme_inherit(args)
  
  # Get density_fill colours (inherits from theme internally)
  if(.all_na(density_fill) | .empty(density_fill)){
    .args <- formalArgs(".cyto_plot_density_fill")
    args[["density_fill"]] <- do.call(".cyto_plot_density_fill",
                                      args[.args])
  }

  # Number of overlays
  ovn <- length(x) - 1
  
  # Calculate the mean maximum y value for kernel densities
  if(args[["density_modal"]]){
    y_range <- 100
  }else{
    y_range<- mean(unlist(lapply(x, function(d){
      if(!.all_na(d)){
        max(d$y) - min(d$y)
      }else{
        NA
      }
    })), na.rm = TRUE)
  }
  ofst <- seq(0,
              ovn * args[["density_stack"]] * y_range,
              args[["density_stack"]] * y_range)
  abline(
    h = ofst,
    col = args[["density_line_col"]],
    lwd = args[["density_line_width"]],
    lty = args[["density_line_type"]],
    xpd = FALSE
  )
  
  # Minimum for each distribution
  mn <- unlist(lapply(x, function(z){
    if(!.all_na(z)){
      min(z$y)
    }else{
      0
    }
  }))
  
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
      args[["density_fill"]],
      args[["density_line_col"]],
      args[["density_line_width"]],
      args[["density_line_type"]])
    
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
      rev(args[["density_fill"]]),
      rev(args[["density_line_col"]]),
      rev(args[["density_line_width"]]),
      rev(args[["density_line_type"]]))
    
  }
  
}