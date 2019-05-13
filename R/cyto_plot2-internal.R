# Internal Function Definition for cyto_plot -----------------------------------

# ARGUMENTS: These internal functions inherit all arguments from cyto_plot as a
# named list. For a full list of supported arguments refer to ?cyto_plot.
# Argument list received has already inherited cyto_plot_theme arguments. All
# arguments in these internal functions must be extracted from and replaced in
# args (e.g. args[["title"]]). All checks are performed at the upper cyto_plot
# layers and all work is performed in the lower internal layers. Missing
# arguments will be replaced with "".

# ARGUMENT CONVENTIONS: Arguments that can be replaced internally are set to
# missing by default and replaced if not assigned NA. Arguments that can be
# either supplied or not are set to NA by default. Stick top using "channels"
# for arguments accepting channel name(s) this make it easier to pass arguments
# through do.call. 

# MISSING ARGUMENTS:
# - title
# - xlab
# - ylab
# - label
# - label_stat

# CHANNELS: "channels" of length 1 in cyto_plot_1d or length 2 for cyto_plot_2d.
# "channels" have already been converted to valid names in the upper cyto_plot
# layer.

# ADDING NEW FEATURES TO CYTO_PLOT: 
# 1. Add the argument to cyto_plot with an appropriate default. This will 
# automatically be passed with all the other arguments to the internal plotting
# functions in a named list. 
# 2. Modify the code in the called cyto_plot internal function (below) to add 
# the new feature. New arguments can be accessed from the argument list by name 
# (e.g. args[["argument_name"]]).

# DO.CALL does not work on internal functions...

#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @noRd
.cyto_plot_1d <- function(x, ...){
  UseMethod(".cyto_plot_1d")
}

#' @importFrom flowCore exprs parameters identifier
#' @importFrom flowWorkspace pData
#' @importFrom graphics plot axis title abline polygon legend par box
#' @importFrom grDevices adjustcolor
#' @importFrom stats density
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @export
.cyto_plot_1d.flowFrame <- function(x, ...) {
  
  # Get current graphics parameters and reset on exit
  pars <- par("mar")
  on.exit(par(pars))
  
  # Assign arguments to function environment (eliminates args subsetting)
  .args_update(...)
  
  # Restrict x to display percentage events
  x <- cyto_sample(x, display)

  # Convert overlay to list of flowFrames
  if(!.all_na(overlay)){
    if(any(inherits(overlay,"flowFrame") |
           inherits(overlay, "flowSet"))){
      overlay <- cyto_convert(overlay, "list of flowFrames")
    }
  }
  
  # Restrict overlay to display percentage events - bypass in gate_draw
  if(!getOption("CytoRSuite_gate_draw") & display != 1){ 
    overlay <- cyto_sample(overlay, display)
  }
  
  # Combine x and overlay into the same list
  if(!.all_na(overlay)){
    fr_list <- c(list(x), overlay)
  }else{
    fr_list <- list(x)
  }
  
  # Name each element of fr_list with identifier - used in legend
  names(fr_list) <- unlist(lapply(fr_list, function(y){
    identifier(y)
  }))
  
  # Get kernel density for each list element
  fr_dens <- lapply(fr_list, function(x){
    
    suppressWarnings(.cyto_density(x,
                                   channel = channels,
                                   smooth = density_smooth,
                                   modal = density_modal))
    
  })
  
  # Number of overlays
  ovn <- length(fr_dens) - 1
  
  # Calculate the mean maximum y value for kernel densities
  if(density_modal){
    y_max <- 100
  }else{
    y_max <- mean(unlist(lapply(fr_dens, function(d){
      if(!.all_na(d)){
        max(d$y)
      }else{
        NA
      }
    })), na.rm = TRUE)
  }
  
  # Stacked distributions require shifting of y values
  shft <- seq(0,
              ovn * density_stack * y_max,
              density_stack * y_max)
  
  # Adjust y values if stacking is been applied
  if(density_stack > 0){  
    # Shift distributions for stacking
    lapply(seq_len(length(fr_dens)), function(z){
      if(!.all_na(fr_dens[[z]])){
        fr_dens[[z]]$y <<- fr_dens[[z]]$y + shft[z]
      }
    })
  }
  
  # XLIM
  if(.all_na(xlim)){
    xlim <- .cyto_range(fr_list,
                        channels = channels[1],
                        limits = limits,
                        plot = TRUE)[,1]
  }
  
  # YLIM 
  if(.all_na(ylim)){
    ylim <- c(0, y_max + ovn * density_stack * y_max)
  }
  
  # TITLE
  title <- .cyto_plot_title(x,
                            channels = channels,
                            overlay = overlay,
                            title = title)
  
  # AXES LABELS
  labs <- .cyto_plot_axes_label(x,
                                channels = channels,
                                xlab = xlab,
                                ylab = ylab,
                                density_modal = density_modal)
  
  # XLAB
  xlab <- labs[[1]]
  
  # YLAB
  ylab <- labs[[2]]
  
  # POPUP
  if(popup){
    cyto_plot_window()
  }

  # LEGEND TEXT - required for setting plot margins
  if(.all_na(legend_text)){
    legend_text <- names(fr_list)
  }
  
  # MARGINS
  .cyto_plot_margins(x,
                     overlay = overlay,
                     legend = legend,
                     legend_text = legend_text,
                     title = title,
                     axes_text = axes_text)

  # EMPTY PLOT - handles margins and axes limits internally
  args <- .args_list()
  .args <- formalArgs("cyto_plot_empty")
  do.call("cyto_plot_empty", 
          args[names(args) %in% .args])

  # DENSITY FILL - inherits theme internally
  if(.all_na(density_fill)){
    density_fill <- .cyto_plot_density_fill(fr_dens,
                                            density_fill = density_fill,
                                            density_cols = density_cols,
                                            density_fill_alpha = density_fill_alpha)
  }

  # DENSITY 
  cyto_plot_density(fr_dens,
                    density_modal = density_modal,
                    density_stack = density_stack,
                    density_cols = density_cols,
                    density_fill = density_fill,
                    density_fill_alpha = density_fill_alpha,
                    density_line_type = density_line_type,
                    density_line_width = density_line_width,
                    density_line_col = density_line_col)
  
  # LEGEND
  if(legend != FALSE){
    .cyto_plot_legend(channels = channels,
                      legend = legend,
                      legend_text = legend_text,
                      legend_text_font = legend_text_font,
                      legend_text_size = legend_text_size,
                      legend_text_col = legend_text_col,
                      legend_line_col = legend_line_col,
                      legend_box_fill = legend_box_fill,
                      legend_point_col = legend_point_col,
                      density_fill = density_fill,
                      density_fill_alpha = density_fill_alpha,
                      density_line_type = density_line_type,
                      density_line_width = density_line_width,
                      density_line_col = density_line_col,
                      point_shape = point_shape,
                      point_size = point_size,
                      point_col = point_col,
                      point_col_alpha = point_col_alpha)
  }
  
  # GATES & LABEL - without overlay
  if (.all_na(overlay)) {
    
    # GATES
    if(!.all_na(gate)){
      gate <- cyto_plot_gate(gate,
                             channels = channels,
                             gate_line_col = gate_line_col,
                             gate_line_width = gate_line_width,
                             gate_line_type = gate_line_type)
    }
    
    # LABEL?
    if(!.all_na(gate) & 
       .empty(label)){
      label<- TRUE # turn on labels if gate
      
    }else if(.all_na(gate) &
             .empty(label)){
      label <- FALSE # turn off labels if no gate
    }
    
    # STAT
    if(!.all_na(gate) &
       .empty(label_stat)){
      label_stat <- "freq"
    }else if(.all_na(gate) &
             .empty(label_stat)){
      label_stat <- NA
    }
    
    # LABELS 
    if (label == TRUE) {
      
      suppressMessages(cyto_plot_label(
        x = fr_list[[1]],
        channels = channels,
        gate = gate,
        trans = axes_trans,
        text_x = label_box_x,
        text_y = label_box_y,
        text = label_text,
        stat = label_stat,
        text_size = label_text_size,
        text_font = label_text_font,
        text_col = label_text_col,
        box_alpha = label_box_alpha
      ))
      
    }

  # GATE & LABEL - with overlay
  } else if (!.all_na(overlay)) {
    
    # Calculate label y positions for stacked overlays
    if(density_stack != 0){
     
      # Need to compute y label positions
      if(.all_na(label_box_y)){
        label_box_y <- unlist(
          lapply(rep(seq(1,length(fr_list)),
                     length.out = length(gate) * length(fr_list),
                     each = length(gate)),
                function(x){
                (0.5 * density_stack * y_max) +
                ((x-1) * density_stack * y_max)
          }))
        
      }
    
      # LABEL?
      if(!.all_na(gate) & 
         .empty(label)){
        label<- TRUE # turn on labels if gate
      }else if(.all_na(gate) &
               .empty(label)){
        label <- FALSE # turn off labels if no gate
      }
      
      # STAT
      if(!.all_na(gate) &
         .empty(label_stat)){
        label_stat <- "freq"
      }else if(.all_na(gate) &
               .empty(label_stat)){
        label_stat <- NA
      }
      
      # Gate density overlays
      .cyto_plot_overlay_gate(
        x = fr_list[[1]],
        channel = channels,
        trans = axes_trans,
        overlay = fr_list[2:length(fr_list)],
        gate = gate,
        density_stack = density_stack,
        density_modal = density_modal,
        label = label,
        label_text = label_text,
        label_stat = label_stat,
        label_text_size = label_text_size,
        label_text_font = label_text_font,
        label_text_col = label_text_col,
        label_box_x = label_box_x,
        label_box_y = label_box_y,
        label_box_alpha = label_box_alpha,
        gate_line_col = gate_line_col,
        gate_line_width = gate_line_width,
        gate_line_type = gate_line_type
      )
      
    }else if(!.all_na(gate) & density_stack == 0){
      message("Gating overlays without stacking is not supported.")
    }
  }
}
  
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @noRd
cyto_plot_1d.flowSet <- function(x, ...){
  
}

#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @noRd
.cyto_plot_2d <- function(x, ...){
  UseMethod(".cyto_plot_2d")
}

#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @noRd
.cyto_plot_2d.flowFrame <- function(x, ...){
  
  # Get current graphics parameters and reset on exit
  pars <- par("mar")
  on.exit(par(pars))
  
  # Assign arguments to function environment (eliminates args subsetting)
  .args_update(...)
  
  # Restrict x to display percentage events
  x <- cyto_sample(x, display)
  
  # Convert overlay to list
  if(!.all_na(overlay) &
     any(inherits(overlay, "flowFrame") |
         inherits(overlay, "flowSet"))){
    overlay <- cyto_convert(overlay, "list of flowFrames")
  }
  
  # Restrict overlay to display percentage events - bypass in gate_draw
  if(!getOption("CytoRSuite_gate_draw") & display != 1){ 
    overlay <- cyto_sample(overlay, display)
  }
  
  # Combine x and overlay into the same list
  if(!.all_na(overlay)){
    fr_list <- c(list(x), overlay)
  }else{
    fr_list <- list(x)
  }
  
  # Name each element of fr_list with identifier - used in legend
  names(fr_list) <- unlist(lapply(fr_list, function(y){
    identifier(y)
  }))
  
  # Number of overlays
  ovn <- length(fr_list) - 1
  
  # TITLE - always name of flowFrame or "Combined Events"
  title <- .cyto_plot_title(x,
                            channels = channels,
                            overlay = overlay,
                            title = title)
  
  # AXES LABELS
  labs <- .cyto_plot_axes_label(x,
                                channels = channels,
                                xlab = xlab,
                                ylab = ylab,
                                density_modal = density_modal)
  
  # XLAB
  xlab <- labs[[1]]
  
  # YLAB
  ylab <- labs[[2]]
  
  # POPUP
  if(popup){
    cyto_plot_window()
  }
  
  # LEGEND TEXT - required for setting plot margins
  if(.all_na(legend_text)){
    legend_text <- names(fr_list)
  }
  
  # MARGINS
  .cyto_plot_margins(x,
                     overlay = overlay,
                     legend = legend,
                     legend_text = legend_text,
                     title = title,
                     axes_text = axes_text)
  
  # EMPTY PLOT - handles margins and axes limits internally
  args <- .args_list()
  .args <- formalArgs("cyto_plot_empty")
  do.call("cyto_plot_empty", 
          args[names(args) %in% .args])
  
  # POINT COL - list
  if(.all_na(point_col)){
    point_col <- .cyto_plot_point_col(fr_list,
                                      channels = channels,
                                      point_col_scale = point_col_scale,
                                      point_cols = point_cols,
                                      point_col = point_col,
                                      point_col_alpha = point_col_alpha)
  }
  
  # POINTS - list of point colours
  cyto_plot_point(fr_list,
                  channels = channels,
                  point_shape = point_shape,
                  point_size = point_size,
                  point_col_scale = point_col_scale,
                  point_cols = point_cols,
                  point_col = point_col,
                  point_col_alpha = point_col_alpha)
  
  # POINT DENSITY COLOUR SCALE
  point_cols <- .cyto_plot_point_cols(point_cols)
  
  # POINT_COL LEGEND - vector (replace density colours with first point_cols)
  point_col <- unlist(lapply(point_col, function(z){
    
    # colours defined for each point
    if(length(z) > 1){
      return(point_cols[1])
    }else{
      return(z)
    }
    
  }))
  
  # LEGEND
  if(legend != FALSE){
    .cyto_plot_legend(channels = channels,
                      legend = legend,
                      legend_text = legend_text,
                      legend_text_font = legend_text_font,
                      legend_text_size = legend_text_size,
                      legend_text_col = legend_text_col,
                      legend_line_col = legend_line_col,
                      legend_box_fill = legend_box_fill,
                      legend_point_col = legend_point_col,
                      density_fill = density_fill,
                      density_fill_alpha = density_fill_alpha,
                      density_line_type = density_line_type,
                      density_line_width = density_line_width,
                      density_line_col = density_line_col,
                      point_shape = point_shape,
                      point_size = point_size,
                      point_col = point_col,
                      point_col_alpha = point_col_alpha)
  }
  
  # GATES
  if(!.all_na(gate)){
    gate <- cyto_plot_gate(gate,
                           channels = channels,
                           gate_line_type = gate_line_type,
                           gate_line_width = gate_line_width,
                           gate_line_col = gate_line_col)
  }
  
  # LABEL
  if(.empty(label)){
    label <- TRUE
  }
  
  # LABEL STAT
  if(!.all_na(gate) &
     .empty(label_stat)){
    label_stat <- "freq"
  }else if(.all_na(gate) &
           .empty(label_stat)){
    label_stat <- NA
  }
  
  # LABELS
  if(!.all_na(gate) & label == TRUE){
    
    cyto_plot_label(fr_list[[1]],
                    channels = channels,
                    trans = axes_trans,
                    text = label_text,
                    gate = gate,
                    stat = label_stat,
                    text_x = label_box_x,
                    text_y = label_box_y,
                    text_font = label_text_font,
                    text_size = label_text_size,
                    text_col = label_text_col,
                    box_alpha = label_box_alpha)
    
  # LABELS WITHOUT GATES  
  }else if(.all_na(gate) &
           .all_na(label_text) &
           label == TRUE){
    
    # label - limited to # layers - arg_split
    mapply(
      function(label_text,
               label_stat,
               label_text_font,
               label_text_size,
               label_text_col,
               label_box_x,
               label_box_y,
               label_box_alpha) {
        if (label_stat == "percent") {
          label_stat <- NA
        }
        suppressMessages(cyto_plot_label(
          x = fr,
          channels = channels,
          gates = gate,
          trans = axes_trans,
          text = label_text,
          stat = label_stat,
          text_x = label_box_x,
          text_y = label_box_y,
          text_size = label_text_size,
          text_font = label_text_font,
          text_col = label_text_col,
          box_alpha = label_box_alpha
        ))
      }, label_text,
      label_stat,
      label_text_font,
      label_text_size,
      label_text_col,
      label_box_x,
      label_box_y,
      label_box_alpha
    )
    
  }
  
}

#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @noRd
.cyto_plot_2d.flowSet <- function(x, ...){
  
}