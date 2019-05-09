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
  
  # Pull down arguments to named list - ... includes x as well
  assign("args", ...)
  
  # Restrict x to display percentage events
  args[["x"]] <- cyto_sample(args[["x"]], args[["display"]])

  # Convert overlay to list
  if(!.all_na(args[["overlay"]])){
    args[["overlay"]] <- cyto_convert(args[["overlay"]], "flowFrame list")
  }
  
  # Restrict overlay to display percentage events - bypass in gate_draw
  if(!getOption("CytoRSuite_gate_draw") & args[["display"]] != 1){ 
    args[["overlay"]] <- lapply(args[["overlay"]], function(x){
      cyto_sample(x, args[["display"]])
    })
  }
  
  # Combine x and overlay into the same list
  if(!.all_na(args[["overlay"]])){
    args[["fr_list"]] <- c(list(args[["x"]]), args[["overlay"]])
  }else{
    args[["fr_list"]] <- list(args[["x"]])
  }
  
  # Name each element of fr_list with identifier - used in legend
  names(args[["fr_list"]]) <- unlist(lapply(args[["fr_list"]], function(y){
    identifier(y)
  }))
  
  # Get kernel density for each list element
  args[["fr_dens"]] <- lapply(args[["fr_list"]], function(x){
    
    suppressWarnings(.cyto_density(x,
                                   args[["channels"]],
                                   args[["density_smooth"]],
                                   args[["density_modal"]]))
    
  })
  
  # Number of overlays
  ovn <- length(args[["fr_dens"]]) - 1
  
  # Calculate the mean maximum y value for kernel densities
  if(args[["density_modal"]]){
    y_max <- 100
  }else{
    y_max <- mean(unlist(lapply(args[["fr_dens"]], function(d){
      if(!.all_na(d)){
        max(d$y)
      }else{
        NA
      }
    })), na.rm = TRUE)
  }
  
  # Stacked distributions require shifting of y values
  shft <- seq(0,
              ovn * args[["density_stack"]] * y_max,
              args[["density_stack"]] * y_max)
  
  # Adjust y values if stacking is been applied
  if(args[["density_stack"]] > 0){  
    # Shift distributions for stacking
    lapply(seq_len(length(args[["fr_dens"]])), function(z){
      if(!.all_na(args[["fr_dens"]][[z]])){
        args[["fr_dens"]][[z]]$y <<- args[["fr_dens"]][[z]]$y + shft[z]
      }
    })
  }
  
  # YLIM 
  if(.all_na(args[["ylim"]])){
    args[["ylim"]] <- c(0, y_max + ovn * args[["density_stack"]] * y_max)
  }

  # TITLE
  args[["title"]] <- .cyto_plot_title(args[["x"]],
                            args[["channels"]],
                            args[["overlay"]],
                            args[["title"]])
  
  # AXES LABELS
  labs <- .cyto_plot_axes_label(args[["x"]],
                                args[["channels"]],
                                args[["xlab"]],
                                args[["ylab"]],
                                args[["density_modal"]])
  
  # XLAB
  args[["xlab"]] <- labs[[1]]
  
  # YLAB
  args[["ylab"]] <- labs[[2]]
  
  # POPUP
  if(args[["popup"]]){
    cyto_plot_window()
  }
   
  # LEGEND TEXT - required for setting plot margins
  if(.all_na(args[["legend_text"]])){
    args[["legend_text"]] <- names(args[["fr_list"]])
  }
  
  # MARGINS
  .cyto_plot_margins(args[["x"]],
                     args[["overlay"]],
                     args[["legend"]],
                     args[["legend_text"]],
                     args[["title"]],
                     args[["axes_text"]])
  
  # EMPTY PLOT - handles margins and axes limits internally
  .args <- formalArgs("cyto_plot_empty")
  do.call("cyto_plot_empty", 
          args[names(args) %in% .args])

  # DENSITY FILL - inherits theme internally
  if(.all_na(args[["density_fill"]])){
    args[["density_fill"]] <- .cyto_plot_density_fill(args[["fr_dens"]],
                                                      args[["density_fill"]],
                                                      args[["density_cols"]],
                                                      args[["density_fill_alpha"]])
  }
  
  # DENSITY 
  cyto_plot_density(args[["fr_dens"]],
                    args[["density_modal"]],
                    args[["density_stack"]],
                    args[["density_cols"]],
                    args[["density_fill"]],
                    args[["density_fill_alpha"]],
                    args[["density_line_type"]],
                    args[["density_line_width"]],
                    args[["density_line_col"]])

  # LEGEND
  if(args[["legend"]] != FALSE){
    .cyto_plot_legend(args[["channels"]],
                      args[["legend"]],
                      args[["legend_text"]],
                      args[["legend_text_font"]],
                      args[["legend_text_size"]],
                      args[["legend_text_col"]],
                      args[["legend_line_col"]],
                      args[["legend_box_fill"]],
                      args[["legend_point_col"]],
                      args[["density_fill"]],
                      args[["density_fill_alpha"]],
                      args[["density_line_type"]],
                      args[["density_line_width"]],
                      args[["density_line_col"]],
                      args[["point_shape"]],
                      args[["point_size"]],
                      args[["point_col"]],
                      args[["point_alpha"]])
  }

  # GATES - no overlay
  if (.all_na(args[["overlay"]])) {
    if(!.all_na(args[["gate"]])){
      args[["gate"]] <- cyto_plot_gate(args[["gate"]],
                             channels = args[["channels"]],
                             gate_line_col = args[["gate_line_col"]],
                             gate_line_width = args[["gate_line_width"]],
                             gate_line_type = args[["gate_line_type"]])
    }
    
    # LABELS
    if (!.all_na(args[["gate"]]) & args[["label"]]) {
      
      # Population names missing - show percentage only
      suppressMessages(cyto_plot_label(
        x = args[["fr_list"]][[1]],
        channels = args[["channels"]],
        gate = args[["gate"]],
        trans = args[["axes_trans"]],
        text = args[["label_text"]],
        stat = args[["label_stat"]],
        text_size = args[["label_text_size"]],
        text_font = args[["label_text_font"]],
        text_col = args[["label_text_col"]],
        box_alpha = args[["label_box_alpha"]]
      ))
    }

  } else if (!.all_na(args[["overlay"]]) & 
             args[["density_stack"]] != 0 & 
             !.all_na(args[["gate"]])) {
    
    # Need to compute y label positions
    if(.all_na(args[["label_box_y"]])){
      args[["label_box_y"]] <- unlist(
        lapply(rep(seq(1,length(args[["fr_list"]])),
                   length.out = length(args[["gate"]]) * length(args[["fr_list"]]),
                   each = length(args[["gate"]])),
              function(x){
              (0.5 * args[["density_stack"]] * y_max) +
              ((x-1) * args[["density_stack"]] * y_max)
        }))
    }
    
    .cyto_overlay_gate(
      x = args[["fr_list"]][[1]],
      channel = args[["channels"]],
      trans = args[["axes_trans"]],
      overlay = args[["fr_list"]][2:length(args[["fr_list"]])],
      gate = args[["gate"]],
      density_stack = args[["density_stack"]],
      density_modal = args[["density_modal"]],
      label_text = args[["label_text"]],
      label_stat = args[["label_stat"]],
      label_text_size = args[["label_text_size"]],
      label_text_font = args[["label_text_font"]],
      label_text_col = args[["label_text_col"]],
      label_box_x = args[["label_box_x"]],
      label_box_y = args[["label_box_y"]],
      label_box_alpha = args[["label_box_alpha"]],
      gate_line_col = args[["gate_line_col"]],
      gate_line_width = args[["gate_line_width"]],
      gate_line_type = args[["gate_line_type"]]
    )
  } else if (!.all_na(args[["overlay"]]) & 
             args[["density_stack"]] == 0 & 
             !.all_na(args[["gate"]])) {
    message("Gating overlays without stacking is not supported.")
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
  
  # Pull down arguments to named list - ... includes x as well
  assign("args", ...)
  
  # Restrict x to display percentage events
  args[["x"]] <- cyto_sample(args[["x"]], args[["display"]])
  
  # Convert overlay to list
  if(!.all_na(args[["overlay"]])){
    args[["overlay"]] <- cyto_convert(args[["overlay"]], "flowFrame list")
  }
  
  # Restrict overlay to display percentage events - bypass in gate_draw
  if(!getOption("CytoRSuite_gate_draw") & args[["display"]] != 1){ 
    args[["overlay"]] <- lapply(args[["overlay"]], function(x){
      cyto_sample(x, args[["display"]])
    })
  }
  
  # Combine x and overlay into the same list
  if(!.all_na(args[["overlay"]])){
    args[["fr_list"]] <- c(list(args[["x"]]), args[["overlay"]])
  }else{
    args[["fr_list"]] <- list(args[["x"]])
  }
  
  # Name each element of fr_list with identifier - used in legend
  names(args[["fr_list"]]) <- unlist(lapply(args[["fr_list"]], function(y){
    identifier(y)
  }))
  
  # Number of overlays
  ovn <- length(args[["fr_list"]]) - 1
  
  # TITLE - always name of flowFrame or "Combined Events"
  args[["title"]] <- .cyto_plot_title(args[["x"]],
                                      args[["channels"]],
                                      args[["overlay"]],
                                      args[["title"]])
  
  # AXES LABELS
  labs <- .cyto_plot_axes_label(args[["x"]],
                                args[["channels"]],
                                args[["xlab"]],
                                args[["ylab"]],
                                args[["density_modal"]])
  
  # XLAB
  args[["xlab"]] <- labs[[1]]
  
  # YLAB
  args[["ylab"]] <- labs[[2]]
  
  # POPUP
  if(args[["popup"]]){
    cyto_plot_window()
  }
  
  # LEGEND TEXT - required for setting plot margins
  if(.all_na(args[["legend_text"]])){
    args[["legend_text"]] <- names(args[["fr_list"]])
  }
  
  # MARGINS
  .cyto_plot_margins(args[["x"]],
                     args[["overlay"]],
                     args[["legend"]],
                     args[["legend_text"]],
                     args[["title"]],
                     args[["axes_text"]])
  
  # EMPTY PLOT - handles margins and axes limits internally
  .args <- formalArgs("cyto_plot_empty")
  do.call("cyto_plot_empty", 
          args[names(args) %in% .args])
  
  # POINT COL - list
  if(.all_na(args[["point_col"]])){
    args[["point_col"]] <- .cyto_plot_point_col(args[["fr_list"]],
                                                args[["channels"]],
                                                args[["point_col_scale"]],
                                                args[["point_cols"]],
                                                args[["point_col"]],
                                                args[["point_col_alpha"]])
  }
  
  # POINTS - list of point colours
  cyto_plot_point(args[["fr_list"]],
                  args[["channels"]],
                  args[["point_shape"]],
                  args[["point_size"]],
                  args[["point_col_scale"]],
                  args[["point_cols"]],
                  args[["point_col"]],
                  args[["point_col_alpha"]])
  
  # POINT DENSITY COLOUR SCALE
  args[["point_cols"]] <- .cyto_plot_point_cols(args[["point_cols"]])
  
  # POINT_COL LEGEND - vector (replace density colours with first point_cols)
  args[["point_col"]] <- unlist(lapply(args[["point_col"]], function(z){
    
    # colours defined for each point
    if(length(z) > 1){
      return(args[["point_cols"]][1])
    }else{
      return(z)
    }
    
  }))
  
  # LEGEND
  if(args[["legend"]] != FALSE){
    .cyto_plot_legend(args[["channels"]],
                      args[["legend"]],
                      args[["legend_text"]],
                      args[["legend_text_font"]],
                      args[["legend_text_size"]],
                      args[["legend_text_col"]],
                      args[["legend_line_col"]],
                      args[["legend_box_fill"]],
                      args[["legend_point_col"]],
                      args[["density_fill"]],
                      args[["density_fill_alpha"]],
                      args[["density_line_type"]],
                      args[["density_line_width"]],
                      args[["density_line_col"]],
                      args[["point_shape"]],
                      args[["point_size"]],
                      args[["point_col"]],
                      args[["point_col_alpha"]])
  }
  
  # GATES
  
  # LABELS
  
}

#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @noRd
.cyto_plot_2d.flowSet <- function(x, ...){
  
}