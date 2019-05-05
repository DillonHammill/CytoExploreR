# Internal Function Definition for cyto_plot -----------------------------------

# ARGUMENTS: 
# These internal functions inherit all arguments from cyto_plot as a named list.
# For a full list of supported arguments refer to ?cyto_plot.

# ARGUMENT CONVENTIONS:
# Arguments that can be replaced internally are set to missing by default and
# replaced if not assigned NA. Arguments that can be either supplied or not are
# set to NA by default.

# ADDING NEW FEATURES TO CYTO_PLOT: 
# 1. Add the argument to cyto_plot with an appropriate default. This will 
# automatically be passed with all the other arguments to the internal plotting
# functions in a named list. 
# 2. Modify the code in the called cyto_plot internal function (below) to add 
# the new feature. New arguments can be accessed form the argument list by name 
# (e.g. args[["argument_name"]]).

#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' @export
.cyto_plot_1d_v2 <- function(x, ...){
  UseMethod(".cyto_plot_1d_v2")
}

#' @importFrom flowCore exprs parameters identifier
#' @importFrom flowWorkspace pData
#' @importFrom graphics plot axis title abline polygon legend par box
#' @importFrom grDevices adjustcolor
#' @importFrom stats density
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @export
.cyto_plot_1d_v2.flowFrame <- function(x, ...) {
  
  # Plot attributes commented with capitalisation
  
  # Get current graphics paraneters and reset on exit
  pars <- par("mar")
  on.exit(par(pars))
  
  # Pull down arguments to named list
  args <- ...
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Return channel if is marker supplied
  args[["channel"]] <- cyto_channels_extract(args[["x"]], 
                                             args[["channel"]], 
                                             TRUE)
  
  # Restrict x to display percentage events
  args[["x"]] <- cyto_sample(args[["x"]], args[["display"]])
  
  # Convert overlay to list
  if(!.all_na(args[["overlay"]])){
    args[["overlay"]] <- .cyto_convert(args[["overlay"]], "flowFrame list")
  }
  
  # Restrict overlay to display percentage events - bypass in gate_draw
  if(!getOption("CytoRSuite_gate_draw") & args[["display"]] != 1){ 
    args[["overlay"]] <- lapply(args[["overlay"]], function(x){
      cyto_sample(x, args[["display"]])
    })
  }
  
  # Number of overlays
  ovn <- length(args[["overlay"]])
  
  # Combine x and overlay into the same list
  if(!.all_na(args[["overlay"]])){
    args[["fr_list"]] <- c(list(args[["x"]]), args[["overlay"]])
  }else{
    args[["fr_list"]] <- list(args[["x"]])
  }
  
  # Name each element of fr_list with identifier - used in legend
  names(args[["fr_list"]]) <- unlist(lapply(args[["fr_list"]], function(x){
    identifier(x)
  }))
  
  # Get kernel density for each list element
  args[["fr_dens"]] <- lapply(args[["fr_list"]], function(x){
    
    .cyto_density(args[["x"]],
                  args[["channel"]],
                  args[["density_smooth"]],
                  args[["density_modal"]])
    
  })
  
  # Calculate the mean maximum y value for kernel densities
  if(args[["density_modal"]]){
    y_max <- 100
  }else{
    y_max <- mean(unlist(lapply(args[["fr_dens"]], function(d){
      max(d$y) # max(NA) returns NA
    })), na.rm = TRUE)
  }
  
  # Stacked distributions require shifting of y values
  shft <- seq(0,
              ovn * args[["density_stack"]] * y_max,
              args[["density_stack"]] * y_max)
  
  # Shift distributions for stacking
  lapply(seq_len(length(args[["fr_dens"]])), function(z){
    args[["fr_dens"]][[z]]$y <<- args[["fr_dens"]][[z]]$y + shft[z]
  })
  
  # YLIM - bypass.cyto_plot_axes_limits
  if(.all_na(args[["ylim"]])){
    args[["ylim"]] <- c(0, y_max + ovn * args[["density_stack"]] * y_max)
  }
  
  # TITLE
  args[["title"]] <- .cyto_plot_title(args[["x"]],
                            args[["channel"]],
                            args[["overlay"]],
                            args[["title"]])
  
  # AXES LABELS
  labs <- .cyto_plot_axes_label(args[["x"]],
                                args[["channel"]],
                                args[["xlab"]],
                                args[["ylab"]],
                                args[["density_modal"]])
  
  # XLAB
  args[["xlab"]] <- labs[[1]]
  
  # YLAB
  args[["ylab"]] <- labs[[2]]
  
  # POPUP
  if(popup){
    cyto_plot_window()
  }
   
  # LEGEND TEXT - required for setting plot margins
  if(.empty(args[["legend_text"]])){
    args[["legend_text"]] <- names(args[["fr_list"]])
  }
  
  # MARGINS
  .cyto_plot_margins(args[["x"]],
                     args[["overlay"]],
                     args[["legend"]],
                     args[["legend_text"]],
                     args[["title"]],
                     args[["axes_text"]])
  
  print(par("mar"))
  
  # EMPTY PLOT - handles margins and axes limits internally
  .args <- formalArgs("cyto_plot_empty")
  do.call("cyto_plot_empty", 
          c(args["channel"], args[names(args) %in% .args]))

  # DENSITY FILL - inherits theme internally
  if(.empty(args[["density_fill"]])){
    .args <- formalArgs(".cyto_plot_density_fill")
    args[["density_fill"]] <- do.call(".cyto_plot_density_fill",
                                      c(list("x" = args[["fr_dens"]]),
                                        args[names(args) %in% .args][-1]))
  }
  
  print(args[["density_fill"]])
  
  # DENSITY 
  cyto_plot_density(args[["fr_dens"]],
                    args[["density_cols"]],
                    args[["density_fill"]],
                    args[["density_fill_alpha"]],
                    args[["density_line_type"]],
                    args[["density_line_width"]],
                    args[["density_line_col"]])
  
  # GATES
  
  
  # LABELS
  
  
  # LEGEND TEXT
  
  # LEGEND
  if(args[["legend"]] != FALSE){
    .args <- formalArgs(".cyto_plot_legend")
    .args[1] <- "channel"
    do.call(".cyto_plot_legend",
            args[names(args) %in% .args])
  }
  
}
  







