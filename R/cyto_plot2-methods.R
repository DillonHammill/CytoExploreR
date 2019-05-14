#' @importFrom purrr transpose pmap
#' @export
cyto_plot.GatingSet <- function(){
  
  # Checks ---------------------------------------------------------------------
  
  # Set plot method - determine when dev.off should be called when saving
  if(is.null(getOption("CytoRSuite_cyto_plot_method"))){
    options("CytoRSuite_cyto_plot_method" = "GatingSet")
  }
  
  # Signal a custom plot if layout is FALSE
  if(!is.null(layout)){
    if(all(layout == FALSE)){
      options("CytoRSuite_cyto_plot_custom" = TRUE)
    }
  }
  
  # No parent supplied
  if(missing(parent)){
    stop("Please supply the name of the parent population to plot.")
  }
  
  # Check if channels are supplied
  if (missing(channels)) {
    stop("Supply channel/marker(s) to construct the plot.")
  }
  
  # Graphics parameters --------------------------------------------------------
  
  # Current graphics parameters
  old_pars <- par(c("mar","mfrow"))
  
  # Reset graphics parameters on exit
  on.exit(par(old_pars))
  
  # Extract channels & transformations -----------------------------------------

  # Get valid channel names if markers are supplied
  channels <- cyto_channels_extract(x, channels, plot = TRUE)
  
  # Get transformations
  if (.all_na(axes_trans)) {
    
    # Some transforms found - replace these entries in transList
    trnsfrms <- lapply(channels, function(channel) {
      getTransformations(x[[1]], channel, only.function = FALSE)
    })
    names(trnsfrms) <- channels
    
    # Remove NULL transforms
    trnsfrms[unlist(lapply(trnsfrms, "is.null"))] <- NULL
    
    if (length(trnsfrms) == 0) {
      axes_trans <- NULL
    } else {
      axes_trans <- transformerList(names(trnsfrms), trnsfrms)
    }
    
    axes_trans <- cyto_transform_convert(axes_trans, inverse = FALSE)
  } else {
    axes_trans <- cyto_transform_convert(axes_trans, inverse = FALSE)
  }
  
  # Extract & format data for plotting (list of flowFrame lists) ---------------
  
  # Extract data from GatingSet
  fs <- cyto_extract(x, parent)
  
  # Add data to list and group if necessary - convert groups to flowFrames
  if(!missing(group_by)){
    fs_list <- cyto_group_by(fs, group_by)
    fr_list <- lapply(fr_list, function(z){cyto_convert(z, "flowFrame")})
  }else{
    fr_list <- cyto_convert(fs, "list of flowFrames")
    names(fr_list) <- sampleNames(x)
  }
  
  # Add overlay to list and group if necessary
  if(!.all_na(overlay)){
    
    # overlay must be list of flowFrame lists
    
    # Overlay is the names of populations
    if(is.character(overlay) & all(overlay %in% basename(getNodes(x)))){
      
      # Pull out list of flowSets
      nms <- overlay
      overlay <- lapply(overlay, function(z){
        cyto_extract(x, z)
      })
      names(overlay) <- nms
    }
    
    # Repeat flowFrame once per plot - no grouping
    if(inherits(overlay, "flowFrame")){
      overlay_list <- rep(list(list(overlay)), length(fr_list))
    # Convert to list of flowFrame lists
    }else if(inherits(overlay, "flowSet")){
      
      # Apply same grouping to overlay
      if(!missing(group_by)){
        overlay_list <- cyto_group_by(overlay, group_by)
        overlay_list <- lapply(overlay_list, function(z){
          cyto_convert(z, "flowFrame")
        })
      }else{
        overlay_list <- cyto_convert(overlay, "list of flowFrames")
      }
      
      overlay_list <- lapply(overlay_list, function(z){
          list(z)
      })
    # Convert list of flowSets to list of flowFrame lists
    }else if(inherits(overlay, "list")){
      
      # Only lists of flowSets are supported
      if(!all(unlist(lapply(overlay, function(z){inherits(z, "flowSet")})))){
        stop("Only lists of flowSets are supported for overlay.")
      }
      
      # Apply same grouping to overlay
      if(!missing(group_by)){
        overlay_list <- lapply(overlay, function(z){
          cyto_group_by(z, group_by)})
        overlay_list <- lapply(overlay_list, function(z){
          lapply(z, function(y){
            cyto_convert(y, "flowFrame")
          })
        })
      }else{
        overlay_list <- lapply(overlay, function(z){
          cyto_convert(z, "list of flowFrames")
        })
      }
      
      overlay_list <- overlay_list %>% transpose()
    }
  }
  
  # Combine base layers with overlay into list of flowFrame lists
  nms <- names(fr_list)
  fr_list <- lapply(seq_len(length(fr_list)), function(z){
    c(fr_list[z], overlay_list[[z]])
  })
  names(fr_list) <- nms
  
  # Extract gate objects -------------------------------------------------------
  
  # Allow alias = "" to plot all appropriate gates
  if (!missing(alias)) {
    # Plot all appropriate gates if alias is an empty character string
    if (all(alias == "")) {
      gt <- templateGen(x[[1]])
      gt <- gt[basename(gt$parent) == parent, ]
      alias <- gt$alias[gt$dims == paste(channels, collapse = ",")]
      
      # No gates constructed in the supplied channels
      if (length(alias) == 0) {
        alias <- NA
      }
    }
  }  
  
  # Extract gate objects directly from x
  if(!.all_na(alias)){
    if(missing(group_by)){
      gate <- lapply(nms, function(nm){
        gt <- lapply(alias, function(z){
          getGate(x[[match(nm, pd$group_by)]], z)
        })
        names(gt) <- alias
        return(gt)
      })
    }else{
      gate <- lapply(seq_len(length(x)), function(z) {
        gt <- lapply(alias, function(y) {
          getGate(x[[z]], y)
        })
        names(gt) <- alias
        return(gt)
      })
    }
  }
  
  # Prepare arguments for plotting ---------------------------------------------
  
  # Labels - use alias  if no text supplied
  if (.all_na(label_text)){
    if(!.all_na(alias)) {
      label_text <- alias
    }
  }
  
  # Legend text
  if(missing(legend_text)){
    if(!.all_na(overlay)){
      if(class(overlay) == "character"){
        legend_text <- c(parent, overlay)
      }
    }else{
      legend_text <- paste0("layer", seq_len(length(fr_list[[1]])))
    }
  }
  
  # Title
  if(missing(title)){
    
    # Extract plot title from names of fr_list
    title <- names(fr_list)
    
    # Add parent name to each plot
    title <- unlist(lapply(title, function(z){
      
      # Parent name
      if(parent == "root"){
        pt == "All Events"
      }else{
        pt <- parent
      }
      
      # All samples grouped together
      if(z == "all"){
        z <- "Combine Events"
      }
      
      # Paste together sampleName and parent name
      paste(z, "\n", pt, sep = " ")
      
    }))
    
  }
  
  # X axis breaks and labels - pass through axes_text argument
  if(axes_text[1] == TRUE){
    axes_text_x <- .cyto_plot_axes_text(x,
                                        channels = channels,
                                        axes_trans = axes_trans)[[1]]
  }else{
    axes_text_x <- FALSE
  }
  
  # Y axis breaks and labels - pass through axes_text argument
  if(axes_text[2] == TRUE){
    if(length(channels) == 2){
      axes_text_y <- .cyto_plot_axes_text(x,
                                          channels = channels[2],
                                          axes_trans = axes_trans)[[1]]
    }else{
      axes_text_y <- NA
    }
    
  }else{
    axes_text_y <- FALSE
  }
  axes_text <- list(axes_text_x, axes_text_y)
  
  # X axis limits
  if(.all_na(xlim)){
    xlim <- lapply(fr_list, function(z){
      .cyto_plot_axes_limits(z,
                             channels = channels[1],
                             limits = limits,
                             plot = TRUE)[,1]
    })
    xlim <- c(min(unlist(xlim)), max(unlist(xlim)))
  }
  
  # Y axis limits - 1D y limit calculated downstream
  if(.all_na(ylim) & length(channels) != 1){
    ylim <- lapply(fr_list, function(z){
      .cyto_plot_axes_limits(z,
                             channels = channels[2],
                             limits = limits,
                             plot = TRUE)[,1]
    })
    ylim <- c(min(unlist(ylim)), max(unlist(ylim)))
  }
  
  # Pop-up window
  if(popup == TRUE){
    cyto_plot_window()
  }
  
  # Layout - missing/off/supplied
  if(missing(layout)){
    
    # Layout dimensions
    layout <- .cyto_plot_layout(fr_list,
                                layout = layout,
                                density_stack = density_stack,
                                density_layers = density_layers)
    
  }else if(all(layout == FALSE) | .all_na(layout)){
    
    # Use current dimensions
    layout <- par("mfrow")
    
  }
  par("mfrow" = layout)
  
  # Number of plots per page
  np <- par("mfrow")[1] * par("mfrow")[2]
  
  # Repeat arguments as required -----------------------------------------------
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Pass relevant arguments to .cyto_plot_args_split
  .args <- formalArgs(".cyto_plot_args_split")
  
  # 1D density distributions
  if(length(channels) == 1){
    
    # No overlay
    if(.all_na(overlay)){
      
      # No stacking -  separate panels
      if(all(density_stack == 0)){
        args <- .cyto_plot_args_split(args[names(args) %in% .args],
                                      channels = channels,
                                      n = length(fr_list[[1]]) * length(fr_list),
                                      plots = length(fr_list), 
                                      layers = 1,
                                      gates = length(gate[[1]]))
        
      # Stacking - multiple layers per panel
      }else{
        
        args <- .cyto_plot_args_split(args[names(args) %in% .args],
                                      channels = channels,
                                      n = length(fr_list[[1]]) * 
                                        length(fr_list),
                                      plots = ceiling(length(fr_list) / 
                                                        density_layers), 
                                      layers = density_layers,
                                      gates = length(gate[[1]]))
        
      }
 
    # Overlay  
    }else if(!.all_na(overlay)){
      
      args <- .cyto_plot_args_split(args[names(args) %in% .args],
                                    channels = channels,
                                    n = length(fr_list[[1]]) * length(fr_list),
                                    plots = length(fr_list), 
                                    layers = length(fr_list[[1]]),
                                    gates = length(gate[[1]]))
      
    }
    
  # 2D scatter plots
  }else if (length(channels) == 2) {
    
    args <- .cyto_plot_args_split(args[names(args) %in% .args],
                                  channels = channels,
                                  n = length(fr_list[[1]]) * length(fr_list),
                                  plots = length(fr_list), 
                                  layers = length(fr_list[[1]]),
                                  gates = length(gate[[1]]))
    
  }
  
  # Update arguments
  .args_update(args)
  
  # Calls to cyto_plot internal methods ----------------------------------------
  
  # 1D density distributions
  if(length(channels) == 1){
    
  
  # 2D scatter plots
  }else if(length(channels) == 2){
    
  }
  
  # Record and/or save ---------------------------------------------------------
  
  # Turn off graphics device for saving
  if (getOption("CytoRSuite_cyto_plot_save")) {
    if (inherits(x, getOption("CytoRSuite_cyto_plot_method"))) {
      if (!getOption("CytoRSuite_cyto_plot_custom")) {
        
        # Close graphics device
        dev.off()
      }
      
      # Reset CytoRSuite_cyto_plot_save
      options("CytoRSuite_cyto_plot_save" = FALSE)
      
      # Reset CytoRSuite_cyto_plot_method
      options("cytoRSuite_cyto_plot_method" = NULL)
      
      # Reset custom plotting flag
      options("CytoRSuite_cyto_plot_custom")
    }
  }
  
  # Record plot for assignment
  invisible(recordPlot())
  
}