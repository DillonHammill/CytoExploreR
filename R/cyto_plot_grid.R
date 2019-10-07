#' Grid Layout for cyto_plot
#'
#' \code{cyto_plot_grid} is wrapper for cyto_plot which provides a grid layout
#' without excess white space. \code{cyto_plot_grid} is particularly useful for
#' visualising data grouped using \code{group_by}.
#'
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' @param ... additional arguments passed to \code{cyto_plot}.
#'
#' @importFrom flowCore parameters flowFrame
#' @importFrom flowWorkspace pData
#' @importFrom graphics par mtext text
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @name cyto_plot_grid
NULL

#' @noRd
#' @export
cyto_plot_grid <- function(x, ...){
  UseMethod("cyto_plot_grid")
}

#' @rdname cyto_plot_grid
#' @export
cyto_plot_grid.flowSet <- function(x,
                                   channels = NULL,
                                   axes_trans = NULL,
                                   group_by = FALSE,
                                   overlay = NULL,
                                   gate = NA,
                                   display = 1,
                                   layout = NULL,
                                   popup = FALSE,
                                   xlab,
                                   ylab,
                                   title = c(NA,NA),
                                   format = "row",
                                   header = NULL,
                                   panel_x_label = NULL,
                                   panel_y_label = NULL,
                                   axes_label_text_font = 2,
                                   axes_label_text_size = 1,
                                   axes_label_text_col = "black",
                                   title_text_font = c(2,2),
                                   title_text_size = c(1,1),
                                   title_text_col = c("black","black"),
                                   header_text_font = 2,
                                   header_text_size = 1.5,
                                   header_text_col = "black",
                                   panel_label_font = c(1,1),
                                   panel_label_size = c(1,1),
                                   panel_label_col = c("black","black"), ...) {
  
  # cyto_plot_grid is designed to plot multiple groups only (1-3 variables)
  if(all(group_by[1] %in% c(TRUE,"all"))){
    stop("Use cyto_plot when grouping all samples into the same group.")
  }else if(length(group_by) > 3){
    stop("cyto_plot_grid supports a maximum of 3 grouping variables.")
  }
  
  # Signal cyto_plot_grid is being used
  options("cyto_plot_grid" = TRUE)
  
  # Set plot method
  if(is.null(getOption("cyto_plot_method"))){
    options("cyto_plot_method" = "grid/flowSet")
  }
  
  # Channels missing
  if(missing(channels)){
    stop("Supply channel/marker(s) to construct the plot")
  }
  
  # Assign x to fs
  fs <- x
  smp <- length(fs)
  
  # Extract pData
  pd <- cyto_details(fs)
  
  # Check channels
  channels <- cyto_channels_extract(fs, 
                                    channels, 
                                    plot = TRUE)
  
  # No grouping
  if(all(group_by == FALSE)){
    
    # Set layout
    layout <- .cyto_plot_layout(x = fs,
                                layout = layout,
                                density_stack = density_stack,
                                density_layers = density_layers)
    
    # Grouping  
  }else{
    
    # Set layout
    if(length(group_by) == 1){
      
      if(is.null(layout)){
        
        # Layout in single row
        if(format == "row"){
          layout <- c(1, nlevels(as.factor(pd[,group_by])))
          #Layout in a single column
        }else if(format == "column"){
          layout <- c(nlevels(as.factor(pd[,group_by])), 1)
        }
        
      }
      
    }else{
      
      # Layout uses first 2 variables - 3rd indicates # of sheets
      if(is.null(layout)){
        layout <- c(nlevels(as.factor(pd[,group_by[2]])),
                    nlevels(as.factor(pd[,group_by[1]])))
      }
      
    }
    
  }
  
  # Grid dimensions
  nrow <- layout[1]
  ncol <- layout[2]
  
  # Number of plots per sheet
  np <- prod(layout)
  
  # Total number of plots (np*ns)
  if(all(group_by == FALSE)){
    tnp <- length(fs)
  }else{
    tnp <- prod(unlist(
      lapply(group_by, function(x){
        nlevels(as.factor(pd[,x]))
      })
    ))
  }
  
  # Number of sheets
  ns <- ceiling(tnp/np)

  # Each level of group_by[3] is on a new sheet
  if(length(group_by) == 3){
    ns <- ns * nlevels(as.factor(pd[,group_by[3]]))
  }
  
  # Use layout to determine which plots require axes
  axes_text <- rep(list(c(FALSE,FALSE)), np * ns)
  
  # Grid matrix per sheet
  grd <- matrix(seq_len(ncol * nrow), ncol = ncol, nrow = nrow, byrow = TRUE)
  grd <- lapply(seq_len(ns), function(x){
    grd * ns
  })
  
  # Plots requiring Y axis
  axes_y <- LAPPLY(grd, function(x){
    x[,1]
  })
  lapply(axes_y, function(x){
    axes_text[[x]][2] <<- TRUE
  })
  
  # Plots requiring X axis
  axes_x <- LAPPLY(grd, function(x){
    x[nrow(x),]
  })
  lapply(axes_x, function(y){
    axes_text[[y]][1] <<- TRUE
  })
  
  # Extract channel info
  fr.data <- pData(parameters(fs[[1]]))
  fr.channels <- BiocGenerics::colnames(fs[[1]])
  
  # X axis title
  if (missing(xlab)) {
    if (!is.na(fr.data$desc[which(fr.channels == channels[1])])) {
      xlab <- paste(fr.data$desc[which(fr.channels == channels[1])],
                    channels[1],
                    sep = " "
      )
    } else if (is.na(fr.data$desc[which(fr.channels == channels[1])])) {
      xlab <- paste(channels[1], sep = " ")
    }
  }
  
  # Y axis title
  if (missing(ylab)) {
    if (!is.na(fr.data$desc[which(fr.channels == channels[2])])) {
      ylab <- paste(fr.data$desc[which(fr.channels == channels[2])],
                    channels[2],
                    sep = " "
      )
    } else if (is.na(fr.data$desc[which(fr.channels == channels[2])])) {
      ylab <- paste(channels[2], sep = " ")
    }
  }
  
  # Pop-up window
  cyto_plot_new(popup = popup)
  
  # Set plot layout
  par(mfrow = layout)
  
  # Remove marginal space around each plot
  par(mar = c(0,0,0,0))
  
  # Add space for labels, titles and headers
  if(is.null(header)){
    par(oma = c(6,6,5,5))
  }else{
    par(oma = c(6,6,8,5))
  }
  
  # Create a list containing all groups (without missing levels)
  if(all(group_by == FALSE)){
    grp_exp <- rep(list(NA), tnp)
  }else if(length(group_by) == 1){
    grp_exp <- rep(list(NA), nlevels(as.factor(pd[,group_by])))
    names(grp_exp) <- levels(as.factor(pd[,group_by]))
  }else if(length(group_by) == 2){
    grp_exp <- rep(list(NA), nlevels(as.factor(pd[,group_by[1]])) *
                     nlevels(as.factor(pd[,group_by[2]])))
    names(grp_exp) <- LAPPLY(
      levels(as.factor(pd[,group_by[2]])),
      function(x){
        paste(levels(as.factor(pd[,group_by[1]])), x)
      })
  }else if(length(group_by) == 3){
    grp_exp <- rep(list(NA), nlevels(as.factor(pd[,group_by[1]])) *
                     nlevels(as.factor(pd[,group_by[2]])) *
                     nlevels(as.factor(pd[,group_by[3]])))
    names(grp_exp) <- LAPPLY(
      levels(as.factor(pd[group_by[3]])),
      function(y){
        paste(LAPPLY(
          levels(as.factor(pd$group_by[2])),
          function(x){
            paste(levels(as.factor(pd$group_by[1])), x)
          }), y)
      })
    
  }
  
  # Group samples into list of named flowFrames
  if(all(group_by == FALSE)){
    fr.list <- lapply(seq_len(length(fs)), function(x){
      fs[[x]]
    })
  }else{
    fr.list <- cyto_group_by(fs,
                             group_by = group_by)
    fr.list <- lapply(fr.list, function(z){
      cyto_convert(z, "flowFrame")
    })
  }
  
  # Group overlay into list of named flowFrame lists
  if(!is.null(overlay)){
    overlay <- .cyto_overlay_check(x = fs,
                                   overlay = overlay,
                                   display = NULL)
    
    if(all(group_by != FALSE)){
      overlay <- .cyto_overlay_merge(x = fs,
                                     overlay = overlay,
                                     group_by = group_by,
                                     display = NULL)
    }
  }else{
    overlay <- NA
  }
  
  # Replace elements of grp_exp with matching fr.list elements
  if(all(group_by == FALSE)){
    grp_exp[seq_len(length(fs))] <- fr.list
  }else{
    ind <- match(names(fr.list), names(grp_exp))
    grp_exp[ind] <- fr.list
  }
  
  # Replace NA with grp_exp with empty flowFrame
  if(any(is.na(grp_exp))){
    
    # Make flowFrame with same range as fs[[1]] - prevent axes issues
    mn <- LAPPLY(
      colnames(fs[[1]]), 
      function(x){
        min(exprs(fs[[1]])[,x])
      }
    )
    mx <- LAPPLY(
      colnames(fs[[1]]),
      function(x){
        max(exprs(fs[[1]])[,x])
      }
    )
    
    empty_flowFrame <- matrix(c(mn,mx),
                              ncol = BiocGenerics::ncol(fs[[1]]),
                              nrow = 2,
                              byrow = TRUE)
    colnames(empty_flowFrame) <- colnames(fs[[1]])
    empty_flowFrame <- flowCore::flowFrame(exprs = empty_flowFrame)
    grp_exp[is.na(grp_exp)] <- rep(list(empty_flowFrame), 
                                   length(which(is.na(grp_exp))))
    
  }
  
  # Logical indicating which panels require X/Y labels
  panel_label <- rep(list(c(FALSE,FALSE)), length(grp_exp))
  
  # Plots requiring y panel label
  panel_label_y <- LAPPLY(grd, function(x){
    x[,ncol(x)]
  })
  lapply(panel_label_y, function(y){
    panel_label[[y]][2] <<- TRUE
  })
  
  # Plots requiring x axis
  panel_label_x <- LAPPLY(grd, function(x){
    x[1,]
  })
  lapply(panel_label_x, function(y){
    panel_label[[y]][1] <<- TRUE
  })
  
  # Single group_by panel label 
  if(length(group_by) == 1 & all(group_by != FALSE)){
    
    if(format == "row"){
      panel_label <- rep(list(c(TRUE,FALSE)),
                         nlevels(as.factor(pd[,group_by])))
      panel_label[[length(panel_label)]][2] <- TRUE
    }else if(format == "column"){
      panel_label <- rep(list(c(FALSE,TRUE)),
                         nlevels(as.factor(pd[,group_by])))
      panel_label[[1]][1] <- TRUE
    }
    
  }
  
  # X panel labels
  if(missing(panel_x_label)){
    if(all(group_by == FALSE)){
      panel_x_label <- NA
    }else if(length(group_by) == 1){
      if(format == "row"){
        panel_x_label <- levels(as.factor(pd[,group_by]))
      }else if(format == "column"){
        panel_x_label <- NA
      }
    }else{
      panel_x_label <- levels(as.factor(pd[,group_by[1]]))
    }
  }
  
  # Y panel labels
  if(missing(panel_y_label)){
    if(all(group_by == FALSE)){
      panel_y_label <- NA
    }else if(length(group_by) == 1){
      if(format == "row"){
        panel_y_label <- NA
      }else if(format == "column"){
        panel_y_label <- levels(as.factor(pd[,group_by]))
      }
    }else{
      panel_y_label <- levels(as.factor(pd[,group_by[2]]))
    }
  }
  
  # Header
  if(is.null(header)){
    
    # Use 3rd variable for plot headers
    if(length(group_by) == 3){
      header <- rep(levels(as.factor(pd[,group_by[3]])), each = ns)
    }
    
  }else if(length(header) != ns){
    
    header <- rep(header, each = ns)
    
  }
  
  # cyto_plot arguments - rep length(grp_exp) - SKIP
  if(1 == 2){
      ovn <- length(overlay) + 1
  args <- .args_list()
  args <- .cyto_plot_args_split(args,
    channels = channels,
    n = ovn * np,
    plots = np,
    layers = ovn,
    gates = length(gate))
  }
  
  # run through each element of grp_exp - empty plot if NA
  cnt <- 0
  sht <- 0
  mapply(function(x, 
                  overlay,
                  axes_text, 
                  panel_label){
    
    # Update counter
    cnt <<- cnt + 1
    
    # Update sheet number
    if(x %in% seq_len(ns)*np){
      sht <<- sht + 1
    }
    
    # Replace NA with grp_exp with empty flowFrame
    if(is.na(grp_exp[[x]])){
      
      # Make flowFrame with same range as fs[[1]] - prevent axes issues
      mn <- LAPPLY(
        colnames(fs[[1]]), 
        function(x){
          min(exprs(fs[[1]])[,x])
        }
      )
      mx <- LAPPLY(
        colnames(fs[[1]]),
        function(x){
          max(exprs(fs[[1]])[,x])
        }
      )
      
      empty_flowFrame <- matrix(c(mn,mx),
                                ncol = BiocGenerics::ncol(fs[[1]]),
                                nrow = 2,
                                byrow = TRUE)
      colnames(empty_flowFrame) <- colnames(fs[[1]])
      empty_flowFrame <- flowCore::flowFrame(exprs = empty_flowFrame)
      
      # Call to .cyto_plot_empty
      .cyto_plot_empty(empty_flowFrame,
                       channels = channels,
                       axes_trans = axes_trans,
                       axes_text = axes_text,
                       overlay = NA,
                       xlim = xlim,
                       ylim = ylim,
                       ...)
      
    }else{
      # Call to cyto_plot flowFrame method
      cyto_plot(grp_exp[[x]],
                channels = channels,
                axes_trans = axes_trans,
                overlay = overlay,
                gate = gate,
                title = NA,
                axes_text = axes_text, ...)
    }
    
    # Add X panel labels where necessary
    if(panel_label[1] & !is.na(panel_x_label)){
      
      mtext(panel_x_label[x/sht],
            outer = FALSE,
            side = 3,
            line = 0.5,
            font = panel_label_font[1],
            cex = panel_label_size[1],
            col = panel_label_col[1])
      
    }
    
    # Add Y panel labels where necessary
    if(panel_label[2] & !is.na(panel_y_label)){
      
      # mtext does not support text rotation - need to use text for y labels
      text(x = par("usr")[2] + 0.07*(par("usr")[2] - par("usr")[1]),
           y = mean(par("usr")[3:4]),
           labels = panel_y_label[x/(ncol*sht)],
           font = panel_label_font[2],
           cex = 1.5*panel_label_size[2],
           col = panel_label_col[2],
           srt = 270,
           xpd = NA)
      
    }
    
    # Add axes labels, panel labels, titles and headers when plot is full
    if(cnt %% prod(layout) == 0){
      
      # X Axis label
      mtext(xlab, 
            outer = TRUE, 
            side = 1, 
            line = 3.5,
            font = axes_label_text_font,
            cex = axes_label_text_size,
            col = axes_label_text_col)
      
      # Y Axis label
      mtext(ylab, 
            outer = TRUE, 
            side = 2, 
            line = 3.5,
            font = axes_label_text_font,
            cex = axes_label_text_size,
            col = axes_label_text_col)
      
      # Add X Title
      if(is.na(title[1]) & all(group_by != FALSE)){
        title[1] <- group_by[1]
      }
      
      if(!is.na(title[1])){
        mtext(title[1], 
              outer = TRUE, 
              side = 3, 
              line = 2.5,
              font = title_text_font[1],
              cex = title_text_size[1],
              col = title_text_col[1])
      }
      
      # Add Y Title
      if(is.na(title[2]) & all(group_by != FALSE)){
        title[2] <- group_by[2]
      }
      
      if(!is.na(title[2])){
        
        # mtext does not support text rotation - need to use text for y labels
        text(x = par("usr")[2] + 0.2*(par("usr")[2] - par("usr")[1]),
             y = 2*par("usr")[4], # HALF GRID SIZE
             labels = title[2],
             font = title_text_font[2],
             cex = 1.8*title_text_size[2],
             col = title_text_col[2],
             srt = 270,
             xpd = NA)
        
      }
      
      # Add Header
      if(!is.null(header)){
        mtext(header[sht], 
              outer = TRUE, 
              side = 3, 
              line = 5.5,
              font = header_text_font,
              cex = header_text_size,
              col = header_text_col)
      }
      
      if(popup & cnt != ns*np){
        .cyto_plot_window()
      }
      
    }
    
  }, seq_len(length(grp_exp)),
  overlay,
  axes_text,
  panel_label
  )
  
  # Reset cyto_plot_grid option
  options("cyto_plot_grid" = FALSE)
  
  # Reset parameters to default
  par(mar = c(5.1,4.1,4.1,2.1))
  par(oma = c(0,0,0,0))
  
  # Turn off graphics device for saving
  if(getOption("cyto_plot_save")){
    
    if(inherits(x, 
                basename(getOption("cyto_plot_method")))){
      
      # Close graphics device
      dev.off()
      
      # Reset cyto_plot_save
      options("cyto_plot_save" = FALSE)
      
      # Reset cyto_plot_method
      options("cyto_plot_method" = NULL)
      
    }
    
  }
  
} 

#' @rdname cyto_plot_grid
#' @export
cyto_plot_grid.GatingSet <- function(x,
                                     parent,
                                     alias,
                                     channels,
                                     group_by,
                                     header,
                                     title_text_x,
                                     title_text_y, ...) {
  
  # Signal cyto_plot _grid is being used
  options("cyto_plot_grid" = TRUE)
  
  # Reset cyto_plot_grid option
  options("cyto_plot_grid" = FALSE)
  
}
