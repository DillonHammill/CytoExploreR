#'  Grid Layout for cyto_plot
#'
#'  \code{cyto_plot_grid} is wrapper for cyto_plot which provides a grid layout
#'  without excess white space. \code{cyto_plot_grid} is particularly useful for
#'  visualising data grouped using \code{group_by}.
#'
#'  @param x object of class \code{flowSet} or \code{GatingSet}.
#'  @param ... additional method-specific arguments for \code{cyto_plot_grid}.
#'
#'  @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#'  @export
setGeneric(
  name = "cyto_plot_grid",
  def = function(x, ...) {
    standardGeneric("cyto_plot_grid")
  }
)

#'  Grid layout for cyto_plot - flowSet Method
#'
#'  \code{cyto_plot_grid} is wrapper for cyto_plot which provides a grid layout
#'  without excess white space. \code{cyto_plot_grid} is particularly useful for
#'  visualising and facetting data grouped using \code{group_by}.
#'
#'  @param x object of class \code{flowSet}.
#'  @param ... additional method-specific arguments for \code{cyto_plot}.
#'  @param format indicates whether to plot in a "row" or "column" when a single
#'    variable is supplied.
#'
#'  @importFrom flowCore parameters flowFrame
#'  @importFrom flowWorkspace pData
#'  @importfrom graphics par mtext text
#'
#'  @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#'  @export
setMethod(
  cyto_plot_grid,
  signature = "flowSet",
  definition = function(x,
                        channels = NULL,
                        axes_trans = NULL,
                        group_by = FALSE,
                        overlay = NULL,
                        gate = NA,
                        limits = "machine",
                        display = NULL,
                        layout = NULL,
                        popup = FALSE,
                        xlim = NULL,
                        ylim = NULL,
                        title = c(NA,NA),
                        xlab,
                        ylab,
                        format = "row",
                        header = NULL,
                        panel_x_label = NULL,
                        panel_y_label = NULL,
                        density_modal = TRUE,
                        density_smooth = 1.5,
                        density_stack = 0,
                        density_fill = NA,
                        density_fill_alpha = 1,
                        density_line_type = 1,
                        density_line_width = 1,
                        density_line_col = "black",
                        point_shape = ".",
                        point_col = NA,
                        point_alpha = NA,
                        contour_lines = 15,
                        contour_line_type = 1,
                        contour_line_width = 1,
                        contour_line_col = "black",
                        axes_text_font = 1,
                        axes_text_size = 1,
                        axes_text_col = "black",
                        axes_label_text_font = 1,
                        axes_label_text_size = 1.1,
                        axes_label_text_col = "black",
                        title_text_font = c(2,2),
                        title_text_size = c(1.2,1.2),
                        title_text_col = c("black","black"),
                        legend = FALSE,
                        gate_line_type = 1,
                        gate_line_width = 2.5,
                        gate_line_col = "red",
                        label = TRUE,
                        label_text = NA,
                        label_stat = "percent",
                        label_text_font = 2,
                        label_text_size = 1,
                        label_text_col = "black",
                        label_box_x = NA,
                        label_box_y = NA,
                        label_box_alpha = 0.6,
                        border_line_type = 1,
                        border_line_width = 1,
                        border_line_col = "black",
                        header_text_font = 2,
                        header_text_size = 1.5,
                        header_text_col = "black",
                        panel_label_font = c(2,2),
                        panel_label_size = c(1,1),
                        panel_label_col = c("black","black"),
                         ...) {
    
    # cyto_plot_grid is designed to plot multiple groups only (1-3 variables)
    if(all(group_by[1] %in% c(TRUE,"all"))){
      stop("Use cyto_plot when grouping all samples into the same group.")
    }else if(length(group_by) > 3){
      stop("cyto_plot_grid supports a maximum of 3 grouping variables.")
    }
    
    # Signal cyto_plot_grid is being used
    options("CytoRSuite_cyto_plot_grid" = TRUE)
    
    # Set plot method
    if(is.null(getOption("CytoRSuite_cyto_plot_method"))){
      options("CytoRSuite_cyto_plot_method" = "grid/flowSet")
    }
    
    # Channels missing
    if(missing(channels)){
      stop("Supply channel/marker(s) to construct the plot")
    }
    
    # Assign x to fs
    fs <- x
    smp <- length(fs)
    
    # Extract pData
    pd <- pData(fs)
    
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
          layout <- c(nlevels(as.factor(pd[,group_by[1]])),
                      nlevels(as.factor(pd[,group_by[2]])))
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
    grd <- matrix(seq_len(ncol*nrow), ncol = ncol, nrow=nrow, byrow = TRUE)
    grd <- lapply(seq_len(ns), function(x){
      grd*ns
    })
    
    # Plots requiring Y axis
    axes_y <- unlist(lapply(grd, function(x){
      x[,1]
    }))
    lapply(axes_y, function(x){
      axes_text[[x]][2] <<- TRUE
    })

    # Plots requiring X axis
    axes_x <- unlist(lapply(grd, function(x){
      x[nrow(x),]
    }))
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
    
    # Pop-up window?
    if(popup){
      .cyto_plot_window()
    }
    
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
      names(grp_exp) <- unlist(lapply(
        levels(as.factor(pd[,group_by[2]])),
        function(x){
          paste(levels(as.factor(pd[,group_by[1]])), x)
        }))
    }else if(length(group_by) == 3){
      grp_exp <- rep(list(NA), nlevels(as.factor(pd[,group_by[1]])) *
                       nlevels(as.factor(pd[,group_by[2]])) *
                       nlevels(as.factor(pd[,group_by[3]])))
      names(grp_exp) <- unlist(lapply(
        levels(as.factor(pd[group_by[3]])),
        function(y){
          paste(unlist(lapply(
            levels(as.factor(pd$group_by[2])),
            function(x){
              paste(levels(as.factor(pd$group_by[1])), x)
            })), y)
        }))
      
    }
    
    # Group samples into list of named flowFrames
    if(all(group_by == FALSE)){
      fr.lst <- lapply(seq_len(length(fs)), function(x){
        fs[[x]]
      })
    }else{
      fr.lst <- .cyto_merge(x = fs,
                            group_by = group_by,
                            display = NULL)
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
    
    # Replace elements of grp_exp with matching fr.lst elements
    if(all(group_by == FALSE)){
      grp_exp[seq_len(length(fs))] <- fr.lst
    }else{
      ind <- match(names(fr.lst), names(grp_exp))
      grp_exp[ind] <- fr.lst
    }
    
    # Replace NA with grp_exp with empty flowFrame
    if(any(is.na(grp_exp))){
      
      # Make flowFrame with same range as fs[[1]] - prevent axes issues
      mn <- unlist(lapply(
        colnames(fs[[1]]), 
        function(x){
          min(exprs(fs[[1]])[,x])
          }
        ))
      mx <- unlist(lapply(
        colnames(fs[[1]]),
        function(x){
          max(exprs(fs[[1]])[,x])
        }
      ))
      
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
    panel_label_y <- unlist(lapply(grd, function(x){
      x[,ncol(x)]
    }))
    lapply(panel_label_y, function(y){
      panel_label[[y]][2] <<- TRUE
    })
    
    # Plots requiring x axis
    panel_label_x <- unlist(lapply(grd, function(x){
      x[1,]
    }))
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
    
    # cyto_plot arguments - rep length(grp_exp)
    args <- as.list(environment())
    args <- .arg_split(args[-c(match(
      c("fr.lst",
        "grp_exp",
        "fr.channels",
        "fr.data",
        "panel_label_x",
        "panel_label_y",
        "axes_x",
        "axes_y",
        "grd",
        "ns",
        "tnp",
        "np",
        "ncol",
        "nrow",
        "pd",
        "smp",
        "fs",
        "x",
        "channels",
        "group_by",
        "layout",
        "popup",
        "title",
        "xlab", 
        "ylab",
        "format",
        "header",
        "panel_x_label",
        "panel_y_label",
        "title_text_font",
        "title_text_size",
        "title_text_col",
        "header_text_font",
        "header_text_size",
        "header_text_col",
        "panel_label_font",
        "panel_label_size",
        "panel_label_col"
      ), names(args)))],
      channels = channels,
      n,
      plots,
      layers,
      gates)
    
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
      
      # Overlay
      if(is.na(overlay)){
        overlay <- NULL
      }
      
      # Replace NA with grp_exp with empty flowFrame
      if(is.na(grp_exp[[x]])){
        
        # Make flowFrame with same range as fs[[1]] - prevent axes issues
        mn <- unlist(lapply(
          colnames(fs[[1]]), 
          function(x){
            min(exprs(fs[[1]])[,x])
          }
        ))
        mx <- unlist(lapply(
          colnames(fs[[1]]),
          function(x){
            max(exprs(fs[[1]])[,x])
          }
        ))
        
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
                         overlay = NULL,
                         xlim = xlim,
                         ylim = ylim,
                         limits = limits,
                         axes_text = axes_text,
                         axes_text_font = axes_text_font,
                         axes_text_size = axes_text_size,
                         axes_text_col = axes_text_col,
                         axes_label_text_font = axes_label_text_font,
                         axes_label_text_size = axes_label_text_size,
                         axes_label_text_col = axes_label_text_col,
                         border_line_type = border_line_type,
                         border_line_width = border_line_width,
                         border_line_col = border_line_col)
        
      }else{
        # Call to cyto_plot flowFrame method
        cyto_plot(grp_exp[[x]],
                  channels = channels,
                  overlay = overlay,
                  legend = FALSE,
                  density_stack = density_stack,
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
        text(x = par("usr")[3] + 0.07*(par("usr")[3] - par("usr")[1]),
             y = mean(par("usr")[c(2,4)]),
             labels = panel_y_label[x/(ncol*sht)],
             font = panel_label_font[2],
             cex = 1.25*panel_label_size[2],
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
          text(x = par("usr")[3] + 0.18*(par("usr")[3] - par("usr")[1]),
               y = par("usr")[4],
               labels = title[2],
               font = title_text_font[2],
               cex = 1.25*title_text_size[2],
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
    options("CytoRSuite_cyto_plot_grid" = FALSE)
    
    # Reset parameters to default
    par(mar = c(5.1,4.1,4.1,2.1))
    par(oma = c(0,0,0,0))
    
    # Turn off graphics device for saving
    if(getOption("CytoRSuite_cyto_plot_save")){
      
      if(inherits(x, 
                  basename(getOption("CytoRSuite_cyto_plot_method")))){
        
        # Close graphics device
        dev.off()
        
        # Reset CytoRSuite_cyto_plot_save
        options("CytoRSuite_cyto_plot_save" = FALSE)
        
        # Reset CytoRSuite_cyto_plot_method
        options("cytoRSuite_cyto_plot_method" = NULL)
        
      }
      
    }
    
  })  

#' Grid layout for cyto_plot - GatingSet Method
#' 
#'  \code{cyto_plot_grid} is wrapper for cyto_plot which provides a grid layout
#'  without excess white space. \code{cyto_plot_grid} is particularly useful for
#'  visualising data grouped using \code{group_by}.
#'  
#'  @param x object of class \code{GatingSet}.
#'  @param ... additional method-specific arguments for cyto_plot.
#'  
#'  @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'  
#'  @export
setMethod(
  cyto_plot_grid,
  signature = "GatingSet",
  definition = function(x,
                        parent,
                        alias,
                        channels,
                        group_by,
                        header,
                        title_text_x,
                        title_text_y, ...) {
            
    # Signal cyto_plot _grid is being used
    options("CytoRSuite_cyto_plot_grid" = TRUE)
    
    # Reset cyto_plot_grid option
    options("CytoRSuite_cyto_plot_grid" = FALSE)
    
  })
