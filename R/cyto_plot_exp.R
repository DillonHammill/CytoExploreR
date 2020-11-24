## CYTO_PLOT_EXPLORE -----------------------------------------------------------

#' @importFrom graphics par
#' @export
cyto_plot_exp <- function(x,
                              parent = "root",
                              overlay = NA,
                              channels_x = NULL,
                              channels_y = NULL,
                              popup = FALSE,
                              layout = NULL,
                              merge_by = NULL,
                              axes_trans = NA,
                              axes_limits = "machine",
                              title = NA,
                              header,
                              header_text_size = 1,
                              header_text_font = 2,
                              header_text_col = "black",
                              page = "samples",
                              ...) {
  
  # CYTO_PLOT_EXIT -------------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- .par()
  
  # NEW PLOT METHOD
  if(is.null(getOption("cyto_plot_method"))) {
    # CYTO_PLOT_METHOD
    options("cyto_plot_method" = "explore")
    # CYTO_PLOT_EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_method" = NULL)
      options("cyto_plot_par" = NULL)
    })
  }
  
  # EXTRACT & PREPARE DATA -----------------------------------------------------
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformer_extract(x)
  }
  
  # CYTOFRAME/CYTOSET
  x <- cyto_extract(x, parent = parent)
  
  # OVERLAY
  if(!.all_na(overlay)) {
    if(is.character(overlay)) {
      
    }
  }
  
  # MERGE_BY
  if(!is.null(merge_by)) {
    if(cyto_class(x, "flowFrame")) {
      cf_list <- cyto_list(x)
    } else {
      cf_list <- cyto_merge_by(x, 
                               merge_by = merge_by)
    }
  } else {
    cf_list <- cyto_list(x)
  }
  
  # OVERLAY
  
  # CHANNELS -------------------------------------------------------------------
  
  # CHANNELS_X - FLUORESCENT CHANNELS
  if(is.null(channels_x)) {
    channels_x <- cyto_fluor_channels(cf_list)
  }
  
  # CHANNELS_Y - CHANNELS
  if(is.null(channels_y)) {
    channels_y <- cyto_channels(cf_list)
  }
  
  # PLOT ARGUMENTS -------------------------------------------------------------
  
  # HEADER
  if(missing(header)) {
    if(page == "samples") {
      header <- channels_y
    } else if(page == "channels") {
      header <- cyto_names(cf_list)
    }
  }
  
  # LAYOUT
  if(page == "samples") {
    layout <- .cyto_plot_layout(cf_list, layout)
  } else if(page == "channels") {
    layout <- .cyto_plot_layout(channels_x, layout)
  }
  
  # OUTER MARGINS
  if(!.all_na(header)) {
    oma <- c(0, 0, 3, 0)
  }
  
  # PLOT CONSTRUCTION ----------------------------------------------------------
  
  # CYTO_PLOT_NEW
  if(getOption("cyto_plot_method") == "explore" &
     !getOption("cyto_plot_custom")) {
    cyto_plot_new(popup,
                  layout = layout, 
                  oma = oma)
  }
  
  # PAGE - SAMPLES
  if(page == "samples") {

  # PAGE - CHANNELS  
  } else if(page == "channels") {
    print("YASS")
    # PLOTS
    p <- lapply(seq_along(cf_list), function(z){
      # CYTOFRAME
      cf <- cf_list[[z]]
      print(cf)
      # OVERLAY
      
      p <- lapply(seq_along(channels_y), function(y) {
        # Y CHANNEL
        chan_y <- channels_y[y]
        print(chan_y)
        p <- lapply(seq_along(channels_x), function(w){
          # X CHANNEL
          chan_x <- channels_x[w]
          print(chan_x)
          # CONSTRUCT 1D PLOT
          if(chan_x == chan_y) {
            cyto_plot(cf,
                      channels = chan_x,
                      axes_trans = axes_trans,
                      axes_limits = axes_limits,
                      title = title,
                      legend = FALSE,
                      ...)
            # CONSTRUCT 2D PLOT 
          } else {
            cyto_plot(cf,
                      channels = c(chan_x, chan_y),
                      axes_trans = axes_trans,
                      axes_limits = axes_limits,
                      title = title,
                      legend = FALSE,
                      ...)
          }
          # HEADER
          if(par("page") | w == length(channels_x)) {
            # HEADER
            if(!.all_na(header)) {
              mtext(header[z],
                    outer = TRUE,
                    cex = header_text_size,
                    font = header_text_font,
                    col = header_text_col)
            }
            # RECORD
            p <- cyto_plot_record()
            # NEW DEVICE
            if(w != length(channels_x)) {
              cyto_plot_new()
            }
          } else {
            p <- NULL
          }
          return(p)
        })
        # REMOVE EMPTY RECORD
        p <- p[!LAPPLY(p, "is.null")]
        names(p) <- channels_y[y]
        # NEW DEVICE
        cyto_plot_new()
        # RETURN RECORDED PLOTS
        return(p)
      })
      return(p)
    })
    names(p) <- cyto_names(cf_list)
  }
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if(getOption("cyto_plot_save") &
     getOption("cyto_plot_method") == "explore") {
    # CLOSE GRAPHICS DEVICE
    dev.off()
    # RESET CYTO_PLOT_SAVE
    options("cyto_plot_save" = FALSE)
  }
  
  # RETURN RECORDED PLOTS
  invisible(p)
  
}