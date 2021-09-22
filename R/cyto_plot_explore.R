## CYTO_PLOT_EXPLORE -----------------------------------------------------------

#' Explore cytometry data in a series of bivariate scatterplots
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @export
cyto_plot_explore <- function(x,
                              parent = "root",
                              channel_x = NULL,
                              channels_y = NULL,
                              select = NULL,
                              merge_by = "name",
                              order = "channels",
                              layout,
                              header,
                              ...) {
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "explore")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # PREPARE PLOT PARAMETERS ----------------------------------------------------
  
  # X CHANNELS - DEFAULT TO FLUORESCENT CHANNELS
  if(is.null(channels_x)) {
    channels_x <- cyto_fluor_channels(x)
  }
  
  # Y CHANNELS - DEFAULT ALL CHANNELS
  if(is.null(channels_y)) {
    channels_y <- cyto_channels(x)
  }
  
  # EXTRACT NAMES OF VARIABLES FOR MERGING
  if(cyto_class(merge_by, "list", TRUE)) {
    vars <- names(merge_by)
  } else {
    vars <- merge_by
  }
  
  # GROUPS
  grps <- cyto_groups(
    x,
    select = select,
    group_by = group_by,
    details = TRUE
  )
  
  # LAYOUT 
  if(missing(layout)) {
    # SWITCH LAYOUTS FOR SINGLE PLOTS
    if(length(channels_x) == 1) {
      order <- "groups"
    } else if(length(grps) == 1) {
      order <- "channels"
    }
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      layout <- .cyto_plot_layout(channels_y)
    # GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(grps)
    }
  }
  
  # TOTAL PLOTS PER PAGE - BASED ON LAYOUT
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES
  if(grepl("^c", order, ignore.case = TRUE)) {
    n <- length(grps)
    pg <- ceiling(length(channels_y)/np)
    tpg <- n * pg
  } else {
    n <- length(channels)
    pg <- ceiling(length(grps)/np)
    tpg <- n * pg
  }
  
  # PREPARE HEADERS
  if(missing(header)) {
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      header <- rep(
        names(grps),
        each = pg
      )
    # GROUP ORDER
    } else {
      header <- rep(
        cyto_markers_extract(
          x,
          channels = channels,
          append = TRUE
        ),
        each = pg
      )
    }
  # SUPPLIED HEADERS
  } else {
    # NO HEADERS
    if(.all_na(header)) {
      header <- rep(NA, length.out = tpg)
    } else {
      stop(
        paste0(
          "Supply a header for each ",
          ifelse(grepl("^c", order, ignore.case = TRUE),
                 "group!",
                 "channel!")
        )
      )
    }
  }
  
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # CALL CYTO_PLOT - CHANNEL ORDER
  if(grepl("^c", order, ignore.case = TRUE)) {
    # CONSTRUCT & RECORD PLOTS
    plots <- lapply(seq_along(grps), function(z){
      # HEADER COUNTER
      cnt <- (z - 1) * pg
      # RECORDED PLOTS PER GROUP - LOOP THROUGH CHANNELS_X
      lapply(seq_along(channels_x), function(w){
        # X CHANNEL
        x_chan <- channels_x[w]
        # LOOP THROUGH CHANNELS_Y
        lapply(seq_along(channels_y), function(v){
          # Y CHANNEL
          y_chan <- channels_y[v]
          # CONSTRUCT PLOT
          p <- cyto_plot(
            x,
            parent = parent,
            select = grps[[z]][, "name"],
            channels = c(x_chan, y_chan),
            overlay = overlay,
            layout = layout,
            header = header[cnt + ceiling(w/np)],
            page = if(v == length(channels_y)) {
              TRUE
            } else {
              FALSE
            },
            ...
          )
        })
      })
    })
  # CALL CYTO_PLOT - GROUP ORDER
  } else {
    # CONSTRUCT & RECORD PLOTS
    plots <- lapply(seq_along(channels_x), function(z){
      # HEADER COUNTER
      cnt <- (z - 1) * pg
      # X CHANNEL
      x_chan <- channels_x[z]
      # LOOP THROUGH CHANNELS_Y
      plots <- lapply(seq_along(channels_y), function(w){
        # Y CHANNEL
        y_chan <- channels_y[w]
        # CONTRUCT PLOT
        p <- cyto_plot(
          x,
          parent = parent,
          select = select,
          channels = c(x_chan, y_chan),
          overlay = overlay,
          layout = layout,
          header = header[cnt + ceiling(w/np)],
          page = if(w == length(channels_y)) {
            TRUE
          } else {
            FALSE
          },
          ...
        )
        # PREPARE RECORDED PLOTS
        p[LAPPLY(p, "is.null")] <- NULL
        return(p)
      })
    })
  }
  
  # FORMAT RECORDED PLOTS
  plots[LAPPLY(plots, "is.null")] <- NULL
  if(length(plots) == 0) {
    plots <- NULL
  }
  
  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  return(plots)
  
}