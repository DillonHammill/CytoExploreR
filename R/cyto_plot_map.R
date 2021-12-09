## CYTO_PLOT_MAP ---------------------------------------------------------------

#' Plot dimension reduced maps in all dimensions
#' 
#' @export
cyto_plot_map <- function(x,
                          parent = "root",
                          channels = NULL,
                          select = NULL,
                          merge_by = "name",
                          overlay = NA,
                          order = "channels",
                          point_col = NA,
                          point_col_scale = NA,
                          layout,
                          header,
                          title,
                          ...) {
  
  # CYTO_PLOT_COMPLETE ----------------------------------------------------------
  
  # CYTO_PLOT METHOD & EXIT
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "map")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }

  # PREPARE PLOT PARAMETERS ----------------------------------------------------
  
  # GROUPS
  grps <- cyto_groups(
    x,
    select = select,
    group_by = merge_by,
    details = TRUE
  )
  
  # GROUPING VARIABLES
  if(cyto_class(merge_by, "list", TRUE)) {
    vars <- names(merge_by)
  } else {
    vars <- merge_by
  }
  
  # POINT_COL - DEFAULT
  if(.all_na(point_col)) {
    # DEFAULT - ALL FLUORESCENT  CHANNELS
    point_col <- cyto_fluor_channels(x)
  # POINT_COL SUPPLIED
  } else {
    point_col <- cyto_channels_extract(x, point_col)
  }
  
  # POINT_COL_SCALE - LIST - DIFFERENT COLOUR SCALE PER CHANNEL
  if(cyto_class(point_col_scale, "list", TRUE)) {
    point_col_scale <- rep(point_col_scale, length(point_col))
    names(point_col_scale) <- point_col
  # POINT_COL_SCALE - VECTOR/FUNCTION - SAME COLOUR SCALE PER CHANNEL
  } else {
    point_col_scale <- structure(
      lapply(seq_along(point_col), function(z){
        return(point_col_scale)
      }),
      names = point_col
    )
  }
  
  # LAYOUT
  if(missing(layout)) {
    if(length(point_col) == 1) {
      order <- "groups"
    } else if(length(grps) == 1) {
      order <- "channels"
    }
    # CHANNEL ORDER
    if(grepl("^c", order, ignore.case = TRUE)) {
      layout <- .cyto_plot_layout(point_col)
    # GROUP ORDER
    } else {
      layout <- .cyto_plot_layout(grps)
    }
  }
  
  # TOTAL PLOTS PER PAGE
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(unique(unlist(layout)))
  }
  
  # PAGES - CHANNELS ON SAME PAGE
  if(grepl("^c", order, ignore.case = TRUE)) {
    # NUMBER OF GROUPS TO PLOT
    n <- length(grps)
    # PAGES PER GROUP
    pg <- ceiling(length(grps)/np)
    # TOTAL PAGES
    tpg <- n * pg
  # PAGES - GROUPS ON SAME PAGE
  } else {
    # NUMBER OF GROUPS TO PLOT
    n <- length(point_col)
    # PAGES PER GROUP
    pg <- ceiling(length(grps)/np)
    # TOTAL PAGES
    tpg <- n * pg
  }
  
  # DEFAULT HEADERS
  if(missing(header)) {
    # CHANNEL ORDER - GROUPS AS HEADERS
    if(grepl("^c", order, ignore.case = TRUE)) {
      header <- rep(
        names(grps),
        each = pg
      )
    # GROUP ORDER - CHANNELS AS HEADERS
    } else {
      header <- rep(
        cyto_markers_extract(
          x,
          channels = point_col,
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
  
  # TITLES - DEFAULT
  if(missing(title)) {
    # CHANNEL ORDER - POINT_COL HEADERS
    if(grepl("^c", order, ignore.case = TRUE)) {
      title <- rep(point_col, each = length(grps))
    # GROUP ORDER
    } else {
      title <- rep("", length.out = length(grps) * length(point_col))
    }
  # TITLES SUPPLIED
  } else {
    title <- rep(title, length.out = length(grps) * length(point_col))
  }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # CALL CYTO_PLOT - CHANNEL ORDER
  if(grepl("^c", order, ignore.case = TRUE)) {
    plots <- structure(
      lapply(seq_along(grps), function(z){
        # HEADER COUNTER
        cnt <- (z - 1) * pg
        # RECORD PLOTS PER GROUP
        p <- structure(
          lapply(seq_along(point_col), function(w){
            # CONSTRUCT PLOT
            cyto_plot(
              x,
              parent = parent,
              select = grps[[z]][, "name"],
              merge_by = vars,
              channels = channels,
              overlay = overlay,
              layout = layout,
              header = header[cnt + ceiling(w/np)],
              title = title[z * length(point_col) + seq_len(length(point_col))],
              page = if(w == length(point_col)) {
                TRUE
              } else {
                FALSE
              },
              point_col = point_col[w],
              point_col_scale <- point_col_scale[[w]],
              ...
            )
          }),
          names = point_col
        )
        # PREPARE RECORDED PLOTS
        p[LAPPLY(p, "is.null")] <- NULL
        p <- lapply(p, `[[`, 1)
        return(p)
      }),
      names = names(grps)
    )
  # CALL CYTO_PLOT - GROUP ORDER
  } else {
    plots <- structure(
      lapply(seq_along(point_col), function(z){
        # HEADER COUNTER
        cnt <- (z - 1) * pg
        # RECORD PLOTS PER CHANNEL
        p <- cyto_plot(
          x,
          parent = parent,
          select = select,
          merge_by = merge_by,
          channels = channels,
          overlay = overlay,
          layout = layout,
          header = header[cnt + seq_len(pg)],
          title = title[z * length(grps) + seq_len(length(grps))],
          page = TRUE,
          point_col = point_col[z],
          point_col_scale <- point_col_scale[[z]],
          ...
        )
        # PREPARE RECORDED PLOTS
        p[LAPPLY(p, "is.null")] <- NULL
        return(p)
      }),
      names = point_col
    )
  }
  
  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}
