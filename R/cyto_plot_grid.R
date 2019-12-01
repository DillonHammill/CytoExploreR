## CYTO_PLOT_GRID --------------------------------------------------------------

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
#' @importFrom purrr transpose
#' @importFrom magrittr %>%
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
                                   channels,
                                   axes_trans = NA,
                                   group_by = "name",
                                   overlay = NA,
                                   gate = NA,
                                   popup = FALSE,
                                   layout = NULL,
                                   xlab = "",
                                   ylab = "",
                                   limits = "auto",
                                   xlim = NA,
                                   ylim = NA,
                                   title = c(NA,NA),
                                   format = "row",
                                   header = NULL,
                                   panel_x_label = NULL,
                                   panel_y_label = NULL,
                                   negate = FALSE,
                                   density_modal = TRUE,
                                   label,
                                   label_text = NA,
                                   label_stat,
                                   label_position = "auto",
                                   label_text_x = NA,
                                   label_text_y = NA,
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
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- .par(c("oma","mfrow","mar"))
  
  # CHECKS ---------------------------------------------------------------------
  
  # GROUPING VARIABLES (1-3)
  if(all(group_by[1] == "all")){
    stop("Use cyto_plot when grouping all samples into the same group.")
  }else if(length(group_by) > 3){
    stop("cyto_plot_grid supports a maximum of 3 grouping variables.")
  }
  
  # SIGNAL CYTO_PLOT_GRID
  options("cyto_plot_grid" = TRUE)
  
  # METHOD & RESET
  if(is.null(getOption("cyto_plot_method"))){
    # SET PLOT METHOD
    options("cyto_plot_method" = "grid/flowSet")
    # RESET PLOT METHOD & GRAPHICAL PARAMETERS ON EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_method" = NULL)
      options("cyto_plot_grid" = FALSE)
    })
  }else{
    # RESET GRAPHICAL PARAMETERS ON EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_grid" = FALSE)})
  }
  
  # CHANNELS
  if(missing(channels)){
    stop("Supply channel/marker(s) to construct the plot")
  }
  
  # ORGANISE DATA ---------------------------------------------------------------
  
  # ASSIGN X TO FS
  fs <- x
  
  # SAMPLES
  smp <- length(fs)
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(fs)
  
  # GROUP_BY FACTORS - DROP MISSING LEVELS
  lapply(group_by, function(z){
    # GROUP_BY VARIABLES AS FACTORS (LEVELS MAY BE MANUALLY SET)
    if(!is(pd[, z], "factor")){
      pd[, z] <<- factor(pd[, z], levels = unique(pd[, z]))
    }
    # DROP MISSING FACTOR LEVELS
    pd[, z] <<- droplevels(pd[, z])
  })
  
  # CHANNELS
  channels <- cyto_channels_extract(fs, 
                                    channels, 
                                    plot = TRUE)
  
  # LAYOUT/GRID/EXPECTED GROUPS ------------------------------------------------
  
  # EXPECTED GROUPS - NO GROUPING
  if(all(group_by == "name")){
    # LAYOUT - ALL SAMPLES SAME SHEET
    if(is.null(layout)){
      layout <- .cyto_plot_layout(x = fs,
                                  layout = NULL)
    }
    # SHEETS
    sheets <- ceiling(smp/prod(layout))
    # TOTAL PANELS
    panels <- prod(layout) * sheets
    # EXPECTED GROUPS 
    expected_groups <- rep(list(NA), panels)
    names(expected_groups) <- rep(c(cyto_names(fs),
                                    paste0("Blank",
                                           seq_len(panels-smp))),
                                  length.out = panels)
  # EXPECTED GROUPS - ONE GROUPING VARIABLE
  }else if(length(group_by) == 1){
    # LAYOUT MUST BE A SINGLE ROW OR COLUMN
    if(!is.null(layout)){
      if(!any(layout == 1)){
        message("'layout' must be a single row or column.")
        layout <- NULL
      }
    }
    # LAYOUT - ROW/COLUMN
    if(is.null(layout)){
      # ROW LAYOUT
      if(grepl("row", format)){
        layout <- c(1, nlevels(pd[, group_by]))
      # COLUMN LAYOUT
      }else if(grepl("column", format)){
        layout <- c(nlevels(pd[, group_by]), 1)
      }
    }
    # GROUPS
    groups <- nlevels(pd[, group_by])
    # SHEETS
    sheets <- ceiling(groups/prod(layout))
    # PANELS
    panels <- prod(layout) * sheets
    # EXPECTED GROUPS
    expected_groups <- rep(list(NA), panels)
    names(expected_groups) <- rep(c(levels(pd[, group_by]),
                                    paste0("Blank",
                                           seq_len(panels-groups))),
                                  length.out = panels)
  # EXPECTED GROUPS - TWO GROUPING VARIABLES  
  }else if(length(group_by) == 2){
    # LAYOUT
    if(is.null(layout)){
      layout <- c(nlevels(pd[,group_by[2]]),
                  nlevels(pd[,group_by[1]]))
    }
    # GROUPS
    groups <- nlevels(pd[, group_by[1]]) * nlevels(pd[, group_by[2]])
    # SHEETS
    sheets <- ceiling(groups/prod(layout))
    # PANELS
    panels <- prod(layout) * sheets
    # EXPECTED GROUPS
    expected_groups <- rep(list(NA), panels)
    names(expected_groups) <- LAPPLY(
      levels(pd[, group_by[2]]),
      function(x){
        paste(levels(pd[, group_by[1]]), x)
      })
  # EXPECTED GROUPS - THREE GROUPING VARIABLES  
  }else if(length(group_by) == 3){
    # LAYOUT - THIRD VARIABLE INDICATES SHEETS
    if(is.null(layout)){
      layout <- c(nlevels(pd[,group_by[2]]),
                nlevels(pd[,group_by[1]]))
    }
    # GROUPS
    groups <- nlevels(pd[, group_by[1]]) * nlevels(pd[, group_by[2]]) *
      nlevels(pd[, group_by[3]])
    # SHEETS
    sheets <- ceiling(groups/prod(layout))
    # PANELS
    panels <- prod(layout) * sheets
    # EXPECT FULL GRID
    expected_groups <- rep(list(NA), panels)
    names(expected_groups) <- LAPPLY(
      levels(pd[, group_by[3]]),
      function(y){
        paste(LAPPLY(
          levels(pd[, group_by[2]]),
          function(x){
            paste(levels(pd[, group_by[1]]), x)
          }), y)
      })
  }
  
  # GRID DIMENSIONS PER SHEET
  grid <- split(seq_len(panels),
                rep(seq_len(sheets), 
                    each = prod(layout)))
  grid <- lapply(grid, function(z){
    matrix(z,
           nrow = layout[1],
           ncol = layout[2],
           byrow = TRUE)
  })

  # PREPARE DATA FOR CYTO_PLOT CALLS -------------------------------------------
  
  # GROUPING - LIST OF MERGED FLOWFRAMES
  fr_list <- cyto_merge_by(fs,
                           merge_by = group_by)
  
  # TRANSPOSE TO LIST OF FLOWFRAME LISTS
  fr_list <- list(fr_list) %>% transpose()
  
  # OVERLAY
  if (!.all_na(overlay)) {
    # REPEAT FLOWFRAME PER PLOT
    if (inherits(overlay, "flowFrame")) {
      overlay_list <- rep(list(list(overlay)), length(fr_list))
      # FLOWSET TO LIST OF FLOWFRAMES
    } else if (inherits(overlay, "flowSet")) {
      # GROUPING
      overlay_list <- cyto_merge_by(overlay,
                                    merge_by = group_by)
      # LIST OF FLOWFRAME LISTS
      overlay_list <- lapply(overlay_list, function(z) {
        list(z)
      })
      # LIST OF FLOWSETS TO LIST OF FLOWFRAME LISTS
    } else if (inherits(overlay, "list")) {
      # ALLOW LIST OF FLOWFRAMES OF LENGTH FR_LIST
      if (all(LAPPLY(unlist(overlay), function(z) {
        inherits(z, "flowFrame")
      }))) {
        # SAME LENGTH AS FR_LIST
        if (length(overlay) != length(fr_list)) {
          stop(
            paste(
              "'overlay' must be a list of flowFrame lists -",
              "one flowFrame list per plot."
            )
          )
        }
        # NO GROUPING APPLIED
        overlay_list <- overlay
        # LIST OF FLOWSETS
      } else if (all(LAPPLY(overlay, function(z) {
        inherits(z, "flowSet")
      }))) {
        # GROUPING
        overlay_list <- lapply(overlay, function(z){
          cyto_merge_by(z,
                        merge_by = group_by)
        })
        overlay_list <- overlay_list %>% transpose()
        # OVERLAY NOT SUPPORTED
      } else {
        stop(
          paste(
            "'overlay' should be either a flowFrame, flowSet, list of flowFrame",
            "lists or list of flowSet lists."
          )
        )
      }
    }
    # NAME OVERLAY
    names(overlay_list) <- names(fr_list)
    
    # COMBINE FR_LIST & OVERLAY_LIST - LIST OF FLOWFRAME LISTS
    fr_list <- lapply(seq_along(fr_list), function(z){
      c(fr_list[[z]], overlay_list[[z]])
    })
    names(fr_list) <- names(overlay_list)
    
  }
  
  # REPLACE ELEMENTS OF expected_groups WITH MATCHING FR_LIST ELEMENTS
  ind <- match(names(fr_list), names(expected_groups))
  expected_groups[ind] <- fr_list
  
  # REPLACE NA IN expected_groups WITH EMPTY FLOWFRAME
  if(any(LAPPLY(expected_groups, ".all_na"))){
    
    # EXPECTED NUMBER OF LAYERS PER PLOT
    N <- max(LAPPLY(expected_groups, "length"))
    
    # MISSING FACTOR GROUP LEVELS
    empty_group <- which(LAPPLY(expected_groups, ".all_na"))

    # VALID GROUP
    valid_group <- expected_groups[[which(!LAPPLY(expected_groups, ".all_na"))[1]]]
    
    # REPLACE NA WITH LIST OF EMPTY FLOWFRAMES
    lapply(empty_group, function(z){
      # CREATE EMPTY FLOWFRAME
      empty_flowFrame <- cyto_empty(name = names(expected_groups)[z],
                                    channels = cyto_channels(valid_group[[1]]))
      # REPLACE NA WITH LIST OF EMPTY FLOWFRAMES
      expected_groups[[z]] <<- rep(list(empty_flowFrame), N)
    })

  }  

  # PREPARE CYTO_PLOT ARGUMENTS ------------------------------------------------
  
  # REMOVAL NEGATIVE FSC/SSC EVENTS - POINT_COL SCALE ISSUE
  lapply(seq_len(length(channels)), function(z) {
    if (grepl("FSC", channels[z], ignore.case = TRUE) |
        grepl("SSC", channels[z], ignore.case = TRUE)) {
      fr_list <<- lapply(fr_list, function(y) {
        # LIST OF FLOWFRAMES
        lapply(y, function(w) {
          if(nrow(w@exprs) > 0){
            if (min(range(w, type = "data")[, channels[z]]) < 0) {
              coords <- matrix(c(0, Inf), ncol = 1, nrow = 2)
              rownames(coords) <- c("min", "max")
              colnames(coords) <- channels[z]
              nonDebris <- rectangleGate(.gate = coords)
              Subset(w, nonDebris)
            } else {
              return(w)
            }
          }else{
            return(w)
          }
        })
      })
    }
  })
  
  # GATES - LIST OF GATE OBJECT LISTS
  if (!.all_na(gate)) {
    # REPEAT GATE OBJECTS PER PLOT
    if (any(is(gate) %in% c(
      "rectangleGate",
      "polygonGate",
      "ellipsoidGate",
      "quadGate",
      "filters"
    ))) {
      gate <- list(gate)
    }
    # LIST OF GATE OBJECTS - REPEAT PER PLOT
    if (all(LAPPLY(gate, function(z) {
      any(is(z) %in% c(
        "rectangleGate",
        "polygonGate",
        "ellipsoidGate",
        "quadGate",
        "filters"
      ))
    }))) {
      gate <- rep(list(gate), length.out = length(fr_list))
    }
  }
  
  # XLIM
  if (.all_na(xlim)) {
    xlim <- .cyto_range(fr_list,
                        channels = channels[1],
                        limits = limits
    )[, channels[1]]
    # XLIM MANUALLY SUPPLIED
  } else {
    if (!.all_na(axes_trans)) {
      if (channels[1] %in% names(axes_trans)) {
        xlim <- axes_trans[[channels[1]]]$transform(xlim)
      }
    }
  }
  
  # YLIM - 1D CALCULATED LATER
  if (.all_na(ylim)) {
    # 2D PLOT
    if (length(channels) == 2) {
      ylim <- .cyto_range(fr_list,
                          channels = channels[2],
                          limits = limits
      )[, channels[2]]
    }
    # YLIM MANUALLY SUPPLIED
  } else {
    # 2D PLOT
    if (length(channels) == 2) {
      if (!.all_na(axes_trans)) {
        if (channels[2] %in% names(axes_trans)) {
          ylim <- axes_trans[[channels[2]]]$transform(ylim)
        }
      }
    }
  }
  
  # SETUP AXES_TEXT
  axes_text <- rep(list(c(FALSE,FALSE)), panel)
  
  # PLOTS REQUIRING Y AXIS
  axes_y <- LAPPLY(grid, function(x){
    x[, 1]
  })
  lapply(axes_y, function(x){
    axes_text[[x]][2] <<- TRUE
  })
  
  # PLOTS REQUIRING X AXIS
  axes_x <- LAPPLY(grid, function(x){
    x[nrow(x),]
  })
  lapply(axes_x, function(y){
    axes_text[[y]][1] <<- TRUE
  })
  
  # AXES LABELS
  axes_labels <- .cyto_plot_axes_label(x = fs[[1]],
                                       channels = channels,
                                       xlab = xlab,
                                       ylab = ylab,
                                       density_modal = density_modal)
  xlab <- axes_labels[[1]]
  ylab <- axes_labels[[2]]
  
  # PANELS TO LABEL
  panel_label <- rep(list(c(FALSE,FALSE)), panels)
  
  # Y PANEL LABEL
  panel_label_y <- LAPPLY(grid, function(x){
    x[, ncol(x)]
  })
  lapply(panel_label_y, function(y){
    # PANEL Y LABELS REQUIRE MULTIPLE GROUPING VARIABLES OR COLUMN FORMAT
    if(length(group_by) > 1 |
       (length(group_by) == 1 & grepl("column", format))){
      panel_label[[y]][2] <<- TRUE
    }
  })
  
  # X PANEL LABEL
  panel_label_x <- LAPPLY(grid, function(x){
    x[1, ]
  })
  lapply(panel_label_x, function(y){
    # PANEL X LABELS ALWAYS REQUIRED UNLESS COLUMN FORMAT & SINGLE VARIABLE
    if(!length(group_by) == 1 & grepl("column", format)){
      panel_label[[y]][1] <<- TRUE
    }
  })
  
  # X PANEL LABELS
  if(is.null(panel_x_label)){
    if(all(group_by == "name")){
      panel_x_label <- NA
    }else if(length(group_by) == 1){
      if(format == "row"){
        panel_x_label <- levels(pd[,group_by])
      }else if(format == "column"){
        panel_x_label <- NA
      }
    }else{
      panel_x_label <- levels(pd[,group_by[1]])
    }
  }
  
  # Y PANEL LABELS
  if(is.null(panel_y_label)){
    if(all(group_by == "name")){
      panel_y_label <- NA
    }else if(length(group_by) == 1){
      if(format == "row"){
        panel_y_label <- NA
      }else if(format == "column"){
        panel_y_label <- levels(pd[,group_by])
      }
    }else{
      panel_y_label <- levels(pd[,group_by[2]])
    }
  }
  
  # HEADER
  if(is.null(header)){
    # 3RD VARIABLE AS PLOT HEADER
    if(length(group_by) == 3){
      header <- rep(levels(pd[,group_by[3]]), 
                    each = ceiling(nlevels(pd[, group_by[3]]) / sheets))
    }
  }else if(length(header) != sheets){
    header <- rep(header, times = sheets)
  }  
  
  # RESET SAVED LABEL CO-ORDINATES ---------------------------------------------
  
  # PREVIOUS CALL
  previous_call <- getOption("cyto_plot_call")
  
  # CURRENT CALL (gh does not match)
  current_call <- list("x" = fs,
                       "channels" = channels,
                       "overlay" = overlay,
                       "group_by" = group_by,
                       "limits" = limits,
                       "gate" = gate,
                       "negate" = negate,
                       "label" = label,
                       "label_text" = label_text,
                       "label_stat" = label_stat,
                       "label_position" = label_position,
                       "label_text_x" = label_text_x,
                       "label_text_y" = label_text_y,
                       "density_modal" = density_modal,
                       "layout" = layout)
  
  # PREVIOUS CALL BELONGS TO DIFFERENT METHOD
  if(!all(names(previous_call) %in% names(current_call))){
    previous_call <- NULL
  }
  
  # UPDATE CYTO_PLOT_CALL
  options("cyto_plot_call" = current_call)
  
  # RESET SAVED LABEL CO-ORDINATES - MATCHING CALLS / NO CYTO_PLOT_SAVE
  if (!isTRUE(all.equal(previous_call, current_call)) |
      getOption("cyto_plot_save") == FALSE) {
    # RESET SAVED LABEL CO-ORDINATES
    options("cyto_plot_label_coords" = NULL)
    # RESET CYTO_PLOT_MATCH
    options("cyto_plot_match" = NULL)
  }
  
  # REPEAT CYTO_PLOT ARGUMENTS -------------------------------------------------
  
  # PULL DOWN ARGUMENTS
  args <- .args_list()
  
  # REPEAT & SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args)
  
  # PLOTTING ARGUMENTS
  cyto_plot_args <- formalArgs("cyto_plot.flowFrame")
  
  # REMOVE UNWANTED ARGUMENTS
  cyto_plot_args <- cyto_plot_args[-which(
    LAPPLY(c("x",
             "channels",
             "overlay",
             "axes_trans",
             "popup",
             "limits",
             "title",
             "axes_label"), function(z){
               grepl(z, cyto_plot_args)
             }))]
  
  # RESTRICT ARGS TO CYTO_PLOT_ARGS ONLY
  args <- args[names(args) %in% cyto_plot_args]
  
  # SETUP PLOTTING SPACE -------------------------------------------------------
  
  # POPUP
  cyto_plot_new(popup = popup)
  
  # SET LAYOUT
  par(mfrow = layout)
  
  # REMOVE MARGINAL SPACE AROUND EACH PLOT
  par(mar = c(0,0,0,0))
  
  # SPACE FOR LABELS, TITLES & HEADERS
  if(is.null(header)){
    par(oma = c(6,6,5,6)) # used to c(6,6,5,5)
  }else{
    par(oma = c(6,6,8,6)) # used to be c(6,6,8,5)
  }
  
  # PLOT CONTRUCTION -----------------------------------------------------------
  
  # PLOT PER expected_groups ELEMENT
  cnt <- 0
  sht <- 0
  plots <- lapply(seq_along(expected_groups), function(z){
    
    print(expected_groups[[z]])
    
    # COUNTER
    cnt <<- cnt + 1
    
    # SHEET NUMBER
    if(z %in% seq_len(ns)*np){
      sht <<- sht + 1
    }
    
    # OVERLAY
    if(length(expected_groups[[z]]) > 1){
      x <- expected_groups[[z]][[1]]
      overlay <- expected_groups[[z]][-1]
    }else{
      x <- expected_groups[[z]][[1]]
      overlay <- NA
    }
    
    # ARGUMENTS
    cyto_plot_args <- lapply(args, function(y){
      if(is(y, "list")){
        y[[z]]
      }else{
        y[z]
      }
    })
    
    # CYTO_PLOT FLOWFRAME METHOD
    do.call("cyto_plot",
            c(cyto_plot_args,
              list("x" = x,
                   "channels" = channels,
                   "overlay" = overlay,
                   "axes_trans" = axes_trans,
                   "popup" = FALSE,
                   "limits" = limits,
                   "title" = NA)))

    # PLOT LIMITS
    lims <- par("usr")
    xmin <- lims[1]
    xmax <- lims[2]
    xrng <- xmax - xmin
    ymin <- lims[3]
    ymax <- lims[4]
    
    # PANEL X LABELS
    if(panel_label[1] & !is.na(panel_x_label)){
      mtext(panel_x_label[z/sht],
            outer = FALSE,
            side = 3,
            line = 0.5,
            font = panel_label_font[1],
            cex = panel_label_size[1],
            col = panel_label_col[1])
    }
    
    # PANEL Y LABELS
    if(panel_label[2] & !is.na(panel_y_label)){
      # mtext does not support text rotation - need to use text for y labels
      text(x = xmax + 0.07*xrng,
           y = mean(c(ymin, ymax)),
           labels = panel_y_label[z/(ncol*sht)],
           font = panel_label_font[2],
           cex = 1.5*panel_label_size[2],
           col = panel_label_col[2],
           srt = 270,
           xpd = NA)
    }
    
    # AXES LABELS, PANEL LABELS, TITLES and HEADERS - FULL SHEET
    if(cnt %% prod(layout) == 0){
      
      # X AXIS LABEL
      mtext(xlab, 
            outer = TRUE, 
            side = 1, 
            line = 3.5,
            font = axes_label_text_font,
            cex = axes_label_text_size,
            col = axes_label_text_col)
      
      # Y AXIS LABEL
      mtext(ylab, 
            outer = TRUE, 
            side = 2, 
            line = 3.5,
            font = axes_label_text_font,
            cex = axes_label_text_size,
            col = axes_label_text_col)
      
      # X TITLE
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
      
      # Y TITLE
      if(is.na(title[2]) & all(group_by != FALSE)){
        title[2] <- group_by[2]
      }
      
      if(!is.na(title[2])){
        # mtext does not support text rotation - need to use text for y labels
        text(x = xmax + 0.2*xrng,
             y = 2*ymax, # HALF GRID SIZE
             labels = title[2],
             font = title_text_font[2],
             cex = 1.8*title_text_size[2],
             col = title_text_col[2],
             srt = 270,
             xpd = NA)
      }
      
      # HEADER
      if(!is.null(header)){
        mtext(header[sht], 
              outer = TRUE, 
              side = 3, 
              line = 5.5,
              font = header_text_font,
              cex = header_text_size,
              col = header_text_col)
      }
      
      # RECORD PLOT
      plot <- cyto_plot_record()
      
      # POPUP
      if(popup & cnt != ns*np){
        cyto_plot_new(popup)
      }
      
      # RETURN RECORDED PLOT
      return(plot)
      
    }
    
  })
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE
  if(getOption("cyto_plot_save")){
    if(inherits(x, basename(getOption("cyto_plot_method")))){
      
      # CLOSE GRAPHICS DEVICE
      dev.off()
      
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN RECORDED PLOT
  invisible(plots)
  
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
