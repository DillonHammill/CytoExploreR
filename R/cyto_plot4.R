# CYTO_PLOT --------------------------------------------------------------------

#' cyto_plot4
#'
#' @importFrom graphics par
#' @importFrom magrittr %>%
#' @importFrom purrr transpose
#'
#' @name cyto_plot4
NULL

#' @noRd
#' @export
cyto_plot4 <- function(x, ...) {
  UseMethod("cyto_plot4")
}

#' @rdname cyto_plot4
#' @export
cyto_plot4.GatingSet <- function(){
  
}

#' @rdname cyto_plot4
#' @export
cyto_plot4.GatingHierarchy <- function(x,
                                       parent,
                                       alias = NA,
                                       channels,
                                       axes_trans = NA,
                                       overlay = NA,
                                       gate = NA,
                                       limits = "data",
                                       display = 25000,
                                       popup = FALSE,
                                       xlim = NA,
                                       ylim = NA,
                                       xlab,
                                       ylab,
                                       title,
                                       negate,
                                       density_modal = TRUE,
                                       density_smooth = 0.6,
                                       density_stack = 0,
                                       density_cols = NA,
                                       density_fill = NA,
                                       density_fill_alpha = 1,
                                       density_line_type = 1,
                                       density_line_width = 1,
                                       density_line_col = "black",
                                       point_shape = ".",
                                       point_size = 2,
                                       point_col_scale = NA,
                                       point_cols = NA,
                                       point_col = NA,
                                       point_col_alpha = 1,
                                       contour_lines = 0,
                                       contour_line_type = 1,
                                       contour_line_width = 1,
                                       contour_line_col = "black",
                                       contour_line_alpha = 1,
                                       axes_text = c(TRUE, TRUE),
                                       axes_text_font = 1,
                                       axes_text_size = 1,
                                       axes_text_col = "black",
                                       axes_label_text_font = 1,
                                       axes_label_text_size = 1.1,
                                       axes_label_text_col = "black",
                                       title_text_font = 2,
                                       title_text_size = 1.1,
                                       title_text_col = "black",
                                       legend = FALSE,
                                       legend_text = NA,
                                       legend_text_font = 1,
                                       legend_text_size = 1,
                                       legend_text_col = "black",
                                       legend_line_type = NA,
                                       legend_line_width = NA,
                                       legend_line_col = NA,
                                       legend_box_fill = NA,
                                       legend_point_col = NA,
                                       gate_line_type = 1,
                                       gate_line_width = 2.5,
                                       gate_line_col = "red",
                                       gate_fill = "white",
                                       gate_fill_alpha = 0,
                                       label,
                                       label_text,
                                       label_stat,
                                       label_position = "auto",
                                       label_text_x = NA,
                                       label_text_y = NA,
                                       label_text_font = 2,
                                       label_text_size = 1,
                                       label_text_col = "black",
                                       label_fill = "white",
                                       label_fill_alpha = 0.6,
                                       border_line_type = 1,
                                       border_line_width = 1,
                                       border_line_col = "black",
                                       border_fill = "white",
                                       border_fill_alpha = 1, ...){
  
  # GATINGHIERARCHY METHOD - CALLS FLOWFRAME METHOD
  
  # CHECKS ---------------------------------------------------------------------
  
  # PARENT
  if(missing(parent)){
    stop("Supply the name of the 'parent' population to plot.")
  }
  
  # CHANNELS
  if(missing(channels)){
    stop("Supply the channel/marker(s) to construct the plot.")
  }
  
  # GATINGHIERACHY - gh (x available for flowframe method call)
  gh <- x
  
  # TRANSFORMATIONS
  if(.all_na(axes_trans)){
    # TRANSFORMERLIST FROM GATINGHIERARCHY
    axes_trans <- gh@transformation
    if(length(axes_trans) == 0){
      axes_trans <- NA
    }
  }
  
  # PREPARE DATA & ARGUMENTS ---------------------------------------------------
  
  # EXTRACT PARENT POPULATION
  x <- cyto_extract(gh, parent)
  
  # OVERLAY - POPULATION NAMES
  if(!.all_na(overlay)){
    # POPULATION NAMES TO OVERLAY
    if(is.character(overlay)){
      if(all(overlay %in% basename(cyto_nodes(gh)))){
        # ADD OVERLAY TO LABEL_TEXT
        if(missing(label_text)){
          if(!.all_na(alias)){
            if(length(channels) == 1){
              label_text <- c(alias, overlay)
            }else if(length(channels) == 2){
              label_text <- alias
            }
          }
        }
        # EXTRACT POPULATIONS
        nms <- overlay
        overlay <- lapply(overlay, function(z){
          cyto_extract(gh, z)
        })
        names(overlay) <- nms
      }
    # OVERLAY - POPULATIONS
    }else{
      # ADD OVERLAY TO LABEL_TEXT
      if(missing(label_text)){
        if(!.all_na(alias)){
          label_text <- rep(alias, length.out = length(overlay) + 1)
        }
      }
    }
  # NO OVERLAY
  }else if(.all_na(overlay)){
    # LABEL TEXT
    if(missing(label_text)){
      if(!.all_na(alias)){
        label_text <- NA
      }else{
        label_text <- alias
      }
    }
  }
  
  # PREPARE GATES --------------------------------------------------------------
  
  # EMPTY ALIAS
  if(.empty(alias)){
    # ALL GATES IN SUPPLIED CHANNELS
    if(all(alias == "")){
      gt <- templateGen(gh)
      gt <- gt[basename(gt$parent) == parent, ]
      # 2D PLOT - BOTH CHANNELS MATCH
      if(length(channels) == 2){
        alias <- gt$alias[gt$dims == paste(channels, collapse = ",")]
      # 1D PLOT - ONE CHANNEL MATCH
      }else if(length(channels) == 1){
        ind <- lapply(gt$dims, function(z){
          grep(channels, z)
        })
        ind <- LAPPLY(ind, "length") != 0
        alias <- gt$alias[ind]
      }
      # NO ALIAS IN SUPPLIED CHANNELS
      if(length(alias) == 0){
        alias <- NA
      }
    }
  }
  
  # EXTRACT GATE OBJECTS
  if(!.all_na(alias)){
    # LIST OF GATE OBJECTS
    gate <- lapply(alias, function(z){
      getGate(gh, paste(parent, z, sep = "/"))
    })
    # WATCH OUT FOR POPULATIONS SHARING THE SAME GATE
    gate <- unique(gate)
  }
  
  # SUPPORT NEGATED BOOLEAN FILTERS
  
  # SUPPORT NEGATED GATES
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # TITLE
  if(missing(title)){
    # SAMPLENAME
    title <- cyto_names(x)
    # PARENT
    title <- LAPPLY(title, function(z){
      if(parent == "root"){
        pt <- "All Events"
      }else{
        pt <- parent
      }
      # 1D PLOT - STACKED NO OVERLAY - LACK SAMPLENAMES
      if(length(channels) == 1 &
         .all_na(overlay) &
         density_stack != 0){
        pt
      # 1D PLOT - STACKED OVERLAY - SAMPLENAMES ONLY
      }else if(length(channels) == 1 &
               !.all_na(overlay) &
               density_stack != 0){
        z
      # PASTE SAMPLNAME & PARENT  
      }else{
        paste(z, pt, sep = "\n")
      }
      # Stacked 1D plots lack sampleNames in titles if no overlay
      if (length(channels) == 1 &
          .all_na(overlay) &
          density_stack != 0) {
        pt
        # Stacked with overlay display sampleNames only
      } else if (length(channels) == 1 &
                 !.all_na(overlay) &
                 density_stack != 0) {
        z
        # Paste together sampleName and parent name
      } else {
        paste(z, "\n", pt, sep = " ")
      }
    })
  }
  
  # NEGATE
  if(missing(negate)){
    negate <- FALSE
  }
  
  # CALL CYTO_PLOT FLOWFRAME METHOD --------------------------------------------
  
  # PULL DOWN ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT FLOWFRAME ARGUMENTS
  ARGS <- formalArgs("cyto_plot4.flowFrame")
  
  # RESTRCT ARGUMENTS
  args <- args[ARGS]
  
  # CALL FLOWFRAME METHOD
  do.call("cyto_plot4.flowFrame", args)
  
}

#' @rdname cyto_plot4
#' @export
cyto_plot4.flowSet <- function(x,
                               channels,
                               axes_trans = NA,
                               group_by = NA,
                               overlay = NA,
                               gate = NA,
                               limits = "data",
                               display = 1,
                               layout,
                               popup = FALSE,
                               xlim = NA,
                               ylim = NA,
                               xlab,
                               ylab,
                               title,
                               negate = FALSE,
                               density_modal = TRUE,
                               density_smooth = 0.6,
                               density_stack = 0,
                               density_layers = NA,
                               density_cols = NA,
                               density_fill = NA,
                               density_fill_alpha = 1,
                               density_line_type = 1,
                               density_line_width = 1,
                               density_line_col = "black",
                               point_shape = ".",
                               point_size = 2,
                               point_col_scale = NA,
                               point_cols = NA,
                               point_col = NA,
                               point_col_alpha = 1,
                               contour_lines = 0,
                               contour_line_type = 1,
                               contour_line_width = 1,
                               contour_line_col = "black",
                               contour_line_alpha = 1,
                               axes_text = c(TRUE, TRUE),
                               axes_text_font = 1,
                               axes_text_size = 1,
                               axes_text_col = "black",
                               axes_label_text_font = 1,
                               axes_label_text_size = 1.1,
                               axes_label_text_col = "black",
                               title_text_font = 2,
                               title_text_size = 1.1,
                               title_text_col = "black",
                               legend = FALSE,
                               legend_text = NA,
                               legend_text_font = 1,
                               legend_text_size = 1,
                               legend_text_col = "black",
                               legend_line_type = NA,
                               legend_line_width = NA,
                               legend_line_col = NA,
                               legend_box_fill = NA,
                               legend_point_col = NA,
                               gate_line_type = 1,
                               gate_line_width = 2.5,
                               gate_line_col = "red",
                               gate_fill = "white",
                               gate_fill_alpha = 0,
                               label,
                               label_text = NA,
                               label_stat,
                               label_position = "auto",
                               label_text_x = NA,
                               label_text_y = NA,
                               label_text_font = 2,
                               label_text_size = 0.8,
                               label_text_col = "black",
                               label_fill = "white",
                               label_fill_alpha = 0.6,
                               border_line_type = 1,
                               border_line_width = 1,
                               border_line_col = "black",
                               border_fill = "white",
                               border_fill_alpha = 1, ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # METHOD
  if(is.null(getOption("cyto_plot_method"))){
    options("cyto_plot_method" = "flowSet")
  }
  
  # CUSTOM PLOT - LAYOUT FALSE
  if(!missing(layout)){
    if(all(layout == FALSE)){
      options("cyto_plot_custom" = TRUE)
    }
  }
  
  # CHANNELS
  if(missing(channels)){
    stop("Supply channel/marker(s) to construct the plot.")
  }else{
    channels <- cyto_channels_extract(x, channels = channels, plot = TRUE)
  }
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- par(c("mar", "mfrow"))
  
  # RESET ON EXIT
  on.exit(par(old_pars))
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS - MISSING -> EMPTY ""
  .args_update(args)
  
  # CYTO_PLOT_SAVE -------------------------------------------------------------
  
  # POPUP
  if(getOption("cyto_plot_save") == TRUE){
    popup <- FALSE
  }
  
  # SAMPLE PREPARATION - BASE LAYERS & OVERLAYS --------------------------------
  
  # DATA TO LIST & GROUP - CONVERT GROUPS TO FLOWFRAMES
  if(!.all_na(group_by)){
    fs_list <- cyto_group_by(x, group_by)
    fr_list <- lapply(fs_list, function(z){
      cyto_convert(z, "flowFrame")
    })
  }else{
    fr_list <- cyto_convert(x, "list of flowFrames")
    names(fr_list) <- cyto_names(x)
  }
  
  # OVERLAY LIST & GROUP
  if(!.all_na(overlay)){
    # REPEAT FLOWFRAME PER PLOT
    if(inherits(overlay, "flowFrame")){
      overlay_list <- rep(list(list(overlay)), length(fr_list))
    # FLOWSET TO LIST OF FLOWFRAMES
    }else if(inherits(overlay, "flowSet")){
      # GROUPING
      if(!.all_na(group_by)){
        overlay_list <- cyto_group_by(overlay, group_by)
        overlay_list <- lapply(overlay_list, function(z){
          cyto_convert(z, "flowFrame")
        })
      # NO GROUPING
      }else{
        overlay_list <- cyto_convert(overlay, "list of flowFrames")
      }
      # LIST OF FLOWFRAME LISTS
      overlay_list <- lapply(overlay_list, function(z){
        list(z)
      })
    # LIST OF FLOWSETS TO LIST OF FLOWFRAME LISTS
    }else if(inherits(overlay, "list")){
      # ALLOW LIST OF FLOWFRAMES OF LENGTH FR_LIST
      if(all(LAPPLY(unlist(overlay), function(z){
        inherits(z, "flowFrame")
      }))) {
        # SAME LENGTH AS FR_LIST
        if(length(overlay) != length(fr_list)){
          stop(
            paste(
              "'overlay' must be  alist of flowFrame lists -",
              "one flowFrame list per plot."
            ))
        }
        # NO GROUPING APPLIED
        overlay_list <- overlay
      }
    # LIST OF FLOWSETS
    }else if(all(LAPPLY(overlay, function(z){
      inherits(z, "flowSet")
    }))){
      # GROUPING
      if(!.all_na(group_by)){
        overlay_list <- lapply(overlay, function(z){
          cyto_group_by(z, group_by)
        })
        overlay_list <- lapply(overlay_list, function(z){
          lapply(z, function(y){
            cyto_convert(y, "flowFrame")
          })
        })
      # NO GROUPING
      }else{
        overlay_list <- lapply(overlay, function(z){
          cyto_convert(z, "list of flowFrames")
        })
      }
      overlay_list <- overlay_list %>% transpose()
    # OVERLAY NOT SUPPORTED
    }else{
      stop(
        paste(
          "'overlay' should be either a flowFrame, flowSet, list of flowFrame",
          "lists or list of flwoSet lists."
        )
      )
    }
  # NO OVERLAY
  }else if(.all_na(overlay)){
    # SAMPLENAMES
    NMS <- names(fr_list)
    # LIST OF FLOWFRAME LISTS
    fr_list <- lapply(seq_len(length(fr_list)), function(z){
      lst <- list(fr_list[[z]])
      names(lst) <- NMS[z]
      return(lst)
    })
  }
  
  # COMBINE BASE LAYERS WITH OVERLAY - LIST OF FLOWFRAME LISTS -----------------
  
  # OVERLAY
  if(!.all_na(overlay)){
    NMS <- names(fr_list)
    fr_list <- lapply(seq_len(length(fr_list)), function(z){
      c(fr_list[z], overlay_list[[z]])
    })
    names(fr_list) <- NMS
  }
  
  # DENSITY_LAYERS -------------------------------------------------------------
  
  # SUPPORT SPLIT STACKED SAMPLES (NO OVERLAY)
  if(length(channels) == 1 & .all_na(overlay)){
    # CONVERT LIST OF INDIVIDUAL FLOWFRAME LISTS TO LIST(LIST OF FLOWFRAMES)
    fr_list <- list(unlist(fr_list))
    # LEGEND_TEXT - SAMPLENAMES
    if(.all_na(legend_text)){
      legend_text <- NMS
    }
    # DENSITY_LAYERS
    if(!.all_na(density_layers)){
      # SAME # LAYERS PER PLOT
      if(length(fr_list[[1]]) %% density_layers != 0){
        stop("Each plot must have the same number of layers!")
      }
      # INDICES
      ind <- rep(seq_len(length(fr_list[[1]])),
                 each = density_layers,
                 length.out = length(fr_list[[1]]))
      # SPLITTING BY DENSITY_LAYERS
      fr_list <- lapply(unique(ind), function(z){
        fr_list[[1]][ind == z]
      })
    }
  }
  
  # GATES PREPARATION ----------------------------------------------------------
  
  # LIST OF GATE OBJECT LISTS
  if(!.all_na(gate)){
    # REPEAT GATE OBJECTS PER PLOT
    if(any(is(gate) %in% c("rectangleGate",
                           "polygonGate",
                           "ellipsoidGate",
                           "quadGate",
                           "filters"))){
      gate <- list(gate)
    }
    # LIST OF GATE OBJECTS - REPEAT PER PLOT
    if(all(LAPPLY(gate, function(z){any(is(z) %in% c("rectangleGate",
                                                     "polygonGate",
                                                     "ellipsoidGate",
                                                     "quadGate",
                                                     "filters"))}))){
      gate <- rep(list(gate), length.out = length(fr_list))
    }
  }
  
  # GATES STACKED DENSITY - NO OVERLAY
  if(length(channels) == 1 & .all_na(overlay) & density_stack != 0){
    # USE FIRST SET OF GATES
    gate <- list(gate[[1]])
  }
  
  # ARGUMENT PREPARATION -------------------------------------------------------
  
  # XLIM
  if(.all_na(xlim)){
    xlim <- .cyto_range(fr_list,
                        channels = channels[1],
                        limits = limits)[, channels[1]]
  # XLIM MANUALLY SUPPLIED
  }else{
    if(!.all_na(axes_trans)){
      if(channels[1] %in% names(axes_trans)){
        xlim <- axes_trans[[channels[1]]]$transform(xlim)
      }
    }
  }
  
  # YLIM - 1D CALCULATED LATER
  if(.all_na(ylim)){
    # 2D PLOT
    if(length(channels) == 2){
      ylim <- .cyto_range(fr_list,
                          channels = channels[2],
                          limits = limits)[, channels[2]]
    }
  # YLIM MANUALLY SUPPLIED
  }else{
    # 2D PLOT
    if(length(channels) == 2){
      if(!.all_na(axes_trans)){
        if(channels[2] %in% names(axes_trans)){
          ylim <- axes_trans[[channels[2]]]$transform(ylim)
        }
      }
    }
  }
  
  # X AXIS BREAKS & LABELS
  if(axes_text[1] == TRUE){
    axes_text_x <- .cyto_plot_axes_text(x[[1]],
                                        channels = channels[1],
                                        axes_trans = axes_trans,
                                        axes_range = list(xlim, ylim),
                                        limits = limits)[[1]]
  }else{
    axes_text_x <- FALSE
  }
  
  # Y AXIS BREAKS & LABELS
  if(axes_text[2] == TRUE){
    if(length(channels) == 2){
      axes_text_y <- .cyto_plot_axes_text(x[[1]],
                                          channels = channels[2],
                                          axes_trans = axes_trans,
                                          axes_range = list(xlim, ylim),
                                          limits = limits)[[1]]
    }else if(length(channels) == 1){
      axes_text_y <- NA
    }
  }else{
    axes_text_y <- FALSE
  }
  
  # AXES_TEXT
  axes_text <- list(axes_text_x, axes_text_y)
  axes_text <- rep(axes_text, length(fr_list))
  
  # GRAPHICS DEVICE
  if(popup == TRUE){
    cyto_plot_new(popup)
  }
  
  # LAYOUT
  if(.empty(layout)){
    # LAYOUT DIMENSIONS
    layout <- .cyto_plot_layout(fr_list,
                                layout = layout,
                                density_stack = density_stack,
                                density_layers = density_layers)
  }else if(all(layout == FALSE) | .all_na(layout)){
    # CURRENT DIMENSIONS
    layout <- par("mfrow")
  }
  
  # SET LAYOUT
  par("mfrow" = layout)
  np <- layout[1] * layout[2]
  
  # CYTO_PLOT_CALL -------------------------------------------------------------
  
  # PREVIOUS CALL
  previous_call <- getOption("cyto_plot_call")
  
  
  # CURRENT ARGUMENTS
  args <- .args_list()
  
  # CURRENT CALL - ARGUMENTS INFLUENCING LABEL LOCATIONS
  current_call <- args[c("x",
                         "channels",
                         "overlay",
                         "group_by",
                         "limits",
                         "gate",
                         "negate",
                         "label",
                         "label_text",
                         "label_stat",
                         "label_position",
                         "label_text_x",
                         "label_text_y",
                         "density_modal",
                         "density_stack",
                         "density_layers")]
  
  # UPDATE CYTO_PLOT_CALL
  options("cyto_plot_call" = current_call)
  
  # RESET SAVED LABEL CO-ORDINATES - MATCHING CALLS / NO CYTO_PLOT_SAVE
  if(!isTRUE(all.equal(previous_call, current_call)) |
     getOption("cyto_plot_save") == FALSE){
    # RESET SAVED LABEL CO-ORDINATES
    options("cyto_plot_label_coords" = NULL)
    # RESET CYTO_PLOT_MATCH
    options("cyto_plot_match" = NULL)
  }

  # REPEAT & SPLIT ARGUMENTS ---------------------------------------------------
  
  # REPEAT & SPLIT ARGUMENTS
  args <- .cyto_plot_args_split(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # CALL CYTO_PLOT FLOWFRAME METHOD --------------------------------------------
  
  # PASS ARGUMENTS TO CYTO_PLOT FLOWFRAME METHOD
  cnt <- 0
  mapply(
    function(x,
             gate,
             limits,
             display,
             xlab,
             ylab,
             title,
             title_text_font,
             title_text_size,
             title_text_col,
             density_modal,
             density_smooth,
             density_stack,
             density_fill,
             density_fill_alpha,
             density_line_type,
             density_line_width,
             density_line_col,
             point_shape,
             point_size,
             point_col,
             point_col_alpha,
             contour_lines,
             contour_line_type,
             contour_line_width,
             contour_line_col,
             contour_line_alpha,
             axes_text,
             axes_text_font,
             axes_text_size,
             axes_text_col,
             axes_label_text_font,
             axes_label_text_size,
             axes_label_text_col,
             legend,
             legend_text,
             legend_text_font,
             legend_text_size,
             legend_text_col,
             legend_line_type,
             legend_line_width,
             legend_line_col,
             legend_box_fill,
             legend_point_col,
             gate_line_type,
             gate_line_width,
             gate_line_col,
             gate_fill,
             gate_fill_alpha,
             label,
             label_text,
             label_stat,
             label_position,
             label_text_x,
             label_text_y,
             label_text_font,
             label_text_size,
             label_text_col,
             label_fill,
             label_fill_alpha,
             border_line_type,
             border_line_width,
             border_line_col,
             border_fill,
             border_fill_alpha){
      
      # PLOT COUNTER
      cnt <<- cnt + 1
      
      # OVERLAY
      if(length(x) > 1){
        overlay <- x[seq(2, length(x), 1)]
      }else{
        overlay <- NA
      }
      
      # CALL CYTO_PLOT FLOWFRAME METHOD
      cyto_plot4(x[[1]],
                 channels = channels,
                 overlay = overlay,
                 gate = gate,
                 axes_trans = axes_trans,
                 limits = limits,
                 display = display,
                 popup = FALSE,
                 xlim = xlim,
                 ylim = ylim,
                 xlab = xlab,
                 ylab = ylab,
                 title = title,
                 negate = negate,
                 title_text_font = title_text_font,
                 title_text_size = title_text_size,
                 title_text_col = title_text_col,
                 density_modal = density_modal,
                 density_smooth = density_smooth,
                 density_stack = density_stack,
                 density_cols = density_cols,
                 density_fill = density_fill,
                 density_fill_alpha = density_fill_alpha,
                 density_line_type = density_line_type,
                 density_line_width = density_line_width,
                 density_line_col = density_line_col,
                 point_shape = point_shape,
                 point_size = point_size,
                 point_col_scale = point_col_scale,
                 point_cols = point_cols,
                 point_col = point_col,
                 point_col_alpha = point_col_alpha,
                 contour_lines = contour_lines,
                 contour_line_type = contour_line_type,
                 contour_line_width = contour_line_width,
                 contour_line_col = contour_line_col,
                 contour_line_alpha = contour_line_alpha,
                 axes_text = axes_text,
                 axes_text_font = axes_text_font,
                 axes_text_size = axes_text_size,
                 axes_text_col = axes_text_col,
                 axes_label_text_font = axes_label_text_font,
                 axes_label_text_size = axes_label_text_size,
                 axes_label_text_col = axes_label_text_col,
                 legend = legend,
                 legend_text = legend_text,
                 legend_text_font = legend_text_font,
                 legend_text_size = legend_text_size,
                 legend_text_col = legend_text_col,
                 legend_line_type = legend_line_type,
                 legend_line_width = legend_line_width,
                 legend_line_col = legend_line_col,
                 legend_box_fill = legend_box_fill,
                 legend_point_col = legend_point_col,
                 gate_line_type = gate_line_type,
                 gate_line_width = gate_line_width,
                 gate_line_col = gate_line_col,
                 gate_fill = gate_fill,
                 gate_fill_alpha = gate_fill_alpha,
                 label = label,
                 label_text = label_text,
                 label_stat = label_stat,
                 label_position = label_position,
                 label_text_x = label_text_x,
                 label_text_y = label_text_y,
                 label_text_font = label_text_font,
                 label_text_size = label_text_size,
                 label_text_col = label_text_col,
                 label_fill = label_fill,
                 label_fill_alpha = label_fill_alpha,
                 border_line_type = border_line_type,
                 border_line_width = border_line_width,
                 border_line_col = border_line_col,
                 border_fill = border_fill,
                 border_fill_alpha = border_fill_alpha)
      
      # NEW PLOT PAGE
      if(popup == TRUE &
         cnt %% np == 0 &
         length(fr_list) > cnt){
        cyto_plot_new(popup = popup)
        par("mfrow" = layout)
      }
      
    },
    fr_list,
    gate,
    limits,
    display,
    xlab,
    ylab,
    title,
    title_text_font,
    title_text_size,
    title_text_col,
    density_modal,
    density_smooth,
    density_stack,
    density_fill,
    density_fill_alpha,
    density_line_type,
    density_line_width,
    density_line_col,
    point_shape,
    point_size,
    point_col,
    point_col_alpha,
    contour_lines,
    contour_line_type,
    contour_line_width,
    contour_line_col,
    contour_line_alpha,
    axes_text,
    axes_text_font,
    axes_text_size,
    axes_text_col,
    axes_label_text_font,
    axes_label_text_size,
    axes_label_text_col,
    legend,
    legend_text,
    legend_text_font,
    legend_text_size,
    legend_text_col,
    legend_line_type,
    legend_line_width,
    legend_line_col,
    legend_box_fill,
    legend_point_col,
    gate_line_type,
    gate_line_width,
    gate_line_col,
    gate_fill,
    gate_fill_alpha,
    label,
    label_text,
    label_stat,
    label_position,
    label_text_x,
    label_text_y,
    label_text_font,
    label_text_size,
    label_text_col,
    label_fill,
    label_fill_alpha,
    border_line_type,
    border_line_width,
    border_line_col,
    border_fill,
    border_fill_alpha)
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save") == TRUE) {
    if (inherits(x, getOption("cyto_plot_method"))) {
      if (!getOption("cyto_plot_custom")) {
        # CLOSE GRAPHICS DEVICE
        dev.off()
      }
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # REST CYTO_PLOT_METHOD
  if(getOption("cyto_plot_method") == "flowSet"){
    options("cyto_plot_method" = NULL)
  }
  
  # RETURN RECORDED PLOT
  invisible(recordPlot())
  
}

#' @rdname cyto_plot4
#' @export
cyto_plot4.flowFrame <- function(x,
                                 channels,
                                 axes_trans = NA,
                                 overlay = NA,
                                 gate = NA,
                                 limits = "data",
                                 display = 25000,
                                 popup = FALSE,
                                 xlim = NA,
                                 ylim = NA,
                                 xlab,
                                 ylab,
                                 title,
                                 negate = FALSE,
                                 density_modal = TRUE,
                                 density_smooth = 0.6,
                                 density_stack = 0,
                                 density_cols = NA,
                                 density_fill = NA,
                                 density_fill_alpha = 1,
                                 density_line_type = 1,
                                 density_line_width = 1,
                                 density_line_col = "black",
                                 point_shape = ".",
                                 point_size = 2,
                                 point_col_scale = NA,
                                 point_cols = NA,
                                 point_col = NA,
                                 point_col_alpha = 1,
                                 contour_lines = 0,
                                 contour_line_type = 1,
                                 contour_line_width = 1,
                                 contour_line_col = "black",
                                 contour_line_alpha = 1,
                                 axes_text = c(TRUE, TRUE),
                                 axes_text_font = 1,
                                 axes_text_size = 1,
                                 axes_text_col = "black",
                                 axes_label_text_font = 1,
                                 axes_label_text_size = 1.1,
                                 axes_label_text_col = "black",
                                 title_text_font = 2,
                                 title_text_size = 1.1,
                                 title_text_col = "black",
                                 legend = FALSE,
                                 legend_text = NA,
                                 legend_text_font = 1,
                                 legend_text_size = 1,
                                 legend_text_col = "black",
                                 legend_line_type = NA,
                                 legend_line_width = NA,
                                 legend_line_col = NA,
                                 legend_box_fill = NA,
                                 legend_point_col = NA,
                                 gate_line_type = 1,
                                 gate_line_width = 2.5,
                                 gate_line_col = "red",
                                 gate_fill = "white",
                                 gate_fill_alpha = 0,
                                 label,
                                 label_text,
                                 label_stat,
                                 label_position = "auto",
                                 label_text_x = NA,
                                 label_text_y = NA,
                                 label_text_font = 2,
                                 label_text_size = 1,
                                 label_text_col = "black",
                                 label_fill = "white",
                                 label_fill_alpha = 0.6,
                                 border_line_type = 1,
                                 border_line_width = 1,
                                 border_line_col = "black",
                                 border_fill = "white",
                                 border_fill_alpha = 1, ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "flowFrame")
  }

  # CHANNELS
  if (missing(channels)) {
    stop("Supply channel/marker(s) to construct the plot.")
  } else {
    channels <- cyto_channels_extract(x, channels, plot = TRUE)
  }

  # AXES_TRANS
  if (!.all_na(axes_trans)) {
    if (inherits(axes_trans, "transformList")) {
      axes_trans <- NA
      message("Supply a transformerList object to axes_trans to transform axes.")
    }
  }

  # CYTO_PLOT_THEME ------------------------------------------------------------

  # ARGUMENTS
  args <- .args_list()

  # THEME
  args <- .cyto_plot_theme_inherit(args)

  # UPDATE ARGUMENTS - MISSING -> EMPTY ""
  .args_update(args)

  # CYTO_PLOT_SAVE -------------------------------------------------------------

  # POPUP
  if (getOption("cyto_plot_save") == TRUE) {
    popup <- FALSE
  }

  # SAMPLE PREPARATION - LIST OF FLOWFRAMES ------------------------------------

  # Add flowFrame to list
  fr_list <- list(x)
  names(fr_list) <- cyto_names(x)

  # Add overlay to list and group if necessary
  if (!.all_na(overlay)) {
    # overlay must be list of flowFrames
    # flowFrame overlay added to list
    if (inherits(overlay, "flowFrame")) {
      overlay_list <- list(overlay)
      # flowSet overlay convert to list of flowFrames
    } else if (inherits(overlay, "flowSet")) {
      overlay_list <- cyto_convert(overlay, "list of flowFrames")
      # flowFrame list overlay as is - flowSet list overlay use overlay[[1]]
    } else if (inherits(overlay, "list")) {
      # overlay should be list of flowFrames
      if (all(LAPPLY(overlay, function(z) {
        inherits(z, "flowFrame")
      }))) {
        overlay_list <- overlay
        # overlay list of flowSets - use first fs convert to list of flowFrames
      } else if (all(LAPPLY(overlay, function(z) {
        inherits(z, "flowSet")
      }))) {
        overlay <- overlay[[1]]
        overlay_list <- cyto_convert(overlay, "list of flowFrames")
        # overlay not supported
      } else {
        stop(paste(
          "'overlay' should be either the names of the populations to",
          "overlay, a flowFrame, a flowSet or a list of flowFrames."
        ))
      }
    }
  }

  # Combine base layers with overlay into list of flowFrames
  if (!.all_na(overlay)) {
    fr_list <- c(fr_list, overlay_list)
  }

  # SAMPLES
  SMP <- length(fr_list)

  # OVERLAYS
  OVN <- SMP - 1
  
  # SAMPLE PREPARATION ---------------------------------------------------------

  # SAMPLING - SET SEED
  if (display != 1) {
    fr_list <- cyto_sample(fr_list, display = display, seed = 56)
  }

  # REMOVAL NEGATIVE FSC/SSC EVENTS - POINT_COL SCALE ISSUE
  lapply(seq_len(length(channels)), function(z) {
    if (grepl("FSC", channels[z], ignore.case = TRUE) |
      grepl("SSC", channels[z], ignore.case = TRUE)) {
      fr_list <<- lapply(fr_list, function(y) {
        if (min(range(y, type = "data")[, channels[z]]) < 0) {
          coords <- matrix(c(0, Inf), ncol = 1, nrow = 2)
          rownames(coords) <- c("min", "max")
          colnames(coords) <- channels[z]
          nonDebris <- rectangleGate(.gate = coords)
          Subset(y, nonDebris)
        } else {
          return(y)
        }
      })
    }
  })

  # DENSITY PREPARATION --------------------------------------------------------

  # DENSITY
  if (length(channels) == 1) {
    # COMPUTE STACKED KERNEL DENSITY
    fr_dens_list <- suppressMessages(
      .cyto_density(fr_list,
        channel = channels,
        smooth = density_smooth,
        modal = density_modal,
        stack = density_stack
      )
    )
  } else {
    fr_dens_list <- NA
  }
  
  # GATES PREPARATION ----------------------------------------------------------

  # LIST OF GATE OBJECTS
  if (!.all_na(gate)) {
    gate <- cyto_gate_prepare(gate, channels)
  }

  # POPULATIONS ----------------------------------------------------------------

  # POPULATIONS PER LAYER
  NP <- .cyto_gate_count(gate, negate = negate)

  # TOTAL POPULATIONS
  TNP <- NP * SMP
  TNP_split <- split(seq_len(TNP), rep(seq_len(SMP), each = NP))

  # ARGUMENT PREPARATION -------------------------------------------------------

  # LABEL_TEXT
  if (all(LAPPLY(label_text, ".empty"))) {
    label_text <- rep(NA, TNP)
  } else {
    label_text <- rep(c(label_text, rep(NA, TNP)), length.out = TNP)
  }

  # LABEL_STAT
  # 1D PLOT NO STACK
  if (length(channels) == 1 & density_stack == 0) {
    # LABEL_STAT MISSING
    if (all(LAPPLY(label_stat, ".empty"))) {
      # GATE - FREQ STAT
      if (!.all_na(gate)) {
        # LABEL_STAT - BASE LAYER ONLY
        label_stat <- c(
          rep("freq", NP),
          rep(NA, TNP - NP)
        )
        # NO GATE - NO STAT
      } else {
        # LABEL_STAT REMOVED
        label_stat <- rep(NA, TNP)
      }
      # LABEL_STAT SUPPLIED - FILL WITH NA
    } else {
      # GATE - BASE LAYER ONLY
      if (!.all_na(gate)) {
        if (length(label_stat) == 1) {
          label_stat <- rep(label_stat, length.out = NP)
        }
        label_stat <- rep(c(
          label_stat,
          rep(NA, length.out = TNP)
        ),
        length.out = TNP
        )
        # NO GATE
      } else {
        # LABEL EACH LAYER
        label_stat <- rep(label_stat, length.out = TNP)
      }
    }
    # 1D PLOT STACK
  } else if (length(channels) == 1 & density_stack != 0) {
    # LABEL_STAT MISSING
    if (all(LAPPLY(label_stat, ".empty"))) {
      # GATE - FREQ STAT
      if (!.all_na(gate)) {
        # LABEL_STAT - ALL LAYERS
        label_stat <- rep("freq", length.out = TNP)
        # NO GATE
      } else {
        # LABEL_STAT REMOVED
        label_stat <- rep(NA, length.out = TNP)
      }
      # LABEL_STAT SUPPLIED - FILL WITH NA
    } else {
      label_stat <- rep(label_stat, length.out = TNP)
    }
    # 2D PLOT
  } else if (length(channels) == 2) {
    # LABEL_STAT MISSING
    if (all(LAPPLY(label_stat, ".empty"))) {
      # GATE - FREQ STAT
      if (!.all_na(gate)) {
        # LABEL_STAT - BASE LAYER ONLY
        label_stat <- c(
          rep("freq", NP),
          rep(NA, TNP - NP)
        )
        # NO GATE - NO STAT
      } else {
        # LABEL_STAT REMOVED
        label_stat <- rep(NA, length.out = TNP)
      }
      # LABEL_STAT SUPPLIED- FILL WITH NA
    } else {
      # GATE - BASE LAYER ONLY
      if (!.all_na(gate)) {
        if (length(label_stat) == 1) {
          label_stat <- rep(label_stat, length.out = NP)
        }
        label_stat <- rep(c(
          label_stat,
          rep(NA, length.out = TNP)
        ),
        length.out = TNP
        )
        # NO GATE
      } else {
        # LABEL EACH LAYER
        label_stat <- rep(label_stat, length.out = TNP)
      }
    }
  }

  # LABEL
  if (all(LAPPLY(label, ".empty"))) {
    # TURN LABELS ON
    if (!.all_na(c(label_text, label_stat))) {
      label <- TRUE
      # TURN LABELS OFF
    } else {
      label <- FALSE
    }
  }
  
  # LEGEND_TEXT
  if(.all_na(legend_text)){
    legend_text <- cyto_names(fr_list)
  }
  
  # TRANSFORM LABEL_TEXT_X
  if(!.all_na(label_text_x)){
    if(!.all_na(axes_trans)){
      if(channels[1] %in% names(axes_trans)){
        ind <- which(!is.na(label_text_x))
        lapply(ind, function(z){
          label_text_x[z] <<- axes_trans[[channels[1]]]$transform(label_text_x[z])
        })
      }
    }
  }
  
  # TRANSFORM LABEL_TEXT_Y
  if(!.all_na(label_text_y)){
    if(length(channels) == 2){
      if(!.all_na(axes_trans)){
        if(channels[2] %in% names(axes_trans)){
          ind <- which(!is.na(label_text_y))
          lapply(ind, function(z){
            label_text_y[z] <- axes_trans[[channels[2]]]$transform(label_text_y[z])
          })
        }
      }
    }
  }
  
  # POPULATIONS TO LABEL -------------------------------------------------------

  # LIST OF POPULATIONS - NEEDED FOR POSITION & STATISTICS
  if (label == TRUE) {
    pops <- LAPPLY(seq_len(SMP), function(z) {
      # INDICES
      ind <- TNP_split[[z]]
      # LAYER LABELLED - APPLY GATES
      if (!.all_na(label_text[ind]) | !.all_na(label_stat[ind])) {
        # GET ALL POPULATIONS (NOT JUST LABELLED ONES - CAN BE IMPROVED?)
        POPS <- .cyto_label_pops(fr_list[[z]],
          gate = gate,
          negate = negate
        )
        # LAYER NOT LABELLED
      } else {
        POPS <- rep(list(NA), NP)
      }
      return(POPS)
    })
  }

  # COMPUTE LABEL STATISTICS ---------------------------------------------------

  # STATISTICS
  label_stat <- .cyto_label_stat(fr_list,
    pops = pops,
    channels = channels,
    axes_trans = axes_trans,
    label_stat = label_stat,
    density_smooth = density_smooth
  )

  # COMBINE LABEL_TEXT & LABEL_STAT
  label_text <- .cyto_label_text(
    label_text,
    label_stat
  )
  
  # ARGUMENT SPLITTING ---------------------------------------------------------

  # ARGUMENTS
  args <- .args_list()

  # REPEAT ARGUMENTS
  args <- .cyto_plot_args_split(args)

  # UPDATE ARGUMENTS
  .args_update(args) 
  
  # CYTO_PLOT_CALL & LABEL RESET -----------------------------------------------

  # PREVIOUS CALL
  if(getOption("cyto_plot_method") == "flowFrame"){
    previous_call <- getOption("cyto_plot_call")
  }else{
    previous_call <- getOption("cyto_plot_match")
  }

  # SAVE CURRENT CALL - FLOWFRAME METHOD
  if(getOption("cyto_plot_method") == "flowFrame"){
    # 1D PLOT
    if(length(channels) == 1){
      current_call <- args[c("fr_list",
                             "fr_dens_list",
                             "channels",
                             "label",
                             "label_position",
                             "label_text_x",
                             "label_text_y",
                             "negate")]
    # 2D PLOT
    }else if(length(channels) == 2){
      current_call <- args[c("fr_list",
                             "channels",
                             "label",
                             "label_position",
                             "label_text_x",
                             "label_text_y",
                             "negate")]
    }
  # SAVE CURRENT CALL - FLOWSET METHOD  
  }else{
    # 1D PLOT
    if(length(channels) == 1){
      current_call <- args[c("channels",
                             "gate",
                             "label_position",
                             names(fr_dens_list),
                             "negate")]
    # 2D PLOT
    }else if(length(channels) == 2){
      current_call <- args[c("channels",
                             "gate",
                             "label_position",
                             "negate")]
    }
  }

  # RESET SAVED LABEL CO-ORDINATES - FLOWFRAME METHOD ONLY
  if (getOption("cyto_plot_method") == "flowFrame"){
    # RESET SAVED COORDS - NEW CALL OR NO CYTO_PLOT_SAVE
    if(isTRUE(all.equal(previous_call, current_call)) |
      getOption("cyto_plot_save") == FALSE) {
    options("cyto_plot_label_coords" = NULL)
    }
  }
  
  # GRAPHICS DEVICE ------------------------------------------------------------

  # PREPARE GRAPHICS DEVICE
  cyto_plot_new(popup)

  # PLOT CONSTRUCTION ----------------------------------------------------------

  # PULL DOWN ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_EMPTY
  .cyto_plot_empty(args)

  # DENSITY LAYERS
  if (length(channels) == 1) {
    .cyto_plot_density(args)
  }

  # POINT & CONTOUR LAYERS
  if (length(channels) == 2) {
    .cyto_plot_point(args)
  }
  
  # LABEL CO-ORDINATE INHERITANCE
  if(label == TRUE){
    # FLOWFRAME METHOD
    if(getOption("cyto_plot_method") == "flowFrame"){
      # CYTO_PLOT_SAVE & SAME CALL
      if (getOption("cyto_plot_save") &
        isTRUE(all.equal(previous_call, current_call))) {
        # SAVED LABEL CO-ORDINATES
        saved_label_coords <- getOption("cyto_plot_label_coords")
        if (!is.null(saved_label_coords)) {
          label_text_x <- saved_label_coords[, "x"]
          label_text_y <- saved_label_coords[, "y"]
        }
      # COMPUTE LABEL CO-ORDINATES
      } else {
        # COMPUTE OFFSET CO-ORDINATES
        if (label_position == "auto") {
          label_text_xy <- .cyto_label_coords(args)
          label_text_x <- label_text_xy[, "x"]
          label_text_y <- label_text_xy[, "y"]
        }
     }
    # FLOWSET METHOD
    }else{
      # SAME CALL - INHERIT CO-ORDINATES
      if(!is.null(getOption("cyto_plot_match")) &
         any(LAPPLY(getOption("cyto_plot_match"), function(z){
        identical(z, current_call)
      }))){
        # INHERIT SAVED CO-ORDINATES
        ind <- which(LAPPLY(getOption("cyto_plot_match"), function(z){
          identical(z, current_call)
        }))[1]
        # SAVED LABEL CO-ORDINATES 
        saved_label_coords <- getOption("cyto_plot_label_coords")[[ind]]
        if (!is.null(saved_label_coords)) {
          label_text_x <- saved_label_coords[, "x"]
          label_text_y <- saved_label_coords[, "y"]
        }
      # COMPUTE LABEL CO-ORDINATES
      }else{
        # COMPUTE OFFSET CO-ORDINATES
        if (label_position == "auto") {
          label_text_xy <- .cyto_label_coords(args)
          label_text_x <- label_text_xy[, "x"]
          label_text_y <- label_text_xy[, "y"]
        }
      }
    }
  }

  # PULL DOWN UPDATED ARGUMENTS
  args <- .args_list()
  
  # GATES & LABELS
  if (!.all_na(gate)) {
    # PLOT GATE & ASSOCIATED LABELS
    label_text_xy <- .cyto_plot_gate(args)
    # LABELS
  } else {
    if (label == TRUE) {
      label_text_xy <- .cyto_plot_label(args)
    }
  }

  # SAVE LABEL CO-ORDINATES ----------------------------------------------------

  # SAVE TO CYTO_PLOT_LABEL_COORDS - INACTIVE CYTO_PLOT_SAVE
  if(label == TRUE){
    # FLOWFRAME METHOD
    if(getOption("cyto_plot_method") == "flowFrame" &
       getOption("cyto_plot_save") == FALSE){
      options("cyto_plot_label_coords" = label_text_xy)
    # FLOWSET METHOD
    }else{
      # SAVE CO-ORDINATES - CYTO_PLOT_SAVE INACTIVE
      if(getOption("cyto_plot_save") == FALSE){
        if(is.null(getOption("cyto_plot_label_coords"))){
          options("cyto_plot_label_coords" = list(label_text_xy))
        }else{
          options("cyto_plot_label_coords" = c(getOption("cyto_plot_label_coords"),
                                             list(label_text_xy)))
        }
      }
    }
  }

  # SAVE PLOT CALL -------------------------------------------------------------
  
  # FLOWFRAME METHOD
  if(getOption("cyto_plot_method") == "flowFrame" &
     getOption("cyto_plot_save") == FALSE){
    # UPDATE CYTO_PLOT_CALL
    options("cyto_plot_call" = current_call)
  # FLOWSET METHOD
  }else{
    # UPDATE CYTO_PLOT_MATCH - CYTO_PLOT_SAVE INACTIVE
    if(getOption("cyto_plot_save") == FALSE){
      if(is.null(previous_call)){
        options("cyto_plot_match" = list(current_call))
      }else{
        options("cyto_plot_match" = c(previous_call, list(current_call)))
      }
    }
  }

  # RECORD/SAVE ----------------------------------------------------------------

  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save") == TRUE) {
    if (inherits(x, getOption("cyto_plot_method"))) {
      if (!getOption("cyto_plot_custom")) {
        # CLOSE GRAPHICS DEVICE
        dev.off()
      }
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RESET CYTO_PLOT_METHOD
  if(getOption("cyto_plot_method") == "flowFrame"){
    options("cyto_plot_method" = NULL)
  }
  
  # RETURN RECORDED PLOT
  invisible(recordPlot())
}
