## CYTO_PLOT_GATING_SCHEME -----------------------------------------------------

#' Plot Cytometry Gating Strategies
#'
#' \code{cyto_plot_gating_scheme} automatically plots the entire gating scheme
#' and has full support for gate tracking and back-gating through
#' \code{gate_track} and \code{back_gate}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param gatingTemplate name of the gatingTemplate csv file used to gate
#'   \code{x}. If not supplied the gating scheme will be obtained directly from
#'   the GatingSet.
#' @param group_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to "name" by default to plot each sample
#'   separately.
#' @param overlay names of the populations to overlay for back-gating, set to NA
#'   by default to turn off overlays. For simplicity, overlay can be set to
#'   \code{"descendants"} or \code{"children"} to overlay descendant or child
#'   populations within each plot.
#' @param gate_track logical indicating whether gate colour should be tracked
#'   throughout gating scheme, set to TRUE by default.
#' @param header character string to use as the header for the plot layout, set
#'   to "Gating Scheme" by default.
#' @param popup logical indicating whether the gating scheme should be plotted
#'   in a pop-up window, set to FALSE by default.
#' @param layout a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param point_col colour of points in 2D plots set to NA to use default
#'   red-blue colour scale. Control the colour of overlays by supplying multiple
#'   colours to this argument (e.g. c("blue","red")).
#' @param label_text_size numeric to control the size of text in the plot
#'   labels, set to 0.8 by default.
#' @param gate_line_col colour to use for gate borders, set to "red" by default.
#' @param legend logical indicating whether a legend should be included when an
#'   overlay is supplied.
#' @param legend_text vector of character strings to use for legend when an
#'   overlay is supplied.
#' @param legend_text_size character expansion for legend text, set to 1.2 by
#'   default.
#' @param header_text_font numeric to control the font of the header, set to 2
#'   by default for bold font.
#' @param header_text_size numeric to control the size of the header text, set
#'   to 1 by default.
#' @param header_text_col colour to use for header text, set to "black" by
#'   default.
#' @param border_line_width line width for plot border, set to 3 when gate_track
#'   is TRUE.
#' @param border_line_col line colour for plot border, set to "black" by
#'   default.
#' @param title_text_col colour to use for title text, set to "black" by
#'   default.
#' @param ... extra arguments passed to \code{\link{cyto_plot}}.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom graphics plot.new legend
#' @importFrom flowWorkspace gh_pop_get_descendants gh_pop_is_negated
#'   gh_pop_get_gate
#' @importFrom openCyto gt_gating
#' @importFrom methods is
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- cyto_compensate(gs)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs)
#'
#' # Gate using cyto_gate_draw
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # Gating scheme
#' cyto_plot_gating_scheme(gs[[4]])
#'
#' # Back-gating
#' cyto_plot_gating_scheme(gs[[32]],
#'   back_gate = TRUE
#' )
#'
#' # Gate-tracking
#' cyto_plot_gating_scheme(gs[1:2],
#'   gate_track = TRUE
#' )
#' @name cyto_plot_gating_scheme
NULL

#' @noRd
#' @export
cyto_plot_gating_scheme <- function(x, ...) {
  UseMethod("cyto_plot_gating_scheme")
}

#' @rdname cyto_plot_gating_scheme
#' @export
cyto_plot_gating_scheme.GatingSet <- function(x,
                                              gatingTemplate = NULL,
                                              group_by = "name",
                                              overlay = NA,
                                              gate_track = FALSE,
                                              header,
                                              popup = FALSE,
                                              layout = NULL,
                                              point_col = NA,
                                              label_text_size = 0.8,
                                              gate_line_col = "red",
                                              legend = TRUE,
                                              legend_text = NULL,
                                              legend_text_size = 1.2,
                                              header_text_font = 2,
                                              header_text_size = 1,
                                              header_text_col = "black",
                                              border_line_width = NULL,
                                              border_line_col = "black",
                                              title_text_col = "black", ...) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- .par(c("oma", "mfrow"))
  
  # CHECKS ---------------------------------------------------------------------
  
  # METHOD & RESET
  if (is.null(getOption("cyto_plot_method"))) {
    # SET PLOT METHOD
    options("cyto_plot_method" = "Gating/GatingSet")
    # RESET PLOT METHOD & GRAPHICAL PARAMETERS ON EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_method" = NULL)
    })
  } else {
    # RESET GRAPHICAL PARAMETERS ON EXIT
    on.exit(par(old_pars))
  }
  
  # THEME ----------------------------------------------------------------------
  
  # PULL DOWN CYTO_PLOT_THEME
  cyto_plot_theme <- getOption("cyto_plot_theme")
  
  # APPLY GATINGTEMPLATE -------------------------------------------------------
  
  # GATINGTEMPLATE SUPPLIED
  if (!is.null(gatingTemplate)) {
    cyto_gatingTemplate_apply(x, gatingTemplate)
  }
  
  # PREPARE NODES --------------------------------------------------------------
  
  # NODES PATHS (ALL NODES)
  nodes_full <- cyto_nodes(x, path = "full")
  nodes_auto <- cyto_nodes(x, path = "auto")
  
  # PREPARE OVERLAY ------------------------------------------------------------
  
  # OVERLAY - FULL PATH
  if(!.all_na(overlay)){
    # DESCENDANTS
    if(all(unique(overlay) %in% c("descendants", "children"))){
      overlay_nodes <- nodes_full
      # SPECIFIC
    }else{
      overlay_nodes <- cyto_nodes_convert(x,
                                          nodes = overlay,
                                          path = "full")
    }
    # NO OVERLAY
  }else{
    overlay_nodes <- NULL
  }
  
  # PREPARE GATINGTEMPLATE -----------------------------------------------------
  
  # GATING SCHEME FROM GATINGHIERARCHY
  gatingTemplate <- gh_generate_template(x[[1]])
  # EXTRACT CHANNELS FROM GATINGTEMPLATE
  chans <- strsplit(gatingTemplate$dims, ",", fixed = TRUE)
  # WATCH FOR ENTRIES WITHOUT CHANNELS
  chans <- lapply(chans, function(z) {
    if (length(z) == 0) {
      return(c("", ""))
    } else if (length(z) == 1) {
      return(c(z, ""))
    } else {
      return(z)
    }
  })
  chans <- do.call("rbind", chans)
  
  # EXTRACT GATING SCHEME INFORMATION
  gatingTemplate <- data.frame(
    parent = gatingTemplate$parent,
    parent_auto = cyto_nodes_convert(x,
                                     nodes = gatingTemplate$parent,
                                     path = "auto"
    ),
    alias = gatingTemplate$alias,
    alias_auto = LAPPLY(seq_len(nrow(gatingTemplate)), function(z) {
      cyto_nodes_convert(x,
                         nodes = gatingTemplate$alias[z],
                         anchor = gatingTemplate$parent[z],
                         path = "auto"
      )
    }),
    alias_full = LAPPLY(seq_len(nrow(gatingTemplate)), function(z) {
      cyto_nodes_convert(x,
                         nodes = gatingTemplate$alias[z],
                         anchor = gatingTemplate$parent[z],
                         path = "full"
      )
    }),
    xchannel = chans[, 1],
    ychannel = chans[, 2],
    stringsAsFactors = FALSE
  )
  
  # PARENTS (PLOT ARGUMENTS - BORDERS & TITLES) - FULL PATH
  parent_nodes <- unique(gatingTemplate$parent)
  
  # ALIASES (GATE ARGUMENTS - GATE_LINE_COL & GATE_LINE_WIDTH) - FULL PATH
  alias_nodes <- unique(paste0(
    gatingTemplate[gatingTemplate$parent %in% parent_nodes & 
                     !.empty(xchannel) & 
                     !.empty(ychannel), c("parent", "/", "alias")]))
  
  # NODE CHARACTERISTICS -------------------------------------------------------
  
  # NODE TYPES
  node_types <- LAPPLY(nodes_full, function(node_full){
    node_type <- c()
    if(node_full %in% parent_nodes){
      node_type <- c(node_type, "parent")
    }
    if(node_full %in% alias_nodes){
      node_type <- c(node_type, "alias")
    }
    if(node_full %in% overlay_nodes){
      node_type <- c(node_type, "overlay")
    }
    return(paste0(node_type, sep = "/"))
  })
  
  # NODE ARGUMENTS - SET NA AS DEFAULT
  node_arguments <- data.frame(
    "node_full" = nodes_full,
    "node_auto" = nodes_auto,
    "node_type" = node_types, 
    "point_col" = rep(NA, length(nodes_full)),
    "gate_line_col" = rep(NA, length(nodes_full)),
    "density_fill" = rep(NA, length(nodes_full)),
    "border_line_col" = rep(NA, length(nodes_full)),
    "title_text_col" = rep(NA, length(nodes_full)),
    stringsAsFactors = FALSE
  )
  
  # PREPARE NODE_ARGUMENTS -----------------------------------------------------
  
  # BORDER_LINE_COL/TITLE_TEXT_COL (PARENT ARGUMENTS)
  
  # GATE_LINE_COL/GATE_LINE_WIDTH (ALIAS ARGUMENTS)
  
  # POINT_COL/DENSITY_FILL (OVERLAY ARGUMENTS)
  
  # NODE COLOURS
  if (.all_na(point_col)) {
    # OVERLAYS
    if (!.all_na(overlay)) {
      if (!is.null(cyto_plot_theme)) {
        if ("point_cols" %in% names(cyto_plot_theme)) {
          point_cols <- cyto_plot_theme[["point_cols"]]
        } else {
          point_cols <- .cyto_plot_colour_palette(type = "point_cols")
        }
      } else {
        point_cols <- .cyto_plot_colour_palette(type = "point_cols")
      }
      # COLOUR PALETTE
      point_col_palette <- colorRampPalette(point_cols)
      # NODE COLOURS
      node_cols <- point_col_palette(nodes_count)
      # NO OVERLAY
    } else {
      node_cols <- rep(NA, nodes_count)
    }
  } else {
    node_cols <- rep(point_col, length.out = nodes_count)
  }
  
  # GATE TRACKING
  if (gate_track == TRUE) {
    node_details[, "gate_line_col"] <- node_cols
    node_details[, "border_line_col"] <- node_cols
  }
  
  # GATE TRACKING AND/OR BACK GATING
  if (gate_track == TRUE | .all_na(overlay)) {
    node_details[, "title_text_col"] <- node_cols
  }
  
  # BORDER_LINE_WIDTH
  if (missing(border_line_width)) {
    if (gate_track == TRUE) {
      border_line_width <- 3
    } else {
      border_line_width <- 1
    }
  }
  
  # RESET SAVED LABEL CO-ORDINATES ---------------------------------------------
  
  # PREVIOUS CALL
  previous_call <- getOption("cyto_plot_call")
  
  # CURRENT CALL (gh does not match)
  current_call <- list(
    "sample" = cyto_names(x),
    "gts" = gts,
    "overlay" = overlay,
    "gate_track" = gate_track
  )
  
  # PREVIOUS CALL BELONGS TO DIFFERENT METHOD
  if (!all(names(previous_call) %in% names(current_call))) {
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
  
  # PREPARE SAMPLES ------------------------------------------------------------
  
  # GROUP_BY
  gs_list <- cyto_group_by(x, group_by)
  
  # HEADER
  if (missing(header)) {
    header <- names(gs_list)
  }
  
  # PLOT CONSTRUCTION ----------------------------------------------------------
  
  # GATING SCHEME FOR EACH GATINGSET
  plots <- lapply(seq_along(gs_list), function(z) {
    
    # GatingSet
    gs <- gs_list[[z]]
    
    # HEADER
    if (!is.null(header)) {
      header <- header[[z]]
    }
    
    # POPUP
    if (popup == TRUE) {
      cyto_plot_new(popup)
    }
    
    # PLOT NUMBER (WATCH OUT BOOLEAN POPULATIONS)
    np <- nrow(unique(gts[
      !LAPPLY(paste(gts$xchannel, gts$ychannel), ".empty"),
      c("parent", "xchannel", "ychannel")
    ]))
    
    # LAYOUT (PARENT NUMBER)
    if (is.null(layout)) {
      if (legend == FALSE) {
        layout <- c(
          n2mfrow(np)[2],
          n2mfrow(np)[1]
        )
      } else {
        if (!.all_na(overlay)) {
          layout <- c(
            n2mfrow(np)[2],
            n2mfrow(np)[1]
          )
        } else {
          layout <- c(
            n2mfrow(np + 1)[2],
            n2mfrow(np + 1)[1]
          )
        }
      }
      par(mfrow = layout)
    } else {
      if (layout[1] != FALSE) {
        par(mfrow = layout)
      }
    }
    
    # HEADER SPACE
    if (!is.null(header)) {
      par(oma = c(0, 0, 3, 0))
    }
    
    # PLOT GATING SCHEME
    mapply(function(parents) {
      
      # PARENT ENTRIES
      gt <- gts[gts$parent == parents, ]
      
      # REMOVE BOOLEAN ENTRIES
      gt_gates <- gt[!LAPPLY(paste(gt$xchannel, gt$ychannel), ".empty"), ]
      
      # Parent may have gates in different channels
      # construct a plot for each channel set
      lapply(
        seq(1, nrow(unique(gt_gates[, c("xchannel", "ychannel")]))),
        function(x) {
          parent <- as.character(parents)
          xchannel <- as.character(unique(gt_gates[, "xchannel"])[x])
          ychannel <- as.character(unique(gt_gates[, "ychannel"])[x])
          channels <- c(xchannel, ychannel)
          # RETAIN BOOLEAN POPULATIONS
          gt_split <- unique(as.vector(gt[gt$parent == parents &
                                            gt$xchannel %in% c(xchannel, "") &
                                            gt$ychannel %in% c(ychannel, ""), ]))
          alias <- gt_split[, "alias"]
          alias_auto <- gt_split[, "alias_auto"]
          alias_full <- gt_split[, "alias_full"]
          
          print(alias_full)
          
          # HISTOGRAMS
          if (channels[1] == channels[2] |
              (!.empty(channels[1]) & .empty(channels[2]))) {
            channels <- channels[1]
          }
          
          # OVERLAYS
          if(!.all_na(overlay)){
            # DESCENDANTS
            if(all(unique(overlay) == "descendants")){
              overlay <- unique(c(alias_full, LAPPLY(alias_full, function(z){
                gh_pop_get_descendants(gs[[1]], z, path = "full")
              })))
              # CHILDREN
            }else if(all(unique(overlay) == "children")){
              overlay <- unique(c(alias_full, LAPPLY(alias_full, function(z){
                gh_pop_get_children(gs[[1]], z, path = "full")
              })))
              # SPECIFIC POPULATIONS 
            }else{
              overlay <- unique(LAPPLY(alias_full, function(z){
                LAPPLY(overlay, function(y){
                  if(y %in% c(alias_full, 
                              gh_pop_get_descendants(gs[[1]], 
                                                     z,
                                                     path = "full"))){
                    return(y)
                  }else{
                    return(NULL)
                  }
                })
              }))
            }
          }
          
          # POINT_COL
          point_col <- node_details[, "point_col"][match(parent, nodes_full)]
          if (!.all_na(overlay)) {
            point_col <- c(
              point_col,
              node_details[, "point_col"][match(
                overlay[!is.na(overlay)], nodes_full
              )]
            )
          }
          
          # DENSITY_FILL
          density_fill <- node_details[, "density_fill"][match(
            c(parent, overlay[!is.na(overlay)]),
            nodes_full
          )]
          if (!.all_na(overlay)) {
            density_fill <- c(
              density_fill,
              node_details[, "density_fill"][match(
                overlay[!is.na(overlay)], nodes_full
              )]
            )
          }
          
          # GATE_LINE_COL
          gate_line_col <- node_details[, "gate_line_col"][match(
            alias_full, nodes_full
          )]
          
          # BORDER_LINE_COL
          border_line_col <- node_details[, "border_line_col"]
          border_line_col <- border_line_col[match(parent, nodes_full)]
          
          # TITLE_TEXT_COL
          title_text_col <- node_details[, "title_text_col"]
          title_text_col <- title_text_col[match(parent, nodes_full)]
          
          # SHORT PATHS
          parent <- cyto_nodes_convert(gs,
                                       nodes = parent,
                                       path = "auto")
          if(!.all_na(overlay)){
            overlay <- cyto_nodes_convert(gs,
                                          nodes = overlay,
                                          path = "auto")
          }
          
          # TITLE
          if (parent == "root") {
            title <- "All Events"
          } else {
            title <- parent
          }
          
          # CYTO_PLOT (TURN OFF LAYOUT)
          cyto_plot(gs,
                    parent = parent,
                    group_by = "all",
                    alias = alias,
                    overlay = overlay,
                    channels = channels,
                    title = title,
                    legend = FALSE,
                    layout = FALSE,
                    point_col = point_col,
                    density_fill = density_fill,
                    gate_line_col = gate_line_col,
                    border_line_col = border_line_col,
                    border_line_width = border_line_width,
                    label_text_size = label_text_size,
                    ...
          )
        }
      )
    }, parents)
    
    # HEADER
    if (!is.null(header)) {
      mtext(header,
            outer = TRUE,
            font = header_text_font,
            cex = header_text_size,
            col = header_text_col
      )
    }
    
    # LEGEND
    if (!.all_na(overlay)) {
      
      # DESCENDANTS OR CHILDREN
      if(all(unique(overlay) %in% c("descendants", "children"))){
        overlay <- nodes_auto
      }else{
        overlay <- cyto_nodes_convert(x,
                                      nodes = overlay,
                                      path = "auto"
        )
      }
      
      # LEGEND
      if (legend == TRUE) {
        
        # LEGEND_TEXT
        if (is.null(legend_text)) {
          if (class(overlay) == "character") {
            legend_text <- overlay
          } else {
            stop("Please supply vector of names to use in the legend.")
          }
        }
        
        # DUMMY PLOT - LEGEND
        plot.new()
        
        # LEGEND
        legend("center",
               legend = legend_text,
               fill = node_details[, "point_col"][match(overlay, nodes_auto)],
               bty = "n",
               cex = legend_text_size,
               x.intersp = 0.5
        )
      }
    }
    # RECORD PLOT
    cyto_plot_record()
  })
  names(plots) <- names(gs_list)
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE
  if (getOption("cyto_plot_save")) {
    if (is(x, basename(getOption("cyto_plot_method")))) {
      
      # CLOSE GRAPHICS DEVICE
      dev.off()
      
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN LIST OF RECORDED PLOTS
  invisible(plots)
}

#' @rdname cyto_plot_gating_scheme
#' @export
cyto_plot_gating_scheme.GatingHierarchy <- function(x,
                                                    gatingTemplate = NULL,
                                                    back_gate = FALSE,
                                                    gate_track = FALSE,
                                                    show_all = FALSE,
                                                    header,
                                                    popup = FALSE,
                                                    layout,
                                                    point_col = NA,
                                                    label_text_size = 0.8,
                                                    gate_line_col = "red",
                                                    legend = TRUE,
                                                    legend_text,
                                                    legend_text_size = 1.2,
                                                    header_text_font = 2,
                                                    header_text_size = 1,
                                                    header_text_col = "black",
                                                    border_line_width,
                                                    border_line_col = "black",
                                                    title_text_col = "black",
                                                    ...) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # CURRENT PARAMETERS
  old_pars <- .par(c("oma", "mfrow"))
  
  # CHECKS ---------------------------------------------------------------------
  
  # METHOD & RESET
  if (is.null(getOption("cyto_plot_method"))) {
    # SET PLOT METHOD
    options("cyto_plot_method" = "Gating/GatingHierarchy")
    # RESET PLOT METHOD & GRAPHICAL PARAMETERS ON EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_method" = NULL)
    })
  } else {
    # RESET GRAPHICAL PARAMETERS ON EXIT
    on.exit(par(old_pars))
  }
  
  # GATINHIERARCHY
  gh <- x
  
  # GATINGTEMPLATE SUPPLIED
  if (!is.null(gatingTemplate)) {
    cyto_gatingTemplate_apply(gh, gatingTemplate)
  }
  
  # POPULATIONS
  pops <- cyto_nodes(gh, path = "auto")
  
  # NUMBER POPULATIONS
  nodes_count <- length(pops)
  
  # THEME ----------------------------------------------------------------------
  
  # PULL DOWN CYTO_PLOT_THEME
  cyto_plot_theme <- getOption("cyto_plot_theme")
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # BACK_GATE
  if (back_gate[1] == TRUE) {
    back_gate <- "all"
  }
  
  # BACK GATING
  if (back_gate[1] != FALSE) {
    if (back_gate[1] == "all") {
      overlay <- pops[-1]
    } else {
      overlay <- back_gate
    }
  } else {
    overlay <- NA
  }
  
  # POPULATION COLOURS - NOT SPECIFIED
  if (.all_na(point_col)) {
    # BACK GATE OR GATE TRACK (SELECT FROM POINT_COLS)
    if (back_gate[1] != FALSE | gate_track == TRUE) {
      if (!is.null(cyto_plot_theme)) {
        if ("point_cols" %in% names(cyto_plot_theme)) {
          point_cols <- cyto_plot_theme[["point_cols"]]
        }
      } else {
        point_cols <- .cyto_plot_colour_palette(type = "point_cols")
      }
      # COLOUR PALETTE
      cols <- colorRampPalette(point_cols)
      # POPULATION COLOURS
      node_cols <- cols(nodes_count)
    } else {
      node_cols <- rep(NA, nodes_count)
    }
    # POPULATION COLOURS - MANUALLY SUPPLIED
  } else {
    node_cols <- rep(point_col, length.out = nodes_count)
  }
  
  # POPULATION COLOURS (DATA FRAME)
  node_details <- data.frame(
    "node" = pops,
    "point_col" = node_cols,
    "gate_line_col" = rep(gate_line_col, nodes_count),
    "density_fill" = node_cols,
    "border_line_col" = rep(border_line_col, nodes_count),
    "title_text_col" = rep(title_text_col, nodes_count),
    stringsAsFactors = FALSE
  )
  
  # GATE TRACKING
  if (gate_track == TRUE) {
    node_details[, "gate_line_col"] <- node_cols
    node_details[, "border_line_col"] <- node_cols
  }
  
  # GATE TRACKING AND/OR BACK GATING
  if (gate_track == TRUE | back_gate[1] != FALSE) {
    node_details[, "title_text_col"] <- node_cols
  }
  
  # BORDER_LINE_WIDTH
  if (missing(border_line_width)) {
    if (gate_track == TRUE) {
      border_line_width <- 3
    } else {
      border_line_width <- 1
    }
  }
  
  # HEADER
  if (missing(header)) {
    header <- cyto_names(gh)
  }
  
  # GATING SCHEME --------------------------------------------------------------
  
  # GATING SCHEME FROM GATINGHIERARCHY
  gt <- gh_generate_template(gh)
  # EXTRACT CHANNELS FROM GATINGTEMPLATE
  chans <- strsplit(gt$dims, ",", fixed = TRUE)
  # WATCH FOR ENTRIES WITHOUT CHANNELS
  chans <- lapply(chans, function(z) {
    if (length(z) == 0) {
      return(c("", ""))
    } else if (length(z) == 1) {
      return(c(z, ""))
    } else {
      return(z)
    }
  })
  chans <- do.call("rbind", chans)
  # EXTRACT GATING SCHEME INFORMATION
  gts <- data.frame(
    parent = basename(gt$parent),
    alias = gt$alias,
    xchannel = chans[, 1],
    ychannel = chans[, 2],
    stringsAsFactors = FALSE
  )
  
  # UNIQUE PARENTS
  parents <- unique(gts$parent)
  
  # RESET SAVED LABEL CO-ORDINATES ---------------------------------------------
  
  # PREVIOUS CALL
  previous_call <- getOption("cyto_plot_call")
  
  # CURRENT CALL (gh does not match)
  current_call <- list(
    "sample" = cyto_names(gh),
    "gts" = gts,
    "back_gate" = back_gate,
    "gate_track" = gate_track
  )
  
  # PREVIOUS CALL BELONGS TO DIFFERENT METHOD
  if (!all(names(previous_call) %in% names(current_call))) {
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
  
  # PLOT LAYOUT ----------------------------------------------------------------
  
  # nodeUP
  if (popup == TRUE) {
    cyto_plot_new(popup)
  }
  
  # PLOT NUMBER (WATCH OUT BOOLEAN POPULATIONS)
  np <- nrow(unique(gts[
    !LAPPLY(paste(gts$xchannel, gts$ychannel), ".empty"),
    c("parent", "xchannel", "ychannel")
  ]))
  
  # LAYOUT (PARENT NUMBER)
  if (missing(layout)) {
    if (legend == FALSE) {
      layout <- c(
        n2mfrow(np)[2],
        n2mfrow(np)[1]
      )
    } else {
      if (back_gate == FALSE) {
        layout <- c(
          n2mfrow(np)[2],
          n2mfrow(np)[1]
        )
      } else {
        layout <- c(
          n2mfrow(np + 1)[2],
          n2mfrow(np + 1)[1]
        )
      }
    }
    par(mfrow = layout)
  } else {
    if (layout[1] != FALSE) {
      par(mfrow = layout)
    }
  }
  
  if (!is.null(header)) {
    par(oma = c(0, 0, 3, 0))
  }
  
  # PLOT CONSTRUCTION ----------------------------------------------------------
  
  # NEW PLOT PER PARENT
  mapply(function(parents) {
    
    # PARENT ENTRIES
    gt <- gts[gts$parent == parents, ]
    
    # REMOVE BOOLEAN ENTRIES
    gt_gates <- gt[!LAPPLY(paste(gt$xchannel, gt$ychannel), ".empty"), ]
    
    # Parent may have gates in different channels
    # construct a plot for each channel set (watch for boolean populations)
    lapply(
      seq(1, nrow(unique(gt_gates[, c("xchannel", "ychannel")]))),
      function(x) {
        parent <- as.character(parents)
        xchannel <- as.character(unique(gt_gates[, "xchannel"])[x])
        ychannel <- as.character(unique(gt_gates[, "ychannel"])[x])
        channels <- c(xchannel, ychannel)
        
        # RETAIN BOOLEAN POPULATIONS
        alias <- as.vector(gt[gt$parent == parents &
                                gt$xchannel %in% c(xchannel, "") &
                                gt$ychannel %in% c(ychannel, ""), "alias"])
        
        # BOOLEAN ORDER
        bool_gate <- LAPPLY(alias, function(z) {
          is(gh_pop_get_gate(gh, paste0(parent, "/", z)), "booleanFilter")
        })
        if (any(bool_gate)) {
          alias <- c(alias[-which(bool_gate)], alias[which(bool_gate)])
        }
        
        # Histograms - stacked & gated
        if (channels[1] == channels[2] |
            (!.empty(channels[1]) & .empty(channels[2]))) {
          channels <- channels[1]
        }
        
        # BACK_GATE
        if (back_gate[1] != FALSE) {
          if (back_gate[1] == "all") {
            if (!show_all) {
              overlay <- c(alias, LAPPLY(
                seq_along(alias),
                function(x) {
                  gh_pop_get_descendants(gh,
                                         alias[x],
                                         path = "auto"
                  )
                }
              ))
            }
          } else {
            if (!show_all) {
              if (any(back_gate %in%
                      c(alias, LAPPLY(seq_along(alias), function(x) {
                        gh_pop_get_descendants(gh,
                                               alias[x],
                                               path = "auto"
                        )
                      })))) {
                ind <- c(
                  alias,
                  LAPPLY(
                    seq_along(alias),
                    function(x) {
                      gh_pop_get_descendants(gh,
                                             alias[x],
                                             path = "auto"
                      )
                    }
                  )
                )
                overlay <- back_gate[back_gate %in% ind]
              } else {
                overlay <- NA
              }
            }
          }
        }
        
        # POINT_COL
        point_col <- node_details[, "point_col"]
        point_col <- point_col[match(c(parent, overlay), nodes)]
        
        # GATE_LINE_COL
        gate_line_col <- node_details[, "gate_line_col"]
        gate_line_col <- gate_line_col[match(alias, node_details[, "node"])]
        
        # DENSITY_FILL
        density_fill <- node_details[, "density_fill"]
        density_fill <- density_fill[match(
          c(parent, overlay),
          pops
        )]
        
        # BORDER_LINE_COL
        border_line_col <- node_details[, "border_line_col"]
        border_line_col <- border_line_col[match(c(parent, overlay), pops)]
        
        # TITLE_TEXT_COL
        title_text_col <- node_details[, "title_text_col"]
        title_text_col <- title_text_col[match(c(parent, overlay), pops)]
        
        # TITLE
        if (parent == "root") {
          title <- "All Events"
        } else {
          title <- parent
        }
        
        # CYTO_PLOT
        cyto_plot(gh,
                  parent = parent,
                  alias = alias,
                  overlay = overlay,
                  channels = channels,
                  title = title,
                  legend = FALSE,
                  point_col = point_col,
                  density_fill = density_fill,
                  gate_line_col = gate_line_col,
                  border_line_col = border_line_col,
                  border_line_width = border_line_width,
                  title_text_col = title_text_col,
                  label_text_size = label_text_size, ...
        )
      }
    )
  }, parents)
  
  # HEADER
  if (!is.null(header)) {
    mtext(header,
          outer = TRUE,
          font = header_text_font,
          cex = header_text_size,
          col = header_text_col
    )
  }
  
  # LEGEND
  if (!.all_na(overlay)) {
    if (legend == TRUE) {
      
      # LEGEND_TEXT
      if (missing(legend_text)) {
        if (class(overlay) == "character") {
          legend_text <- overlay
        } else {
          stop("Please supply vector of names to use in the legend.")
        }
      }
      
      # DUMMY PLOT - LEGEND
      plot.new()
      
      # LEGEND
      legend("center",
             legend = legend_text,
             fill = node_details[, "point_col"][match(overlay, pops)],
             bty = "n",
             cex = legend_text_size,
             x.intersp = 0.5
      )
    }
  }
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE
  if (getOption("cyto_plot_save")) {
    if (is(x, basename(getOption("cyto_plot_method")))) {
      
      # CLOSE GRAPHICS DEVICE
      dev.off()
      
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN RECORDED PLOT
  invisible(cyto_plot_record())
}
