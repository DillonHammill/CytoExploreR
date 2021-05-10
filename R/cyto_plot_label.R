## CYTO_PLOT_LABEL -------------------------------------------------------------

#' Add labels to cyto_plot
#'
#' Interactively label existing plots with boxed labels containing text and
#' statistics. \code{cyto_plot_label} can only label the last plot and can only
#' add labels on a per sample basis (i.e. single layer only). Locations for each
#' label must be manually selected. More complex labelling can be achieved
#' directly through \code{cyto_plot}.
#'
#' @param x object of class \code{flowFrame}.
#' @param parent name of the parent population to extract when a
#'   \code{GatingHierarchy} object is supplied.
#' @param channels a vector indicating the fluorescent channel(s) to be used for
#'   gating.
#' @param alias name of teh gated populations to label when a GatingHierarchy
#'   object is supplied. This is equivalent to the \code{gate} argument for
#'   \code{flowFrame} objects.
#' @param gate object of class
#'   \code{\link[flowCore:rectangleGate-class]{rectangleGate}},
#'   \code{\link[flowCore:polygonGate-class]{polygonGate}},
#'   \code{\link[flowCore:ellipsoidGate-class]{ellipsoidGate}},
#'   \code{\link[flowCore:quadGate-class]{quadGate}}, \code{"list"} or
#'   \code{\link[flowCore:filters-class]{filters}}.
#' @param label_text character string to include in the label above the
#'   statistic (e.g. population name(s)).
#' @param label_stat indicates the type of statistic to include in the label,
#'   can be either code{"count"}, \code{"median"}, \code{"mean"}, \code{"mode"},
#'   \code{"geo mean"} or \code{"CV"}. Only count and percent statistics are
#'   supported for 2D plots.
#' @param label_text_x vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param label_text_y vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param trans object of class \code{\link[flowWorkspace]{transformerList}}
#'   which was used to transform the fluorescent channels of the supplied
#'   flowFrame.
#' @param negate logical indicating whether a label should be included for the
#'   neagted population (i.e. events outside the gate). Set to FALSE by default
#'   to only calculate statistics for events within the gate.
#' @param label_text_font integer [1,2,3,4] passed to \code{text} to alter the
#'   font, set to \code{2} by default for a bold font.
#' @param label_text_size numeric character expansion used to control the size
#'   of the text in the labels, set to \code{0.8} by default. See \code{?text}
#'   for details.
#' @param label_text_col specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param label_fill fill colour to use for labels, set to "white" by default.
#' @param label_fill_alpha numeric [0,1] controls the transparency of the fill
#'   colour, set to \code{0.6} by default.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param hist_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust kernel density for mode
#'   calculation.
#' @param ... not in use.
#'
#' @return add a boxed text label to an existing cyto_plot.
#'
#' @importFrom flowCore Subset
#' @importFrom methods is
#' @importFrom flowWorkspace gh_pop_get_gate gh_pop_is_negated
#' @importFrom openCyto gh_generate_template
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(fs))
#' gs <- transform(gs, trans)
#'
#' # Gate using gate_draw
#' gt_gating(Activation_gatingTemplate, gs)
#'
#' # Plot
#' cyto_plot(gs[[32]],
#'   parent = "CD4 T Cells",
#'   channels = c("CD69")
#' )
#'
#' # Label - median fluorescent intensity
#' cyto_plot_label(cyto_data_extract(gs, "CD4 T Cells")[[1]][[32]],
#'   channels = "CD69",
#'   label_stat = "median",
#'   label_text = "MedFI",
#'   trans = trans,
#'   label_text_x = 3,
#'   label_text_y = 50
#' )
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @name cyto_plot_label
NULL

#' @noRd
#' @export
cyto_plot_label <- function(x, ...){
  UseMethod("cyto_plot_label")
}

#' @rdname cyto_plot_label
#' @export
cyto_plot_label.GatingHierarchy <- function(x,
                                            parent,
                                            alias = NA,
                                            channels = NULL,
                                            display = 1,
                                            gate = NA,
                                            negate = NA,
                                            label_text,
                                            label_stat,
                                            label_text_x = NA,
                                            label_text_y = NA,
                                            label_text_font = 2,
                                            label_text_size = 0.8,
                                            label_text_col = "black",
                                            label_fill = "white",
                                            label_fill_alpha = 0.6,
                                            hist_smooth = 1, ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # PARENT
  if(missing(parent)){
    stop("Supply the name of the parental population to label.")
  }
  
  # CHANNELS
  if(is.null(channels)){
    stop("Supply the names of the channel(s) used to construct the plot.")
  }
  
  # EXTRACT DATA FROM GATINGHIERARCHY ------------------------------------------
  
  # PARENTAL FLOWFRAME
  fr <- cyto_data_extract(x, parent)[[1]]
  
  # TRANSFORMATIONS
  trans <- cyto_transformers_extract(x)
  
  # EXTRACT GATE(S) ------------------------------------------------------------
  
  # EMPTY ALIAS - EXTRACT ALIAS IN SUPPLIED CHANNEL(S)
  if (.empty(alias)) {
    # Plot all appropriate gates if alias is an empty character string
    if (all(alias == "")) {
      gt <- gh_generate_template(x)
      gt <- gt[basename(gt$parent) == parent, ]
      
      # Match both channels for 2D plots
      if (length(channels) == 2) {
        alias <- gt$alias[gt$dims == paste(channels, collapse = ",")]
        
        # At least 1 channel match
      } else if (length(channels) == 1) {
        ind <- lapply(gt$dims, function(z) {
          grep(channels, z)
        })
        ind <- LAPPLY(ind, "length") != 0
        alias <- gt$alias[ind]
      }
      
      # No gates constructed in the supplied channels
      if (length(alias) == 0) {
        alias <- NA
      }
    }
  }
  
  # ALIAS MULTIPLE POPULATIONS - QUADGATE
  if(!.all_na(alias)){
    alias <- LAPPLY(alias, function(z){
      if(grepl(",", z)){
        strsplit(z, ",")[[1]]
      }
    })
  }
  
  # EXTRACT GATE(S)
  if (!.all_na(alias)) {
    # POPULATIONS
    PNS <- length(alias)
    # EXTRACT GATE
    gate <- LAPPLY(seq_len(PNS), function(z) {
      # GATE
      gt <- gh_pop_get_gate(x, paste(parent, alias[z], sep = "/"))
      # BOOLEANFILTER - SUPPORT NEGATED & FILTERS
      if(is(gt, "booleanFilter")){
        # BOOLEAN LOGIC -  !ALIAS1&!ALIAS2&!ALIAS3
        bool <- gt@deparse
        # POPULATIONS SEPARATED BY & (!ALIAS 1 !ALIAS2 ...)
        bool_split <- strsplit(bool, "&")[[1]]
        # SUPPORTED BOOLEANFILTER
        if(all(grepl("!", bool_split)) &
           all(substring(bool_split, 2) %in% cyto_nodes(x, path  = "auto"))){
          # UPDATE ALIAS - MUST CONTAIN REFERENCE POPULATIONS
          alias <<- c(substring(bool_split, 2),
                      alias[z])
          # GATE
          gt <- LAPPLY(alias[-length(alias)], function(z){
            gh_pop_get_gate(x, paste(parent, z, spe = "/"))
          })
          # SET NEGATE TO TRUE
          negate <<- TRUE
        }else{
          message("Skipping unsupported booolean gates...")
          gt <- NA
        }
      }
      return(gt)
    })
    names(gate) <- alias
  }
  
  # LIST OF UNIQUE GATE OBJECTS
  if(!.all_na(gate)){
    cyto_gate_prepare(gate, channels)
  }
  
  # NEGATE
  if(!.all_na(gate)){
    # ALL NEGATED GATES - SET NEGATE TO TRUE
    if(all(lapply(gate, function(z){gh_pop_is_negated(z)})))
      # NEGATE NOT SPECIFIED
      if(missing(negate)){
        negate <- TRUE
      }
  }
  
  # PREPARE LABEL ARGUMENTS ----------------------------------------------------
  
  # LABEL_TEXT
  if(missing(label_text)){
    # GATE FILTERID - ALIAS & FILTERID SHOULD MATCH
    if(!.all_na(gate)){
      label_text <- LAPPLY(gate, function(z){z@filterId})
      # NO LABEL_TEXT
    }else{
      label_text <- NA
    }
  }
  
  # LABEL_STAT
  if(missing(label_stat)){
    # GATE
    if(!.all_na(gate)){
      label_stat <- rep("freq", length(alias))
      # NO GATE
    }else{
      label_stat <- NA
    }
  }
  
  # CALL FLOWFRAME METHOD ------------------------------------------------------
  
  # RETAIN LABEL CO-ORDINATES
  label_text_xy <- cyto_plot_label(x = fr,
                                   channels = channels,
                                   trans = trans,
                                   display = display,
                                   gate = gate,
                                   negate = negate,
                                   label_text = label_text,
                                   label_stat = label_stat,
                                   label_text_x = label_text_x,
                                   label_text_y = label_text_y,
                                   label_text_font = label_text_font,
                                   label_text_size = label_text_size,
                                   label_text_col = label_text_col,
                                   label_fill = label_fill,
                                   label_fill_alpha = label_fill_alpha,
                                   hist_smooth = hist_smooth)
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)
  
}

#' @rdname cyto_plot_label
#' @export
cyto_plot_label.flowFrame <- function(x,
                                      channels = NULL,
                                      trans = NA,
                                      display = 1,
                                      gate = NA,
                                      negate = FALSE,
                                      label_text = NA,
                                      label_stat = NA,
                                      label_text_x = NA,
                                      label_text_y = NA,
                                      label_text_font = 2,
                                      label_text_size = 0.8,
                                      label_text_col = "black",
                                      label_fill = "white",
                                      label_fill_alpha = 0.6,
                                      hist_smooth = 1, ...){
  
  # CHECKS ---------------------------------------------------------------------
  
  # FLOWFRAME
  if(!is(x, "flowFrame")){
    stop("'x' must be a flowFrame object.")
  }
  
  # SAMPLING
  if(display != 1){
    x <- cyto_sample(x, display = display, seed = 56)
  }
  
  # STATISTIC
  if(!.all_na(label_stat)){
    # VALID STATISTICS
    label_stat <- LAPPLY(label_stat, ".cyto_stat_check")
    # SUPPORTED STATISTICS
    if(length(channels) == 2 & 
       !all(label_stat %in% c("count","freq"))){
      message("Only frequency and count statistics are supported in 2D plots")
      # REPLACE WITH NA
      label_stat[!label_stat %in% c("count", "freq")] <- NA
    }
  }
  
  # CYTO_PLOT_THEME ------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # PREPARE GATES --------------------------------------------------------------
  
  # LIST OF GATE OBJECTS
  if(!.all_na(gate)){
    gate <- cyto_gate_prepare(gate, channels)
  }
  
  # POPULATIONS TO LABEL -------------------------------------------------------
  
  # POPULATIONS
  PNS <- .cyto_gate_count(gate)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # REPEAT ARGUMENTS - GATE
  if(!.all_na(gate)){
    
    # REPEAT LABEL ARGS PER POPULATION
    args[["label_text"]] <- rep(c(label_text, rep(NA, PNS)), 
                                length.out = PNS)
    label_args <- formalArgs(cyto_plot_labeller)
    label_args <- label_args[grepl("label_", label_args)]
    label_args <- label_args[-match("label_text", label_args)]
    lapply(label_args, function(z){
      args[[z]] <<- rep(args[[z]], length.out = PNS)
    })
    
    # REPEAT ARGUMENTS - WITHOUT GATE
  }else{
    
    # REPEAT LABEL ARGS PER LABEL_TEXT
    label_args <- formalArgs(cyto_plot_labeller)
    label_args <- label_args[grepl("label_", label_args)]
    label_args <- label_args[-match("label_text", label_args)]
    lapply(label_args, function(z){
      args[[z]] <<- rep(args[[z]], length.out = length(label_text))
    })
    
  }
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # STATISTICS PER POPULATION --------------------------------------------------
  
  # POPULATIONS
  pops <- .cyto_label_pops(list(x),
                           gate = gate,
                           negate = negate)
  
  # COMPUTE STATISTICS
  label_stat <- .cyto_label_stat(list(x), 
                                 pops = pops,
                                 channels = channels,
                                 label_stat = label_stat,
                                 axes_trans = trans,
                                 gate = gate,
                                 hist_smooth = hist_smooth)
  
  # PREPARE LABEL_TEXT ---------------------------------------------------------
  
  # MERGE LABEL_TEXT & LABEL_STAT
  label_text <- .cyto_label_text(label_text,
                                 label_stat)
  
  
  # LABEL CONSTRUCTION ---------------------------------------------------------
  
  # MANUAL CO-ORDINATES - SUPPLIED OR SELECTED
  label_text_xy <- mapply(
    function(label_text,
             label_text_x,
             label_text_y,
             label_text_font,
             label_text_size,
             label_text_col,
             label_fill,
             label_fill_alpha) {
      
      # SELECT MISSING CO-ORDINATES
      if (any(.all_na(c(label_text_x, label_text_y)))) {
        if(grepl("\n", label_text)){
          txt <- paste(strsplit(label_text, "\n")[[1]], collapse = " ")
        }else{
          txt <- label_text
        }
        message(paste("Select a location on the plot the position the",
                      txt, "label."))
        label_text_xy <- locator(n = 1)
        label_text_x <- label_text_xy[[1]]
        label_text_y <- label_text_xy[[2]]
      }
      
      # PLOT LABELS
      .boxed.labels(
        x = label_text_x,
        y = label_text_y,
        labels = label_text,
        border = FALSE,
        font = label_text_font,
        cex = label_text_size,
        col = label_text_col,
        bg = label_fill,
        alpha.bg = label_fill_alpha
      )
      
      # RETURN LABEL CO-ORDINATES
      return(c(label_text_x, label_text_y))
      
    },
    label_text,
    label_text_x,
    label_text_y,
    label_text_font,
    label_text_size,
    label_text_col,
    label_fill,
    label_fill_alpha,
    SIMPLIFY = FALSE
  )
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL CO-ORDINATES IN MATRIX
  label_text_xy <- do.call("rbind", label_text_xy)
  colnames(label_text_xy) <- c("x","y")
  label_text_xy <- as.matrix(label_text_xy)
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)
  
}

## CYTO_PLOT_LABELLER ----------------------------------------------------------

#' Add labels to existing cyto_plot
#'
#' Convenient labelling function to add prepared text labels to an existing
#' cyto_plot.
#'
#' @param label_text character string to include in the label.
#' @param label_text_x vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param label_text_y vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param label_text_font integer [1,2,3,4] passed to \code{text} to alter the
#'   font, set to \code{2} by default for a bold font.
#' @param label_text_size numeric character expansion used to control the size
#'   of the text in the labels, set to \code{0.8} by default. See \code{?text}
#'   for details.
#' @param label_text_col specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param label_fill fill colour to use for labels, set to "white" by default.
#' @param label_fill_alpha numeric [0,1] controls the transparency of the fill
#'   colour, set to \code{0.6} by default.
#' @param ... not in use.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @export
cyto_plot_labeller <- function(label_text = NA,
                               label_text_x = NA,
                               label_text_y = NA,
                               label_text_font = 2,
                               label_text_size = 0.8,
                               label_text_col = "black",
                               label_fill = "white",
                               label_fill_alpha = 0.6,
                               ...){
  
  # LABEL CONSTRUCTION ---------------------------------------------------------
  
  # MANUAL CO-ORDINATES - SUPPLIED OR SELECTED
  label_text_xy <- mapply(
    function(label_text,
             label_text_x,
             label_text_y,
             label_text_font,
             label_text_size,
             label_text_col,
             label_fill,
             label_fill_alpha) {
      
      # MISSING LABEL
      if(!.all_na(label_text)){
        
        # SELECT MISSING CO-ORDINATES
        if (any(.all_na(c(label_text_x, label_text_y)))) {
          if(grepl("\n", label_text)){
            txt <- paste(strsplit(label_text, "\n")[[1]], collapse = " ")
          }else{
            txt <- label_text
          }
          message(paste("Select a location on the plot the position the",
                        txt, "label."))
          label_text_xy <- locator(n = 1)
          label_text_x <- label_text_xy[[1]]
          label_text_y <- label_text_xy[[2]]
        }
        
        # PLOT LABELS
        .boxed.labels(
          x = label_text_x,
          y = label_text_y,
          labels = label_text,
          border = FALSE,
          font = label_text_font,
          cex = label_text_size,
          col = label_text_col,
          bg = label_fill,
          alpha.bg = label_fill_alpha
        )
        
      }
      
      # RETURN LABEL CO-ORDINATES
      return(c(label_text_x, label_text_y))
      
    },
    label_text,
    label_text_x,
    label_text_y,
    label_text_font,
    label_text_size,
    label_text_col,
    label_fill,
    label_fill_alpha,
    SIMPLIFY = FALSE
  )
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL CO-ORDINATES IN MATRIX
  label_text_xy <- do.call("rbind", label_text_xy)
  colnames(label_text_xy) <- c("x","y")
  label_text_xy <- as.matrix(label_text_xy)
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)
  
}
