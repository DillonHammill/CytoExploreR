## CYTO_PLOT INTERNAL FUNCTIONS ------------------------------------------------

# These functions are used within cyto_plot() to prepare arguments for use in
# cyto_plot().

## BUILD -----------------------------------------------------------------------

#' Build cyto_plot from arguments
#' 
#' @param args list of arguments prepared by cyto_plot.
#' 
#' @noRd
.cyto_plot_build <- function(args){
  
  # ARGUMENTS - CYTO_PLOT CLASS
  class(args) <- "cyto_plot"
  
  # POPUP
  if(getOption("cyto_plot_method") == "flowFrame") {
    cyto_plot_new(args$popup)
  }
  
  # CYTO_PLOT_EMPTY
  cyto_plot_empty(args)
  
  # HISTOGRAMS
  if(length(args$channels) == 1) {
    cyto_plot_hist(args)
    # POINTS & CONTOURS
  } else {
    cyto_plot_point(args)
  }
  
  # AUTOMATIC LABEL POSITION - BYPASS ON SAVING - MUST BE AFTER CYTO_PLOT_EMPTY
  if(!getOption("cyto_plot_save") &
     args$label == TRUE &
     args$label_position == "auto") {
    label_text_xy <- do.call(".cyto_plot_label_coords", args)
    args$label_text_x <- label_text_xy[, "x"]
    args$label_text_y <- label_text_xy[, "y"]
  }
  
  # GATES & LABELS
  if (!.all_na(args$gate)) {
    # PLOT GATE & ASSOCIATED LABELS
    label_text_xy <- .cyto_plot_gate_label(args)
    args$label_text_x <- label_text_xy[, "x"]
    args$label_text_y <- label_text_xy[, "y"]
    # LABELS
  } else {
    if (args$label == TRUE) {
      label_text_xy <- do.call("cyto_plot_labeller", args)
      args$label_text_x <- label_text_xy[, "x"]
      args$label_text_y <- label_text_xy[, "y"]
    }
  }
  
  return(args)
  
}

## DATA ------------------------------------------------------------------------

#' Prepare data for cyto_plot
#' @return list of singular cytosets for each plot
#' @importFrom flowWorkspace cytoset
#' @importFrom purrr transpose
#' @noRd
.cyto_plot_data <- function(x,
                            overlay = NULL,
                            merge_by = "name",
                            select = NULL,
                            display = 50000,
                            barcode = FALSE,
                            seed = 56) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CYTOSETS ONLY
  if(!cyto_class(x, "flowSet")) {
    stop("cyto_plot only supports cytoset objects!")
  }
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(x, select)
  }
  
  # GROUP
  x <- cyto_group_by(x, merge_by)
  
  # OVERLAY
  if(!is.null(overlay)) {
    # OVERLAY LIST
    if(!cyto_class(overlay, "list", TRUE)) {
      overlay <-  list(overlay)
    }
    # OVERLAY - SELECT & GROUP
    overlay <- lapply(seq_along(overlay), function(z){
      cs <- overlay[[z]]
      # CHECK
      if(!cyto_class(cs, "flowSet")) {
        stop("cyto_plot only supports cytosets!")
      }
      # SELECT - IF POSSIBLE
      cs <- tryCatch(
        cyto_select(cs, select),
        error = function(e){return(cs)}
      )
      # GROUP - ONLY IF SAME GROUPS
      if(setequal(cyto_groups(cs, merge_by), names(x))) {
        cs <- cyto_group_by(cs, merge_by)
      } else {
        cs <- structure(
          rep(list(cs), length(x)),
          names  = names(x)
        )
      }
      # LAYER
      return(cs)
    })
    # COMBINE & TRANSPOSE
    x <- transpose(c(list(x), overlay))
  }
  
  # SAMPLE & MERGE
  x <- structure(
    # LOOP THROUGH EACH PLOT
    lapply(seq_along(x), function(z){
      # PLOT LAYERS - LIST OF CYTOSETS
      cs_list <- x[[z]]
      # GROUP NAME
      grp <- names(x)[z]
      # CHECK EACH LAYER FOR EVENTS IN PREVIOUS LAYERS - CANNOT STORE GROUP NAME
      i <- lapply(seq_along(cs_list) , function(w){
        # CYTOSET
        cs <- cs_list[[w]]
        # IDENTIFIERS
        ids <- cyto_names(cs)
        # CYTOFRAME
        structure(
          lapply(seq_along(cs), function(q) {
            # IDENTIFIER
            id <- ids[q]
            l <- list(
              "total" = nrow(cs[[q]]),
              "layer" = NULL,
              "match" = NULL,
              "sample" = NULL
            )
            for(v in seq_len(w-1)) {
              # FIRST LAYER
              if(v == 0){
                break()
                # SUBSEQUENT LAYERS
              } else {
                # SAMPLE ID MATCH PREVIOUS LAYER
                if(id %in% cyto_names(cs_list[[v]])) {
                  # MATCHING EVENT IDS?
                  m <- sum(
                    cyto_exprs(cs[[q]], "Event-ID") %in%
                      cyto_exprs(cs_list[[v]][[match(id, ids)]], "Event-ID")
                  )
                  # ANY MATCHING EVENTS
                  if(m > 0) {
                    # UPDATE LIST
                    l$layer <- v
                    l$match <- m
                    break()
                  }
                }
              }
            }
            return(l)
          }),
          names = ids
        )
      })
      
      # SAMPLE - SAVE TO CS_SUB_LIST - NEED CS_LIST CANNOT UPDATE IN PLACE
      cs_sub_list <- list()
      lapply(seq_along(cs_list), function(w){
      
        # CYTOSET
        cs <- cs_list[[w]]
        # IDENTIFIERS
        ids <- cyto_names(cs)
        
        # ALL CYTOFRAMES CONTAIN NEW EVENTS
        if(is.null(LAPPLY(i[[w]], `[[`, "layer"))) {
          # COMPUTE SAMPLE SIZES
          n <- cyto_sample_n(cs, 
                             display = display)
          for(id in ids) {
            i[[w]][[id]]$sample <<- n[id]
          }
          cs_sub_list[[w]] <<- cyto_sample(cs,
                                           display = n,
                                           seed = seed)
        # SOME/ALL CYTOFRAMES CONTAIN EVENTS ON PLOT
        } else {
          # FIND WHICH CYTOFRAMES HAVE EVENTS ON PLOT
          m <- LAPPLY(i[[w]], function(s){
            !is.null(s$layer)
          })
          # STORE CYTOFRAMES IN LIST
          cf_list <- list()
          # CYTOFRAMES WITH EVENTS ON PLOT
          for(id in names(m[m])) {
            # MATCHING LAYER
            l <- i[[w]][[id]]$layer
            # ORIGINAL SIZE OF REFERNCE LAYER
            t <- i[[l]][[id]]$total
            # SAMPLE SIZE OF REFERENCE LAYER
            s <- i[[l]][[id]]$sample
            # CYTOFRAME SUBSET OF PREVIOUS LAYER
            if(i[[w]][[id]]$match == nrow(cs[[id]])){
              # USE EVENT IDs FROM REFERENCE LAYER
              cf_list[[id]] <- cs[[id]][
                cyto_exprs(cs[[id]], "Event-ID") %in% 
                  cyto_exprs(cs_sub_list[[l]][[id]], "Event-ID") 
                ,]
              i[[w]][[id]]$sample <<- nrow(cf_list[[id]])
            # CYTOFRAME OVERLAP WITH PREVIOUS LAYER
            } else {
              # EVENT IDs FROM LAYER SAMPLE
              e <- cyto_exprs(cs[[id]], "Event-ID")[
                cyto_exprs(cs[[id]], "Event-ID") %in%
                  cyto_exprs(cs_sub_list[[l]][[id]], "Event-ID")
              ]
              # NON-OVERLAPPING PORTION OF CYTOFRAME
              cf_new <- cs[[id]][
                !cyto_exprs(cs[[id]], "Event-ID") %in%
                  cyto_exprs(cs_list[[l]][[id]], "Event-ID"),
              ]
              n <- nrow(cyto_exprs(cf_new))
              # SAMPLE NON-OVERLAPPING PORTION OF CYTOFRAME
              e <- c(e,
                     cyto_exprs(cf_new)[
                       sample(
                         nrow(cf_new),
                         round(
                           (n/i[[l]][[id]]$total)*(i[[l]][[id]]$sample)
                           )    
                         )
                       , "Event-ID"])
              # COMPLETE SAMPLE CYTOFRAME
              cf_list[[id]] <- cs[[id]][
                cyto_exprs(cs[[id]], "Event-ID") %in% e
              , ]
              # STORE SAMPLE SIZE
              i[[w]][[id]]$sample <<- length(e)
            }
          }
          # CYTOFRAMES WITHOUT EVENTS ON PLOT
          # NEW LAYER - ALL NEW CYTOFRAMES
          if(length(m[!m]) == length(cs)) {
            # COMPUTE SAMPLE SIZES
            n <- cyto_sample_n(cs, 
                               display = display)
            # SAMPLE
            for(id in names(m[!m])) {
              i[[w]][[id]]$sample <<- n[id]
              cf_list[[id]] <- cyto_sample(cs[[id]],
                                           display = n[id],
                                           seed = seed)
            }
          # SOME NEW CYTOFRAME(S)
          } else if(length(m[!m]) != 0 & length(m[!m]) < length(cs)) {
            # RATIO SAMPLING TO COMPUTE N FOR NEW CYTOFRAMES
            n <- round(
              sum(LAPPLY(i[[w]][names(m[m])], `[[`, "sample")) *
              (sum(LAPPLY(i[[w]][names(m[!m])], `[[`, "total")) /
                 sum(LAPPLY(i[[w]][names(m[m])], `[[`, "total")))
            )
            # COMPUTE SAMPLE SIZES FOR NEW CYTOFRAMES
            n <- cyto_sample_n(cs[names(m[!m])],
                               display = n)
            # SAMPLE
            for(id in names(m[!m])) {
              i[[w]][[id]]$sample <<- n[id]
              cf_list[[id]] <- cyto_sample(cs[[id]],
                                           display = n[id])
            }
          }
          # ADD CYTOFRAMES TO NEW CYTOSET
          cs_sub_list[[w]] <<- cytoset(cf_list)
        }
      })
      
      # COERCE
      lapply(seq_along(cs_sub_list), function(r){
        if(length(cs_sub_list[[r]]) > 1) {
          cytoset(
            structure(
              list(
                as(cs_sub_list[[r]], "cytoframe")
              ),
              names = grp
            )
          )
        } else {
          return(cs_sub_list[[r]])
        }
      })
    }),
    names = names(x)
  )
  
  # DATA
  return(x)
  
}

## AXES ------------------------------------------------------------------------

## AXES LIMITS ----

#' Compute axes limits for cyto_plot
#'
#' The lower limit is always set to zero unless there is data below this limit.
#' The upper limit is determined by the axes_limits argument.
#'
#' @param x cytometry object(s) which require range calculation.
#' @param parent name of parent population to extract for range calculation.
#' @param channels name(s) of channel(s).
#' @param axes_limits either "auto", "data" or "machine". "auto" use data limits but
#'   always includes zero.
#' @param plot logical indicating whether a check should be performed for
#'   channel length.
#' @param buffer fraction indicating the amount of buffering to be added on top
#'   of the upper and lower limit, set to 0.03 by default.
#' @param anchor logical indicating if the lower limit should be anchored to
#'   zero if the data range is above zero, set to TRUE by default.
#'
#' @return vector containing minimum and maximum values.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_axes_limits <- function(x,
                                   parent = NULL,
                                   channels = NA,
                                   axes_limits = "auto",
                                   plot = FALSE,
                                   buffer = 0.04,
                                   anchor = TRUE){
  
  # COMPATIBILITY --------------------------------------------------------------
  
  # FLOWCORE USES INSTRUMENT
  if(axes_limits == "instrument"){
    axes_limits <- "machine"
  }
  
  # DATA -----------------------------------------------------------------------
  
  # LIST CYTOFRAMES/CYTOSETS
  if(cyto_class(x, "list")){
    x <- unlist(x) 
    x <- lapply(x, function(z){
      # GATINGHIERARCHY/GATINGSET
      if(cyto_class(z, "GatingSet")){
        z <- cyto_data_extract(z, 
                               parent = parent)[[1]]
      }
      return(z)
    })
  # GATINGHIERARCHY/GATINGSET
  }else if(cyto_class(x, "GatingSet")){
    x <- list(cyto_data_extract(x, 
                                parent = parent)[[1]])
  # CYTOFRAME/CYTOSET
  } else {
    x <- list(x)
  }
  
  # CHANNELS -------------------------------------------------------------------
  
  # MARKERS TO CHANNELS
  if(!.all_na(channels)){
    channels <- cyto_channels_extract(x[[1]], channels, plot)
  }else{
    channels <- cyto_channels(x[[1]])
  }
  
  # TIME CHANNEL ALWAYS USE DATA LIMITS
  if(grepl("^Time", channels, ignore.case = TRUE)){
    axes_limits <- "data"
  } 
  
  # RANGE ----------------------------------------------------------------------
  
  # DATA RANGE
  data_range <- lapply(x, function(z){
    if(cyto_class(z, "flowFrame")){
      if(nrow(z) == 0){
        type <- "instrument"
      }else{
        type <- "data"
      }
      rng <- suppressWarnings(
        range(z[, channels],
              type = type)
      )
    }else if(cyto_class(z, "flowSet")){
      rng <- suppressWarnings(
        cyto_apply(z, function(y){
          if(nrow(y) == 0){
            type <- "instrument"
          }else{
            type <- "data"
          }
          range(y, type = type)
        }, 
        input = "cytoframe", 
        channels = channels,
        inverse = FALSE)
      )
      rng <- do.call("rbind", rng)
    }
    return(rng)
  })
  data_range <- do.call("rbind", data_range)
  
  # MIN/MAX DATA RANGE
  data_range <- lapply(seq_len(ncol(data_range)), function(z){
    c(min(data_range[, z]), max(data_range[, z]))
  })
  data_range <- do.call("cbind", data_range)
  colnames(data_range) <- channels
  rownames(data_range) <- c("min", "max")
  
  # MACHINE RANGE
  if(axes_limits == "machine"){
    machine_range <- lapply(x, function(z){
      if(cyto_class(z, "flowFrame")){
        rng <- suppressWarnings(
          range(z,
                type = "instrument")[, channels, drop = FALSE]
        )
      }else if(cyto_class(z, "flowSet")){
        rng <- suppressWarnings(
          cyto_apply(z,
                     "range",
                     type = "instrument",
                     input = "cytoframe",
                     channels = channels,
                     inverse = FALSE)
        )
        rng <- do.call("rbind", rng)
      }
      return(rng)
    })
    machine_range <- do.call("rbind", machine_range)
    
    # MIN/MAX DATA RANGE
    machine_range <- lapply(seq_len(ncol(machine_range)), function(z){
      c(min(machine_range[, z]), max(machine_range[, z]))
    })
    machine_range <- do.call("cbind", machine_range)
    colnames(machine_range) <- channels
    rownames(machine_range) <- c("min", "max")
    
    # REPLACE MAX DATA RANGE 
    data_range["max", ] <- machine_range["max", ]
  }
  
  # REPLACE LOWER LIMIT IF > 0 - AUTO
  if(axes_limits != "data" & anchor == TRUE){
    if(any(data_range[1,] > 0)){
      data_range[1, data_range[1,] > 0] <- 0
    }
  }
  
  # BUFFER ---------------------------------------------------------------------
  
  # ADD 2% BUFFER EITHER SIDE - 4% TOTAL IN PLOT
  if(buffer != 0){
    # BUFFER
    for(z in channels) {
      # RANGE
      RNG <- data_range[, z][2] - data_range[, z][1]
      # ADD BUFFER
      data_range[, z][1] <- data_range[, z][1] - buffer * RNG
      data_range[, z][2] <- data_range[, z][2] + buffer * RNG
    }
  }
  
  return(data_range)
  
}

## AXES TEXT ----

#' Axes ticks and text
#'
#' @param x list of cytoframes/cytosets.
#' @param channels name(s) of the channel(s) used to construct the plot.
#' @param axes_trans transformerList.
#' @param axes_range named list of axes limits for each each axis (i.e.
#'   list(xlim,ylim)).
#' @param axes_limits either "auto", "data" or "machine".
#'
#' @return list containing axis labels and breaks.
#'
#' @importFrom methods is
#' @importFrom grDevices .axisPars
#'
#' @noRd
.cyto_plot_axes_text <- function(x,
                                 channels,
                                 axes_trans = NA,
                                 axes_range = list(NA, NA),
                                 axes_limits = "data") {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNELS
  channels <- rep(c(channels, rep(NA, 2)), length.out = 2) # x & y
  
  # AXES_TRANS
  if(!.all_na(axes_trans) & !cyto_class(axes_trans, "transformerList", TRUE)) {
    stop("'axes_trans' must be an object of class transformerList.")
  }
  
  # PAD SUPPLIED AXES RANGES ---------------------------------------------------
  
  # AXES RANGE - 4% BUFFER
  axes_range <- structure(
    lapply(axes_range, function(z){
      if(!.all_na(z)) {
        # RANGE
        rng <- z[2] - z[1]
        pad <- (rng - rng / 1.04) / 2
        return(c(z[1] - pad,
                 z[2] + pad))
      } else {
        return(z)
      }
    }), 
    names = names(axes_range))
  
  # AXES TICKS & LABELS --------------------------------------------------------
  
  # LOOP THROUGH CHANNELS
  axes_text <- lapply(channels, function(z){
    # 1D Y AXIS
    if(.all_na(z)) {
      return(NA)
    }
    # LINEAR CHANNEL
    if(.all_na(axes_trans) | !z %in% names(axes_trans)) {
      # COMPUTE AXP
      axp <- unlist(.axisPars(axes_range[[z]]))
      # COMPUTE AXES LABELS
      axis_breaks <- seq(axp[1], axp[2], (axp[2] - axp[1])/axp[3])
      # AXES MINOR INTERVAL
      axis_interval <- (axis_breaks[2] - axis_breaks[1])/axp[3]
      # COMPUTE MAJOR TICK LOCATIONS
      axis_ticks_min <- axp[1] - floor(
        (axp[1] - axes_range[[z]][1])/axis_interval
      ) * axis_interval
      axis_ticks_max <- axp[2] + floor(
        (axes_range[[z]][2] - axp[2])/axis_interval
      ) * axis_interval
      axis_ticks <- seq(axis_ticks_min, axis_ticks_max, axis_interval)
      axis_labels <- unlist(lapply(axis_ticks, function(z){
        if(z %in% axis_breaks){
          return(z)
        } else {
          return("")
        }
      }))
      axis_text <- list("label" = axis_labels,
                        "at" = axis_ticks)
      # TRANSFORMED CHANNEL  
    } else {
      # AXIS TICKS
      axis_ticks <- c(
        sort(LAPPLY(
          1 * 10^seq_len(9),
          function(z) {
            -seq(90, 10, -10) * z
          }
        )),
        seq(-9, 9, 1),
        LAPPLY(
          1 * 10^seq_len(9),
          function(z) {
            seq(10, 90, 10) * z
          }
        )
      )
      # AXES LABELS
      axis_labels <- lapply(axis_ticks, function(z) {
        if (z != 0) {
          pwr <- log10(abs(z))
          if (z < 0) {
            pwr <- -pwr
          }
        }
        if (z == 0) {
          quote(0)
        } else if (pwr == 0) {
          quote("")
        } else if (abs(pwr) %% 1 == 0) {
          substitute(10^pwr)
        } else {
          quote("")
        }
      })
      axis_labels <- do.call("expression", axis_labels)
      # AXES TEXT
      trans_func <- axes_trans[[z]]$transform
      inv_func <- axes_trans[[z]]$inverse
      # AXES RANGE ON LINEAR SCALE
      if (!.all_na(axes_range[[z]])) {
        rng <- inv_func(axes_range[[z]])
      } else {
        rng <- inv_func(.cyto_plot_axes_limits(x,
                                               channels = z,
                                               axes_limits = axes_limits
        )[, z])
      }
      # RESTRICT AXES TICKS & LABELS BY RANGE
      ind <- which(axis_ticks > rng[1] & axis_ticks < rng[2])
      axis_ticks <- axis_ticks[ind]
      axis_labels <- axis_labels[ind]
      # BREAKS - TRANSFORMED SCALE
      axes_breaks <- signif(trans_func(axis_ticks))
      # BREAKS & LABELS
      return(list("label" = axis_labels, "at" = axes_breaks))
    }
  })
  names(axes_text) <- channels
  
  # RETURN PREPARED AXES_TEXT
  return(axes_text)
  
}

## AXES LABELS ----

#' Get axes titles for cyto_plot
#'
#' @param x list of cytoframes/cytosets.
#' @param channels used to construct the plot.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param hist_stat "percent", "count" or "density".
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_axes_label <- function(x,
                                  channels,
                                  xlab,
                                  ylab,
                                  hist_stat = "count") {
  
  # CHANNELS
  channels <- c(channels, NA)[1:2]
  
  # MARKERS
  markers <- cyto_markers(x[[1]])
  
  # AXES_LABELS
  cnt <- 0
  axes_labels <- lapply(channels, function(z) {
    assign("cnt", cnt + 1, envir = parent.frame(2))
    # LABELS
    if (cnt == 1) {
      lab <- xlab
    } else {
      lab <- ylab
    }
    # XLAB OR SCATTER YLAB
    if (!.all_na(z)) {
      # XLAB/YLAB
      if (missing(lab) | .empty(lab)) {
        # MARKER ASSIGNED - DIFFERENT FROM CHANNEL
        ind <- match_ind(z, names(markers))
        if (length(ind) > 0 & !z %in% markers) {
          return(paste(markers[ind], z))
        } else {
          return(z)
        }
      } else if (.all_na(lab)) {
        return(NA)
      }
      # HISTOGRAM YLAB
    } else {
      if (missing(lab) | .empty(lab)) {
        # COUNT
        if (hist_stat == "count") {
          return("Count")
          # PERCENT
        } else if (hist_stat == "percent") {
          return("% of Mode")
          # DENSITY
        } else if (hist_stat == "density") {
          return("Density")
        }
      } else if (.all_na(lab)) {
        return(NA)
      }
    }
  })
  names(axes_labels) <- c("xlab", "ylab")
  
  # RETURN AXES LABELS
  return(axes_labels)
}

## LAYOUT/MARGINS --------------------------------------------------------------

## LAYOUT ----

#' Set plot layout
#'
#' @param x list of  cytoframe/cytoset lists.
#' @param layout grid dimensions c(nr, nc), NA or FALSE.
#'
#' @importFrom grDevices n2mfrow
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_layout <- function(x,
                              layout = NULL) {
  
  # LAYOUT
  if (is.null(layout) | .empty(layout)) {
    if(is.vector(x) & length(x) == 1) {
      layout <- rev(n2mfrow(x))
    }else {
      layout <- rev(n2mfrow(length(x)))
    }
  }
  return(layout)
  
}

## MARGINS ----

#' Set plot margins
#'
#' @param x list of flowFrames or density objects to plot.
#' @param legend logical indicating whether a legend should be included in the
#'   plot.
#' @param title if NULL remove excess space above plot.
#' @param axes_text vector of logicals indicating whether the x and y axes
#'   should be included on the plot.
#' @param margins a vector of length 4 to control the margins around the bottom,
#'   left, top and right of the plot, set to c(NA, NA, NA, NA) by default to let
#'   `cyto_plot` compute optimal margins.
#' @param point_col required to make space for point_col_scale legend.
#'
#' @importFrom methods is
#'
#' @noRd
.cyto_plot_margins <- function(x,
                               channels = NULL,
                               legend = FALSE,
                               legend_text = NA,
                               legend_text_size = 1,
                               title,
                               axes_text = list(TRUE, TRUE),
                               margins = c(NA, NA, NA, NA),
                               point_col = NA) {
  
  # ARGUMENTS
  args <- .args_list()
  
  # MARGINS
  if (is.null(margins) | .empty(margins) | .all_na(margins)) {
    
    # DEFAULT
    mar <- c(5.1, 5.1, 4.1, 2.1)
    
    # POINT_COL_SCALE
    if(length(channels) == 2 & 
       (any(is.na(point_col)) | 
        any(point_col %in% c(cyto_channels(x), cyto_markers(x))))) {
      mar[4] <- mar[4] + 1
    }
    
    # LEGEND TEXT
    if (length(x) > 1 &
        legend != FALSE &
        !.all_na(legend_text)) {
      mar[4] <- 7 + max(nchar(legend_text)) * 0.32 * mean(legend_text_size)
    }
    
    # TITLE
    if (.all_na(title)) {
      mar[3] <- 2.1
    }
    
    # X AXIS
    if (!all(is(axes_text[[1]], "list"))) {
      if (.all_na(axes_text[[1]])) {
        # NA == FALSE returns NA not T/F
      } else if (all(axes_text[[1]] == FALSE)) {
        mar[1] <- 4.1
      }
    }
    
    # Y AXIS
    if (!all(is(axes_text[[2]], "list"))) {
      if (.all_na(axes_text[[2]])) {
        # NA == FALSE return NA not T/F
      } else if (all(axes_text[[2]] == FALSE)) {
        mar[2] <- 4.1
      }
    }
  } else {
    if (length(margins) != 4) {
      stop("'margins' must be a vector with 4 elements.")
    }
    mar <- margins
  }
  
  # SET MAR PARAMETER
  par("mar" = mar)
}

## LABELS ----------------------------------------------------------------------

## LABEL_TEXT & LABEL_STAT ARGUMENTS ----

#' Prepare label_text and label_stat arguments
#' @noRd
.cyto_plot_label_args <- function(x,
                                  channels = NULL,
                                  label = "",
                                  label_text = "",
                                  label_stat = "",
                                  gate = NA,
                                  negate = FALSE,
                                  hist_stack = 0,
                                  ...) {
  
  
  if(any(is.function(label_stat))) {
    stop("label_stat must be a character string!")
  }
  
  # POPULATIONS PER LAYER
  NP <- .cyto_gate_count(gate, negate = negate)
  
  # TOTAL POPULATIONS
  SMP <- length(x)
  TNP <- NP * SMP
  TNP_split <- split(seq_len(TNP), rep(seq_len(SMP), each = NP))
  
  # LABEL_TEXT
  if (all(LAPPLY(label_text, ".empty"))) {
    label_text <- rep(NA, TNP)
  } else {
    label_text <- rep(c(label_text, rep(NA, TNP)), length.out = TNP)
  }
  
  # LABEL_STAT FILLED WITH EMPTY
  if (!all(LAPPLY(label_stat, ".empty")) &
      any(LAPPLY(label_stat, ".empty"))) {
    label_stat[LAPPLY(label_stat, ".empty")] <- NA
  }
  
  # LABEL_STAT
  # 1D PLOT NO STACK
  if (length(channels) == 1 & hist_stack == 0) {
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
  } else if (length(channels) == 1 & hist_stack != 0) {
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
  
  # RETURN LABEL ARGUMENTS
  return(list("label" = label,
              "label_text" = label_text,
              "label_stat" = label_stat))
  
}

## LABEL POPULATIONS ----

#' Get a list of gated populations to label
#'
#' @param x list of cytoframes/cytosets per layer to be gated
#' @param gate list of gate objects to apply to each element of x, gates only
#'   required for base layer.
#' @param negate logical indicating whether negated population should be
#'   included.
#' @param label_stat required to identify which layers require gating to get
#'   statistics
#' @param ... not in use.
#'
#' @return list of flowFrames lists
#'
#' @importFrom flowWorkspace flowSet_to_cytoset
#' @importFrom flowCore Subset split quadGate
#' @importFrom methods is
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_pops <- function(x,
                                  gate,
                                  negate = FALSE,
                                  label_stat = NA,
                                  ...) {
  
  # FLOWFRAME LIST
  if(cyto_class(x, c("flowFrame", "flowSet"))){
    x <- list(x)
  }
  
  # NO GATES -------------------------------------------------------------------
  
  # RETURN X
  if (.all_na(gate)) {
    pops <- lapply(seq_len(length(x)), function(z){
      x[z]
    })
    return(pops)
  }
  
  # PREPARE GATES --------------------------------------------------------------
  
  # LIST OF GATE OBJECT LISTS
  if(cyto_class(gate, "list") &
     all(LAPPLY(gate, "cyto_class", "list", TRUE))){
    # USE BASE LAYER GATES
    gate <- gate[[1]]
  }
  
  # LIST OF GATE OBJECTS
  if (cyto_class(gate, "list", TRUE)) {
    if (all(LAPPLY(gate, 
                   "cyto_class", 
                   c("rectangleGate",
                     "polygonGate",
                     "ellipsoidGate",
                     "quadGate",
                     "filters"), 
                   TRUE))) {
      gate <- unlist(gate)
    }
  } else if (cyto_class(gate, "filters", TRUE)) {
    gate <- unlist(gate)
  } else if (cyto_class(gate, 
                        c("rectangleGate",
                          "polygonGate",
                          "ellipsoidGate",
                          "quadGate",
                          "filters"), 
                        TRUE)) {
    gate <- list(gate)
  }
  
  # List of RectangleGates to QuadGate (MUST BE 4 RECTANGLES)
  if (length(gate) == 4 & all(LAPPLY(gate, function(z) {
    cyto_class(z, "rectangleGate") & any(grepl("quad", names(attributes(z))))
  }))) {
    # CHANNELS
    chans <- as.character(parameters(gate[[1]]))
    quad_order <- LAPPLY(gate, function(z) {
      z@filterId
    })
    gate <- list(.cyto_gate_quad_convert(gate, channels = chans))
  }
  
  # GATES ----------------------------------------------------------------------
  
  # NEGATED GATE - QUADGATES EXCLUDED
  if (negate == TRUE & !any(LAPPLY(gate, function(z) {
    is(z, "quadGate")
  }))) {
    if (length(gate) > 1) {
      gate <- c(gate, list(do.call("|", gate)))
    } else {
      gate <- c(gate, gate)
    }
  }
  
  # LABEL_STAT -----------------------------------------------------------------
  
  # LAYERS
  layers <- length(x)
  
  # LABELS PER LAYER (SAME LENGTH AS GATE WHICH INCLUDES NEGATE)
  labels_per_layer <- split(label_stat,
                            rep(seq_len(layers), 
                                each = length(label_stat)/layers))
  
  # POPULATIONS ----------------------------------------------------------------
  
  # CAREFUL CANNOT NEGATE INDIVIDUAL QUADRANTS EITHER
  
  # ARGUMENTS
  args <- .args_list()
  
  # GATING PER LAYER (list of pops per layer)
  pops <- lapply(seq_along(x), function(z){
    cs <- x[[z]]
    pops_to_label <- labels_per_layer[[z]]
    gated <- c()
    gated_pops <- lapply(seq_len(length(gate)), function(w) {
      # NO GATE OR NO STAT
      if (.all_na(gate[[w]]) | is.na(pops_to_label[[w]])) {
        assign("gated", c(gated, FALSE), envir = parent.frame(2))
        res <- cs
      }else{
        assign("gated", c(gated, TRUE), envir = parent.frame(2))
      }
      # NEGATED POPULATION
      if (negate == TRUE & w == length(gate)) {
        res <- split(cs, gate[[w]])[[2]]
        # *** CYTOSET CONVERSION ***
        if(cyto_class(res, "flowSet", TRUE)) {
          res <- flowSet_to_cytoset(res)
        }
        # GATED POPULATIONS
      } else {
        # QUADGATES RETURN MULTIPLE POPULATIONS
        if (is(gate[[w]], "quadGate")) {
          if ("quad_order" %in% names(args)) {
            quads <- unlist(strsplit(gate[[w]]@filterId, "\\|"))
            res <- split(cs, gate[[w]])[c(2, 1, 3, 4)][match(
              quad_order,
              quads
            )] # FIX ORDER

          } else {
            res <- split(cs, gate[[w]])[c(2, 1, 3, 4)] # FIX ORDER
          }
          # *** CYTOSET CONVERSION ***
          res <- structure(
            lapply(res, function(q){
              if(cyto_class(q, "flowSet", TRUE)){
                flowSet_to_cytoset(q)
              }
            }),
            names = names(res)
          )
          # SINGLE POPULATIONS
        } else {
          # RECTANGLE BELONGS TO QUADGATE
          if (is(gate[[w]], "rectangleGate") &
              any(grepl("quad", names(attributes(gate[[w]]))))) {
            q <- names(attributes(gate[[w]])[["quadrants"]])
            coords <- .cyto_gate_coords(gate[w],
                                        channels = as.character(
                                          parameters(gate[[w]]))
            )
            chans <- colnames(coords)
            coords <- lapply(colnames(coords), function(y) {
              unique(coords[, y][is.finite(coords[, y])])
            })
            names(coords) <- chans
            qg <- quadGate(
              filterId = paste(q, collapse = "|"),
              .gate = coords
            )
            p <- split(cs, qg)[c(2, 1, 3, 4)] # FIX ORDER
            names(p) <- q
            # *** CYTOSET CONVERSION ***
            p <- structure(
              lapply(p, function(b){
                if(cyto_class(b, "flowSet", TRUE)){
                  flowSet_to_cytoset(b)
                }
              }),
              names = names(p)
            )
            p[[match(gate[[w]]@filterId, names(p))]]
          } else {
            Subset(cs, gate[[w]])
          }
        }
      }
      return(res)
    })
    # SIGNAL IF GATE PRESENT
    if(any(gated == TRUE)){
      gated[gated == TRUE] <- "gated"
    }
    if(any(gated == FALSE)){
      gated[gated == FALSE] <- "gated"
    }
    names(gated_pops) <- gated
    return(gated_pops)
  })
  
  # RETURN LIST OF GATED POPULATIONS
  return(pops)
}

## LABEL STATISTICS ----

#' Compute and prepare statistics for labels
#'
#' @param x list of parental flowFrame objects.
#' @param pops list of population lists to label.
#' @param channels vector of channels used to construct the plot.
#' @param label_stat names of statistics to include in labels, supplied per
#'   layer.
#'
#' @return computed statistics to include in labels.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_stat <- function(x,
                                  pops,
                                  channels,
                                  axes_trans = NA,
                                  label_stat,
                                  gate = NA,
                                  hist_smooth = 1,
                                  hist_bins = 256,
                                  ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # NO LABEL_STAT
  if (.all_na(label_stat)) {
    return(label_stat)
  }
  
  # FLOWFRAME LIST
  if(is(x, "flowFrame")){
    x <- list(x)
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # LABEL_STAT
  label_stat <- split(label_stat,
                      rep(seq_len(length(x)),
                          each  = length(pops[[1]])))
  
  # COMPUTE STAT AGAINST BASE LAYER (NO GATES)
  if(.all_na(gate)){
    x <- rep(x[1], length(x))
  }
  
  # COMPUTE STATISTICS ---------------------------------------------------------
  
  # STATISTICS
  LABEL_STAT <- LAPPLY(seq_len(length(x)), function(z) {
    # LABEL_STAT
    ST <- lapply(seq_len(length(pops[[z]])), function(y) {
      # STATISTIC SUPPLIED
      if (!.all_na(label_stat[[z]][y])) {
        # FREQUENCY
        if(grepl(paste0("freq", "$"), 
                 label_stat[[z]][y], ignore.case = TRUE)) {
          res <- cyto_stats_compute(x[[z]],
                                    channels = channels,
                                    gate = pops[[z]][[y]],
                                    stat = label_stat[[z]][y],
                                    inverse = TRUE,
                                    trans = axes_trans,
                                    round = 2,
                                    input = "matrix",
                                    smooth = hist_smooth,
                                    bins = hist_bins,
                                    markers = FALSE)
          # OTHER STATISTICS - GATED POPS
        } else {
          res <- cyto_stats_compute(pops[[z]][[y]],
                                    channels = channels,
                                    stat = label_stat[[z]][y],
                                    inverse = TRUE,
                                    trans = axes_trans,
                                    round = 2,
                                    input = "matrix",
                                    smooth = hist_smooth,
                                    bins = hist_bins,
                                    markers = FALSE)
        }
        # DROP NAME COLUMN
        res <- res[, -match("name", colnames(res)), drop = FALSE]
        # EXTRACT STATISTICS - EITHER SINGLE OR DOUBLE
        res <- res[1, , drop = TRUE]
        # ROUND 2 DECIMAL PLACES
        if(!grepl("count$", label_stat[[z]][y])) {
          res <- .round(res)
        }
        # STATSTICS REQUIRE %
        if(any(LAPPLY(c("freq",
                        "percent",
                        "cv",
                        "rcv"), function(w){
                          grepl(paste0(w, "$"), 
                                label_stat[[z]][y], 
                                ignore.case = TRUE)
                        }))) {
          res <- paste(res, "%")
        }
        # STATISTIC PER CHANNEL - 2D PLOT
        if(length(res) > 1) {
          res <- paste("x =", res[1], "\n", "y =", res[2])
        }
      } else {
        res <- NA
      }
      return(res)
    })
    return(ST)
  })
  
  # RETURN COMPUTED STATISTICS -------------------------------------------------
  return(LABEL_STAT)
}

## PREPARE LABEL_TEXT ----

#' Prepare text for labels to include statistics
#'
#' @param label_text text to include in the labels.
#' @param label_stat vector of computed statistics to include in the labels.
#'
#' @return vector of finalised labels incorporating both label_text and
#'   label_stat.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_text <- function(label_text,
                                  label_stat,
                                  ...) {
  
  # MERGE LABEL_TEXT & LABEL_STAT
  for(z in seq_len(length(label_text))) {
    # NO LABEL_TEXT & NO LABEL_STAT
    if (.all_na(label_text[z]) & .all_na(label_stat[z])) {
      label_text[z] <- NA
      # NO LABEL_TEXT & LABEL_STAT
    } else if (.all_na(label_text[z]) & !.all_na(label_stat[z])) {
      label_text[z] <- label_stat[z]
      # LABEL_TEXT & NO LABEL_STAT
    } else if (!.all_na(label_text[z]) & .all_na(label_stat[z])) {
      label_text[z] <- label_text[z]
      # LABEL_TEXT & LABEL_STAT
    } else if (!.all_na(label_text[z]) & !.all_na(label_stat[z])) {
      label_text[z] <- paste(label_text[z], label_stat[z], sep = "\n")
    }
  }
  
  # RETURN LABEL_TEXT
  return(label_text)
}

## LABEL CO-ORDINATES ----

#' Compute offset label co-ordinates
#'
#' Used internally within cyto_plot to compute offset label co-ordinates. Only
#' called if label_position and label are set to TRUE.
#'
#' @param args list of named cyto_plot arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' @noRd
.cyto_plot_label_coords <- function(x,
                                    channels = NULL,
                                    pops = NULL,
                                    axes_limits = c(NA, NA),
                                    gate = NA,
                                    label_text = NA,
                                    label_text_x = NA,
                                    label_text_y = NA,
                                    label_text_size = 1.8,
                                    hist_stack = 0.5,
                                    hist_layers = NA,
                                    hist_smooth = 1,
                                    hist_bins = 256,
                                    d = NULL,
                                    ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # X - CYTO_PLOT ARGUMENTS
  if(class(x) == "cyto_plot") {
    .args_update(x)
  }
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # PLOT LIMITS
  lims <- par("usr")
  
  # X LIMITS
  xmin <- lims[1]
  xmax <- lims[2]
  xrng <- xmax - xmin
  xpad <- (xrng - xrng / 1.04) / 2 # 2% BUFFER
  xmin <- xmin + xpad
  xmax <- xmax - xpad
  xrng <- xmax - xmin
  
  # Y LIMITS
  ymin <- lims[3]
  ymax <- lims[4]
  yrng <- ymax - ymin
  ypad <- (yrng - yrng / 1.04) / 2 # 2% BUFFER
  ymin <- ymin + ypad
  ymax <- ymax - ypad
  yrng <- ymax - ymin
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # LAYERS & POPULATIONS
  L <- length(pops)
  GNP <- length(pops[[1]])
  GC <- .cyto_gate_count(gate, total = TRUE)
  
  # PREPARE LABEL ARGUMENTS PER LAYER
  label_text <- split(label_text, rep(seq_len(L), each = GNP))
  label_text_x <- split(label_text_x, rep(seq_len(L), each = GNP))
  label_text_y <- split(label_text_y, rep(seq_len(L), each = GNP))
  label_text_size <- split(label_text_size, rep(seq_len(L), each = GNP))
  
  # hist_stack --------------------------------------------------------------
  
  # 1D PLOT
  if (length(channels) == 1) {
    # STACKING
    if (hist_stack != 0 &
        ifelse(.all_na(hist_layers), TRUE, hist_layers != 1)) {
      stk <- LAPPLY(names(d), function(z){
        as.numeric(unlist(strsplit(z, "-"))[1])
      })
      y_coords <- stk + 0.5 * (
        as.numeric(unlist(strsplit(names(d)[2], "-"))[1]) - 
          as.numeric(unlist(strsplit(names(d)[1], "-"))[1])
      )
      # REPEAT (Y_COORDS/LAYER)
      y_coords <- rep(y_coords, each = GNP)
      # NO STACKING
    } else {
      # DEFAULT Y COORD - 50% Y RANGE
      y_coords <- rep(0.5 * yrng, L * GNP)
    }
    y_coords <- split(y_coords, rep(seq_len(L), each = GNP))
  }
  
  # COMPUTE LABEL CO-ORDINATES -------------------------------------------------
  
  # GATE CENTERS - (IGNORE SUPPLIED LABEL COORDS)
  if (!.all_na(gate)) {
    gate_centers <- .cyto_gate_center(gate,
                                      channels = channels,
                                      text_x = rep(NA, GC),
                                      text_y = rep(NA, GC)
    )
  }
  
  # LABEL_TEXT_XY - MATRIX
  label_text_xy <- lapply(seq_len(L), function(z) {
    # TEMPORARY STORAGE VECTORS
    text_x <- c()
    text_y <- c()
    # COMPUTE LABEL CO-ORDINATES
    for(y in seq_len(GNP)) {
      # LABEL
      if (!.all_na(label_text[[z]][y])) {
        # ID PLOT - CENTER OF RANGE
        if (length(channels) == 1) {
          # GATE
          if (!.all_na(gate)) {
            # GATED POPULATION
            if (y <= nrow(gate_centers)) {
              # X COORD - GATE CENTER
              if (.all_na(label_text_x[[z]][y])) {
                text_x[y] <- gate_centers[y, "x"]
              } else {
                text_x[y] <- label_text_x[[z]][y]
              }
              # Y COORD - STACKS/LIMITS
              if (.all_na(label_text_y[[z]][y])) {
                text_y[y] <- y_coords[[z]][y]
              }else if(!.all_na(label_text_y[[z]][y])) {
                text_y[y] <- label_text_y[[z]][y]
              }
              # NEGATED POPULATION
            } else if (y > nrow(gate_centers)) {
              # X COORD - RANGE CENTER
              if (.all_na(label_text_x[[z]][y])) {
                # NO EVENTS
                if (cyto_stat_count(pops[[z]][[y]]) == 0) {
                  # RANGE CENTER - PLOT LIMITS
                  text_x[y] <- mean(c(xmin, xmax))
                } else {
                  # RANGE
                  rng <- suppressMessages(
                    .cyto_plot_axes_limits(pops[[z]][[y]],
                                           channels = channels[1],
                                           axes_limits = axes_limits,
                                           buffer = 0,
                                           anchor = FALSE
                    )[, channels]
                  )
                  # UNIMODAL - 50% RANGE
                  if (abs(diff(rng)) <= 0.6 * xrng) {
                    text_x[y] <- quantile(rng, 0.5)
                    # UMULTIMODAL - 56% RANGE
                  } else {
                    text_x[y] <- quantile(rng, 0.56)
                  }
                }
                # X COORD MANUALLY SUPPLIED
              } else if (!.all_na(label_text_x[[z]][y])) {
                text_x[y] <- label_text_x[[z]][y]
              }
              # Y COORD - STACKS/LIMITS
              if (.all_na(label_text_y[[z]][y])) {
                text_y[y] <- y_coords[[z]][y]
                # Y COORD MANUALLY SUPPLIED
              } else if (!.all_na(label_text_y[[z]][y])) {
                text_y[y] <- label_text_y[[z]][y]
              }
            }
            # NO GATE
          } else if (.all_na(gate)) {
            # X COORD - RANGE CENTER
            if (.all_na(label_text_x[[z]][y])) {
              # NO EVENTS
              if (cyto_stat_count(pops[[z]][[y]]) == 0) {
                # RANGE CENTER - PLOT LIMITS
                text_x[y] <- mean(c(xmin, xmax))
              } else {
                # RANGE
                rng <- suppressMessages(
                  .cyto_plot_axes_limits(pops[[z]][[y]],
                                         channels = channels[1],
                                         axes_limits = axes_limits,
                                         buffer = 0,
                                         anchor = FALSE
                  )[, channels]
                )
                # UNIMODAL - 50% RANGE
                if (abs(diff(rng)) <= 0.6 * xrng) {
                  text_x[y] <- quantile(rng, 0.5)
                  # MULTIMODAL - 56% RANGE
                } else {
                  text_x[y] <- quantile(rng, 0.56)
                }
              }
              # X COORD MANUALLY SUPPLIED
            } else if (!.all_na(label_text_x[[z]][y])) {
              text_x[y] <- label_text_x[[z]][y]
            }
            # Y COORD - STACKS/LIMITS
            if (.all_na(label_text_y[[z]][y])) {
              text_y[y] <- y_coords[[z]][y]
              # Y COORD MANUALLY SUPPLIED
            } else if (!.all_na(label_text_y[[z]][y])) {
              text_y[y] <- label_text_y[[z]][y]
            }
          }
          # 2D PLOT - MODE
        } else if (length(channels) == 2) {
          # GATE
          if (!.all_na(gate)) {
            # GATED POPULATION
            if (y <= nrow(gate_centers)) {
              # X COORD - GATE CENTER
              if (.all_na(label_text_x[[z]][y])) {
                text_x[y] <- gate_centers[y, "x"]
              } else {
                text_x[y] <- label_text_x[[z]][y]
              }
              # Y COORD - GATE CENTER
              if (.all_na(label_text_y[[z]][y])) {
                text_y[y] <- gate_centers[y, "y"]
              } else {
                text_y[y] <- label_text_y[[z]][y]
              }
              # NEGATED POPULATION
            } else if (y > nrow(gate_centers)) {
              # X COORD - MODE/RANGE CENTER
              if (.all_na(label_text_x[[z]][y])) {
                # NO EVENTS
                if (cyto_stat_count(pops[[z]][[y]]) < 2) {
                  # RANGE CENTER
                  text_x[y] <- mean(c(xmin, xmax))
                } else {
                  # MODE
                  text_x[y] <- cyto_apply(pops[[z]][[y]],
                                          "cyto_stat_mode",
                                          input = "matrix",
                                          channels = channels[1],
                                          hist_smooth = hist_smooth,
                                          hist_bins = hist_bins,
                                          inverse = FALSE)[, 1]
                }
                # X COORD MANUALLY SUPPLIED
              } else if (!.all_na(label_text_x[[z]][y])) {
                text_x[y] <- label_text_x[[z]][y]
              }
              # Y COORD - MODE/RANGE CENTER
              if (.all_na(label_text_y[[z]][y])) {
                # NO EVENTS
                if (cyto_stat_count(pops[[z]][[y]]) == 0) {
                  # RANGE CENTER
                  text_y[y] <- mean(c(ymin, ymax))
                } else {
                  # MODE
                  text_y[y] <- cyto_apply(pops[[z]][[y]],
                                          "cyto_stat_mode",
                                          input = "matrix",
                                          channels = channels[2],
                                          hist_smooth = hist_smooth,
                                          hist_bins = hist_bins,
                                          inverse = FALSE)[, 1]
                }
                # Y COORD MANUALLY SUPPLIED
              } else if (!.all_na(label_text_y[[z]][y])) {
                text_y[y] <- label_text_y[[z]][y]
              }
            }
            # NO GATE
          } else if (.all_na(gate)) {
            # X COORD - MODE/RANGE CENTER
            if (.all_na(label_text_x[[z]][y])) {
              # NO EVENTS
              if (cyto_stat_count(pops[[z]][[y]]) < 2) {
                # RANGE CENTER
                text_x[y] <- mean(c(xmin, xmax))
              } else {
                # MODE
                text_x[y] <- cyto_apply(pops[[z]][[y]],
                                        "cyto_stat_mode",
                                        input = "matrix",
                                        channels = channels[1],
                                        hist_smooth = hist_smooth,
                                        hist_bins = hist_bins,
                                        inverse = FALSE)[, 1]
              }
              # X COORD SUPPLIED MANUALLY
            } else if (!.all_na(label_text_x[[z]][y])) {
              text_x[y] <- label_text_x[[z]][y]
            }
            # Y COORD - MODE
            if (.all_na(label_text_y[[z]][y])) {
              # NO EVENTS
              if (cyto_stat_count(pops[[z]][[y]]) < 2) {
                text_y[y] <- mean(c(ymin, ymax))
              } else {
                # MODE
                text_y[y] <- cyto_apply(pops[[z]][[y]],
                                        "cyto_stat_mode",
                                        input = "matrix",
                                        channels = channels[2],
                                        hist_smooth = hist_smooth,
                                        hist_bins = hist_bins,
                                        inverse = FALSE)[, 1]
              }
              # Y COORD MANUALLY SUPPLIED
            } else if (!.all_na(label_text_y[[z]][y])) {
              text_y[y] <- label_text_y[[z]][y]
            }
          }
        }
        # NO LABEL
      } else if (.all_na(label_text[[z]][y])) {
        text_x[y] <- NA
        text_y[y] <- NA
      }
    }
    
    # MATRIX
    text_xy <- matrix(c(text_x, text_y),
                      ncol = 2,
                      nrow = GNP,
                      byrow = FALSE
    )
    colnames(text_xy) <- c("x", "y")
    return(text_xy)
  })
  
  # UPDATE LABEL_COORDS
  label_text_x <- lapply(label_text_xy, function(z){
    unlist(z[, "x", drop = TRUE])
  })
  label_text_y <- lapply(label_text_xy, function(z){
    unlist(z[, "y", drop = TRUE])
  })
  
  # OFFSET LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL DIMENSIONS
  label_dims <- lapply(seq_len(L), function(z){
    lapply(seq_len(GNP), function(y){
      # COMPUTE LABEL DIMENSIONS
      if (!.all_na(label_text[[z]][y])) {
        .cyto_plot_label_dims(
          label_text = label_text[[z]][y],
          label_text_x = label_text_x[[z]][y],
          label_text_y = label_text_y[[z]][y],
          label_text_size = label_text_size[[z]][y]
        )
      } else {
        matrix(rep(NA, 4),
               ncol = 2,
               dimnames = list(NULL, c("x", "y"))
        )
      }
    })
  })
  
  # OFFSET BY LAYER
  if (length(channels) == 1 & hist_stack != 0) {
    # OFFSET PER LAYER
    for(z in seq_len(L)) {
      # LABEL OVERLAP
      if (.cyto_plot_label_overlap(label_dims[[z]])) {
        # LABEL HEIGHT - OFFSETTING
        label_height <- max(LAPPLY(label_dims[[z]], function(y) {
          max(y[, "y"]) - min(y[, "y"])
        }), na.rm = TRUE)
        # LABEL HEIGHT BUFFERING
        label_height <- 1.18 * label_height
        # Y COORDS TO OFFSET
        text_y <- label_text_y[[z]]
        # OFFSET Y CO-ORDINATES (EXCLUDE NA)
        label_text_y[[z]][!is.na(text_y)] <- tryCatch(
          {
            .suppress_all_messages(
              .spread.labels(text_y[!is.na(text_y)],
                             mindiff = label_height,
                             min = ymin,
                             max = ymax
              )
            )
          },
          error = function(e) {
            return(text_y[!is.na(text_y)])
          }
        )
      }
    }
    # OFFSET ALL LABELS
  } else {
    # COLLAPSE LABEL_DIMS
    label_dims <- unlist(label_dims, recursive = FALSE)
    # LABEL OVERLAP
    if (.cyto_plot_label_overlap(label_dims)) {
      # LABEL HEIGHT - OFFSETTING
      label_height <- max(LAPPLY(label_dims, function(y) {
        max(y[, "y"]) - min(y[, "y"])
      }), na.rm = TRUE)
      # LABEL HEIGHT BUFFERING
      label_height <- 1.18 * label_height
      # OFFSET Y CO-ORDINATES (EXCLUDE NA)
      text_y <- unlist(label_text_y)
      text_y[!is.na(text_y)] <-
        .suppress_all_messages(
          .spread.labels(text_y[!is.na(text_y)],
                         mindiff = label_height,
                         min = ymin,
                         max = ymax
          )
        )
      # UPDATE LABEL_TEXT_Y
      label_text_y <- split(text_y,
                            rep(seq_len(L), each = GNP))
    }
  }
  
  # RETURN COMPUTED LABEL CO-ORDINATES -----------------------------------------
  
  # LABEL CO-ORDINATE MATRIX
  label_text_xy <- matrix(c(unlist(label_text_x), 
                            unlist(label_text_y)),
                          ncol = 2,
                          byrow = FALSE
  )
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
}

## LABEL DIMENSIONS ----

#' Compute label dimensions
#'
#' Labels should already contain statistic as well. Co-ordiante for each label
#' must already be computed.
#'
#' @importFrom graphics strwidth strheight
#'
#' @return upper left and bottom right x and y co-ordinates of labels
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_dims <- function(label_text,
                                  label_text_x,
                                  label_text_y = NA,
                                  label_text_size = 1,
                                  ...) {
  
  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # RESET CEX & XPD
  old_pars <- par(c("cex", "xpd"))
  on.exit(par(old_pars))
  
  # SET CEX & XPD
  par(cex = label_text_size)
  par(xpd = TRUE)
  
  # LABEL DIMENSIONS -----------------------------------------------------------
  
  # ARGUMENTS
  xpad <- 1.2
  ypad <- 1.2
  adj <- 0.5
  
  # BOX ADJUSTMENT
  box_adj <- adj + (xpad - 1) * label_text_size * (0.5 - adj)
  
  # BOX DIMENSIONS
  lwidths <- strwidth(label_text)
  rwidths <- lwidths * (1 - box_adj)
  lwidths <- lwidths * box_adj
  bheights <- theights <- strheight(label_text) * 0.5
  
  # BOX X COORDS
  xr <- label_text_x - lwidths * xpad
  xl <- label_text_x + lwidths * xpad
  
  # BOX Y COORDS
  yb <- label_text_y - bheights * ypad
  yt <- label_text_y + theights * ypad
  
  # LABEL DIMENSIONS MATRIX ----------------------------------------------------
  
  # MATRIX - TOP LEFT THEN BOTTOM RIGHT
  coords <- matrix(c(
    min(c(xl, xr)),
    max(c(yb, yt)),
    max(c(xl, xr)),
    min(c(yb, yt))
  ),
  ncol = 2,
  byrow = TRUE
  )
  colnames(coords) <- c("x", "y")
  
  # RETURN LABEL DIMENSIONS ----------------------------------------------------
  return(coords)
}

## LABEL OVERLAP ----

#' Check if any labels are overlapping.
#'
#' @param x list of label dimensions.
#'
#' @return TRUE or FALSE based on whether any overlapping labels are found.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label_overlap <- function(x) {
  
  # For each rectangle in x
  overlaps <- LAPPLY(seq_len(length(x)), function(y) {
    
    # Check if other rectangles overlap
    LAPPLY(seq_len(length(x))[-y], function(z) {
      
      # Co-ordinates of reference label
      x1 <- x[[y]][, "x"]
      y1 <- x[[y]][, "y"]
      
      # Co-ordinates of comparison label
      x2 <- x[[z]][, "x"]
      y2 <- x[[z]][, "y"]
      
      # MISSING LABELS - NO OVERLAP
      if (any(is.na(c(x1, x2, y1, y2)))) {
        return(FALSE)
      }
      
      # X co-ordinates are overlapping
      if (min(x2) >= min(x1) & min(x2) <= max(x1) |
          max(x2) >= min(x1) & max(x2) <= max(x1)) {
        # Y co-ordinates are also overlapping
        if (min(y2) >= min(y1) & min(y2) <= max(y1) |
            max(y2) >= min(y1) & max(y2) <= max(y1)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
      
      # Non-overlapping x and y co-ordinates
      return(FALSE)
    })
  })
  
  # RETURN TRUE OR FALSE
  return(any(overlaps))
}

## GATE LABELS -----------------------------------------------------------------

#' Add labels to gates as they are plotted
#' 
#' @param args list of cyto_plot arguments
#'
#' @noRd
.cyto_plot_gate_label <- function(args) {
  
  # GENERAL --------------------------------------------------------------------
  
  # SAMPLES
  SMP <- length(args$x)
  
  # POPULATIONS PER LAYER
  NP <- length(args$pops[[1]])
  
  # TOTAL POPULATIONS
  TNP <- length(args$x) * NP
  
  # GATES
  NG <- length(args$gate)
  
  # POPULATIONS PER GATE
  P <- .cyto_gate_count(args$gate,
                        negate = FALSE,
                        total = FALSE)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # CYTO_PLOT_GATE & CYTO_PLOT_LABELLER ARGUMENTS
  gate_args <- formalArgs("cyto_plot_gate.list")
  label_args <- formalArgs("cyto_plot_labeller")
  
  # RESTRICT ARGUMENTS
  args <- args[names(args) %in% c(gate_args, "label", label_args)]
  
  # SPLIT GATE_FILL ARGUMENTS BY POPULATIONS
  for(z in gate_args[which(grepl("gate_fill", gate_args))]){
    if(z %in% names(args)) {
      args[[z]] <- split(args[[z]], rep(seq_len(NG), times = P))
    }
  }
  
  # SPLIT LABEL ARGUMENTS BY POPULATIONS PER LAYER
  for(z in label_args) {
    if(z %in% names(args)) {
      args[[z]] <- rep(args[[z]], length.out = TNP)
      args[[z]] <- split(args[[z]], rep(seq_len(SMP), each = NP))
      # RE-ARRANGE LABEL COORDS PER GATE
      args[[z]] <- lapply(seq_len(NP), function(y){
        LAPPLY(args[[z]], `[[`, y)
      })
    }
  }
  
  # GATE ARGUMENTS
  gate_args <- args[names(args) %in% gate_args]
  
  # LABEL ARGUMENTS
  label_args <- args[names(args) %in% label_args]
  
  # GATE & ASSOCIATED LABELS ---------------------------------------------------
  
  # GATES & LABELS
  label_text_xy <- lapply(seq_len(NP), function(z) {
    # GATED POPULATION(S) - GATE & LABEL(S)
    if (z <= NG) {
      # PLOT GATE
      do.call(
        "cyto_plot_gate",
        c(
          list("channels" = gate_args$channels),
          lapply(
            gate_args[!grepl("channels", names(gate_args))],
            `[[`, z
          )
        )
      )
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (args$label == TRUE &
          !.all_na(args$label_text[[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(args$label_text[z])),
                          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
      # NEGATED POPULATION(S) - LABEL(s) ONLY
    } else {
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (args$label == TRUE &
          !.all_na(args$label_text[[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(args$label_text[[z]])),
                          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
    }
    return(text_xy)
  })
  label_text_xy <- do.call("rbind", label_text_xy)
  
  # RE-ARRANGE LABEL ARGUMENTS -------------------------------------------------
  
  # UPDATE LABEL_TEXT_X & LABEL_TEXT_Y
  args$label_text_x <- label_text_xy[, "x"]
  args$label_text_y <- label_text_xy[, "y"]
  
  # REVERT LABEL_TEXT_X & LABEL_TEXT_Y TO ORIGINAL FORMAT
  if (SMP > 1) {
    args$label_text_x <- LAPPLY(seq_len(SMP), function(z) {
      args$label_text_x[names(args$label_text_x) == z]
    })
    args$label_text_y <- LAPPLY(seq_len(SMP), function(z) {
      args$label_text_y[names(args$label_text_y) == z]
    })
  }
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL_COORDS MATRIX
  label_text_xy <- matrix(c(args$label_text_x, args$label_text_y),
                          ncol = 2,
                          byrow = FALSE
  )
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
  
}

## LEGENDS ---------------------------------------------------------------------

## LEGEND ----

#' Create a legend for cyto_plot
#'
#' \code{.cyto_plot_margins} will handle setting the plot margins to make space
#' for the legend.
#'
#' @param x list of flowFrame objects to include in the plot.
#' @param channels name of the channels or markers to be used to construct the
#'   plot.
#' @param legend logical indicating whether a legend should be included for
#'   plots including overlays, set to FALSE by default.
#' @param legend_text vector of labels to use for the legend.
#' @param legend_text_font numeric indicating the font to use for legend text,
#'   set to 2 for bold font by default. See \code{\link[graphics:par]{?par}}
#'   font for details.
#' @param legend_text_size character expansion for legend text, set to 1 by
#'   default.
#' @param legend_text_col colour to use for legend text, set to "black by
#'   default.
#' @param legend_line_col vector of line colours to use for legend.
#' @param legend_box_fill vector of fill colours to use for legend.
#' @param legend_point_col vector of colours to use for points in legend.
#' @param hist_cols vector colours to draw from when selecting density fill
#'   colours if none are supplied to hist_fill.
#' @param hist_fill colour(s) used to fill polygons.
#' @param hist_fill_alpha numeric [0,1] used to control fill transparency,
#'   set to 1 by default to remove transparency.
#' @param hist_line_type line type(s) to use for border(s), set to solid
#'   lines by default.
#' @param hist_line_width line width for border.
#' @param hist_line_col colour(s) for border line, set to "black" by default.
#' @param point_shape point character to use for points, set to "." by default
#'   to maximise plotting speed.
#' @param point_size numeric specifying the degree of character expansion for
#'   points, set to 2 by default.
#' @param point_col_scale vector of colours to use for density gradient.
#' @param point_cols vector colours to draw from when selecting colours for
#'   points if none are supplied to point_col.
#' @param point_col colours to use for points, set to NA by default to blue-red
#'   density colour scale.
#' @param point_alpha numeric [0,1] used to control colour transparency, set to
#'   1 by default to remove transparency.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom graphics legend strheight
#' @importFrom grDevices adjustcolor
#'
#' @noRd
.cyto_plot_legend <- function(x,
                              channels,
                              legend = "fill",
                              legend_text = NA,
                              legend_text_font = 1,
                              legend_text_size = 1,
                              legend_text_col = "black",
                              legend_line_type = NA,
                              legend_line_width = NA,
                              legend_line_col = NA,
                              legend_box_fill = NA,
                              legend_point_col = NA,
                              hist_cols = NA,
                              hist_fill = NA,
                              hist_fill_alpha = 1,
                              hist_line_type = 1,
                              hist_line_width = 1,
                              hist_line_col = "black",
                              point_shape = ".",
                              point_size = 2,
                              point_col_scale = NA,
                              point_cols = NA,
                              point_col = NA,
                              point_col_alpha = 1) {
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # ARGUMENTS
  args <- .args_list()
  
  # CYTO_PLOT_THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # LEGEND_TEXT ----------------------------------------------------------------
  
  # Estimate legend height using strheight
  lgnd <- paste(legend_text, collapse = " \n ")
  lgnd_height <- strheight(lgnd,
                           cex = legend_text_size,
                           font = legend_text_font
  )
  
  # LEGEND POSITION ------------------------------------------------------------
  
  # Calculate y center of plot
  cnt <- par("usr")[3] + (par("usr")[4] - par("usr")[3]) / 2
  
  # Legend for 1D density distributions
  if (length(channels) == 1) {
    
    # Set default legend type to fill
    if (legend == TRUE) {
      legend <- "fill"
    }
    
    # Reverse legend text order for legend
    legend_text <- rev(legend_text)
    
    # Line legend
    if (legend == "line") {
      
      # Revert to hist_line_col if no colours supplied
      if (.all_na(legend_line_col)) {
        legend_line_col <- hist_line_col
      }
      
      # Revert to hist_line_type if not specified
      if (.all_na(legend_line_type)) {
        legend_line_type <- hist_line_type
      }
      
      # Revert to hist_line_width if not specified
      if (.all_na(legend_line_width)) {
        legend_line_width <- hist_line_width
      }
      
      # Construct legend
      legend(
        x = 1.07 * par("usr")[2],
        y = cnt + 0.52 * lgnd_height,
        legend = legend_text,
        text.font = rev(legend_text_font),
        cex = legend_text_size,
        text.col = rev(legend_text_col),
        col = rev(legend_line_col),
        lty = rev(legend_line_type),
        lwd = rev(legend_line_width),
        xpd = TRUE,
        bty = "n",
        x.intersp = 0.5
      )
      # Fill legend
    } else if (legend == "fill") {
      
      # COLOURS
      hist_fill <- .cyto_plot_hist_fill(x,
                                        hist_fill = hist_fill,
                                        hist_cols = hist_cols,
                                        hist_fill_alpha = 1
      )
      
      # Revert to hist_fill if no legend fill colours supplied
      if (.all_na(legend_box_fill)) {
        legend_box_fill <- hist_fill
      }
      # Alpha adjust colours if suppplied directly to legend_box_fill
      if (!.all_na(legend_box_fill) &
          !all(hist_fill_alpha == 1)) {
        legend_box_fill <- mapply(
          function(legend_box_fill,
                   hist_fill_alpha) {
            adjustcolor(legend_box_fill, hist_fill_alpha)
          }, legend_box_fill, hist_fill_alpha
        )
      }
      
      # Construct legend
      legend(
        x = 1.07 * par("usr")[2],
        y = cnt + 0.52 * lgnd_height,
        legend = legend_text,
        fill = rev(legend_box_fill),
        xpd = TRUE,
        bty = "n",
        x.intersp = 0.5,
        cex = legend_text_size,
        text.col = rev(legend_text_col),
        text.font = rev(legend_text_font)
      )
    }
    
    # Legend for 2D scatter plot
  } else if (length(channels) == 2) {
    
    # CYTO_PLOT_POINT_COL_SCALE
    point_col_scale <- .cyto_plot_point_col_scale(point_col_scale)
    
    # Prepare point_col - alpha adjust later
    point_col <- .cyto_plot_point_col(x,
                                      channels = channels,
                                      point_col_scale = point_col_scale,
                                      point_cols = point_cols,
                                      point_col = point_col,
                                      point_col_alpha = 1
    )
    
    # Prepare point col - use first density colour
    point_col <- LAPPLY(point_col, function(z) {
      if (length(z) > 1) {
        return(point_col_scale[1])
      } else {
        return(z)
      }
    })
    
    # Revert to point_col if no legend point cols supplied
    if (.all_na(legend_point_col)) {
      legend_point_col <- point_col
    }
    
    # Alpha adjust colours supplied directly to legend_point_col
    if (!.all_na(legend_point_col) &
        !all(point_col_alpha == 1)) {
      legend_point_col <- mapply(function(col, alpha) {
        adjustcolor(col, alpha)
      }, legend_point_col, point_col_alpha)
    }
    
    legend(
      x = 1.08 * par("usr")[2],
      y = cnt + 0.6 * lgnd_height,
      legend = rev(legend_text),
      col = rev(legend_point_col),
      pch = rev(point_shape),
      pt.cex = rev(2 * point_size),
      xpd = TRUE,
      bty = "n",
      x.intersp = 0.7,
      cex = legend_text_size,
      text.col = rev(legend_text_col),
      text.font = rev(legend_text_font)
    )
  }
}

## TITLE -----------------------------------------------------------------------

#' Title for cyto_plot
#'
#' @param x flowFrame object.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_title <- function(x,
                             channels,
                             overlay = NA,
                             title = "") {
  
  # x can be a list
  if (class(x) == "list") {
    if (length(x) > 1) {
      overlay <- x[2:length(x)]
      x <- x[[1]]
    }
  }
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Update arguments
  .args_update(args)
  
  # 1D density distributions
  if (length(channels) == 1) {
    
    # missing/empty replace with valid title
    if (.empty(title)) {
      
      # stacked/overlays lack a title
      if (.all_na(overlay)) {
        title <- cyto_names(x)
      } else {
        title <- NA
      }
      
      # NA will remove title in cyto_plot_empty
    } else if (.all_na(title)) {
      title <- NA
    }
    
    # 2D scatterplots
  } else if (length(channels) == 2) {
    
    # missing title replaced with sample name
    if (.empty(title)) {
      title <- cyto_names(x)
      # NA will remove title in cyto_plot_empty
    } else if (.all_na(title)) {
      title <- NA
    }
  }
  
  return(title)
}

## HEADER ----------------------------------------------------------------------

#' Add header to plot
#' 
#' @importFrom graphics mtext
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_header <- function(header,
                              header_text_font = 2,
                              header_text_size = 1,
                              header_text_col = "black") {
  
  # HEADER - OMA MUST BE PRESET - c(0,0,3,0)
  mtext(header, 
        outer = TRUE, 
        font = header_text_font,
        cex = header_text_size,
        col = header_text_col)
  
}

## POINTS ----------------------------------------------------------------------

## POINT COLOUR SCALE ----

#' Get density gradient colours for cyto_plot
#'
#' @param point_col_scale vector of ordered colours to use for point density
#'   colour scale.
#'
#' @return a list of colorRampPalette functions to be used in densCols.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_point_col_scale <- function(point_col_scale = NA) {
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Use default colour scale
  if (.all_na(args[["point_col_scale"]])) {
    args[["point_col_scale"]] <- .cyto_plot_colour_palette("point_col_scale")
  }
  
  return(args[["point_col_scale"]])
}

## POINT COLOUR ----

#' Get point colours for cyto_plot
#'
#' @param x list of flowFrames.
#' @param channels used to construct the plot.
#' @param point_col_scale vector of colours to use for density gradient.
#' @param point_cols vector colours to select from when choosing a colour for
#'   each layer in x.
#' @param point_col vector of length x indicating colours to use for each layer.
#'   If NA set to default density gradient.
#' @param point_col_alpha transparency to use for point colours.
#'
#' @importFrom grDevices densCols colorRampPalette adjustcolor
#' @importFrom flowCore exprs
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_point_col <- function(x,
                                 channels,
                                 point_col_scale,
                                 point_cols,
                                 point_col,
                                 point_col_alpha = 1) {
  
  # Expected number of colours
  SMP <- length(x)
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme - possibly remove?
  args <- .cyto_plot_theme_inherit(args)
  
  # Update arguments
  .args_update(args)
  
  # No colours supplied for density gradient
  if (.all_na(point_col_scale)) {
    point_col_scale <- .cyto_plot_colour_palette(type = "point_col_scale")
  }
  
  # Make colorRampPalette
  if (class(point_col_scale) != "function") {
    col_scale <- colorRampPalette(point_col_scale)
  } else {
    col_scale <- point_col_scale
  }
  
  # No colours supplied for selection
  if (.all_na(point_cols)) {
    point_cols <- .cyto_plot_colour_palette(type = "point_cols")
  }
  
  # Make colorRampPalette
  if (class(point_cols) != "function") {
    cols <- colorRampPalette(point_cols)
  } else {
    cols <- point_cols
  }
  
  # Repeat point_col arguments SMP times
  point_col <- rep(point_col, length.out = SMP)
  point_col_alpha <- rep(point_col_alpha, length.out = SMP)
  
  # Convert point_col to list
  if (!is(point_col, "list")) {
    point_col <- lapply(seq(1, SMP), function(z) {
      point_col[z]
    })
  }
  
  # First layer contains density gradient if no other colour is designated
  if (all(LAPPLY(point_col, ".all_na"))) {
    
    # Extract data
    fr_exprs <- exprs(x[[1]])[, channels]
    
    # Too few events for density computation
    if (!is.null(nrow(fr_exprs))) {
      if (nrow(fr_exprs) >= 2) {
        # Get density colour for each point
        point_col[[1]] <- suppressWarnings(
          densCols(fr_exprs,
                   colramp = col_scale
          )
        )
      }
    } else {
      point_col[[1]] <- point_col_scale[1]
    }
  }
  
  # Remaining colours are selected one per layer from point_cols
  if (any(LAPPLY(point_col, ".all_na"))) {
    
    # Number of layers missing colours
    n <- length(point_col[LAPPLY(point_col, ".all_na")])
    
    # Pull colours out of point_cols
    clrs <- cols(n)
    
    # Replace NA values in point_col with selected colours
    point_col[LAPPLY(point_col, ".all_na")] <- clrs
  }
  
  # RANGE CALIBRATION
  cyto_cal <- .cyto_calibrate_recall()
  
  # 1D COLOUR SCALE
  point_col <- lapply(point_col, function(z) {
    if (length(z) == 1) {
      # NAME OF CHANNEL/MARKER
      if (z %in% c(
        cyto_channels(x[[1]]),
        cyto_markers(x[[1]])
      )) {
        # CONVERT TO CHANNEL
        z <- cyto_channels_extract(
          x[[1]],
          z
        )
        # MATRIX
        fr_exprs <- exprs(x[[1]])
        # CALIBRATION
        if (!is.null(cyto_cal)) {
          if (z %in% colnames(cyto_cal)) {
            cyto_range <- c(
              min(cyto_cal[, z]),
              max(cyto_cal[, z])
            )
          } else {
            cyto_range <- c(
              min(fr_exprs[, z]),
              max(fr_exprs[, z])
            )
          }
        } else {
          cyto_range <- c(
            min(fr_exprs[, z]),
            max(fr_exprs[, z])
          )
        }
        # RESCALE
        rescale <- (fr_exprs[, z] - cyto_range[1]) /
          (cyto_range[2] - cyto_range[1])
        rescale[rescale > 1] <- 1
        rescale[rescale < 0] <- 0
        # POINT_COLOUR_SCALE
        col_scale <- grDevices::colorRamp(point_col_scale)
        # POINT COLOURS
        col <- col_scale(rescale)
        col <- grDevices::rgb(col[, 1],
                              col[, 2],
                              col[, 3],
                              maxColorValue = 255
        )
        return(col)
        # NAME OF A COLOUR
      } else {
        return(z)
      }
    } else {
      return(z)
    }
  })
  
  # Adjust colors by point_fill_alpha - REMOVE CHECK FOR ALPHA != 1
  for(z in seq_len(SMP)) {
    point_col[[z]] <- adjustcolor(point_col[[z]], point_col_alpha[z])
  }
  
  return(point_col)
}

## HISTOGRAMS ------------------------------------------------------------------

## DENSITY ----

#' Prepare histograms for cyto_plot
#' @param x list of cytoframes.
#' @noRd
.cyto_plot_hist <- function(x,
                            channels = NULL,
                            hist_stat = "count",
                            hist_smooth = 1,
                            hist_bins = 256,
                            hist_stack = 0,
                            ...) {
  
  # CHANNELS
  channels <- cyto_channels_extract(x[[1]], channels, plot = TRUE)
  
  # BANDWIDTH - BIN INSTRUMENT RANGE
  rng <- range(x[[1]][, channels], type = "instrument")
  hist_bandwidth <- LAPPLY(colnames(rng), function(z){
    (rng["max", z] - rng["min", z])/hist_bins
  })
  names(hist_bandwidth) <- channels
  # KERNEL DENSITY 
  d <- suppressWarnings(
    cyto_apply(x, 
               "cyto_stat_density",
               input = 2,
               channels = channels,
               stat = hist_stat,
               bandwidth = hist_bandwidth,
               smooth = hist_smooth,
               simplify = FALSE)
  )
  d <- lapply(d, `[[`, 1)
  # STACKING
  if(hist_stack != 0) {
    hist_heights <- round(LAPPLY(d, function(D){
      if(.all_na(D)){
        return(NA)
      }else{
        D$range[2]
      }
    }), 2)
    hist_stack <- max(hist_heights, na.rm = TRUE) * hist_stack
    hist_levels <- c(0, seq_along(x)) * hist_stack
    d <- lapply(seq_along(d), function(z){
      D <- d[[z]]
      if(z > 1 & !.all_na(D)){
        D$y <- D$y + hist_levels[z]
      }
      return(D)
    })
    # STORE RANGES IN NAMES
    names(d) <- LAPPLY(seq_along(d), function(z){
      if(is.na(hist_heights[z])) {
        paste(hist_levels[z], hist_levels[z], sep = "-")
      } else {
        paste(hist_levels[z], hist_levels[z] + hist_heights[z], sep = "-")
      }
    })
  } else {
    # STORE RANGES IN NAMES
    names(d) <- LAPPLY(d, function(z){
      if(.all_na(z)) {
        paste(c(0,0), collapse = "-")
      } else {
        paste(z$range, collapse = "-")
      }
    })
  }
  return(d)
}

## HISTOGRAM FILL ----

#' Get density fill colours for cyto_plot
#'
#' @param x list of flowFrame or density objects.
#' @param density_fill vector of colours to use for each layer.
#' @param density_cols vector of colls to use to select density_fill colours.
#'
#' @importFrom grDevices adjustcolor colorRampPalette
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_hist_fill <- function(x,
                                 hist_fill = NA,
                                 hist_cols = NA,
                                 hist_fill_alpha = 1) {
  
  # INHERIT CYTO_PLOT_THEME ----------------------------------------------------
  
  # Pull down arguments to named list
  args <- .args_list()
  
  # Inherit arguments from cyto_plot_theme
  args <- .cyto_plot_theme_inherit(args)
  
  # Update arguments
  .args_update(args)
  
  # GENERAL --------------------------------------------------------------------
  
  # Expected number of colours
  SMP <- length(x)
  
  # DENSITY_FILL ---------------------------------------------------------------
  
  # No hist_cols supplied
  if (.all_na(hist_cols)) {
    hist_cols <- .cyto_plot_colour_palette(type = "hist_cols")
  }
  
  # Make colorRampPalette
  if (class(hist_cols) != "function") {
    cols <- colorRampPalette(hist_cols)
  } else {
    cols <- hist_cols
  }
  
  # No colours supplied to hist_fill either
  if (.all_na(hist_fill)) {
    
    # Pull out a single colour per layer
    hist_fill <- cols(SMP)
    
    # Colours supplied manually to hist_fill
  } else {
    
    # Too few colours supplied - pull others from cols
    if (length(hist_fill) < SMP) {
      hist_fill <- c(
        hist_fill,
        cols(SMP - length(hist_fill))
      )
      
      # Too many colours supplied
    } else if (length(hist_fill) > SMP) {
      hist_fill <- hist_fill[seq_len(SMP)]
    }
  }
  
  # Adjust colors by hist_fill_alpha
  hist_fill <- mapply(function(hist_fill, hist_fill_alpha) {
    if (hist_fill_alpha != 1) {
      adjustcolor(hist_fill, hist_fill_alpha)
    } else {
      hist_fill
    }
  }, hist_fill, hist_fill_alpha, USE.NAMES = FALSE)
  
  return(hist_fill)
}

## THEME -----------------------------------------------------------------------

## INHERIT THEME ----

#' Inherit cyto_plot_theme arguments
#'
#' @param x list of named cyto_plot arguments.
#'
#' @return updated list of named arguments if cyto_plot_theme has been set.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_theme_inherit <- function(x) {
  
  # extract cyto_plot_theme arguments
  args <- getOption("cyto_plot_theme")
  
  if (!is.null(args)) {
    for(y in names(args)) {
      x[[y]] <- args[[y]]
    }
  }
  
  return(x)
}

## PALETTES --------------------------------------------------------------------

#' cyto_plot colour palette
#'
#' @param type indicates whether to return the "point_cols", "point_col_scale"
#'   or "hist_cols" colour palette.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_colour_palette <- function(type = "point_cols") {
  # POINT COLOUR PALETTE
  if (type == "point_cols") {
    pal <- c(
      "grey25",
      "bisque4",
      "brown1",
      "red",
      "darkred",
      "chocolate",
      "orange",
      "yellow",
      "yellowgreen",
      "green",
      "limegreen",
      "turquoise",
      "aquamarine",
      "cyan",
      "cornflowerblue",
      "blue",
      "blueviolet",
      "purple4",
      "purple",
      "magenta",
      "deeppink"
    )
    
    # POINT COLOUR SCALE
  } else if (type == "point_col_scale") {
    # NEW VIRIDIS PALETTE
    # viridis_pal <- viridis(9, option = "D")
    # plasma_pal <- viridis(12, option = "C")
    # viridis_pal <- viridis_pal[-9]
    # plasma_pal <- plasma_pal[-c(1:4)]
    # custom_pal <- c(viridis_pal, rev(plasma_pal))
    
    pal <- c("#440154",
             "#472D7B",
             "#3B528B",
             "#2C728E",
             "#21908C",
             "#27AD81",
             "#5DC863",
             "#AADC32",
             "#F0F921",
             "#FCD225",
             "#FDAD32",
             "#F58C46",
             "#E76F5A",
             "#D5546E",
             "#C03A83",
             "#A62098")
    
    # DENSITY COLOUR PALETTE
  } else if (type == "hist_cols") {
    pal <- c(
      "grey50",
      "bisque4",
      "brown1",
      "red",
      "darkred",
      "chocolate",
      "orange",
      "yellow",
      "yellowgreen",
      "green",
      "limegreen",
      "turquoise",
      "aquamarine",
      "cyan",
      "cornflowerblue",
      "blue",
      "blueviolet",
      "purple4",
      "purple",
      "magenta",
      "deeppink"
    )
  }
  
  return(pal)
}
