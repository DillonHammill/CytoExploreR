## CYTO_PLOT_COMPENSATION ------------------------------------------------------

#' @importFrom flowWorkspace cytoset
#' @noRd
cyto_plot_comp <- function(x,
                           parent = NULL,
                           channels = NULL,
                           axes_trans = NA,
                           axes_limits = "machine",
                           title = NA,
                           hist_fill = c("blue", "grey"),
                           hist_fill_alpha = 0.5,
                           header = NULL,
                           header_text_size = 1,
                           header_text_font = 2,
                           header_text_col = "black",
                           ...) {

  # CYTO_PLOT_EXIT -------------------------------------------------------------

  # CURRENT SET PARAMETERS
  old_pars <- .par()

  # NEW PLOT METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    # CYTO_PLOT_METHOD
    options("cyto_plot_method" = "compensation")
    # CYTO_PLOT_EXIT
    on.exit({
      par(old_pars)
      options("cyto_plot_method" = NULL)
      options("cyto_plot_par" = NULL)
    })
  }

  # CHECKS ---------------------------------------------------------------------

  # EXPERIMENT DETAILS
  pd <- cyto_details(x)

  # CHANNELS
  if (is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  }

  # TRANSFORMATIONS
  axes_trans <- cyto_transformers_extract(x)

  # CHANNEL_MATCH --------------------------------------------------------------

  # CYTO_DETAILS - NO CHANNEL
  if (!"channel" %in% colnames(pd)) {
    # CHANNEL_MATCH SUPPLIED
    if (!is.null(channel_match)) {
      # CHANNEL NAMES OR FILE
      if (is.character(channel_match)) {
        # CHANNEL/MARKER
        if (all(channel_match %in% c(
          cyto_channels(x),
          cyto_markers(x)
        ))) {
          channel_match <- cbind(pd, "channel" = channel_match)
          # FILE NAME
        } else {
          channel_match <- read_from_csv(channel_match)
        }
        # DATA.FRAME OR MATRIX
      } else {
        if (cyto_class(channel_match, c(
          "data.frame",
          "matrix",
          "tibble"
        ))) {
          if (!all(c("name", "channel") %in% colnames(channel_match))) {
            stop("'channel_match' must contain columns 'name' and 'channel'.")
          }
        }
      }
      # SORT CHANNEL_MATCH
      channel_match <- channel_match[match_ind(
        cyto_names(x),
        channel_match$name
      ), ]
      # UPDATE CYTO_DETAILS
      pd$channel <- channel_match$channel
      # NO CHANNEL_MATCH
    } else if (is.null(channel_match)) {
      # INTERACTIVELY SELECT CHANNEL
      channel_match <- data.frame(
        "name" = cyto_names(x),
        "channel" = cyto_channel_select(x)
      )
      pd$channel <- channel_match$channel
    }
  } else {
    channel_match <- pd
  }

  # PARENTS --------------------------------------------------------------------

  # PARENT MISSING
  if (is.null(parent)) {
    if (!"parent" %in% colnames(pd)) {
      if ("parent" %in% colnames(channel_match)) {
        parent <- channel_match[, "parent"]
        pd[, "parent"] <- parent
      } else {
        nodes <- cyto_nodes(x, path = "auto")
        parent <- rep(nodes[length(nodes)], length(x))
        pd[, "parent"] <- parent
      }
    }
  } else {
    parent <- rep(parent, length.out = length(x))
    pd[, "parent"] <- parent
  }

  # PREPARE DATA ---------------------------------------------------------------

  # UPDATE CYTO_DETAILS
  cyto_details(x) <- pd

  # ISOLATE UNSTAINED CONTROL
  if (any(grepl("unstained", pd$channel, ignore.case = TRUE))) {
    ind <- which(grepl("unstained", pd$channel, ignore.case = TRUE))[1]
    NIL <- x[ind]
    x <- x[-ind]
    pd <- cyto_details(x)
    neg_pops <- lapply(seq_along(x), function(z) {
      cyto_extract(
        NIL,
        pd$parent[z],
        copy = TRUE
      )
    })
    names(neg_pops) <- cyto_names(x)
  } else {
    neg_pops <- structure(rep(list(NA, length.out = length(pos_pops))), 
                          names = names(pos_pops))
  }

  # POSITIVE POPULATIONS
  pos_pops <- lapply(seq_along(x), function(z) {
    cyto_extract(x[[z]],
      parent = pd$parent[z],
      copy = TRUE
    )
  })
  names(pos_pops) <- cyto_names(x)

  # COMPENSATE & TRANSFORM -----------------------------------------------------

  # COMPENSATE & TRANSFORM
  if (compensate) {
    # LINEARISE & COMPENSATE
    pos_pops <- cyto_apply(pos_pops,
      "cyto_compensate",
      spillover = spillover,
      trans = axes_trans,
      inverse = TRUE,
      simplify = FALSE,
      copy = FALSE
    )
    if (!is.na(neg_pops[[1]])) {
      neg_pops <- cyto_apply(neg_pops,
        "cyto_compensate",
        spillover = spillover,
        trans = axes_trans,
        inverse = TRUE,
        simplify = FALSE,
        copy = FALSE
      )
    }
    # TRANSFORMERS
    if (.all_na(axes_trans)) {
      axes_trans <- cyto_transformer_biex(cytoset(pos_pops),
        plot = FALSE
      )
    }
    # TRANSFORM
    pos_pops <- structure(
      lapply(pos_pops, function(cf) {
        cyto_transform(cf,
          trans = axes_trans,
          inverse = FALSE,
          plot = FALSE
        )
      }),
      names = names(pos_pops)
    )
    if (!is.na(neg_pops[[1]])) {
      neg_pops <- structure(
        lapply(neg_pops, function(cf) {
          cyto_transform(cf,
            trans = axes_trans,
            inverse = FALSE,
            plot = FALSE
          )
        }),
        names = names(neg_pops)
      )
    }
  # TRANSFORM ONLY
  } else {
    # TRANSFORMERS
    if (.all_na(axes_trans)) {
      axes_trans <- cyto_transformer_biex(cytoset(pos_pops),
        plot = FALSE
      )
      # TRANSFORM
      pos_pops <- structure(
        lapply(pos_pops, function(cf) {
          cyto_transform(cf,
            trans = axes_trans,
            inverse = FALSE,
            plot = FALSE
          )
        }),
        names = names(pos_pops)
      )
      if (!is.na(neg_pops[[1]])) {
        neg_pops <- structure(
          lapply(neg_pops, function(cf) {
            cyto_transform(cf,
              trans = axes_trans,
              inverse = FALSE,
              plot = FALSE
            )
          }),
          names = names(neg_pops)
        )
      }
    }
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # TITLE
  if(is.null(header)) {
    header <- names(pos_pops)
  }
  
  # LAYOUT
  layout <- .cyto_plot_layout(channels, layout)
  
  # OUTER MARGINS
  if(!.all_na(header)) {
    outer_margins <- c(0, 0, 3, 0)
  }
  
  # PLOT CONSTRUCTION ----------------------------------------------------------
  
  # OPEN NEW DEVICE
  cyto_plot_new(popup,
                layout,
                outer_margins)
  
  # PLOTS
  plots <- lapply(seq_along(pos_pops), function(z){
    # POPULATIONS
    pos_pop <- pos_pops[[z]]
    neg_pop <- neg_pops[[z]]
    # CHANNEL
    channel <- pd$channel[pd$name == names(pos_pops)[z]]
    # CHANNELS
    p <- lapply(seq_along(channels), function(y) {
      # MATCHING CHANNEL - HISTOGRAM
      if(y == channel) {
        cyto_plot(pos_pop,
                  overlay = neg_pop,
                  channels = y,
                  axes_trans = axes_trans,
                  axes_limits = axes_limits,
                  title = title,
                  legend = FALSE,
                  hist_stack = hist_stack,
                  hist_fill = hist_fill,
                  hist_fill_alpha = hist_fill_alpha,
                  ...)
      # NON-MATCHING CHANNEL - SCATTERPLOT  
      } else {
        cyto_plot(pos_pop,
                  overlay = neg_pop,
                  channels = c(channels, y),
                  axes_trans = axes_trans,
                  axes_limits = axes_limits,
                  legend = FALSE,
                  title = title,
                  ...)
      }
      # HEADER & NEW DEVICE
      if(par("page") | y == length(channels)) {
        # HEADER
        if(!.all_na(header)) {
          mtext(header[z],
                outer = TRUE,
                cex = header_text_size,
                font = header_text_font,
                col = header_text_col)
        }
        # NEW DEVICE
        if(z != length(pos_pops)) {
          cyto_plot_new()
        }
      }
    })
    
  })
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE FOR SAVING
  if (getOption("cyto_plot_save")) {
    if (getOption("cyto_plot_method") == "compensation") {
      # CLOSE GRAPHICS DEVICE
      dev.off()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
    }
  }
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}
