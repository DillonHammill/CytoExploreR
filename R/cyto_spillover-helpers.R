## CYTO_SPILLOVER_POPS ---------------------------------------------------------

#' Return list of negative and positive populations using single stain controls
#' @param ... additional arguments passed to \code{cyto_plot}.
#' @importFrom flowWorkspace flowSet_to_cytoset
#' @return list of negative and positive flowSets
#' @noRd
.cyto_spillover_pops <- function(x,
                                 parent = NULL,
                                 axes_trans = NULL,
                                 channel_match = NULL,
                                 axes_limits = "machine",
                                 ...) {
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # PREPARE CHANNEL_MATCH VARIABLE (MARKERS TO CHANNELS)
  if (any(grepl("channel", colnames(pd), ignore.case = TRUE))) {
    # MARKERS TO CHANNELS
    ind <- which(grepl("channel", colnames(pd), ignore.case = TRUE))
    pd[, "channel"] <- LAPPLY(pd[, ind], function(z) {
      if (!grepl("unstained", z, ignore.case = TRUE)) {
        return(cyto_channels_extract(x, z))
      } else {
        return(z)
      }
    })
  }
  
  # CHANNELS
  channels <- cyto_fluor_channels(x)
  
  # TRANSFORMATIONS
  axes_trans <- cyto_transformers_extract(x)
  
  # PREPARE CHANNEL_MATCH ------------------------------------------------------
  
  # CHANNEL MATCH MISSING
  if (!"channel" %in% colnames(pd)) {
    # TRY CHANNEL_MATCH
    if (is.null(channel_match)) {
      pd$channel <- paste(cyto_channel_select(x))
    } else {
      if (is(channel_match, "data.frame") |
          is(channel_match, "matrix") |
          is(channel_match, "tibble")) {
        if (!all(c("name", "channel") %in% colnames(channel_match))) {
          stop("channel_match must contain columns 'name' and 'channel'.")
        }
        cm <- channel_match
        chans <- cm$channel[match(cyto_names(x), rownames(cm))]
        pd$channel <- paste(chans)
      } else {
        cm <- read_from_csv(channel_match)
        chans <- cm$channel[match(cyto_names(x), rownames(cm))]
        pd$channel <- paste(chans)
      }
    }
  }
  
  # PREPARE PARENTS ------------------------------------------------------------
  
  # PARENT PER CONTROL
  if (is(x, "GatingHierarchy") | is(x, "GatingSet")) {
    if (is.null(parent)) {
      if (!"parent" %in% colnames(pd)) {
        if (!is.null(channel_match)) {
          if ("parent" %in% colnames(cm)) {
            parent <- cm[, "parent"]
            pd[, "parent"] <- parent
          } else {
            nodes <- cyto_nodes(x, path = "auto")
            parent <- rep(nodes[length(nodes)], length(x))
            pd[, "parent"] <- parent
          }
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
  } else {
    pd[, "parent"] <- rep("root", length.out = length(x))
  }
  
  # UPDATE CHANNEL_MATCH FILE & CYTO_DETAILS -----------------------------------
  
  # SAVE CHANNEL_MATCH
  if (is.null(channel_match)) {
    channel_match <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Compensation-Channels.csv"
    )
  }
  
  # WRITE TO CSV
  write_to_csv(pd,
               channel_match)
  
  # CYTO_DETAILS - MAC ORDERING ISSUE
  lapply(colnames(pd), function(z) {
    cyto_details(x)[, z] <<- pd[, z]
  })
  
  # SAMPLE NAME COLUMN (AVOID INDEXING USING "NAME")
  pd_name <- colnames(pd)[which(LAPPLY(colnames(pd), function(z) {
    all(cyto_names(x) %in% pd[, z])
  }))]
  
  # REMOVE EXCESS CONTROLS -----------------------------------------------------
  
  # MULTIPLE CONTROLS PER CHANNEL (BYPASS UNSTAINED)
  if (length(unique(pd[, "channel"])) < length(pd[, pd_name])) {
    chans <- unique(pd[, "channel"])
    lapply(chans, function(z) {
      # RESTRICT CYTO_DETAILS
      pd_chunk <- pd[pd$channel == z, ]
      # MULTIPLE CONTROLS
      if (nrow(pd_chunk) > 1) {
        # BYPASS UNSTAINED
        if (!grepl("unstained", z, ignore.case = TRUE)) {
          # MESSAGE
          message(paste0(
            "Selecting the control with best signal for the ",
            z, " channel."
          ))
          # EXTRACT DATA
          fr_list <- lapply(seq_len(nrow(pd_chunk)), function(z) {
            cyto_extract(
              x[[pd_chunk$name[z]]],
              pd_chunk$parent[z],
              copy = TRUE
            )
          })
          names(fr_list) <- cyto_names(fr_list)
          fs <- as(fr_list, "flowSet")
          # CALCULATE MEDFI
          MEDFI <- suppressMessages(
            cyto_stats_compute(fs,
                               channels = z,
                               stat = "median",
                               trans = axes_trans
            )
          )
          # MAXIMUM SIGNAL
          max_MEDFI <- max(MEDFI[, ncol(MEDFI)])
          ind <- which(MEDFI[, ncol(MEDFI)] != max_MEDFI)
          remove_names <- MEDFI[ind, pd_name, drop = TRUE]
          # REMOVE MISSING FACTOR LEVELS
          if (is.factor(remove_names)) {
            droplevels(remove_names)
          }
          # REMOVE SAMPLES - LOW SIGNAL
          x <<- x[-match(remove_names, cyto_details(x)[, pd_name])]
          # REMOVE MISSING FACTOR LEVELS
          if (is.factor(cyto_details(x)[, pd_name])) {
            cyto_details(x)[, pd_name] <<- droplevels(
              cyto_details(x)[, pd_name]
            )
          }
        }
      }
    })
  }
  
  # RESTRICT CYTO_DETAILS
  pd <- cyto_details(x)
  
  # SPILLOVER POPULATIONS ------------------------------------------------------
  
  # UNIVERSAL REFERENCE
  if (any(grepl("unstained", pd[, "channel"], ignore.case = TRUE))) {
    
    # REMOVE EXCESS UNSTAINED CONTROLS - SELECT FIRST INSTANCE
    if (length(which(grepl("unstained", pd[, "channel"],
                           ignore.case = TRUE
    ))) > 1) {
      message("Removing excess unstained controls...")
      x <- x[-which(grepl("unstained", pd[, "channel"],
                          ignore.case = TRUE
      ))[-1]]
      pd <- cyto_details(x)
    }
    
    # EXTRACT UNSTAINED POPULATIONS - LIST OF FLOWFRAMES
    NIL_data <- x[which(grepl("unstained", pd[, "channel"],
                              ignore.case = TRUE
    ))]
    NIL <- lapply(unique(parent)[!is.na(unique(parent))], function(z) {
      # EITHER FLOWFRAME OR GATINGHIERARCHY
      cyto_extract(NIL_data[[1]],
                   z,
                   copy = TRUE
      )
    })
    names(NIL) <- unique(parent)
    
    # EXTRACT STAINED POPULATIONS - LIST OF FLOWFRAMES
    POS_gs <- x[-which(grepl("unstained", pd[, "channel"],
                             ignore.case = TRUE
    ))]
    POS <- lapply(seq_along(POS_gs), function(z) {
      cyto_extract(
        POS_gs[[z]],
        pd[, "parent"][match(cyto_names(POS_gs[[z]]), pd[, pd_name])],
        copy = TRUE
      )
    })
    names(POS) <- cyto_names(POS_gs)
    
    # RESTRICT CYTO_DETAILS
    pd <- pd[pd[, pd_name] != cyto_names(NIL_data), ]
    
    # NAMES
    nms <- names(POS)
    
    # SAMPLES
    smp <- length(POS)
    
    # GATE POSITIVE SIGNAL
    pos_pops <- list()
    neg_pops <- list()
    pops <- lapply(seq_len(smp), function(z) {
      
      # Extract flowFrame
      fr <- POS[[z]]
      
      # Parent
      parent <- pd[pd[, pd_name] == cyto_names(fr), "parent"]
      
      # Channel
      chan <- pd$channel[z]
      
      # Reference
      ref <- NIL[[parent]]
      
      # Plot
      cyto_plot(ref,
                channels = chan,
                overlay = fr,
                density_stack = 0,
                axes_trans = axes_trans,
                popup = TRUE,
                density_fill = c("red", "dodgerblue"),
                legend = FALSE,
                density_fill_alpha = 0.6,
                title = nms[z],
                axes_limits = axes_limits, ...
      )
      
      # cyto_gate_draw on each flowFrame using interval gate on selected channel
      gt <- cyto_gate_draw(
        x = fr,
        alias = paste0(chan, "+"),
        channels = chan,
        type = "interval",
        plot = FALSE
      )
      fr <- Subset(fr, gt[[1]])
      
      # SAVE GATED POPULATIONS
      pos_pops[[z]] <<- fr
      neg_pops[[z]] <<- ref
    })
    
    # INTERNAL REFERENCE - POSITIVE & NEGATIVE WITHIN CONTROL
  } else if (!any(grepl("unstained", pd[, "channel"], ignore.case = TRUE))) {
    
    # EXTRACT POPULATIONS
    pops <- lapply(seq_along(x), function(z) {
      cyto_extract(
        x[[z]],
        pd[, "parent"][match(cyto_names(x[[z]]), pd[, pd_name])],
        copy = TRUE
      )
    })
    names(pops) <- cyto_names(x)
    
    # NAMES
    nms <- names(pops)
    
    # SAMPLES
    smp <- length(pops)
    
    # NEGATIVE POPULATIONS
    neg_pops <- list()
    
    # POSITIVE POPULATIONS
    pos_pops <- list()
    
    # GATE NEGATIVE & POSITIVE SIGNAL PER CONTROL
    lapply(seq_len(smp), function(z) {
      
      # CONTROL
      fr <- pops[[z]]
      
      # CHANNEL
      chan <- pd$channel[z]
      
      # PLOT
      cyto_plot(fr,
                channels = chan,
                density_stack = 0,
                axes_trans = axes_trans,
                popup = TRUE,
                density_fill = "dodgerblue",
                legend = FALSE,
                density_fill_alpha = 0.6,
                title = nms[z],
                axes_limits = axes_limits, ...
      )
      
      # GATE NEGATIVE POPULATION
      gt <- cyto_gate_draw(
        x = fr,
        alias = paste0(chan, "-"),
        channels = chan,
        type = "interval",
        plot = FALSE
      )
      neg_pop <- Subset(fr, gt[[1]])
      
      # GaATE POSITIVE POPULATION
      gt <- cyto_gate_draw(
        x = fr,
        alias = paste0(chan, "+"),
        channels = chan,
        type = "interval",
        plot = FALSE
      )
      pos_pop <- Subset(fr, gt[[1]])
      
      neg_pops[[z]] <<- neg_pop
      pos_pops[[z]] <<- pos_pop
    })
  }
  
  # Add names to pops lists
  names(neg_pops) <- nms
  names(pos_pops) <- nms
  
  # Convert neg_pops and pos_pops to flowSets
  neg_pops <- flowSet_to_cytoset(flowSet(neg_pops))
  pos_pops <- flowSet_to_cytoset(flowSet(pos_pops))
  
  # cytoset required to retain details
  cyto_details(pos_pops) <- pd[pd$name %in% cyto_names(pos_pops), ]
  
  # Return list of negative and positive flowSets
  return(list(
    "negative" = neg_pops,
    "positive" = pos_pops
  ))
}
