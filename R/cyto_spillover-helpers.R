## CYTO_SPILLOVER_POPS ---------------------------------------------------------

#' Return list of negative and positive populations using single stain controls
#' @param ... additional arguments passed to \code{cyto_plot}.
#' @importFrom flowWorkspace flowSet_to_cytoset cytoset
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
      if (cyto_class(channel_match, c("matrix", "data.frame", "tibble"))) {
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
  if (cyto_class(x, "GatingSet")) {
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
  
  # INDEX BY ROWNAMES
  
  # REMOVE EXCESS CONTROLS -----------------------------------------------------
  
  # MULTIPLE CONTROLS PER CHANNEL (BYPASS UNSTAINED)
  if (length(unique(pd[, "channel"])) < length(rownames(pd))) {
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
          cf_list <- lapply(seq_len(nrow(pd_chunk)), function(z) {
            cyto_data_extract(
              x[rownames(pd_chunk)[z]], # subset by rownames not "name"
              pd_chunk$parent[z],
              copy = TRUE
            )[[1]][[1]]
          })
          names(cf_list) <- rownames(pd_chunk)
          cs <- cytoset(cf_list)
          # CALCULATE MEDFI
          MEDFI <- suppressMessages(
            cyto_stats_compute(cs,
                               channels = z,
                               stat = "median",
                               trans = axes_trans
            )
          )
          # MAXIMUM SIGNAL - CYTO_STATS_COMPUTE DROPS ROWNAMES
          max_MEDFI <- max(MEDFI[, ncol(MEDFI)])
          ind <- which(MEDFI[, ncol(MEDFI)] != max_MEDFI)
          remove_names <- rownames(cyto_details(cs))[ind]
          # # REMOVE MISSING FACTOR LEVELS
          # if (is.factor(remove_names)) {
          #   droplevels(remove_names)              # rownames can't be factor
          # }
          # REMOVE SAMPLES - LOW SIGNAL
          x <<- x[-match(remove_names, 
                         rownames(cyto_details(x)))]
          # # REMOVE MISSING FACTOR LEVELS
          # if (is.factor(rownames(cyto_details(x)))) {
          #   rownames(cyto_details(x)) <<- droplevels(
          #     rownames(cyto_details(x))
          #   )
          # }
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
      # CYTOSET
      cyto_data_extract(NIL_data[1],
                        z,
                        copy = TRUE)[[1]]
    })
    names(NIL) <- unique(parent)[!is.na(unique(parent))]
    
    # EXTRACT STAINED POPULATIONS - LIST OF FLOWFRAMES
    POS_gs <- x[-which(grepl("unstained", pd[, "channel"],
                             ignore.case = TRUE
    ))]
    POS <- lapply(seq_along(POS_gs), function(z) {
      # CYTOSET
      cyto_data_extract(
        POS_gs[z],
        pd[, "parent"][match(rownames(cyto_details(POS_gs))[z], rownames(pd))],
        copy = TRUE
      )[[1]]
    })
    names(POS) <- cyto_names(POS_gs)
    
    # RESTRICT CYTO_DETAILS
    pd <- pd[rownames(pd) != rownames(cyto_details(NIL_data)), ]
    
    # NAMES
    nms <- names(POS)
    
    # SAMPLES
    smp <- length(POS)
    
    # GATE POSITIVE SIGNAL
    pos_pops <- list()
    neg_pops <- list()
    pops <- lapply(seq_len(smp), function(z) {
      
      # CYTOSET
      cs <- POS[[z]]
      
      # Parent
      parent <- pd[rownames(pd) == rownames(cyto_details(cs)), "parent"]
      
      # Channel
      chan <- pd$channel[z]
      
      # Reference
      ref <- NIL[[parent]]
      
      # Plot
      cyto_plot(ref,
                channels = chan,
                overlay = cs,
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
        x = cs,
        alias = paste0(chan, "+"),
        channels = chan,
        type = "interval",
        plot = FALSE
      )
      cs <- Subset(cs, gt[[1]])
      
      # SAVE GATED POPULATIONS - CYTOFRAMES
      pos_pops[[z]] <<- cs[[1]]
      neg_pops[[z]] <<- ref[[1]]
    })
    
    
    
    # INTERNAL REFERENCE - POSITIVE & NEGATIVE WITHIN CONTROL
  } else if (!any(grepl("unstained", pd[, "channel"], ignore.case = TRUE))) {
    
    # EXTRACT POPULATIONS
    pops <- lapply(seq_along(x), function(z) {
      cyto_data_extract(
        x[z],
        pd[, "parent"][match(rownames(cyto_details(x[z])), rownames(pd))],
        copy = TRUE
      )
    })
    names(pops) <- cyto_names(x) # cytoset/GatingHierarchy/GatingSet
    
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
      cs <- pops[[z]]
      
      # CHANNEL
      chan <- pd$channel[z]
      
      # PLOT
      cyto_plot(cs,
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
        x = cs,
        alias = paste0(chan, "-"),
        channels = chan,
        type = "interval",
        plot = FALSE
      )
      neg_pop <- Subset(cs, gt[[1]])
      
      # GATE POSITIVE POPULATION
      gt <- cyto_gate_draw(
        x = cs,
        alias = paste0(chan, "+"),
        channels = chan,
        type = "interval",
        plot = FALSE
      )
      pos_pop <- Subset(cs, gt[[1]])
      
      # SAVE GATED POPULATIONS - CYTOFRAMES
      neg_pops[[z]] <<- neg_pop[[1]]
      pos_pops[[z]] <<- pos_pop[[1]]
    })
  }
  
  # Add names to pops lists
  names(neg_pops) <- nms
  names(pos_pops) <- nms
  
  # Convert neg_pops and pos_pops to cytosets
  neg_pops <- cytoset(neg_pops)
  pos_pops <- cytoset(pos_pops)
  
  # cytoset required to retain details
  cyto_details(pos_pops) <- pd[rownames(pd) %in% 
                                 rownames(cyto_details(pos_pops)), ]
  
  # Return list of negative and positive cytosets
  return(list(
    "negative" = neg_pops,
    "positive" = pos_pops
  ))
}
