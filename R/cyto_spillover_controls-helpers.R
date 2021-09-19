## CYTO_SPILL_ANNOTATE ---------------------------------------------------------

#' Annotate compensation controls with channels and parents
#' @noRd
.cyto_spillover_controls_annotate <- function(x,
                                              parent = "root",
                                              channel_match = NULL) {
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # CHANNEL_MATCH NOT SUPPLIED - MAY BE PRESENT
  if(is.null(channel_match)) {
    # SEARCH FOR CHANNEL MATCH FILE
    channel_match <- tryCatch(
      list.files()[grepl(
        "channel-match.csv$",
        list.files(),
        ignore.case = TRUE
      )],
      error = function(e) {
        return(NULL)
      }
    )
    if(length(channel_match) == 0){
      channel_match <- NULL
    }
  }
  
  # CHANNEL_MATCH SUPPLIED
  if(!is.null(channel_match)) {
    # CHANNEL_MATCH - VECTOR OF CHANNELS OR FILE NAME
    if(is.null(dim(channel_match))) {
      # CHANNEL_MATCH FILE
      if(all(file_exists(file_ext_append(channel_match)))) {
        # FIND CORRECT CHANNEL_MATCH FILE
        channel_match <- lapply(channel_match, function(z){
          # READ FILE
          cm <- read_from_csv(z)
          # CHECK COLUMNS - CYTO_CHANNEL_MATCH()
          if(!all(c("name", "channel") %in% colnames(cm))) {
            return(NULL)
            # VALID CHANNEL MATCH
          } else {
            # CHANNEL_MATCH - MUST CONTAIN ALL SAMPLES
            if(length(cyto_match(x, cm[, "name"]) == nrow(cm))) {
              return(cm)
            } else {
              return(NULL)
            }
          }
        })
        channel_match <- channel_match[!LAPPLY(channel_match, "is.null")][[1]]
        # CHANNELS/MARKERS
      } else {
        channel_match <- pd
        channel_match$channel <- cyto_channels_extract(x, 
                                                       channels = channel_match, 
                                                       skip = "Unstained")
      }
      # CHANNEL_MATCH - ARRAY
    } else {
      # CHECK COLNAMES
      if(!all(c("name", "channel") %in% colnames(channel_match))) {
        warning(
          "Use cyto_channel_match() to create 'channel_match' file!"
        )
        channel_match <- NULL
      }
      # CHECK SAMPLES - CREATE FROM SCRATCH IF ARE SAMPLES MISSING
      if(length(cyto_match(x, channel_match[, "name"]) != nrow(channel_match))) {
        warning(
          paste0(
            "'channel_match' does not cover all samples in this ",
            cyto_class(x), "!"
          )
        )
        channel_match <- NULL
      }
    }
  }
  
  # ANNOTATE CHANNELS
  if(!"channel" %in% colnames(pd)) {
    # CREATE NEW CHANNEL_MATCH FILE
    if(is.null(channel_match)) {
      channel_match <- cyto_channel_match(x)
    }
    # UPDATE CHANNELS IN EXPERIMENT DETAILS
    pd <- cbind(pd, 
                "channel" = channel_match[
                  cyto_match(x, channel_match[, "name"]), "channel"])
    cyto_details(x) <- pd
  }
  
  # ANNOTATE PARENT
  if(!any(grepl("parent", colnames(pd), ignore.case = TRUE))) {
    pd$parent <- rep(parent, length.out = length(x))
    cyto_details(x) <- pd
  }
  
  return(x)
  
}

## CYTO_SPILLOVER_CONTROLS -----------------------------------------------------

#' Prepare compensation controls for cyto_spillover_compute
#' @importFrom flowWorkspace cytoset
#' @noRd
.cyto_spillover_controls_prepare <- function(x,
                                             parent = "root",
                                             type = "roca",
                                             channel_match = NULL,
                                             axes_trans = NA) {
  
  # ANNOTATE CONTROLS
  x <- .cyto_spillover_controls_annotate(x,
                                         parent = parent,
                                         channel_match = channel_match)
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  #ISOLATE UNSTAINED CONTROL
  if(any(grepl("Unstained", pd[, "channel"], ignore.case = TRUE))) {
    ind <- cyto_match(x, channel = "Unstained")
    x_NIL <- x[ind]
    x <- x[-ind]
  } else {
    x_NIL <- NULL
  }
  
  # EXTRACT CYTOSET COPY
  cs <- cytoset(
    structure(
      lapply(seq_along(x), function(z){
        cyto_data_extract(
          x[z],
          parent = pd$parent[z],
          format = "cytoframe",
          copy = TRUE,
          trans = axes_trans,
          inverse = if(grepl("^r", type, ignore.case = TRUE)) {
            TRUE         # LINEAR DATA REQUIRED
          } else {
            FALSE        # TRANSFORMED DATA REQUIRED
          }
        )[[1]][[1]]
      }),
      names = cyto_names(x)
    )
  )
  
  # UPDATE EXPERIMENT DETAILS
  cyto_details(cs) <- pd[cyto_match(cs, rownames(pd)), ]
  
  # UPDATED EXPERIMENT DETAILS
  pd <- cyto_details(cs)

  # IDENTIFY HIGHEST QUALITY CONTROLS
  if(any(duplicated(pd$channel))) {
    # SAMPLES TO KEEP
    cyto_ind <- c()
    lapply(unique(pd$channel), function(z){
      # MULTIPLE UNSTAINED CONTROLS - USE MIN AVERAGE CV
      if(grepl("Unstained", z, ignore.case = TRUE)) {
        # COMPUTE CV IN ALL CHANNELS
        cv <- cyto_stats_compute(
          cs[pd$channel == z],
          channels = pd[, "channel"][pd[, "channel"] != z],
          stat = "cv",
          inverse = TRUE,
          trans = axes_trans, 
          details = FALSE,
          markers = FALSE,
          format = "wide"
        )[, z]
        cyto_ind <<- c(
          cyto_ind, 
          which.min(
            rowMeans(cv[, pd[, "channel"][pd[, "channel"] != z]])
          )
        )
        # MULTIPLE STAINED CONTROLS
      } else {
        # COMPUTE CHANNEL MEDFI
        cyto_ind <<- c(
          cyto_ind,
          which.max(
            cyto_stats_compute(
              cs[pd$channel == z],
              channels = z,
              stat = "median",
              inverse = TRUE,
              trans = axes_trans, 
              details = FALSE,
              markers = FALSE,
              format = "wide"
            )[, z]
          ))
      }
    })
    cs <- cs[cyto_ind]
  }

  # UPDATE EXPERIMENT DETAILS
  pd <- cyto_details(cs)
  
  # ROCA METHOD
  if(grepl("^r", type, ignore.case = TRUE)) {
    return(cs)
  # BAGWELL METHOD
  } else {
    # EXTRACT PARENTS FOR UNSTAINED CONTROL
    if(!is.null(x_NIL)) {
      x_NIL <- cytoset(
        structure(
          lapply(pd$channel, function(z){
            cyto_data_extract(
              x_NIL,
              parent = pd[, "parent"][pd[, "channel"] == z],
              format = "cytoset",
              copy = TRUE,
              trans = NA,
              inverse = FALSE
            )[[1]][[1]]
          }),
          names = cyto_names(cs)
        )
      )
    }
    # LIST OF NEGATIVE & POSITIVE EVENTS PER CHANNEL
    return(
      structure(
        lapply(cyto_names(cs), function(z){
          list(
            "-" = if(!is.null(x_NIL)){
                cyto_select(x_NIL, z)
              } else {
                NULL
              },
            "+" = cyto_select(cs, z)
          )
        }),
        names = pd[, "channel"]
      )
    )
  }
  
}

## CYTO_SPILLOVER_POPS ---------------------------------------------------------

#' Return list of negative and positive populations using single stain controls
#' @param ... additional arguments passed to \code{cyto_plot}.
#' @importFrom flowWorkspace flowSet_to_cytoset cytoset
#' @return list of negative and positive flowSets
#' @noRd
.cyto_spillover_controls_pops <- function(x,
                                          axes_trans = NA,
                                          axes_limits = "machine",
                                          ...) {
  
  # GATE POPULATIONS PER CHANNEL
  structure(
    lapply(names(x), function(z) {
      # REFERENCE CYTOSETS
      neg_events <- x[[z]][["-"]]
      pos_events <- x[[z]][["+"]]
      # CYTO_PLOT
      cyto_plot(
        if(is.null(neg_events)){
          pos_events
        } else {
          neg_events
        },
        channels = z,
        overlay = if(!is.null(neg_events)){
          pos_events
        } else {
          NA
        },
        hist_stack = 0,
        hist_fill = if(is.null(neg_events)){
          "dodgerblue"
        } else {
          c("red", "dodgerblue")
        },
        hist_fill_alpha = 0.6,
        title = cyto_names(pos_events),
        axes_limits = axes_limits, 
        axes_trans = axes_trans,
        legend = FALSE,
        ...
      )
      # GATE NEGATIVE POPULATION
      if(is.null(neg_events)) {
        neg_gt <- cyto_gate_draw(
          x = pos_events,
          alias = paste0(z, "-"),
          channels = z,
          type = "interval",
          plot = FALSE
        )
      } else {
        neg_gt <- NULL
      }
      # GATE POSITIVE POPULATION
      pos_gt <- cyto_gate_draw(
        x = pos_events,
        alias = paste0(z, "+"),
        channels = z,
        type = "interval",
        plot = FALSE
      )
      # RETURN GATED POSITIVE/NEGATIVE POPULATIONS
      return(
        list(
          "-" = if(!is.null(neg_gt)) {
            cyto_gate_apply(
              pos_events,
              neg_gt
            )[[1]][[1]] # LIST OF POPULATIONS PER GATE
          } else {
            neg_events
          },
          "+" = cyto_gate_apply(
            pos_events,
            pos_gt
          )[[1]][[1]] # LIST OF POPULATIONS PER GATE
        )
      )
    }),
    names = names(x)
  )
  
}
