## CYTO_SPILLOVER_SPREAD_COMPUTE -----------------------------------------------

#' Compute total spillover spreading matrix
#' 
#' @export
cyto_spillover_spread_compute <- function(x,
                                          parent = "root",
                                          select = NULL,
                                          channels = NULL,
                                          save_as = NULL,
                                          compensated = FALSE,
                                          spillover = NULL,
                                          spillover_spread = NULL,
                                          axes_trans = NA,
                                          ...) {
  
  # BACKWARDS COMPATIBILITY ----------------------------------------------------
  
  # SPILLOVER_SPREAD
  if(!is.null(spillover_spread)) {
    warning(
      paste0(
        "'spillover_spread' is now deprecated in favour of 'save_as' - ",
        "please use that argument instead."
      )
    )
    save_as <- spillover_spread
  }
  
  # CHANNEL_MATCH DEPRECATED
  if("channel_match" %in% names(list(...))) {
    stop(
      paste0(
        "'channel_match' is no longer supported. CytoExploreR will ",
        "automatically create and search for a file called ",
        "'Compensation-Details.csv'."
      )
    )
  }
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # SELECT SAMPLES
  x <- cyto_select(
    x,
    select
  )
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    # EXCLUDE HEIGHT/WIDTH PARAMETERS
    channels <- channels[!grepl("-H$|-W$", channels, ignore.case = TRUE)]
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # MATCH CHANNELS
  pd <- cyto_channel_match(
    x,
    channels = channels
  )
  pd <- pd[match(rownames(cyto_details(x)), rownames(pd)), , drop = FALSE]
  
  # CHECK COMPENSATION DETAILS
  lapply(
    seq_len(nrow(pd)),
    function(z) {
      # CHANNEL
      if(!pd$channel[z] %in% c("Unstained",
                               "unstained",
                               cyto_channels(x))) {
        stop(
          paste0(
            pd$channel[z], 
            " is not a valid channel for this ", 
            cyto_class(x), 
            "!"
          )
        )
      }
      # PARENT
      if(cyto_class(x, "GatingSet")) {
        cyto_nodes_convert(x, pd$parent[z]) # ERRORS IF MISSING
      }
    }
  )
  
  # MISSING VARIABLES
  vars <- c("name", "group", "parent", "channel")
  if(!all(vars %in% colnames(pd))) {
    stop(
      paste0(
        "cyto_details(x) is missing required variables: ",
        paste0(
          vars[!vars %in% colnames(pd)],
          collapse = " & "
        )
      )
    )
  }
  
  
  
  
  
  
  
  
  
}
