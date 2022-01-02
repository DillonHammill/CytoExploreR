## CYTO_CHANNEL_MATCH ----------------------------------------------------------

#' Matches each compensation control to a fluorescent channel
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels names of the possible channels or markers to be matched to
#'   each control, set to all area fluorescent parameters by default.
#' @param save_as name of a CSV file to which the channel matching should be
#'   written for downstream use, set to \code{"Compensation-Details.csv"}
#'   prefixed with the date by default. Users can set custom file names here,
#'   but the file name should contain \code{"Compensation-Details.csv"} in order
#'   to be automatically detected by CytoExploreR within
#'   \code{cyto_spillover_compute()}, \code{cyto_spillover_edit()},
#'   \code{cyto_spillover_spread()} and \code{cyto_plot_compensation()}.
#' @param strip logical indicating whether overlapping characters in file names
#'   should be stripped prior to matching, set to TRUE by default for more
#'   accurate matching.
#' @param ignore.case logical indicating whether all case insensitive matches
#'   should be found, set to TRUE by default.
#' @param ... additional arguments passed to \code{link{grep}} when performing
#'   character matching.
#'
#' @return a data.frame written to a CSV file containing information matching
#'   each file name to a channel. This channel matching information is also
#'   automatically added to the \code{cyto_details()} of the supplied samples
#'   where it can be easily accessed by CytoExploreR downstream.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Compensation GatingSet
#' gs <- GatingSet(Compensation)
#'
#' # Channel matching
#' cyto_channel_match(gs)
#'
#' @export
cyto_channel_match2 <- function(x,
                                channels = NULL,
                                save_as = NULL,
                                strip = TRUE,
                                ignore.case = TRUE,
                                ...) {
  
  # TODO: GROUP & PARENT COLUMNS TOO?
  
  # CYTOFRAMES NOT SUPPORTED
  if(cyto_class(x, "flowFrame")) {
    stop(
      paste0(
        "cyto_channel_match() only supports objects of class cytoset, ",
        "GatingHierarchy or GatingSet!"
      )
    )
  }
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # CREATE CHANNEL MATCH OUTPUT
  channel_match <- cbind(
    pd[, "name", drop = FALSE],
    "channel" = NA
  )
  
  # NAMES
  file_names <- rownames(pd)
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    channels <- channels[!grepl("-H$|-W$", channels, ignore.case = TRUE)]
  } else {
    channels <- cyto_channels_extract(x, channels = channels)
  }
  
  # LOCATE UNSTAINED CONTROL(S)
  unst_ind <- grep("Unst|NIL", file_names, ignore.case = TRUE)
  
  # ANNOTATE UNSTAINED CONTROLS - REMOVE FROM FILE NAMES
  if(length(unst_ind) > 0) {
    channel_match$channel[unst_ind] <- "Unstained"
    file_names <- file_names[-unst_ind]
  }
  
  # MATCH FILENAMES TO MARKERS/CHANNELS ----------------------------------------
  
  # ADDITIONAL FILE NAMES TO MATCH
  if(length(file_names) > 0) {
    # REPLACE ABNORMAL CHARACTERS
    file_names_strip <- gsub(",2f,", "/", file_names, ignore.case = TRUE)
    # REMOVE WHITESPACE & SPECIAL CHARACTERS
    file_names_strip <- gsub("[^[:alnum:]]", "", file_names_strip)
    # STRIP OVERLAPPING CHARACTERS - NO PADDING REQUIRED
    if(strip) {
      file_names_strip <- .cyto_string_strip(file_names_strip)
    }
    # COMBINE MARKERS & CHANNELS FOR MATCHING - "MARKER CHANNEL"
    channels_strip <- unlist(
      lapply(
        seq_along(channels),
        function(z) {
          # CHANNEL
          channel <- channels[z]
          # MARKER(S) ASSIGNED
          if(!is.na(names(channel)) & names(channel) != "NA") {
            # HANDLE MULTIPLE MARKERS PER CHANNEL
            markers <- strsplit(names(channel), "\\||\\/")[[1]]
            # STORE CHANNEL NAME
            channel <- rep(channel, length.out = length(markers))
            # APPEND CHANNEL NAME
            channel <- structure(
              paste(
                markers,
                channel
              ),
              names = channel
            )
            # NO MARKER
          } else {
            names(channel) <- channel # STORE CHANNEL NAME 
          }
          return(channel)
        }
      )
    )
    # MATCH FILENAMES TO MARKER/CHANNEL COMBOS
    file_channel_match <- structure(
      lapply(
        seq_along(file_names_strip),
        function(z) {
          # SPLIT FILE NAME INTO INDIVIDUAL CHARACTERS
          file_name_split <- strsplit(file_names_strip[[z]], "")[[1]]
          # LOOP THROUGH MARKER/CHANNEL COMBOS
          channel_match_opts <- structure(
            lapply(
              channels_strip,
              function(channel) {
                # SPLIT ALPHNUMERIC CHARACTERS - NO SUFFIX
                channel_split <- strsplit(
                  gsub(
                    "[^[:alnum:]]",
                    "",
                    gsub(
                      "-A$|-H$|-W$",
                      "",
                      channel
                    ),
                  ),
                  ""
                )[[1]]
                # MATCHES PER CHARACTER
                channel_split_match <- structure(
                  lapply(
                    channel_split,
                    function(char) {
                      grep(
                        char, 
                        file_name_split, 
                        ignore.case = ignore.case, 
                        ...
                      )
                    }
                  ),
                  names = channel_split
                )
                # NUMBER OF CHARACTERS
                n <- length(channel_split_match)
                # STORE MATCHES
                cm <- list()
                # STARTING LETTER INDEX
                for(i in seq_len(n)) {
                  # TERMINATE IF LARGER MATCH CANNOT BE LOCATED
                  if(length(cm) > 0) {
                    if((n - i) < max(unlist(lapply(cm, "sum")))) {
                      break()
                    }
                  }
                  # STORE MATCHING CHARACTERS
                  m <- c()
                  # LETTER MATCH 
                  for(v in channel_split_match[[i]]) {
                    # ADD LETTER AS MATCH
                    m <- structure(
                      c(1),
                      names = names(channel_split_match[i])
                    )
                    # DESCEND TREE - SEARCH FOR INDEX +1 - ALLOW DELETIONS ONLY
                    if(i + 1 <= n) {
                      for(q in (i + 1):n) {
                        # EMPTY MATCH
                        if(length(channel_split_match[[q]]) == 0) {
                          # ALLOW EMPTY MATCH
                          m <- c(
                            m, 
                            structure(
                              0,
                              names = names(channel_split_match[q])
                            )
                          )
                        } else {
                          # TODO: MAKE MATCHING MORE SENSITIVE
                          # MATCH LOCATED
                          if(any(channel_split_match[[q]] > v)) {
                            v <- channel_split_match[[q]][
                              channel_split_match[[q]] > v
                            ][1]
                            m <- c(
                              m,
                              structure(
                                1,
                                names = names(channel_split_match[q])
                              )
                            )
                          # NO MATCH FOUND
                          } else {
                            break()
                          }
                        }
                      }
                    } else {
                      break()
                    }
                  }
                  # STORE VALID MATCHESD - MULTIPLE CHARACTERS
                  if(sum(m) > 1) { 
                    cm <- c(cm, list(m))
                  }
                }
                # RETURN LONGEST MATCH
                if(length(cm) == 0) {
                  return(0)
                } else {
                  return(
                    max(
                      unlist(
                        lapply(
                          cm,
                          "sum"
                        )
                      )
                    )
                  )
                }
              }
            ),
            names = channels
          )
          # UPDATE CHANNEL_MATCH
          if(sum(unlist(channel_match_opts)) > 0) {
            ind <- which(
              channel_match_opts == max(unlist(channel_match_opts))
            )
            # MULTIPLE CHANNEL MATCHES - CHOOSE SHORTEST CHANNEL OPTION
            if(length(ind) > 0) {
              ind <- ind[which.min(nchar(channels[ind]))]
              # CANNOT ASSIGN CHANNELS WITH LENGTH - AMBIGUOUS
              if(length(ind) == 1) {
                channel_match[
                  match(file_names[z], rownames(channel_match)),
                  "channel"
                ] <<- channels[ind] # NEED TO STORE ORIGINAL CHANNELS
              }
            }
          }
        }
      ),
      names = file_names
    )
  }
  
  # INTERACTIVE EDITING & EXPORT -----------------------------------------------
  
  # UNASSIGNED CHANNELS - RESORT TO CHANNEL STATISTICS
  if(any(is.na(channel_match$channel))) {
    warning(
      paste0(
        "The following samples could not be matched to a channel by name:",
        paste0(
          channel_match$name[is.na(channel_match$channel)],
          collapse = "\n"
        ),
        "\n",
        "CytoExploreR will make an educated guess for these samples using ",
        "the intensities in the fluorescent channels."
      )
    )
    # TODO: USE INTENSITIES TO CHOOSE CHANNELS
  }
  
  # INTERACTIVE CHANNEL MATCHING
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    rn <- rownames(channel_match)
    rownames(channel_match) <- NULL
    channel_match <- data_edit(
      channel_match,
      logo = CytoExploreR_logo(),
      title = "Channel Match Editor",
      row_edit = FALSE,
      col_readonly = "name",
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane"
    )
    rownames(channel_match) <- rn
  }

  # UPDATE DETAILS IN SAMPLES
  cyto_details(x)$channel <- channel_match$channel

  # SAVE_AS
  if(!.all_na(save_as)) {
    # DEFAULT FILE NAME
    if(is.null(save_as)) {
      save_as <- paste0(
        format(Sys.Date(), "%d%m%y"),
        "-", "Compensation-Details.csv"
      )
    }
    # WRITE TO CSV
    write_to_csv(
      channel_match,
      save_as,
      row.names = TRUE
    )
  }
  
  # RETURN CHANNEL MATCHING
  return(channel_match)
  
}

#' Helper function to remove overlapping string fragments
#' @noRd
.cyto_string_strip <- function(x,
                               pad = FALSE) {
  # MULTIPLE STRINGS REQUIRED
  if(length(x) == 1) {
    return(x)
  }
  # MAX SEARCH DEPTH
  depth <- max(nchar(x))
  # LEFT SEARCH
  for(i in 1:depth) {
    if(length(unique(substring(x, 1, 1))) == 1) {
      x <- gsub("^.", "", x)
    } else {
      break()
    }
  }
  # RIGHT SEARCH
  for(i in 1:depth) {
    frag <- unlist(
      lapply(
        x,
        function(z) {
          substring(
            z,
            nchar(z),
            nchar(z)
          )
        }
      )
    )
    if(length(unique(frag)) == 1) {
      x <- gsub(".$", "", x)
    } else {
      break()
    }
  }
  # PADDING - SAME WIDTH
  if(pad) {
    x <- format(x, width = max(nchar(x)))
  }
  # REMAINDER
  return(x)
}
