## CYTO_PLOT_COMPENSATION ------------------------------------------------------

#' Visualise compensation of spillover in all channels
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing data for
#'   single colour compensation controls.
#' @param parent names of the parent populations to extract from each control
#'   for plotting.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details.
#' @param overlay can be either \code{"unstained"}, \code{"compensated"},
#'   \code{"both"} or \code{"none"} to allow control over which data is overlaid
#'   onto the plots, set to \code{"both"} by default.
#' @param channels names of the channels or markers in which comensation should
#'   be visualised, set to all area fluorescentce parameter by default.
#' @param spillover a matrix or name of a CSV file containing the spillover
#'   coefficients that should be applied to the data when \code{compensate} is
#'   TRUE.
#' @param compensated logical required for \code{cytoset} objects to indicate
#'   whether the data has already been compensated prior to plotting, set to
#'   FALSE by default.
#' @param channel_match for internal use only.
#' @param axes_trans object of class \code{transformerList} containing the
#'   definitions of the transformations that have been applied to the data.
#'   Transformers will be automatically extracted from \code{GatingSet} objects
#'   but must be supplied manually for \code{cytoset} objects. If no transformer
#'   definitions are extracted or supplied, new transformers of the logicle type
#'   will be generated internally and applied to the data for better
#'   visualisation.
#' @param axes_limits options include \code{"auto"}, \code{"trim"},
#'   \code{"data"} or \code{"machine"} to use optimised, lower trimmed (remove
#'   lowest 1% of events), data or machine limits respectively. Set to
#'   \code{"auto"} by default to use optimised axes ranges. Fine control over
#'   axes limits can be obtained by altering the \code{xlim} and \code{ylim}
#'   arguments.
#' @param layout a vector of the length 2 of form \code{c(#rows, #columns)} or a
#'   matrix indicating the dimensions of the grid for plotting.
#' @param header a vector of text to include at the top of each plot layout, set
#'   to the names of each single colour control by default.
#' @param point_col vector of colours to use for layers in scatter plots, set to
#'   \code{c("magenta", "blue", "grey40")} by default.
#' @param point_col_alpha numeric [0,1] to control point colour transparency in
#'   2-D scatter plots, set to \code{0.5} by default.
#' @param hist_fill vector of colours to use for layers in 1-D histograms, by
#'   default matches \code{point_col} for consistency.
#' @param hist_fill_alpha numeric [0,1] used to control histogram fill colour
#'   transparency, set to \code{0.5} by default.
#' @param lines logical indicating whether a robust linear model should be
#'   fitted to the data and added to the plots, set to TRUE by default.
#' @param line_type type of line to use for fitted robust linear models, set to
#'   1 by default for solid lines.
#' @param line_width width of lines to use for fitted robust linear models, set
#'   to 2 by default.
#' @param line_col colour to use for fitted robust linear models, set to
#'   \code{"red"} by default.
#' @param line_col_alpha numeric [0,1] to control transparency of fitted robust
#'   linear models, set to 1 by default to use solid colours.
#' @param text logical indicating whether the plots should be annotated with
#'   their respective spillover coefficients in \code{spillover}, set to TRUE by
#'   default.
#' @param text_font font to use for spillover coefficient text, set to 1 by
#'   default for plain font.
#' @param text_size numeric to control the size of the spillover coefficient
#'   text, set to 1.2 by default.
#' @param text_col colour to use for spillover coefficient text, set to
#'   \code{"black"} by default.
#' @param text_col_alpha numeric [0,1] to control transparency of spillover
#'   coefficient text, set to 1 by default to use solid colours.
#' @param header_text_font numeric to control the font of the header text, set
#'   to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param header_text_size numeric to control the size of the header text, set
#'   to 1 by default.
#' @param header_text_col colour to use for the header text, set to "black" by
#'   default.
#' @param ... additional arguments passed to \code{\link{cyto_plot}}.
#'
#' @importFrom flowWorkspace cytoset
#' @importFrom graphics abline text
#' @importFrom MASS rlm
#' @importFrom grDevices adjustcolor
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @return a list of recorded plots for each single colour control.
#'
#' @export
cyto_plot_compensation <- function(x,
                                   parent = "root",
                                   select = NULL,
                                   overlay = "both",
                                   channels = NULL,
                                   spillover = NULL,
                                   compensated = FALSE,
                                   channel_match = NA,
                                   axes_trans = NA,
                                   axes_limits = "machine",
                                   layout,
                                   header,
                                   point_col = c("magenta", 
                                                 "blue", 
                                                 "grey40"),
                                   point_col_alpha = 0.5,
                                   hist_fill = c("magenta", 
                                                 "blue", 
                                                 "grey40"),
                                   hist_fill_alpha = 0.5,
                                   lines = TRUE,
                                   line_type = 1,
                                   line_width = 2,
                                   line_col = "red",
                                   line_col_alpha = 1,
                                   text = TRUE,
                                   text_font = 1,
                                   text_size = 1.2,
                                   text_col = "black",
                                   text_col_alpha = 1,
                                   header_text_font = 2,
                                   header_text_size = 1,
                                   header_text_col = "black",
                                   ...) {
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # SIGNAL CALL TO CYTO_PLOT & RESET
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "compensation")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # PREPARE UNCOMPENSATED & COMPENSATED DATA -----------------------------------
  
  # SELECT & COPY
  x <- cyto_copy(
    cyto_select(
      x,
      select
    )
  )
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    channels <- channels[!grepl("-H$|-W$", channels, ignore.case = TRUE)]
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # OVERLAY
  if(!is.character(overlay)) {
    stop(
      "'overlay' must be either 'none', 'compensated', 'unstained' or 'both'!"
    )
  } else {
    if(any(grepl("^b", overlay, ignore.case = TRUE))) {
      overlay <- c("unstained", "compensated")
    }
  }
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # APPLIED TRANSFORMERS
  if(!.all_na(axes_trans)) {
    attributes(axes_trans)$applied <- TRUE
  }
  
  # DEFAULT TRANSFROMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_define(
      x,
      parent = "root",
      channels = channels,
      type = "biex",
      plot = FALSE
    )
    attributes(axes_trans)$applied <- FALSE
  }
  
  # APPLIED SPILLOVER MATRICES
  if(cyto_class(x, "GatingSet")) {
    spill <- cyto_spillover_extract(x)
  } else {
    if(compensated) {
      spill <- cyto_spillover_extract(x)
    } else {
      spill <- NULL
    }
  }
  
  # CHANNEL MATCH - CYTO_SPILLOVER_EDIT() BYPASS
  if(!.all_na(channel_match)) {
    warning(
      "'channel_match' is now reserved for internal use only!"
    )
    pd <- channel_match
  } else {
    pd <- cyto_channel_match(
      x,
      channels = channels
    )
  }
  pd <- pd[match(rownames(cyto_details(x)), rownames(pd)), , drop = FALSE]
  
  # UPDATE EXPERIMENT DETAILS
  cyto_details(x) <- pd
  
  # UNCOMPENSATED DATA PER GROUP
  x_uncomp <- cytoset(
    structure(
      lapply(
        seq_along(x), 
        function(z) {
          # PARENT POPULATION
          if("parent" %in% colnames(pd)) {
            pop <- pd$parent[z]
          } else {
            pop <- parent
          }
          print(pop)
          # EXTRACT LINEAR DATA
          cs <- cyto_data_extract(
            x[z],
            parent = pop,
            copy = FALSE,
            trans = if(attributes(axes_trans)$applied) {
              axes_trans
            } else {
              NA
            },
            inverse = TRUE
          )[[1]]
          # DECOMPENSATE
          if(!is.null(spill)) {
            cs <- cyto_compensate(
              cs,
              spillover = spill[[z]],
              remove = TRUE,
              copy = FALSE
            )
          }
          # RE-APPLY TRANSFORMERS
          cs <- cyto_transform(
            cs,
            trans = axes_trans,
            copy = FALSE,
            plot = FALSE,
            quiet = TRUE
          )
          # RETURN CYTOFRAME
          return(cs[[1]])
          
        }
      ),
      names = cyto_names(x)
    )
  )
  
  # PREPARE SPILLOVER MATRICES
  spillover <- .cyto_spillover_prepare(
    x_uncomp,
    spillover = spillover
  )
  
  # COMPENSATED DATA
  if(any(grepl("^c", overlay, ignore.case = TRUE))) {
    # INVERSE TRANSFORMATIONS
    x_comp <- cyto_transform(
      x_uncomp,
      trans = axes_trans,
      inverse = TRUE,
      copy = TRUE,
      plot = FALSE,
      quiet = TRUE
    )
    # APPLY COMPENSATION
    x_comp <- cyto_compensate(
      x_comp,
      spillover = spillover, # NULL - STORED SPILLOVER MATRIX
      copy = FALSE
    )
    # RE-APPLY TRANSFORMATIONS
    x_comp <- cyto_transform(
      x_comp,
      trans = axes_trans,
      inverse = FALSE,
      copy = FALSE,
      plot = FALSE,
      quiet = TRUE
    )
  } else {
    x_comp <- NULL
  }
  
  # EXPERIMENT DETAILS WITHOUT UNSTAINED
  pd_comp <- pd[
    !grepl("Unstained", pd$channel, ignore.case = TRUE), , drop = FALSE
  ]
  
  # PREPARE GRAPHICS DEVICE ----------------------------------------------------
  
  # LAYOUT
  if(missing(layout)) {
    layout <- .cyto_plot_layout(channels)
  }
  
  # PLOTS PER PAGE
  if(is.null(dim(layout))) {
    np <- prod(layout)
  } else {
    np <- length(
      unique(
        unlist(
          layout
        )
      )
    )
  }
  
  # SAMPLES
  n <- nrow(pd_comp)
  
  # PAGES PER SAMPLE
  pg <- ceiling(length(channels)/np)
  
  # TOTAL PAGES
  tpg <- n * pg
  
  # HEADER
  if(missing(header)) {
    header <- rep(
      cyto_names(
        cyto_select(
          x_uncomp,
          pd_comp$name
        )
      ),
      each = pg
    )
  } else {
    header <- rep(
      header,
      length.out = nrow(pd_comp) * pg
    )
  }
  
  # REPEAT HEADER ARGUMENTS
  header_text_font <- rep(header_text_font, length.out = length(header))
  header_text_size <- rep(header_text_size, length.out = length(header))
  header_text_col <- rep(header_text_col, length.out = length(header))
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
    
  # PLOTS PER SAMPLE
  plots <- structure(
    # PAGE PER SAMPLE
    lapply(
      seq_along(pd_comp$name),
      function(z) {
        # SAMPLE INDEX
        x_ind <- cyto_match(
          x_uncomp,
          pd_comp$name[z],
          exact = TRUE
        )
        # GROUP
        grp <- pd_comp$group[z]
        # SEARCH FOR UNSTAINED CONTROL WITHIN GROUP
        unst_ind <- which(
          pd$group == grp & grepl("Unstained", pd$channel, ignore.case = TRUE)
        )
        # UNSTAINED CONTROL LOCATED
        if(length(unst_ind) > 0 & 
           any(grepl("^u", overlay, ignore.case = TRUE))) {
          NIL <- cyto_select(
            x_uncomp,
            pd$name[unst_ind[1]] # USE FIRST UNSTAINED CONTROL
          )
        } else {
          NIL <- NULL
        }
        # HEADER COUNTER
        cnt <- (z - 1) * pg
        # X CHANNEL
        xchan <- pd_comp$channel[z]
        # PLOT FOR EACH CHANNEL
        p <- structure(
          lapply(
            seq_along(channels),
            function(w) {
              # Y CHANNEL
              ychan <- channels[w]
              # OVERLAY - COMPENSATED | UNSTAINED
              if(is.null(x_comp) & is.null(NIL)) {
                overlay <- NA
              } else {
                overlay <- list(
                  x_comp[x_ind],
                  NIL
                )
                overlay[LAPPLY(overlay, "is.null")] <- NULL
              }
              # PLOT - HISTOGRAM | SCATTER
              cyto_plot(
                x_uncomp[x_ind],
                channels = if(xchan == ychan) {
                  xchan
                } else {
                  c(xchan, ychan)
                },
                overlay = overlay,
                axes_trans = axes_trans,
                axes_limits = axes_limits,
                point_col = point_col,
                point_col_alpha = point_col_alpha,
                hist_fill = hist_fill,
                hist_fill_alpha = hist_fill_alpha,
                header = "", # HEADERS MUST BE MANUALLY ADDED
                title = NA,
                layout = layout,
                page = FALSE, # MANUAL PAGING REQUIRED
                ...
              )
              # SPILLOVER COEFFICIENTS | MODELS
              if(xchan != ychan) {
                # LINES
                if(lines) {
                  # ROBUST LINEAR MODEL
                  rlm <- cyto_apply(
                    if(is.null(x_comp)) {
                      x_uncomp[x_ind]
                    } else {
                      x_comp[x_ind]
                    },
                    function(q, chans = NULL){
                      rlm(q[, chans[2]] ~ q[, chans[1]])
                    },
                    input = "matrix",
                    channels = c(xchan, ychan),
                    copy = FALSE,
                    simplify = FALSE,
                    chans = c(xchan, ychan)
                  )[[1]]
                  # PLOT LINEAR MODEL
                  abline(
                    rlm,
                    lty = line_type,
                    lwd = line_width,
                    col = adjustcolor(line_col, line_col_alpha),
                  )
                }
                # SPILLOVER COEEFICIENTS
                if(!is.null(spillover[[x_ind]]) & text & 
                   all(c(xchan, ychan) %in% colnames(spillover[[x_ind]]))) {
                  text(
                    labels = paste0(
                      "spill = ", 
                      round(
                        spillover[[x_ind]][
                          match_ind(xchan, 
                                    colnames(spillover[[x_ind]])), 
                          match_ind(ychan, 
                                    colnames(spillover[[x_ind]]))
                        ] * 100, 2), 
                      "%"
                    ),
                    x = mean(.par("usr")[[1]][1:2]),
                    y = .par("usr")[[1]][3] + 0.9 * diff(.par("usr")[[1]][3:4]),
                    cex = text_size,
                    font = text_font,
                    col = adjustcolor(text_col, text_col_alpha)
                  )
                }
              }
              # PAGE & HEADER & RECORD
              rec <- NULL
              if(.par("page")[[1]] | w == length(channels)) {
                # HEADER INDEX
                ind <- cnt + ceiling(w/np)
                # HEADER
                if(!.all_na(header[ind])) {
                  .cyto_plot_header(
                    header[ind],
                    header_text_font = header_text_font[ind],
                    header_text_size = header_text_size[ind],
                    header_text_col = header_text_col[ind]
                  )
                }
                # NEW PAGE
                cyto_plot_new_page()
                # RECORD
                rec <- cyto_plot_record()
              }
              return(rec)
            }
          ),
          names = channels
        )
      }
    ),
    names = pd_comp$name
  )
  
  # RECORDED PLOTS -------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}