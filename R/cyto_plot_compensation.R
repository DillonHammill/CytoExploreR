## CYTO_PLOT_COMPENSATION ------------------------------------------------------

#' Visualise compensation of spillover in all channels
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing data for
#'   single colour compensation controls.
#' @param parent names of the parent populations to extract from each control
#'   for plotting.
#' @param channel_match a vector of channels one per control or the name of a
#'   CSV file which indicates the name of the channel associated with each
#'   control. If \code{channel_match} is not supplied or the channels are not
#'   annotated in \code{cyto_details(x)} the user will be asked to interactively
#'   enter this information prior to plotting.
#' @param spillover a matrix or name of a CSV file containing the spillover
#'   coefficients that should be applied to the data when \code{compensate} is
#'   TRUE.
#' @param compensate logical required for \code{cytoset} objects to indicate
#'   whether the data should be compensated prior to plotting, set to TRUE by
#'   default. This argument should be set to FALSE if the \code{spillover} has
#'   already been applied to the data.
#' @param overlay can be either \code{"unstained"}, \code{"compensated"},
#'   \code{"both"} or \code{"none"} to allow control over which data is overlaid
#'   onto the plots, set to \code{"both"} by default.
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
                                   channel_match = NULL,
                                   spillover = NULL, 
                                   compensate = TRUE,
                                   overlay = "both",
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
                                   ...) {
  
  # TODO: SUPPORT COMPENSATED DATA AND MATRIX MANUALLY SUPPLIED
  # REMOVE CURRENT MATRIX & APPLY NEW MATRIX
  
  # CYTO_PLOT_COMPLETE ---------------------------------------------------------
  
  # SIGNAL CALL TO CYTO_PLOT & RESET
  if(is.null(cyto_option("cyto_plot_method"))) {
    # SET CYTO_PLOT_METHOD
    cyto_option("cyto_plot_method", "comp")
    # CYTO_PLOT_EXIT
    on.exit({
      cyto_plot_complete()
    })
  }
  
  # PREPARE COMPENSATED & UNCOMPENSATED DATA -----------------------------------
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # APPLIED TRANSFORMERS
  if(!.all_na(axes_trans)){
    attributes(axes_trans)$applied <- TRUE
  }
  
  # APPLIED SPILLOVER MATRICES
  if(cyto_class(x, "GatingSet")) {
    spill <- cyto_spillover_extract(x)[[1]]
  } else {
    spill <- NULL
  }
  
  # OVERLAY - COMPENSATED DATA
  if(cyto_class(overlay, "cytoset")) {
    x_comp <- .cyto_spillover_controls_prepare(
      overlay,
      type = "bagwell",
      channel_match = channel_match,
      axes_trans = axes_trans
    )
    overlay <- "compensated"
  } else {
    x_comp <- NULL
  }
  
  # FORMAT OVERLAY
  if(overlay == "both") {
    overlay <- c("compensated", "unstained")
  } else if(.all_na(overlay)) {
    overlay <- "none"
  }
  
  # PREPARE DATA - GATINGSET CANNOT BE INVERSE TRANSFORMED
  x_uncomp <- .cyto_spillover_controls_prepare(
    x,
    parent = parent,
    type = "bagwell",
    channel_match = channel_match,
    axes_trans = axes_trans
  )
  
  # DEFAULT TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_define(x,
                                           parent = parent[[1]],
                                           type = "logicle",
                                           plot = FALSE)
    attributes(axes_trans)$applied <- FALSE
  }
  
  # DATA IS UNCOMPENSATED
  if(is.null(spill)) {
    # USE SUPPLIED SPILLOVER MATRIX
    if(!is.null(spillover)) {
      if(is.null(dim(spillover))) {
        spill <- read_from_csv(spillover)
      } else {
        spill <- spillover
      }
    }
    # PREPARE COMPENSATED DATA - BYPASS OVERLAY - MANUALLY SUPPLIED
    if(is.null(x_comp) & "compensated" %in% overlay) {
      # SPILLOVER REQUIRED
      if(is.null(spill)) {
        stop(
          "Supply a matrix/file to 'spillover' to overlay compensated data!"
        )
      # SPILLOVER SUPPLIED
      } else {
        # PREPARE COMPENSATED DATA
        x_comp <- structure(
          lapply(seq_along(x_uncomp), function(z){
            structure(
              lapply(x_uncomp[[z]], function(w){
                # COPY DATA - UNLINK FROM X_UNCOMP
                w <- cyto_copy(w)
                # UNCOMPENSATED DATA IS TRANSFORMED
                if(attributes(axes_trans)$applied) {
                  # INVERSE TRANSFORMATIONS
                  w <- cyto_transform(w,
                                      trans = axes_trans,
                                      channels = names(x_uncomp),
                                      inverse = TRUE,
                                      copy = FALSE,
                                      plot = FALSE,
                                      quiet = TRUE)
                }
                # COMPENSATE
                w <- cyto_compensate(w,
                                     spillover = spill,
                                     remove = FALSE,
                                     copy = FALSE)
                # APPLY TRANSFORMATIONS
                w <- cyto_transform(w,
                                    trans = axes_trans,
                                    channels = names(x_uncomp),
                                    inverse = FALSE,
                                    copy = FALSE,
                                    plot = FALSE,
                                    quiet = TRUE)
                return(w)
              }),
              names = names(x_uncomp[[z]])
            )
          }),
          names = names(x_uncomp)
        )
      }
    }
  # DATA IS COMPENSATED
  } else {
    # ASSIGN TO COMPENSATED DATA
    if(is.null(x_comp)) {
      x_comp <- x_uncomp
    }
    # PREPARE UNCOMPENSATED DATA
    x_uncomp <- structure(
      lapply(seq_along(x_comp), function(z){
        structure(
          lapply(x_comp[[z]], function(w){
            # COPY DATA - UNLINK FROM X_COMP
            w <- cyto_copy(w)
            # COMPENSATED DATA IS TRANSFORMED
            if(attributes(axes_trans)$applied) {
              # INVERSE TRANSFORMATIONS
              w <- cyto_transform(w,
                                  trans = axes_trans,
                                  channels = names(x_uncomp),
                                  inverse = TRUE,
                                  copy = FALSE,
                                  plot = FALSE,
                                  quiet = TRUE)
            }
            # DECOMPENSATE
            w <- cyto_compensate(w,
                                 spillover = spill,
                                 remove = TRUE,
                                 copy = FALSE)
            # APPLY TRANSFORMATIONS
            w <- cyto_transform(w,
                                trans = axes_trans,
                                channels = names(x_uncomp),
                                inverse = FALSE,
                                copy = FALSE,
                                plot = FALSE,
                                quiet = TRUE)
            return(w)
          }),
          names = names(x_comp[[z]])
        )
      }),
      names = names(x_comp)
    )
  }
  
  # COMPENSATED DATA NOT REQUIRED
  if(!"compensated" %in% overlay) {
    x_comp <- NA
  }
  
  # COMBINE DATA
  x <- structure(
    lapply(seq_along(x_uncomp), function(z){
      res <- list()
      # UNCOMPENSATED
      res[["uncompensated"]] <- x_uncomp[[z]][["+"]]
      # COMPENSATED
      if(!.all_na(x_comp)) {
        res[["compensated"]] <- x_comp[[z]][["+"]]
      } else {
        res[["compensated"]] <- NA
      }
      # UNSTAINED - UNCOMEPNSATED
      if(!"compensated" %in% overlay) {
        if(!is.null(x_uncomp[[z]][["-"]])) {
          res[["unstained"]] <- x_uncomp[[z]][["-"]]
        } else {
          res[["unstained"]] < NA
        }
      # UNSTAINED COMPENSATED
      } else {
        if(!is.null(x_comp[[z]][["-"]])) {
          res[["unstained"]] <- x_comp[[z]][["-"]]
        } else {
          res[["unstained"]] < NA
        }
      }
      return(res)
    }),
    names = names(x_uncomp)
  )
  
  # PREPARE GRAPHICS DEVICE ----------------------------------------------------
  
  # LAYOUT
  if(missing(layout)) {
    layout <- .cyto_plot_layout(x)
  }
  
  # HEADER
  if(missing(header)) {
    header <- LAPPLY(x, function(z){
      cyto_names(z[["uncompensated"]])
    })
  }
  header <- rep(header, length.out = length(x))
  names(header) <- names(x)
  
  # HEADER SPACE
  if(!.all_na(header)) {
    oma <- c(0,0,3,0)
  } else {
    oma <- .par("oma")[[1]]
  }

  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # PLOTS PER CONTROL
  plots <- structure(
    # PAGE PER CONTROL
    lapply(names(x), function(z){
      # PLOT FOR EACH CHANNEL
      p <- lapply(names(x), function(v){
        # HISTOGRAM/SCATTER
        cyto_plot(
          x[[z]][["uncompensated"]],
          channels = if(z == v) {
            z
          } else {
            c(z, v)
          },
          overlay = if(!.all_na(overlay)) {
            x[[z]][overlay]
          } else {
            NA
          },
          axes_trans = axes_trans,
          axes_limits = axes_limits,
          point_col = point_col,
          point_col_alpha = point_col_alpha,
          hist_fill = hist_fill,
          hist_fill_alpha = hist_fill_alpha,
          header = header[z],
          title = NA,
          layout = layout,
          page = if(v == names(x)[length(x)]) {
            TRUE
          } else {
            FALSE
          },
          ...
        )
        # SPILLOVER COEFFICIENTS/LINES
        if(z != v) {
          # LINES
          if(lines) {
            # ROBUST LINEAR MODEL
            rlm <- cyto_apply(
              if("compensated" %in% overlay) {
                x[[z]][["compensated"]]
              } else {
                x[[z]][["uncompensated"]]
              },
              function(q, chans = NULL){
                rlm(q[, chans[2]] ~ q[, chans[1]])
              },
              input = "matrix",
              channels = c(z, v),
              copy = FALSE,
              simplify = FALSE,
              chans = c(z, v)
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
          if(!is.null(spill) & text) {
            text(
              labels = paste0("sp = ", round(spill[z, v]*100, 2), "%"),
              x = mean(.par("usr")[[1]][1:2]),
              y = .par("usr")[[1]][3] + 0.9 * diff(.par("usr")[[1]][3:4]),
              cex = text_size,
              font = text_font,
              col = adjustcolor(text_col, text_col_alpha)
            )
          }
        }
        # RECORD PLOTS
        if(.par("page")[[1]]) {
          return(cyto_plot_record())
        } else {
          return(NULL)
        }
      })
      # REMOVE NULL PLOTS
      p[LAPPLY(p, "is.null")] <- NULL
      return(p)
    }),
    names = LAPPLY(x, function(w){
      cyto_names(w[["uncompensated"]])
    })
  )
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # RETURN RECORDED PLOTS
  invisible(plots)
  
}
