## CYTO_PLOT_PROFILE -----------------------------------------------------------

#' Plot expression profile in all channels
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels a vector channels to use to construct the plots, set to all
#'   channels by default.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied.
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which
#'   was used to transform the channels of the supplied flowFrame.
#'   \code{cyto_plot} does not support in-line transformations and as such the
#'   transformations should be applied to the data prior to plotting. The
#'   transformerList is used internally to ensure that the axes on the
#'   constructed plots are appropriately labelled.
#' @param group_by a vector of experiment variables to sort and merge samples
#'   into groups prior to plotting, set to NA by default to prevent merging.
#'   To merge all samples set this argument to \code{"all"}.
#' @param select designates which samples will be plotted and used for
#'   determining the best location to set the drawn gate(s). Filtering steps
#'   should be comma separated and wrapped in a list. Refer to
#'   \code{\link{cyto_select}}.
#' @param layout a vector of the length 2 indicating the dimensions of the grid
#'   for plotting \code{c(#rows, #columns)}.
#' @param header character string to include in the plot header, set to the
#'   sample names by deafult. The header can be removed by setting this argument
#'   to NA.
#' @param header_text_font numeric indicating the font to use for the header,
#'   set to 2 by default for bold font.
#' @param header_text_size numeric to control the text size of the header, set
#'   to 1 by deafult.
#' @param density_stack numeric [0,1] indicating the degree of offset for 1-D
#'   density distributions with overlay, set to 0.5 by default.
#' @param ... additional arguments passed to \code{\link{cyto_plot}}.
#'
#' @importFrom grDevices n2mfrow recordPlot
#' @importFrom graphics par mtext
#' @importFrom methods is
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @seealso \code{\link{cyto_plot}}
#'
#' @examples
#' 
#' # Load in CytoExploreR to access data
#' library(CytoExploreRData)
#'
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#'
#' # Apply compensation
#' gs <- cyto_compensate(gs, fs[[1]]@description$SPILL)
#'
#' # Transform fluorescent channels
#' gs <- cyto_transform(gs, trans_type = "logicle")
#'
#' # Apply gates
#' gt <- Activation_gatingTemplate
#' gt_gating(gt, gs)
#'
#' # Plot expression profile of T Cells in all channels
#' cyto_plot_profile(gs[1:9],
#'   parent = "T Cells"
#' )
#' 
#' # Group samples by Treatment & select highest OVAConc
#'  cyto_plot_profile(gs[1:9],
#'   parent = "CD4 T Cells",
#'   group_by = "Treatment",
#'   select = list("OVAConc" = 500),
#' )
#' 
#' @name cyto_plot_profile
NULL

#' @noRd
#' @export
cyto_plot_profile <- function(x, ...) {
  UseMethod("cyto_plot_profile")
}

#' @rdname cyto_plot_profile
#' @export
cyto_plot_profile.GatingSet <- function(x,
                                        parent = NULL,
                                        channels = NULL,
                                        group_by = NA,
                                        select = NULL,
                                        axes_trans = NA,
                                        layout = NULL,
                                        header = NULL,
                                        header_text_font = 2,
                                        header_text_size = 1,
                                        ...) {

  # CHECKS ---------------------------------------------------------------------

  # PLOT METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "profile/GatingSet")
  }

  # PARENT
  if (is.null(parent)) {
    stop("Please supply the name of the parent population to plot.")
  }

  # TRANSFORMATIONS
  axes_trans <- cyto_transformer_extract(x)
  
  # PREPARE DATA ---------------------------------------------------------------

  # EXTRACT PARENT
  fs <- cyto_extract(x, parent)

  # Plots
  p <- cyto_plot_profile(
    fs,
    channels = channels,
    group_by = group_by,
    select = select,
    axes_trans = axes_trans,
    layout = layout,
    header = header,
    header_text_font = header_text_font,
    header_text_size = header_text_size,
    ...
  )

  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save")) {
    if (is(x, basename(getOption("cyto_plot_method")))) {
      # CLOSE GRAPHICS DEVICE
      dev.off()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
      # RESET CYTO_PLOT_METHOD
      options("cyto_plot_method" = NULL)
    }
  }
  
  # RETURN RECORDED PLOTS
  invisible(p)
  
}

#' @rdname cyto_plot_profile
#' @export
cyto_plot_profile.GatingHierarchy <- function(x,
                                             parent = NULL,
                                             channels = NULL,
                                             axes_trans = NA, 
                                             header = NULL,
                                             header_text_font = 2,
                                             header_text_size = 1, ...) {

  # CHECKS ---------------------------------------------------------------------

  # PLOT METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "profile/GatingHierarchy")
  }

  # PARENT
  if (is.null(parent)) {
    stop("Please supply the name of the parent population to plot.")
  }

  # TRANSFORMATIONS
  axes_trans <- cyto_transformer_extract(x)
  
  # PREPARE DATA ---------------------------------------------------------------

  # EXTRACT PARENT
  fr <- cyto_extract(x, parent)

  # Plots
  p <- cyto_plot_profile(
    fr,
    channels = channels,
    axes_trans = axes_trans,
    header = header,
    header_text_font = header_text_font,
    header_text_size = header_text_size, ...
  )
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save")) {
    if (is(x, basename(getOption("cyto_plot_method")))) {
      # CLOSE GRAPHICS DEVICE
      dev.off()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
      # RESET CYTO_PLOT_METHOD
      options("cyto_plot_method" = NULL)
    }
  }
  
  # RETURN RECORDED PLOTS
  invisible(p)
  
}

#' @rdname cyto_plot_profile
#' @export
cyto_plot_profile.flowSet <- function(x,
                                      channels = NULL,
                                      group_by = NA,
                                      select = NA,
                                      axes_trans = NULL,
                                      layout = NULL,
                                      header = NULL,
                                      header_text_font = 2,
                                      header_text_size = 1,
                                      density_stack = 0.5, ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # PLOT METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "profile/flowSet")
  }
  
  # PREPARE SAMPLES ------------------------------------------------------------
  
  # GROUP_BY
  if(!.all_na(group_by)){
    fs_list <- cyto_group_by(x, group_by)
  }else{
    fs_list <- list(x)
  }
  
  # SELECT
  if(!.all_na(select)){
    fs_list <- lapply(fs_list, function(z) {
       tryCatch(cyto_select(z, select), error = function(e) {
        z
      })
    })
  }

  # COLLAPSE GROUPS - FLOWFRAME PER GROUP
  if(!.all_na(group_by)){
    fr_list <- lapply(fs_list, function(z){
      cyto_convert(z, "flowFrame")
    })
    names(fr_list) <- group_by
  # FLOWSET TO FLOWFRAME LIST
  }else{
    fr_list <- cyto_convert(fs_list[[1]], "list of flowFrames")
    names(fr_list) <- cyto_names(x)
  }
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # TITLE
  if (is.null(header)) {
    header <- "Expression Profile"
  }

  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # CYTO_PLOT_PROFILE
  if (length(fr_list) == 1) {
    p <- cyto_plot_profile(
      x = fr_list[[1]],
      channels = channels,
      axes_trans = axes_trans,
      layout = layout,
      header = header,
      header_text_font = header_text_font,
      header_text_size = header_text_size, ...
    )
  } else {
    p <- cyto_plot_profile(
      x = fr_list[[1]],
      overlay = fr_list[2:length(fr_list)],
      channels = channels,
      axes_trans = axes_trans,
      layout = layout,
      header = header,
      header_text_font = header_text_font,
      header_text_size = header_text_size,
      density_stack = density_stack, ...
    )
  }
  
  # RECORD/SAVE ----------------------------------------------------------------
  
  # TURN OFF GRAPHICES DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save")) {
    if (is(x, basename(getOption("cyto_plot_method")))) {
      # CLOSE GRAPHICS DEVICE
      dev.off()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
      # RESET CYTO_PLOT_METHOD
      options("cyto_plot_method" = NULL)
    }
  }
  
  # RETURN RECORDED PLOT
  invisible(p)
  
}

#' @rdname cyto_plot_profile
#' @export
cyto_plot_profile.flowFrame <- function(x,
                                        channels = NULL,
                                        axes_trans = NA,
                                        layout = NULL,
                                        header = NULL,
                                        header_text_font = 2,
                                        header_text_size = 1,
                                        ...) {

  # CHECKS ---------------------------------------------------------------------
  
  # PLOT METHOD
  if (is.null(getOption("cyto_plot_method"))) {
    options("cyto_plot_method" = "profile/flowFrame")
  }

  # CHANNELS
  if (is.null(channels)) {
    channels <- cyto_channels(x, exclude = "Time")
  }

  # PREPARE CHANNELS
  channels <- cyto_channels_extract(x, channels = channels, plot = FALSE)

  # GRAPHICAL PARAMETERS -------------------------------------------------------
  
  # OLD PARAMETERS
  old_pars <- .par(c("mfrow","oma"))
  
  # RESET PARAMETERS ON EXIT
  on.exit(par(old_pars))
  
  # PREPARE PLOTTING SAPCE -----------------------------------------------------
  
  # LAYOUT DIMENSIONS
  if (is.null(layout)) {
    layout <- c(n2mfrow(length(channels))[2], n2mfrow(length(channels))[1])
    par(mfrow = layout)
  } else if (!is.null(layout)) {
    if (layout[1] == FALSE) {
      # DO NOTHING
    } else {
      par(mfrow = layout)
    }
  }

  # HEADER SPACE
  if (!.all_na(header)) {
    par(oma = c(0, 0, 3, 0))
  }
  
  # CONSTRUCT PLOTS ------------------------------------------------------------
  
  # NUMBER OF PLOTS
  NP <- length(channels)
  NP <- split(seq_len(NP), ceiling(seq_len(NP)/prod(layout)))
  NP <- LAPPLY(NP, "max")
  
  # HEADER SAMPLENAMES
  if (is.null(header)) {
    header <- cyto_names(x)
  } 
  
  # CONSTRUCT PLOTS
  p <- lapply(seq_len(length(channels)), function(z) {
    
    # PLOT
    cyto_plot(x,
      channels = channels[z],
      axes_trans = axes_trans,
      title = NA, ...
    )
    
    # HEADER
    if(z %in% NP){
      # ADD HEADER
      if (!.all_na(header)) {
        mtext(header, 
              outer = TRUE, 
              font = header_text_font,
              cex = header_text_size)
      }
      # RECORD PLOT
      cyto_plot_record()
    }
  })
  
  # RECORD/SAVE ----------------------------------------------------------------

  # TURN OFF GRAPHICS DEVICE - CYTO_PLOT_SAVE
  if (getOption("cyto_plot_save")) {
    if (is(x, basename(getOption("cyto_plot_method")))) {
      # CLOSE GARPHICS DEVICE
      dev.off()
      # RESET CYTO_PLOT_SAVE
      options("cyto_plot_save" = FALSE)
      # RESET CYTO_PLOT_METHOD
      options("cyto_plot_method" = NULL)
    }
  }
  
  # RETURN RECORDED PLOTS
  invisible(p)
}
