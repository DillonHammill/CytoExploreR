## CYTO_PLOT_LABEL -------------------------------------------------------------

#' Add labels to cyto_plot
#'
#' \code{cyto_plot_label()} prepares the data in the same way as
#' \code{cyto_plot()} by allows uses to interactively add labels to plots that
#' have already been created.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the population to plot when a \code{GatingHierarchy} or
#'   \code{GatingSet} object is supplied, set to the \code{"root"} node by
#'   default.
#' @param alias name of the gated population(s) to gated in the plot when a
#'   \code{GatingHierarchy} or \code{GatingSet} object is supplied. Setting
#'   \code{alias} to "" will automatically plot any gates constructed in the
#'   supplied channels. \code{alias} is equivalent to the \code{gate} argument
#'   for \code{cytoset} objects.
#' @param channels name of the channel(s) or marker(s) to be used to construct
#'   the plot.
#' @param axes_trans object of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} which was used
#'   to transform the channels of the supplied data.  The transformerList is
#'   used internally to ensure that the data is converted to a linear scale
#'   (when required) prior to computing statistics for labels.
#' @param events numeric to control the number or percentage of events to
#'   display. Values [0,1] indicate the percentage of events to display (i.e.
#'   value of 1 will display all events), whilst values larger than 1 indicate
#'   the number of events to display. The default value for \code{events} is set
#'   to 50000 to display 50000 events only.
#' @param negate logical indicating whether a label should be included for the
#'   negated population when gate objects are supplied, set to FALSE by default.
#'   Setting \code{negate = TRUE} will result in the creation of a boolean
#'   filter that contains all the events outside the gates supplied to
#'   \code{alias} or \code{gate}. If such a boolean gate exists in the
#'   GatingSet/GatingHierarchy it will automatically be extracted when
#'   \code{alias = ""}. In order to prevent plotting of these boolean gates,
#'   users will need to explicitly pass the names of the gates they want to
#'   display to \code{alias}.
#' @param gate gate objects to be apply to data prior to computing statistics
#'   for labels, can be either objects of class \code{rectangleGate},
#'   \code{polygonGate}, \code{ellipsoidGate}, \code{quadGate} or
#'   \code{filters}. Lists of these supported gate objects are also supported.
#' @param overlay name(s) of the populations to overlay or a \code{cytoset},
#'   \code{list of cytosets} or \code{list of cytoset lists} containing
#'   populations to be overlaid onto the plot(s). This argument can be set to
#'   "children" or "descendants" when a \code{GatingSet} or
#'   \code{GatingHierarchy} to overlay all respective nodes.
#' @param merge_by a vector of pData variables to sort and merge samples into
#'   groups prior to plotting, set to NA by default to prevent merging. To
#'   merge all samples set this argument to \code{TRUE} or \code{"all"}.
#' @param label_text vector of population names to use in the labels. To exclude
#'   the population names set this argument to NA.
#' @param label_stat indicates the type of statistic to include in the plot
#'   labels, can be \code{"percent"}, \code{"count"}, \code{"mean"},
#'   \code{"median"}, \code{"mode"} or \code{"geo mean"}, set to
#'   \code{"percent"} for gated data or \code{NA} to exclude statistics for
#'   un-gated data. Currently, only \code{"percent"} and \code{"count"} are
#'   supported for 2-D scatter plots.
#' @param label_text_x vector of x co-ordinate(s) to manually adjust the
#'   position plot label(s) on the plot. To interactively position labels set
#'   either \code{label_text_x} or \code{label_text_y} to "select".
#' @param label_text_y vector of y co-ordinate(s) to manually adjust the
#'   position plot label(s) on the plot. To interactively position labels set
#'   either \code{label_text_x} or \code{label_text_y} to "select".
#' @param label_text_font numeric to control the font of text in plot labels,
#'   set to 2 for bold font by default. See \code{\link[graphics:par]{font}} for
#'   alternatives.
#' @param label_text_size numeric to control the size of text in the plot
#'   labels, set to 1 by default.
#' @param label_text_col colour(s) to use for text in plot labels, set to
#'   \code{"black"} by default.
#' @param label_text_col_alpha numeric [0, 1] to control the transparency of the
#'   text colour, set to 1 by default to remove transparency.
#' @param label_fill fill colour(s) to use for labels, set to "white" by
#'   default.
#' @param label_fill_alpha numeric to control background fill transparency of
#'   label, set to 0.6 by default to introduce some transparency.
#' @param hist_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}} to adjust the smoothness of the kernel
#'   density for histograms, set to \code{1} by default. Only values greater or
#'   equal to 1 are supported.
#' @param ... not in use.
#'
#' @return a matrix containing the xy co-ordinates of the labels.
#'
#' @importFrom grDevices adjustcolor
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @examples
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#'
#' # Call cyto_plot()
#' cyto_plot(gs[1], parent = "root", channels = c("FSC-A", "SSC-A"))
#'
#' # Add label after plotting
#' cyto_plot_label(
#'   gs[1],
#'   parent = "root",
#'   channels = c("FSC-A", "SSC-A"),
#'   label_stat = "count",
#'   label_text = "count",
#'   label_text_x = 150000,
#'   label_text_y = 150000
#' )
#'
#' @export
cyto_plot_label <- function(x,
                            parent = "root",
                            alias = NA,
                            channels = NULL,
                            axes_trans = NA,
                            events = 50000,
                            negate = FALSE,
                            gate = NA,
                            overlay = NA,
                            merge_by = "name",
                            label_text,
                            label_stat,
                            label_text_x = NA,
                            label_text_y = NA,
                            label_text_font = 2,
                            label_text_size = 0.8,
                            label_text_col = "black",
                            label_text_col_alpha = 1,
                            label_fill = "white",
                            label_fill_alpha = 0.6,
                            hist_smooth = 1,
                            ...) {
  
  # ARGUMENTS ------------------------------------------------------------------
  
  # PULL DOWN ARGUMENTS
  args <- .args_list(...)
  
  # INHERIT THEME
  args <- .cyto_plot_theme_inherit(args)
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # TRANSFORMERS
  if(.all_na(args$axes_trans)) {
    args$axes_trans <- cyto_transformers_extract(args$x)
  }
  
  # CONVERT CHANNELS
  args$channels <- cyto_channels_extract(args$x[[1]], args$channels)
  
  # GATE - .CYTO_PLOT_GATES() -> LIST OF GATE OBJECTS
  args$gate <- cyto_func_execute(".cyto_plot_gates", args)[[1]]
  
  # DATA - .CYTO_PLOT_DATA() -> LIST OF MERGED CYTOSETS PER LAYER
  args$x <- cyto_func_execute(".cyto_plot_data", args)[[1]]
  
  # PREPARE GATES
  if(!.all_na(args$gate)) {
    args$gate <- cyto_gate_prepare(args$gate, args$channels)
  }
  
  # PREPARE LABEL ARGUMENTS
  args$label <- TRUE
  args[c("label", "label_text", "label_stat")] <- cyto_func_execute(
    ".cyto_plot_label_args",
    args
  )
  
  # COMPUTE STATISTICS FOR POPULATIONS -----------------------------------------
  
  # POPULATIONS TO LABEL
  args$pops <- cyto_func_execute(
    ".cyto_plot_label_pops",
    args
  )
  
  # COMPUTE STATISTICS
  args$label_stat <- cyto_func_execute(
    ".cyto_plot_label_stat",
    args
  )
  
  # PREPARE LABELS
  args$label_text <- cyto_func_execute(
    ".cyto_plot_label_text",
    args
  )
  
  # UPDATE ARGUMENTS
  .args_update(args)
  
  # LABEL CONSTRUCTION ---------------------------------------------------------
  
  # LABELS ADDED MANUALLY
  label_text_xy <- mapply(
    function(label_text,
             label_text_x,
             label_text_y,
             label_text_font,
             label_text_size,
             label_text_col,
             label_text_col_alpha,
             label_fill,
             label_fill_alpha) {
      
      # SELECT MISSING CO-ORDINATES
      if (any(.all_na(c(label_text_x, label_text_y)))) {
        if(grepl("\n", label_text)){
          txt <- paste(strsplit(label_text, "\n")[[1]], collapse = " ")
        }else{
          txt <- label_text
        }
        message(paste("Select a location on the plot the position the",
                      txt, "label."))
        label_text_xy <- locator(n = 1)
        label_text_x <- label_text_xy[[1]]
        label_text_y <- label_text_xy[[2]]
      }
      
      # PLOT LABELS
      .boxed.labels(
        x = label_text_x,
        y = label_text_y,
        labels = label_text,
        border = FALSE,
        font = label_text_font,
        cex = label_text_size,
        col = adjustcolor(label_text_col, label_text_col_alpha),
        bg = label_fill,
        alpha.bg = label_fill_alpha
      )
      
      # RETURN LABEL CO-ORDINATES
      return(c(label_text_x, label_text_y))
      
    },
    label_text,
    label_text_x,
    label_text_y,
    label_text_font,
    label_text_size,
    label_text_col,
    label_text_col_alpha,
    label_fill,
    label_fill_alpha,
    SIMPLIFY = FALSE
  )
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL CO-ORDINATES IN MATRIX
  label_text_xy <- do.call("rbind", label_text_xy)
  colnames(label_text_xy) <- c("x","y")
  label_text_xy <- as.matrix(label_text_xy)
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)
  
}

## CYTO_PLOT_LABELLER ----------------------------------------------------------

#' Add labels to existing cyto_plot
#'
#' Convenient labelling function to add prepared text labels to an existing
#' cyto_plot.
#'
#' @param label_text character string to include in the label.
#' @param label_text_x vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param label_text_y vector containing the x co-ordinates for the plot labels.
#'   Label positions can be interactively selected if no co-ordinates are
#'   manually supplied.
#' @param label_text_font integer [1,2,3,4] passed to \code{text} to alter the
#'   font, set to \code{2} by default for a bold font.
#' @param label_text_size numeric character expansion used to control the size
#'   of the text in the labels, set to \code{0.8} by default. See \code{?text}
#'   for details.
#' @param label_text_col specify text colour in label for each gate, defaults to
#'   \code{"black"} for all gates.
#' @param label_text_col_alpha numeric [0, 1] to control the transparency of the
#'   text colour, set to 1 by default to remove transparency.
#' @param label_fill fill colour to use for labels, set to "white" by default.
#' @param label_fill_alpha numeric [0,1] controls the transparency of the fill
#'   colour, set to \code{0.6} by default.
#' @param ... not in use.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom grDevices adjustcolor
#'
#' @export
cyto_plot_labeller <- function(label_text = NA,
                               label_text_x = NA,
                               label_text_y = NA,
                               label_text_font = 2,
                               label_text_size = 0.8,
                               label_text_col = "black",
                               label_text_col_alpha = 1,
                               label_fill = "white",
                               label_fill_alpha = 0.6,
                               ...){
  
  # LABEL CONSTRUCTION ---------------------------------------------------------
  
  # MANUAL CO-ORDINATES - SUPPLIED OR SELECTED
  label_text_xy <- mapply(
    function(label_text,
             label_text_x,
             label_text_y,
             label_text_font,
             label_text_size,
             label_text_col,
             label_text_col_alpha,
             label_fill,
             label_fill_alpha) {
      
      # MISSING LABEL
      if(!.all_na(label_text)){
        
        # SELECT MISSING CO-ORDINATES
        if (any(.all_na(c(label_text_x, label_text_y)))) {
          if(grepl("\n", label_text)){
            txt <- paste(strsplit(label_text, "\n")[[1]], collapse = " ")
          }else{
            txt <- label_text
          }
          message(
            paste(
              "\n Select a location on the plot the position the",
              txt,
              "label."
            )
          )
          label_text_xy <- locator(n = 1)
          label_text_x <- label_text_xy[[1]]
          label_text_y <- label_text_xy[[2]]
        }
        
        # PLOT LABELS
        .boxed.labels(
          x = label_text_x,
          y = label_text_y,
          labels = label_text,
          border = FALSE,
          font = label_text_font,
          cex = label_text_size,
          col = adjustcolor(label_text_col, label_text_col_alpha),
          bg = label_fill,
          alpha.bg = label_fill_alpha
        )
        
      }
      
      # RETURN LABEL CO-ORDINATES
      return(c(label_text_x, label_text_y))
      
    },
    label_text,
    label_text_x,
    label_text_y,
    label_text_font,
    label_text_size,
    label_text_col,
    label_text_col_alpha,
    label_fill,
    label_fill_alpha,
    SIMPLIFY = FALSE
  )
  
  # RETURN LABEL CO-ORDINATES --------------------------------------------------
  
  # LABEL CO-ORDINATES IN MATRIX
  label_text_xy <- do.call("rbind", label_text_xy)
  colnames(label_text_xy) <- c("x","y")
  label_text_xy <- as.matrix(label_text_xy)
  
  # RETURN LABEL CO-ORDINATES
  invisible(label_text_xy)
  
}
