## CYTO_SPILLOVER_EDIT  -------------------------------------------------------

#' Interactively Edit Spillover Matrices in Real-Time
#'
#' \code{cyto_spillover_edit} provides an interactive shiny interface for
#' editing fluorescent spillover matrices.
#'
#' \code{cyto_spillover_edit} takes on either a
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' compensation controls and/or samples. It is recommended that samples be
#' pre-gated based on FSC and SSC parameters to obtain a homogeneous population
#' for calculation of fluorescent spillover. The compensation controls should
#' also be transformed prior to using \code{cyto_spillover_edit}.
#'
#' Users begin by selecting the unstained control and a stained control from
#' dropdown menus of sample names. \code{cyto_spillover_edit} leverages
#' \code{cyto_plot} to plot the stained sample and overlay the unstained control
#' in black. Users should then select the channel associated with the selected
#' control on the \code{x axis} and go through all other channels on the \code{y
#' axis}.
#'
#' The displayed spillover matrix is extracted directly from the
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spillover argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time.
#'
#' To aid in selection of an appropriate spillover value, the median fluorescent
#' intensity of the unstained control is indicated by a red line and median
#' fluorescent intensity of the stained control is tracked with a purple line.
#' These features can be turned off by de-selecting the check boxes. Changes to
#' the spillover matrix are automatically saved to a csv file called
#' \code{"date-Spillover-Matrix.csv"} in the case where the \code{spillover} is
#' not specified or to the same name as the specified \code{spillover}.
#'
#' @param x an object of class \code{flowSet} or \code{GatingSet}.
#' @param parent name of the parent population to plot when a \code{GatingSet}
#'   object is supplied.
#' @param channel_match name of csv file matching the name of each sample to a
#'   fluorescent channel. The \code{channel_match} file must contain the columns
#'   "name" and "channel". The \code{channel_match} file not required to use
#'   \code{cyto_spillover_edit} but is used internally to automatically select
#'   channels associated with the selcted samples.
#' @param spillover name of a square spillover matrix csv file or spillover
#'   matrix to edit. Setting \code{spillover} to NULL (the default) will result
#'   in extraction of the spillover matrix directly from the supplied samples
#'   (i.e. edit the spillover matrix constructed on the cytometer). Similarly,
#'   if the supplied file name does not exist the spillover matrix from the
#'   first sample will be used and the edited matrix will be saved to the
#'   specified file.
#' @param axes_trans an object of class \code{transformerList} containing
#'   transformers to used to transform the fluorescent channels of the samples
#'   for visualisation.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges.
#' @param display numeric passed to \code{cyto_plot} to control the number of
#'   events to be displayed in the plots, set to 2000 events by default.
#' @param point_size integer passed to \code{cyto_plot} to control the size of
#'   the points in all plots, set to 3 by default.
#' @param axes_text_size numeric pasedd to \code{cyto_plot} to control the size
#'   of axes text, set to 1.7 by default.
#' @param axes_label_text_size numeric passed to \code{cyto_plot} to control the
#'   text size of axes labels, set to 2 by default.
#' @param title_text_size numeric passed to \code{cyto_plot} to control the text
#'   size of titles above each plot, set to 2 by default.
#' @param header_text_size numeric passed to \code{cyto_plot_compensation} to
#'   control size of the header text, set to 1.5 by default.
#' @param viewer logical indicating whether the spillover matrix editor should
#'   be launched in the RStudio viewer pane, set to FALSE by default.
#' @param ... additional arguments passed to \code{cyto_plot}.
#'
#' @return edited spillover matrix and save to designated \code{spillover} csv
#'   file. Saved filename defaults to \code{date-Spillover-Matrix.csv} is not
#'   specified.
#'
#' @importFrom flowCore compensate fsApply exprs Subset each_col
#' @importFrom utils read.csv write.csv
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarPanel selectInput
#'   checkboxInput actionButton mainPanel plotOutput reactiveValues observe
#'   eventReactive renderPlot tabsetPanel tabPanel sidebarLayout fluidRow
#'   updateSelectInput onStop stopApp runApp updateCheckboxInput paneViewer
#' @importFrom rhandsontable rhandsontable rHandsontableOutput hot_to_r
#'   renderRHandsontable hot_cols hot_rows
#' @importFrom shinythemes shinytheme
#' @importFrom magrittr %>%
#' @importFrom methods as is
#' @importFrom stats median
#' @importFrom graphics lines layout
#' @importFrom tools file_ext
#' @importFrom flowWorkspace gs_cyto_data
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_spillover_compute}}
#' @seealso \code{\link{cyto_plot_compensation}}
#' @seealso \code{\link{cyto_plot}}
#'
#' @name cyto_spillover_edit
NULL

#' @noRd
#' @export
cyto_spillover_edit <- function(x, ...) {
  UseMethod("cyto_spillover_edit")
}

#' @rdname cyto_spillover_edit
#' @export
cyto_spillover_edit.GatingSet <- function(x,
                                          parent = NULL,
                                          channel_match = NULL,
                                          spillover = NULL,
                                          axes_trans = NULL,
                                          axes_limits = "machine",
                                          display = 2000,
                                          point_size = 3,
                                          axes_text_size = 1.7,
                                          axes_label_text_size = 2,
                                          title_text_size = 2,
                                          header_text_size = 1.5,
                                          viewer = FALSE,
                                          ...) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # COPY
  gs <- cyto_copy(x)

  # SAMPLE NAMES
  nms <- cyto_names(gs)

  # EXPERIMENT DETAILS
  pd <- cyto_details(gs)

  # CHANNELS
  channels <- cyto_fluor_channels(gs)

  # TRANSFORMATIONS
  axes_trans <- cyto_transformer_extract(gs)

  # COMPENSATION
  comp <- cyto_spillover_extract(gs)

  # GS SHOULD NOT BE COMPENSATED (TRANSFORMED LATER IF REQUIRED)
  if(!is.null(comp)){
    cs <- cyto_extract(gs, "root")
    # UNTRANSFORMED - REMOVE COMPENSATION
    if(.all_na(axes_trans)){
      cyto_compensate(cs, 
                      spillover = comp, 
                      remove = TRUE)
    # TRANSFORMED - REMOVE COMPENSATION
    }else{
      # INVERSE TRANSFORM
      cyto_transform(cs,
                     trans = axes_trans,
                     inverse = TRUE,
                     plot = FALSE)
      # DECOMPENSATE
      cyto_compensate(cs,
                      spillover = comp,
                      remove = TRUE)
      # TRANSFORM
      cyto_transform(cs,
                     trans = axes_trans,
                     plot = FALSE)
    }
  }
  
  # INVERSE TRANSFORMATIONS - GS_LINEAR
  if (.all_na(axes_trans)) {
    gs_linear <- cyto_copy(gs)
    # REVERSE COMPENSATION
    if(!is.null(comp)){
      cs_linear <- cyto_extract(gs_linear, "root")
      cyto_compensate(cs_linear, 
                      spillover = comp, 
                      remove = TRUE)
    }
    axes_trans <- cyto_transformer_biex(gs_linear,
      channels = channels,
      plot = FALSE
    )
    gs <- cyto_transform(gs,
      trans = axes_trans,
      plot = FALSE
    )
  } else {
    gs_linear <- cyto_copy(gs)
    # INVERSE TRANSFORMATIONS
    if(any(channels %in% names(axes_trans))){
      cs_linear <- cyto_extract(gs_linear, "root")
      trans <- axes_trans[match_ind(channels, names(axes_trans))]
      trans <- cyto_transformer_combine(trans)
      cyto_transform(cs_linear,
        trans = trans,
        inverse = TRUE,
        plot = FALSE
      )
    }
    # REVERSE COMPENSATION
    if(!is.null(comp)){
      cs_linear <- cyto_extract(gs_linear, "root")
      cyto_compensate(cs_linear,
                      spillover = comp,
                      remove = TRUE)
    }
  }

  # PREPARE CHANNEL_MATCH ----------------------------------------------------

  # PREPARE CHANNEL_MATCH VARIABLE (MARKERS TO CHANNELS)
  if (any(grepl("channel", colnames(pd), ignore.case = TRUE))) {
    # MARKERS TO CHANNELS
    ind <- which(grepl("channel", colnames(pd), ignore.case = TRUE))
    pd[, "channel"] <- LAPPLY(pd[, ind], function(z) {
      if (!grepl("unstained", z, ignore.case = TRUE)) {
        return(cyto_channels_extract(gs, z))
      } else {
        return(z)
      }
    })
    # CHANNEL_MATCH MISSING IN EXPERIMENT DETAILS
  } else {
    # CHANNEL_MATCH SUPPLIED
    if (!is.null(channel_match)) {
      # CHANNEL_MATCH OBJECT
      if (is(channel_match, "data.frame") |
        is(channel_match, "matrix") |
        is(channel_match, "tibble")) {
        if (!all(c("name", "channel") %in% colnames(channel_match))) {
          stop("channel_match must contain columns 'name' and 'channel'.")
        }
      } else {
        # File extension
        channel_match <- file_ext_append(channel_match, ".csv")
        # File exists
        if(file.exists(channel_match)){
          channel_match <- read.csv(channel_match,
                                    header = TRUE,
                                    row.names = 1,
                                    stringsAsFactors = FALSE)
        # file does not exist
        }else{
          stop(paste(channel_match, 
                     "does not exist or lacks required permissions."))
        }
      }
      # Names of samples don't match any listed in channel_match
      if (!any(nms %in% rownames(channel_match))) {
        # Add NA channel selections to pd
        pd$channel <- rep(NA, nrow(pd))
        # channel_match contains selections for supplied samples
      } else {
        # Add NA channel selections to pd
        pd$channel <- rep(NA, nrow(pd))
        # Replace with channel selections from channel_match
        lapply(nms, function(z) {
          # Update channel selction if covered by channel_match file
          if (z %in% rownames(channel_match)) {
            chn <- channel_match$channel[match(z, rownames(channel_match))]
            pd[match(z, pd$name), "channel"] <<- chn
          }
        })
      }
      # CHANNELS NOT SPECIFIED
    } else {
      pd$channel <- rep(NA, nrow(pd))
    }
  }

  # PREPARE PARENTS ------------------------------------------------------------

  # PARENT MISSING
  if (is.null(parent)) {
    if (!"parent" %in% colnames(pd)) {
      if (!is.null(channel_match)) {
        if ("parent" %in% colnames(channel_match)) {
          parent <- channel_match[, "parent"]
          pd[, "parent"] <- channel_match[, "parent"]
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

  # PREPARE SPILLOVER MATRIX ---------------------------------------------------

  # Spillover matrix supplied
  if (!is.null(spillover)) {
    # spillover is a data.frame or matrix or tibble
    if (is(spillover, "matrix") |
      is(spillover, "data.frame")) {
      # spill should be a matrix
      spill <- spillover
      spill <- as.matrix(spill)
      # Save edited spillover matrix to date-Spillover-Matrix.csv
      spillover <- paste0(
        format(Sys.Date(), "%d%m%y"),
        "-", "Spillover-Matrix.csv"
      )
      # spillover is the name of csv file
    } else {
      # File extension missing
      spillover <- file_ext_append(spillover, ".csv")
      # File does not exist
      if(!file.exists(spillover)){
        # Use first spillover matrix
        spill <- cyto_spillover_extract(
          cyto_extract(gs)
          )[[1]]
      # Use matrix from file
      }else{
        spill <- read.csv(spillover,
                          header = TRUE,
                          row.names = 1)
        spill <- as.matrix(spill)
      }
    }
    # No spillover matrix supplied
  } else {
    spillover <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Spillover-Matrix.csv"
    )
    # Use spillover matrix attached to first sample
    spill <- cyto_spillover_extract(
      cyto_extract(gs)
      )[[1]]
  }
  colnames(spill) <- channels
  rownames(spill) <- channels

  # PREPARE DEFAULT SELECTIONS -------------------------------------------------

  # UNSTAINED SAMPLE
  if (any(grepl("unstained", pd$channel, ignore.case = TRUE))) {
    unstained_initial <- pd$name[grepl("unstained",
      pd$channel,
      ignore.case = TRUE
    )][1]
  } else {
    unstained_initial <- NA
  }

  # SAMPLE
  if (!.all_na(pd$channel)) {
    # UNSTAINED INCLUDED
    if (any(grepl("unstained", pd$channel))) {
      editor_initial_sample <- pd$name[!grepl("unstained",
        pd$channel,
        ignore.case = TRUE
      )][1]
      # NO UNSTAINED - USE FIRST SAMPLE
    } else {
      editor_initial_sample <- pd$name[1]
    }
  } else {
    editor_initial_sample <- pd$name[1]
  }

  # X CHANNEL
  if (!.all_na(pd$channel[pd$name == editor_initial_sample])) {
    editor_initial_xchannel <- pd$channel[pd$name == editor_initial_sample]
  } else{
    editor_initial_xchannel <- channels[1]
  }

  # SHINY SPILLOVER MATRIX EDITOR ----------------------------------------------

  # APPLICATION
  app <- shinyApp(
    ui <- fluidPage(
      theme = shinytheme("yeti"),
      titlePanel("CytoExploreR Spillover Matrix Editor"),
      tabsetPanel(
        tabPanel("Editor",
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              selectInput(
                inputId = "editor_unstained",
                label = "Select Unstained Control:",
                choices = c("No Unstained Control", nms),
                selected = unstained_initial
              ),
              selectInput(
                inputId = "editor_sample",
                label = "Select Sample:",
                choices = nms,
                selected = editor_initial_sample
              ),
              selectInput(
                inputId = "editor_xchannel",
                label = "X Axis:",
                choices = channels,
                selected = editor_initial_xchannel
              ),
              selectInput(
                inputId = "editor_ychannel",
                label = "Y Axis:",
                choices = channels,
                selected = channels[2]
              ),
              checkboxInput(
                inputId = "editor_unstained_overlay",
                label = "Overlay Unstained Control",
                value = TRUE
              ),
              checkboxInput(
                inputId = "editor_unstained_median",
                label = "Unstained Control Median",
                value = TRUE
              ),
              checkboxInput(
                inputId = "editor_median_tracker",
                label = "Median Tracker",
                value = TRUE
              ),
              actionButton("editor_save_button", "Save")
            ),
            mainPanel(
              rHandsontableOutput("spillover_matrix",
                height = "300px"
              ),
              plotOutput("editor_plots",
                height = "400px",
                width = "80%"
              )
            )
          )
        ),
        tabPanel("Plots",
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              selectInput(
                inputId = "plots_unstained",
                label = "Select Unstained Control:",
                choices = c("No Unstained Control", nms),
                selected = unstained_initial
              ),
              selectInput(
                inputId = "plots_sample",
                label = "Select Sample:",
                choices = nms,
                selected = editor_initial_sample
              ),
              selectInput(
                inputId = "plots_xchannel",
                label = "Channel",
                choices = channels,
                selected = editor_initial_xchannel
              ),
              checkboxInput(
                inputId = "plots_unstained_overlay",
                label = "Overlay Unstained Control",
                value = TRUE
              ),
              checkboxInput(
                inputId = "plots_uncompensated_underlay",
                label = "Underlay Uncompensated Control",
                value = TRUE
              ),
              checkboxInput(
                inputId = "plots_compensated_overlay",
                label = "Overlay Compensated Control",
                value = TRUE
              )
            ),
            mainPanel(
              fluidRow(
                plotOutput("plots", height = "700px", width = "100%")
              )
            )
          )
        )
      )
    ),
    # Shiny application server
    server <- function(input, output, session) {
      values <- reactiveValues()

      observe({
        if (!is.null(input$spillover_matrix)) {
          spill <- hot_to_r(input$spillover_matrix)
          rownames(spill) <- colnames(spill)
          values$spill <- spill
        } else {
          values$spill <- spill * 100
        }
      })

      output$spillover_matrix <- renderRHandsontable({
        rhandsontable(values$spill,
          rowHeaderWidth = 105,
          readOnly = FALSE
        ) %>%
          hot_cols(
            type = "numeric",
            colWidths = 105,
            format = "0.000",
            halign = "htCenter",
            renderer = "
                    function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);
                    if(value < 0 ){
                    td.style.background = 'lightblue';
                    } else if (value == 0 ){
                    td.style.background = 'white';
                    } else if (value > 0 & value <= 10) {
                    td.style.background = 'lightgreen';
                    } else if (value > 10 & value <= 25){
                    td.style.background = 'yellow';
                    } else if (value > 25 & value <= 50){
                    td.style.background = 'orange';
                    } else if (value > 50 & value < 100){
                    td.style.background = 'red';
                    } else if (value == 100){
                    td.style.background = 'darkgrey';
                    } else if (value > 100){
                    td.style.background = 'violet';
                    }
                    }"
          ) %>%
          hot_rows(rowHeights = 20)
      })

      # Update sample selection in plots tab & x channel selection (match)
      observe({
        gh <- input$editor_sample
        # Update sample selection in Plots tab
        updateSelectInput(session,
          "plots_sample",
          selected = gh
        )
        # Use channel_match to set xchannel
        xchan <- pd$channel[match(gh, pd$name)]
        # Only update x channel if not NA
        if (!.all_na(xchan)) {
          updateSelectInput(session,
            "editor_xchannel",
            selected = xchan
          )
          updateSelectInput(session,
            "plots_xchannel",
            selected = xchan
          )
        }
      })

      # Re-apply compensation and transformations
      gs_comp <- eventReactive(values$spill, {
        # COPY
        gs_linear_copy <- cyto_copy(gs_linear)
        # Data must be LINEAR at this step
        gs_linear_comp <- compensate(gs_linear_copy, values$spill / 100)
        # Get transformed data
        trans <- axes_trans[match_ind(channels, names(axes_trans))]
        trans <- cyto_transformer_combine(trans)
        gs_trans <- cyto_transform(gs_linear_comp,
          trans = trans,
          plot = FALSE
        )
        return(gs_trans)
      })

      # Turn off unstained aspects if no unstained control is supplied
      observe({
        unst <- input$editor_unstained
        # Turn off unstained components if unstained is NA
        if (unst == "No Unstained Control") {
          updateCheckboxInput(session,
            "editor_unstained_overlay",
            value = FALSE
          )
          updateCheckboxInput(session,
            "editor_unstained_median",
            value = FALSE
          )
          # Turn on unstained components
        } else {
          updateCheckboxInput(session,
            "editor_unstained_overlay",
            value = TRUE
          )
          updateCheckboxInput(session,
            "editor_unstained_median",
            value = TRUE
          )
        }
        # Update Plots tab unstained control based on editor selection
        updateSelectInput(session,
          "plots_unstained",
          selected = unst
        )
      })

      # Update channel selection on plots tab based on editor selection
      observe({
        xchan <- input$editor_xchannel
        updateSelectInput(session,
          "plots_xchannel",
          selected = xchan
        )
      })

      # Update sample & channel selection in editor based on Plots selection
      observe({
        smp <- input$plots_sample
        updateSelectInput(session,
          "editor_sample",
          selected = smp
        )
        xchan <- input$plots_xchannel
        updateSelectInput(session,
          "editor_xchannel",
          selected = xchan
        )
      })

      # Update overlay unstained check box based on unstained selection (Plots)
      observe({
        unst <- input$plots_unstained
        # Remove unstained overlay
        if (unst == "No Unstained Control") {
          updateCheckboxInput(session,
            "plots_unstained_overlay",
            value = FALSE
          )
          # Turn on unstained overlay
        } else {
          updateCheckboxInput(session,
            "plots_unstained_overlay",
            value = TRUE
          )
        }
        # Update unstained control selection in editor tab
        updateSelectInput(session,
          "editor_unstained",
          selected = unst
        )
      })

      output$editor_plots <- renderPlot({

        # Set up plot layout
        layout(matrix(c(1, 1, 2, 2, 1, 1, 3, 3),
          byrow = TRUE,
          ncol = 4
        ))

        # No unstained control - no unstained overlay or unstained median
        if (input$editor_unstained == "No Unstained Control") {

          # Plots
          cyto_plot(gs_comp()[[input$editor_sample]],
            channels = c(
              input$editor_xchannel,
              input$editor_ychannel
            ),
            parent = pd[pd$name == input$editor_sample, "parent"],
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = input$editor_sample,
            point_size = point_size,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size, ...
          )

          # Median tracker
          if (input$editor_median_tracker == TRUE) {

            # MEDIAN TRACKER
            pop <- cyto_extract(
              gs_comp()[[input$editor_sample]],
              pd[pd$name == input$editor_sample, "parent"]
            )
            .cyto_median_tracker(
              pop,
              c(
                input$editor_xchannel,
                input$editor_ychannel
              )
            )
          }

          # Density distribution in associated channel
          cyto_plot(gs_comp()[[input$editor_sample]],
            channels = input$editor_xchannel,
            parent = pd[pd$name == input$editor_sample, "parent"],
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = NA,
            ylab = "Density",
            density_fill = "white",
            density_line_col = "blue",
            density_line_width = 2,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size
          )

          # Density distribution in other channel
          cyto_plot(gs_comp()[[input$editor_sample]],
            channels = input$editor_ychannel,
            parent = pd[pd$name == input$editor_sample, "parent"],
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = NA,
            ylab = "Density",
            density_fill = "white",
            density_line_col = "blue",
            density_line_width = 2,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size
          )

          # Add MedFI label in other channel
          cyto_plot_label(gs_comp()[[input$editor_sample]],
            channels = input$editor_ychannel,
            parent = pd[pd$name == input$editor_sample, "parent"],
            display = display,
            label_text = "MedFI",
            label_stat = "median",
            label_text_x = 0.75 * par("usr")[2],
            label_text_y = 30,
            label_text_col = "red",
            label_text_size = 1.1
          )

          # Unstained control supplied
        } else if (input$editor_unstained != "No Unstained Control") {

          # NO UNSTAINED OVERLAY
          if (input$editor_unstained_overlay == FALSE) {

            # Plot
            cyto_plot(gs_comp()[[input$editor_sample]],
              channels = c(
                input$editor_xchannel,
                input$editor_ychannel
              ),
              parent = pd[pd$name == input$editor_sample, "parent"],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = input$editor_sample,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size, ...
            )

            # Median line through unstained control
            if (input$editor_unstained_median == TRUE) {
              # PULL OUT NIL - RAW DATA
              NIL <- cyto_extract(gs_comp()[[input$editor_unstained]],
                pd[pd$name == input$editor_sample, "parent"],
                raw = TRUE
              )[[1]]
              # COMPUTE MEDFI IN ALL CHANNELS
              medFI <- lapply(channels, function(z) {
                median(NIL[, z])
              })
              medFI <- do.call("cbind", medFI)
              colnames(medFI) <- channels
              cutoff <- medFI[1, input$editor_ychannel]
              abline(h = cutoff, col = "red", lwd = 2)
            }

            # Median tracker through stained control
            if (input$editor_median_tracker == TRUE) {

              # MEDIAN TRACKER
              pop <- cyto_extract(
                gs_comp()[[input$editor_sample]],
                pd[pd$name == input$editor_sample, "parent"]
              )
              .cyto_median_tracker(
                pop,
                c(
                  input$editor_xchannel,
                  input$editor_ychannel
                )
              )
            }

            # Density distribution in associated channel
            cyto_plot(gs_comp()[[input$editor_sample]],
              channels = input$editor_xchannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = NA,
              ylab = "Density",
              density_stack = 0,
              density_fill = "white",
              density_fill_alpha = 0,
              density_line_col = "blue",
              density_line_width = 2,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = 1.25 * title_text_size
            )

            # Density distribution in other channel
            cyto_plot(gs_comp()[[input$editor_sample]],
              channels = input$editor_ychannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = NA,
              ylab = "Density",
              density_stack = 0,
              density_fill = "white",
              density_fill_alpha = 0,
              density_line_col = "blue",
              density_line_width = 2,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = 1.25 * title_text_size
            )

            # Add label for Stained median
            cyto_plot_label(gs_comp()[[input$editor_sample]],
              channels = input$editor_ychannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              display = display,
              label_text = "MedFI",
              label_stat = "median",
              label_text_x = 0.75 * par("usr")[2],
              label_text_y = 30,
              label_text_col = "blue",
              label_text_size = 1.1
            )

            # UNSTAINED OVERLAY
          } else if (input$editor_unstained_overlay == TRUE) {

            # Plot
            cyto_plot(gs_comp()[[input$editor_sample]],
              channels = c(
                input$editor_xchannel,
                input$editor_ychannel
              ),
              parent = pd[pd$name == input$editor_sample, "parent"],
              overlay = cyto_extract(
                gs_comp()[[input$editor_unstained]],
                pd[pd$name == input$editor_sample, "parent"]
              ),
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = input$editor_sample,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size, ...
            )

            # Median line through unstained control
            if (input$editor_unstained_median == TRUE) {
              # PULL OUT NIL - RAW DATA
              NIL <- cyto_extract(gs_comp()[[input$editor_unstained]],
                pd[pd$name == input$editor_sample, "parent"],
                raw = TRUE
              )[[1]]
              # COMPUTE MEDFI IN ALL CHANNELS
              medFI <- lapply(channels, function(z) {
                median(NIL[, z])
              })
              medFI <- do.call("cbind", medFI)
              colnames(medFI) <- channels
              cutoff <- medFI[1, input$editor_ychannel]
              abline(h = cutoff, col = "red", lwd = 2)
            }

            # MEDIAN TRACKER
            if (input$editor_median_tracker == TRUE) {

              # MEDIAN TRACKER
              pop <- cyto_extract(
                gs_comp()[[input$editor_sample]],
                pd[pd$name == input$editor_sample, "parent"]
              )
              .cyto_median_tracker(
                pop,
                c(
                  input$editor_xchannel,
                  input$editor_ychannel
                )
              )
            }

            # Density distribution in associated channel - unstained overlay
            cyto_plot(gs_comp()[[input$editor_unstained]],
              channels = input$editor_xchannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              overlay = cyto_extract(
                gs_comp()[[input$editor_sample]],
                pd[pd$name == input$editor_sample, "parent"]
              ),
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = NA,
              ylab = "Density",
              density_stack = 0,
              density_fill = c("grey70", "white"),
              density_fill_alpha = c(1, 0),
              density_line_col = c("black", "blue"),
              density_line_width = 2,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = 1.25 * title_text_size
            )

            # Density distribution in other channel - unstained overlay
            cyto_plot(gs_comp()[[input$editor_unstained]],
              channels = input$editor_ychannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              overlay = cyto_extract(
                gs_comp()[[input$editor_sample]],
                pd[pd$name == input$editor_sample, "parent"]
              ),
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = NA,
              ylab = "Density",
              density_stack = 0,
              density_fill = c("grey70", "white"),
              density_fill_alpha = c(1, 0),
              density_line_col = c("black", "blue"),
              density_line_width = 2,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = 1.25 * title_text_size
            )

            # Add label for unstained median
            cyto_plot_label(gs_comp()[[input$editor_unstained]],
              channels = input$editor_ychannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              display = display,
              label_text = "MedFI",
              label_stat = "median",
              label_text_x = 0.75 * par("usr")[2],
              label_text_y = 80,
              label_text_col = "grey40",
              label_text_size = 1.1
            )

            # Add label for Stained median
            cyto_plot_label(gs_comp()[[input$editor_sample]],
              channels = input$editor_ychannel,
              parent = pd[pd$name == input$editor_sample, "parent"],
              display = display,
              label_text = "MedFI",
              label_stat = "median",
              label_text_x = 0.75 * par("usr")[2],
              label_text_y = 30,
              label_text_col = "blue",
              label_text_size = 1.1
            )
          }
        }
      })

      # Plots tab uses cyto_plot_compensation
      output$plots <- renderPlot({

        # Compensation plots - fill colours are reversed internally
        if (input$plots_uncompensated_underlay == TRUE) {
          if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(gs[[input$plots_sample]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent", ],
              axes_trans = axes_trans,
              overlay = cyto_extract(
                gs_comp()[[input$plots_sample]],
                pd[pd$name == input$plots_sample, "parent", ]
              ),
              axes_limits = axes_limits,
              display = display,
              point_col = c("blue", "red"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("blue", "red")
            )
          } else if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(gs[[input$plots_sample]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent"],
              axes_trans = axes_trans,
              overlay = list(
                cyto_extract(
                  gs_comp()[[input$plots_sample]],
                  pd[pd$name == input$plots_sample, "parent"]
                ),
                cyto_extract(
                  gs_comp()[[input$plots_unstained]],
                  pd[pd$name == input$plots_sample, "parent"]
                )
              ),
              axes_limits = axes_limits,
              display = display,
              point_col = c("blue", "red", "black"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("grey", "red", "blue")
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(gs[[input$plots_sample]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent"],
              axes_trans = axes_trans,
              overlay = cyto_extract(
                gs_comp()[[input$plots_unstained]],
                pd[pd$name == input$plots_sample, "parent"]
              ),
              axes_limits = axes_limits,
              display = display,
              point_col = c("blue", "black"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("grey", "blue")
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(gs[[input$plots_sample]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent"],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              point_col = "blue",
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = "blue"
            )
          }
        } else if (input$plots_uncompensated_underlay == FALSE) {
          if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(gs_comp()[[input$plots_sample]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent"],
              axes_trans = axes_trans,
              overlay = cyto_extract(
                gs_comp()[[input$plots_unstained]],
                pd[pd$name == input$plots_sample, "parent"]
              ),
              axes_limits = axes_limits,
              display = display,
              point_col = c("red", "black"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("grey", "red")
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(gs_comp()[[input$plots_unstained]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent"],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              point_col = "black",
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = "grey"
            )
          } else if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(gs_comp()[[input$plots_sample]],
              channel = input$plots_xchannel,
              parent = pd[pd$name == input$plots_sample, "parent"],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              point_col = "red",
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = "red"
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == FALSE) {
            # No data to plot
          }
        }
      })

      # Save edited spillover matrix to csv file
      observe({
        input$saveBtn
        class(values$spill) <- "numeric"
        spill.mat <- values$spill / 100
        write.csv(spill.mat, spillover)
      })

      # Return edited matrix on application close
      onStop(function() {
        spill.mat <- read.csv(spillover, 
                              header = TRUE, 
                              row.names = 1)
        colnames(spill.mat) <- rownames(spill.mat)
        stopApp(spill.mat)
      })
    }
  )

  # Run the shiny application
  if(viewer){
    sp <- runApp(app, 
                 launch.browser = paneViewer(),
                 quiet = TRUE)
  }else{
    sp <- runApp(app,
                 quiet = TRUE)
  }

  # RETURN UPDATED SPILLOVER MATRIX
  return(sp)
}

#' @rdname cyto_spillover_edit
#' @export
cyto_spillover_edit.flowSet <- function(x,
                                        channel_match = NULL,
                                        spillover = NULL,
                                        axes_trans = NULL,
                                        axes_limits = "machine",
                                        display = 2000,
                                        point_size = 3,
                                        axes_text_size = 1.7,
                                        axes_label_text_size = 2,
                                        title_text_size = 2,
                                        header_text_size = 1.5,
                                        viewer = FALSE,
                                        ...) {

  # Assign x to fs
  fs <- x

  # Copy fs (prevent transforming underlying data)
  fs <- cyto_copy(fs)

  # Extract sample names
  nms <- cyto_names(fs)

  # Extract fluorescent channels
  channels <- cyto_fluor_channels(fs)

  # Extract cyto details to pd
  pd <- cyto_details(fs)

  # Inverse transformations
  if (.all_na(axes_trans)) {
    fs_linear <- cyto_copy(fs)
    axes_trans <- cyto_transformer_biex(fs,
      channels = channels,
      plot = FALSE
    )
    fs <- cyto_transform(fs,
      trans = axes_trans,
      plot = FALSE
    )
  } else {
    trans <- axes_trans[match_ind(channels, names(axes_trans))]
    trans <- cyto_transformer_combine(trans)
    fs_linear <- cyto_transform(cyto_copy(fs),
      trans = trans,
      inverse = TRUE,
      plot = FALSE
    )
  }

  # Channel match file supplied
  if (!is.null(channel_match)) {
    # channel_match is a data.frame or matrix or tibble
    if (is(channel_match, "data.frame") |
      is(channel_match, "matrix")) {
      # channel_match must contain "name" and "channel" columns
      if (!any(grepl("name", colnames(channel_match), ignore.case = TRUE)) |
        !any(grepl("channel", colnames(channel_match), ignore.case = TRUE))) {
        stop("'channel_match' must contain the columns 'name' and 'channel'.")
      }
      # channel_match is the name of a csv file
    } else {
      # channel_match must contain csv file extension
      channel_match <- file_ext_append(channel_match, ".csv")
      # file exists
      if(file.exists(channel_match)){
        channel_match <- read.csv(channel_match,
                                  header = TRUE,
                                  row.names = 1,
                                  stringsAsFactors = FALSE)
      # file does not exist
      }else{
        stop(paste(channel_match, 
                   "does not exist or lacks required permissions."))
      }
    }

    # Names of samples don't match any listed in channel_match
    if (!any(nms %in% rownames(channel_match))) {
      # Add NA channel selections to pd
      pd$channel <- rep(NA, nrow(pd))
      # channel_match contains selections for supplied samples
    } else {
      # Add NA channel selections to pd
      pd$channel <- rep(NA, nrow(pd))
      # Replace with channel selections from channel_match
      lapply(nms, function(z) {
        # Update channel selction if covered by channel_match file
        if (z %in% rownames(channel_match)) {
          chn <- channel_match$channel[match(z, rownames(channel_match))]
          pd[match(z, pd$name), "channel"] <<- chn
        }
      })
    }

    # No channel match file supplied
  } else if (is.null(channel_match)) {
    # Add NA channel selections to pd
    pd$channel <- rep(NA, nrow(pd))
  }

  # Spillover matrix supplied
  if (!is.null(spillover)) {
    # spillover is a data.frame or matrix or tibble
    if (is(spillover, "matrix") |
      is(spillover, "data.frame")) {
      # spill should be a matrix
      spill <- spillover
      spill <- as.matrix(spill)
      # Save edited spillover matrix to date-Spillover-Matrix.csv
      spillover <- paste0(
        format(Sys.Date(), "%d%m%y"),
        "-", "Spillover-Matrix.csv"
      )
      # spillover is the name of csv file
    } else {
      # File extension missing
      spillover <- file_ext_append(spillover, ".csv")
      # File does not exist
      if(!file.exists(spillover)){
        # Use first spillover matrix
        spill <- cyto_spillover_extract(fs)[[1]]
        # Use matrix from file
      }else{
        spill <- read.csv(spillover,
                          header = TRUE,
                          row.names = 1)
        spill <- as.matrix(spill)
      }
    }
    # No spillover matrix supplied
  } else {
    spillover <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Spillover-Matrix.csv"
    )
    # Use spillover matrix attached to first sample
    spill <- cyto_spillover_extract(fs)[[1]]
  }
  colnames(spill) <- channels
  rownames(spill) <- channels

  # Unstained supplied?
  if (any(grepl("Unstained", pd$channel, ignore.case = TRUE))) {
    unstained_initial <- pd$name[grepl("Unstained",
      pd$channel,
      ignore.case = TRUE
    )][1]
  } else {
    unstained_initial <- NA
  }

  # Selected sample - unstained control supplied
  if (!.all_na(pd$channel)) {
    # Unstained control included
    if (any(grepl("Unstained", pd$channel))) {
      editor_initial_sample <- pd$name[!grepl("Unstained",
        pd$channel,
        ignore.case = TRUE
      )][1]
      # No unstained control - use first sample
    } else {
      editor_initial_sample <- pd$name[1]
    }
  } else {
    editor_initial_sample <- pd$name[1]
  }

  # X channel selection
  if (!.all_na(pd$channel[pd$name == editor_initial_sample])) {
    editor_initial_xchannel <- pd$channel[pd$name == editor_initial_sample]
  } else {
    editor_initial_xchannel <- channels[1]
  }

  # Shiny application
  app <- shinyApp(
    ui <- fluidPage(
      theme = shinytheme("yeti"),
      titlePanel("CytoExploreR Spillover Matrix Editor"),
      tabsetPanel(
        tabPanel("Editor",
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              selectInput(
                inputId = "editor_unstained",
                label = "Select Unstained Control:",
                choices = c("No Unstained Control", nms),
                selected = unstained_initial
              ),
              selectInput(
                inputId = "editor_sample",
                label = "Select Sample:",
                choices = nms,
                selected = editor_initial_sample
              ),
              selectInput(
                inputId = "editor_xchannel",
                label = "X Axis:",
                choices = channels,
                selected = editor_initial_xchannel
              ),
              selectInput(
                inputId = "editor_ychannel",
                label = "Y Axis:",
                choices = channels,
                selected = channels[2]
              ),
              checkboxInput(
                inputId = "editor_unstained_overlay",
                label = "Overlay Unstained Control",
                value = TRUE
              ),
              checkboxInput(
                inputId = "editor_unstained_median",
                label = "Unstained Control Median",
                value = TRUE
              ),
              checkboxInput(
                inputId = "editor_median_tracker",
                label = "Median Tracker",
                value = TRUE
              ),
              actionButton("editor_save_button", "Save")
            ),
            mainPanel(
              rHandsontableOutput("spillover_matrix", height = "300px"),
              plotOutput("editor_plots", height = "400px", width = "80%")
            )
          )
        ),
        tabPanel("Plots",
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              selectInput(
                inputId = "plots_unstained",
                label = "Select Unstained Control:",
                choices = c("No Unstained Control", nms),
                selected = unstained_initial
              ),
              selectInput(
                inputId = "plots_sample",
                label = "Select Sample:",
                choices = nms,
                selected = editor_initial_sample
              ),
              selectInput(
                inputId = "plots_xchannel",
                label = "Channel",
                choices = channels,
                selected = editor_initial_xchannel
              ),
              checkboxInput(
                inputId = "plots_unstained_overlay",
                label = "Overlay Unstained Control",
                value = TRUE
              ),
              checkboxInput(
                inputId = "plots_uncompensated_underlay",
                label = "Underlay Uncompensated Control",
                value = TRUE
              ),
              checkboxInput(
                inputId = "plots_compensated_overlay",
                label = "Overlay Compensated Control",
                value = TRUE
              )
            ),
            mainPanel(
              fluidRow(
                plotOutput("plots", height = "700px", width = "100%")
              )
            )
          )
        )
      )
    ),
    # Shiny application server
    server <- function(input, output, session) {
      values <- reactiveValues()

      observe({
        if (!is.null(input$spillover_matrix)) {
          spill <- hot_to_r(input$spillover_matrix)
          rownames(spill) <- colnames(spill)
          values$spill <- spill
        } else {
          values$spill <- spill * 100
        }
      })

      output$spillover_matrix <- renderRHandsontable({
        rhandsontable(values$spill,
          rowHeaderWidth = 105,
          readOnly = FALSE
        ) %>%
          hot_cols(
            type = "numeric",
            colWidths = 105,
            format = "0.000",
            halign = "htCenter",
            renderer = "
                    function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);
                    if(value < 0 ){
                    td.style.background = 'lightblue';
                    } else if (value == 0 ){
                    td.style.background = 'white';
                    } else if (value > 0 & value <= 10) {
                    td.style.background = 'lightgreen';
                    } else if (value > 10 & value <= 25){
                    td.style.background = 'yellow';
                    } else if (value > 25 & value <= 50){
                    td.style.background = 'orange';
                    } else if (value > 50 & value < 100){
                    td.style.background = 'red';
                    } else if (value == 100){
                    td.style.background = 'darkgrey';
                    } else if (value > 100){
                    td.style.background = 'violet';
                    }
                    }"
          ) %>%
          hot_rows(rowHeights = 20)
      })

      # Update sample selection in plots tab & x channel selection (match)
      observe({
        fr <- input$editor_sample
        # Update sample selection in Plots tab
        updateSelectInput(session, "plots_sample", selected = fr)
        # Use channel_match to set xchannel
        xchan <- pd$channel[match(fr, pd$name)]
        # Only update x channel if not NA
        if (!.all_na(xchan)) {
          updateSelectInput(session, "editor_xchannel", selected = xchan)
          updateSelectInput(session, "plots_xchannel", selected = xchan)
        }
      })

      # Re-apply compensation and transformations
      fs.comp <- eventReactive(values$spill, {
        # Make a copy
        fs_copy <- cyto_copy(fs_linear)
        # Data must be LINEAR at this step
        fs_comp <- compensate(fs_copy, values$spill / 100)
        # Get transformed data
        trans <- axes_trans[match_ind(channels, names(axes_trans))]
        trans <- cyto_transformer_combine(trans)
        fs_trans <- cyto_transform(fs_comp,
          trans = trans,
          plot = FALSE
        )
        return(fs_trans)
      })

      # Turn off unstained aspects if no unstained control is supplied
      observe({
        unst <- input$editor_unstained
        # Turn off unstained components if unstained is NA
        if (unst == "No Unstained Control") {
          updateCheckboxInput(session, "editor_unstained_overlay", value = FALSE)
          updateCheckboxInput(session, "editor_unstained_median", value = FALSE)
          # Turn on unstained components
        } else {
          updateCheckboxInput(session, "editor_unstained_overlay", value = TRUE)
          updateCheckboxInput(session, "editor_unstained_median", value = TRUE)
        }
        # Update Plots tab unstained control based on editor selection
        updateSelectInput(session, "plots_unstained", selected = unst)
      })

      # Update channel selection on plots tab based on editor selection
      observe({
        xchan <- input$editor_xchannel
        updateSelectInput(session, "plots_xchannel", selected = xchan)
      })

      # Update sample & channel selection in editor based on Plots selection
      observe({
        smp <- input$plots_sample
        updateSelectInput(session, "editor_sample", selected = smp)
        xchan <- input$plots_xchannel
        updateSelectInput(session, "editor_xchannel", selected = xchan)
      })

      # Update overlay unstained check box based on unstained selection (Plots)
      observe({
        unst <- input$plots_unstained
        # Remove unstained overlay
        if (unst == "No Unstained Control") {
          updateCheckboxInput(session, "plots_unstained_overlay", value = FALSE)
          # Turn on unstained overlay
        } else {
          updateCheckboxInput(session, "plots_unstained_overlay", value = TRUE)
        }
        # Update unstained control selection in editor tab
        updateSelectInput(session, "editor_unstained", selected = unst)
      })

      output$editor_plots <- renderPlot({

        # Set up plot layout
        layout(matrix(c(1, 1, 2, 2, 1, 1, 3, 3), byrow = TRUE, ncol = 4))

        # No unstained control - no unstained overlay or unstained median
        if (input$editor_unstained == "No Unstained Control") {

          # Plots
          cyto_plot(fs.comp()[[input$editor_sample]],
            channels = c(
              input$editor_xchannel,
              input$editor_ychannel
            ),
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = input$editor_sample,
            point_size = point_size,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size, ...
          )

          # Median tracker
          if (input$editor_median_tracker == TRUE) {

            # MEDIAN TRACKER
            .cyto_median_tracker(
              fs.comp()[[input$editor_sample]],
              c(
                input$editor_xchannel,
                input$editor_ychannel
              )
            )
          }

          # Density distribution in associated channel
          cyto_plot(fs.comp()[[input$editor_sample]],
            channels = input$editor_xchannel,
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = NA,
            ylab = "Density",
            density_fill = "white",
            density_line_col = "blue",
            density_line_width = 2,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size
          )

          # Density distribution in other channel
          cyto_plot(fs.comp()[[input$editor_sample]],
            channels = input$editor_ychannel,
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = NA,
            ylab = "Density",
            density_fill = "white",
            density_line_col = "blue",
            density_line_width = 2,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size
          )

          # Add MedFI label in other channel
          cyto_plot_label(fs.comp()[[input$editor_sample]],
            channels = input$editor_ychannel,
            trans = axes_trans,
            display = display,
            label_text = "MedFI",
            label_stat = "median",
            label_text_x = 0.75 * par("usr")[2],
            label_text_y = 30,
            label_text_col = "red",
            label_text_size = 1.1
          )

          # Unstained control supplied
        } else if (input$editor_unstained != "No Unstained Control") {

          # Unstained Overlay
          if (input$editor_unstained_overlay == FALSE) {

            # Plot
            cyto_plot(fs.comp()[[input$editor_sample]],
              channels = c(
                input$editor_xchannel,
                input$editor_ychannel
              ),
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = input$editor_sample,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size, ...
            )

            # Median line through unstained control
            if (input$editor_unstained_median == TRUE) {
              medians <- fsApply(fs.comp(), each_col, "median")[
                input$editor_unstained,
                channels
              ]
              MFI <- data.frame("Channel" = channels, "Median" = medians)

              cutoff <- MFI[match(input$editor_ychannel, MFI$Channel), ]

              abline(h = cutoff[2], col = "red", lwd = 2)
            }

            # Median tracker through stained control
            if (input$editor_median_tracker == TRUE) {

              # MEDIAN TRACKER
              .cyto_median_tracker(
                fs.comp()[[input$editor_sample]],
                c(
                  input$editor_xchannel,
                  input$editor_ychannel
                )
              )
            }

            # Density distribution in associated channel - unstained overlay
            cyto_plot(fs.comp()[[input$editor_sample]],
              channels = input$editor_xchannel,
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = NA,
              ylab = "Density",
              density_stack = 0,
              density_fill = "white",
              density_fill_alpha = 0,
              density_line_col = "blue",
              density_line_width = 2,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = 1.25 * title_text_size
            )

            # Add label for Stained median
            cyto_plot_label(fs.comp()[[input$editor_sample]],
              channels = input$editor_ychannel,
              trans = axes_trans,
              display = display,
              label_text = "MedFI",
              label_stat = "median",
              label_text_x = 0.75 * par("usr")[2],
              label_text_y = 30,
              label_text_col = "blue",
              label_text_size = 1.1
            )

            # No unstained overlay
          } else if (input$editor_unstained_overlay == TRUE) {

            # Plot
            cyto_plot(fs.comp()[[input$editor_sample]],
              channels = c(
                input$editor_xchannel,
                input$editor_ychannel
              ),
              overlay = fs.comp()[[input$editor_unstained]],
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              title = input$editor_sample,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size, ...
            )

            # Median line through unstained control
            if (input$editor_unstained_median == TRUE) {
              medians <- fsApply(
                fs.comp(),
                each_col,
                median
              )[input$editor_unstained, channels]
              MFI <- data.frame("Channel" = channels, "Median" = medians)

              cutoff <- MFI[match(input$editor_ychannel, MFI$Channel), ]

              abline(h = cutoff[2], col = "red", lwd = 3)
            }

            # Median tracker through stained control
            if (input$editor_median_tracker == TRUE) {

              # MEDIAN TRACKER
              .cyto_median_tracker(
                fs.comp()[[input$editor_sample]],
                c(
                  input$editor_xchannel,
                  input$editor_ychannel
                )
              )
            }
          }

          # Density distribution in associated channel - unstained overlay
          cyto_plot(fs.comp()[[input$editor_unstained]],
            channels = input$editor_xchannel,
            overlay = fs.comp()[[input$editor_sample]],
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = NA,
            ylab = "Density",
            density_stack = 0,
            density_fill = c("grey70", "white"),
            density_fill_alpha = c(1, 0),
            density_line_col = c("black", "blue"),
            density_line_width = 2,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = 1.25 * title_text_size
          )

          # Density distribution in other channel - unstained overlay
          cyto_plot(fs.comp()[[input$editor_unstained]],
            channels = input$editor_ychannel,
            overlay = fs.comp()[[input$editor_sample]],
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            display = display,
            title = NA,
            ylab = "Density",
            density_stack = 0,
            density_fill = c("grey70", "white"),
            density_fill_alpha = c(1, 0),
            density_line_col = c("black", "blue"),
            density_line_width = 2,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = 1.25 * title_text_size
          )

          # Add label for unstained median
          cyto_plot_label(fs.comp()[[input$editor_unstained]],
            channels = input$editor_ychannel,
            trans = axes_trans,
            display = display,
            label_text = "MedFI",
            label_stat = "median",
            label_text_x = 0.75 * par("usr")[2],
            label_text_y = 80,
            label_text_col = "grey40",
            label_text_size = 1.1
          )

          # Add label for Stained median
          cyto_plot_label(fs.comp()[[input$editor_sample]],
            channels = input$editor_ychannel,
            trans = axes_trans,
            display = display,
            label_text = "MedFI",
            label_stat = "median",
            label_text_x = 0.75 * par("usr")[2],
            label_text_y = 30,
            label_text_col = "blue",
            label_text_size = 1.1
          )
        }
      })

      # Plots tab uses cyto_plot_compensation
      output$plots <- renderPlot({

        # Compensation plots - fill colours are reversed internally
        if (input$plots_uncompensated_underlay == TRUE) {
          if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              overlay = fs.comp()[[input$plots_sample]],
              axes_limits = axes_limits,
              display = display,
              point_col = c("blue", "red"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("blue", "red")
            )
          } else if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              overlay = list(
                fs.comp()[[input$plots_sample]],
                fs.comp()[[input$plots_unstained]]
              ),
              axes_limits = axes_limits,
              display = display,
              point_col = c("blue", "red", "black"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("grey", "red", "blue")
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              overlay = fs.comp()[[input$plots_unstained]],
              axes_limits = axes_limits,
              display = display,
              point_col = c("blue", "black"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("grey", "blue")
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              point_col = "blue",
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = "blue"
            )
          }
        } else if (input$plots_uncompensated_underlay == FALSE) {
          if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs.comp()[[input$plots_sample]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              overlay = fs.comp()[[input$plots_unstained]],
              axes_limits = axes_limits,
              display = display,
              point_col = c("red", "black"),
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = c("grey", "red")
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs.comp()[[input$plots_unstained]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              point_col = "black",
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = "grey"
            )
          } else if (input$plots_compensated_overlay == TRUE &
            input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs.comp()[[input$plots_sample]],
              channel = input$plots_xchannel,
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              display = display,
              point_col = "red",
              point_alpha = 0.6,
              header = input$plots_sample,
              title = NA,
              point_size = point_size,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              header_text_size = header_text_size,
              density_fill = "red"
            )
          } else if (input$plots_compensated_overlay == FALSE &
            input$plots_unstained_overlay == FALSE) {
            # No data to plot
          }
        }
      })

      # Save edited spillover matrix to csv file
      observe({
        input$saveBtn
        class(values$spill) <- "numeric"
        spill.mat <- values$spill / 100
        write.csv(spill.mat, spillover)
      })

      # Return edited matrix on application cloase
      onStop(function() {
        spill.mat <- read.csv(spillover, 
                              header = TRUE, 
                              row.names = 1)
        colnames(spill.mat) <- rownames(spill.mat)
        stopApp(spill.mat)
      })
    }
  )

  # Run the shiny application
  if(viewer){
    sp <- runApp(app, 
                 launch.browser = paneViewer(),
                 quiet = TRUE)
  }else{
    sp <- runApp(app,
                 quiet = TRUE)
  }

  # Return updated spillover matrix
  return(sp)
}

#' Median Tracker
#' Add median tracker to plot
#' @param x object of class flowFrame
#' @param channels channels used to construct the plot
#' @importFrom stats loess median predict
#' @importFrom graphics par
#' @importFrom graphics lines
#' @noRd
.cyto_median_tracker <- function(x,
                                 channels = NULL) {

  # RAW DATA
  raw_data <- cyto_extract(x, raw = TRUE)[[1]]
  raw_data <- raw_data[, channels]
  # SORT RAW DATA
  raw_data <- raw_data[order(raw_data[, channels[1]]), ]
  raw_data <- as.data.frame(raw_data)
  # SPLIT PLOT RANGE INTO CHUNKS
  xlim <- par("usr")
  # NUMBER OF CHUNKS
  n <- 25
  chunks <- seq(
    from = floor(xlim[1]),
    to = ceiling(xlim[2]),
    by = (ceiling(xlim[2]) - floor(xlim[1])) / n
  )
  # SPLIT SORTED RAW DATA INTO CHUNKS
  raw_data_chunks <- lapply(seq_len(n - 1), function(z) {
    if (z < n - 1) {
      raw_data[raw_data[, channels[1]] >= chunks[z] &
        raw_data[, channels[1]] < chunks[z + 1], ]
    } else {
      raw_data[raw_data[, channels[1]] >= chunks[z] &
        raw_data[, channels[1]] <= chunks[z + 1], ]
    }
  })
  # COMPUTE MEDIANS IN BOTH CHANNELS PER CHUNK
  xmedians <- LAPPLY(
    raw_data_chunks,
    function(z) {
      # ONLY COMPUTE IF SUFFICIENT EVENTS
      if (nrow(z) > 30) {
        median(z[, channels[1]])
      } else {
        NA
      }
    }
  )
  xmedians <- xmedians[!is.na(xmedians)]
  ymedians <- LAPPLY(
    raw_data_chunks,
    function(z) {
      # ONLY COMPUTE IF SUFFICIENT EVENTS
      if (nrow(z) > 30) {
        median(z[, channels[2]])
      } else {
        NA
      }
    }
  )
  ymedians <- ymedians[!is.na(ymedians)]
  # PREPARE MEDFIs
  medians <- data.frame(xmedians, ymedians)
  colnames(medians) <- c(
    channels[1],
    channels[2]
  )
  vals_x <- medians[, channels[1]]
  vals_y <- medians[, channels[2]]
  loessMod <- loess(vals_y ~ vals_x,
    data = medians,
    span = 0.9
  )
  loessMod <- predict(loessMod)
  # ADD LINE TO PLOT
  lines(medians[, channels[1]],
    loessMod,
    col = "purple2",
    lwd = 3
  )
}
