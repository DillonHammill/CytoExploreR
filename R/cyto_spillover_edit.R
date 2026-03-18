## CYTO_SPILLOVER_EDIT ---------------------------------------------------------

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
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}} when a \code{flowSet} or
#'   \code{GatingSet} is supplied. Refer to \code{\link{cyto_select}} for more
#'   details.
#' @param channels names of the channels or markers in which compensation should
#'   be visualised, set to all area fluorescence parameters by default.
#'   \code{channels} can be ued to restrict the list of parameters that are
#'   displayed within the spillover editor.
#' @param channel_match logical indicating whether a call should be made to
#'   \code{cyto_channel_match()} to automatically detect an appropriate parental
#'   population for each control and to match each sample with a fluorescent
#'   channel. \code{channel_match} is set to TRUE by default, users can set this
#'   argument to FALSE if they supply samples stained with multiple antibodies.
#' @param spillover name of a square spillover matrix csv file or spillover
#'   matrix to edit. Setting \code{spill} to NULL (the default) will result in
#'   extraction of the spillover matrix generated on the cytometer which is
#'   attached to the samples. Similarly, if the supplied spillover matrix csv
#'   file does not exist, the spillover matrix attached to the first sample will
#'   be used and the edited spillover matrix will be saved to the specified
#'   file.
#' @param compensated logical required when a \code{cytoset} is supplied to
#'   indicate whether the supplied data has been compensated prior to passing it
#'   to \code{cyto_spillover_edit()}, set to FALSE by default.
#' @param save_as name of a csv file to which the edited spillover matrix should
#'   be written, set to \code{Spillover-Matrix.csv} prefixed with the date by
#'   default.
#' @param axes_trans an object of class \code{transformerList} containing
#'   transformers to used to transform the fluorescent channels of the samples
#'   for visualisation.
#' @param axes_limits options include \code{"auto"}, \code{"data"} or
#'   \code{"machine"} to use optimised, data or machine limits respectively. Set
#'   to \code{"machine"} by default to use entire axes ranges.
#' @param events numeric passed to \code{cyto_plot} to control the number of
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
#'   file. Saved filename defaults to \code{date-Spillover-Matrix.csv} if not
#'   specified.
#'
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarPanel selectInput
#'   checkboxInput actionButton mainPanel plotOutput reactiveValues observe
#'   eventReactive renderImage tabsetPanel tabPanel sidebarLayout fluidRow
#'   updateSelectInput onStop stopApp runApp updateCheckboxInput paneViewer icon
#'   span img NS reactive moduleServer observeEvent column
#'   updateCheckboxGroupInput checkboxGroupInput renderPlot tags tagList HTML
#' @importFrom rhandsontable rhandsontable rHandsontableOutput hot_to_r
#'   renderRHandsontable hot_cols hot_rows
#' @importFrom bslib bs_theme
#' @importFrom rhandsontable %>%
#' @importFrom stats median
#' @importFrom graphics lines layout
#' @importFrom flowWorkspace gs_cyto_data
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_spillover_compute}}
#' @seealso \code{\link{cyto_plot_compensation}}
#' @seealso \code{\link{cyto_plot}}
#'
#' @export
cyto_spillover_edit <- function(x,
                                select = NULL,
                                channels = NULL,
                                channel_match = TRUE,
                                spillover = NULL,
                                compensated = FALSE,
                                save_as = NULL,
                                axes_trans = NA,
                                axes_limits ="machine",
                                events = 2000,
                                point_size = 4,
                                axes_text_size = 1.7,
                                axes_label_text_size = 2,
                                title_text_size = 1.5,
                                header_text_size = 1.5,
                                viewer = FALSE,
                                ...) {
  
  # ENFORCE HIGH-PERFORMANCE RASTERIZATION -------------------------------------
  options(shiny.useragg = TRUE)
  
  # PREPARE DATA ---------------------------------------------------------------
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    stop("'x' must be either a cytoset or GatingSet object.")
  }
  
  # COPY SELECTED DATA
  x <- cyto_copy(
    cyto_select(
      x, 
      select
    )
  )
  
  # CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
    channels <- channels[
      !grepl(
        "\\-H$|\\-W$", 
        channels,
        ignore.case = TRUE
      )
    ]
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # SAMPLE NAMES
  nms <- cyto_names(x)
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  pd_order <- rownames(pd)
  
  # INTERACTIVE CHANNEL_MATCH
  if(isTRUE(channel_match) || is.character(channel_match)) {
    channel_match <- if(is.character(channel_match)){
      channel_match 
    }else {
      NULL
    }
    pd_new <- cyto_channel_match(
      x,
      channels = channels,
      file = channel_match
    )
    pd <- pd_new[match(pd_order, rownames(pd_new)), , drop = FALSE]
  # BYPASS CHANNEL_MATCH
  } else {
    if(!"group" %in% colnames(pd)){
      pd$group   <- rep(NA, nrow(pd))
    }
    if(!"parent" %in% colnames(pd)) {
      pd$parent  <- rep("root", nrow(pd))
    }
    if(!"channel" %in% colnames(pd)){
      pd$channel <- rep(NA, nrow(pd))
    }
  }
  
  # UPDATE EXPERIMENT DETAILS
  cyto_details(x) <- pd
  
  # TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_extract(x)
  }
  
  # APPLIED TRANSFORMERS
  if(!.all_na(axes_trans)) {
    attributes(axes_trans)$applied <- TRUE
  }
  
  # DEFAULT TRANSFORMERS
  if(.all_na(axes_trans)) {
    axes_trans <- cyto_transformers_define(
      x,
      parent = "root",
      channels = channels,
      type = "asinhx",
      plot = FALSE
    )
    attributes(axes_trans)$applied <- FALSE
  }
  
  # APPLIED SPILLOVER MATRICES
  if(cyto_class(x, "GatingSet")) {
    spill <- cyto_spillover_extract(x)
  } else {
    # CYTOSET COMPENSATION FLAG
    if(compensated) {
      spill <- cyto_spillover_extract(x)
    } else {
      spill <- NULL
    }
  }
  
  # DEFAULT SPILLOVER MATRIX TO EDIT
  if(is.null(spillover)) {
    # ATTEMPT TO EXTRACT SPILLOVER MATRIX FROM CYTOSET
    spillover <- cyto_spillover_extract(
      cyto_data_extract(
        x,
        parent = "root",
        format = "cytoset",
        copy = FALSE
      )[[1]]
    )
    # NO SPILLOVER MATRIX FOUND
    if(!is.null(spillover)) {
      spillover <- spillover[!ulapply(spillover, "is.null")][[1]]
      # TEMPLATE SPILLOVER MATRIX
    } else {
      spillover <- matrix(
        0,
        nrow = length(channels),
        ncol = length(channels),
        dimnames = list(
          channels,
          channels
        )
      )
      diag(spillover) <- 1
    }
    # PREPARE SUPPLIED SPILLOVER MATRIX
  } else {
    spillover <- .cyto_spillover_prepare(
      x,
      spillover = spillover
    )[[1]]
  }
  
  # Pre-computed constant: axes_trans and channels do not change after startup.
  trans_ind           <- if(any(channels %in% names(axes_trans))) match_ind(channels, names(axes_trans)) else NULL
  axes_trans_combined <- if(!is.null(trans_ind)) cyto_transformers_combine(axes_trans[trans_ind]) else NULL

  # X -> LINEAR SCALE
  if(attributes(axes_trans)$applied & !is.null(axes_trans_combined)) {
    x <- cyto_transform(
      x,
      trans = axes_trans_combined,
      inverse = TRUE,
      copy = TRUE,
      plot = FALSE,
      quiet = TRUE
    )
  }
  # X -> DECOMPENSATE
  if(!is.null(spill)) {
    x <- cyto_compensate(
      x,
      spillover = spill,
      remove = TRUE
    )
  }
  # X - LINEAR UNCOMPENSATED
  x_linear <- cyto_copy(x)
  # X - TRANSFORMED UNCOMPENSATED
  if(!is.null(axes_trans_combined)) {
    x <- cyto_transform(
      x,
      trans = axes_trans_combined,
      inverse = FALSE,
      copy = FALSE,
      plot = FALSE,
      quiet = TRUE
    )
  }
  
  # SHINY DEFAULTS -------------------------------------------------------------
  
  # UNSTAINED SAMPLE
  if(any(grepl("unstained", pd$channel, ignore.case = TRUE))) {
    NIL_select <- pd$name[grepl("Unstained", pd$channel, ignore.case = TRUE)][1]
  } else {
    NIL_select <- "None"
  }
  
  # SAMPLE
  if (!.all_na(pd$channel)) {
    # AVOID UNSTAINED CONTROL
    if (any(grepl("Unstained", pd$channel))) {
      ID_select <- pd$name[!grepl("Unstained",pd$channel,ignore.case = TRUE)][1]
      # USE FIRST SAMPLE
    } else {
      ID_select <- pd$name[1]
    }
  } else {
    ID_select <- pd$name[1]
  }
  
  # X CHANNEL - pre-compute initial value from ID_select so the channel
  # dropdown is seeded correctly on first flush, before ID() resolves
  xchan_select <- {
    ch <- pd$channel[pd$name == ID_select]
    if(length(ch) > 0 && !is.na(ch[1])) ch[1] else NULL
  }
  
  # EDITOR OPTIONS
  if (any(grepl("Unstained", pd$channel))) {
    editor_opts_select <- c("tracker", "line", "overlay")
  }else{
    editor_opts_select <- "tracker"
  }
  
  # PLOTS OPTIONS
  plots_opts_select <- c("uncompensated", "compensated")
  if (any(grepl("Unstained", pd$channel))) {
    plots_opts_select <- c("unstained",
                           plots_opts_select)
  }
  
  # SAVE_AS
  if(is.null(save_as)) {
    # USE SPILLOVER FILE NAME
    if(is.character(spillover)) {
      save_as <- spillover
    } else {
      save_as <- cyto_file_name(
        paste0(
          format(Sys.Date(), "%d%m%y"),
          "-", "Spillover-Matrix.csv"
        )
      )
    }
  }
  # SHINY APPLICATION ----------------------------------------------------------
  app <- shinyApp(
    ui <- fluidPage(
      theme = bs_theme(bootswatch = "yeti", version = 3),
      tags$head(tags$style(HTML("
        .shiny-plot-output { position: relative; }
        .shiny-plot-output.recalculating { opacity: 0.4; }
        .shiny-plot-output.recalculating::after {
          content: '';
          display: block;
          position: absolute;
          top: 50%;
          left: 50%;
          width: 50px;
          height: 50px;
          margin: -25px 0 0 -25px;
          border-radius: 50%;
          border: 5px solid #ccc;
          border-top-color: #2196F3;
          animation: cyto-spin 0.8s linear infinite;
        }
        @keyframes cyto-spin { to { transform: rotate(360deg); } }
      "))),
      titlePanel(span(img(src = CytoExploreR_logo(), width = 35), "CytoExploreR Spillover Matrix Editor")),
      tabsetPanel(id = "tabs",
        tabPanel("Editor", fluid = TRUE, sidebarLayout(
          sidebarPanel(width = 3,
                       cytoSelectUI("editor_select_unst", label = "Select unstained sample:"),
                       cytoSelectUI("editor_select", label = "Select sample:"),
                       nodeSelectUI("editor_node_select", label = "Select parent:"),
                       channelSelectUI("editor_xchannel", label = "X axis:"),
                       channelSelectUI("editor_ychannel", label = "Y axis:"),
                       optionsUI("editor_options", label = NULL, selected = editor_opts_select,
                                 choiceNames = list("Overlay unstained control", "Unstained control median", "Median tracker"),
                                 choiceValues = list("overlay", "line", "tracker")),
                       spillSaveUI("editor_save"),
                       actionButton("close", "Close")
          ),
          mainPanel(width = 9, spillEditUI("editor_spill"), editPlotUI("editor_plot"))
        )),
        tabPanel("Plots", fluid = TRUE, sidebarLayout(
          sidebarPanel(width = 3,
                       cytoSelectUI("plots_select_unst", label = "Select unstained sample:"),
                       cytoSelectUI("plots_select", label = "Select sample:"),
                       nodeSelectUI("plots_node_select", label = "Select parent:"),
                       channelSelectUI("plots_xchannel", label = "Channel:"),
                       optionsUI("plots_options", label = NULL, selected = plots_opts_select,
                                 choiceNames = list("Overlay unstained control", "Overlay compensated data", "Fit robust linear models"),
                                 choiceValues = list("unstained", "compensated", "models"))
          ),
          mainPanel(width = 9, compPlotUI("plots"))
        ))
      )
    ),
    
    server <- function(input, output, session) {
      
      editor_parent <- nodeSelectServer("editor_node_select", data = reactive(x), selected = reactive("root"))
      NIL <- cytoSelectServer("editor_select_unst", data = reactive(x), choices = "None", selected = reactive(NIL_select))
      ID <- cytoSelectServer("editor_select", data = reactive(x), selected = reactive(ID_select))

      # Auto-populate x channel from pData channel column for selected sample.
      # Falls back to xchan_select (pre-computed from ID_select) before ID()
      # resolves so channelSelectServer gets a correct initial value on the
      # first flush, preventing a spurious extra render with the wrong channel.
      selected_xchan <- reactive({
        id <- tryCatch(ID(), error = function(e) NULL)
        if(is.null(id)) return(xchan_select)
        ch <- pd$channel[pd$name == id]
        if(length(ch) > 0 && !is.na(ch[1])) ch[1] else NULL
      })
      xchan <- channelSelectServer("editor_xchannel", data = reactive(x), selected = selected_xchan)
      ychan <- channelSelectServer("editor_ychannel", data = reactive(x))
      editor_opts <- optionsServer("editor_options", selected = reactive(editor_opts_select))
      
      NIL_linear_subset <- reactive({
        req(NIL(), NIL() != "None")
        cyto_data_extract(x_linear, select = NIL(), parent = editor_parent(), events = events, copy = TRUE)[[1]]
      })

      ID_linear_subset <- reactive({
        req(ID())
        cyto_data_extract(x_linear, select = ID(), parent = editor_parent(), events = events, copy = TRUE)[[1]]
      })
      
      spill_edit <- spillEditServer("editor_spill", data = reactive(x), spill = reactive(spillover), xchan = xchan, ychan = ychan)
      
      NIL_compensated <- reactive({
        req(NIL_linear_subset(), spill_edit())
        cyto_compensate(NIL_linear_subset(), spillover = spill_edit())
      })
      
      ID_compensated <- reactive({
        req(ID_linear_subset(), spill_edit())
        cyto_compensate(ID_linear_subset(), spillover = spill_edit())
      })
      
      NIL_comp_trans <- reactive({
        req(NIL_compensated())
        if(!is.null(axes_trans_combined)) {
          cyto_transform(NIL_compensated(), trans = axes_trans_combined,
                         inverse = FALSE, copy = FALSE, plot = FALSE, quiet = TRUE)
        } else NIL_compensated()
      })

      ID_comp_trans <- reactive({
        req(ID_compensated())
        if(!is.null(axes_trans_combined)) {
          cyto_transform(ID_compensated(), trans = axes_trans_combined,
                         inverse = FALSE, copy = FALSE, plot = FALSE, quiet = TRUE)
        } else ID_compensated()
      })
      
      editPlotServer("editor_plot",
                     ID_comp_trans = ID_comp_trans, NIL_comp_trans = NIL_comp_trans,
                     NIL = NIL,
                     parent = editor_parent, opts = editor_opts, xchan = xchan, ychan = ychan,
                     axes_trans = axes_trans, axes_limits = axes_limits, events = events,
                     point_size = point_size, axes_text_size = axes_text_size,
                     axes_label_text_size = axes_label_text_size, title_text_size = title_text_size, ...)
      
      plots_NIL <- cytoSelectServer("plots_select_unst", data = reactive(x), choices = "None", selected = NIL)
      plots_ID <- cytoSelectServer("plots_select", data = reactive(x), selected = ID)
      plots_parent <- nodeSelectServer("plots_node_select", data = reactive(x), selected = editor_parent)
      plots_xchan <- channelSelectServer("plots_xchannel", data = reactive(x), selected = xchan)
      plots_opts <- optionsServer("plots_options", selected = reactive(plots_opts_select))

      compPlotServer("plots",
                     ID_comp_trans = ID_comp_trans, NIL_comp_trans = NIL_comp_trans,
                     opts = plots_opts, xchan = plots_xchan, channels = reactive(channels),
                     spillover = spill_edit, axes_trans = axes_trans, axes_limits = axes_limits,
                     events = events, point_size = point_size, axes_text_size = axes_text_size,
                     axes_label_text_size = axes_label_text_size, title_text_size = title_text_size,
                     active_tab = reactive(input$tabs), ...)
      
      spillSaveServer("editor_save", spill = spill_edit, save_as = save_as)
      
      observeEvent(input$close, {
        write_to_csv(spill_edit(), save_as)
        stopApp(read_from_csv(save_as))
      })
      
      onStop(function() stopApp(read_from_csv(save_as)))
    }
  )
  
  if(viewer) runApp(app, launch.browser = paneViewer(), quiet = TRUE) else runApp(app, quiet = TRUE)
}

# UI & SERVER MODULES ----------------------------------------------------------

#' @noRd
cytoSelectUI <- function(id, label = NULL, ...) {
  selectInput(NS(id, "select"), label = label, choices = NULL, ...)
}

#' @noRd
cytoSelectServer <- function(id, data = reactive(NULL), choices = NULL, selected = reactive(NULL)) {
  moduleServer(id, function(input, output, session){
    pending_sel <- reactiveVal(NULL)
    observe({
      sel <- selected()
      pending_sel(sel)
      updateSelectInput(session, "select", choices = c(cyto_names(data()), choices), selected = sel)
    })
    observeEvent(input$select, {
      ps <- pending_sel()
      if(is.null(ps) || input$select != ps) {
        pending_sel(input$select)
      }
    }, ignoreInit = TRUE)
    return(reactive({ if(!is.null(pending_sel())) pending_sel() else input$select }))
  })
}

#' @noRd
channelSelectUI <- function(id, label = NULL, ...) {
  fluidRow(
    column(10, style = "padding-right: 5px;", selectInput(NS(id, "select"), label = label, choices = NULL)),
    column(1, style = "padding-left: 4px; padding-right: 0px; margin-top: 26px; margin-left: 0px; margin-right: 0px;",
           actionButton(NS(id, "down"), label = NULL, icon = icon("arrow-down"), style = "background-color: #FF0000; padding: 4px;")),
    column(1, style = "padding-left: 0px; padding-right: 2px; margin-top: 26px; margin-left: 0px; margin-right: 0px;",
           actionButton(NS(id, "up"), label = NULL, icon = icon("arrow-up"), style = "background-color: #33CC33; padding: 4px;"))
  )
}

#' @noRd
channelSelectServer <- function(id, data = reactive(NULL), selected = reactive(NULL), exclude = reactive(NULL), ...) {
  moduleServer(id, function(input, output, session) {

    channels_excl <- reactive({
      req(data())
      chans <- unname(cyto_fluor_channels(data()))
      if(!.empty(exclude(), null = TRUE)) chans <- chans[-match(exclude(), chans)]
      chans
    })

    # Server-side effective selection — updated synchronously when selected()
    # changes so downstream reactives see the new channel in the same flush
    # cycle as the sample switch, eliminating an intermediate render with the
    # wrong channel before the browser round-trip completes.
    pending_sel <- reactiveVal(NULL)

    observe({
      chans <- channels_excl()
      sel <- selected()
      if(!is.null(sel)) {
        effective <- if(is.na(sel)) chans[1] else sel
        pending_sel(effective)
        updateSelectInput(session, "select", choices = chans, selected = effective)
      } else {
        # Seed pending_sel with the first channel so downstream reactives have
        # a valid value before the browser's first round-trip populates
        # input$select.  Use isolate() so this read does not add pending_sel
        # as a reactive dependency of this observer (which would re-trigger it
        # every time pending_sel changes and create a loop).
        if(is.null(isolate(pending_sel()))) {
          pending_sel(chans[1])
        }
        updateSelectInput(session, "select", choices = chans)
      }
    })

    # User manually changed the dropdown — sync pending_sel only when the new
    # value differs from what was last set programmatically, so that the
    # updateSelectInput round-trip is not treated as a user action.
    observeEvent(input$select, {
      ps <- pending_sel()
      if(is.null(ps) || input$select != ps) {
        pending_sel(input$select)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$up, {
      req(input$select, channels_excl())
      ind <- match(input$select, channels_excl())
      new_ind <- if(ind < length(channels_excl())) ind + 1 else 1
      new_chan <- channels_excl()[new_ind]
      pending_sel(new_chan)
      updateSelectInput(session, "select", choices = channels_excl(), selected = new_chan)
    })

    observeEvent(input$down, {
      req(input$select, channels_excl())
      ind <- match(input$select, channels_excl())
      new_ind <- if(ind > 1) ind - 1 else length(channels_excl())
      new_chan <- channels_excl()[new_ind]
      pending_sel(new_chan)
      updateSelectInput(session, "select", choices = channels_excl(), selected = new_chan)
    })

    return(reactive({ if(is.null(pending_sel())) input$select else pending_sel() }))
  })
}

#' @noRd
nodeSelectUI <- function(id, label = NULL, ...) {
  selectInput(NS(id, "select"), label = label, choices = "root", ...)
}

#' @noRd
nodeSelectServer <- function(id, data = reactive(NULL), choices = NULL, selected = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    observe({
      if(cyto_class(data(), "GatingSet")) {
        updateSelectInput(session, "select", choices = cyto_nodes(data(), path = "auto"), selected = selected())
      }
    })
    return(reactive({ req(input$select); input$select }))
  })
}

#' @noRd
spillEditUI <- function(id, height = "300px", ...) {
  tagList(
    tags$script(HTML(paste0(
      "if (!window._cytoHighlight) window._cytoHighlight = {};",
      "Shiny.addCustomMessageHandler('selectHotCell_", id, "', function(msg) {",
      "  window._cytoHighlight[msg.id] = { row: msg.row, col: msg.col };",
      "  (function tryRender(n) {",
      "    var el = document.getElementById(msg.id);",
      "    if (el) {",
      "      var hw = HTMLWidgets.getInstance(el);",
      "      if (hw && hw.hot) {",
      "        hw.hot.scrollViewportTo(msg.row, msg.col);",
      "        hw.hot.render();",
      "        return;",
      "      }",
      "    }",
      "    if (n > 0) { setTimeout(function() { tryRender(n - 1); }, 100); }",
      "  })(10);",
      "});"
    ))),
    rHandsontableOutput(NS(id, "spill"), height = height, width = "99%", ...)
  )
}

#' @noRd
spillEditServer <- function(id, data = reactive(NULL), spill = reactive(NULL), xchan = reactive(NULL), ychan = reactive(NULL), ...) {
  moduleServer(id, function(input, output, session){
    
    # source_version increments only on external spill changes, not user edits.
    # renderRHandsontable depends solely on source_version so user edits do not
    # trigger a full table re-render.
    values <- reactiveValues(spill = NULL, source_version = 0L)

    spill_mat <- reactive({
      if(.empty(spill(), null = TRUE)) {
        if(!.empty(data(), null = TRUE)) {
          sp <- cyto_spillover_extract(cyto_data_extract(data(), parent = "root", copy = FALSE)[[1]])[[1]]
        } else sp <- NULL
      } else {
        if(is.character(spill())) {
          sp <- if(file_ext(spill()) %in% "mtx") read_from_mtx(spill()) else read_from_csv(spill(), data.table = FALSE)
        } else sp <- spill()
      }
      if(!.empty(sp, null = TRUE)) {
        sp <- sp * 100
        rownames(sp) <- colnames(sp)
      }
      return(sp)
    })

    observe({
      values$spill <- spill_mat()
      values$source_version <- isolate(values$source_version) + 1L
    })

    # Push cell selection to the browser without triggering a full table re-render.
    # values$spill is read with isolate() so user edits to the table do not
    # re-fire this observer — hot.render() called from JS can otherwise
    # propagate back to input$spill -> values$spill -> here, creating a loop.
    # values$source_version is read non-isolated so the observer re-fires when
    # the table first loads (spill was NULL when xchan/ychan initially resolved).
    # The row/col indices depend only on matrix dimension (channel names), which
    # are fixed for the lifetime of the editor session.
    observe({
      req(xchan(), ychan())
      values$source_version
      sp <- isolate(values$spill)
      req(sp)
      session$sendCustomMessage(
        paste0("selectHotCell_", id),
        list(
          id  = session$ns("spill"),
          row = match(xchan(), rownames(sp)) - 1L,
          col = match(ychan(), colnames(sp)) - 1L
        )
      )
    })

    output$spill <- renderRHandsontable({
      values$source_version
      sp <- isolate(values$spill)
      req(sp)
      hot_id <- session$ns("spill")
      rhandsontable(sp, rowHeaderWidth = 105, readOnly = FALSE, manualColumnResize = TRUE) %>%
        hot_cols(type = "numeric", colWidths = 105, format = "0.000", halign = "htCenter",
                 renderer = paste0("
                   function (instance, td, row, col, prop, value, cellProperties) {
                     Handsontable.renderers.TextRenderer.apply(this, arguments);
                     var hl = window._cytoHighlight && window._cytoHighlight['", hot_id, "'];
                     if (hl && row === hl.row && col === hl.col) {
                       td.style.border = 'solid'; td.style.borderWidth = '3px'; td.style.borderColor = 'black';
                     } else {
                       td.style.border = '';
                     }
                     if(value < 0) td.style.background = 'lightblue';
                     else if (value == 0) td.style.background = 'white';
                     else if (value > 0 && value <= 10) td.style.background = 'lightgreen';
                     else if (value > 10 && value <= 25) td.style.background = 'yellow';
                     else if (value > 25 && value <= 50) td.style.background = 'orange';
                     else if (value > 50 && value < 100) td.style.background = 'red';
                     else if (value == 100) td.style.background = 'darkgrey';
                     else if (value > 100) td.style.background = 'violet';
                   }")) %>%
        hot_rows(rowHeights = 20)
    })
    
    observeEvent(input$spill, {
      sp <- hot_to_r(input$spill)
      rownames(sp) <- colnames(sp)
      if(any(is.na(sp))) {
        if(is.null(values$spill)) sp[is.na(sp)] <- 0
        else sp[is.na(sp)] <- values$spill[is.na(sp)]
      }
      # Skip update if value is unchanged — prevents spurious re-renders when
      # rhandsontable sends back its initial value after the first browser render.
      if(is.null(values$spill) || !isTRUE(all.equal(sp, values$spill, check.attributes = FALSE))) {
        values$spill <- sp
      }
    })
    
    spill_reactive <- reactive({
      req(values$spill)
      return(values$spill / 100)
    })
    
    return(spill_reactive %>% debounce(400))
  })
}

#' @noRd
optionsUI <- function(id, ...) {
  checkboxGroupInput(NS(id, "options"), ...)
}

#' @noRd
optionsServer <- function(id, selected = reactive(NULL), ...) {
  moduleServer(id, function(input, output, session){
    observe({
      if(!.empty(selected(), null = TRUE)) {
        updateCheckboxGroupInput(session, "options", selected = selected())
      }
    })
    return(reactive({ input$options }))
  })
}

#' @noRd
spillSaveUI <- function(id, ...){
  actionButton(NS(id, "save"), "Save", ...)
}

#' @noRd
spillSaveServer <- function(id, spill = reactive(NULL), save_as = NULL, ...) {
  moduleServer(id, function(input, output, session){
    observeEvent(input$save, { write_to_csv(spill(), save_as) })
    return(spill)
  })
}

#' @noRd
editPlotUI <- function(id, ...) {
  plotOutput(NS(id, "plot"), width = "70%", ...)
}

#' @noRd
editPlotServer <- function(id, ID_comp_trans = reactive(NULL), NIL_comp_trans = reactive(NULL),
                           NIL = reactive(NULL),
                           opts = reactive(NULL), xchan = reactive(NULL), ychan = reactive(NULL),
                           axes_trans = NA, axes_limits = "machine", events = 2000,
                           point_size = 3, axes_text_size = 1.7, axes_label_text_size = 2,
                           title_text_size = 1.5, ...) {
  moduleServer(id, function(input, output, session){

    # Cache nil_data as reactive — avoid tryCatch inside renderPlot
    nil_data <- reactive({
      tryCatch(NIL_comp_trans(), error = function(e) NULL)
    })

    # Pre-compute overlay flag — only invalidates when nil_data or opts change
    use_overlay <- reactive({
      !is.null(nil_data()) && "overlay" %in% opts()
    })

    # Cache NIL medians — only recomputes when nd changes, not on channel navigation
    nil_medians <- reactive({
      nd <- nil_data()
      req(nd)
      cyto_apply(nd, "cyto_stat_median", input = "matrix", inverse = FALSE, copy = FALSE)
    })

    output$plot <- renderPlot({
      req(ID_comp_trans(), xchan(), ychan())

      nd <- nil_data()

      # If an unstained control is selected, gate on its data being ready.
      # Without this, the plot renders once with nd=NULL (before NIL()'s
      # browser round-trip completes) and again immediately after NIL_comp_trans
      # resolves, producing an unnecessary double render on load.
      nil_name <- tryCatch(NIL(), error = function(e) NULL)
      if (!is.null(nil_name) && nil_name != "None") {
        req(nd)
      }

      overlay <- use_overlay()

      cyto_plot_custom(layout = matrix(c(1, 1, 2, 2, 1, 1, 3, 3), byrow = TRUE, ncol = 4), popup = FALSE)
      
      suppressPrint(
        cyto_plot(
          ID_comp_trans(),
          channels = c(xchan(), ychan()),
          overlay = if(overlay) nd else NA,
          axes_trans = axes_trans, axes_limits = axes_limits, events = events,
          title = cyto_names(ID_comp_trans()), point_size = point_size,
          axes_text_size = axes_text_size, axes_label_text_size = axes_label_text_size,
          title_text_size = title_text_size, key_text_size = 2, key_title_text_size = 2,
          popup = FALSE, margins = c(4,5,6,6.5), ...
        )
      )
      
      usr <- .par("usr")[[1]]
      
      if(!is.null(nd) && "line" %in% opts()) {
        medFI <- tryCatch(nil_medians(), error = function(e) NULL)
        if(!is.null(medFI)) abline(h = medFI[1, ychan()], col = "red", lwd = 2)
      }
      
      if("tracker" %in% opts()) {
        .cyto_median_tracker(ID_comp_trans(), c(xchan(), ychan()))
      }
      
      for(i in seq_len(2)) {
        chan <- c(xchan(), ychan())[i]
        lims <- if(i == 1) usr[1:2] else usr[3:4]
        
        if(chan %in% names(axes_trans)) {
          label_text_x <- .cyto_transform(min(lims) + 0.90 * diff(lims), trans = axes_trans, channel = chan, inverse = TRUE)
        } else {
          label_text_x <- min(lims) + 0.90 * diff(lims)
        }
        
        suppressPrint(
          cyto_plot(
            if(overlay) nd else ID_comp_trans(),
            channels = chan,
            overlay = if(overlay) ID_comp_trans() else NA,
            axes_trans = axes_trans, axes_limits = axes_limits, events = events,
            title = NA, ylab = "Density",
            hist_fill = if(overlay) c("grey70", "white") else "white",
            hist_fill_alpha = if(overlay) c(1, 0) else 1,
            hist_line_col = if(overlay) c("black", "blue") else "blue",
            hist_line_width = 2, axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size, title_text_size = title_text_size,
            popup = FALSE,
            label_text = if(chan == ychan()) (if(overlay) c("MedFI", "MedFI") else c("MedFI")) else NA,
            label_stat = if(chan == ychan()) (if(overlay) c("median", "median") else c("median")) else NA,
            label_text_x = if(chan == ychan()) (if(overlay) rep(label_text_x, 2) else label_text_x) else NA,
            label_text_y = if(chan == ychan()) (if(overlay) c(80, 30) else c(50)) else NA,
            label_text_col = if(chan == ychan()) (if(overlay) c("grey40", "blue") else c("blue")) else "black",
            label_text_size = 1.1, margins = c(5,6,2,4), ...
          )
        )
      }
      cyto_plot_complete()
    })
  })
}

#' @noRd
compPlotUI <- function(id, ...) {
  plotOutput(NS(id, "cyto_plot_comp"), height = "800px", ...)
}

#' @noRd
compPlotServer <- function(id, ID_comp_trans = reactive(NULL), NIL_comp_trans = reactive(NULL),
                           opts = reactive(NULL), xchan = reactive(NULL), channels = reactive(NULL),
                           spillover = reactive(NULL), axes_trans = NA, axes_limits = "machine",
                           events = 2000, point_size = 3, axes_text_size = 1.7,
                           axes_label_text_size = 2, title_text_size = 1.5,
                           active_tab = reactive(NULL), ...) {
  moduleServer(id, function(input, output, session){

    layout_dims <- reactive({ c(ceiling(length(channels())/4), 4) })

    # Cache nil_data as reactive — avoid tryCatch inside renderPlot
    nil_data <- reactive({
      tryCatch(NIL_comp_trans(), error = function(e) NULL)
    })

    # Cache cytoset assembly — avoid constructing new cytoset every render
    comp_cs <- reactive({
      req(ID_comp_trans(), xchan())
      nd <- nil_data()
      if(!is.null(nd)) {
        cf_list <- structure(list(ID_comp_trans()[[1]], nd[[1]]),
                             names = c(cyto_names(ID_comp_trans()), cyto_names(nd)))
        cs <- cytoset(cf_list[!ulapply(cf_list, "is.null")])
      } else {
        cs <- ID_comp_trans()
      }
      pd <- cyto_details(cs)
      pd$channel <- c(xchan(), "unstained")[seq_along(cs)]
      cyto_details(cs) <- pd
      cs
    })

    # Pre-compute overlay colours — only recalculated when opts() changes
    plot_colours <- reactive({
      o <- opts()
      has_unst <- "unstained" %in% o
      has_comp <- "compensated" %in% o
      cols <- if(has_unst && has_comp) {
        c("magenta", "blue", "grey40")
      } else if(has_unst) {
        c("magenta", "grey40")
      } else if(has_comp) {
        c("magenta", "blue")
      } else {
        c("magenta")
      }
      overlay <- if(!has_unst && !has_comp) "none" else o[o %in% c("unstained", "compensated")]
      list(cols = cols, overlay = overlay)
    })

    output$cyto_plot_comp <- renderPlot({
      req(comp_cs(), channels(), isTRUE(active_tab() == "Plots"))
      pc <- plot_colours()

      suppressPrint(
        suppressWarnings(
          cyto_plot_compensation(
            comp_cs(), channels = channels(), channel_match = cyto_details(comp_cs()),
            overlay = pc$overlay,
            spillover = spillover(), compensated = TRUE, axes_trans = axes_trans,
            axes_limits = axes_limits, events = events, point_size = point_size,
            point_col = pc$cols, hist_fill = pc$cols,
            lines = "models" %in% opts(), text = TRUE, text_size = 1.5,
            axes_text_size = 1.7, axes_label_text_size = 2, title_text_size = 2,
            layout = layout_dims(), popup = FALSE, ...
          )
        )
      )
    }, height = reactive({ 250 * layout_dims()[1] }))
  })
}

# ------------------------------------------------------------------------------
# OPTIMIZED MEDIAN TRACKER (C++ INTEGRATION)
# ------------------------------------------------------------------------------

#' Median Tracker (Accelerated via Rcpp)
#' Add median tracker to plot utilizing C++ binning and calculation
#' @param x object of class flowFrame
#' @param channels channels used to construct the plot
#' @noRd
.cyto_median_tracker <- function(x, channels = NULL) {
  
  # Extract numeric matrix directly
  raw_data <- cyto_data_extract(
    x, 
    format = "matrix", 
    channels = channels, 
    copy = FALSE
  )[[1]][[1]]
  
  if(nrow(raw_data) > 0) {
    # Execute C++ binning and median algorithm
    medians <- binned_median_cpp(
      x = raw_data[, channels[1]], 
      y = raw_data[, channels[2]], 
      n_bins = 25
    )
    
    # Filter NA outputs originating from insufficient bin counts (n < 30)
    medians <- medians[!is.na(medians$x) & !is.na(medians$y), ]
    
    # Proceed strictly if sufficient degrees of freedom exist for LOESS
    if(nrow(medians) > 3) {
      loessMod <- suppressWarnings(
        loess(
          y ~ x, 
          data = medians, 
          span = 0.9
        )
      )
      
      lines(
        medians$x, 
        predict(loessMod), 
        col = "purple2", 
        lwd = 3
      )
    }
  }
}