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
#'   updateCheckboxGroupInput checkboxGroupInput
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
  
  # PREPARE DATA ---------------------------------------------------------------
  
  # CYTOSET/GATINGSET
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    stop(
      "'x' must be either a cytoset or GatingSet object!"
    )
  }
  
  # COPY
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
  
  # INTERACTIVE CHANNEL_MATCH
  if(isTRUE(channel_match) | is.character(channel_match)) {
    pd <- cyto_channel_match(
      x,
      channels = channels,
      file = if(isTRUE(channel_match)) {
        NULL
      } else{
        channel_match
      }
    )
    pd <- pd[match(rownames(cyto_details(x)), rownames(pd)), , drop = FALSE]
  # BYPASS CHANNEL_MATCH
  } else {
    # BYPASS CYTO_CHANNEL_MATCH TO ALLOW FULL STAINED SAMPLES
    lapply(
      c("group", "parent", "channel"),
      function(z) {
        if(!z %in% colnames(pd)) {
          if(z == "parent") {
            pd[, z] <<- rep("root", nrow(pd))
          } else {
            pd[, z] <<- rep(NA, nrow(pd))
          }
        }
      }
    )
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
      type = "biex",
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
      spillover <- spillover[!LAPPLY(spillover, "is.null")][[1]]
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
  
  # X -> LINEAR SCALE
  if(attributes(axes_trans)$applied & any(channels %in% names(axes_trans))) {
    x <- cyto_transform(
      x,
      trans = cyto_transformers_combine(
        axes_trans[match_ind(channels, names(axes_trans))]
      ),
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
  if(any(channels %in% names(axes_trans))) {
    x <- cyto_transform(
      x,
      trans = cyto_transformers_combine(
        axes_trans[match_ind(channels, names(axes_trans))]
      ),
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
  
  # X CHANNEL
  xchan_select <- NULL
  
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
  
  # APPLICATION
  app <- shinyApp(
    # USER INTERFACE
    ui <- fluidPage(
      theme = bs_theme(
        bootswatch = "yeti",
        version = 3
      ),
      titlePanel(
        span(
          img(
            src = CytoExploreR_logo(),
            width = 35
          ), 
          "CytoExploreR Spillover Matrix Editor"
        )
      ),
      tabsetPanel(
        tabPanel(
          "Editor",
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              width = 3,
              cytoSelectUI("editor_select_unst",
                           label = "Select unstained sample:"),
              cytoSelectUI("editor_select",
                           label = "Select sample:"),
              nodeSelectUI("editor_node_select",
                           label = "Select parent:"),
              channelSelectUI("editor_xchannel",
                              label = "X axis:"),
              channelSelectUI("editor_ychannel",
                              label = "Y axis:"),
              optionsUI("editor_options",
                        label = NULL,
                        selected = editor_opts_select,
                        choiceNames = list(
                          "Overlay unstained control",
                          "Unstained control median",
                          "Median tracker"
                        ),
                        choiceValues = list(
                          "overlay",
                          "line",
                          "tracker"
                        )),
              spillSaveUI("editor_save")
            ),
            mainPanel(
              width = 9,
              spillEditUI(
                "editor_spill"
              ),
              editPlotUI(
                "editor_plot"
              )
            )
          )
        ),
        tabPanel(
          "Plots",
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              width = 3,
              cytoSelectUI(
                "plots_select_unst",
                label = "Select unstained sample:"
              ),
              cytoSelectUI(
                "plots_select",
                label = "Select sample:"
              ),
              nodeSelectUI(
                "plots_node_select",
                label = "Select parent:"
              ),
              channelSelectUI(
                "plots_xchannel",
                label = "Channel:"
              ),
              optionsUI(
                "plots_options",
                 label = NULL,
                 selected = plots_opts_select,
                 choiceNames = list(
                   "Overlay unstained control",
                   "Overlay compensated data",
                   "Fit robust linear models"
                 ),
                 choiceValues = list(
                   "unstained",
                   "compensated",
                   "models"
                 )
              )
            ),
            mainPanel(
              width = 9,
              compPlotUI(
                "plots"
              )
            )
          )
        )
      )
    ),
    # SERVER
    server <- function(input, 
                       output, 
                       session) {
      
      # STORAGE ----------------------------------------------------------------
      
      # VALUES
      values <- reactiveValues(
        parent = "root",
        spill = spillover,
        NIL_select = NIL_select,
        NIL_comp_trans = NULL,
        ID_select = ID_select,
        ID_comp_trans = NULL,
        xchan_select = xchan_select,
        editor_opts_select = editor_opts_select,
        plots_opts_select = plots_opts_select
      )
      
      # PARENT -----------------------------------------------------------------
      
      # PARENT DROPDOWN
      editor_parent <- nodeSelectServer(
        "editor_node_select",
        data = reactive(x),
        selected = reactive(values$parent)
      )
      
      # PARENT MANUALLY SELECTED
      observe({
        values$parent <- editor_parent()
      })
      
      # PARENT INHERITED FROM SELECTED SAMPLE
      observe({
        if(!.empty(ID(), null = TRUE)) {
          values$parent <- cyto_details(
            cyto_select(
              x,
              ID(),
              exact = TRUE
            )
          )$parent
        }
      })
      
      # UNSTAINED CONTROL ------------------------------------------------------
      
      # NAME - UNSTAINED CONTROL
      NIL <- cytoSelectServer(
        "editor_select_unst", 
        data = reactive(x),
        choices = "None",
        selected = reactive(values$NIL_select)
      )
      
      # NIL UPDATE
      observe({
        if(!.empty(NIL(), null = TRUE)) {
          values$NIL_select <- NIL()
        }
      })
      
      # UNSTAINED CONTROL - COMPENSATED & TRANSFORMED
      observe({
        # UNSTAINED CONTROL SELECTED  
        if(!.empty(values$NIL_select, null = TRUE) & 
           values$NIL_select != "None") {
          # UNSTAINED CONTROL - LINEAR DATA
          NIL_linear <- cyto_data_extract(
            x_linear,
            select = values$NIL_select,
            parent = values$parent,
            copy = TRUE
          )[[1]]
          # COMPENSATE
          if(!is.null(values$spill)) {
            NIL_linear <- cyto_compensate(
              NIL_linear,
              spillover = values$spill
            )
          }
          # TRANSFORM
          if(any(channels %in% names(axes_trans))) {
            values$NIL_comp_trans <- cyto_transform(
              NIL_linear,
              trans = cyto_transformers_combine(
                axes_trans[match_ind(channels, names(axes_trans))]
              ),
              inverse = FALSE,
              copy = FALSE,
              plot = FALSE,
              quiet = TRUE
            )
          }
        # NO UNSTAINED CONTROL SELECTED
        } else {
          values$NIL_comp_trans <- NULL
        }
      })
      
      # SELECTED CONTROL -------------------------------------------------------
      
      # NAME - SELECTED CONTROL
      ID <- cytoSelectServer(
        "editor_select",
        data = reactive(x),
        selected = reactive(values$ID_select)
      )
      
      # ID UPDATE
      observe({
        if(!.empty(ID(), null = TRUE)) {
          values$ID_select <- ID()
        }
      })
      
      # SELECTED CONTROL - COMPENSATED & TRANSFORMED
      observe({
        # SAMPLE SELECTED
        if(!.empty(values$ID_select, null = TRUE)) {
          # # PARENT
          # values$parent <- cyto_details(x_linear[values$ID_select])[, "parent"]
          # X CHANNEL
          values$xchan_select <- cyto_details(
            cyto_select(
              x,
              values$ID_select,
              exact = TRUE
            )
          )$channel
          # LINEAR DATA
          ID_linear <- cyto_data_extract(
            x_linear,
            select = values$ID_select,
            parent = values$parent,
            copy = TRUE
          )[[1]]
          # COMPENSATE
          if(!is.null(values$spill)) {
            ID_linear <- cyto_compensate(
              ID_linear,
              spillover = values$spill
            )
          }
          # TRANSFORM
          if(any(channels %in% names(axes_trans))) {
            values$ID_comp_trans <- cyto_transform(
              ID_linear,
              trans = cyto_transformers_combine(
                axes_trans[match_ind(channels, names(axes_trans))]
              ),
              inverse = FALSE,
              copy = FALSE,
              plot = FALSE,
              quiet = TRUE
            )
          }
        # NO SAMPLE SELECTED
        } else {
          values$ID_comp_trans <- NULL
        }
      })
      
      # CHANNELS ---------------------------------------------------------------
      
      # X CHANNEL
      xchan <- channelSelectServer(
        "editor_xchannel",
        data = reactive(x),
        selected = reactive(values$xchan_select)
      )
      
      # # X CHANNEL UPDATE
      # observe({
      #   values$xchan_select <- pd$channel[cyto_match(x, values$ID_select)]
      # })
      
      observe({
        values$xchan_select <- xchan()
      })
      
      # Y CHANNEL
      ychan <- channelSelectServer(
        "editor_ychannel",
        data = reactive(x)
        # exclude = xchan # AWESOME BUT TIME CONSUMING
      )
      
      # SPILLOVER MATRIX -------------------------------------------------------
      
      # EDIT SPILLOVER MATRIX
      spill_edit <- spillEditServer(
        "editor_spill",
        data = reactive(x),
        spill = reactive({values$spill}),
        xchan = xchan,
        ychan = ychan
      )
      
      # STORE EDITS TO SPILLOVER MATRIX
      observe({
        values$spill <- spill_edit()
      })
      
      # OPTIONS ----------------------------------------------------------------
      
      # EDITOR OPTIONS
      editor_opts <- optionsServer(
        "editor_options",
        selected = reactive(values$editor_opts_select)
      )
      
      # EDITOR OPTIONS UPDATE
      observe({
        if(.empty(values$NIL_comp_trans, null = TRUE)) {
          values$editor_opts_select <- values$editor_opts_select[
            !values$editor_opts_select %in% c("overlay", "line")
          ]
        } else {
          values$editor_opts_select <- unique(
            c(values$editor_opts_select,
              "overlay",
              "line")
          )
        }
      })
      
      # SCATTER/HISTOGRAMS -----------------------------------------------------
      
      # EDITOR PLOT
      editPlotServer(
        "editor_plot",
        ID_comp_trans = reactive({values$ID_comp_trans}),
        NIL_comp_trans = reactive({values$NIL_comp_trans}),
        parent = reactive({values$parent}),
        opts = editor_opts,
        xchan = reactive({values$xchan_select}),
        ychan = ychan,
        axes_trans = axes_trans,
        axes_limits = axes_limits,
        events = events,
        point_size = point_size,
        axes_text_size = axes_text_size,
        axes_label_text_size = axes_label_text_size,
        title_text_size = title_text_size,
        ...
      )
      
      # PLOTS TAB --------------------------------------------------------------
      
      # NAME - UNSTAINED CONTROL
      plots_NIL <- cytoSelectServer(
        "plots_select_unst", 
        data = reactive(x),
        choices = "None",
        selected = reactive(values$NIL_select)
      )
      
      observe({
        if(!.empty(plots_NIL(), null = TRUE)) {
          values$NIL_select <- plots_NIL()
        }
      })
      
      # NAME - SELECTED CONTROL
      plots_ID <- cytoSelectServer(
        "plots_select",
        data = reactive(x),
        selected = reactive(values$ID_select)
      )
      
      observe({
        if(!.empty(plots_ID(), null = TRUE)) {
          values$ID_select <- plots_ID()
        }
      })
      
      # PARENT DROPDOWN
      plots_parent <- nodeSelectServer(
        "plots_node_select",
        data = reactive(x),
        selected = reactive(values$parent)
      )
      
      # PARENT MANUALLY SELECTED
      observe({
        if(!.empty(plots_parent(), null = TRUE)) {
          values$parent <- plots_parent()
        }
      })
      
      # CHANNEL
      plots_xchan <- channelSelectServer(
        "plots_xchannel",
        data = reactive(x),
        selected = reactive(values$xchan_select)
      )
      
      observe({
        if(!.empty(plots_xchan(), null = TRUE)) {
          values$xchan_select <- plots_xchan()
        }
      })
      
      # PLOTS OPTIONS
      plots_opts <- optionsServer(
        "plots_options",
        selected = reactive({values$plots_opts_select})
      )
      
      # UPDATE PLOT OPTIONS
      observe({
        values$plots_opts_select <- plots_opts()
      })
      
      # COMPENSATION PLOTS
      compPlotServer(
        "plots",
        ID_comp_trans = reactive({values$ID_comp_trans}),
        NIL_comp_trans = reactive({values$NIL_comp_trans}),
        opts = reactive({values$plots_opts_select}),
        xchan = reactive({values$xchan_select}),
        channels = reactive({channels}),
        spillover = reactive({values$spill}),
        axes_trans = axes_trans,
        axes_limits = axes_limits,
        events = events,
        point_size = point_size,
        axes_text_size = axes_text_size,
        axes_label_text_size = axes_label_size,
        title_text_size = title_text_size,
        ...
      )
      
      # SAVE & EXPORT ----------------------------------------------------------
      
      # SAVE SPILLOVER MATRIX
      spill_save <- spillSaveServer(
        "editor_save",
        spill = reactive(values$spill),
        save_as = save_as
      )
      
      # RETURN
      onStop(function() {
        stopApp(read_from_csv(save_as))
      })
    }
  )
  
  # RUN SHINY APPLICATION
  if(viewer) {
    runApp(app,
           launch.browser = paneViewer(),
           quiet = TRUE)
  } else {
    runApp(app,
           quiet = TRUE)
  }
  
}

#' @noRd
cytoSelectUI <- function(id,
                         label = NULL,
                         ...) {
  selectInput(
    NS(id, "select"),
    label = label,
    choices = NULL,
    ...)
}

#' @noRd
cytoSelectServer <- function(id,
                             data = reactive(NULL),
                             choices = NULL,
                             selected = reactive(NULL)) {
  
  moduleServer(id, function(input,
                            output, session){
    
    # NAMESPACE
    ns <- session$ns
    
    # VALUES
    values <- reactiveValues(
      select = NULL
    )
    
    # UPDATE UI OPTIONS
    observe({
      updateSelectInput(
        session,
        "select",
        choices = c(cyto_names(data()), choices),
        selected = selected()
      )
    })
    
    observeEvent(input$select, {
      values$select <- input$select
    })
    
    return(reactive({values$select}))
  })
  
}

#' @noRd
channelSelectUI <- function(id,
                            label = NULL,
                            ...) {
  
  fluidRow(
    column(
      10,
      style = "padding-right: 5px;",
      selectInput(
        NS(id, "select"),
        label = label,
        choices = NULL
      )
    ),
    column(
      1,
      style = paste(
        "padding-left: 4px;",
        "padding-right: 0px;",
        "margin-top: 26px;",
        "margin-left: 0px;",
        "margin-right: 0px;"
      ),
      actionButton(
        NS(id, "down"),
        label = NULL,
        icon = icon("arrow-down"),
        style = paste0(
          "background-color: #FF0000;",
          "padding: 4px;"
        )
      )
    ),
    column(
      1,
      style = paste(
        "padding-left: 0px;",
        "padding-right: 2px;",
        "margin-top: 26px;",
        "margin-left: 0px;",
        "margin-right: 0px;"
      ),
      actionButton(
        NS(id, "up"),
        label = NULL,
        icon = icon("arrow-up"),
        style = paste0(
          "background-color: #33CC33;",
          "padding: 4px"
        )
      )
    )
  )
  
}

#' @noRd
channelSelectServer <- function(id,
                                data = reactive(NULL),
                                selected = reactive(NULL),
                                exclude = reactive(NULL),
                                ...) {
  
  moduleServer(id, function(input, 
                            output, 
                            session) {
    
    # NAMESPACE
    ns <- session$ns
    
    # VALUES
    values <- reactiveValues(
      channels = NULL,
      channels_excl = NULL,
      select = NULL
    )
    
    # CHANNELS
    observe({
      if(!.empty(data(), null = TRUE)) {
        values$channels <- unname(cyto_fluor_channels(data()))
        values$channels_excl <- unname(cyto_fluor_channels(data()))
      }
    })
    
    # EXCLUDE CHANNELS
    observe({
      if(!.empty(values$channels, null = TRUE) & 
         !.empty(exclude(), null = TRUE)) {
        values$channels_excl <- values$channels[
          -match(exclude(), values$channels)
        ]
      }
    })
    
    # UPDATE UI OPTIONS
    observe({
      updateSelectInput(
        session,
        "select",
        choices = values$channels_excl
      )
    })
    
    # STORE SELECTION
    observeEvent(input$select, {
      values$select <- input$select
    })
    
    observe({
      if(!is.null(selected())) {
        if(is.na(selected())) {
          select <- values$channels[1]
        } else {
          select <- selected()
        }
        updateSelectInput(
          session,
          "select",
          selected = selected()
        )
      }
    })
    
    # CHANNEL UP 
    observeEvent(input$up, {
      ind <- match(values$select, values$channels_excl)
      if(ind < length(values$channels_excl)) {
        ind <- ind + 1
      } else {
        ind <- 1
      }
      updateSelectInput(
        session,
        "select",
        choices = values$channels_excl,
        selected = values$channels_excl[ind]
      )
    })
    
    # CHANNEL DOWN
    observeEvent(input$down, {
      ind <- match(values$select, values$channels_excl)
      if(ind > 1) {
        ind <- ind - 1
      } else {
        ind <- length(values$channels_excl)
      }
      updateSelectInput(
        session,
        "select",
        choices = values$channels_excl,
        selected = values$channels_excl[ind]
      )
    })
    
    return(
      reactive({
        values$select
      })
    )
    
  })
  
}

#' @noRd
nodeSelectUI <- function(id,
                         label = NULL,
                         ...) {
  
  selectInput(
    NS(id, "select"),
    label = label,
    choices = "root",
    ...
  )
  
}

#' @noRd
nodeSelectServer <- function(id,
                             data = reactive(NULL),
                             choices = NULL,
                             selected = reactive(NULL)) {
  
  moduleServer(
    id, 
    function(input, output, session) {
    
    # NAMESPACE
    ns <- session$ns
    
    # VALUES
    values <- reactiveValues(
      select = NULL
    )
    
    # UPDATE UI OPTIONS
    observe({
      if(cyto_class(data(), "GatingSet")) {
        updateSelectInput(
          session,
          "select",
          choices = cyto_nodes(data(), path = "auto"),
          selected = selected()
        )
      }
    })
    
    observeEvent(input$select, {
      values$select <- input$select
    })
    
    return(reactive({values$select}))
  })
  
}

#' @noRd
spillEditUI <- function(id,
                        height = "300px",
                        ...) {
  # SPILLOVER MATRIX
  rHandsontableOutput(NS(id, "spill"),
                      height = height,
                      width = "99%",
                      ...)
  
}

#' @noRd
spillEditServer <- function(id,
                            data = reactive(NULL),
                            spill = reactive(NULL),
                            xchan = reactive(NULL),
                            ychan = reactive(NULL),
                            ...) {
  
  moduleServer(id, function(input, 
                            output, 
                            session){
    
    # VALUES
    values <- reactiveValues(
      spill = NULL,
      index = c(-1, -1),
      rows = NULL,
      cols = NULL
    )
    
    # PREPARE SPILLOVER MATRIX
    spill_mat <- reactive({
      if(.empty(spill(), null = TRUE)) {
        # EXTRACT SPILLOVER MATRIX - (FIRST SAMPLE)
        if(!.empty(data(), null = TRUE)) {
          sp <- cyto_spillover_extract(
            cyto_data_extract(
              data(),
              parent = "root",
              copy = FALSE
            )[[1]]
          )[[1]]
        } else {
          sp <- NULL
        }
      } else {
        # SPILLOVER CSV FILENAME
        if(is.character(spill())) {
          sp <- read_from_csv(spill(),
                              data.table = FALSE)
          # SPILLOVER MATRIX SUPPLIED
        } else {
          sp <- spill()
        }
      }
      # FORMAT SPILLOVER MATRIX
      if(!.empty(sp, null = TRUE)) {
        # PERCENT
        sp <- sp * 100
        # SQUARE MATRIX
        rownames(sp) <- colnames(sp)
      }
      return(sp)
    })
    # LOAD SPILLOVER MATRIX
    observe({
      values$spill <- spill_mat()
    })
    # HIGHLIGHT
    observe({
      if(!.empty(xchan(), null = TRUE) & !.empty(ychan(), null = TRUE)) {
        if(!.empty(values$spill, null = TRUE)) {
          values$index <- c(match(xchan(), rownames(values$spill)) - 1,
                            match(ychan(), colnames(values$spill)) - 1)
        }
      }
    })
    # RHANDSONTABLE
    output$spill <- renderRHandsontable({
      # SPILLOVER MATRIX
      if(!.empty(values$spill, null = TRUE)) {
        rhandsontable(
          values$spill,
          rowHeaderWidth = 105,
          readOnly = FALSE,
          manualColumnResize = TRUE,
          row_highlight = values$index[1],
          col_highlight = values$index[2]
        ) %>%
          hot_cols(
            type = "numeric",
            colWidths = 105,
            format = "0.000",
            halign = "htCenter",
            renderer = paste0(
              "function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);

                    if(instance.params) {
                      hrows = instance.params.row_highlight;
                      hrows = hrows instanceof Array ? hrows : [hrows];
                      hcols = instance.params.col_highlight;
                      hcols = hcols instanceof Array ? hcols : [hcols];
              
                      if (hcols.includes(col) && hrows.includes(row)) {
                        td.style.border = 'solid';
                        td.style.borderWidth = '3px';
                        td.style.borderColor = 'black';
                      }
                    }
              
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
            )
          ) %>%
          hot_rows(rowHeights = 20)
      }
    })
    # STORE SPILLOVER MATRIX EDITS
    observeEvent(input$spill, {
      sp <- hot_to_r(input$spill)
      rownames(sp) <- colnames(sp)
      if(any(is.na(sp))) {
        if(is.null(values$spill)) {
          sp[is.na(sp)] <- 0
        } else {
          # NO RE-RENDER
          sp[is.na(sp)] <- values$spill[is.na(sp)]
        }
      }
      values$spill <- sp
    })
    # RETURN EDITED SPILLOVER MATRIX
    return(
      reactive({
        if(is.null(values$spill)) {
          return(values$spill)
        } else {
          return(values$spill/100)
        }
      })
    )
  })
  
}

#' @noRd
optionsUI <- function(id,
                      ...) {
  
  checkboxGroupInput(
    NS(id, "options"),
    ...
  )
  
}

#' @noRd
optionsServer <- function(id, 
                          selected = reactive(NULL),
                          ...) {
  
  moduleServer(id, function(input,
                            output, session){
    
    # VALUES
    values <- reactiveValues(
      options = NULL
    )
    # UPDATE OPTIONS
    observeEvent(input$options, {
      values$options <- input$options
    })
    observe({
      if(!.empty(selected(), null = TRUE)) {
        updateCheckboxGroupInput(
          session,
          "options",
          selected = selected()
        )
      }
    })
    # RETURN OPTIONS
    return(reactive({values$options}))
    
  })
  
}

#' @noRd
spillSaveUI <- function(id, 
                        ...){
  
  actionButton(
    NS(id, "save"),
    "Save",
    ...
  )
  
}

#' @noRd
spillSaveServer <- function(id,
                            spill = reactive(NULL),
                            save_as = NULL,
                            ...) {
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    observe({
      input$save
      write_to_csv(spill(), save_as)
    })
    
    return(reactive({spill}))
    
  })
  
}

#' @noRd
editPlotUI <- function(id,
                       ...) {
  
  plotOutput(
    NS(id, "plot"),
    width = "70%",
    ...
  )
  
}

#' @noRd
editPlotServer <- function(id,
                           ID_comp_trans = reactive(NULL),
                           NIL_comp_trans = reactive(NULL),
                           opts = reactive(NULL),
                           xchan = reactive(NULL),
                           ychan = reactive(NULL),
                           axes_trans = NA,
                           axes_limits = "machine",
                           events = 2000,
                           point_size = 3,
                           axes_text_size = 1.7,
                           axes_label_text_size = 2,
                           title_text_size = 1.5,
                           ...) {
  
  moduleServer(id, function(input,
                            output,
                            session){

    # PLOTS
    output$plot <- renderPlot({
      # BYPASS PLOTS FOR MISSING DATA
      if((!.empty(ID_comp_trans(), null = TRUE) | 
          (!.empty(ID_comp_trans(), null = TRUE) &
           !.empty(NIL_comp_trans(), null = TRUE))) &
         (!.empty(xchan(), null = TRUE) & !.empty(ychan(), null = TRUE))) {
        # OVERLAY
        overlay <- FALSE
        if(!.empty(NIL_comp_trans(), null = TRUE) & 
           "overlay" %in% opts()){
          overlay <- TRUE
        }
        # CUSTOM LAYOUT
        cyto_plot_custom(
          layout = matrix(
            c(1, 1, 2, 2, 1, 1, 3, 3),
            byrow = TRUE,
            ncol = 4
          ),
          popup = FALSE
        )
        # SCATTERPLOT - NO PROGRESS BAR
        suppressPrint(
          cyto_plot(
            ID_comp_trans(),
            channels = c(xchan(), ychan()),
            overlay = if(overlay) {
              NIL_comp_trans()
            } else {
              NA
            },
            axes_trans = axes_trans,
            axes_limits = axes_limits,
            events = events,
            title = cyto_names(ID_comp_trans()),
            point_size = point_size,
            axes_text_size = axes_text_size,
            axes_label_text_size = axes_label_text_size,
            title_text_size = title_text_size,
            key_text_size = 2,
            key_title_text_size = 2,
            popup = FALSE,
            margins = c(4,5,6,6.5),
            ...
          )
        )
        # STORE PLOT LIMITS
        usr <- .par("usr")[[1]]
        # UNSTAINED MEDIAN
        if(!.empty(NIL_comp_trans(), null = TRUE) & "line" %in% opts()) {
          medFI <- cyto_apply(
            NIL_comp_trans(),
            "cyto_stat_median",
            input = "matrix",
            inverse = FALSE,
            copy = FALSE
          )
          abline(h = medFI[1, ychan()],
                 col = "red",
                 lwd = 2)
        }
        # MEDIAN TRACKER
        if("tracker" %in% opts()) {
          .cyto_median_tracker(
            ID_comp_trans(),
            c(
              xchan(),
              ychan()
            )
          )
        }
        # HISTOGRAMS - CORRECT & OTHER CHANNEL
        for(i in seq_len(2)) {
          # CHANNEL
          chan <- c(xchan(), ychan())[i]
          # COMPUTE LABEL TEXT LOCATIONS
          if(chan %in% names(axes_trans)) {
            # AXES_LIMITS FROM 2D PLOT
            if(i == 1) {
              lims <- usr[1:2]
            } else {
              lims <- usr[3:4]
            }
            # LABEL_TEXT_X ON LINEAR SCALE
            label_text_x <- .cyto_transform(
              min(lims) + 0.90 * diff(lims),
              trans = axes_trans,
              channel = chan,
              inverse = TRUE
            )
          # NO TRANSFORMERS - LINEAR SCALE
          } else {
            label_text_x <- min(lims) + 
              0.90 * diff(lims)
          }
          # CONSTRUCT PLOT - NO PROGRESS BAR
          suppressPrint(
            cyto_plot(
              if (overlay) {
                NIL_comp_trans()
              } else {
                ID_comp_trans()
              },
              channels = chan,
              overlay = if (overlay) {
                ID_comp_trans()
              } else {
                NA
              },
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              events = events,
              title = NA,
              ylab = "Density",
              hist_fill = if (overlay){
                c("grey70", "white")
              }else {
                "white" 
              },
              hist_fill_alpha = if (overlay){
                c(1, 0)
              }else {
                1
              },
              hist_line_col = if (overlay){
                c("black", "blue")
              } else {
                "blue"
              }, 
              hist_line_width = 2,
              axes_text_size = axes_text_size,
              axes_label_text_size = axes_label_text_size,
              title_text_size = title_text_size,
              popup = FALSE,
              label_text = if (chan == ychan()){
                if (overlay) {
                  c("MedFI", "MedFI")
                } else {
                  c("MedFI")
                }
              } else {
                NA
              },
              label_stat = if (chan == ychan()){
                if (overlay) {
                  c("median", "median")
                } else {
                  c("median")
                }
              } else {
                NA
              },
              label_text_x = if (chan == ychan()) {
                if(overlay) {
                  rep(
                    label_text_x,
                    2
                  )
                } else {
                  label_text_x
                }
              } else {
                NA
              },
              label_text_y = if (chan == ychan()) {
                if(overlay) {
                  c(80, 30)
                } else {
                  c(50)
                }
              } else {
                NA
              },
              label_text_col = if (chan == ychan()){
                if(overlay){
                  c("grey40", "blue")
                } else {
                  c("blue")
                }
              } else {
                "black"
              },
              label_text_size = 1.1,
              margins = c(5,6,2,4),
              ...
            )
          )
        }
        # RESET - CYTO_PLOT_COMPLETE
        cyto_plot_complete()
      }
    })
  })
}

#' @noRd
compPlotUI <- function(id,
                       ...) {
  
  plotOutput(
    NS(id, "cyto_plot_comp"),
    height = "800px",
    ...
  )
  
}

#' @importFrom flowWorkspace cytoset
#' @noRd
compPlotServer <- function(id,
                           ID_comp_trans = reactive(NULL),
                           NIL_comp_trans = reactive(NULL),
                           opts = reactive(NULL),
                           xchan = reactive(NULL),
                           channels = reactive(NULL),
                           spillover = reactive(NULL),
                           axes_trans = NA,
                           axes_limits = "machine",
                           events = 2000,
                           point_size = 3,
                           axes_text_size = 1.7,
                           axes_label_text_size = 2,
                           title_text_size = 1.5,
                           ...) {
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    # VALUES
    values <- reactiveValues(
      layout = NULL
    )
    
    # CUSTOM LAYOUT
    observe({
      values$layout <- c(ceiling(length(channels())/4), 4)
    })
    
    # PLOTS
    output$cyto_plot_comp <- renderPlot({
      # BYPASS PLOTS FOR MISSING DATA
      if((!.empty(ID_comp_trans(), null = TRUE) | 
          (!.empty(ID_comp_trans(), null = TRUE) &
           !.empty(NIL_comp_trans(), null = TRUE))) &
         (!.empty(xchan(), null = TRUE) & !.empty(channels(), null = TRUE))) {
        # PREPARE DATA
        if(!.empty(ID_comp_trans(), null = TRUE)) {
          # COMBINE UNSTAINED
          if(!.empty(ID_comp_trans(), null = TRUE)) {
            cf_list <- structure(
              list(
                ID_comp_trans()[[1]],
                NIL_comp_trans()[[1]]
              ),
              names = c(
                cyto_names(ID_comp_trans()),
                cyto_names(NIL_comp_trans())
              )
            )
            cf_list <- cf_list[!LAPPLY(cf_list, "is.null")]
            cs <- cytoset(
              cf_list
            )
          # STAINED ONLY
          } else {
            cs <- ID_comp_trans()
          }
        }
        # UPDATE EXPERIMENT DETAILS
        pd <- cyto_details(cs)
        pd$channel <- c(xchan(), "unstained")[seq_along(cs)]
        # CYTO_PLOT_COMPENSATION
        suppressPrint(
          suppressWarnings(
            cyto_plot_compensation(
              cs,
              channels = channels(),
              channel_match = pd,
              overlay = if(!any(c("unstained", "compensated") %in% opts())){
                "none"
              } else {
                opts()[opts() %in% c("unstained", "compensated")]
              },
              spillover = spillover(),
              compensated = TRUE,
              axes_trans = axes_trans,
              axes_limits = axes_limits,
              events = events,
              point_size = point_size,
              point_col = if(all(c("unstained", "compensated") %in% opts())) {
                c("magenta", "blue", "grey40")
              } else if("unstained" %in% opts()) {
                c("magenta", "grey40")
              } else if("compensated" %in% opts()) {
                c("magenta", "blue")
              } else {
                c("magenta")
              },
              hist_fill = if(all(c("unstained", "compensated") %in% opts())) {
                c("magenta", "blue", "grey40")
              } else if("unstained" %in% opts()) {
                c("magenta", "grey40")
              } else if("compensated" %in% opts()) {
                c("magenta", "blue")
              } else {
                c("magenta")
              },
              lines = if("models" %in% opts()) {
                TRUE
              } else {
                FALSE
              },
              text = TRUE,
              text_size = 1.5,
              axes_text_size = 1.7,
              axes_label_text_size = 2,
              title_text_size = 2,
              layout = values$layout,
              popup = FALSE,
              ...
            )
          )
        )
      }
    },
    height = reactive({250 * values$layout[1]})
    )
  })
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
  raw_data <- cyto_data_extract(
    x, 
    format = "matrix",
    channels = channels,
    copy = FALSE
  )[[1]][[1]]
  
  # BYPASS EMPTY SAMPLES
  if(nrow(raw_data) > 0) {
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
    raw_data_chunks <- lapply(
      seq_len(n - 1), 
      function(z) {
        if (z < n - 1) {
          raw_data[raw_data[, channels[1]] >= chunks[z] &
                   raw_data[, channels[1]] < chunks[z + 1], ]
        } else {
          raw_data[raw_data[, channels[1]] >= chunks[z] &
                   raw_data[, channels[1]] <= chunks[z + 1], ]
        }
      }
    )
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
    loessMod <- suppressWarnings(
      loess(
        vals_y ~ vals_x,
        data = medians,
        span = 0.9
      )
    )
    loessMod <- predict(loessMod)
    # ADD LINE TO PLOT
    lines(
      medians[, channels[1]],
      loessMod,
      col = "purple2",
      lwd = 3
    )
  }
  
}
