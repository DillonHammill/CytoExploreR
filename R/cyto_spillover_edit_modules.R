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
#' @param parent name of the parent population to plot when a \code{GatingSet}
#'   object is supplied.
#' @param channel_match name of csv file matching the name of each sample to a
#'   fluorescent channel. The \code{channel_match} file must contain the columns
#'   "name" and "channel". The \code{channel_match} file not required to use
#'   \code{cyto_spillover_edit} but is used internally to automatically select
#'   channels associated with the selected samples.
#' @param spill name of a square spillover matrix csv file or spillover matrix
#'   to edit. Setting \code{spill} to NULL (the default) will result in
#'   extraction of the spillover matrix generated on the cytometer which is
#'   attached to the samples. Similarly, if the supplied spillover matrix csv
#'   file does not exist, the spillover matrix attached to the first sample will
#'   be used and the edited spillover matrix will be saved to the specified
#'   file.
#' @param spillover added for backwards compatibility with older versions of
#'   CytoExploreR, this performs the same operations as \code{spill} but will be
#'   eventually phased out in favour of the more concise \code{spill} argument.
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
#' @importFrom shinythemes shinytheme
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
                                parent = "root",
                                channel_match = NULL,
                                spill = NULL,
                                spillover = NULL,
                                save_as = NULL,
                                axes_trans = NA,
                                axes_limits ="machine",
                                events = 2000,
                                point_size = 4,
                                axes_text_size = 1.7,
                                axes_label_text_size = 2,
                                title_text_size = 2,
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
  x <- cyto_copy(x)
  
  # CHANNELS
  channels <- cyto_fluor_channels(x)
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # SAMPLE NAMES
  nms <- cyto_names(x)
  
  # CHANNEL_MATCH - EXPERIMENT DETAILS
  if(any(grepl("channel", colnames(pd), ignore.case = TRUE))) {
    # MARKERS TO CHANNELS
    ind <- which(grepl("channel", colnames(pd), ignore.case = TRUE))
    pd[, "channel"] <- LAPPLY(pd[, ind], function(z) {
      if (!grepl("unstained", z, ignore.case = TRUE)) {
        return(cyto_channels_extract(x, z))
      } else {
        return(z)
      }
    })
  # CHANNEL_MATCH - MANUALLY SUPPLIED
  } else {
    if(!is.null(channel_match)) {
      # CHANNEL MATCH - MATRIX/DATA.FRAME
      if(cyto_class(channel_match, c("data.frame",
                                     "matrix",
                                     "tibble"))) {
        if(!all(c("name", "channel") %in% colnames(channel_match))) {
          stop(
            "'channel_match' must contain columns 'name' and 'channel'!"
          )
        }
      # CHANNEL_MATCH FILENAME
      } else {
        channel_match <- read_from_csv(channel_match)
      }
      # UPDATE EXPERIMENT DETAILS - CHANNEL SELECTION
      pd$channel <- rep(NA, nrow(pd))
      lapply(nms, function(z) {
        if (z %in% channel_match$name) {
          chn <- channel_match$channel[match(z, channel_match$name)]
          pd[match(z, pd$name), "channel"] <<- chn
        }
      })
    # CHANNEL_MATCH - NOT SUPPLIED
    } else {
      # LEAVE ALONE
      # pd$channel = NULL
      # pd$channel[c(1:4)] = NULL
    }
  }
  
  # PARENT  
  if (is.null(parent)) {
    if (!"parent" %in% colnames(pd)) {
      if (!is.null(channel_match)) {
        if ("parent" %in% colnames(channel_match)) {
          pd[match(channel_match[, "name"],rownames(pd)), "parent"] <- 
            channel_match[, "parent"]
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
  cyto_details(x) <- pd
  
  # UNCOMPENSATED/TRANSFORMED & LINEAR GATINGSET
  if(cyto_class(x, "GatingSet")) {
    # TRANSFORMERS
    axes_trans <- cyto_transformers_extract(x)
    cs <- cyto_data_extract(
      x, 
      parent = "root",
      copy = FALSE
    )[[1]]
    # COMPENSATION
    comp <- cyto_spillover_extract(x)
    # COMPENSATED GATINGSET
    if(!is.null(comp)) {
      # REVERSE TRANSFORMATIONS
      if(!.all_na(axes_trans)) {
        trans <- axes_trans[match_ind(channels, names(axes_trans))]
        trans <- cyto_transformers_combine(trans)
        cs <- cyto_transform(cs, 
                             trans = trans, 
                             inverse = TRUE, 
                             plot = FALSE,
                             quiet = TRUE)
      }
      # REMOVE COMPENSATION
      cs <- cyto_compensate(cs, 
                            spillover = comp,
                            remove = TRUE)
      # GET DEFAULT TRANSFORMERS
      if (.all_na(axes_trans)) {
        axes_trans <- cyto_transformers_define(x,
                                               channels = channels,
                                               type = "biex",
                                               plot = FALSE)
      }
      # RE-APPLY TRANSFORMATIONS
      trans <- axes_trans[match_ind(channels, names(axes_trans))]
      trans <- cyto_transformers_combine(trans)
      cs <- cyto_transform(cs,
                           trans = trans,
                           plot = FALSE,
                           quiet = TRUE)
      # REPLACE DATA
      gs_cyto_data(x) <- cs
    # UNCOMPENSATED GATINGSET
    } else {
      # UNTRANSFOMED GATINGSET
      if (.all_na(axes_trans)) {
        # GET DEFAULT TRANSFORMERS
        axes_trans <- cyto_transformers_define(x,
                                               channels = channels,
                                               type = "biex",
                                               plot = FALSE)
        # RE-APPLY TRANSFORMATIONS
        trans <- axes_trans[match_ind(channels, names(axes_trans))]
        trans <- cyto_transformers_combine(trans)
        cs <- cyto_transform(cs,
                             trans = trans,
                             plot = FALSE,
                             quiet = TRUE)
        # REPLACE DATA
        gs_cyto_data(x) <- cs
      }
    }
    # LINEAR UNTRANSFORMED GATINGSET
    if (!.all_na(axes_trans)) {
      # GS TRANSFORMED & UNCOMPENSATED
      x_linear <- cyto_copy(x)
      cs <- cyto_data_extract(x_linear, "root")[[1]]
      # INVERSE TRANSFORMATIONS
      trans <- axes_trans[match_ind(channels, names(axes_trans))]
      trans <- cyto_transformers_combine(trans)
      cs <- cyto_transform(cs,
                           trans = trans,
                           inverse = TRUE,
                           plot = FALSE,
                           quiet = TRUE)
      # REPLACE DATA
      gs_cyto_data(x_linear) <- cs
    }
  # UNCOMPENSATED/TRANSFORMED & LINEAR CYTOSET
  } else {
    # LINEAR DATA
    if(.all_na(axes_trans)) {
      x_linear <- cyto_copy(x)
      axes_trans <- cyto_transformers_define(x,
                                             channels = channels,
                                             type = "biex",
                                             plot = FALSE)
      x <- cyto_transform(x,
                          trans = axes_trans,
                          plot = FALSE,
                          quiet = TRUE)
    # TRANSFORMED DATA
    } else {
      trans <- axes_trans[match_ind(channels, names(axes_trans))]
      trans <- cyto_transformers_combine(trans)
      x_linear <- cyto_transform(cyto_copy(x),
                                 trans = trans,
                                 inverse = TRUE,
                                 plot = FALSE,
                                 quiet = TRUE)
    }
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
  
  # SPILL - BACKWARDS COMPATIBILITY
  if(!is.null(spillover)) {
    spill <- spillover
    message(
      paste0(
        "'spillover' is now deprecated in favour of the more concise 'spill' ",
        "argument. Please pass your spillover matrix to 'spill' instead."
      )
   )
  }
  
  # SAVE_AS
  if(is.character(spill)) {
      save_as <- spill
  } else {
    save_as <- paste0(
      format(Sys.Date(), "%d%m%y"),
      "-", "Spillover-Matrix.csv"
    )
  }
  
  # SHINY APPLICATION ----------------------------------------------------------
  
  # APPLICATION
  app <- shinyApp(
    # USER INTERFACE
    ui <- fluidPage(
      theme = shinytheme("yeti"),
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
              dataSelectUI("editor_select_unst",
                           label = "Select unstained sample:"),
              dataSelectUI("editor_select",
                           label = "Select sample:"),
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
              dataSelectUI("plots_select_unst",
                           label = "Select unstained sample:"),
              dataSelectUI("plots_select",
                           label = "Select sample:"),
              channelSelectUI("plots_xchannel",
                              label = "Channel:"),
              optionsUI("plots_options",
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
                        ))
            ),
            mainPanel(
              width = 9,
              # compPlotUI(
              #   "plots"
              # )
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
        spill = NULL,
        NIL_select = NIL_select,
        NIL_comp_trans = NULL,
        ID_select = ID_select,
        ID_comp_trans = NULL,
        xchan_select = xchan_select,
        editor_opts_select = editor_opts_select,
        plots_opts_select = plots_opts_select
      )
      
      # UNSTAINED CONTROL ------------------------------------------------------
      
      # NAME - UNSTAINED CONTROL
      NIL <- dataSelectServer(
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
          # TRANSFORMERS
          trans <- axes_trans[match_ind(channels, names(axes_trans))]
          trans <- cyto_transformers_combine(trans)
          # COMPENSATE & TRANSFORM
          values$NIL_comp_trans <- cyto_transform(
            if(!is.null(values$spill)){
              cyto_compensate(
                NIL_linear,
                values$spill
              )
            } else {
              NIL_linear
            },
            trans = trans,
            plot = FALSE,
            quiet = TRUE
          )
        } else {
          values$NIL_comp_trans <- NULL
        }
      })
      
      # SELECTED CONTROL -------------------------------------------------------
      
      # NAME - SELECTED CONTROL
      ID <- dataSelectServer(
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
        if(!.empty(values$ID_select, null = TRUE)) {
          # PARENT
          values$parent <- cyto_details(x_linear[values$ID_select])[, "parent"]
          # X CHANNEL
          values$xchan_select <- cyto_details(x_linear)$channel[
            cyto_match(x_linear, values$ID_select)
          ]
          # LINEAR DATA
          ID_linear <- cyto_data_extract(
            x_linear,
            select = values$ID_select,
            parent = values$parent,
            copy = TRUE
          )[[1]]
          # TRANSFORMERS
          trans <- axes_trans[match_ind(channels, names(axes_trans))]
          trans <- cyto_transformers_combine(trans)
          # COMPENSATE & TRANSFORM
          values$ID_comp_trans <- cyto_transform(
            if(!is.null(values$spill)){
              cyto_compensate(
                ID_linear,
                values$spill
              )
            } else {
              ID_linear
            },
            trans = trans,
            plot = FALSE,
            quiet = TRUE
          )
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
        spill = reactive(spill),
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
        xchan = xchan,
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
      plots_NIL <- dataSelectServer(
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
      plots_ID <- dataSelectServer(
        "plots_select",
        data = reactive(x),
        selected = reactive(values$ID_select)
      )
      
      observe({
        if(!.empty(plots_ID(), null = TRUE)) {
          values$ID_select <- plots_ID()
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
        selected = reactive(values$plots_opts_select)
      )
      
      # COMPENSATION PLOTS
      # compPlotServer(
      #   
      # )
      
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
dataSelectUI <- function(id,
                         label = NULL,
                         ...) {
  selectInput(
    NS(id, "select"),
    label = label,
    choices = NULL,
    ...)
}

#' @noRd
dataSelectServer <- function(id,
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
        values$channels <- cyto_fluor_channels(data())
        values$channels_excl <- cyto_fluor_channels(data())
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
      updateSelectInput(
        session,
        "select",
        selected = selected()
      )
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
      write_to_csv(spill(),
                   save_as)
    })
    
    return(reactive({spill}))
    
  })
  
}

#' @noRd
editPlotUI <- function(id,
                       ...) {
  
  plotOutput(
    NS(id, "plot"),
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
                           title_text_size = 2,
                           ...) {
  
  moduleServer(id, function(input,
                            output,
                            session){

    # PLOTS
    output$plot <- renderImage({
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
        # SAVE IMAGE
        temp <- paste0(tempdir(),
                       .Platform$file.sep,
                       "cyto_spillover_edit.png")
        cyto_plot_save(
          temp,
          height = 6,
          width = 12,
          res = 300
        )
        # CUSTOM LAYOUT
        cyto_plot_custom(
          layout = matrix(
            c(1, 1, 2, 2, 1, 1, 3, 3),
            byrow = TRUE,
            ncol = 4
          ),
          popup = FALSE
        )
        # SCATTERPLOT
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
          ...
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
          # CONSTRCUT PLOT
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
            ...
          )
        }
        # COMPLETE
        cyto_plot_complete()
        # RENDER IMAGE
        list(src = temp,
             height = 400)
      } else {
        # SAVE BLANK IMAGE
        temp <- paste0(tempdir(),
                       .Platform$file.sep,
                       "cyto_spillover_edit.png")
        cyto_plot_save(temp,
                       height = 6,
                       width = 12,
                       res = 300)
        cyto_plot_complete()
        list(src = temp,
             height = 400)
      }
    }, deleteFile = TRUE)
  })
}

#' @noRd
compPlotUI <- function(id,
                       ...) {
  
  plotOutput(
    NS(id, "cyto_plot_comp"),
    ...
  )
  
}

#' @noRd
compPlotServer <- function(id,
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
                           title_text_size = 2,
                           ...) {
  
  moduleServer(id, function(input,
                            output,
                            session){
    
    # PLOTS
    output$cyto_plot_comp <- renderImage({
      # BYPASS PLOTS FOR MISSING DATA
      if((!.empty(ID_comp_trans(), null = TRUE) | 
          (!.empty(ID_comp_trans(), null = TRUE) &
           !.empty(NIL_comp_trans(), null = TRUE))) &
         (!.empty(xchan(), null = TRUE) & !.empty(ychan(), null = TRUE))) {
        # SAVE IMAGE
        temp <- paste0(tempdir(),
                       .Platform$file.sep,
                       "cyto_spillover_edit.png")
        cyto_plot_save(
          temp,
          height = 6,
          width = 12,
          res = 300
        )
        # CYTO_PLOT_COMPENSATION
        
        # COMPLETE
        cyto_plot_complete()
        # RENDER IMAGE
        list(src = temp,
             height = 400)
      } else {
        # SAVE BLANK IMAGE
        temp <- paste0(tempdir(),
                       .Platform$file.sep,
                       "cyto_spillover_edit.png")
        cyto_plot_save(temp,
                       height = 6,
                       width = 12,
                       res = 300)
        cyto_plot_complete()
        list(src = temp,
             height = 400)
      }
    }, deleteFile = TRUE)
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
  raw_data <- cyto_data_extract(x, 
                                format = "matrix",
                                copy = FALSE)[[1]][[1]]
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
  loessMod <- suppressWarnings(
    loess(vals_y ~ vals_x,
          data = medians,
          span = 0.9)
  )
  loessMod <- predict(loessMod)
  # ADD LINE TO PLOT
  lines(medians[, channels[1]],
        loessMod,
        col = "purple2",
        lwd = 3
  )
}
