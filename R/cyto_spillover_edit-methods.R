#' Edit Spillover Matrix
#'
#' Edit spillover matrices in real-time using a shiny interface.
#'
#' \code{cyto_spillover_edit} provides an interactive shiny interface for
#' editing fluorescent spillover matrices. \code{cyto_spillover_edit} takes on
#' either a \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' untransformed single stain compensation controls and a universal unstained
#' control. It is recommended that samples be pre-gated based on FSC and SSC
#' parameters to obtain a homogeneous population for calculation of fluorescent
#' spillover. Users begin by selecting the unstained control and a stained
#' control from dropdown menus of sample names. \code{cyto_spillover_edit}
#' leverages \code{cyto_plot} to plot the stained sample and overlay the
#' unstained control in black. Users should then select the channel associated
#' with the selected control on the \code{x axis} and go through all other
#' channels on the \code{y axis}. The displayed spillover matrix is extracted
#' directly from the \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spillover argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time. To aid in
#' selection of an appropriate spillover value, the median fluorescent intensity
#' of the unstained control is indicated by a red line and median fluorescent
#' intensity of the stained control is tracked with a purple line. These
#' features can be turned off by de-selecting the check boxes. Changes to the
#' spillover matrix are automatically saved to a csv file called
#' \code{"Spillover-Matrix.csv"} in the case where the \code{spillover} is not
#' specified or to the same name as the specified \code{spillover}.
#' \code{cyto_spillover_edit} has methods for both
#' \code{\link[flowCore:flowSet-class]{flowSet}} and
#' \code{\link[flowCore:flowSet-class]{flowSet}} objects refer to their
#' respective help pages for more information.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_spillover_edit,flowSet-method}}
#' @seealso \code{\link{cyto_spillover_edit,GatingSet-method}}
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#'
#' @export
cyto_spillover_edit <- function(x, ...) {
  UseMethod("cyto_spillover_edit")
}

#' Edit Spillover Matrix - flowSet Method
#'
#' Edit spillover matrices in real-time using a shiny interface.
#'
#' \code{cyto_spillover_edit} provides an interactive shiny interface for
#' editing fluorescent spillover matrices. \code{cyto_spillover_edit} takes on
#' either a \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' untransformed single stain compensation controls and a universal unstained
#' control. It is recommended that samples be pre-gated based on FSC and SSC
#' parameters to obtain a homogeneous population for calculation of fluorescent
#' spillover. Users begin by selecting the unstained control and a stained
#' control from dropdown menus of sample names. \code{cyto_spillover_edit}
#' leverages \code{cyto_plot} to plot the stained sample and overlay the
#' unstained control in black. Users should then select the channel associated
#' with the selected control on the \code{x axis} and go through all other
#' channels on the \code{y axis}. The displayed spillover matrix is extracted
#' directly from the \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spillover argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time. To aid in
#' selection of an appropriate spillover value, the median fluorescent intensity
#' of the unstained control is indicated by a red line and median fluorescent
#' intensity of the stained control is tracked with a purple line. These
#' features can be turned off by de-selecting the check boxes. Changes to the
#' spillover matrix are automatically saved to a csv file called
#' \code{"Spillover-Matrix.csv"} in the case where the \code{spillover} is not
#' specified or to the same name as the specified \code{spillover}.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param channel_match name of .csv file containing the names of the samples in
#'   a column called "name" and their matching channel in a column called
#'   "channel". \code{cyto_spillover_edit} will the guide you through the
#'   channel selection process and generate a channel match file called
#'   "Compensation-Channels.csv" automatically. If you already have a complete
#'   channel_match and would like to bypass the channel selection process,
#'   simply pass the name of the channel_match to this argument (e.g.
#'   "Compensation-Channels.csv").
#' @param spillover name of spillover matrix csv file including .csv file
#'   extension to use as a starting point for editing. If \code{spillover} is
#'   not supplied the spillover matrix will be extracted directly from the
#'   \code{\link[flowCore:flowSet-class]{flowSet}} and the edited matrix saved
#'   to \code{"Spillover-Matrix.csv"}.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines. Set to all events by default.
#' @param axes_trans object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the supplied flowSet.
#' @param ... additional arguments passed to
#'   \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return save edited spillover matrix to .csv file named
#'   "Spillover-matrix.csv" or spillover.
#'
#' @importFrom flowWorkspace sampleNames pData
#' @importFrom flowCore compensate fsApply sampleFilter exprs Subset each_col
#' @importFrom utils read.csv write.csv
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarPanel selectInput
#'   checkboxInput actionButton mainPanel plotOutput reactiveValues observe
#'   eventReactive renderPlot tabsetPanel tabPanel sidebarLayout fluidRow
#'   updateSelectInput
#' @importFrom rhandsontable rhandsontable rHandsontableOutput hot_to_r
#'   renderRHandsontable hot_cols hot_rows
#' @importFrom shinythemes shinytheme
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom stats median loess predict
#' @importFrom graphics lines layout
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_spillover_edit,GatingSet-method}}
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#'
#' # Load in compensation controls
#' fs <- Compensation
#' gs <- GatingSet(fs)
#'
#' # Gate single cells using gate_draw
#' gt <- Compensation_gatingTemplate
#' gating(gt, gs)
#'
#' # Channel match file
#' cmfile <- system.file("extdata",
#'   "Compensation-Channels.csv",
#'   package = "CytoRSuiteData"
#' )
#'
#' # Edit spillover matrix
#' cyto_spillover_edit(getData(gs, "Single Cells"),
#'   channel_match = cmfile
#' )
#' }
#' @export
cyto_spillover_edit.flowSet <- function(x,
                                        channel_match = NULL,
                                        spillover = NULL,
                                        display = 1,
                                        axes_trans = NULL, ...) {
  
  # Assign x to fs
  fs <- x
  
  # Extract sample names
  nms <- cyto_names(fs)
  
  # Extract fluorescent channels
  channels <- cyto_fluor_channels(fs)
  
  # Extract pData information
  pd <- cyto_details(fs)
  
  # Select a fluorescent channel for each compensation control
  if (is.null(channel_match)) {
    pd$channel <- paste(cyto_channel_select(fs))
    write.csv(pd, "Compensation-Channels.csv", row.names = FALSE)
  } else {
    if (inherits(channel_match, "data.frame") |
        inherits(channel_match, "matrix") |
        inherits(channel_match, "tibble")) {
      if (!all(c("name", "channel") %in% colnames(channel_match))) {
        stop("channel_match should contain columns 'name' and 'channel'.")
      }
      cm <- channel_match
    } else {
      if (getOption("CytoRSuite_wd_check") == TRUE) {
        if (file_wd_check(channel_match)) {
          cm <- read.csv(channel_match, header = TRUE, row.names = 1)
        } else {
          stop(paste(channel_match, "is not in this working directory."))
        }
      } else {
        cm <- read.csv(channel_match, header = TRUE, row.names = 1)
      }
    }
    chans <- cm$channel[match(cyto_names(fs), row.names(cm))]
    pd$channel <- paste(chans)
  }
  
  # Read in spillover matrix to object spill
  if (!is.null(spillover)) {
    if (inherits(spillover, "matrix") |
        inherits(spillover, "data.frame") |
        inherits(spillover, "tibble")) {
      spill <- spillover
      spillover <- "Spillover-Matrix.csv"
    } else {
      if (getOption("CytoRSuite_wd_check") == TRUE) {
        if (file_wd_check(spillover)) {
          spill <- read.csv(spillover, header = TRUE, row.names = 1)
          spill <- as.matrix(spill)
        } else {
          stop(paste(spillover, "is not in this working directory."))
        }
      } else {
        spill <- read.csv(spillover, header = TRUE, row.names = 1)
        spill <- as.matrix(spill)
      }
    }
  } else {
    spillover <- "Spillover-Matrix.csv"
    spill <- fs[[1]]@description$SPILL
  }
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  shinyApp(
    ui <- fluidPage(
      theme = shinytheme("yeti"),
      
      titlePanel("CytoRSuite Spillover Matrix Editor"),
      
      tabsetPanel(
        tabPanel("Editor",
                 fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel(
                     selectInput(
                       inputId = "Unstained",
                       label = "Select Unstained Control:",
                       choices = cyto_names(fs),
                       selected = pd$name[match("Unstained", pd$channel)]
                     ),
                     selectInput(
                       inputId = "flowFrame",
                       label = "Select Sample:",
                       choices = cyto_names(fs)
                     ),
                     selectInput(
                       inputId = "xchannel",
                       label = "X Axis:",
                       choices = channels
                     ),
                     selectInput(
                       inputId = "ychannel",
                       label = "Y Axis:",
                       choices = channels
                     ),
                     checkboxInput(
                       inputId = "NIL",
                       label = "Overlay Unstained Control",
                       value = TRUE
                     ),
                     checkboxInput(
                       inputId = "median",
                       label = "Unstained Control Median",
                       value = TRUE
                     ),
                     checkboxInput(
                       inputId = "trace",
                       label = "Median Tracker",
                       value = TRUE
                     ),
                     actionButton("saveBtn", "Save")
                   ),
                   mainPanel(
                     rHandsontableOutput("spillover", height = "300px"),
                     plotOutput("plot", height = "400px", width = "80%")
                   )
                 )
        ),
        tabPanel("Plots",
                 fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel(
                     selectInput(
                       inputId = "Unstained2",
                       label = "Select Unstained Control:",
                       choices = cyto_names(fs),
                       selected = pd$name[match("Unstained", pd$channel)]
                     ),
                     selectInput(
                       inputId = "flowFrame2",
                       label = "Select Sample:",
                       choices = cyto_names(fs)
                     ),
                     checkboxInput(
                       inputId = "NIL2",
                       label = "Overlay Unstained Control",
                       value = TRUE
                     ),
                     checkboxInput(
                       inputId = "Uncomp",
                       label = "Underlay Uncompensated Control",
                       value = TRUE
                     ),
                     checkboxInput(
                       inputId = "Comp",
                       label = "Overlay Compensated Control",
                       value = TRUE
                     )
                   ),
                   mainPanel(fluidRow(
                     plotOutput("plots", height = "700px", width = "100%")
                   ))
                 )
        )
      )
    ),
    server <- function(input, output, session) {
      values <- reactiveValues()
      
      observe({
        if (!is.null(input$spillover)) {
          spill <- hot_to_r(input$spillover)
          rownames(spill) <- colnames(spill)
          values$spill <- spill
        } else {
          values$spill <- spill * 100
        }
      })
      
      
      output$spillover <- renderRHandsontable({
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
      
      # Update xchannel selection based on selected control - spillover
      observe({
        fr <- input$flowFrame
        xchan <- pd$channel[match(input$flowFrame, pd$name)]
        updateSelectInput(session, "xchannel", selected = xchan)
        updateSelectInput(session, "flowFrame2", selected = fr)
      })
      
      # Apply compensation after each edit
      fs.comp <- eventReactive(values$spill, {
        fs <- compensate(fs, values$spill / 100)
        
        # Get trnsformed data
        fs <- cyto_transform(fs, axes_trans)
        
        return(fs)
      })
      
      output$plot <- renderPlot({
        layout(rbind(c(1, 1, 2, 2), c(1, 1, 3, 3)))
        
        # Plot
        if (input$NIL == FALSE) {
          cyto_plot(fs.comp()[[input$flowFrame]],
                    channels = c(input$xchannel, input$ychannel),
                    display = display,
                    axes_trans = axes_trans,
                    title = input$flowFrame,
                    point_size = 3,
                    limits = "machine",
                    axes_text_size = 1.7,
                    axes_label_text_size = 2,
                    title_text_size = 2, ...
          )
          
          if (input$median == TRUE) {
            medians <- fsApply(fs.comp(), each_col, "median")[
              input$Unstained,
              channels
              ]
            MFI <- data.frame("Channel" = channels, "Median" = medians)
            
            cutoff <- MFI[match(input$ychannel, MFI$Channel), ]
            
            abline(h = cutoff[2], col = "red", lwd = 2)
          }
          
          if (input$trace == TRUE) {
            cells <- exprs(fs.comp()[[input$flowFrame]])
            cells <- cells[order(cells[, input$xchannel]), ]
            cells <- as.data.frame(cells)
            
            n <- nrow(cells)
            splt <- seq(1, 20, 1)
            r <- ceiling(nrow(cells) / 20)
            cuts <- splt * r
            cells <- split(cells, cumsum(seq_len(nrow(cells)) %in% cuts))
            
            xmedians <- lapply(cells, function(x) median(x[, input$xchannel]))
            ymedians <- lapply(cells, function(x) median(x[, input$ychannel]))
            
            medians <- data.frame(unlist(xmedians), unlist(ymedians))
            colnames(medians) <- c(input$xchannel, input$ychannel)
            
            vals.y <- medians[, input$ychannel]
            vals.x <- medians[, input$xchannel]
            loessMod <- loess(vals.y ~ vals.x,
                              data = medians, span = 0.9
            )
            loessMod <- predict(loessMod)
            
            lines(medians[, input$xchannel],
                  loessMod,
                  col = "magenta3",
                  lwd = 3
            )
          }
        } else if (input$NIL == TRUE) {
          cyto_plot(fs.comp()[[input$flowFrame]],
                    channels = c(input$xchannel, input$ychannel),
                    overlay = fs.comp()[[input$Unstained]],
                    display = display,
                    axes_trans = axes_trans,
                    title = input$flowFrame,
                    point_size = 3,
                    limits = "machine",
                    axes_text_size = 1.7,
                    axes_label_text_size = 2,
                    title_text_size = 2, ...
          )
          
          if (input$median == TRUE) {
            medians <- fsApply(
              fs.comp(),
              each_col,
              median
            )[input$Unstained, channels]
            MFI <- data.frame("Channel" = channels, "Median" = medians)
            
            cutoff <- MFI[match(input$ychannel, MFI$Channel), ]
            
            abline(h = cutoff[2], col = "red", lwd = 3)
          }
          
          if (input$trace == TRUE) {
            cells <- exprs(fs.comp()[[input$flowFrame]])
            cells <- cells[order(cells[, input$xchannel]), ]
            cells <- as.data.frame(cells)
            
            n <- nrow(cells)
            splt <- seq(1, 20, 1)
            r <- ceiling(nrow(cells) / 20)
            cuts <- splt * r
            cells <- split(cells, cumsum(seq_len(nrow(cells)) %in% cuts))
            
            xmedians <- lapply(cells, function(x) median(x[, input$xchannel]))
            ymedians <- lapply(cells, function(x) median(x[, input$ychannel]))
            
            medians <- data.frame(unlist(xmedians), unlist(ymedians))
            colnames(medians) <- c(input$xchannel, input$ychannel)
            
            y <- medians[, input$ychannel]
            x <- medians[, input$xchannel]
            loessMod <- loess(y ~ x, data = medians, span = 0.9)
            loessMod <- predict(loessMod)
            
            lines(medians[, input$xchannel],
                  loessMod,
                  col = "purple",
                  lwd = 3
            )
          }
        }
        
        # Density distribution in channel
        cyto_plot(fs.comp()[[input$Unstained]],
                  channels = input$xchannel,
                  overlay = fs.comp()[[input$flowFrame]],
                  title = NA,
                  ylab = "Density",
                  density_fill = c("grey70", "white"),
                  density_fill_alpha = c(1, 0),
                  density_line_col = c("black", "blue"),
                  density_line_width = 2,
                  axes_trans = axes_trans,
                  axes_text_size = 1.7,
                  axes_label_text_size = 2,
                  title_text_size = 2.5
        )
        
        # Get median statistics for plot labels
        Nil.median <- median(
          exprs(
            .getRawData(fs.comp()[[input$Unstained]], axes_trans)
          )[, input$ychannel]
        )
        stain.median <- median(
          exprs(
            .getRawData(fs.comp()[[input$flowFrame]], axes_trans)
          )[, input$ychannel]
        )
        
        # Density distribution in selected channel
        cyto_plot(fs.comp()[[input$Unstained]],
                  channels = input$ychannel,
                  overlay = fs.comp()[[input$flowFrame]],
                  title = NA,
                  ylab = "Density",
                  density_fill = c("grey70", "white"),
                  density_fill_alpha = c(1, 0),
                  density_line_col = c("black", "red"),
                  density_line_width = 2,
                  axes_trans = axes_trans,
                  axes_text_size = 1.7,
                  axes_label_text_size = 2,
                  title_text_size = 2.5
        )
        
        # Add label for unstained median
        cyto_plot_label(fs.comp()[[input$Unstained]],
                        trans = axes_trans,
                        channels = input$ychannel,
                        text = "MedFI",
                        stat = "median",
                        text_x = 3.7,
                        text_y = 80,
                        text_col = "grey40",
                        text_size = 1.1
        )
        
        # Add label for Stained median
        cyto_plot_label(fs.comp()[[input$flowFrame]],
                        trans = axes_trans,
                        channels = input$ychannel,
                        text = "MedFI",
                        stat = "median",
                        text_x = 3.7,
                        text_y = 30,
                        text_col = "red",
                        text_size = 1.1
        )
      })
      
      output$plots <- renderPlot({
        xchan <- pd$channel[match(input$flowFrame2, pd$name)]
        if (input$Uncomp == TRUE) {
          if (input$Comp == TRUE & input$NIL2 == FALSE) {
            cyto_plot_compensation(fs[[input$flowFrame2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$flowFrame2]],
                                   point_col = c("blue", "red"),
                                   display = display, point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          } else if (input$Comp == TRUE & input$NIL2 == TRUE) {
            cyto_plot_compensation(fs[[input$flowFrame2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   overlay = list(
                                     fs.comp()[[input$flowFrame2]],
                                     fs.comp()[[input$Unstained2]]
                                   ),
                                   point_col = c("blue", "red", "black"),
                                   display = display,
                                   point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          } else if (input$Comp == FALSE & input$NIL2 == TRUE) {
            cyto_plot_compensation(fs[[input$flowFrame2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$Unstained2]],
                                   point_col = c("blue", "black"),
                                   display = display,
                                   point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          } else if (input$Comp == FALSE & input$NIL2 == FALSE) {
            cyto_plot_compensation(fs[[input$flowFrame2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   point_col = "blue",
                                   display = display,
                                   point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          }
        } else if (input$Uncomp == FALSE) {
          if (input$Comp == TRUE & input$NIL2 == TRUE) {
            cyto_plot_compensation(fs.comp()[[input$flowFrame2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$Unstained2]],
                                   point_col = c("blue", "black"),
                                   display = display,
                                   point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          } else if (input$Comp == FALSE & input$NIL2 == TRUE) {
            cyto_plot_compensation(fs.comp()[[input$Unstained2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   point_col = "black",
                                   display = display,
                                   point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          } else if (input$Comp == TRUE & input$NIL2 == FALSE) {
            cyto_plot_compensation(fs.comp()[[input$flowFrame2]],
                                   channel = xchan,
                                   axes_trans = axes_trans,
                                   point_col = "red",
                                   display = display,
                                   point_alpha = 0.6,
                                   header = input$flowFrame2,
                                   title = NA,
                                   point_size = 3,
                                   axes_text_size = 1.7,
                                   axes_label_text_size = 2,
                                   title_text_size = 2,
                                   header_text_size = 1.5
            )
          } else if (input$Comp == FALSE & input$NIL2 == FALSE) {
            
          }
        }
      })
      
      observe({
        input$saveBtn
        class(values$spill) <- "numeric"
        spill.mat <- values$spill / 100
        write.csv(spill.mat, spillover)
      })
    }
  )
}

#' Edit Spillover Matrix GatingSet Method
#'
#' Edit spillover matrices in real-time using a shiny interface.
#'
#' \code{cyto_spillover_edit} provides an interactive shiny interface for
#' editing fluorescent spillover matrices. \code{cyto_spillover_edit} takes on
#' either a \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' untransformed single stain compensation controls and a universal unstained
#' control. It is recommended that samples be pre-gated based on FSC and SSC
#' parameters to obtain a homogeneous population for calculation of fluorescent
#' spillover. Users begin by selecting the unstained control and a stained
#' control from dropdown menus of sample names. \code{cyto_spillover_edit}
#' leverages \code{cyto_plot} to plot the stained sample and overlay the
#' unstained control in black. Users should then select the channel associated
#' with the selected control on the \code{x axis} and go through all other
#' channels on the \code{y axis}. The displayed spillover matrix is extracted
#' directly from the \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spillover argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time. To aid in
#' selection of an appropriate spillover value, the median fluorescent intensity
#' of the unstained control is indicated by a red line and median fluorescent
#' intensity of the stained control is tracked with a purple line. These
#' features can be turned off by de-selecting the check boxes. Changes to the
#' spillover matrix are automatically saved to a csv file called
#' \code{"Spillover-Matrix.csv"} in the case where the \code{spillover} is not
#' specified or to the same name as the specified \code{spillover}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the pre-gated population to be plotted (e.g. "Single
#'   Cells").
#' @param channel_match name of .csv file containing the names of the samples in
#'   a column called "name" and their matching channel in a column called
#'   "channel". \code{cyto_spillover_edit} will the guide you through the
#'   channel selection process and generate a channel match file called
#'   "Compensation-Channels.csv" automatically. If you already have a complete
#'   channel_match and would like to bypass the channel selection process,
#'   simply pass the name of the channel_match to this argument (e.g.
#'   "Compensation-Channels.csv").
#' @param spillover name of spillover matrix csv file including .csv file
#'   extension to use as a starting point for editing. If \code{spillover} is
#'   not supplied the spillover matrix will be extracted directly from the
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and the edited
#'   matrix saved to \code{"Spillover-Matrix.csv"}.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines. Set to all events by default.
#' @param axes_trans object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the flowSet.
#' @param ... additional arguments passed to
#'   \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return save edited spillover matrix to .csv file named
#'   "Spillover-matrix.csv" or spillover.
#'
#' @importFrom flowWorkspace sampleNames pData
#' @importFrom flowCore compensate fsApply sampleFilter exprs Subset each_col
#' @importFrom utils read.csv write.csv
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarPanel selectInput
#'   checkboxInput actionButton mainPanel plotOutput reactiveValues observe
#'   eventReactive renderPlot tabsetPanel tabPanel sidebarLayout fluidRow
#'   updateSelectInput
#' @importFrom rhandsontable rhandsontable rHandsontableOutput hot_to_r
#'   renderRHandsontable hot_cols hot_rows
#' @importFrom shinythemes shinytheme
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom stats median loess predict
#' @importFrom graphics lines
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_spillover_edit,flowSet-method}}
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#'
#' # Load in compensation controls
#' fs <- Compensation
#' gs <- GatingSet(fs)
#'
#' # Gate single cells using gate_draw
#' gt <- Compensation_gatingTemplate
#' gating(gt, gs)
#'
#' # Channel match file
#' cmfile <- system.file("extdata",
#'   "Compensation-Channels.csv",
#'   package = "CytoRSuiteData"
#' )
#'
#' # Edit spillover matrix
#' cyto_spillover_edit(gs,
#'   parent = "Single Cells",
#'   channel_match = cmfile
#' )
#' }
#' @export
cyto_spillover_edit.GatingSet <- function(x,
                                          parent = NULL,
                                          channel_match = NULL,
                                          spillover = NULL,
                                          display = 1,
                                          axes_trans = NULL, ...) {
  
  # Assign x to gs
  gs <- x
  
  # Extract sample names
  nms <- cyto_names(gs)
  
  # Extract channels
  channels <- cyto_fluor_channels(gs)
  
  # Extract population for downstream plotting
  if (!is.null(parent)) {
    fs <- cyto_extract(gs, parent)
  } else if (is.null(parent)) {
    fs <- cyto_extract(gs, cyto_nodes(gs)[length(cyto_nodes(gs))])
  }
  
  # Make call to cyto_spillover_edit flowSet method
  cyto_spillover_edit(
    x = fs,
    channel_match = channel_match,
    spillover = spillover,
    axes_trans = axes_trans,
    display = display, ...
  )
}
