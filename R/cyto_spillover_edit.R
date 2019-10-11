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
#' for calculation of fluorescent spillover. The data should not be transformed
#' prior to using \code{cyto_spillover_edit} as transformations will be applied
#' internally when required. Users can supply their own custom transformers to
#' \code{axes_trans} or the default biexponential transformers will be used.
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
#'   (i.e. edit the spillover matrix constructed on the cytometer).
#' @param axes_trans an object of class \code{transformerList} containing
#'   transformers to used internally to transform the fluorescent channels of
#'   the samples for visualisation. \code{cyto_spillover_edit} expects
#'   un-transformed data and will automatically resort to using biexponential
#'   transformers internally if no transformers are supplied to this argument.
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
#'   updateSelectInput onStop stopApp runApp updateCheckboxInput
#' @importFrom rhandsontable rhandsontable rHandsontableOutput hot_to_r
#'   renderRHandsontable hot_cols hot_rows
#' @importFrom shinythemes shinytheme
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom stats median loess predict
#' @importFrom graphics lines layout
#' @importFrom tools file_ext
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
                                           display = 2000,
                                           point_size = 3,
                                           axes_text_size = 1.7,
                                           axes_label_text_size = 2,
                                           title_text_size = 2,
                                           header_text_size = 1.5, ...) {
  
  # Assign x to gs
  gs <- x
  
  # Extract sample names
  nms <- cyto_names(gs)
  
  # Extract channels
  channels <- cyto_fluor_channels(gs)
  
  # Message channels should not be transformed
  trans <- gs[[1]]@transformation
  
  # Extract population for downstream plotting
  if (!is.null(parent)) {
    fs <- cyto_extract(gs, parent)
  } else if (is.null(parent)) {
    fs <- cyto_extract(gs, cyto_nodes(gs)[length(cyto_nodes(gs))])
  }  

  # Extract transformers from GatingSet
  trans <- x[[1]]@transformation
  
  # Reverse any applied transformations to get LINEAR data
  if(any(channels %in% names(trans))){
    # Extract transformations that have been applied
    trans <- cyto_transformer_combine(trans[names(trans) %in% channels])
    # Inverse transformations
    fs <- cyto_transform(fs, 
                         trans,
                         inverse = TRUE,
                         plot = FALSE)
  }
  
  # Get complete transformerList
  axes_trans <- .cyto_transformer_complete(x, axes_trans)
  
  # Make call to cyto_spillover_edit flowSet method
  spill <- cyto_spillover_edit(
    x = fs,
    channel_match = channel_match,
    spillover = spillover,
    axes_trans = axes_trans,
    display = display,
    point_size = point_size,
    axes_text_size = axes_text_size,
    axes_label_text_size = axes_label_text_size,
    title_text_size = title_text_size,
    header_text_size = header_text_size, ...)
  
  return(spill)
  
}

#' @rdname cyto_spillover_edit
#' @export
cyto_spillover_edit.flowSet <- function(x,
                                         channel_match = NULL,
                                         spillover = NULL,
                                         axes_trans = NULL,
                                         display = 2000,
                                         point_size = 3,
                                         axes_text_size = 1.7,
                                         axes_label_text_size = 2,
                                         title_text_size = 2,
                                         header_text_size = 1.5, ...) {
  
  # Assign x to fs
  fs <- x
  
  # Extract sample names
  nms <- cyto_names(fs)
  
  # Extract fluorescent channels
  channels <- cyto_fluor_channels(fs)
  
  # Extract cyto details to pd
  pd <- cyto_details(fs)
  
  # Channel match file supplied
  if(!is.null(channel_match)){
  
    # channel_match is a data.frame or matrix or tibble
    if(inherits(channel_match, "data.frame") |
       inherits(channel_match, "matrix")){
      # channel_match must contain "name" and "channel" columns
      if(!any(grepl("name", colnames(channel_match), ignore.case = TRUE)) |
         !any(grepl("channel", colnames(channel_match), ignore.case = TRUE))) {
        stop("'channel_match' must contain the columns 'name' and 'channel'.")
      }
    # channel_match is the name of a csv file
    }else{
      # channel_match muct contain csv file extension
      if(!file_ext(channel_match) == "csv"){
        channel_match <- paste0(channel_match, ".csv")
      }
      # Working directory check
      if(getOption("CytoExploreR_wd_check") == TRUE){
        # channel_match exists in working directory
        if(file_wd_check(channel_match)){
          channel_match <- read.csv(channel_match, 
                                    header = TRUE, 
                                    row.names = 1,
                                    stringsAsFactors = FALSE)
        }else{
          stop(paste(channel_match, " is not in this working directory."))
        }
      # Bypass working directory checks
      }else{
        channel_match <- read.csv(channel_match, 
                                  header = TRUE, 
                                  row.names = 1,
                                  stringsAsFactors = FALSE)
      }
    }
    
    # Names of samples don't match any listed in channel_match
    if(!any(nms %in% row.names(channel_match))){
      # Add NA channel selections to pd
      pd$channel <- rep(NA, nrow(pd))
    # channel_match contains selections for supplied samples
    }else{
      # Add NA channel selections to pd
      pd$channel <- rep(NA, nrow(pd))
      # Replace with channel selections from channel_match
      lapply(nms, function(z){
        # Update channel selction if covered by channel_match file
        if(z %in% row.names(channel_match)){
          chn <- channel_match$channel[match(z, row.names(channel_match))]
          pd[match(z, pd$name),"channel"] <<- chn
        }
      })
    }
    
  # No channel match file supplied  
  }else if(is.null(channel_match)){
    # Add NA channel selections to pd
    pd$channel <- rep(NA, nrow(pd))
  }
  
  # Spillover matrix supplied
  if (!is.null(spillover)) {
    # spillover is a data.frame or matrix or tibble
    if (inherits(spillover, "matrix") |
        inherits(spillover, "data.frame")) {
      # spill should be a matrix
      spill <- spillover
      spill <- as.matrix(spill)
      # Save edited spillover matrix to date-Spillover-Matrix.csv
      spillover <- paste0(format(Sys.Date(), "%d%m%y"), 
                          "-", "Spillover-Matrix.csv")
    # spillover is the name of csv file
    } else {
      # File extension missing
      if(!file_ext(spillover) == "csv"){
        spillover <- paste0(spillover, ".csv")
      }
      # Working directory checks
      if (getOption("CytoExploreR_wd_check") == TRUE) {
        if (file_wd_check(spillover)) {
          spill <- read.csv(spillover, header = TRUE, row.names = 1)
          spill <- as.matrix(spill)
        } else {
          stop(paste(spillover, "is not in this working directory."))
        }
      # Bypass working directory checks 
      } else {
        spill <- read.csv(spillover, header = TRUE, row.names = 1)
        spill <- as.matrix(spill)
      }
    }
    
  # No spillover matrix supplied
  } else {
    spillover <- paste0(format(Sys.Date(), "%d%m%y"), 
                        "-", "Spillover-Matrix.csv")
    # Use spillover matrix attached to first sample
    spill <- fs[[1]]@description$SPILL
  }
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  # Unstained supplied?
  if(any(grepl("Unstained", pd$channel, ignore.case = TRUE))){
    unstained_initial <- pd$name[grepl("Unstained", 
                                       pd$channel, 
                                       ignore.case = TRUE)][1]
  }else{
    unstained_initial <- NA
  }
  
  # Get complete transformerList
  axes_trans <- .cyto_transformer_complete(fs, axes_trans)
  
  # Selected sample - unstained control supplied
  if(!.all_na(pd$channel)){
    # Unstained control included
    if(any(grepl("Unstaied", pd$channel))){
    editor_initial_sample <- pd$name[!grepl("Unstained", 
                                            pd$channel, 
                                            ignore.case = TRUE)][1]
    # No unstained control - use first sample
    }else{
      editor_initial_sample <- pd$name[1]
    }
  }else{
    editor_initial_sample <- pd$name[1]
  }
  
  # X channel selection
  if(!.all_na(pd$channel[pd$name == editor_initial_sample])){
    editor_initial_xchannel <- pd$channel[pd$name == editor_initial_sample]
  }else{
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
        if(!is.null(input$spillover_matrix)) {
          spill <- hot_to_r(input$spillover_matrix)
          rownames(spill) <- colnames(spill)
          values$spill <- spill
        }else{
          values$spill <- spill * 100
        }
      })
      
      output$spillover_matrix <- renderRHandsontable({
        rhandsontable(values$spill,
                      rowHeaderWidth = 105,
                      readOnly = FALSE) %>%
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
        if(!.all_na(xchan)){
          updateSelectInput(session, "editor_xchannel", selected = xchan)
          updateSelectInput(session, "plots_xchannel", selected = xchan)
        }
      })
      
      # Re-apply compensation and transformations
      fs.comp <- eventReactive(values$spill, {
        # Data must be LINEAR at this step
        fs <- compensate(fs, values$spill / 100)
        
        # Get transformed data
        fs <- cyto_transform(fs, axes_trans, plot = FALSE)
        
        return(fs)
      })
      
      # Turn off unstained aspects if no unstained control is supplied
      observe({
        unst <- input$editor_unstained
        # Turn off unstained components if unstained is NA
        if(unst == "No Unstained Control"){
          updateCheckboxInput(session,"editor_unstained_overlay",value = FALSE)
          updateCheckboxInput(session, "editor_unstained_median", value = FALSE)
        # Turn on unstained components
        }else{
          updateCheckboxInput(session,"editor_unstained_overlay",value = TRUE)
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
        if(unst == "No Unstained Control"){
          updateCheckboxInput(session, "plots_unstained_overlay", value = FALSE)
        # Turn on unstained overlay
        }else{
          updateCheckboxInput(session, "plots_unstained_overlay", value = TRUE)
        }
        # Update unstained control selection in editor tab
        updateSelectInput(session, "editor_unstained", selected = unst)
      })
      
      output$editor_plots <- renderPlot({
        
        # Set up plot layout
        layout(matrix(c(1,1,2,2,1,1,3,3), byrow = TRUE, ncol = 4))
        
        # No unstained control - no unstained overlay or unstained median
        if(input$editor_unstained == "No Unstained Control") {
          
          # Plots
          cyto_plot(fs.comp()[[input$editor_sample]],
                    channels = c(input$editor_xchannel,
                                 input$editor_ychannel),
                    axes_trans = axes_trans,
                    display = display,
                    title = input$editor_sample,
                    point_size = point_size,
                    axes_text_size = axes_text_size,
                    axes_label_text_size = axes_label_text_size,
                    title_text_size = title_text_size, ...)
          
          # Median tracker
          if(input$editor_median_tracker == TRUE){
            
            # Calculate medians for loess fitting
            cells <- exprs(fs.comp()[[input$editor_sample]])
            cells <- cells[order(cells[, input$editor_xchannel]), ]
            cells <- as.data.frame(cells)
            
            n <- nrow(cells)
            splt <- seq(1, 20, 1)
            r <- ceiling(nrow(cells) / 20)
            cuts <- splt * r
            cells <- split(cells, cumsum(seq_len(nrow(cells)) %in% cuts))
            
            xmedians <- lapply(cells, 
                               function(x){ 
                                 median(x[, input$editor_xchannel])
                                 })
            ymedians <- lapply(cells, 
                               function(x){ 
                                 median(x[, input$editor_ychannel])
                                 })
            
            medians <- data.frame(unlist(xmedians), unlist(ymedians))
            colnames(medians) <- c(input$editor_xchannel, 
                                   input$editor_ychannel)
            
            vals.y <- medians[, input$editor_ychannel]
            vals.x <- medians[, input$editor_xchannel]
            loessMod <- loess(vals.y ~ vals.x,
                              data = medians, span = 0.9
            )
            loessMod <- predict(loessMod)
            
            lines(medians[, input$editor_xchannel],
                  loessMod,
                  col = "purple2",
                  lwd = 3
            )
            
          }

          # Density distribution in associated channel
          cyto_plot(fs.comp()[[input$editor_sample]],
                    channels = input$editor_xchannel,
                    axes_trans = axes_trans,
                    display = display,
                    title = NA,
                    ylab = "Density",
                    density_fill = "white",
                    density_line_col = "blue",
                    density_line_width = 2,
                    axes_text_size = axes_text_size,
                    axes_label_text_size = axes_label_text_size,
                    title_text_size = title_text_size)

          # Density distribution in other channel
          cyto_plot(fs.comp()[[input$editor_sample]],
                    channels = input$editor_ychannel,
                    axes_trans = axes_trans,
                    display = display,
                    title = NA,
                    ylab = "Density",
                    density_fill = "white",
                    density_line_col = "blue",
                    density_line_width = 2,
                    axes_text_size = axes_text_size,
                    axes_label_text_size = axes_label_text_size,
                    title_text_size = title_text_size)

          # Add MedFI label in other channel
          cyto_plot_label(fs.comp()[[input$editor_sample]],
                           channels = input$editor_ychannel,
                           trans = axes_trans,
                           display = display,
                           label_text = "MedFI",
                           label_stat = "median",
                           label_text_x = 0.75*par("usr")[2],
                           label_text_y = 30,
                           label_text_col = "red",
                           label_text_size = 1.1
          )
          
        # Unstained control supplied  
        }else if(input$editor_unstained != "No Unstained Control"){
          
          # Unstained Overlay
          if(input$editor_unstained_overlay == FALSE){
            
            # Plot
            cyto_plot(fs.comp()[[input$editor_sample]],
                      channels = c(input$editor_xchannel, 
                                   input$editor_ychannel),
                      axes_trans = axes_trans,
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
              cells <- exprs(fs.comp()[[input$editor_sample]])
              cells <- cells[order(cells[, input$editor_xchannel]), ]
              cells <- as.data.frame(cells)
              
              n <- nrow(cells)
              splt <- seq(1, 20, 1)
              r <- ceiling(nrow(cells) / 20)
              cuts <- splt * r
              cells <- split(cells, cumsum(seq_len(nrow(cells)) %in% cuts))
              
              xmedians <- lapply(cells, 
                                 function(x){ 
                                   median(x[, input$editor_xchannel])
                                   })
              ymedians <- lapply(cells, 
                                 function(x){
                                   median(x[, input$editor_ychannel])
                                   })
              
              medians <- data.frame(unlist(xmedians), unlist(ymedians))
              colnames(medians) <- c(input$editor_xchannel, 
                                     input$editor_ychannel)
              
              vals.y <- medians[, input$editor_ychannel]
              vals.x <- medians[, input$editor_xchannel]
              loessMod <- loess(vals.y ~ vals.x,
                                data = medians, span = 0.9
              )
              loessMod <- predict(loessMod)
              
              lines(medians[, input$editor_xchannel],
                    loessMod,
                    col = "purple2",
                    lwd = 3
              )
              
            }
            
          # No unstained overlay  
          }else if(input$editor_unstained_overlay == TRUE){
            
            # Plot
            cyto_plot(fs.comp()[[input$editor_sample]],
                      channels = c(input$editor_xchannel, 
                                   input$editor_ychannel),
                      overlay = fs.comp()[[input$editor_unstained]],
                      axes_trans = axes_trans,
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
              cells <- exprs(fs.comp()[[input$editor_sample]])
              cells <- cells[order(cells[, input$editor_xchannel]), ]
              cells <- as.data.frame(cells)
              
              n <- nrow(cells)
              splt <- seq(1, 20, 1)
              r <- ceiling(nrow(cells) / 20)
              cuts <- splt * r
              cells <- split(cells, cumsum(seq_len(nrow(cells)) %in% cuts))
              
              xmedians <- lapply(cells, function(x){
                median(x[, input$editor_xchannel])
                })
              ymedians <- lapply(cells, function(x){
                median(x[, input$editor_ychannel])
                })
              
              medians <- data.frame(unlist(xmedians), unlist(ymedians))
              colnames(medians) <- c(input$editor_xchannel, 
                                     input$editor_ychannel)
              
              y <- medians[, input$editor_ychannel]
              x <- medians[, input$editor_xchannel]
              loessMod <- loess(y ~ x, data = medians, span = 0.9)
              loessMod <- predict(loessMod)
              
              lines(medians[, input$editor_xchannel],
                    loessMod,
                    col = "purple2",
                    lwd = 3
              )
            }
            
          }
          
          # Density distribution in associated channel - unstained overlay
          cyto_plot(fs.comp()[[input$editor_unstained]],
                    channels = input$editor_xchannel,
                    overlay = fs.comp()[[input$editor_sample]],
                    axes_trans = axes_trans,
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
                    title_text_size = 1.25 * title_text_size)
          
          # Density distribution in other channel - unstained overlay
          cyto_plot(fs.comp()[[input$editor_unstained]],
                    channels = input$editor_ychannel,
                    overlay = fs.comp()[[input$editor_sample]],
                    axes_trans = axes_trans,
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
                           label_text_x = 0.75*par("usr")[2],
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
                           label_text_x = 0.75*par("usr")[2],
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
                                   density_fill = c("blue","red")
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
                                   density_fill = c("grey", "red","blue")
            )
          } else if (input$plots_compensated_overlay == FALSE & 
                     input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$plots_unstained]],
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
                                   density_fill = c("grey","blue")
            )
          } else if (input$plots_compensated_overlay == FALSE & 
                     input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
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
      onStop(function(){
        spill.mat <- read.csv(spillover, header = TRUE, row.names = 1)
        colnames(spill.mat) <- rownames(spill.mat)
        stopApp(spill.mat)
      })
      
    }
    
  )
  
  # Run the shiny application
  sp <- runApp(app)
  
  # Return updated spillover matrix
  return(sp)
  
}
  
