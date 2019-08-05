# CYTO_SPILLOVER_EDIT ----------------------------------------------------------

#' Interactively Edit Spillover Matrices in Real-Time
#' 
#' #' @return save edited spillover matrix to .csv file named
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
#' @importFrom tools file_ext
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @name cyto_spillover_edit2
NULL

#' @noRd
#' @export
cyto_spillover_edit2 <- function(x, ...) {
  UseMethod("cyto_spillover_edit2")
}

#' @rdname cyto_spillover_edit2
#' @export
cyto_spillover_edit2.flowSet <- function(x,
                                         channel_match = NULL,
                                         spillover = NULL,
                                         axes_trans = NULL,
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
      if(getOption("CytoRSuite_wd_check") == TRUE){
        # channel_match exists in working directory
        if(file_wd_check(channel_match)){
          channel_match <- read.csv(channel_match, header = TRUE, row.names = 1)
        }else{
          stop(paste(channel_match, " is not in this working directory."))
        }
      # Bypass working directory checks
      }else{
        channel_match <- read.csv(channel_match, header = TRUE, row.names = 1)
      }
    }
    
    # Add channel selections to pd
    ind <- match(cyto_names(fs),row.names(channel_match))
    chans <- channel_match$channel[ind]
    pd$channel <- paste(chans)
    
  # No channel match file supplied  
  }else if(!is.null(channel_match)){
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
      if (getOption("CytoRSuite_wd_check") == TRUE) {
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
                                       ignore.case = TRUE)[1]]
  }else{
    unstained_initial <- NA
  }
  
  # Transformations - data untransformed but transformers can be supplied
  if(!is.null(axes_trans)){
    if(!inherits(axes_trans, "transformerList")){
      stop("'axes_trans' must be an object of class transformerList.")
    }
    message("cyto_spillover_edit expects untransformed samples.")
  # Transformations missing resort to logicle transformation (flow cytometry)
  }else{
    message("Applying logicle transformations to visualise the data.")
    axes_trans <- cyto_transformer_logicle(fs, plot = FALSE)
  }
  
  # Shiny application
  shinyApp(
    ui <- fluidPage(
      theme = shinytheme("yeti"),
      titlePanel("CytoEXploreR Spillover Matrix Editor"),
      tabsetPanel(
        tabPanel("Editor",
                 fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel(
                     selectInput(
                       inputId = "editor_unstained",
                       label = "Select Unstained Control:",
                       choices = c(nms, NA),
                       selected = unstained_initial
                     ),
                     selectInput(
                       inputId = "editor_sample",
                       label = "Select Sample:",
                       choices = nms
                     ),
                     selectInput(
                       inputId = "editor_xchannel",
                       label = "X Axis:",
                       choices = channels
                     ),
                     selectInput(
                       inputId = "editor_ychannel",
                       label = "Y Axis:",
                       choices = channels
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
                       choices = c(nms, NA),
                       selected = unstained_initial
                     ),
                     selectInput(
                       inputId = "plots_sample",
                       label = "Select Sample:",
                       choices = nms
                     ),
                     selectInput(
                       inputId = "plots_xchannel",
                       label = "Channel",
                       choices = channels
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
        xchan <- pd$channel[match(input$editor_sample, pd$name)]
        # Only update x channel if not NA
        if(!.all_na(xchan)){
          updateSelectInput(session, editor_xchannel, selected = xchan)
        }
        # Update sample selection in Plots tab 
        updateSelectInput(session, "plots_sample", selected = fr)
      })
      
      # Re-apply compensation and transformations
      fs.comp <- eventReactive(values$spill, {
        # Data must be LINEAR at this step
        fs <- compensate(fs, values$spill / 100)
        
        # Get trnsformed data
        fs <- cyto_transform(fs, axes_trans, plot = FALSE)
        
        return(fs)
      })
      
      # Turn off unstained aspects if no unstained control is supplied
      observe({
        unst <- input$editor_unstained
        # Turn off unstained components if unstained is NA
        if(.all_na(unst)){
          # remove unstained overlay in editor
          updateCheckboxInput(session,"editor_unstained_overlay",value = FALSE)
          updateCheckboxInput(session, "editor_unstained_median", value = FALSE)
        }
        # Update Plots tab unstained control based on editor selection
        updateSelectInput(session, "plots_unstained", selected = unst)
      })
      
      # Update channel selection on plots tab based on editor selection
      observe({
        xchan <- input$editor_xchannel
        updateSelectInput(session, "plots_xchannel", selected = xchan)
      })
      
      # Update overlay unstained check box based on unstained selection (Plots)
      observe({
        unst <- input$plots_unstained
        if(.all_na(unst)){
          updateCheckboxInput(session, "plots_unstained_overlay", value = FALSE)
        }
      })
      
      output$editor_plots <- renderPlot({
        
        layout(rbind(1, 1, 2, 2), c(1, 1, 3, 3))
        
        # No unstained control - no unstained overlay or unstained median
        if(.all_na(input$editor_unstained)) {
          
          # Plots
          cyto_plot(fs.comp()[[input$editor_sample]],
                    channels = c(input$editor_xchannel,
                                 input$editor_ychannel),
                    axes_trans = axes_trans,
                    title = input$editor_sample,
                    point_size = point_size,
                    axes_text_size = axes_text_size,
                    axes_label_text_size = axes_label_text_size,
                    title_text_size = title_text_size)
          
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
            colnames(medians) <- c(input$editor_xchannel, input$editor_ychannel)
            
            vals.y <- medians[, input$editor_ychannel]
            vals.x <- medians[, input$editor_xchannel]
            loessMod <- loess(vals.y ~ vals.x,
                              data = medians, span = 0.9
            )
            loessMod <- predict(loessMod)
            
            lines(medians[, input$editor_xchannel],
                  loessMod,
                  col = "cyan3",
                  lwd = 3
            )
            
          }
          
          # Density distribution in associated channel
          cyto_plot(fs.comp()[[input$editor_sample]],
                    channels = input$editor_xchannel,
                    axes_trans = axes_trans,
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
                    title = NA,
                    ylab = "Density",
                    density_fill = "white",
                    density_line_col = "blue",
                    density_line_width = 2,
                    axes_text_size = axes_text_size,
                    axes_label_text_size = axes_label_text_size,
                    title_text_size = title_text_size)
          
          # Add MedFI label in other channel
          cyto_plot_label2(fs.comp()[[input$editor_sample]],
                           trans = axes_trans,
                           channels = input$editor_ychannel,
                           label_text = "MedFI",
                           label_stat = "median",
                           label_text_x = 3.7,
                           label_text_y = 30,
                           label_text_col = "red",
                           label_text_size = 1.1
          )
          
        # Unstained control supplied  
        }else if(!.all_na(input$editor_unstained)){
          
          # Unstained Overlay
          if(input$editor_unstained_overlay){
            
            # Plot
            cyto_plot(fs.comp()[[input$editor_sample]],
                      channels = c(input$editor_xchannel, 
                                   input$editor_ychannel),
                      axes_trans = axes_trans,
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
                    col = "cyan3",
                    lwd = 3
              )
              
            }
            
          # No unstained overlay  
          }else if(!input$editor_unstained_overlay){
            
            # Plot
            cyto_plot(fs.comp()[[input$editor_sample]],
                      channels = c(input$editor_xchannel, 
                                   input$editor_ychannel),
                      overlay = fs.comp()[[input$editor_unstained]],
                      axes_trans = axes_trans,
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
                    col = "cyan3",
                    lwd = 3
              )
            }
            
          }
          
          # Density distribution in associated channel - unstained overlay
          cyto_plot(fs.comp()[[input$editor_unstained]],
                    channels = input$editor_xchannel,
                    overlay = fs.comp()[[input$editor_sample]],
                    axes_trans = axes_trans,
                    title = NA,
                    ylab = "Density",
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
                    title = NA,
                    ylab = "Density",
                    density_fill = c("grey70", "white"),
                    density_fill_alpha = c(1, 0),
                    density_line_col = c("black", "red"),
                    density_line_width = 2,
                    axes_trans = axes_trans,
                    axes_text_size = axes_text_size,
                    axes_label_text_size = axes_label_text_size,
                    title_text_size = 1.25 * title_text_size
          )
          
          # Add label for unstained median
          cyto_plot_label2(fs.comp()[[input$editor_unstained]],
                           trans = axes_trans,
                           channels = input$editor_ychannel,
                           label_text = "MedFI",
                           label_stat = "median",
                           label_text_x = 3.7,
                           label_text_y = 80,
                           label_text_col = "grey40",
                           label_text_size = 1.1
          )
          
          # Add label for Stained median
          cyto_plot_label2(fs.comp()[[input$editor_unstained]],
                           trans = axes_trans,
                           channels = input$editor_ychannel,
                           label_text = "MedFI",
                           label_stat = "median",
                           label_text_x = 3.7,
                           label_text_y = 30,
                           label_text_col = "red",
                           label_text_size = 1.1
          )
        }
        
      })
      
      # Plots tab uses cyto_plot_compensation
      output$plots <- renderPlot({
        
        # Compensation plots
        if (input$plots_uncompensated_underlay == TRUE) {
          if (input$plots_compensated_overlay == TRUE & 
              input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$plots_sample]],
                                   point_col = c("blue", "red"),
                                   point_alpha = 0.6,
                                   header = input$plots_sample,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
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
                                   point_col = c("blue", "red", "black"),
                                   point_alpha = 0.6,
                                   header = input$plots_samples,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
            )
          } else if (input$plots_compensated_overlay == FALSE & 
                     input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$plots_unstained]],
                                   point_col = c("blue", "black"),
                                   point_alpha = 0.6,
                                   header = input$plots_sample,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
            )
          } else if (input$plots_compensated_overlay == FALSE & 
                     input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   point_col = "blue",
                                   point_alpha = 0.6,
                                   header = input$plots_sample,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
            )
          }
        } else if (input$plots_uncompensated_underlay == FALSE) {
          if (input$plots_compensated_overlay == TRUE & 
              input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs.comp()[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   overlay = fs.comp()[[input$plots_unstained]],
                                   point_col = c("blue", "black"),
                                   point_alpha = 0.6,
                                   header = input$plots_sample,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
            )
          } else if (input$plots_compensated_overlay == FALSE & 
                     input$plots_unstained_overlay == TRUE) {
            cyto_plot_compensation(fs.comp()[[input$plots_unstained]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   point_col = "black",
                                   point_alpha = 0.6,
                                   header = input$plots_sample,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
            )
          } else if (input$plots_compensated_overlay == TRUE & 
                     input$plots_unstained_overlay == FALSE) {
            cyto_plot_compensation(fs.comp()[[input$plots_sample]],
                                   channel = input$plots_xchannel,
                                   axes_trans = axes_trans,
                                   point_col = "red",
                                   point_alpha = 0.6,
                                   header = input$plots_sample,
                                   title = NA,
                                   point_size = point_size,
                                   axes_text_size = axes_text_size,
                                   axes_label_text_size = axes_label_text_size,
                                   title_text_size = title_text_size,
                                   header_text_size = header_text_size
            )
          } else if (input$plots_compensated_overlay == FALSE & 
                     input$plots_unstained_overlay == FALSE) {
            # No data to plot
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

#' @rdname cyto_spillover_edit2
#' @export
cyto_spillover_edit2.GatingSet <- function(x,
                                           parent = NULL,
                                           channel_match = NULL,
                                           spillover = NULL,
                                           axes_trans = NULL,
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
  
  # Extract population for downstream plotting
  if (!is.null(parent)) {
    fs <- cyto_extract(gs, parent)
  } else if (is.null(parent)) {
    fs <- cyto_extract(gs, cyto_nodes(gs)[length(cyto_nodes(gs))])
  }
  
  # Transformations - data untransformed but transformers can be supplied
  if(!is.null(axes_trans)){
    if(!inherits(axes_trans, "transformerList")){
      stop("'axes_trans' must be an object of class transformerList.")
    }
    message("cyto_spillover_edit expects untransformed samples.")
    # Transformations missing resort to logicle transformation (flow cytometry)
  }else{
    message("Applying logicle transformations to visualise the data.")
    axes_trans <- cyto_transformer_logicle(fs, plot = FALSE)
  }
  
  # Make call to cyto_spillover_edit flowSet method
  cyto_spillover_edit(
    x = fs,
    channel_match = channel_match,
    spillover = spillover,
    axes_trans = axes_trans,
    point_size = point_size,
    axes_text_size = axes_text_size,
    axes_label_text_size = axes_label_text_size,
    title_text_size = title_text_size,
    header_text_size = header_text_size, ...)
  
}
