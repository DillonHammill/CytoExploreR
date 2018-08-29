#' Edit Spillover Matrix
#' 
#' Edit spillover matrix using single stain compensation controls and an unstained control.
#' 
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(name="editSpillover",
           def=function(x, ...){standardGeneric("editSpillover")}
)


#' Edit Spillover Matrix flowSet Method
#' 
#' @param x object of class \code{flowSet}.
#' @param spfile name of spillover matrix csv file including .csv file extension to use as a starting point for editing.
#' @param subSample numeric indicating the number of events to plot, set to 5000 events by default.
#' 
#' @return write edited spillover matrix to csv file called \code{"Spillover Matrix.csv"} for later use.
#' 
#' @export
setMethod(editSpillover, signature = "flowSet", definition = function(x, spfile = NULL, subSample = 5000, ...){
  
  require(shiny)
  require(shinythemes)
  require(rhandsontable)
  require(ggcyto)
  
  fs <- x
  nms <- sampleNames(fs)
  channels <- getChannels(fs)
  
  # Read in spillover matrix to object spill
  if(!is.null(spfile)){
    
    spill <- read.csv(spfile, header = TRUE, row.names = 1)
    spill <- as.matrix(spill)
    
  }else{
    
    spfile <- "Spillover Matrix.csv"
    spill <- fs[[1]]@description$SPILL
    
  }
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  # Rhandsontable does handle decimal points if none are supplied in the dataset - if matrix is empty edit first value in second column to 0.01
  if(all(spill %in% c(0,1))){
    
    spill[1,2] <- 0.0001
    
  }
  
  shinyApp(
    ui <- fluidPage(
      
      theme = shinytheme("yeti"),
      
      titlePanel("cytoSuite Spillover Matrix Editor"),
      
      sidebarPanel(
        selectInput(inputId = "Unstained", label = "Select Unstained Control:", choices = sampleNames(fs)),
        selectInput(inputId = "flowFrame", label = "Select Sample:", choices = sampleNames(fs)),
        selectInput(inputId = "xchannel", label = "X Axis:", choices = channels, selected = channels[1]),
        selectInput(inputId = "ychannel", label = "Y Axis:", choices = channels, selected = channels[2]),
        checkboxInput(inputId = "NIL", label = "Overlay Unstained Control", value = TRUE),
        checkboxInput(inputId = "median", label = "Unstained Control Median", value = TRUE),
        checkboxInput(inputId = "trace", label = "Median Tracker", value = TRUE),
        actionButton("saveBtn", "Save")
      ),
      
      mainPanel(
        rHandsontableOutput("spillover"),
        plotOutput("plot", height = "500px", width = "50%")
      )
    ),
    
    server = function(input, output, session){
      
      values <- reactiveValues()
      
      observe({
        
        if(!is.null(input$spillover)){
          
          spill <- hot_to_r(input$spillover)
          rownames(spill) <- colnames(spill)
          values$spill <- spill
          
        }else{
          
          values$spill <- spill*100
          
        }
        
      })    
      
      
      output$spillover <- renderRHandsontable({
        
        
        rhandsontable(values$spill, rowHeaderWidth = 105, readOnly = FALSE) %>% hot_cols(type = "numeric", colWidths = 105, format = "0.000", halign = "htCenter", renderer = "
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
                                                                                         }") %>% hot_rows(rowHeights = 20)
      
      })
      
      fs.comp <- eventReactive(values$spill, {
        
        fs <- compensate(fs, values$spill/100)
        
        return(fs)
        
      })
      
      
      output$plot <- renderPlot({
        
        p <- ggcyto(fs.comp()[[input$flowFrame]], max_nrow_to_plot = subSample, aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)))
        p <- p + geom_hex(bins = 150)
        
        if(input$NIL == TRUE){
          p <- p + geom_point(data = Subset(fs.comp()[[input$Unstained]], sampleFilter(size = subSample)), color = "black", size = 1.5, alpha = 0.6)
        }
        
        if(input$median == TRUE){      
          
          medians <- fsApply(fs.comp(), each_col, median)[input$Unstained,channels]
          MFI <- data.frame("Channel" = channels, "Median" = medians)
          
          cutoff <- MFI[match(input$ychannel, MFI$Channel),]
          
          p <- p + geom_hline(aes(yintercept = Median), color = "red", size = 1.2, data = cutoff)
        }
        
        if(input$trace == TRUE){
          cells <- exprs(fs.comp()[[input$flowFrame]])
          cells <- cells[order(cells[,input$xchannel]),]
          cells <- as.data.frame(cells)
          
          n <- nrow(cells)
          splt <- seq(1,20,1)
          r <- ceiling(nrow(cells)/20)
          cuts <- splt*r
          cells <- split(cells, cumsum(1:nrow(cells) %in% cuts))
          
          xmedians <- lapply(cells, function(x) median(x[,input$xchannel]))
          ymedians <- lapply(cells, function(x) median(x[,input$ychannel]))
          
          medians <- data.frame(unlist(xmedians), unlist(ymedians))
          colnames(medians) <- c(input$xchannel,input$ychannel)
          
          p <- p + geom_smooth(method = "loess", se = FALSE, mapping = aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)), data = medians, color = "magenta", size = 1.2)
          
        }
        
        # if transformation applied use axis_inverse_trans 0 < range <= 5
        # Ranges
        xrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel])
        yrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel])
        
        # X Axis
        if(xrange > -5 & xrange <= 10){
          
          p <- p + scale_x_continuous(limits = c(-4,5))
          
        }else{
          
          p <- p + scale_x_logicle(limits = c(-10000,250000))
          
        }
        
        # Y Axis
        if(yrange > -5 & yrange <= 10){
          
          p <- p + scale_y_continuous(limits = c(-4,5))
          
        }else{
          
          p <- p + scale_y_logicle(limits = c(-10000,250000))
          
        }
        
        p <- p + facet_null()
        p <- p + ggtitle(paste(input$flowFrame),"\n")
        p <- p + xlab(paste("\n",input$xchannel))
        p <- p + ylab(paste(input$ychannel,"\n"))
        p <- p + editSpillover_theme()
        print(p)
        
      })
      
      observe({
        input$saveBtn
        class(values$spill) <- "numeric"
        spill.mat <- values$spill/100
        write.csv(spill.mat, spfile)
      })
      
    })
  
  
})

#' Edit Spillover Matrix GatingSet Method
#'
#' @param x object of class \code{GatingSet}.
#' @param parent name of the population to use for plotting (defaults to "root").
#' @param spfile name of spillover matrix csv file including .csv file extension to use as a starting point for editing.
#' @param subSample numeric indicating the number of events to plot, set to 5000 events by default.
#' 
#' @return write edited spillover matrix to csv file called \code{"Spillover Matrix.csv"} for later use.
#' 
#' @export
setMethod(editSpillover, signature = "GatingSet", definition = function(x, parent = "root", spfile = NULL, subSample = 5000, ...){
  
  require(shiny)
  require(shinythemes)
  require(rhandsontable)
  require(ggcyto)
  
  gs <- x
  nms <- sampleNames(gs)
  channels <- getChannels(gs)
  
  # Extract population for downstream plotting
  fs <- getData(gs, parent)
  
  # Read in spillover matrix to object spill
  if(!is.null(spfile)){
    
    spill <- read.csv(spfile, header = TRUE, row.names = 1)
    spill <- as.matrix(spill)
    
  }else{
    
    spfile <- "Spillover Matrix.csv"
    spill <- fs[[1]]@description$SPILL
    
  }
  colnames(spill) <- channels
  rownames(spill) <- channels
  
  # Rhandsontable does handle decimal points if none are supplied in the dataset - if matrix is empty edit first value in second column to 0.01
  if(all(spill %in% c(0,1))){
    
    spill[1,2] <- 0.0001
    
  }
  
  shinyApp(
    ui <- fluidPage(
      
      theme = shinytheme("yeti"),
      
      titlePanel("cytoSuite Spillover Matrix Editor"),
      
      sidebarPanel(
        selectInput(inputId = "Unstained", label = "Select Unstained Control:", choices = sampleNames(fs)),
        selectInput(inputId = "flowFrame", label = "Select Sample:", choices = sampleNames(fs)),
        selectInput(inputId = "xchannel", label = "X Axis:", choices = channels, selected = channels[1]),
        selectInput(inputId = "ychannel", label = "Y Axis:", choices = channels, selected = channels[2]),
        checkboxInput(inputId = "NIL", label = "Overlay Unstained Control", value = TRUE),
        checkboxInput(inputId = "median", label = "Unstained Control Median", value = TRUE),
        checkboxInput(inputId = "trace", label = "Median Tracker", value = TRUE),
        actionButton("saveBtn", "Save")
      ),
      
      mainPanel(
        rHandsontableOutput("spillover"),
        plotOutput("plot", height = "500px", width = "50%")
      )
    ),
    
    server = function(input, output, session){
      
      values <- reactiveValues()
      
      observe({
        
        if(!is.null(input$spillover)){
          
          spill <- hot_to_r(input$spillover)
          rownames(spill) <- colnames(spill)
          values$spill <- spill
          
        }else{
          
          values$spill <- spill*100
          
        }
        
      })    
      
      
      output$spillover <- renderRHandsontable({
        
        
        rhandsontable(values$spill, rowHeaderWidth = 105, readOnly = FALSE) %>% hot_cols(type = "numeric", colWidths = 105, format = "0.000", halign = "htCenter", renderer = "
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
                                                                                         }") %>% hot_rows(rowHeights = 20)
      
      })
      
      fs.comp <- eventReactive(values$spill, {
        
        fs <- compensate(fs, values$spill/100)
        
        return(fs)
        
      })
      
      
      output$plot <- renderPlot({
        
        p <- ggcyto(fs.comp()[[input$flowFrame]], max_nrow_to_plot = subSample, aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)))
        p <- p + geom_hex(bins = 150)
        
        if(input$NIL == TRUE){
          p <- p + geom_point(data = Subset(fs.comp()[[input$Unstained]], sampleFilter(size = subSample)), color = "black", size = 1.5, alpha = 0.6)
        }
        
        if(input$median == TRUE){      
          
          medians <- fsApply(fs.comp(), each_col, median)[input$Unstained,channels]
          MFI <- data.frame("Channel" = channels, "Median" = medians)
          
          cutoff <- MFI[match(input$ychannel, MFI$Channel),]
          
          p <- p + geom_hline(aes(yintercept = Median), color = "red", size = 1.2, data = cutoff)
        }
        
        if(input$trace == TRUE){
          cells <- exprs(fs.comp()[[input$flowFrame]])
          cells <- cells[order(cells[,input$xchannel]),]
          cells <- as.data.frame(cells)
          
          n <- nrow(cells)
          splt <- seq(1,20,1)
          r <- ceiling(nrow(cells)/20)
          cuts <- splt*r
          cells <- split(cells, cumsum(1:nrow(cells) %in% cuts))
          
          xmedians <- lapply(cells, function(x) median(x[,input$xchannel]))
          ymedians <- lapply(cells, function(x) median(x[,input$ychannel]))
          
          medians <- data.frame(unlist(xmedians), unlist(ymedians))
          colnames(medians) <- c(input$xchannel,input$ychannel)
          
          p <- p + geom_smooth(method = "loess", se = FALSE, mapping = aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)), data = medians, color = "magenta", size = 1.2)
          
        }
        
        # if transformation applied use axis_inverse_trans 0 < range <= 5
        # Ranges
        xrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel])
        yrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel])
        
        # X Axis
        if(xrange > -5 & xrange <= 10){
          
          p <- p + scale_x_continuous(limits = c(-4,5))
          
        }else{
          
          p <- p + scale_x_logicle(limits = c(-10000,250000))
          
        }
        
        # Y Axis
        if(yrange > -5 & yrange <= 10){
          
          p <- p + scale_y_continuous(limits = c(-4,5))
          
        }else{
          
          p <- p + scale_y_logicle(limits = c(-10000,250000))
          
        }
        
        p <- p + facet_null()
        p <- p + ggtitle(paste(input$flowFrame),"\n")
        p <- p + xlab(paste("\n",input$xchannel))
        p <- p + ylab(paste(input$ychannel,"\n"))
        p <- p + editSpillover_theme()
        print(p)
        
      })
      
      observe({
        input$saveBtn
        class(values$spill) <- "numeric"
        spill.mat <- values$spill/100
        write.csv(spill.mat, spfile)
      })
      
    })
  
  
})