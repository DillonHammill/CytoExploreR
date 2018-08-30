#' Edit Spillover Matrix
#' 
#' \code{editSpillover} provides an interactive shiny interface for editing fluorescent spillover matrices. \code{editSpillover} takes on either a \code{flowSet} 
#' or \code{GatingSet} containing untransformed single stain compensation controls and a universal unstained control. It is recommended that samples be pre-gated based on FSC
#' and SSC parameters to obtain a homogeneous population for calculation of fluorescent spillover. Users begin by selecting the unstained control and a stained control
#' from dropdown menus of sample names. \code{editSpillover} leverages \code{ggcyto} to plot the stained sample and overlay the unstained control in black. Users should
#' then select the channel associated with the selected control on the \code{x axis} and go through all other channels on the \code{y axis}. The displayed spillover
#' matrix is extracted directly from the \code{flowSet} or \code{GatingSet} unless another spillover matrix is supplied through the spfile argument. To edit the spillover
#' matrix simply modify the appropriate cell in the the table. The new spillover matrix will be re-applied to the samples with each edit and automatically re-plotted
#' so you can track changes in real time. To add in selection of an appropriate spillover value, the median fluorescent intensity of the unstained control is indicated by
#' a red line and median fluorescent intensity of the stained control is tracked with a blue line. These features can be turned off by de-selecting the check boxes. Changes
#' to the spillover matrix are automatically saved to a csv file called \code{"Spillover Matrix.csv"} in the case where the \code{spfile} is not specified or to the same 
#' name as the specified \code{spfile}. \code{editSpillover} has methods for both \code{flowSet} and \code{GatingSet} objects 
#' refer to their respective help pages for more information - ?`editSpillover,flowSet-method` or ?`editSpillover,GatingSet-method`.
#' 
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(name = "editSpillover",
           def = function(x, ...){standardGeneric("editSpillover")}
)


#' Edit Spillover Matrix - flowSet Method
#' 
#' \code{editSpillover} provides an interactive shiny interface for editing fluorescent spillover matrices. \code{editSpillover} takes on either a \code{flowSet} 
#' or \code{GatingSet} containing untransformed single stain compensation controls and a universal unstained control. It is recommended that samples be pre-gated based on FSC
#' and SSC parameters to obtain a homogeneous population for calculation of fluorescent spillover. Users begin by selecting the unstained control and a stained control
#' from dropdown menus of sample names. \code{editSpillover} leverages \code{ggcyto} to plot the stained sample and overlay the unstained control in black. Users should
#' then select the channel associated with the selected control on the \code{x axis} and go through all other channels on the \code{y axis}. The displayed spillover
#' matrix is extracted directly from the \code{flowSet} or \code{GatingSet} unless another spillover matrix is supplied through the spfile argument. To edit the spillover
#' matrix simply modify the appropriate cell in the the table. The new spillover matrix will be re-applied to the samples with each edit and automatically re-plotted
#' so you can track changes in real time. To add in selection of an appropriate spillover value, the median fluorescent intensity of the unstained control is indicated by
#' a red line and median fluorescent intensity of the stained control is tracked with a blue line. These features can be turned off by de-selecting the check boxes. Changes
#' to the spillover matrix are automatically saved to a csv file called \code{"Spillover Matrix.csv"} in the case where the \code{spfile} is not specified or to the same 
#' name as the specified \code{spfile}. 
#'
#' @param x object of class \code{flowSet}.
#' @param spfile name of spillover matrix csv file including .csv file extension to use as a starting point for editing. If \code{spfile} is not supplied
#' the spillover matrix will be extracted directly from the \code{flowSet} or \code{GatingSet}.
#' @param subSample numeric indicating the number of events to plot, set to 5000 events by default.
#' @param ... additional arguments (not used).
#' 
#' @return save edited spillover matrix to .csv file named "Spillover Matrix.csv" or spfile.
#' 
#' @importFrom flowWorkspace sampleNames
#' @importFrom flowCore compensate fsApply sampleFilter exprs Subset each_col
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
          
          p <- p + geom_smooth(method = "loess", se = FALSE, mapping = aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)), data = medians, color = "cyan", size = 1.2)
          
        }
        
        # if transformation applied use axis_inverse_trans -5 < range <= 10
        # Ranges
        xrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel])
        yrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel])
        
        # X Axis
        if(xrange > -5 & xrange <= 10){
          
          p <- p + scale_x_continuous(limits = c(-4,5.5))
          
        }else{
          
          p <- p + scale_x_logicle(limits = c(-10000,262143))
          
        }
        
        # Y Axis
        if(yrange > -5 & yrange <= 10){
          
          p <- p + scale_y_continuous(limits = c(-4,5.5))
          
        }else{
          
          p <- p + scale_y_logicle(limits = c(-10000,262143))
          
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
#' \code{editSpillover} provides an interactive shiny interface for editing fluorescent spillover matrices. \code{editSpillover} takes on either a \code{flowSet} 
#' or \code{GatingSet} containing untransformed single stain compensation controls and a universal unstained control. It is recommended that samples be pre-gated based on FSC
#' and SSC parameters to obtain a homogeneous population for calculation of fluorescent spillover. Users begin by selecting the unstained control and a stained control
#' from dropdown menus of sample names. \code{editSpillover} leverages \code{ggcyto} to plot the stained sample and overlay the unstained control in black. Users should
#' then select the channel associated with the selected control on the \code{x axis} and go through all other channels on the \code{y axis}. The displayed spillover
#' matrix is extracted directly from the \code{flowSet} or \code{GatingSet} unless another spillover matrix is supplied through the spfile argument. To edit the spillover
#' matrix simply modify the appropriate cell in the the table. The new spillover matrix will be re-applied to the samples with each edit and automatically re-plotted
#' so you can track changes in real time. To add in selection of an appropriate spillover value, the median fluorescent intensity of the unstained control is indicated by
#' a red line and median fluorescent intensity of the stained control is tracked with a blue line. These features can be turned off by de-selecting the check boxes. Changes
#' to the spillover matrix are automatically saved to a csv file called \code{"Spillover Matrix.csv"} in the case where the \code{spfile} is not specified or to the same 
#' name as the specified \code{spfile}. 
#'
#' @param x object of class \code{GatingSet}.
#' @param name of the pre-gated population to use for downstream calculations, set to the last node of the GatingSet by default (e.g. "Single Cells").
#' @param spfile name of spillover matrix csv file including .csv file extension to use as a starting point for editing. If \code{spfile} is not supplied
#' the spillover matrix will be extracted directly from the \code{flowSet} or \code{GatingSet}.
#' @param subSample numeric indicating the number of events to plot, set to 5000 events by default.
#' @param ... additional arguments (not used).
#' 
#' @return save edited spillover matrix to .csv file named "Spillover Matrix.csv" or spfile.
#' 
#' @importFrom flowWorkspace sampleNames getData
#' @importFrom flowCore compensate fsApply sampleFilter exprs Subset each_col
#' 
#' @export
setMethod(editSpillover, signature = "GatingSet", definition = function(x, alias = NULL, spfile = NULL, subSample = 5000, ...){
  
  require(shiny)
  require(shinythemes)
  require(rhandsontable)
  require(ggcyto)
  
  gs <- x
  nms <- sampleNames(gs)
  channels <- getChannels(gs)
  
  # Extract population for downstream plotting
  if(!is.null(alias)){
    
    fs <- getData(gs, alias)
    
  }else if(is.null(alias)){
    
    fs <- getData(gs, getNodes(gs)[length(getNodes(gs))])
    
  }
  
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
          
          p <- p + geom_smooth(method = "loess", se = FALSE, mapping = aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)), data = medians, color = "cyan", size = 1.2)
          
        }
        
        # if transformation applied use axis_inverse_trans -5 < range <= 10
        # Ranges
        xrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$xchannel])
        yrange <- max(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel]) - min(exprs(fs.comp()[[input$flowFrame]])[,input$ychannel])
        
        # X Axis
        if(xrange > -5 & xrange <= 10){
          
          p <- p + scale_x_continuous(limits = c(-4,5.5))
          
        }else{
          
          p <- p + scale_x_logicle(limits = c(-10000,262143))
          
        }
        
        # Y Axis
        if(yrange > -5 & yrange <= 10){
          
          p <- p + scale_y_continuous(limits = c(-4,5.5))
          
        }else{
          
          p <- p + scale_y_logicle(limits = c(-10000,262143))
          
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