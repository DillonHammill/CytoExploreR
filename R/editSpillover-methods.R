#' Edit Spillover Matrix
#'
#' Edit spillover matrices in real-time using a shiny interface.
#'
#' \code{editSpillover} provides an interactive shiny interface for editing
#' fluorescent spillover matrices. \code{editSpillover} takes on either a
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' untransformed single stain compensation controls and a universal unstained
#' control. It is recommended that samples be pre-gated based on FSC and SSC
#' parameters to obtain a homogeneous population for calculation of fluorescent
#' spillover. Users begin by selecting the unstained control and a stained
#' control from dropdown menus of sample names. \code{editSpillover} leverages
#' \code{plotCyto} to plot the stained sample and overlay the unstained control
#' in black. Users should then select the channel associated with the selected
#' control on the \code{x axis} and go through all other channels on the \code{y
#' axis}. The displayed spillover matrix is extracted directly from the
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spfile argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time. To aid in
#' selection of an appropriate spillover value, the median fluorescent intensity
#' of the unstained control is indicated by a red line and median fluorescent
#' intensity of the stained control is tracked with a purple line. These
#' features can be turned off by de-selecting the check boxes. Changes to the
#' spillover matrix are automatically saved to a csv file called
#' \code{"Spillover Matrix.csv"} in the case where the \code{spfile} is not
#' specified or to the same name as the specified \code{spfile}.
#' \code{editSpillover} has methods for both
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
#' @seealso \code{\link{editSpillover,flowSet-method}}
#' @seealso \code{\link{editSpillover,GatingSet-method}}
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#'
#' @export
setGeneric(name = "editSpillover",
           def = function(x, ...){standardGeneric("editSpillover")}
)


#' Edit Spillover Matrix - flowSet Method
#'
#' Edit spillover matrices in real-time using a shiny interface.
#'
#' \code{editSpillover} provides an interactive shiny interface for editing
#' fluorescent spillover matrices. \code{editSpillover} takes on either a
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' untransformed single stain compensation controls and a universal unstained
#' control. It is recommended that samples be pre-gated based on FSC and SSC
#' parameters to obtain a homogeneous population for calculation of fluorescent
#' spillover. Users begin by selecting the unstained control and a stained
#' control from dropdown menus of sample names. \code{editSpillover} leverages
#' \code{plotCyto} to plot the stained sample and overlay the unstained control
#' in black. Users should then select the channel associated with the selected
#' control on the \code{x axis} and go through all other channels on the \code{y
#' axis}. The displayed spillover matrix is extracted directly from the
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spfile argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time. To aid in
#' selection of an appropriate spillover value, the median fluorescent intensity
#' of the unstained control is indicated by a red line and median fluorescent
#' intensity of the stained control is tracked with a purple line. These
#' features can be turned off by de-selecting the check boxes. Changes to the
#' spillover matrix are automatically saved to a csv file called
#' \code{"Spillover Matrix.csv"} in the case where the \code{spfile} is not
#' specified or to the same name as the specified \code{spfile}.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param spfile name of spillover matrix csv file including .csv file extension
#'   to use as a starting point for editing. If \code{spfile} is not supplied
#'   the spillover matrix will be extracted directly from the
#'   \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param subSample numeric indicating the number of events to plot, set to 5000
#'   events by default.
#' @param transList \code{\link[flowCore:transformList-class]{transformList}}
#'   object generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the flowSet.
#' @param ... additional arguments passed to \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return save edited spillover matrix to .csv file named "Spillover
#'   Matrix.csv" or spfile.
#'
#' @importFrom flowWorkspace sampleNames pData
#' @importFrom flowCore compensate fsApply sampleFilter exprs Subset each_col
#' @importFrom utils read.csv write.csv
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarPanel selectInput
#'   checkboxInput actionButton mainPanel plotOutput reactiveValues observe
#'   eventReactive renderPlot
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
#' @seealso \code{\link{editSpillover,GatingSet-method}}
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#'
#' @export
setMethod(editSpillover, signature = "flowSet", definition = function(x, spfile = NULL, subSample = 5000, transList = NULL, ...){
  
  require(shiny)
  require(shinythemes)
  require(rhandsontable)
  
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
  
  # Rhandsontable does handle decimal points if none are supplied in the dataset - if matrix is empty edit first value in second column to 0.0001
  if(all(spill %in% c(0,1))){
    
    spill[1,2] <- 0.0001
    
  }
  
  shinyApp(
    ui <- fluidPage(
      
      theme = shinytheme("yeti"),
      
      titlePanel("cytoRSuite Spillover Matrix Editor"),
      
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
      
      # Apply compensation after each edit
      fs.comp <- eventReactive(values$spill, {
        
        fs <- compensate(fs, values$spill/100)
        
        # Apply logicle transformation to fluorescent channels
        if(!is.null(transList)){
          
          if(class(transList)[1] != "transformList"){
            
            stop("Please supply a transList object of class transformList.")
            
          }else{
            
            fs <- transform(fs, transList)
            
          }
          
        }else{
          
          transList <<- estimateLogicle(as(fs,"flowFrame"), channels)
          fs <- transform(fs, transList)
          
        }
        
        return(fs)
        
      })
      
      output$plot <- renderPlot({
        
        # Axes Limits
        fr.exprs <- exprs(fs.comp()[[input$flowFrame]])
        xrange <- range(fr.exprs[, input$xchannel])
        yrange <- range(fr.exprs[, input$ychannel])
        
        # X Axis Limits
        if(xrange[1] > 0){
          
          xmin <- -0.5
          
        }else if(xrange[1] <= 0){
          
          xmin <- xrange[1]
        }
        xlim <- c(xmin,4.5)
        
        
        # Y Axis Limits
        if(yrange[1] > 0){
          
          ymin <- -0.5
          
        }else if(yrange[1] <= 0){
          
          ymin <- yrange[1]
        } 
        ylim <- c(ymin,4.5)
        
        # Plot
        if(input$NIL == FALSE){
          
          plotCyto(fs.comp()[[input$flowFrame]], channels = c(input$xchannel,input$ychannel), subSample = subSample, transList = transList, main = input$flowFrame, cex.pts = 3, xlim = xlim, ylim = ylim, ...)
          
          if(input$median == TRUE){
            
            medians <- fsApply(fs.comp(), each_col, "median")[input$Unstained,channels]
            MFI <- data.frame("Channel" = channels, "Median" = medians)
            
            cutoff <- MFI[match(input$ychannel, MFI$Channel),]
            
            abline(h = cutoff[2], col = "red", lwd = 2)
            
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
            
            loessMod <- loess(medians[,input$ychannel] ~ medians[,input$xchannel], data = medians, span = 0.9)
            loessMod <- predict(loessMod)
            
            lines(medians[,input$xchannel], loessMod, col = "magenta3", lwd = 3)
            
          }
          
        }else if(input$NIL == TRUE){
          
          plotCyto(fs.comp()[[input$flowFrame]], channels = c(input$xchannel,input$ychannel), overlay = fs.comp()[[input$Unstained]], subSample = subSample, transList = transList, main = input$flowFrame, cex.pts = 3, xlim = xlim, ylim = ylim, ...)
          
          if(input$median == TRUE){
            
            medians <- fsApply(fs.comp(), each_col, median)[input$Unstained,channels]
            MFI <- data.frame("Channel" = channels, "Median" = medians)
            
            cutoff <- MFI[match(input$ychannel, MFI$Channel),]
            
            abline(h = cutoff[2], col = "red", lwd = 2)
            
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
            
            loessMod <- loess(medians[,input$ychannel] ~ medians[,input$xchannel], data = medians, span = 0.9)
            loessMod <- predict(loessMod)
            
            lines(medians[,input$xchannel], loessMod, col = "magenta3", lwd = 3)
            
          }
          
        }
        
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
#' Edit spillover matrices in real-time using a shiny interface.
#'
#' \code{editSpillover} provides an interactive shiny interface for editing
#' fluorescent spillover matrices. \code{editSpillover} takes on either a
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} containing
#' untransformed single stain compensation controls and a universal unstained
#' control. It is recommended that samples be pre-gated based on FSC and SSC
#' parameters to obtain a homogeneous population for calculation of fluorescent
#' spillover. Users begin by selecting the unstained control and a stained
#' control from dropdown menus of sample names. \code{editSpillover} leverages
#' \code{plotCyto} to plot the stained sample and overlay the unstained control
#' in black. Users should then select the channel associated with the selected
#' control on the \code{x axis} and go through all other channels on the \code{y
#' axis}. The displayed spillover matrix is extracted directly from the
#' \code{\link[flowCore:flowSet-class]{flowSet}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} unless another
#' spillover matrix is supplied through the spfile argument. To edit the
#' spillover matrix simply modify the appropriate cell in the the table. The new
#' spillover matrix will be re-applied to the samples with each edit and
#' automatically re-plotted so you can track changes in real-time. To aid in
#' selection of an appropriate spillover value, the median fluorescent intensity
#' of the unstained control is indicated by a red line and median fluorescent
#' intensity of the stained control is tracked with a purple line. These
#' features can be turned off by de-selecting the check boxes. Changes to the
#' spillover matrix are automatically saved to a csv file called
#' \code{"Spillover Matrix.csv"} in the case where the \code{spfile} is not
#' specified or to the same name as the specified \code{spfile}.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the pre-gated population to be plotted (e.g. "Single
#'   Cells").
#' @param spfile name of spillover matrix csv file including .csv file extension
#'   to use as a starting point for editing. If \code{spfile} is not supplied
#'   the spillover matrix will be extracted directly from the
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param subSample numeric indicating the number of events to plot, set to 5000
#'   events by default.
#' @param transList object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace:transformerList]{transformerList}}
#'   generated by \code{\link[flowCore:logicleTransform]{estimateLogicle}} which
#'   was used to transform the fluorescent channels of the flowSet.
#' @param ... additional arguments passed to \code{\link{plotCyto,flowFrame-method}}.
#'
#' @return save edited spillover matrix to .csv file named "Spillover
#'   Matrix.csv" or spfile.
#'
#' @importFrom flowWorkspace sampleNames pData
#' @importFrom flowCore compensate fsApply sampleFilter exprs Subset each_col
#' @importFrom utils read.csv write.csv
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarPanel selectInput
#'   checkboxInput actionButton mainPanel plotOutput reactiveValues observe
#'   eventReactive renderPlot
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
#' @seealso \code{\link{editSpillover,flowSet-method}}
#' @seealso \code{\link{plotCyto1d,flowFrame-method}}
#' @seealso \code{\link{plotCyto2d,flowFrame-method}}
#'
#' @export
setMethod(editSpillover, signature = "GatingSet", definition = function(x, parent = NULL, spfile = NULL, subSample = 5000, transList = NULL, ...){
  
  require(shiny)
  require(shinythemes)
  require(rhandsontable)
  
  gs <- x
  nms <- sampleNames(gs)
  channels <- getChannels(gs)
  
  # Transformations
  if(is.null(transList)){
    
    trnsfrms <- lapply(channels, function(channel){getTransformations(gs[[1]], channel, only.function = FALSE)})
    names(trnsfrms) <- channels
  
    # Remove NULL transforms
    trnsfrms[sapply(trnsfrms, is.null)] <- NULL
  
    if(length(trnsfrms) == 0){
    
      transList <- NULL
    
    }
    
  }
  
  # Convert to transformerList
  if(!is.null(transList)){
    
    transList <- transformerList(names(transList), transList)
    
  }
  
  # Extract population for downstream plotting
  if(!is.null(parent)){
    
    fs <- getData(gs, parent)
    
  }else if(is.null(parent)){
    
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
      
      titlePanel("cytoRSuite Spillover Matrix Editor"),
      
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
      
      # Apply compensation after each edit
      fs.comp <- eventReactive(values$spill, {
        
        fs <- compensate(fs, values$spill/100)
        
        # Apply logicle transformation to fluorescent channels
        if(!is.null(transList)){
          
          if(class(transList)[1] == "transformerList"){
            
            transList <<- transformList(names(transList),lapply(transList, `[[`, "transform"))
            
          }else if(class(transList)[1] == "transformList"){
            
            transList <<- transList
            fs <- transform(fs, transList)
            
          }else{
            
            stop("Supplied transList should be of class transformList or transFormerList.")
            
          }
          
        }else{
          
          transList <<- estimateLogicle(as(fs,"flowFrame"), channels)
          fs <- transform(fs, transList)
          
        }
        
        return(fs)
        
      })
      
      output$plot <- renderPlot({
        
        # Axes Limits
        fr.exprs <- exprs(fs.comp()[[input$flowFrame]])
        xrange <- range(fr.exprs[, input$xchannel])
        yrange <- range(fr.exprs[, input$ychannel])
        
        # X Axis Limits
        if(xrange[1] > 0){
          
          xmin <- -0.5
          
        }else if(xrange[1] <= 0){
          
          xmin <- xrange[1]
        }
        xlim <- c(xmin,4.5)
        
        
        # Y Axis Limits
        if(yrange[1] > 0){
          
          ymin <- -0.5
          
        }else if(yrange[1] <= 0){
          
          ymin <- yrange[1]
        } 
        ylim <- c(ymin,4.5)
        
        # Plot
        if(input$NIL == FALSE){
          
          plotCyto(fs.comp()[[input$flowFrame]], channels = c(input$xchannel,input$ychannel), subSample = subSample, transList = transList, main = input$flowFrame, cex.pts = 3, xlim = xlim, ylim = ylim, ...)
          
          if(input$median == TRUE){
            
            medians <- fsApply(fs.comp(), each_col, "median")[input$Unstained,channels]
            MFI <- data.frame("Channel" = channels, "Median" = medians)
            
            cutoff <- MFI[match(input$ychannel, MFI$Channel),]
            
            abline(h = cutoff[2], col = "red", lwd = 2)
            
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
            
            loessMod <- loess(medians[,input$ychannel] ~ medians[,input$xchannel], data = medians, span = 0.9)
            loessMod <- predict(loessMod)
            
            lines(medians[,input$xchannel], loessMod, col = "magenta3", lwd = 3)
            
          }
          
        }else if(input$NIL == TRUE){
          
          plotCyto(fs.comp()[[input$flowFrame]], channels = c(input$xchannel,input$ychannel), overlay = fs.comp()[[input$Unstained]], subSample = subSample, transList = transList, main = input$flowFrame, cex.pts = 3, xlim = xlim, ylim = ylim, ...)
          
          if(input$median == TRUE){
            
            medians <- fsApply(fs.comp(), each_col, median)[input$Unstained,channels]
            MFI <- data.frame("Channel" = channels, "Median" = medians)
            
            cutoff <- MFI[match(input$ychannel, MFI$Channel),]
            
            abline(h = cutoff[2], col = "red", lwd = 2)
            
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
            
            loessMod <- loess(medians[,input$ychannel] ~ medians[,input$xchannel], data = medians, span = 0.9)
            loessMod <- predict(loessMod)
            
            lines(medians[,input$xchannel], loessMod, col = "magenta3", lwd = 3)
            
          }
          
        }
        
      })
      
      observe({
        input$saveBtn
        class(values$spill) <- "numeric"
        spill.mat <- values$spill/100
        write.csv(spill.mat, spfile)
      })
      
    })
  
})