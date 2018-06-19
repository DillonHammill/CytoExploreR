#' Edit Spillover Matrix
#' 
#' @param x object of class \code{flowSet} or \code{GatingSet}.
#' 
#' @export
setGeneric(name="editSpillover",
           def=function(x, ...){standardGeneric("editSpillover")}
)

#' Edit Spillover Matrix flowSet
#' 
#' @param x an object of class \code{flowSet}.
#' @param spfile name of spillover csv file to be edited (defaults to NULL). Spillover extracted from samples if not supplied.
#' @param gtfile name of gatingTemplate csv file to apply to flowSet prior to plotting (defaults to NULL). If not supplied samples will be gated FSC-A/SSC-A and
#' SSC-W/SSC-H to gate single cells. 
#' 
#' @export
setMethod(editSpillover, signature = "flowSet", definition = function(x, spfile = NULL, gtfile = NULL){
  
  require(shiny)
  require(shinythemes)
  require(rhandsontable)
  require(ggcyto)
  
  # Rename x to fs
  fs <- x
  
  # Add flowSet to GatingSet
  gs <- GatingSet(fs)
  
  # Transform fluorescent channels
  channels <- getChannels(gs)
  trans <- estimateLogicle(gs[[1]], channels)
  
  # Gate GatingSet using template if provided
  if(!is.null(gtfile)){
    
  gt <- gatingTemplate(gtfile)
  gating(gt, gs)
  
  }else{ # Construct gates manually
    
  gs <- drawGate(gs, parent = "root", alias = "Cells", channels = c("FSC-A","SSC-A"), gate_type = "polygon", file = "Compensation gatingTemplate.csv")
  gs <- drawGate(gs, parent = "Cells", alias = "Single Cells", channels = c("SSC-W","SSC-H"), gate_type = "polygon", template = "Template",file = "Compensation gatingTemplate.csv")
    
  }
  
  # Extract population for downstream plotting
  fs <- getData(gs, "Single Cells")
  
  # Read in spillover matrix to object spill
  if(!is.null(spfile)){
    
    spill <- read.csv(spfile)
    
  }else{
    
    spill <- fs[[1]]@description$SPILL
    
  }
  colnames(spill) <- channels
  rownames(spill) <- channels

  shinyApp(
  ui <- fluidPage(
    
    theme = shinytheme("yeti"),
    
    titlePanel("Edit Spillover Matrix"),
    
    sidebarPanel(
      selectInput(inputId = "Unstained", label = "Select Unstained Control:", choices = sampleNames(fs)),
      selectInput(inputId = "flowFrame", label = "Select sample:", choices = sampleNames(fs)),
      selectInput(inputId = "xchannel", label = "x axis:", choices = channels, selected = channels[1]),
      selectInput(inputId = "ychannel", label = "y axis:", choices = channels, selected = channels[2]),
      checkboxInput(inputId = "NIL", label = "Overlay Unstained Control", value = TRUE),
      checkboxInput(inputId = "median", label = "Unstained Control Median", value = TRUE),
      checkboxInput(inputId = "trace", label = "Median Tracker", value = TRUE),
      actionButton("saveBtn", "Save")
    ),
    
    mainPanel(
      rHandsontableOutput("spillover"),
      plotOutput("plot", height = "500px", width = "70%")
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
      
      
      rhandsontable(values$spill, rowHeaderWidth = 105, readOnly = FALSE) %>% hot_cols(type = "numeric", colWidths = 105, format = "0.00",
                                                                                       halign = "htCenter") %>% hot_rows(rowHeights = 20)
      
    })
    
    fs.comp <- eventReactive(values$spill, {
      
      gs <- GatingSet(fs)
      gs <- compensate(gs, values$spill/100)
      gs <- transform(gs, trans)
      
      return(gs)
      
    })
    
    
    output$plot <- renderPlot({
      
      p <- ggcyto(fs.comp()@data[[input$flowFrame]], max_nrow_to_plot = 10000, subset = "root", aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)))
      p <- p + geom_point(alpha = 0.4)
      
      if(input$NIL == TRUE){
        p <- p + geom_point(data = Subset(getData(fs.comp(),"root")[[input$Unstained]], sampleFilter(size = 10000)), alpha = 0.4,color = "red")
      }
      
      if(input$median == TRUE){      
        
        medians <- fsApply(fs.comp()@data, each_col, median)[input$Unstained,channels]
        MFI <- data.frame("Channel" = channels, "Median" = medians)
        
        cutoff <- MFI[match(input$ychannel, MFI$Channel),]
        
        p <- p + geom_hline(aes(yintercept = Median), color = "green", size = 1.2, data = cutoff)
      }
      
      if(input$trace == TRUE){
        cells <- exprs(getData(fs.comp(),"root")[[input$flowFrame]])
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
        
        p <- p + geom_line(mapping = aes_(x = as.name(input$xchannel), y = as.name(input$ychannel)), data = medians, color = "turquoise2", size = 1.2)
        
      }

      p <- p + scale_x_continuous(limits = c(NA,5))
      p <- p + scale_y_continuous(limits = c(NA,5))
      p <- p + axis_x_inverse_trans()
      p <- p + axis_y_inverse_trans()
      p <- p + facet_null()
      p <- p + ggtitle(paste(getData(fs.comp(),"root")[[input$flowFrame]]@description$GUID,"\n"))
      p <- p + xlab(paste("\n",input$xchannel))
      p <- p + ylab(paste(input$ychannel,"\n"))
      p <- p + editSpillover_theme()
      print(p)
      
    })
    
    observe({
      input$saveBtn
      class(values$spill) <- "numeric"
      spill.mat <- values$spill/100
      write.csv(spill.mat, "Spillover Matrix.csv")
    })
    
  })

  
})
