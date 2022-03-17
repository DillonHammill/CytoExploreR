## RHANDSONTABLE DATA EDITOR ---------------------------------------------------

#' Interactive data editor
#'
#' @param x object of class matrix or data.frame with colnames specified.
#' @param title title to include in above the table.
#' @param type can be either "editor", "menu" or "selector", set to "editor" by
#'   default. The editor option allows complete editing of the data_matrix with
#'   text input. The menu option converts the entries in the last column to
#'   checkboxes for logical inputs. The selector option converts entries in the
#'   last column to dropdown menus that contain choices supplied to the
#'   \code{options} argument.
#' @param options vector of options to use in dropdown menus in the last column.
#' @param save_as name of a csv file to which the edited table should be saved.
#' @param viewer logical indicating whether the table editor should be launched
#'   in the RStudio viewer pane, set to TRUE by default. Make sure to use the
#'   "Save & Close" to close the data editor or else the changes will not be
#'   saved.
#' @param logo path to image to include a logo in data editor title panel, set
#'   to CytoExploreR by default.
#'
#' @importFrom shiny shinyApp fluidPage titlePanel mainPanel runApp onStop img
#'   div paneViewer observeEvent
#' @importFrom bslib bs_theme
#' @importFrom rhandsontable rHandsontableOutput renderRHandsontable
#'   rhandsontable hot_col hot_to_r
#' @importFrom utils read.csv write.csv
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
data_editor <- function(x,
                        title = "Data Editor",
                        type = "editor",
                        options = NULL,
                        save_as = NULL,
                        viewer = TRUE,
                        logo = "CytoExploreR") {
  
  # SELECTOR OPTIONS
  if(type == "selector" & is.null(options)){
    stop("options required to use slector type.")
  }else if(type == "selector" & !is.null(options)){
    # EMPTY DEFAULT
    if(!any(LAPPLY(options, ".empty"))){
      options <- c(options, "")
    }
  }
  
  # PULL DOWN ROW NAMES
  row_names <- rownames(x)
  
  # CONVERT TO CHRACTER MATRIX
  if(is(x, "data.frame")){
    x <- as.matrix(x)
  }
  
  # EDITOR DATA
  if(type == "editor"){
    # MOVE COLNAMES INTO MATRIX
    if(is.null(colnames)){
      stop("colnames must be assigned!")
    }else{
      x <- rbind(colnames(x), x)
    }
  # MENU DATA
  }else if(type %in% "menu"){
    x[, ncol(x)] <- rep(NA, nrow(x))
  # SELECTOR DATA
  }else if(type == "selector"){
    x[, ncol(x)] <- rep("", nrow(x))
  }
  
  # CONVERT TO DATA FRAME
  x <- data.frame(x, stringsAsFactors = FALSE)
  
  # CytoExploreR logo from GitHub man folder
  if(logo == "CytoExploreR"){
    logo <- paste0(
      "https://raw.githubusercontent.com/DillonHammill/CytoExploreR",
      "/master/man/figures/logo.png"
    )
  }
  
  # TEMP_FILE
  temp_file <- NULL
  
  # DATA EDITOR
  app <- shinyApp(
    
    # USER INTERFACE
    ui <- fluidPage(
      theme = bs_theme(version = 3, bootswatch = "yeti"),
      titlePanel(div(img(src = logo, width = 100), title)),
      mainPanel(rHandsontableOutput("x")),
      actionButton("save_and_close", "Save & Close")
    ),
    
    # SERVER
    server <- function(input, output, session) {
      
      # VALUES
      values <- reactiveValues()
      
      # DATA EDITS
      observe({
        if (!is.null(input$x)) {
          values[["x"]] <- hot_to_r(input$x)
        } else {
          values[["x"]] <- x
        }
        write.csv(values[["x"]],
                  temp_file,
                  row.names = FALSE)
      })
      
      # TABLE

      output$x <- renderRHandsontable({
        # EDITOR
        if(type == "editor"){
          suppressWarnings(
            rhandsontable(values[["x"]],
                          contextMenu = TRUE,
                          useTypes = FALSE,
                          colHeaders = NULL,
                          rowHeaders = NULL,
                          readOnly = FALSE,
                          halign = "htCenter")
          )
        # MENU
        }else if(type == "menu"){
          suppressWarnings(
            rhandsontable(values[["x"]],
                          contextMenu = TRUE,
                          rowHeaders = NULL) %>%
            hot_col(colnames(x)[-length(colnames(x))],
                    halign = "htCenter",
                    readOnly = TRUE,
                    useTypes = FALSE) %>%
            hot_col(col = colnames(x)[length(colnames(x))],
                    type = "checkbox",
                    halign = "htCenter",
                    readOnly = FALSE,
                    allowInvalid = FALSE) 
          )
        # SELECTOR
        }else if(type == "selector"){
          # SELECTOR
          suppressWarnings(
            rhandsontable(values[["x"]],
                          contextMenu = TRUE,
                          rowHeaders = NULL,
                          manualColumnResize = TRUE) %>%
            hot_col(col = colnames(x)[length(colnames(x))],
                    type = "dropdown",
                    source = options,
                    halign = "htCenter",
                    readOnly = FALSE)
          )
        }

      })

      # MANUAL CLOSE
      observeEvent(input$save_and_close, {
        stopApp({
          dm <- read.csv(temp_file,
                         header = TRUE,
                         stringsAsFactors = FALSE)})
          unlink(temp_file)
          return(dm)
      })
      
    },
    
    # CREATE TEMP FILE
    onStart <- function(){
      temp_file <<- tempfile(fileext = ".csv")
    }
  )

  # RUN DATA EDITOR - INTERACTIVE MODE ONLY
  if(getOption("CytoExploreR_interactive")){
    if (viewer == TRUE) {
      x <- runApp(app,
                  launch.browser = paneViewer(),
                  quiet = TRUE)
    } else {
      x <- runApp(app,
                  quiet = TRUE)
    }
  }

  # EDITOR DATA
  if(type == "editor"){
    # UPDATE COLUMN NAMES
    colnames(x) <- x[1, ]
    # REMOVE COLNAMES FROM MATRIX
    x <- x[-1, , drop = FALSE]
    # CONVERT NUMERIC CHARACTERS
    lapply(seq_len(ncol(x)), function(z){
      # NUMBERS TO NUMERIC
      if(all(grepl("^[0-9 ]+$", x[, z]))){
        x[, z] <<- as.numeric(x[, z])
     }
    })
    # CONVERT EMPTY CHARACTERS TO NA
    lapply(seq_len(ncol(x)), function(z){
      if(any(LAPPLY(x[, z], ".empty"))){
        x[,z][which(LAPPLY(x[, z], ".empty"))] <<- NA
      }
    })
    # ADD BACK ROW NAMES
    rownames(x) <- row_names
  # SELECTOR DATA
  }else if(type == "selector"){
    # CONVERT EMPTY CHARACTERS TO NA
    lapply(seq_len(ncol(x)), function(z){
      if(any(LAPPLY(x[, z], ".empty"))){
        x[,z][which(LAPPLY(x[, z], ".empty"))] <<- NA
      }
    })
  }
  
  # WRITE TO FILE
  if(!is.null(save_as)){
    write.csv(x,
              save_as,
              row.names = FALSE)
  }
  
  # RETURN UPDATED DATA MATRIX
  return(x)
  
}