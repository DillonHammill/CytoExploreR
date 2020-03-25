## RHANDSONTABLE DATA EDITOR ---------------------------------------------------

#' Interactive data editor
#'
#' @param x object of class matrix or data.frame, matrices will be coerced to
#'   data.frames.
#' @param title title to include in above the table.
#' @param menu logical indicating whether the last column should contain
#'   dropdown menus, set to FALSE by default.
#' @param options vector of options to use in dropdpwn menus when menu is TRUE.
#' @param save_as name of a csv file to which the edited table should be saved.
#' @param viewer logical indicating whether the table editor should be launched
#'   in the RStudio viewer pane, set to TRUE by default. Mke sure to use the
#'   "Save & Close" to close the data editor or else the changes will not be
#'   saved.
#'
#' @importFrom shiny shinyApp fluidPage titlePanel mainPanel runApp onStop img
#'   div paneViewer
#' @importFrom shinythemes shinytheme
#' @importFrom rhandsontable rHandsontableOutput renderRHandsontable
#'   rhandsontable hot_col hot_to_r
#' @importFrom utils read.csv write.csv
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
data_editor <- function(x,
                        title = NULL,
                        menu = FALSE,
                        options = NULL,
                        save_as = NULL,
                        viewer = TRUE) {
  
  # DATA FRAME
  if (is(x, "matrix")) {
    x <- as.data.frame(x)
  }
  
  # SAVE_AS
  if (is.null(save_as)) {
    stop("Supply a file name to save_as.")
  }
  
  # WATCH OUT - ASIS
  if(any(LAPPLY(colnames(x), function(z){is(x[, z], "AsIs")}))){
    lapply(colnames(x), function(z){
      if(is(x[, z], "AsIs")){
        x[, z] <<- as.character(x[, z])
      }
    })
  }
  
  # ASSIGN X TO DATA_MATRIX
  data_matrix <- x
  
  # CytoExploreR logo from GitHub man folder
  logo <- paste0(
    "https://raw.githubusercontent.com/DillonHammill/CytoExploreR",
    "/master/man/figures/logo.png"
  )
  
  app <- shinyApp(
    
    # USER INTERFACE
    ui <- fluidPage(
      theme = shinytheme("yeti"),
      titlePanel(div(img(src = logo, width = 100), title)),
      mainPanel(rHandsontableOutput("data_matrix")),
      actionButton("save_and_close", "Save & Close")
    ),
    
    # SERVER
    server <- function(input, output, session) {
      
      # VALUES
      values <- reactiveValues()
      
      # DATA EDITS
      observe({
        if (!is.null(input$data_matrix)) {
          values[["data_matrix"]] <- hot_to_r(input$data_matrix)
        } else {
          values[["data_matrix"]] <- data_matrix
        }
        write.csv(values[["data_matrix"]],
                  save_as,
                  row.names = FALSE
        )
      })
      
      # ADD COLUMNS
      
      
      # TABLE
      output$data_matrix <- renderRHandsontable({
        if (menu == FALSE) {
          suppressWarnings(rhandsontable(values[["data_matrix"]],
                        contextMenu = TRUE,
                        useTypes = FALSE,
                        rowHeaders = NULL # remove row names
          ) %>% hot_cols(readOnly = FALSE,
                         halign = "htCenter"))
        } else if (menu == TRUE) {
          if (is.null(options)) {
            stop("Supply a vector of options for dropdown menus.")
          }
          suppressWarnings(rhandsontable(values[["data_matrix"]],
                                         contextMenu = TRUE,
                                         useTypes = FALSE,
                                         rowHeaders = NULL # remove row names
          ) %>%
            hot_col(colnames(x)[-length(colnames(x))],
                    halign = "htCenter",
                    readOnly = TRUE
            ) %>%
            hot_col(
              col = colnames(x)[length(colnames(x))],
              type = "dropdown",
              source = options,
              halign = "htCenter",
              readOnly = FALSE
            ))
        }
      })
      
      # MANUAL CLOSE
      observeEvent(input$save_and_close, {
        stopApp(read.csv(save_as, header = TRUE, stringsAsFactors = FALSE))
      })
      
      # SESSION ENDS
      session$onSessionEnded(function() {
        stopApp(read.csv(save_as, header = TRUE, stringsAsFactors = FALSE))
      })
    }
  )
  
  # RUN DATA EDITOR
  if (viewer == TRUE) {
    runApp(app, launch.browser = paneViewer())
  } else {
    runApp(app)
  }
}