## CYTO_ARRAY_ANALYSE ----------------------------------------------------------

#' Analyse cytokine bead array data
#'
#' @param x GatingSet.
#' @param array array details.
#' @param model model to use.
#' @param channel detection channel.
#' @param format either long or wide return.
#' @param save_as csv file for results
#' @param popup whether standard curves should plotted in pop-up window, set to
#'   FALSE by default.
#' @param ... additional arguments passed to model.
#'
#' @importFrom tidyr spread
#' @importFrom graphics par mtext
#' @importFrom DataEditR data_edit
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @name cyto_array_analyse
NULL

#' @noRd
cyto_array_analyse <- function(x,
                               ...) {
  UseMethod("cyto_array_analyse")
}

#' @noRd
cyto_array_analyse.GatingSet <- function(x,
                                         array = NULL,
                                         model = LL.5(),
                                         channel = NULL,
                                         format = "wide",
                                         save_as = NULL,
                                         popup = FALSE,
                                         ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CHANNEL
  if(is.null(channel)){
    stop("Supply the name of the detection channel to 'channel'.")
  }else{
    channel <- cyto_channels_extract(x, channel)
  }
  
  # GRAPHICAL PARAMETERS
  old_pars <- .par(c("mfrow", "oma"))
  on.exit(par(old_pars))
  
  # CYTO_DETAILS
  cyto_details <- cyto_details(x)
  
  # ARRAY DETAILS --------------------------------------------------------------
  
  # ARRAY DETAILS - REQUIRED
  if (is.null(array)) {
    message("Update details of the bead array to continue...")
    array <- matrix(c(
      "A1",
      "IL-2",
      "10000",
      "4"
    ),
    ncol = 4,
    dimnames = list(
      NULL,
      c(
        "bead_id",
        "analyte",
        "max_concentration",
        "dilution_factor"
      )
    )
    )
    array <- data_edit(array,
                       title = "Bead Array Details Editor",
                       save_as = paste0(format(Sys.Date(), "%d%m%y"), 
                                        "-", 
                                        "Bead-Array-Details.csv"),
                       logo = CytoExploreR_logo(),
                       write_fun = "write_to_csv",
                       viewer = "pane"
    )
  } else {
    if (is(array, "character")) {
      # EDIT ARRAY DETAILS
      array <- data_editor(array,
                           title = "Bead Array Details Editor",
                           save_as = array,
                           logo = CytoExploreR_logo(),
                           write_fun = "write_to_csv"
      )
    } else {
      # EDIT ARRAY DETAILS
      array <- data_editor(array,
                           title = "Bead Array Details Editor",
                           save_as = paste0(format(Sys.Date(), "%d%m%y"), 
                                            "-", 
                                            "Bead-Array-Details.csv"),
                           logo = CytoExploreR_logo(),
                           write_fun = "write_to_csv")
    }
  }
  
  # DILUTION FACTOR - WHOLE NUMBER
  ind <- which(array[, "dilution_factor"] < 0 )
  if (length(ind) > 0) {
    array[ind, "dilution_factor"] <- 1 / array[ind, "dilution_factor"]
  }
  
  # BEAD MEDFI -----------------------------------------------------------------
  
  # BEAD POPULATION NAMES
  bead_pops <- unique(array[, "bead_id"])
  
  # COMPUTE MEDFI PER BEAD POPULATION
  message(
    paste0("Computing medFI for each bead population in ", channel, "...")
  )
  bead_medFI <- cyto_stats_compute(x,
                                   alias = bead_pops,
                                   stat = "median",
                                   channel = channel,
                                   format = "wide"
  )
  
  # WIDE FORMAT
  if(!all(bead_pops %in% colnames(bead_medFI))){
    bead_medFI <- spread(bead_medFI,
                         key = "Population",
                         value = colnames(bead_medFI)[ncol(bead_medFI)])
  }
  
  # SPLIT BY BEAD POPULATION - RETAIN CYTO_DETAILS
  bead_data <- lapply(seq_along(bead_pops), function(z){
    bd <- bead_medFI[, !colnames(bead_medFI) %in% bead_pops[-z]]
    colnames(bd)[ncol(bd)] <- "medFI"
    return(bd)
  })
  names(bead_data) <- bead_pops
  
  # STANDARDS DETAILS ----------------------------------------------------------
  
  # STANDARD COLUMN
  if (!"standard" %in% colnames(cyto_details)) {
    # EDIT CYTO_DETAILS
    cyto_details <- cbind(
      cyto_details,
      matrix(rep(NA, nrow(cyto_details)),
             ncol = 1,
             dimnames = list(
               NULL,
               "standard"
             )
      )
    )
    message("Annotate standards...")
    cyto_details <- data_editor(cyto_details,
                                title = "Experiment Details Editor",
                                logo = CytoExploreR_logo()
    )
  }
  
  # PREPARE STANDARDS ----------------------------------------------------------
  
  # STANDARD DATA
  message("Preparing standards...")
  std_data <- lapply(bead_pops, function(z){
    # PULL OUT STANDARDS 
    std <- bead_data[[z]][!bead_data[[z]][, "standard"] == "NA" & 
                            !is.na(bead_data[[z]][, "standard"]), ]
    # STANDRAD INDICES - REMOVE LEADING C
    if(all(LAPPLY(std[,"standard"], function(w){
      grepl("C", w, ignore.case = TRUE)
    }))){
      std_ind <- LAPPLY(std[, "standard"], function(w){
        as.numeric(gsub("c", "", w, ignore.case = TRUE))
      })
    }else{
      std_ind <- as.numeric(std[, "standard"])
    }
    # MAXIMUM STANDARD INDEX
    std_ind_max <- max(std_ind)
    # MAXIMUM CONCENTRATION
    max_conc <- array[array$bead_id == z, "max_concentration"]
    # STANDARD CONCENTRATIONS
    std_conc <- LAPPLY(std_ind, function(y){
      # ZERO CONCENTRATION
      if(std_ind_max - y == std_ind_max){
        return(0)
        # COMPUTE CONCENTRATION
      }else{
        return(max_conc / unique(array[array$bead_id == z, "dilution_factor"]) ^ 
                 (std_ind_max - y))
      }
    })
    std_conc <- matrix(std_conc, 
                       ncol = 1, 
                       dimnames = list(NULL, "concentration"))
    std <- cbind(std, std_conc)
    return(std)
  })
  names(std_data) <- paste(bead_pops, 
                           "-", 
                           array[array$bead_id == bead_pops, "analyte"])
  
  # PREPARE STANDARD CURVES ----------------------------------------------------
  
  # PREPARE MODEL
  if(!is(model, "llogistic")){
    stop("'model' must be of class llogistic.")
  }
  
  # POP-UP - CYTO_PLOT_SAVE
  if(getOption("cyto_plot_save")){
    popup <- FALSE
  }
  
  # NEW GRAPHICS DEVICE
  cyto_plot_new(popup)
  
  # LAYOUT
  par(mfrow = c(n2mfrow(length(std_data))[2],
                n2mfrow(length(std_data))[1]))
  
  # HEADER
  par(oma = c(0, 0, 3.5, 0))
  
  # MODELS
  message("Constructing standard curves...")
  bead_mod <- lapply(seq_along(std_data), function(z) {
    # MODEL
    std_mod <- drm(concentration  ~ pseudo_log_trans()$transform(medFI),
                   data = std_data[[z]],
                   fct = model,
                   ...
    )
    # PLOT
    plot(pseudo_log_trans()$transform(std_data[[z]][, "medFI"]),
         pseudo_log_trans()$transform(std_data[[z]][, "concentration"]),
         main = names(std_data)[z],
         ylab = expression(Concentration~log[10]*(pg/mL)),
         xlab = expression(medFI~log[10]*(FI)),
         col = "black",
         pch = 16,
         yaxt = "none")
    axis(2, c(0:4), labels = c(0,
                               expression(10^1),
                               expression(10^2),
                               expression(10^3),
                               expression(10^4)))
    # CURVE
    lines(pseudo_log_trans()$transform(std_data[[z]][, "medFI"]),
          pseudo_log_trans()$transform(predict(std_mod,
                                               newdata = pseudo_log_trans()$transform(std_data[[z]][, "medFI"]))),
          col = "red")
    # PLOT
    # plot(std_mod,
    #   type = "all",
    #   main = names(std_data)[z],
    #   ylab = expression(Concentration~log[10]*(pg/mL)),
    #   xlab = expression(medFI~log[10]*(FI)),
    #   col = "black",
    #   pch = 16,
    #   log = ""
    # )
    # HEADER
    if(z == length(std_data)){
      mtext("Standard Curves",
            outer = TRUE,
            font = 2,
            line = 0.5)
    }
    return(std_mod)
  })
  names(bead_mod) <- paste(bead_pops, 
                           "-", 
                           array[array$bead_id == bead_pops, "analyte"])
  
  # INTERPOLATE FROM STANDARD CURVES -------------------------------------------
  
  # NEW GRAPHICS DEVICE
  cyto_plot_new(popup)
  
  # LAYOUT
  par(mfrow = c(n2mfrow(length(std_data))[2],
                n2mfrow(length(std_data))[1]))
  
  # HEADER
  par(oma = c(0, 0, 3.5, 0))
  
  # PREDICT CONCENTRATIONS
  message("Interpolating concentrations...")
  bead_pred <- lapply(seq_along(bead_data), function(z){
    # REMOVE STANDARDS
    mod_data <- bead_data[[z]][bead_data[[z]][, "standard"] == "NA" | 
                                 is.na(bead_data[[z]][, "standard"]), ]
    # DATA TO INTERPOLATE
    mod_pred <- predict(bead_mod[[z]], 
                        as.matrix(mod_data[, "medFI", drop = FALSE]))
    
    # PLOT MODEL
    plot(bead_mod[[z]],
         type = "all",
         main = names(bead_data)[z],
         xlab = expression(Concentration~log[10]*(pg/mL)),
         ylab = expression(medFI~log[10]*(FI)),
         col = "black",
         pch = 16,
         log = "xy"
    )
    # PLOT INTERPOLATED POINTS
    points(mod_pred,
           as.matrix(mod_data[, "medFI"]),
           col = "blue",
           pch = 16)
    # HEADER
    if(z == length(bead_data)){
      mtext("Standard Curves Interpolation",
            outer = TRUE, 
            font = 2,
            line = 0.5)
    }
    # APPEND CONCENTRATIONS
    mod_conc <- matrix(mod_pred, 
                       ncol = 1, 
                       dimnames = list(NULL, "concentration"))
    mod_data <- cbind(bead_data[[z]], mod_conc)
    return(mod_data)
  })
  names(bead_mod) <- paste(bead_pops, 
                           "-", 
                           array[array$bead_id == bead_pops, "analyte"])
  
  # CYTO_PLOT_SAVE -------------------------------------------------------------
  
  # RESET CYTO_PLOT_SAVE
  if(getOption("cyto_plot_save")){
    # RESET
    options("cyto_plot_save" = FALSE)
    # CLOSE GRAPHICS DEVICE
    dev.off()
  }
  
  # FORMAT & SAVE RESULTS ------------------------------------------------------
  
  
  # SAVE RESULTS TO CSV FILE
  if (!is.null(save_as)) {
    write_to_csv(bead_data, save_as)
  }
  
  # RETURN PREDICTIONS
  return(bead_data)
}

#' @noRd
pseudo_log_trans <- function(sigma = 1, base = 10) {
  trans_new(
    "pseudo_log",
    function(x) asinh(x / (2 * sigma)) / log(base),
    function(x) 2 * sigma * sinh(x * log(base))
  )
}
