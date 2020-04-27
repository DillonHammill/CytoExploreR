## CYTO_PLOT CONSTRUCTORS ------------------------------------------------------

# Collection of internal wrappers used within cyto_plot to construct aspects of
# the plot. This includes wrappers for cyto_plot_empty, cyto_plot_density,
# cyto_plot_point, cyto_plot_gate and cyto_plot_label. These wrappers bypass
# some of the calculations in the exported equivalents as these calculations are
# already performed internally within cyto_plot. All of these functions inherit
# a named list of arguments from cyto_plot. Arguments must be extracted directly
# from the list to prevent R CMD CHECK global binding NOTEs.

## .CYTO_PLOT_EMPTY ------------------------------------------------------------

#' Construct empty cyto_plot
#'
#' @param args named list of cyto_plot arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_empty <- function(args) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REMOVE X ARGUMENT
  args <- args[-match("x", names(args))]

  # RENAME FR_LIST ARGUMENT TO X (DO.CALL() ON CYTO_PLOT_EMPTY)
  names(args)[match("fr_list", names(args))] <- "x"

  # CYTO_PLOT_EMPTY ARGUMENTS
  ARGS <- formalArgs("cyto_plot_empty.list")
  
  # PLOT CONSTRUCTION ----------------------------------------------------------

  # DENSITY SUPPLIED - (PULL OUT YLIM - BYPASS DENSITY CONSTRUCTION)
  if (length(args[["channels"]]) == 1 & 
      !.all_na(args[["fr_dens_list"]])) {

    # YLIM
    if (.all_na(args[["ylim"]])) {
      # EXTRACT YLIM - NAMES
      ymin <- as.numeric(
        unlist(
          strsplit(names(args[["fr_dens_list"]])[1], "-")
          )[1])
      ymax <- as.numeric(
        unlist(
          strsplit(names(args[["fr_dens_list"]])[args[["SMP"]]], "-")
          )[2])
      args[["ylim"]] <- c(ymin, ymax)
    }
    
  }
  
  # DO.CALL CYTO_PLOT_EMPTY
  do.call("cyto_plot_empty.list", args[names(args) %in% ARGS])
  
}

## .CYTO_PLOT_DENSITY ----------------------------------------------------------

#' Add density distributions to an existing cyto_plot
#'
#' @param args named list of cyto_plot_arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_density <- function(args) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REMOVE X ARGUMENT
  args <- args[-match("x", names(args))]

  # RENAME FR_LIST ARGUMENT TO X (DO.CALL() ON CYTO_PLOT_EMPTY)
  names(args)[match("fr_dens_list", names(args))] <- "x"

  # CYTO_PLOT_DENSITY ARGUMENTS
  ARGS <- formalArgs("cyto_plot_density.list")

  # RESTRICT SUPPLIED ARGUMENTS
  args <- args[names(args) %in% ARGS]

  # CALL CYTO_PLOT_DENSITY -----------------------------------------------------
  do.call("cyto_plot_density.list", args)
}

## .CYTO_PLOT_POINT ------------------------------------------------------------

#' Add points to an existing cyto_plot
#'
#' @param args named list of cyto_plot arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_point <- function(args) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # REMOVE X ARGUMENT
  args <- args[-match("x", names(args))]

  # RENAME FR_LIST ARGUMENT TO X (DO.CALL() ON CYTO_PLOT_EMPTY)
  names(args)[match("fr_list", names(args))] <- "x"

  # RESET DISPLAY ARGUMENT - DATA ALREADY SAMPLED
  args["display"] <- 1

  # CYTO_PLOT_DENSITY ARGUMENTS
  ARGS <- formalArgs("cyto_plot_point.list")

  # RESTRICT SUPPLIED ARGUMENTS
  args <- args[names(args) %in% ARGS]

  # CALL CYTO_PLOT_DENSITY -----------------------------------------------------
  do.call("cyto_plot_point.list", args)
}

## .CYTO_PLOT_GATE -------------------------------------------------------------

#' Add gates and labels to an existing cyto_plot
#'
#' @param args named list of cyto_plot arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_plot_gate <- function(args) {

  # GENERAL --------------------------------------------------------------------

  # SAMPLES
  SMP <- args[["SMP"]]
  
  # POPULATIONS PER LAYER
  NP <- args[["NP"]]
  
  # TOTAL POPULATIONS
  TNP <- args[["TNP"]]  
  
  # GATES
  NG <- length(args[["gate"]])

  # POPULATIONS PER GATE
  P <- .cyto_gate_count(args[["gate"]], 
                        negate = FALSE, 
                        total = FALSE)
  
  # LABEL
  label <- args[["label"]]

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # CYTO_PLOT_GATE & CYTO_PLOT_LABELLER ARGUMENTS
  gate_args <- formalArgs("cyto_plot_gate.list")
  label_args <- formalArgs("cyto_plot_labeller")

  # RESTRICT ARGUMENTS
  args <- args[names(args) %in% c(gate_args, label_args)]

  # SPLIT GATE_FILL ARGUMENTS BY POPULATIONS
  lapply(gate_args[which(grepl("gate_fill", gate_args))], function(z) {
    args[[z]] <<- split(args[[z]], rep(seq_len(NG), times = P))
  })

  # SPLIT LABEL ARGUMENTS BY POPULATIONS PER LAYER
  lapply(label_args, function(z) {
    args[[z]] <<- rep(args[[z]], length.out = TNP)
  })
  lapply(label_args, function(z) {
    args[[z]] <<- split(args[[z]], rep(seq_len(SMP), 
                                       each = NP))
  })

  # RE-ARRANGE LABEL COORDS PER GATE
  lapply(label_args, function(z) {
    args[[z]] <<- lapply(seq_len(NP), function(y) {
      LAPPLY(args[[z]], `[[`, y)
    })
  })

  # GATE ARGUMENTS
  gate_args <- args[names(args) %in% gate_args]
  
  # LABEL ARGUMENTS
  label_args <- args[names(args) %in% label_args]
  
  # GATE & ASSOCIATED LABELS ---------------------------------------------------

  # GATES & LABELS
  label_text_xy <- lapply(seq_len(NP), function(z) {
    # GATED POPULATION(S) - GATE & LABEL(S)
    if (z <= NG) {
      # PLOT GATE
      do.call(
        "cyto_plot_gate",
        c(
          list("channels" = gate_args[["channels"]]),
          lapply(
            gate_args[!grepl("channels", names(gate_args))],
            `[[`, z
          )
        )
      )
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (label == TRUE &
          !.all_na(args[["label_text"]][[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(args[["label_text"]][z])),
          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
      # NEGATED POPULATION(S) - LABEL(s) ONLY
    } else {
      # RETAIN MANUALLY SELECTED CO-ORDINATES
      if (label == TRUE &
          !.all_na(args[["label_text"]][[z]])) {
        # ADD LABELS
        text_xy <- do.call("cyto_plot_labeller", lapply(label_args, `[[`, z))
      } else {
        # NO LABELS
        text_xy <- matrix(rep(NA, 2 * length(args[["label_text"]][[z]])),
          ncol = 2
        )
        colnames(text_xy) <- c("x", "y")
      }
    }
    return(text_xy)
  })
  label_text_xy <- do.call("rbind", label_text_xy)

  # RE-ARRANGE LABEL ARGUMENTS -------------------------------------------------

  # UPDATE LABEL_TEXT_X & LABEL_TEXT_Y
  args[["label_text_x"]] <- label_text_xy[, "x"]
  args[["label_text_y"]] <- label_text_xy[, "y"]

  # REVERT LABEL_TEXT_X & LABEL_TEXT_Y TO ORIGINAL FORMAT
  if (SMP > 1) {
    args[["label_text_x"]] <- LAPPLY(seq_len(SMP), function(z) {
      args[["label_text_x"]][names(args[["label_text_x"]]) == z]
    })
    args[["label_text_y"]] <- LAPPLY(seq_len(SMP), function(z) {
      args[["label_text_y"]][names(args[["label_text_y"]]) == z]
    })
  }

  # RETURN LABEL CO-ORDINATES --------------------------------------------------

  # LABEL_COORDS MATRIX
  label_text_xy <- matrix(c(args[["label_text_x"]], args[["label_text_y"]]),
    ncol = 2,
    byrow = FALSE
  )
  colnames(label_text_xy) <- c("x", "y")
  return(label_text_xy)
}

## .CYTO_PLOT_LABEL ------------------------------------------------------------

#' Add labels to an existing cyto_plot
#'
#' @param args named list of cyto_plot arguments.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_plot_label <- function(args) {

  # PREPARE ARGUMENTS ----------------------------------------------------------

  # CYTO_PLOT_LABEL ARGUMENTS
  ARGS <- formalArgs("cyto_plot_labeller")

  # LABEL CONSTRUCTION ---------------------------------------------------------

  # CALL CYTO_PLOT_LABELLER
  label_text_xy <- do.call("cyto_plot_labeller", args[names(args) %in% ARGS])

  # RETURN LABEL CO-ORDINATES
  return(label_text_xy)
}
