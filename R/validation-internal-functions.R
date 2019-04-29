## Function Definition for Internal Checks -------------------------------------

# Validation Functions for gate_draw -------------------------------------------

#' Check Gate Type(s) Supplied to gate_draw.
#'
#' @param type vector indicating the types of gates to construct using
#'   \code{gate_draw}.
#' @param alias names of the populations to be gated.
#'
#' @return Stop gating process if type is incorrect or returns \code{type} as
#'   full lower case name(s). If a single type is supplied for multiple
#'   populations, the same type will be used for all populations.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .cyto_gate_type_check(type = "r", alias = c("A", "B", "C"))
#' @noRd
.cyto_gate_type_check <- function(type, alias) {
  if (all(type %in% c("q", "Q", "quadrant", "Quadrant")) & 
      length(alias) != 4) {
    stop("Supply the names of 4 poulations to alias for quadrant gates.")
  }
  
  gts <- c(
    "polygon",
    "Polygon",
    "p",
    "P",
    "rectangle",
    "Rectangle",
    "r",
    "R",
    "interval",
    "Interval",
    "i",
    "I",
    "threshold",
    "Threshold",
    "t",
    "T",
    "boundary",
    "Boundary",
    "b",
    "B",
    "ellipse",
    "Ellipse",
    "e",
    "E",
    "quadrant",
    "Quadrant",
    "q",
    "Q",
    "web",
    "Web",
    "w",
    "W"
  )
  
  if (!all(type %in% gts)) {
    if (length(type[type %in% gts == FALSE]) >= 2) {
      stop(
        paste(
          paste(type[type %in% gts == FALSE], collapse = " & "),
          "are not valid gate types for gate_draw!"
        )
      )
    } else {
      stop(paste(
        type[type %in% gts == FALSE],
        "is not a valid gate type for gate_draw!"
      ))
    }
  }
  
  type[type %in% c("polygon", "Polygon", "p", "P")] <- "polygon"
  
  type[type %in% c("rectangle", "Rectangle", "r", "R")] <- "rectangle"
  
  type[type %in% c("interval", "Interval", "i", "I")] <- "interval"
  
  type[type %in% c("threshold", "Threshold", "t", "T")] <- "threshold"
  
  type[type %in% c("boundary", "Boundary", "b", "B")] <- "boundary"
  
  type[type %in% c("ellipse", "Ellipse", "e", "E")] <- "ellipse"
  
  type[type %in% c("quadrant", "Quadrant", "q", "Q")] <- "quadrant"
  
  type[type %in% c("web", "Web", "w", "W")] <- "web"
  
  # Repeat type to equal length of alias
  if (length(type) != length(alias) &
      type[1] != "quadrant" &
      type[1] != "web") {
    type <- rep(type, length(alias))
  }
  
  return(type)
}

#' Check Alias Supplied to gate_draw
#'
#' @param alias vector indicating the names of the populations to be gated.
#' @param type vector indicating the type(s) of gate(s) to be constructed.
#'
#' @return Stops the gating process if alias is missing or of the incorrect
#'   length given the gate type.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .cyto_alias_check(alias = c("A", "B", "C", "D"), type = "q")
#' @noRd
.cyto_alias_check <- function(alias = NULL, type) {
  if (is.null(alias)) {
    stop("Supply names of populations to 'alias' for checking.")
  }
  
  if (type[1] == "quadrant" & length(alias) != 4) {
    stop("Supply 4 population names to 'alias' to construct quadrant gates.")
  }
  
  if (length(type) > 1) {
    if (length(alias) != length(type)) {
      stop("Length of alias must be the same length as type for multi-gates.")
    }
  }
}

# Validation Functions for gatingTemplate Objects ------------------------------

#' Check gatingTemplate for Existing Entry
#'
#' @param parent name of the parent population.
#' @param alias name of the population of interest.
#' @param gatingTemplate csv file name of the gatingTemplate.
#'
#' @return stops the gating process if an entry already exists in the
#'   gatingTemplate for the supplied alias.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom utils read.csv
#'
#' @noRd
.cyto_gatingTemplate_check <- function(parent, alias, gatingTemplate) {
  if (inherits(gatingTemplate, "gatingTemplate")) {
    stop("'gatingTemplate' should be the name of the gatingTemplate csv file.")
  } else {
    if (getOption("CytoRSuite_wd_check") == TRUE) {
      if (.file_wd_check(gatingTemplate)) {
        gt <- read.csv(gatingTemplate, header = TRUE)
        
        # Parent and alias entries match file
        if (any(gt$parent %in% parent & gt$alias %in% alias)) {
          message(
            paste(
              paste(gt$alias, collapse = " & "),
              "already exists in", gatingTemplate, "."
            )
          )
          stop("Supply another gatingTemplate or edit gate(s) using gate_edit.")
        }
      }
    }
  }
}

# Validation Functions for Statistical Methods ---------------------------------

#' Check Statistic for cyto_stats_compute
#'
#' @param stat cyto_stats_compute statistic.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_stat_check <- function(stat) {
  if (!stat %in% c(
    "mean",
    "Mean",
    "median",
    "Median",
    "mode",
    "Mode",
    "count",
    "Count",
    "events",
    "Events",
    "percent",
    "Percent",
    "freq",
    "Freq",
    "geo mean",
    "Geo mean",
    "Geo Mean",
    "CV",
    "cv"
  )) {
    stop("Supplied statistic not supported.")
  }
  
  if (stat %in% c("mean", "Mean")) {
    stat <- "mean"
  } else if (stat %in% c("median", "Median")) {
    stat <- "median"
  } else if (stat %in% c("mode", "Mode")) {
    stat <- "mode"
  } else if (stat %in% c("count", "Count","events", "Events")) {
    stat <- "count"
  } else if (stat %in% c("percent", "Percent", "freq", "Freq")) {
    stat <- "freq"
  } else if (stat %in% c("geo mean", "Geo mean", "Geo Mean")) {
    stat <- "geo mean"
  } else if (stat %in% c("cv", "CV")) {
    stat <- "CV"
  }
  
  return(stat)
}
