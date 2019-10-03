# .CYTO_GATE_TYPE_CHECK --------------------------------------------------------

#' Check gate types supplied to cyto_gate_draw.
#'
#' @param type vector indicating the types of gates to construct using
#'   \code{gate_draw}.
#' @param channels names of channels used to construct the plot.
#' @param alias names of the populations to be gated.
#' @param negate logical flag indicating whether the negated population should
#'   be included in the gatingTemplate entry.
#'
#' @return Stop gating process if type is incorrect or returns \code{type} as
#'   full lower case name(s). If a single type is supplied for multiple
#'   populations, the same type will be used for all populations.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_gate_type <- function(type, channels, alias, negate = FALSE) {
  
  # DEFAULT GATE TYPES ---------------------------------------------------------
  
  # NO GATE TYPE SUPPLIED
  if(is.null(type)){
    # 1D PLOT - INTERVAL
    if(length(channels) == 1){
      type <- "interval"
    # 2D PLOT - POLYGON
    }else if(length(channels) == 2){
      type <- "polygon"
    }
  }
  
  # SUPPORTED GATE TYPES -------------------------------------------------------
  
  # GATE TYPES
  ind <- LAPPLY(seq_len(length(type)), function(z){
    if(grepl("rectangle", type[z], ignore.case = TRUE) | 
       (nchar(type[z]) == 1 & grepl("r", type[z], ignore.case = TRUE))){
      type[z] <<- "rectangle"
      return(TRUE)
    }else if(grepl("polygon", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("p", type[z], ignore.case = TRUE))){
      type[z] <<- "polygon"
      return(TRUE)
    }else if(grepl("interval", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("i", type[z], ignore.case = TRUE))){
      type[z] <<- "interval"
      return(TRUE)
    }else if(grepl("threshold", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("t", type[z], ignore.case = TRUE))){
      type[z] <<- "threshold"
      return(TRUE)
    }else if(grepl("boundary", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("b", type[z], ignore.case = TRUE))){
      type[z] <<- "boundary"
      return(TRUE)
    }else if(grepl("ellipse", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("e", type[z], ignore.case = TRUE))){
      type[z] <<- "ellipse"
      return(TRUE)
    }else if(grepl("quadrant", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("q", type[z], ignore.case = TRUE))){
      type[z] <<- "quadrant"
      return(TRUE)
    }else if(grepl("web", type[z], ignore.case = TRUE) |
             (nchar(type[z]) == 1 & grepl("w", type[z], ignore.case = TRUE))){
      type[z] <<- "web"
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  # UNSUPPORTED GATE TYPES
  if(!all(ind == TRUE)){
    if(length(!ind) == 1){
      stop(paste(type[!ind], "is not a valid gate type for cyto_gate_draw."))
    }else{
      stop(paste(type[ind], sep = " & ",
                 "are not valid gate types for cyto_gate_draw"))
    }
  }
  
  # PREPARE GATE TYPE ----------------------------------------------------------
  
  # REPEAT GATE TYPE ALIAS TIMES
  if(length(type) < length(alias)){
    type <- rep(type, length.out = length(alias))
    # ONLY SINGLE QUADRANT CALL
    if(length(grep("quadrant", type)) > 1){
      type <- type[-grep("quadrant", type)[-1]]
    } 
    # ONLY SINGLE WEB CALL
    if(length(grep("web", type)) > 1){
      type <- type[-grep("quadrant", type)[-1]]
    }
  }

  # UNSUPPORTED GATE TYPES 1D PLOTS --------------------------------------------
  
  # INTERVAL - BOUNDARY - THRESHOLD 
  if(length(channels) == 1){
    if(any(type %in% c("rectangle",
                       "polygon",
                       "ellipse",
                       "quadrant",
                       "web"))){
      stop("Supported 1D gate types include interval, boundary and threshold.")
    }
  }
  
  # REMOVE NEGATED GATE TYPE
  if(negate == TRUE){
    type <- type[-length(type)]
  }
  
  # RETURN VALID GATE TYPES ----------------------------------------------------
  
  return(type)
}

# .CYTO_ALIAS ------------------------------------------------------------------

#' Check supplied alias
#'
#' @param alias vector of population names.
#' @param type vector gate types to be used to gate the populations.
#' @param negate logical indicating whether the negated population
#'
#' @return stop gating process if alias is missing or of an incorrect length,
#'   and return alias.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_alias <- function(alias, 
                        type,
                        negate = FALSE){
  
  # CHECKS ---------------------------------------------------------------------
  
  # MISSING ALIAS
  if(is.null(alias)){
    stop("Supply the name(s) of the population(s) to 'alias'.")
  }
  
  # GATES ----------------------------------------------------------------------
  
  # GATED POPULATIONS
  N <- sum(LAPPLY(type, function(z){
    if(grepl("quadrant", z, ignore.case = TRUE) | 
       grepl("q", z, ignore.case = TRUE)){
      n <- 4 
    }else{
      n <- 1
    }
  }))
  
  # NEGATE
  if(negate == TRUE){
    N <- N + 1
  }
  
  # CHECK ALIAS ----------------------------------------------------------------
  
  # ALIAS LENGTH == N
  if(!length(alias) == N){
    stop("Supply a name for each of the population(s) to 'alias'.")
    if(negate == TRUE){
      stop("Have you forgotten to include a name for the negated population?")
    }
  }
  
  return(alias)
  
}