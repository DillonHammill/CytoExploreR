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
  
  # SPLIT GATE TYPES
  split_gate_types <- .split_gate_types()
  
  # MULTI GATE TYPES
  multi_gate_types <- c(split_gate_types, "quadrant")
  
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
  
  # CANNOT USED MIXED GATES FOR SPLIT GATE METHODS
  if(any(type %in% split_gate_types) & !all(type %in% split_gate_types)){
    stop(paste("Mixed gates are not supported for",
         paste(split_gate_types, sep = " & "), "gate types."))
  }
  
  # CANNOT NEGATE MULTI GATE METHODS
  if(negate == TRUE & any(type %in% multi_gate_types)){
    stop(paste("Cannot negate", 
               paste(multi_gate_types, sep = " & "),
               "gate_types."))
  }
  
  # PREPARE GATE TYPE ----------------------------------------------------------
  
  # REPEAT GATE TYPE ALIAS TIMES
  if(length(type) < length(alias)){
    type <- rep(type, length.out = length(alias))
    # ONLY SINGLE MULTI GATE TYPE CALL
    lapply(multi_gate_types, function(z){
      if(length(grep(z, type)) > 1){
        type <<- type[-grep(z, type)[-1]]
      }
    })
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
    type[length(type)] <- NA
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
#'   and return alias split by gate type.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_alias <- function(alias, 
                        type,
                        negate = FALSE,
                        split = FALSE){
  
  # CHECKS ---------------------------------------------------------------------
  
  # MISSING ALIAS
  if(is.null(alias)){
    stop("Supply the name(s) of the population(s) to 'alias'.")
  }
  
  # ALIAS PER GATE TYPE --------------------------------------------------------
  
  # EXPECTED ALIAS LENGTH PER GATE TYPE
  N <- LAPPLY(type, function(z){
    if(z == "quadrant"){
      n <- 4 
    }else if(z %in% .split_gate_types()){  
      n <- length(alias)
    }else{
      n <- 1
    }
    return(n)
  })
  
  # INCLUDE NEGATED POPULATION
  if(negate == TRUE){
    N <- c(N, 1)
  }
  
  # CHECK ALIAS ----------------------------------------------------------------
  
  # ALIAS LENGTH == N
  if(!length(alias) == sum(N)){
    stop("Supply a name for each of the population(s) to 'alias'.")
    if(negate == TRUE){
      stop("Have you forgotten to include a name for the negated population?")
    }
  }  
  
  # SPLIT ALIAS ----------------------------------------------------------------
  
  # PREPARE ALIAS
  alias <- split(alias, rep(seq_len(length(type)), times = N))
  return(alias)
  
}

#' @noRd
.split_gate_types <- function(){
  c("web")
}
