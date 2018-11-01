#' Removes any observation from the given flowFrame object that has values
#' outside the given range for the specified channels
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param min a numeric vector that sets the lower bounds for data filtering
#' @param max a numeric vector that sets the upper bounds for data filtering
#' @param channels \code{character} specifying which channel to operate on
#' @return a \code{flowFrame} object
#' @examples
#' \dontrun{
#'  library(flowClust)
#'  data(rituximab)
#'  # Consider the range of values for FSC.H and SSC.H
#'  summary(rituximab)
#' 
#'  # Truncates any observations with FSC.H outside [100, 950]
#'  rituximab2 <- .truncate_flowframe(rituximab, channels = "FSC.H", min = 100, max = 950)
#'  summary(rituximab2)
#'  # Next, truncates any observations with SSC.H outside [50, 1000]
#'  rituximab3 <- .truncate_flowframe(rituximab2, channels = "SSC.H", min = 50, max = 1000)
#'  summary(rituximab3)
#'
#'  # Instead, truncates both channels at the same time
#'  rituximab4 <- .truncate_flowframe(rituximab, channels = c("FSC.H", "SSC.H"),
#'  min = c(100, 50), max = c(950, 1000))
#'  summary(rituximab4)
#' }
#' 
#' @noRd 
.truncate_flowframe <- function(flow_frame, channels, min = NULL, max = NULL) {
  channels <- as.character(channels)
  num_channels <- length(channels)
  
  # For comparison purposes, we update the min and max values to -Inf and Inf,
  # respectively, if NULL.
  if (is.null(min)) {
    min <- rep(-Inf, num_channels)
  }
  if (is.null(max)) {
    max <- rep(Inf, num_channels)
  }
  
  if (!(num_channels == length(min) && num_channels == length(max))) {
    stop("The lengths of 'min' and 'max' must match the number of 'channels' given.")
  }
  
  gate_coordinates <- lapply(seq_len(num_channels), function(i) {
    c(min[i], max[i])
  })
  names(gate_coordinates) <- channels
  
  truncate_filter <- rectangleGate(gate_coordinates)
  Subset(flow_frame, truncate_filter)
}