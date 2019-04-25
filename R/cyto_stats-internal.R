## Internal Statistical Functions ----------------------------------------------

# In order to calculate statistics on the original linear scale a transformation
# object (either of class transformList or transformerList) must be supplied.
# Otherwise the statistic will be calculated on the transformed scale.

# To add new statistical functions to cyto_stats_compute follow the same format
# as below and add the name of the function to cyto_stats_compute options.

#' Calculate Number of Events
#' 
#' @param x object of class flowFrame.
#'
#' @return data.frame containing the number events.
#'
#' @importFrom flowCore exprs
#' @importFrom tibble tibble
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
.cyto_count <- function(x){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Assign x to fr
  fr <- x
  
  # Extract raw data and calculate mean directly
  res <- tibble("count" = BiocGenerics::nrow(fr))
  
  return(res)
  
}

#' Calculate Arithmetic Mean Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#'
#' @return data.frame containing the mean fluorescent intensity for each of
#'   the supplied channels.
#'
#' @importFrom flowCore exprs
#' @importFrom tibble as_tibble
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_mean <- function(x, 
                       channels = NULL,
                       trans = NULL){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformList/transformerList needed for transformed chans
  if(is.null(trans)){
    message(
      paste(
        "'trans' requires a transformList/transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Assign x to fr
  fr <- x
  
  # Check channels
  channels <- cyto_channel_check(fr, channels, FALSE)
  
  # Inverse transformations
  if(!is.null(trans)){
    inv <- cyto_trans_check(trans = trans, inverse = TRUE)
    fr <- transform(fr, inv)
  }
  
  # Extract raw data and calculate mean directly - colMeans for speed
  res <- colMeans(exprs(fr)[, channels])
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_marker_extract(fr, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

#' Calculate Geometric Mean Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#'
#' @return data.frame containing the geometric mean fluorescent intensity for
#'   each of the supplied channels.
#'
#' @importFrom flowCore exprs
#' @importFrom tibble as_tibble
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_geometric_mean <- function(x, 
                                 channels = NULL,
                                 trans = NULL){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Assign x to fr
  fr <- x
  
  # Check channels
  channels <- cyto_channel_check(fr, channels, FALSE)
  
  # Geometric mean calculate as inverse of arithmetic mean of transformed data
  if(is.null(trans)){
    
    log_exprs <- suppressWarnings(log(exprs(fr)[,channels]))
    geo_mean <- suppressWarnings(exp(colMeans(log_exprs)))
      
    # Zero/negative values present
    if(any(is.nan(geo_mean))){
      stop(
        paste0(
          "Supply transformList/transformerList to calculate geometric mean."
        )
      )
    }
    
  }else if(!is.null(trans)){
    
    # Convert tranform object to transformList
    trans <- cyto_trans_check(trans, inverse = FALSE)
    
    fr_mean <- colMeans(exprs(fr)[,channels])
    
    res <- unlist(lapply(channels, function(x){
      
      # Channel has been transformed
      if(x %in% BiocGenerics::colnames(trans)){
        
        # Inverse transformations
        inv <- cyto_trans_check(trans, inverse = TRUE)
        
        # Inverse transformation on calculated arithmetic mean
        geo_mean <- inv@transforms[[x]]@f(fr_mean[x])
        
      # Channel has not been transformed  
      }else{
        
        geo_mean <- exp(mean(log(exprs(fr)[,x])))
        
      }
      
      return(geo_mean)
      
    }))
    
  }
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_marker_extract(fr, channels)
  res <- as_tibble(res)

  return(res)
  
}

#' Calculate Median Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#'
#' @return data.frame containing the median fluorescent intensity for each of
#'   the supplied channels.
#'
#' @importFrom flowCore exprs
#' @importFrom tibble as_tibble
#' @importFrom robustbase colMedians
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_median <- function(x, 
                        channels = NULL,
                        trans = NULL){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformList/transformerList needed for transformed chans
  if(is.null(trans)){
    message(
      paste(
        "'trans' requires a transformList/transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Assign x to fr
  fr <- x
  
  # Check channels
  channels <- cyto_channel_check(fr, channels, FALSE)
  
  # Inverse transformations
  if(!is.null(trans)){
    inv <- cyto_trans_check(trans = trans, inverse = TRUE)
    fr <- transform(fr, inv)
  }
  
  # Extract raw data and calculate median directly - colMedians for speed
  res <- colMedians(exprs(fr)[,channels])
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_marker_extract(fr, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

#' Calculate Mode Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#' @param density_smooth numeric smoothing factor passed to
#'   \code{\link[stats:density]{density}} to control the degree of smoothing,
#'   set to 1.5 by default.
#'
#' @return data.frame containing the mode fluorescent intensity for each of
#'   the supplied channels.
#'
#' @importFrom flowCore exprs
#' @importFrom tibble as_tibble
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_mode <- function(x, 
                       channels = NULL,
                       trans = NULL,
                       density_smooth = 1.5){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformList/transformerList needed for transformed chans
  if(is.null(trans)){
    message(
      paste(
        "'trans' requires a transformList/transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Assign x to fr
  fr <- x
  
  # Check channels
  channels <- cyto_channel_check(fr, channels, FALSE)
  
  # Inverse transformations
  if(!is.null(trans)){
    inv <- cyto_trans_check(trans = trans, inverse = TRUE)
    fr <- transform(fr, inv)
  }
  
  # Extract raw data and calculate mode directly
  res <- unlist(lapply(channels, function(x){
    d <- .cyto_density(fr, 
                       channel = x, 
                       density_smooth = density_smooth, 
                       modal = FALSE)
    d$x[d$y == max(d$y)]
  }))
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_marker_extract(fr, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

#' Calculate Robust Coefficient of Variation
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#'
#' @return data.frame containing the robust coefficient of variation for each
#'   of the supplied channels.
#'
#' @importFrom flowCore exprs
#' @importFrom tibble as_tibble
#' @importFrom robustbase colMedians
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_CV <- function(x, 
                    channels = NULL,
                    trans = NULL){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformList/transformerList needed for transformed chans
  if(is.null(trans)){
    message(
      paste(
        "'trans' requires a transformList/transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Assign x to fr
  fr <- x
  
  # Check channels
  channels <- cyto_channel_check(fr, channels, FALSE)
  
  # Inverse transformations
  if(!is.null(trans)){
    inv <- cyto_trans_check(trans = trans, inverse = TRUE)
    fr <- transform(fr, inv)
  }
  
  # Extract raw data and calculate CV directly
  fr_median <- colMedians(exprs(fr)[,channels])
  res <- unlist(lapply(channels, function(x){
    md <- fr_median[x]
    rSD <- median(abs(exprs(fr)[,x] - md)) * 1.4826
    rSD/md * 100
  }))
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_marker_extract(fr, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

#' Calculate Kernel Density
#'
#' @param x object of class flowFrame.
#' @param channel channel to calculate kernel density.
#' @param density_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}}.
#' @param modal logical indicating whether densities should be normalised to
#'   mode.
#'
#' @return calculated kernel density for supplied channel.
#'
#' @importFrom flowCore exprs
#' @importFrom stats density
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_density <- function(x,
                          channel = NULL,
                          density_smooth = 1.5,
                          modal = TRUE) {
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Assign x to fr
  fr <- x
  
  # Check channels
  channel <- cyto_channel_check(fr, channel, FALSE)
  
  # Calculate kernel density
  dens <- density(exprs(fr)[,channel],
                  adjust = density_smooth)
  
  # Normalise to mode if required
  if(modal){
    dens$y <- (dens$y / max(dens$y)) * 100
  }
    
  return(dens)
}

#' Calculate combined range of cytometry objects
#'
#' The lower limit is always set to zero unless there is data below this limit.
#' The upper limit is determined by the limits argument.
#'
#' @param x cytometry object(s) which require range calculation.
#' @param channels name(s) of channels.
#' @param limits either "data" or "machine".
#' @param plot logical indicating whether a check should be performed for
#'   channel length.
#' @param buffer FALSE or fraction indcating the amount of buffering to be added
#'   to the lower limit for transformed channels, set to 0.1 by default.
#'
#' @return vector containing minimum and maximum values.
#'
#' @importFrom flowCore flowSet
#' @importFrom flowWorkspace getData
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_range <- function(x,
                        parent = NULL,
                        channels = NULL,
                        limits = "data",
                        plot = FALSE,
                        buffer = 0.1) {
  
  # Time parameter always uses data limits
  if("Time" %in% channels){
    limits <- "data"
  }
  
  # Convert GatingHierarchy/GatingSet to flowFrame/flowSet
  if(inherits(x, "GatingHierarchy") | inherits(x, "GatingSet")){
    x <- getData(x, parent)
  }
  
  # Convert lists to conventional objects
  if(class(x) == "list"){
    
    # Convert list of flowFrames to flowSet
    if(all(unlist(lapply(x,"class")) == "flowFrame")){
      x <- .cyto_convert(x, "flowSet")
    }
    
    # Convert list of flowSets to flowSet
    if(all(unlist(lapply(x, function(y){inherits(y,"flowSet")})))){
      x <- .cyto_convert(x, "flowSet")
    }
    
  }
  
  # Convert markers to channels
  if(!is.null(channels)){
    channels <- cyto_channel_check(x, channels = channels, plot = plot)
  }else{
    channels <- BiocGenerics::colnmaes(x)
  }
  
  # Extract ranges
  if(inherits(x, "flowFrame")){
    
    # Lower bound - use data limits
    rng <- suppressWarnings(range(x, type = "data")[,channels])
    
    # Upper bound
    if(limits == "data"){
      mx <- suppressWarnings(range(x, type = "data")[,channels][2,])
    }else if(limits == "machine"){
      mx <- suppressWarnings(range(x, type = "instrument")[,channels][2,])
    }
    rng[2,] <- mx
    
  }else if(inherits(x, "flowSet")){
    
    # Lower bound - use data limits
    rng <- suppressWarnings(range(x[[1]], type = "data")[,channels])
    
    # Upper bound
    if(limits == "data"){
      mx <- do.call("rbind",fsApply(x, function(fr){
        suppressWarnings(range(fr, type = "data")[,channels][2,])
      }))
      mx <- unlist(lapply(channels, function(chan){
        mx <- mx[,chan]
        mx <- mx[is.finite(mx)]
        max(mx)
      }))
    }else if(limits == "machine"){
      mx <- do.call("rbind", fsApply(x, function(fr){
        suppressWarnings(range(fr, type = "instrument")[,channels][2,])
      }))
      mx <- unlist(lapply(channels, function(chan){
        mx <- mx[,chan]
        mx <- mx[is.finite(mx)]
        max(mx)
      }))
    }
    rng[2,] <- mx
    
  }
  
  # Replace lower data limit if > 0
  if(any(rng[1,] > 0)){
    rng[1, rng[1,] > 0] <- 0
  }
  
  # Add buffer to lower limit if channel is transformed
  if(buffer != FALSE){
    lapply(channels, function(chan){
      if(rng[,chan][2] < 6){
        rng[,chan][1] <<- rng[,chan][1] - buffer*(rng[,chan][2] - rng[,chan][1])
      }
    })
  }
  
  return(rng)
}
