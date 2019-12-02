## Internal Statistical Functions ----------------------------------------------

# In order to calculate statistics on the original linear scale a transformation
# object (either of class transformList or transformerList) must be supplied.
# Otherwise the statistic will be calculated on the transformed scale.

# To add new statistical functions to cyto_stats_compute follow the same format
# as below and add the name of the function to cyto_stats_compute options.

# Gates are applied to flowFrame prior to statistics computation.

## COUNT -----------------------------------------------------------------------

#' Calculate Number of Events
#' 
#' @param x object of class flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#'
#' @return data.frame containing the number events.
#'
#' @importFrom flowCore Subset
#' @importFrom tibble tibble
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
.cyto_count <- function(x,
                        gate = NA){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }
  }
  
  # Gate flowFrame
  if(!.all_na(gate)){
    x <- Subset(x, gate)
  }
  
  # Return count statistic
  res <- tibble("count" = BiocGenerics::nrow(x))
  
  return(res)
  
}

## FREQUENCY -------------------------------------------------------------------

#' Calculate frequency of gated population
#'
#' @param x object of class flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#'
#' @return data.frame containing the population frequency.
#'
#' @importFrom flowCore Subset
#' @importFrom tibble tibble
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_freq <- function(x,
                       gate = NA){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Must supply gate to compute frequency
  if(.all_na(gate)){
    stop("A 'gate' object is required to compute frequency.")
  }
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }else{
      # Gating with rectangleGate objects is slow ...
      if(inherits(gate, "rectangleGate")){
        gate <- as(gate, "polygonGate")
      }
    }
  }
  
  # Extract population
  fr <- Subset(x, gate)
  
  # Compute frequency
  res <- .cyto_count(fr)/.cyto_count(x) * 100
  res <- round(res, 2)
  
  # Return calculated frequency
  res <- tibble("freq" = as.numeric(res))
  
  return(res)
  
}


## MEAN ------------------------------------------------------------------------

#' Calculate Arithmetic Mean Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#'
#' @return data.frame containing the mean fluorescent intensity for each of
#'   the supplied channels.
#'
#' @importFrom flowCore exprs Subset
#' @importFrom tibble as_tibble
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_mean <- function(x, 
                       channels = NULL,
                       trans = NA,
                       gate = NA){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformList/transformerList needed for transformed chans
  if(.all_na(trans)){
    message(
      paste(
        "'trans' requires a transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Check channels
  channels <- cyto_channels_extract(x, channels, FALSE)
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }else{
      # Gating with rectangleGate objects is slow ...
      if(inherits(gate, "rectangleGate")){
        gate <- as(gate, "polygonGate")
      }
    }
  }
  
  
  # Gate flowFrame
  if(!.all_na(gate)){
    x <- Subset(x, gate)
  }
  
  # Get raw data
  if(!.all_na(trans)){
    x <- cyto_transform(x, trans = trans, inverse = TRUE, plot = FALSE)
  }
  
  # Extract raw data and calculate mean directly - colMeans for speed
  res <- colMeans(exprs(x)[, channels, drop = FALSE])
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_markers_extract(x, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

## GEOMETRIC MEAN --------------------------------------------------------------

#' Calculate Geometric Mean Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#'   
#'
#' @return data.frame containing the geometric mean fluorescent intensity for
#'   each of the supplied channels.
#'
#' @importFrom flowCore exprs Subset
#' @importFrom tibble as_tibble
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_geometric_mean <- function(x, 
                                 channels = NULL,
                                 trans = NA,
                                 gate = NA){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }

  # Check channels
  channels <- cyto_channels_extract(x, channels, FALSE)
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }else{
      # Gating with rectangleGate objects is slow ...
      if(inherits(gate, "rectangleGate")){
        gate <- as(gate, "polygonGate")
      }
    }
  }
  
  
  # Gate flowFrame
  if(!.all_na(gate)){
    x <- Subset(x, gate)
  }
  
  # Geometric mean calculate as inverse of arithmetic mean of transformed data
  if(.all_na(trans)){
    
    log_exprs <- suppressWarnings(log(exprs(x)[,channels,drop = FALSE]))
    geo_mean <- suppressWarnings(exp(colMeans(log_exprs)))
      
    # Zero/negative values present
    if(any(is.nan(geo_mean))){
      stop(
        paste0(
          "Supply transformerList to calculate geometric mean."
        )
      )
    }
  
  # transformerList supplied
  }else if(!.all_na(trans)){
    
    # Calculate arithmetic mean on transformed scale
    fr_mean <- colMeans(exprs(x)[,channels, drop = FALSE])
    
    # Inverse result to linear scale for geometric mean
    res <- LAPPLY(channels, function(z){
      
      # Channel has been transformed
      if(z %in% names(trans)){
        
        # Inverse transformations
        inv <- cyto_transform_extract(trans, inverse = TRUE)
        
        # Inverse transformation on calculated arithmetic mean
        geo_mean <- inv@transforms[[z]]@f(fr_mean[z])
        
      # Channel has not been transformed  
      }else{
        
        geo_mean <- exp(mean(log(exprs(x)[,z, drop = FALSE])))
        
      }
      
      return(geo_mean)
      
    })
    
  }
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_markers_extract(x, channels)
  res <- as_tibble(res)

  return(res)
  
}

## MEDIAN ----------------------------------------------------------------------

#' Calculate Median Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#'
#' @return data.frame containing the median fluorescent intensity for each of
#'   the supplied channels.
#'
#' @importFrom flowCore exprs Subset
#' @importFrom tibble as_tibble
#' @importFrom robustbase colMedians
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_median <- function(x, 
                        channels = NULL,
                        trans = NA,
                        gate = NA){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformerList needed for transformed chans
  if(.all_na(trans)){
    message(
      paste(
        "'trans' requires a transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Check channels
  channels <- cyto_channels_extract(x, channels, FALSE)
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }else{
      # Gating with rectangleGate objects is slow ...
      if(inherits(gate, "rectangleGate")){
        gate <- as(gate, "polygonGate")
      }
    }
  }
  
  # Gate flowFrame
  if(!.all_na(gate)){
    x <- Subset(x, gate)
  }
  
  # Get raw data
  if(!.all_na(trans)){
    x <- cyto_transform(x,
                        trans = trans,
                        inverse = TRUE, 
                        plot = FALSE)
  }
  
  # Extract raw data and calculate median directly - colMedians for speed
  res <- colMedians(exprs(x)[,channels,drop = FALSE])
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_markers_extract(x, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

## MODE ------------------------------------------------------------------------

#' Calculate Mode Fluorescent Intensity
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#' @param density_smooth numeric smoothing factor passed to
#'   \code{\link[stats:density]{density}} to control the degree of smoothing,
#'   set to 1.5 by default.
#'
#' @return data.frame containing the mode fluorescent intensity for each of
#'   the supplied channels.
#'
#' @importFrom flowCore exprs Subset
#' @importFrom tibble as_tibble
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_mode <- function(x, 
                       channels = NULL,
                       trans = NA,
                       gate = NA,
                       density_smooth = 0.6){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformerList needed for transformed chans
  if(.all_na(trans)){
    message(
      paste(
        "'trans' requires a transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }
  
  # Check channels
  channels <- cyto_channels_extract(x, channels, FALSE)
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }else{
      # Gating with rectangleGate objects is slow ...
      if(inherits(gate, "rectangleGate")){
        gate <- as(gate, "polygonGate")
      }
    }
  }
  
  # Gate flowFrame
  if(!.all_na(gate)){
    x <- Subset(x, gate)
  }
  
  # Get raw data
  if(!.all_na(trans)){
    x <- cyto_transform(x, trans = trans, inverse = TRUE, plot = FALSE)
  }
  
  # Extract raw data and calculate mode directly
  res <- LAPPLY(channels, function(z){
    d <- .cyto_density(x, 
                       channel = z, 
                       smooth = density_smooth, 
                       modal = FALSE)
    d$x[d$y == max(d$y)]
  })
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_markers_extract(x, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

## COEFFICIENT OF VARIATION ----------------------------------------------------

#' Calculate Robust Coefficient of Variation
#'
#' @param x object of class flowFrame.
#' @param channels name(S) of the channels for which the mean fluorescent
#'   intensity should be calculated.
#' @param trans transformation object used to transform the channels of the
#'   flowFrame.
#' @param gate object of class \code{rectangleGate}, \code{polygonGate} or
#'   \code{ellipsoidGate}.
#'
#'
#' @return data.frame containing the robust coefficient of variation for each
#'   of the supplied channels.
#'
#' @importFrom flowCore exprs Subset
#' @importFrom tibble as_tibble
#' @importFrom robustbase colMedians
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_CV <- function(x, 
                    channels = NULL,
                    trans = NA,
                    gate = NA){
  
  # Throw error for invalid object
  if(!inherits(x, "flowFrame")){
    stop("'x' should be a flowFrame object.")
  }
  
  # Reminder that transformerList needed for transformed chans
  if(.all_na(trans)){
    message(
      paste(
        "'trans' requires a transformerList to calculate",
        "statistics on a linear scale for transformed channels."
      )
    )
  }

  # Check channels
  channels <- cyto_channels_extract(x, channels, FALSE)
  
  # Only single gate objects are supported - calculate stats separately
  if(!.all_na(gate)){
    if(!any(c(inherits(gate, "rectangleGate"),
              inherits(gate, "polygonGate"),
              inherits(gate, "ellipsoidGate")))){
      stop(
        paste("Only rectangleGate, polygonGate and ellipsoidGate objects are",
              "supported.")
      )
    }else{
      # Gating with rectangleGate objects is slow ...
      if(inherits(gate, "rectangleGate")){
        gate <- as(gate, "polygonGate")
      }
    }
  }
  
  
  # Gate flowFrame
  if(!.all_na(gate)){
    x <- Subset(x, gate)
  }
  
  # Get raw data
  if(!.all_na(trans)){
    x <- cyto_transform(x, trans = trans, inverse = TRUE, plot = FALSE)
  }
  
  # Extract raw data and calculate CV directly
  fr_median <- colMedians(exprs(x)[,channels,drop = FALSE])
  res <- LAPPLY(channels, function(z){
    md <- fr_median[z]
    rSD <- median(abs(exprs(x)[,z] - md)) * 1.4826
    rSD/md * 100
  })
  
  # return transposed tibble for cbinding to pData
  res <- t(res)
  colnames(res) <- cyto_markers_extract(x, channels)
  res <- as_tibble(res)
  
  return(res)
  
}

## DENSITY ---------------------------------------------------------------------

#' Calculate Kernel Density
#'
#' @param x object of class flowFrame.
#' @param channel channel to calculate kernel density.
#' @param density_smooth smoothing parameter passed to
#'   \code{\link[stats:density]{density}}.
#' @param modal logical indicating whether densities should be normalised to
#'   mode.
#' @param ... additional arguments passed to
#'   \code{\link[stats:density]{density}}.
#'
#' @return calculated kernel density for supplied channel.
#'
#' @importFrom flowCore exprs
#' @importFrom stats density
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_density <- function(x, ...){
  UseMethod(".cyto_density")
}

#' @noRd
.cyto_density.flowFrame <- function(x,
                                    channel = NULL,
                                    smooth = 0.6,
                                    modal = TRUE, ...){
  
  # Check channels
  channel <- cyto_channels_extract(x, channel, FALSE)
  
  # Calculate kernel density
  if(is.null(length(exprs(x)[,channel]))){
    warning("Insufficient events to compute kernel density.")
    return(NA)
  }else{
    if(length(exprs(x)[,channel]) > 2){
      dens <- suppressWarnings(density(exprs(x)[,channel],
                                       adjust = smooth,
                                       ...))
    }else{
      warning("Insufficient events to compute kernel density.")
      return(NA)
    }
  }
  
  # Normalise to mode if required
  if(modal){
    dens$y <- (dens$y / max(dens$y)) * 100
  }
  
  return(dens)
  
}

#' Save y range to each layer
#' Easy to get ylim & label y co-ordinates
#' @importFrom stats bw.nrd0
#' @noRd
.cyto_density.list <- function(x,
                               channel = NULL,
                               smooth = 0.6,
                               modal = TRUE,
                               stack = 0.5,
                               layers = length(x)) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # LIST OF FLOWFRAMES
  if(!all(LAPPLY(x, "class") == "flowFrame")){
    stop("'x' must be a list of flowFrame objects.")
  }
  
  # LAYERS - EQUAL
  if(!length(x) %% layers == 0){
    stop("Number of flowFrames must be divisible by 'layers'.")
  }
  
  # GENERAL --------------------------------------------------------------------
  
  # SAMPLES
  smp <- length(x)
  
  # OVERLAY
  ovn <- smp - 1
  
  # BANDWIDTH ------------------------------------------------------------------
  
  # POOL DATA
  dt <- LAPPLY(x, function(z){exprs(z)[, channel]})
  
  # RESTRICT TO 100000 EVENTS
  if(length(dt) > 100000){
    dt <- sample(dt, 100000)
  }
  
  # JOINT BANDWIDTH
  bw <- bw.nrd0(dt)
  
  # KERNEL DENSITY -------------------------------------------------------------
  
  # CALL FLOWFRAME METHOD
  fr_dens_list <- lapply(x, function(z){
    suppressWarnings(.cyto_density(z,
                  channel = channel,
                  smooth = smooth,
                  modal = modal,
                  bw = bw))
  })
  names(fr_dens_list) <- rep(paste(0, 
                               max(fr_dens_list[[1]]$y),
                               sep = "-"),
                             length.out = length(fr_dens_list))
  
  # STACKING -------------------------------------------------------------------
  
  # STACKING REQUIRED
  if(stack != 0 & ovn != 0){
    
    # YRANGE PER LAYER
    if(!.all_na(fr_dens_list)){
      y_range <- mean(LAPPLY(fr_dens_list, function(d){
        if(!.all_na(d)){
          max(d$y)
        }else{
          NA
        }
      }), na.rm = TRUE)
    }else{
      y_range <- 100
    }

    # STACKING LEVELS
    stk <- seq(0, 
               smp * stack * y_range,
               stack * y_range)
    stk <- stk[seq_len(layers)]
    
    # LAYERS
    stk <- rep(stk, smp/layers)
    
    # APPLY STACKING 
    lapply(seq_len(smp), function(z){
      if(!.all_na(fr_dens_list[[z]])){
        fr_dens_list[[z]]$y <<- fr_dens_list[[z]]$y + stk[z]
      }
    })
    
    # YLIM PER LAYER
    ymin <- stk[seq_len(smp)]
    ymax <- ymin + y_range
    
    # SAVE YLIM - NAMES
    names(fr_dens_list) <- paste(ymin, ymax, sep = "-")
    
  }
  
  # LAYERS ---------------------------------------------------------------------
  
  # RETURN LIST OF DENSITY LISTS
  if(layers != smp){
    ind <- rep(seq_len(smp/layers), each = layers)
    fr_dens_list <- split(fr_dens_list, ind)
  }
  
  # RETURN STACKED DENSITY -----------------------------------------------------
  
  # LIST OF DENSITY OR LIST OF DENSITY LISTS
  return(fr_dens_list)
  
}

## DENSITY RANGES --------------------------------------------------------------

#' Calculate range of density objects
#' 
#' @param x list of density objects.
#' @param axis either "x" or "y".
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @noRd
.cyto_density_range <- function(x,
                                axis = "x"){
  
  if(!all(LAPPLY(x,"class") == "density")){
    stop("'x should be a list of density objects")
  }
  
  # Range
  rng <- lapply(x, function(d){
    if(axis == "x"){
      range(d$x)
    }else if(axis == "y"){
      range(d$y)
    }
  })
  rng <- do.call(rng, "rbind")
  rng <- c(min(rng[,1]), max(rng[,2]))
  
  return(rng)
}

## RANGE -----------------------------------------------------------------------

#' Calculate combined range of cytometry objects
#'
#' The lower limit is always set to zero unless there is data below this limit.
#' The upper limit is determined by the limits argument.
#'
#' @param x cytometry object(s) which require range calculation.
#' @param parent name of parent population to extract for range calculation.
#' @param channels name(s) of channel(s).
#' @param limits either "auto", "data" or "machine". "auto" use data limits but
#'   always includes zero.
#' @param plot logical indicating whether a check should be performed for
#'   channel length.
#' @param buffer fraction indcating the amount of buffering to be added on top
#'   of the upper and lower limit, set to 0.03 by default.
#' @param anchor logical indicating if the lower limit should be anchored to
#'   zero if the data range is above zero, set to TRUE by default.
#'
#' @return vector containing minimum and maximum values.
#'
#' @importFrom flowCore flowSet fsApply
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_range <- function(x, ...){
  UseMethod(".cyto_range")
}

#' @noRd
.cyto_range.flowFrame <- function(x, 
                                  channels = NA,
                                  limits = "auto",
                                  plot = FALSE,
                                  buffer = 0.03,
                                  anchor = TRUE,
                                  ...){
  
  # flowCore compatibility
  if(limits == "instrument"){
    limits <- "machine"
  }
  
  # Convert markers to channels
  if(!.all_na(channels)){
    channels <- cyto_channels_extract(x, channels, plot)
  }else{
    channels <- cyto_channels(x)
  }
  
  # Time parameter always uses data limits
  if("Time" %in% channels){
    limits <- "data"
  } 
  
  # Lower bound - use data limits
  rng <- suppressWarnings(
    range(x, type = "data")[, channels, drop = FALSE]
  )
  
  # Upper bound
  if(limits %in% c("auto", "data")){
    mx <- suppressWarnings(
      range(x, type = "data")[, channels, drop = FALSE][2,]
    )
  }else if(limits == "machine"){
    mx <- suppressWarnings(
      range(x, type = "instrument")[, channels, drop = FALSE][2,]
    )
  }
  rng[2,] <- mx
  
  # Replace lower data limit if > 0 - AUTO
  if(limits != "data" & anchor == TRUE){
    if(any(rng[1,] > 0)){
      rng[1, rng[1,] > 0] <- 0
    }
  }

  # ADD 2% BUFFER EITHER SIDE - 4% TOTAL IN PLOT
  if(buffer != 0){
    # BUFFER
    lapply(channels, function(z){
      # RANGE
      RNG <- rng[, z][2] - rng[, z][1]
      # ADD BUFFER
      rng[, z][1] <<- rng[, z][1] - buffer * RNG
      rng[, z][2] <<- rng[, z][2] + buffer * RNG
    })
  }
  
  return(rng)
}
 
#' @noRd
.cyto_range.flowSet <- function(x,
                                channels = NA,
                                limits = "machine",
                                plot = FALSE,
                                buffer = 0.03,
                                anchor = TRUE, ...){
  
  # Get ranges for each flowFrame as list
  rng <- fsApply(x, function(z){
    .cyto_range(z,
                channels = channels,
                limits = limits,
                plot = plot,
                buffer = buffer,
                anchor = anchor)
  }, simplify = TRUE)
  
  # Combine results
  rng <- do.call("rbind", rng)
  
  # Compute range for all samples for each channel
  rng <- suppressWarnings(lapply(channels, function(z){
    c(min(rng[,z][is.finite(rng[,z])]),
      max(rng[,z][is.finite(rng[,z])]))
  }))
  rng <- do.call("cbind",rng)
  rownames(rng) <- c("min","max")
  colnames(rng) <- channels
  
  return(rng)
  
}

#' @noRd
.cyto_range.GatingHierarchy <- function(x, 
                                        parent,
                                        channels = NA,
                                        limits = "machine",
                                        plot = FALSE,
                                        buffer = 0.03,
                                        anchor = TRUE, ...){
  
  # Extract data from GatingHierarchy
  x <- cyto_extract(x, parent)
  
  # Make call to flowFrame method
  rng <- .cyto_range(x,
                     channels = channels,
                     limits = limits,
                     plot = plot,
                     buffer = buffer,
                     anchor = anchor)
  
  return(rng)
  
}

#' @noRd
.cyto_range.GatingSet <- function(x,
                                  parent,
                                  channels = NA,
                                  limits = "machine",
                                  plot = FALSE,
                                  buffer = 0.03,
                                  anchor = TRUE, ...){
  
  # Extract data from GatingSet
  x <- cyto_extract(x, parent)
  
  # Make call to flowSet method
  rng <- .cyto_range(x,
                     channels = channels,
                     limits = limits,
                     plot = plot,
                     buffer = buffer,
                     anchor = anchor)
  
  return(rng)
  
}

#' @noRd
.cyto_range.list <- function(x,
                             parent,
                             channels = NA,
                             limits = "machine",
                             plot = FALSE,
                             buffer = 0.04,
                             anchor = TRUE){
  
  # Calculate range for each list element
  rng <- lapply(x, function(z){
    .cyto_range(z,
                parent = parent,
                channels = channels,
                limits = limits,
                plot = plot,
                buffer = buffer,
                anchor = anchor)
  })
  
  # Combine results
  rng <- do.call("rbind", rng)
  
  # Compute range for all samples for each channel
  rng <- suppressWarnings(lapply(channels, function(z){
    c(min(rng[,z][is.finite(rng[,z])]),
      max(rng[,z][is.finite(rng[,z])]))
  }))
  rng <- do.call("cbind",rng)
  rownames(rng) <- c("min","max")
  colnames(rng) <- channels
  
  return(rng)
  
}
