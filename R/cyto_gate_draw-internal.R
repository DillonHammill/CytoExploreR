## INTERNAL GATING FUNCTIONS ---------------------------------------------------

# Type is prepared within cyto_gate_draw so we can dispatch to necessary
# function by name.

# These functions DO NOT make calls to cyto_plot anymore! This is instead
# handled at the cyto_gate_draw level.

#' Create a gate by selecting co-ordinates on a plot
#'
#' @param type indicates the type of gate to create, by default set to
#'   \code{interval} and \code{polygon} for 1D and 2D plots respectively.
#'
#' @return rectangleGate, polygonGate, ellipsoidGate, quadrantGate or list of
#'   polygonGates.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @noRd
.cyto_gate_draw_dispatch <- function(args) {
  
  # DISPATCH - NA TYPE FOR NEGATED GATES
  if(!.all_na(args$type)) {
    return(
      do.call(
        paste0(".cyto_gate_", args$type, "_draw"), 
        args)
      )
  } else {
    return(NULL)
  }

}

#' Add labels for populations gated using cyto_gate_draw
#' 
#' @importFrom graphics par
#' 
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @noRd
.cyto_gate_draw_label <- function(x,
                                  alias = NA,
                                  gate = NA,
                                  channels = NULL,
                                  label = TRUE,
                                  label_text_size = 1,
                                  label_text_font = 2,
                                  label_text_col = "black",
                                  label_text_col_alpha = 1,
                                  label_fill = "white",
                                  label_fill_alpha = 0.6, 
                                  ...) {
  
  # X - MUST BE A CYTOSET
  
  # LABEL GATED POPULATIONS
  if(label == TRUE) {
    # PREPARE GATES
    if(!cyto_class(gate, "list", TRUE)) {
      gate <- list(gate)
    }
    # PULLDOWN ARGUMEMTS
    args <- .args_list(...)
    args <- args[!names(args) %in% c("x", "gate", "label", "channels")]
    # LIST OF POPULATIONS PER GATE
    pops <- cyto_gate_apply(x,
                            gate = gate)
    # REPEAT ARGUMENTS - NOT CHANNELS
    args <- structure(
      lapply(args, function(arg) {
        split(
          rep(arg, length.out = length(unlist(pops))),
          rep(1:length(pops), each = LAPPLY(pops, length))
        )
      }),
      names = names(args)
    )
    # UPDATE ARGUMENTS
    .args_update(args)
    # GATE CENTERS - BYPASS GATE CONVERSION?
    label_text_xy <- tryCatch(
      .cyto_gate_center(gate,
                        channels = channels), 
      error = function(e){
        return(rep(list(NA), length(gate)))
      }
    )
    # LABEL EACH POPULATION
    mapply(function(pop_list,
                    alias,
                    label_text_xy,
                    label_text_size,
                    label_text_font,
                    label_text_col,
                    label_text_col_alpha,
                    label_fill,
                    label_fill_alpha) {
      # LABEL EACH POPULATION IN POP_LIST
      lapply(seq_along(pop_list), function(z) {
        # POPULATION
        pop <- pop_list[[z]]
        # GATE CENTER
        if(!.all_na(label_text_xy)) {
          text_xy <- label_text_xy[z, ]
          # POPULATION CENTER
        } else {
          text_xy <- c()
          # X CO-ORDINATE
          if(nrow(pop[[1]]) < 2) {
            text_xy[1] <- mean(par("usr")[1:2])
          } else {
            text_xy[1] <- cyto_apply(
              pop,
              "cyto_stat_mode",
              channels = channels[1],
              input = "matrix",
              copy = FALSE,
              limits = matrix(
                par("usr")[1:2],
                ncol = 1,
                dimnames = list(
                  c("min", "max"),
                  channels[1]
                )
              )
            )[, 1]
          }
          # Y CO-ORDINATE
          if(length(channels) == 1) {
            text_xy[2] <- mean(par("usr")[3:4])
          } else {
            if(nrow(pop[[1]]) < 2) {
              text_xy[2] <- mean(par("usr")[3:4])
            } else {
              text_xy[2] <- cyto_apply(
                pop,
                "cyto_stat_mode",
                channels = channels[2],
                input = "matrix",
                copy = FALSE,
                limits = matrix(
                  par("usr")[3:4],
                  ncol = 1,
                  dimnames = list(
                    c("min", "max"),
                    channels[2]
                  )
                )
              )[, 1]
            }
          }
        }
        # PREPARE LABEL FOR GATED POPULATION
        stat <- paste(
          .round(
            sum(
              cyto_apply(
                pop,
                "cyto_stat_count",
                input = "matrix",
                copy = FALSE
              )
            ) /
            sum(
              cyto_apply(
                x,
                "cyto_stat_count",
                input = "matrix",
                copy = FALSE
              )
            ) * 100  
          ), 
          "%"  
        )
        # NaN
        stat <- gsub("NaN", "0.00", stat)
        # ALIAS
        if(!.all_na(alias[z])) {
          text <- paste(alias[z], stat, sep = "\n")
        } else {
          text <- stat
        }
        # PLOT LABELS
        cyto_plot_labeller(
          label_text = text,
          label_text_x = text_xy[1],
          label_text_y = text_xy[2],
          label_text_size = label_text_size,
          label_text_font = label_text_font,
          label_text_col = label_text_col,
          label_fill = label_fill,
          label_fill_alpha = label_fill_alpha
        )
      })
    },
    pops,
    alias,
    label_text_xy,
    label_text_size,
    label_text_font,
    label_text_col,
    label_text_col_alpha,
    label_fill,
    label_fill_alpha)
  }
  
}

#' @importFrom flowCore polygonGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator
#' @noRd
.cyto_gate_polygon_draw <- function(alias,
                                    channels,
                                    gate_point_shape = 16,
                                    gate_point_size = 1,
                                    gate_point_col = "red",
                                    gate_point_col_alpha = 1,
                                    gate_line_type = 1,
                                    gate_line_width = 2.5,
                                    gate_line_col = "red",
                                    gate_line_col_alpha = 1, 
                                    ...) {
  
  # HIDE ERROR MESSAGE ON CLOSE
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # INSTRUCTIONS
  message(
    paste(
      "Select at least 3 points to construct a polygon gate around the",
      alias, "population. \n"
    )
  )

  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # COORDS - POINT BY POINT - SET 100 POINT LIMIT
  coords <- list()
  for(z in seq_len(100)){
    # SELECT POINT
    pt <- locator(n = 1,
                  type = "p",
                  pch = gate_point_shape,
                  cex = gate_point_size,
                  col = gate_point_col)
    # POINT SELECTED ADD TO COORDS
    if(!is.null(pt)){
      pt <- do.call("cbind", pt)
    }
    # ADD POINT TO COORDS
    coords[[z]] <- pt
    # JOIN TO PREVIOUS POINT
    if(!is.null(pt) & length(coords) > 1){
      # ADD LINE BETWEEN POINTS
      lines(x = c(coords[[z-1]][, "x"], coords[[z]][, "x"]),
            y = c(coords[[z-1]][, "y"], coords[[z]][, "y"]),
            lty = gate_line_type,
            lwd = gate_line_width,
            col = gate_line_col)
      # CLOSE GATE
    }else if(is.null(pt)){
      if(length(coords) < 3){
        stop("A minimum of 3 points is required to construct a polygon gate.")
      }else{
        # JOIN FIRST & LAST POINTS
        lines(x = c(coords[[1]][, "x"], coords[[z-1]][, "x"]),
              y = c(coords[[1]][, "y"], coords[[z-1]][, "y"]),
              lty = gate_line_type,
              lwd = gate_line_width,
              col = gate_line_col)
        # TERMINATE LOOP
        break()
      }
    }
  }
  coords <- do.call("rbind", coords)
  
  # CHANNELS
  colnames(coords) <- channels
  
  # RETURN GATE OBJECT
  return(
    polygonGate(.gate = coords, filterId = alias)
  )

}

#' @importFrom flowCore rectangleGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator
#' @noRd
.cyto_gate_rectangle_draw <- function(alias,
                                      channels,
                                      gate_point_shape = 16,
                                      gate_point_size = 1,
                                      gate_point_col = "red",
                                      gate_point_col_alpha = 1,
                                      gate_line_type = 1,
                                      gate_line_width = 2.5,
                                      gate_line_col = "red",
                                      gate_line_col_alpha = 1, 
                                      ...) {

  # HIDE ERROR MESSAGE ON CLOSE
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # INSTRUCTIONS
  message(
    paste(
      "Select 2 diagonal points to construct a rectangle gate around the",
      alias, "population. \n"
    )
  )
  
  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # COORDS - POINT BY POINT - SET 2 POINT LIMIT
  coords <- list()
  for(z in seq_len(2)){
    # SELECT POINT
    pt <- locator(n = 1,
                  type = "p",
                  pch = gate_point_shape,
                  cex = gate_point_size,
                  col = gate_point_col)
    # POINT SELECTED ADD TO COORDS
    if(!is.null(pt)){
      pt <- do.call("cbind", pt)
    }
    # ADD POINT TO COORDS
    coords[[z]] <- pt
    # TERMINATE LOOP
    if(is.null(pt)){
      if(length(coords) < 2){
        stop("A minimum of 2 points is required to construct a rectangle gate.")
      }else{
        break()
      }
    }
  }
  coords <- do.call("rbind", coords)
  
  # CHANNELS
  colnames(coords) <- channels
  
  # CONSTRUCT GATE
  gate <- rectangleGate(.gate = coords, filterId = alias)
  
  # PLOT GATE
  cyto_plot_gate(gate,
                 channels = channels,
                 gate_line_type = gate_line_type,
                 gate_line_width = gate_line_width,
                 gate_line_col = gate_line_col
  )
  
  # RETURN GATE OBJECT
  return(
    gate
  )

}

#' @importFrom flowCore rectangleGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator par
#' @noRd
.cyto_gate_interval_draw <- function(alias,
                                     channels,
                                     axis = "x",
                                     gate_point_shape = 16,
                                     gate_point_size = 1,
                                     gate_point_col = "red",
                                     gate_point_col_alpha = 1,
                                     gate_line_type = 1,
                                     gate_line_width = 2.5,
                                     gate_line_col = "red",
                                     gate_line_col_alpha = 1, 
                                     ...) {
  
  # HIDE ERROR MESSAGE ON CLOSE
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # PLOT LIMITS
  par_usr <- par("usr")
  
  # X LIMITS + BUFFER
  par_xlim <- par_usr[1:2]
  par_xrng <- par_xlim[2] - par_xlim[1]
  par_xpad <- (par_xrng - par_xrng / 1.04) / 2
  par_xmin <- par_xlim[1] + 0.5 * par_xpad
  par_xmax <- par_xlim[2] - 0.5 * par_xpad
  
  # Y LIMITS + BUFFER
  par_ylim <- par_usr[3:4]
  par_yrng <- par_ylim[2] - par_ylim[1]
  par_ypad <- (par_yrng - par_yrng / 1.04) / 2
  par_ymin <- par_ylim[1] + 0.5 * par_ypad
  par_ymax <- par_ylim[2] - 0.5 * par_ypad
  
  # INSTRUCTIONS
  message(
    paste(
      "Select the lower and upper bounds of the",
      alias, "population to construct an interval gate. \n"
    )
  )
  
  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # COORDS - POINT BY POINT - SET 2 POINT LIMIT
  coords <- list()
  for(z in seq_len(2)){
    # SELECT POINT
    pt <- locator(n = 1,
                  type = "p",
                  pch = gate_point_shape,
                  cex = gate_point_size,
                  col = gate_point_col)
    # POINT SELECTED ADD TO COORDS
    if(!is.null(pt)){
      pt <- do.call("cbind", pt)
    }
    # ADD POINT TO COORDS
    coords[[z]] <- pt
    # ADD LINES - X COORS + PAR_YMIN & PAR_YMAX
    if(length(channels) == 1 | axis == "x"){
      lines(x = c(coords[[z]][, "x"], coords[[z]][, "x"]),
            y = c(par_ymin, par_ymax),
            lty = gate_line_type,
            lwd = gate_line_width,
            col = gate_line_col)
    }else{
      lines(x = c(par_xmin, par_xmax),
            y = c(coords[[z]][, "y"], coords[[z]][, "y"]),
            lty = gate_line_type,
            lwd = gate_line_width,
            col = gate_line_col)
    }
    # CLOSE GATE
    if(z == 2){
      if(length(channels) == 1 | axis == "x"){
        # BOTTOM LINE
        lines(x = c(coords[[z-1]][, "x"], coords[[z]][, "x"]),
              y = c(par_ymin, par_ymin),
              lty = gate_line_type,
              lwd = gate_line_width,
              col = gate_line_col)
        # TOP LINE
        lines(x = c(coords[[z-1]][, "x"], coords[[z]][, "x"]),
              y = c(par_ymax, par_ymax),
              lty = gate_line_type,
              lwd = gate_line_width,
              col = gate_line_col)
      }else{
        # LOWER LINE
        lines(x = c(par_xmin, par_xmin),
              y = c(coords[[z-1]][, "y"], coords[[z]][, "y"]),
              lty = gate_line_type,
              lwd = gate_line_width,
              col = gate_line_col)
        # UPPER LINE
        lines(x = c(par_xmax, par_xmax),
              y = c(coords[[z-1]][, "y"], coords[[z]][, "y"]),
              lty = gate_line_type,
              lwd = gate_line_width,
              col = gate_line_col)
      }
    }
  }
  coords <- do.call("rbind", coords)
  
  # PREPARE COORDS
  if (axis == "x") {
    if (length(channels) == 1) {
      coords <- data.frame(x = coords[, 1])
      coords <- as.matrix(coords)
      colnames(coords) <- channels[1]
      rownames(coords) <- c("min", "max")
    } else if (length(channels) == 2) {
      coords <- data.frame(x = coords[, 1], y = c(-Inf, Inf))
      coords <- as.matrix(coords)
      colnames(coords) <- channels
      rownames(coords) <- c("min", "max")
    }
    gate <- rectangleGate(.gate = coords, filterId = alias)
  } else if (axis == "y") {
    if (length(channels) == 1) {
      stop("Cannot gate y axis if a single channel is supplied.")
    }
    coords <- data.frame(x = c(-Inf, Inf), y = coords[, 2])
    coords <- as.matrix(coords)
    colnames(coords) <- channels
    rownames(coords) <- c("min", "max")
    
    gate <- rectangleGate(.gate = coords, filterId = alias)
  }
  
  # RETURN GATE OBJECT
  return(gate)
  
}

#' @importFrom flowCore rectangleGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator
#' @noRd
.cyto_gate_threshold_draw <- function(alias,
                                      channels,
                                      gate_point_shape = 16,
                                      gate_point_size = 1,
                                      gate_point_col = "red",
                                      gate_point_col_alpha = 1,
                                      gate_line_type = 1,
                                      gate_line_width = 2.5,
                                      gate_line_col = "red",
                                      gate_line_col_alpha = 1, 
                                      ...) {
  
  # HIDE ERROR MESSAGES
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # INSTRUCTIONS
  message(
    paste(
      "Select the lower bound of the",
      alias, "population to construct a threshold gate. \n"
    )
  )

  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # COORDS - 1 POINT ONLY
  coords <- locator(
    n = 1,
    type = "p",
    pch = gate_point_shape,
    cex = gate_point_size,
    col = gate_point_col
  )
  
  # TIDY GATE COORDS
  if (length(channels) == 1) {
    pts <- data.frame(x = c(coords$x, Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    rownames(pts) <- c("min", "max")
  } else if (length(channels) == 2) {
    pts <- data.frame(x = c(coords$x, Inf), 
                      y = c(coords$y, Inf))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rownames(pts) <- c("min", "max")
  }
  
  # CONSTRUCT GATE
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  # PLOT GATE
  cyto_plot_gate(
    gate = gate,
    channels = channels,
    gate_line_type = gate_line_type,
    gate_line_width = gate_line_width,
    gate_line_col = gate_line_col
  )
  
  # RETURN GATE OBJECT
  return(gate)
  
}

#' @importFrom flowCore rectangleGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator
#' @noRd
.cyto_gate_boundary_draw <- function(alias,
                                     channels,
                                     gate_point_shape = 16,
                                     gate_point_size = 1,
                                     gate_point_col = "red",
                                     gate_point_col_alpha = 1,
                                     gate_line_type = 1,
                                     gate_line_width = 2.5,
                                     gate_line_col = "red",
                                     gate_line_col_alpha = 1, 
                                     ...) {
  
  # HIDE ERROR MESSAGES
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # INSTRUCTIONS
  message(
    paste(
      "Select the upper bound of the",
      alias, "population to construct a boundary gate. \n"
    )
  )
  
  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # GATE COORD 1 POINT
  coords <- locator(
    n = 1,
    type = "p",
    pch = gate_point_shape,
    cex = gate_point_size,
    col = gate_point_col
  )
  
  # TIDY GATE COORDS
  if (length(channels) == 1) {
    pts <- data.frame(x = c(-Inf, coords$x))
    pts <- as.matrix(pts)
    colnames(pts) <- channels[1]
    rownames(pts) <- c("min", "max")
  } else if (length(channels) == 2) {
    pts <- data.frame(x = c(-Inf, coords$x), y = c(-Inf, coords$y))
    pts <- as.matrix(pts)
    colnames(pts) <- channels
    rownames(pts) <- c("min", "max")
  }
  
  # CONSTRUCT GATE
  gate <- rectangleGate(.gate = pts, filterId = alias)
  
  # PLOT GATE
  cyto_plot_gate(
    gate = gate,
    channels = channels,
    gate_line_type = gate_line_type,
    gate_line_width = gate_line_width,
    gate_line_col = gate_line_col
  )
  
  # RETURN GATE OBJECT
  return(gate)
  
}

#' @importFrom flowCore ellipsoidGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator
#' @noRd
.cyto_gate_ellipse_draw <- function(alias,
                                    channels,
                                    gate_point_shape = 16,
                                    gate_point_size = 1,
                                    gate_point_col = "red",
                                    gate_point_col_alpha = 1,
                                    gate_line_type = 1,
                                    gate_line_width = 2.5,
                                    gate_line_col = "red",
                                    gate_line_col_alpha = 1, 
                                    ...) {

  # HIDE ERROR MESSAGES
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # INSTRUCTIONS
  message(
    paste(
      "Select 4 points to define the limits of the",
      alias, "population to construct an ellipsoid gate. \n"
    )
  )
  
  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # COORDS - POINT BY POINT - SET 4 POINT LIMIT
  coords <- list()
  for(z in seq_len(4)){
    # SELECT POINT
    pt <- locator(n = 1,
                  type = "p",
                  pch = gate_point_shape,
                  cex = gate_point_size,
                  col = gate_point_col)
    # POINT SELECTED ADD TO COORDS
    if(!is.null(pt)){
      pt <- do.call("cbind", pt)
    }
    # ADD POINT TO COORDS
    coords[[z]] <- pt
  }
  coords <- do.call("rbind", coords)
  coords <- data.frame(coords)
  colnames(coords) <- channels
  
  # COMPUTE ANTIPODAL CO-ORDINATES ---------------------------------------------
  
  # MAJOR AXIS
  dst <- as.matrix(stats::dist(coords))
  mj.pts <- coords[which(dst == max(dst), arr.ind = TRUE)[1, ], ]
  
  # MINOR AXIS
  mr.pts <- coords[!coords[, 1] %in% mj.pts[, 1] & 
                     !coords[, 2] %in% mj.pts[, 2], ]
  
  # CENTER - MAJOR AXIS
  mj.center <- c(
    (sum(mj.pts[, 1]) / nrow(mj.pts)),
    (sum(mj.pts[, 2]) / nrow(mj.pts))
  )
  
  # CENTROID
  center <- colMeans(coords)
  
  # ADJUST MAJOR AXIS TO CENTROID
  adj <- c((mj.center[1] - center[1]), (mj.center[2] - center[2]))
  mj.pts[, 1] <- mj.pts[, 1] - adj[1]
  mj.pts[, 2] <- mj.pts[, 2] - adj[2]
  
  # UPDATE CO-ORDINATES
  coords <- rbind(mj.pts, mr.pts)
  
  # CONSTRUCT ELLIPSOIDGATE ----------------------------------------------------
  
  # COMPUTE COVARIANCE MATRIX - CANNOT CALL COV FLOWCORE ERROR
  cvm <- .cyto_ellipse_cov(coords)
  
  # CONSTRUCT GATE
  gate <- ellipsoidGate(
    .gate = cvm,
    mean = center,
    filterId = alias,
    distance = 1
  )

  # PLOT GATE
  cyto_plot_gate(
    gate,
    channels = channels,
    gate_line_type = gate_line_type,
    gate_line_width = gate_line_width,
    gate_line_col = gate_line_col
  )

  # RETURN GATE OBJECT
  return(gate)

}

#' @noRd
.cyto_ellipse_cov <- function(x) {

  # CHANNELS
  channels <- colnames(x)
  
  # POINTS ON MAJOR/MINOR AXES
  dst <- as.matrix(stats::dist(x))
  mj.pts <- x[which(dst == max(dst), arr.ind = TRUE)[1, ], ]
  mr.pts <- x[!x[, 1] %in% mj.pts[, 1] &
              !x[, 2] %in% mj.pts[, 2], ]
  
  # CENTER
  center <- colMeans(x)
  
  # MAJOR POINT ABOVE CENTER
  max.pt <- mj.pts[mj.pts[, 2] > center[2], ]
  
  # RADIUS - MAJOR AXIS
  a <- stats::dist(mj.pts) / 2
  
  # RADIUS MINR AXIS
  b <- stats::dist(mr.pts) / 2
  
  # ANGLE BETWEEN CENTER & MAX.PT
  if (max.pt[1] > center[1]) { # angle < pi/2
    mj.pt.ct <- cbind(max.pt[1], center[2])
    colnames(mj.pt.ct) <- channels
    adj <- stats::dist(rbind(center, mj.pt.ct))
    angle <- acos(adj / a) # [-1, 1]
  } else if (max.pt[1] <= center[1]) { # angle >= pi/2
    mj.pt.ct <- cbind(center[1], max.pt[2])
    colnames(mj.pt.ct) <- channels
    opp <- stats::dist(as.matrix(rbind(max.pt, mj.pt.ct)))
    angle <- pi / 2 + asin(opp / a)
  }
  
  # COVARIANCE MATRIX
  cinv <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
  cinv[1, 1] <- (((cos(angle) * cos(angle)) /
                    (a^2)) + ((sin(angle) * sin(angle)) / (b^2)))
  cinv[2, 1] <- sin(angle) * cos(angle) * ((1 / (a^2)) - (1 / (b^2)))
  cinv[1, 2] <- cinv[2, 1]
  cinv[2, 2] <- (((sin(angle) * sin(angle)) / (a^2)) +
                   ((cos(angle) * cos(angle)) / (b^2)))
  cvm <- solve(cinv)
  dimnames(cvm) <- list(channels, channels)
  
  return(cvm)
  
  # CYTOLIB APPROACH:
  #
  # # CHANNELS
  # channels <- colnames(x)
  # 
  # # CENTER
  # center <- colMeans(x)
  # 
  # # INDICES OF FAR RIGHT & FAR LEFT
  # ind <- c(
  #   which.max(x[, 1]),
  #   which.min(x[, 1])
  # )
  # 
  # # FAR RIGHT
  # R <- as.numeric(x[ind[1], ])
  # 
  # # FAR LEFT
  # L <- as.numeric(x[ind[2], ])
  # 
  # # HALF AXIS LENGTH
  # a <- sqrt((L[1] - R[1])^2 + (L[2] - R[2])^2)/2
  # a2 <- a*a
  # 
  # # OTHER AXIS
  # ind <- seq_len(nrow(x))[-ind]
  # V1 <- as.numeric(x[ind[1], ])
  # V2 <- as.numeric(x[ind[2], ])
  # 
  # # HALF AXIS LENGTH
  # b <- sqrt(
  #   (V1[1]-V2[1])^2 + 
  #     (V1[2] - V2[2])^2
  #   )/2
  # b2 <- b*b
  # 
  # # NORMALISE
  # L_norm <- sqrt(L[1]*L[1] + L[2]*L[2])
  # x1 <- L[1]/L_norm
  # y1 <- L[2]/L_norm
  # 
  # V1_norm <- sqrt(V1[1]*V1[1] + V1[2]*V1[2])
  # x2 <- V1[1]/V1_norm
  # y2 <- V1[2]/V1_norm
  # 
  # # COVARIANCE
  # p1 <-  c(x1 * x1 * a2 + x2 * x2 * b2,
  #          x1 * y1 * a2 + x2 * y2 * b2)
  # 
  # p2 <- c(p1[2],
  #         y1 * y1 * a2 + y2 * y2 * b2)
  # 
  # # COVARIANCE MATRIX
  # cov <- matrix(
  #   c(p2, rev(p1)),
  #   ncol = 2,
  #   byrow = FALSE,
  #   dimnames = list(channels, channels)
  # )
  # 
  # return(cov)
  
}

#' @importFrom flowCore quadGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator
#' @noRd
.cyto_gate_quadrant_draw <- function(alias,
                                     channels,
                                     gate_point_shape = 16,
                                     gate_point_size = 1,
                                     gate_point_col = "red",
                                     gate_point_col_alpha = 1,
                                     gate_line_type = 1,
                                     gate_line_width = 2.5,
                                     gate_line_col = "red",
                                     gate_line_col_alpha = 1, 
                                     ...) {
  
  # HIDE ERROR MESSAGES
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE)) 
  
  # ALIAS
  if (is.null(alias)) {
    alias <- c(
      paste0(channels[1], "+", channels[2], "+"),
      paste0(channels[1], "-", channels[2], "+"),
      paste0(channels[1], "+", channels[2], "-"),
      paste0(channels[1], "-", channels[2], "-")
    )
  }
  
  # ALIAS CHECK
  if (!length(alias) == 4) {
    stop("'alias' must contain 4 population names for quadrant gates.")
  }
  
  # INSTRUCTIONS
  message(
    paste("Select the center point to construct quadrant gates. \n")
  )

  # COLOUR - TRANSPARENCY
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha)
  
  # GATE COORDS
  pts <- locator(
    n = 1,
    type = "p",
    pch = gate_point_shape,
    cex = gate_point_size,
    col = gate_point_col
  )
  
  # CO-ORDINATES MATRIX
  pts <- matrix(unlist(pts), ncol = 2)
  
  # CHANNELS
  colnames(pts) <- channels
  
  # QUADGATE CONSTRUCTION
  gate <- quadGate(.gate = pts, 
                   filterId = paste(alias, collapse = "|"))
  
  # PLOT GATE
  cyto_plot_gate(gate, 
                 channels = channels,
                 gate_line_type = gate_line_type,
                 gate_line_width = gate_line_width,
                 gate_line_col = gate_line_col)
  
  # RETURN GATE OBJECT
  return(gate)
  
}

#' @importFrom flowCore polygonGate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics locator par
#' @noRd
.cyto_gate_web_draw <- function(alias,
                                channels,
                                gate_point_shape = 16,
                                gate_point_size = 1,
                                gate_point_col = "red",
                                gate_point_col_alpha = 1,
                                gate_line_type = 1,
                                gate_line_width = 2.5,
                                gate_line_col = "red",
                                gate_line_col_alpha = 1,
                                ...) {     

  # WARNING
  message("Web gates are an experimental feature - use at your own risk!")
  
  # ALIAS
  if(length(alias) == 1) {
    stop(
      "Multiple populations must be supplied to 'alias' to construct web gates!"
    )
  }
  
  # HIDE ERROR MESSAGES
  options("show.error.messages" = FALSE)
  on.exit(options("show.error.messages" = TRUE))
  
  # INSTRUCTIONS
  message("Select the center of the web gate.")
  
  # POINT & LINE COLOURS
  gate_point_col <- adjustcolor(gate_point_col, gate_point_col_alpha)
  gate_line_col <- adjustcolor(gate_line_col, gate_line_col_alpha) 
  
  # GATE CENTER
  center <- locator(
    n = 1,
    type = "p",
    pch = gate_point_shape,
    cex = gate_point_size,
    col = gate_point_col
  )
  
  # User Prompt
  message("Select surrounding co-ordinates on plot edges to draw a web gate.")
  
  # Minimum and maximum limits of plot
  xmin <- round(par("usr")[1], 2)
  xmax <- round(par("usr")[2], 2)
  ymin <- round(par("usr")[3], 2)
  ymax <- round(par("usr")[4], 2)
  
  # Get all gate co-ordinates - c(center, others)
  coords <- lapply(seq_len(length(alias)), function(x) {
    
    # GATE COORDS
    pt <- locator(
      n = 1,
      type = "p",
      pch = gate_point_shape,
      cex = gate_point_size,
      col = gate_point_col
    )
    
    # LINE TO CENTER
    lines(
      x = c(center$x, pt$x),
      y = c(center$y, pt$y),
      lty = gate_line_type,
      lwd = gate_line_width,
      col = gate_line_col
    )
    
    return(c(pt$x, pt$y))
  })
  coords <- as.data.frame(do.call(rbind, coords))
  colnames(coords) <- c("x", "y")
  coords <- rbind(center, coords)
  
  # Determine which quadrants the points are in
  # bottom left anti-clockwise to top left (relative to center)
  quads <- c(0, rep(NA, length(alias)))
  for (i in seq_len(length(coords$x))[-1]) {
    
    # Bottom left Q1
    if (coords[i, ]$x < center$x & coords[i, ]$y <= center$y) {
      quads[i] <- 1
      
      # Bottom right Q2
    } else if (coords[i, ]$x >= center$x & coords[i, ]$y < center$y) {
      quads[i] <- 2
      
      # Top right Q3
    } else if (coords[i, ]$x > center$x & coords[i, ]$y >= center$y) {
      quads[i] <- 3
      
      # Top left Q4
    } else if (coords[i, ]$x <= center$x & coords[i, ]$y > center$y) {
      quads[i] <- 4
    }
  }
  coords[, "Q"] <- quads
  coords <- coords[with(coords, order(coords$Q)), ]
  
  # Push points to plot limits (intersection with plot limits)
  
  # Quadrant 1: find limit intercept and modify point co-ordinates
  if (1 %in% coords$Q) {
    q1 <- coords[coords$Q == 1, ]
    for (x in seq_len(length(q1$Q))) {
      
      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q1[x, "x"], q1[x, "y"]),
        c(xmin, center$y),
        c(xmin, ymin)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q1[x, "x"], q1[x, "y"]),
        c(center$x, ymin),
        c(xmin, ymin)
      )
      
      # Check which axis the point should be pushed onto
      if (vint[2] >= ymin) {
        q1[x, c("x", "y")] <- vint
      } else if (vint[2] < ymin) {
        q1[x, c("x", "y")] <- hint
      }
    }
    coords[coords$Q == 1, ] <- q1
  }
  
  # Quadrant 2: find limit intercept and modify point co-ordinates
  if (2 %in% coords$Q) {
    q2 <- coords[coords$Q == 2, ]
    for (x in seq_len(length(q2$Q))) {
      
      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q2[x, "x"], q2[x, "y"]),
        c(xmax, center$y),
        c(xmax, ymin)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q2[x, "x"], q2[x, "y"]),
        c(center$x, ymin),
        c(xmax, ymin)
      )
      
      # Check which axis the point should be pushed onto
      if (vint[2] >= ymin) {
        q2[x, c("x", "y")] <- vint
      } else if (vint[2] < ymin) {
        q2[x, c("x", "y")] <- hint
      }
    }
    coords[coords$Q == 2, ] <- q2
  }
  
  # Quadrant 3: find limit intercept and modify point co-ordinates
  if (3 %in% coords$Q) {
    q3 <- coords[coords$Q == 3, ]
    for (x in seq_len(length(q3$Q))) {
      
      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q3[x, "x"], q3[x, "y"]),
        c(xmax, ymax),
        c(xmax, center$y)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q3[x, "x"], q3[x, "y"]),
        c(center$x, ymax),
        c(xmax, ymax)
      )
      
      # Check which axis the point should be pushed onto
      if (vint[2] >= ymax) {
        q3[x, c("x", "y")] <- hint
      } else if (vint[2] < ymax) {
        q3[x, c("x", "y")] <- vint
      }
    }
    coords[coords$Q == 3, ] <- q3
  }
  
  # Quadrant 4: find limit intercept and modify point co-ordinates
  if (4 %in% coords$Q) {
    q4 <- coords[coords$Q == 4, ]
    for (x in seq_len(length(q4$Q))) {
      
      # Calculate intersection with horizontal and vertical axes
      vint <- linesIntercept(
        c(center$x, center$y),
        c(q4[x, "x"], q4[x, "y"]),
        c(xmin, ymax),
        c(xmin, center$y)
      )
      hint <- linesIntercept(
        c(center$x, center$y),
        c(q4[x, "x"], q4[x, "y"]),
        c(xmin, ymax),
        c(center$x, ymax)
      )
      
      # Check which axis the point should be pushed onto
      if (vint[2] >= ymax) {
        q4[x, c("x", "y")] <- hint
      } else if (vint[2] < ymax) {
        q4[x, c("x", "y")] <- vint
      }
    }
    coords[coords$Q == 4, ] <- q4
  }
  
  # If multiple points in same quadrant order anticlockwise Q1-Q4
  if (anyDuplicated(coords$Q) != 0) {
    
    # Quadrant 1
    if (1 %in% coords$Q[duplicated(coords$Q)]) {
      
      # Multiple points in Q1 - sort by -y then +x
      q1 <- coords[coords$Q == 1, ]
      q1 <- q1[with(q1, order(-q1$y, q1$x)), ]
      coords[coords$Q == 1, c("x", "y")] <- q1[, c("x", "y")]
    }
    
    # Quadrant 2
    if (2 %in% coords$Q[duplicated(coords$Q)]) {
      
      # Multiple points in Q2 - sort by +x then +y
      q2 <- coords[coords$Q == 2, ]
      q2 <- q2[with(q2, order(q2$x, q2$y)), ]
      coords[coords$Q == 2, c("x", "y")] <- q2[, c("x", "y")]
    }
    
    # Quadrant 3
    if (3 %in% coords$Q[duplicated(coords$Q)]) {
      
      # Multiple points in Q3 - sort by +y then -x
      q3 <- coords[coords$Q == 3, ]
      q3 <- q3[with(q3, order(q3$y, -q3$x)), ]
      coords[coords$Q == 3, c("x", "y")] <- q3[, c("x", "y")]
    }
    
    # Quadrant 4
    if (4 %in% coords$Q[duplicated(coords$Q)]) {
      
      # Multiple points in Q4 - sort by -x then -y
      q4 <- coords[coords$Q == 4, ]
      q4 <- q4[with(q4, order(-q4$x, -q4$y)), ]
      coords[coords$Q == 4, c("x", "y")] <- q4[, c("x", "y")]
    }
  }
  
  # Construct gates using input points
  # Duplicate first point after last point
  coords[(length(coords$Q) + 1), ] <- coords[2, ]
  coords[] <- lapply(coords, round, 4)
  
  # Gate coordinates using input points
  gates <- list()
  for (i in 2:(length(coords$Q) - 1)) {
    gates[[i - 1]] <- rbind(coords[1, ], coords[i, ], coords[i + 1, ])
  }
  
  # Check if a corner lies between the points - add as gate co-ordinate
  # Calculate corner points using min & max values
  Q1 <- c(xmin, ymin, 1)
  Q2 <- c(xmax, ymin, 2)
  Q3 <- c(xmax, ymax, 3)
  Q4 <- c(xmin, ymax, 4)
  Q <- matrix(c(Q1, Q2, Q3, Q4), byrow = TRUE, nrow = 4)
  colnames(Q) <- c("x", "y", "Q")
  Q <- data.frame(Q)
  
  # LAST GATE INHERITS REMAINING CORNERS
  indx <- seq_len(length(alias) - 1)
  
  # Add corners to appropriate gates step-wise
  gates[indx] <- lapply(gates[indx], function(x) {
    
    # DUPLICATION - points in same quadrant
    if (any(duplicated(x$Q))) {
      
      # Quadrant 1
      if (1 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "x"] == xmin & x[3, "x"] != xmin) {
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, x[3, ])
          
          # Remove Q1 from Q
          if (1 %in% Q[, "Q"]) {
            Q <<- Q[-match(1, Q[, "Q"]), ]
          }
        }
      }
      
      # Quadrant 2
      if (2 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "y"] == ymin & x[3, "y"] != ymin) {
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])
          
          # Remove Q2 from Q
          if (2 %in% Q[, "Q"]) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        }
      }
      
      # Quadrant 3
      if (3 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "x"] == xmax & x[3, "x"] != xmax) {
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])
          
          # Remove Q3 from Q
          if (3 %in% Q[, "Q"]) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        }
      }
      
      # Quadrant 4
      if (4 %in% x$Q[duplicated(x$Q)]) {
        if (x[2, "y"] == ymax & x[3, "y"] != ymax) {
          
          # Include Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q4, x[3, ])
          
          # Remove Q4 from Q
          if (4 %in% Q[, "Q"]) {
            Q <<- Q[-match(4, Q[, "Q"]), ]
          }
        }
      }
      
      # ADJACENT - points in adjacent quadrants
    } else if ((x[3, "Q"] - x[2, "Q"]) %in% c(0, 1)) {
      
      # Q1-Q2
      if (x[2, "Q"] == 1 & x[3, "Q"] == 2) {
        if (x[2, "x"] == xmin & x[3, "x"] == xmax) {
          
          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, Q2, x[3, ])
          
          # Remove Q1 and Q2 from Q
          if (any(c(1, 2) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(1, 2), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] == xmin & x[3, "x"] != xmax) {
          
          # Include Q1 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, x[3, ])
          
          # Remove Q1 from Q
          if (any(1 %in% Q[, "Q"])) {
            Q <<- Q[-match(1, Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmin & x[3, "x"] == xmax) {
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])
          
          # Remove Q2 from Q
          if (any(2 %in% Q[, "Q"])) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        }
        
        # Q2-Q3
      } else if (x[2, "Q"] == 2 & x[3, "Q"] == 3) {
        if (x[2, "y"] == ymin & x[3, "y"] == ymax) {
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, x[3, ])
          
          # Remove Q2 and Q3 from Q
          if (any(c(2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] == ymin & x[3, "y"] != ymax) {
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])
          
          # Remove Q2 from Q
          if (any(2 %in% Q[, "Q"])) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        } else if (x[2, "y"] != ymin & x[3, "y"] == ymax) {
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])
          
          # Remove Q3 from Q
          if (any(3 %in% Q[, "Q"])) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        }
        
        # Q3-Q4
      } else if (x[2, "Q"] == 3 & x[3, "Q"] == 4) {
        if (x[2, "x"] == xmax & x[3, "x"] == xmin) {
          
          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, Q4, x[3, ])
          
          # Remove Q3 and Q4 from Q
          if (any(c(3, 4) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(3, 4), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] == xmax & x[3, "x"] != xmin) {
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])
          
          # Remove Q3 from Q
          if (any(3 %in% Q[, "Q"])) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmax & x[3, "x"] == xmin) {
          
          # Include Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q4, x[3, ])
          
          # Remove Q4 from Q
          if (any(4 %in% Q[, "Q"])) {
            Q <<- Q[-match(4, Q[, "Q"]), ]
          }
        }
      }
      
      # SEPARATED - points separated by a quadrant
    } else if (x[3, "Q"] - x[2, "Q"] == 2) {
      
      # Q1-Q3
      if (x[2, "Q"] == 1 & x[3, "Q"] == 3) {
        if (x[2, "x"] == xmin & x[3, "y"] == ymax) {
          
          # Include Q1, Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, Q2, Q3, x[3, ])
          
          # Remove Q1, Q2 and Q3 from Q
          if (any(c(1, 2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(1, 2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] == xmin & x[3, "y"] != ymax) {
          
          # Include Q1 & Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q1, Q2, x[3, ])
          
          # Remove Q1 and Q2 from Q
          if (any(c(1, 2) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(1, 2), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmin & x[3, "y"] == ymax) {
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, x[3, ])
          
          # Remove Q2 and Q3 from Q
          if (any(c(2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "x"] != xmin & x[3, "y"] != ymax) {
          
          # Include Q2 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, x[3, ])
          
          # Remove Q2 from Q
          if (any(2 %in% Q[, "Q"])) {
            Q <<- Q[-match(2, Q[, "Q"]), ]
          }
        }
        
        # Q2-Q4
      } else if (x[2, "Q"] == 2 & x[3, "Q"] == 4) {
        if (x[2, "y"] == ymin & x[3, "x"] == xmin) {
          
          # Include Q2, Q3 & Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, Q4, x[3, ])
          
          # Remove Q2, Q3 and Q4 from Q
          if (any(c(2, 3, 4) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3, 4), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] == ymin & x[3, "x"] != xmin) {
          
          # Include Q2 & Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q2, Q3, x[3, ])
          
          # Remove Q2 and Q3 from Q
          if (any(c(2, 3) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(2, 3), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] != ymin & x[3, "x"] == xmin) {
          
          # Include Q3 & Q4 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, Q4, x[3, ])
          
          # Remove Q3 and Q4
          if (any(c(3, 4) %in% Q[, "Q"])) {
            Q <<- Q[-match(c(3, 4), Q[, "Q"]), ]
          }
        } else if (x[2, "y"] != ymin & x[3, "x"] != xmin) {
          
          # Include Q3 corner in gate
          x <- rbind(x[c(1, 2), ], Q3, x[3, ])
          
          # Remove Q3 from Q
          if (any(3 %in% Q[, "Q"])) {
            Q <<- Q[-match(3, Q[, "Q"]), ]
          }
        }
      }
    }
    
    return(x)
  })
  
  # Last gate inherits remaining corners
  if (nrow(Q) != 0) {
    if (length(which(Q[, "Q"] >= gates[[length(alias)]][2, "Q"])) != 0) {
      g <- Q[which(Q[, "Q"] >= gates[[length(alias)]][2, "Q"]), ]
    }else{
      g <- NULL
    }
    
    if (length(which(Q[, "Q"] < gates[[length(alias)]][2, "Q"])) != 0) {
      r <- Q[which(Q[, "Q"] < gates[[length(alias)]][2, "Q"]), ]
    }else{
      r <- NULL
    }
    
    if (!is.null(g) & !is.null(r)) {
      Q <- rbind(g, r)
    } else if (!is.null(g) & is.null(r)) {
      Q <- g
    } else if (is.null(g) & !is.null(r)) {
      Q <- r
    }
    
    gates[[length(alias)]] <- rbind(
      gates[[length(alias)]][c(1, 2), ],
      Q,
      gates[[length(alias)]][3, ]
    )
  }
  
  # CONSTRUCT GATES
  gates <- lapply(seq(1, length(gates), 1), function(x) {
    coords <- as.matrix(gates[[x]])[, -3]
    colnames(coords) <- channels
    rownames(coords) <- NULL
    
    # CONSTRUCT GATE
    gate <- polygonGate(.gate = coords, filterId = alias[x])
    
    # PLOT GATE
    cyto_plot_gate(gate, 
                   channels = channels,
                   gate_line_type = gate_line_type,
                   gate_line_width = gate_line_width,
                   gate_line_col = gate_line_col)
  })
  names(gates) <- alias
  
  # RETURN GATE OBJECTS
  return(gates)
  
}
