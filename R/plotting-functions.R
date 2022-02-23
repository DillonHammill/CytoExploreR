## .BOXED.LABELS ---------------------------------------------------------------

#' Boxed Labels - Modified plotrix
#'
#' @param x,y  x and y position of the centers of the labels. \code{x} can be a
#'   xy.coords list.
#' @param bg The fill color of the rectangles on which the labels are displayed
#'   (see Details).
#' @param labels Text strings.
#' @param border Whether to draw borders around the rectangles.
#' @param xpad,ypad The proportion of the rectangles to the extent of the text
#'   within.
#' @param srt Rotation of the labels. if 90 or 270 degrees, the box will be
#'   rotated 90 degrees.
#' @param cex Character expansion. See \code{text}.
#' @param adj left/right adjustment. If this is set outside the function, the
#'   box will not be aligned properly.
#' @param xlog Whether the X axis is a log axis.
#' @param ylog Whether the y axis is a log axis.
#' @param alpha.bg Numeric [0,1] controlling the transparency of the background,
#'   set to 0.5 by default.
#' @param ... additional arguments passed to \code{text}.
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#'
#' @importFrom graphics par strwidth strheight rect text
#' @importFrom grDevices col2rgb adjustcolor
#' @importFrom utils modifyList
#'
#' @noRd
.boxed.labels <- function(x,
                          y = NA,
                          labels,
                          bg = ifelse(match(par("bg"), "transparent", 0),
                                      "white", par("bg")
                          ),
                          border = NA,
                          xpad = 1.2,
                          ypad = 1.2,
                          srt = 0,
                          cex = 1,
                          adj = 0.5,
                          xlog = FALSE,
                          ylog = FALSE,
                          alpha.bg = 0.5, ...) {
  border <- NA
  oldpars <- par(c("cex", "xpd"))
  par(cex = cex, xpd = TRUE)
  if(.all_na(y)){
    y <- x
  }
  box.adj <- adj + (xpad - 1) * cex * (0.5 - adj)
  if (srt == 90 || srt == 270) {
    bheights <- strwidth(labels)
    theights <- bheights * (1 - box.adj)
    bheights <- bheights * box.adj
    lwidths <- rwidths <- strheight(labels) * 0.5
  }
  else {
    lwidths <- strwidth(labels)
    rwidths <- lwidths * (1 - box.adj)
    lwidths <- lwidths * box.adj
    bheights <- theights <- strheight(labels) * 0.5
  }
  args <- list(
    x = x, y = y, labels = labels, srt = srt, adj = adj,
    col = ifelse(colSums(col2rgb(bg) * c(1, 1.4, 0.6)) < 350, 
                 "white", 
                 "black")
  )
  args <- modifyList(args, list(...))
  if (xlog) {
    xpad <- xpad * 2
    xr <- exp(log(x) - lwidths * xpad)
    xl <- exp(log(x) + lwidths * xpad)
  }
  else {
    xr <- x - lwidths * xpad
    xl <- x + lwidths * xpad
  }
  if (ylog) {
    ypad <- ypad * 2
    yb <- exp(log(y) - bheights * ypad)
    yt <- exp(log(y) + theights * ypad)
  }
  else {
    yb <- y - bheights * ypad
    yt <- y + theights * ypad
  }
  
  rect(xr,
       yb,
       xl,
       yt,
       col = adjustcolor(col = bg, alpha.f = alpha.bg),
       border = border
  )
  do.call(text, args)
  par(cex = oldpars)
}

## .SPREAD.LABELS --------------------------------------------------------------

#' spread.labs from TeachingDemos
#' 
#' Used internally in cyto_plot to offset overlapping labels in y direction.
#' 
#' @param x x or y co-ordinates to spread out
#' @param mindiff minimum difference between values
#' @param maxiter maximum number of iterations
#' @param stepsize how far to move values with each iteration
#' @param min minimum bound for returned values
#' @param max maximum bound for returned values
#' 
#' @return spread co-ordinates
#' 
#' @author Greg Snow \email{538280@gmail.com}
#' 
#' @noRd
.spread.labels <- function(x, 
                        mindiff, 
                        maxiter=1000, 
                        stepsize=1/10,
                        min=-Inf, 
                        max=Inf) {
  
  unsort <- order(order(x))
  x <- sort(x)
  df <- x[-1] - x[ -length(x) ]
  
  stp <- mindiff * stepsize
  
  i <- 1
  while( any( df < mindiff ) ) {
    tmp <- c( df < mindiff, FALSE )
    if( tmp[1] && (x[1] - stp) < min ) {  # don't move bottom set
      tmp2 <- as.logical( cumprod(tmp) )
      tmp <- tmp & !tmp2
    }
    x[ tmp ] <- x[ tmp ] - stp
    tmp <- c( FALSE, df < mindiff )
    if( tmp[length(tmp)] && (x[length(x)] + stp) > max ) { # don't move top
      tmp2 <- rev( as.logical( cumprod( rev(tmp) ) ) )
      tmp <- tmp & !tmp2
    }
    x[ tmp ] <- x[ tmp] + stp
    
    df <- x[-1] - x[-length(x)]
    i <- i + 1
    if( i > maxiter ) {
      warning("Maximum iterations reached")
      break
    }
  }
  x[unsort]
}

## .line2user ------------------------------------------------------------------

#' @noRd 
.line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}
