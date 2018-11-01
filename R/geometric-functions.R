#' Lines Intersection from retistruct Package
#' 
#' Determine the intersection of two lines L1 and L2 in two dimensions,
#' using the formula described by Weisstein.
#' 
#' @title Determine intersection between two lines 
#' @param P1 vector containing x,y coordinates of one end of L1
#' @param P2 vector containing x,y coordinates of other end of L1
#' @param P3 vector containing x,y coordinates of one end of L2
#' @param P4 vector containing x,y coordinates of other end of L2
#' @param interior.only boolean flag indicating whether only
#' intersections inside L1 and L2 should be returned.
#' @return Vector containing x,y coordinates of intersection of L1
#' and L2.  If L1 and L2 are parallel, this is infinite-valued.  If
#' \code{interior.only} is \code{TRUE}, then when the intersection
#' does not occur between P1 and P2 and P3 and P4, a vector
#' containing \code{NA}s is returned.
#' @source Weisstein, Eric W. "Line-Line Intersection."
#' From MathWorld--A Wolfram Web Resource.
#' \url{http://mathworld.wolfram.com/Line-LineIntersection.html}
#' @author David Sterratt
#' @examples
#' ## Intersection of two intersecting lines
#' linesIntercept(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
#'
#' ## Two lines that don't intersect
#' linesIntercept(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
#' 
#' @noRd
linesIntercept <- function(P1, P2, P3, P4, interior.only=FALSE) {
  P1 <- as.vector(P1)
  P2 <- as.vector(P2)
  P3 <- as.vector(P3)
  P4 <- as.vector(P4)
  
  dx1 <- P1[1] - P2[1]
  dx2 <- P3[1] - P4[1]
  dy1 <- P1[2] - P2[2]
  dy2 <- P3[2] - P4[2]
  
  D <- det(rbind(c(dx1, dy1),
                 c(dx2, dy2)))
  if (D==0) {
    return(c(Inf, Inf))
  }
  D1 <- det(rbind(P1, P2))
  D2 <- det(rbind(P3, P4))
  
  X <- det(rbind(c(D1, dx1),
                 c(D2, dx2)))/D
  Y <- det(rbind(c(D1, dy1),
                 c(D2, dy2)))/D
  
  if (interior.only) {
    ## Compute the fractions of L1 and L2 at which the intersection
    ## occurs
    lambda1 <- -((X-P1[1])*dx1 + (Y-P1[2])*dy1)/(dx1^2 + dy1^2)
    lambda2 <- -((X-P3[1])*dx2 + (Y-P3[2])*dy2)/(dx2^2 + dy2^2)
    if (!((lambda1>0) & (lambda1<1) &
          (lambda2>0) & (lambda2<1))) {
      return(c(NA, NA))
    }
  }
  return(c(X, Y))
}