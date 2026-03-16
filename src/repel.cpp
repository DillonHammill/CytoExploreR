#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//' Repel overlapping labels
//'
//' Iteratively move labels apart to minimize overlap.
//'
//' @param x x-coordinates
//' @param y y-coordinates
//' @param w label widths
//' @param h label heights
//' @param xlim plot x-limits
//' @param ylim plot y-limits
//' @param iter maximum iterations
//' @param force repulsion force factor
//'
//' @noRd
// [[Rcpp::export]]
NumericMatrix cpp_repel(NumericVector x,
                        NumericVector y,
                        NumericVector w,
                        NumericVector h,
                        NumericVector xlim,
                        NumericVector ylim,
                        int iter = 2000,
                        double force = 1.0,
                        NumericVector pad = NumericVector::create(0.0, 0.0)) {

  int n = x.size();
  NumericVector nx(n);
  NumericVector ny(n);

  // Copy initial positions
  for(int i = 0; i < n; ++i) {
    nx[i] = x[i];
    ny[i] = y[i];
  }

  double xmin = xlim[0];
  double xmax = xlim[1];
  double ymin = ylim[0];
  double ymax = ylim[1];

  // Initial temperature (max displacement)
  double max_dim = std::max(xmax - xmin, ymax - ymin);
  double temperature = max_dim / 10.0;

  for(int k = 0; k < iter; ++k) {

    // Cooling
    double temp = temperature * (1.0 - (double)k / (double)iter);

    for(int i = 0; i < n; ++i) {
      double fx = 0;
      double fy = 0;

      // 1. Attraction to data point (spring force)
      double dx_origin = x[i] - nx[i];
      double dy_origin = y[i] - ny[i];

      // Attraction strength
      fx += dx_origin;
      fy += dy_origin;

      // 2. Repulsion from other labels
      for(int j = 0; j < n; ++j) {
        if(i == j) continue;

        double dx = nx[i] - nx[j];
        double dy = ny[i] - ny[j];

        // Effective dimensions with padding
        double wi = w[i] + pad[0];
        double wj = w[j] + pad[0];
        double hi = h[i] + pad[1];
        double hj = h[j] + pad[1];

        double w_sum = (wi + wj) / 2.0;
        double h_sum = (hi + hj) / 2.0;

        if(std::abs(dx) < w_sum && std::abs(dy) < h_sum) {
          // Overlap detected

          // Add jitter if perfectly aligned
          if(dx == 0 && dy == 0) {
             dx = (R::runif(0,1) - 0.5) * 0.01;
             dy = (R::runif(0,1) - 0.5) * 0.01;
          }

          double d2 = dx*dx + dy*dy;
          double d = std::sqrt(d2);

          // Repulsion vector (normalized)
          double ux = dx / d;
          double uy = dy / d;

          // Repulsion strength: stronger force when closer
          double rep = force * 10.0;

          fx += ux * rep;
          fy += uy * rep;
        }
      }

      // Limit movement by temperature
      double f_mag = std::sqrt(fx*fx + fy*fy);
      if(f_mag > 0) {
        double scale = std::min(f_mag, temp) / f_mag;
        nx[i] += fx * scale;
        ny[i] += fy * scale;
      }

      // Clamp to limits
      if(nx[i] < xmin + w[i]/2.0) nx[i] = xmin + w[i]/2.0;
      if(nx[i] > xmax - w[i]/2.0) nx[i] = xmax - w[i]/2.0;
      if(ny[i] < ymin + h[i]/2.0) ny[i] = ymin + h[i]/2.0;
      if(ny[i] > ymax - h[i]/2.0) ny[i] = ymax - h[i]/2.0;
    }
  }

  NumericMatrix out(n, 2);
  for(int i = 0; i < n; ++i) {
    out(i, 0) = nx[i];
    out(i, 1) = ny[i];
  }
  return out;
}
