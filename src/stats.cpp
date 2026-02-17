#include <Rcpp.h>
#include <algorithm> // std::sort, std::unique
#include <cmath>     // std::log std::exp std::pow, std::abs
#include <numeric>   // std::accumulate std::iota
#include <sstream>   // For std::stringstream

using namespace Rcpp;

// --- MEDIAN ---

template <typename T>
double median_inplace(T& x) {
  size_t n = x.size();
  if (n == 0) {
    return R_NaN;
  }
  size_t half = n / 2;
  if (n % 2 == 1) {
    std::nth_element(x.begin(), x.begin() + half, x.end());
    return x[half];
  } else {
    std::nth_element(x.begin(), x.begin() + half - 1, x.end());
    double val1 = x[half - 1];
    double val2 = *std::min_element(x.begin() + half, x.end());
    return (val1 + val2) / 2.0;
  }
}

// NOTE: THIS METHOD SORTS THE INPUT VECTOR
// [[Rcpp::export]]
double rcpp_median_inplace(Rcpp::NumericVector x) {
  auto it = std::partition(x.begin(), x.end(), [](double val) {
    return !R_IsNA(val);
  });
  // Create a view of the non-NA part of the vector
  Rcpp::NumericVector non_na_view(x.begin(), it);
  // COMPUTE MEDIAN WITHOUT COPY
  return median_inplace(non_na_view);
}

// WE COPY THE INPUT HERE TO PREVENT SORTING
// [[Rcpp::export]]
double median_cpp(Rcpp::NumericVector x) {
  Rcpp::NumericVector y = Rcpp::clone(x); // Clone to avoid modifying original
  // Fast NA removal by partitioning
  auto it = std::partition(y.begin(), y.end(), [](double val) {
    return !R_IsNA(val);
  });
  // Create a view of the non-NA part of the vector
  Rcpp::NumericVector non_na_view(y.begin(), it);
  return rcpp_median_inplace(non_na_view);
}

// [[Rcpp::export]]
Rcpp::NumericVector col_median_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericVector out(ncol);
  Rcpp::NumericVector temp_col(nrow); // Allocate ONCE
  for (int j = 0; j < ncol; ++j) {
    temp_col = x.column(j); // Copy data into reusable vector
    out[j] = rcpp_median_inplace(temp_col);
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- GEOMETRIC MEDIAN ---

// Function to calculate Euclidean distance between two vectors
double euclidean_dist(NumericVector p1, NumericVector p2) {
  return sqrt(sum(pow(p1 - p2, 2.0)));
}

// [[Rcpp::export]]
NumericVector geometric_median_cpp(NumericMatrix x, double eps = 1e-5, int maxiter = 100) {
  int n_rows = x.nrow();
  int n_cols = x.ncol();
  CharacterVector cnames = colnames(x); // Get column names from the input matrix
  
  // 1. Initial guess: the centroid (column means)
  NumericVector y = colMeans(x);
  NumericVector y_new(n_cols);
  
  for (int k = 0; k < maxiter; ++k) {
    double den = 0.0; // Denominator for the update formula
    NumericVector num(n_cols); // Numerator vector
    
    // 2. Iterate through each data point to calculate sums
    for (int i = 0; i < n_rows; ++i) {
      NumericVector xi = x(i, _);
      double dist = euclidean_dist(xi, y);
      
      // Handle the case where the estimate coincides with a data point
      if (dist < 1e-9) { // Use a small threshold for floating point comparison
        xi.names() = cnames; // Set names before returning
        return xi;
      }
      
      den += 1.0 / dist;
      num += xi / dist;
    }
    
    y_new = num / den;
    
    // 3. Check for convergence
    if (euclidean_dist(y_new, y) < eps) {
      y_new.names() = cnames; // Set names before returning
      return y_new;
    }
    
    y = y_new; // Update for the next iteration
  }
  
  // Warning if max iterations reached without convergence
  Rcpp::warning("Algorithm did not converge within the specified number of iterations.");
  y.names() = cnames; // Set names before returning
  return y;
}

// --- GEOMETRIC MEAN ---

// [[Rcpp::export]]
double geomean_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) {
    return R_NaN; // Or some other indicator for an empty vector
  }
  double sum_log_x = 0.0;
  for (int i = 0; i < n; ++i) {
    if (x[i] <= 0) {
      Rcpp::warning(
        "Geometric mean is not defined for non-positive numbers. Returning NaN."
      );
      return R_NaN;
    }
    sum_log_x += std::log(x[i]);
  }
  return std::exp(sum_log_x / n);
}

// [[Rcpp::export]]
Rcpp::NumericVector col_geomean_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  Rcpp::NumericVector out(ncol);
  for (int j = 0; j < ncol; ++j) {
    out[j] = geomean_cpp(x.column(j));
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- STANDARD DEVIATION ---

// [[Rcpp::export]]
double sd_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n < 2) {
    return NA_REAL;
  }
  double sum = std::accumulate(x.begin(), x.end(), 0.0);
  double mean = sum / n;
  double sum_sq_diff = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_sq_diff += (x[i] - mean) * (x[i] - mean);
  }
  return std::sqrt(sum_sq_diff / (n - 1));
}

// [[Rcpp::export]]
Rcpp::NumericVector col_sd_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  Rcpp::NumericVector out(ncol);
  for (int j = 0; j < ncol; ++j) {
    out[j] = sd_cpp(x.column(j));
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- MEDIAN ABSOLUTE DEVIATION ---

// [[Rcpp::export]]
double mad_cpp(Rcpp::NumericVector x) {
  double center = median_cpp(x);
  Rcpp::NumericVector abs_devs(x.size());
  for(int i = 0; i < x.size(); ++i) {
    abs_devs[i] = std::abs(x[i] - center);
  }
  return median_cpp(abs_devs) * 1.4826;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_mad_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericVector out(ncol);
  // Allocate temporary vectors ONCE
  Rcpp::NumericVector temp_col(nrow);
  Rcpp::NumericVector abs_devs(nrow);
  for (int j = 0; j < ncol; ++j) {
    temp_col = x.column(j);
    // Need a copy here because the first median calculation modifies the data
    Rcpp::NumericVector col_copy_for_median = Rcpp::clone(temp_col);
    double center = rcpp_median_inplace(col_copy_for_median);
    for (int i = 0; i < nrow; ++i) {
      abs_devs[i] = std::abs(temp_col[i] - center);
    }
    out[j] = rcpp_median_inplace(abs_devs) * 1.4826;
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- ROBUST STANDARD DEVIATION ---

// [[Rcpp::export]]
double rsd_cpp(Rcpp::NumericVector x) {
  double center = median_cpp(x);
  Rcpp::NumericVector abs_devs(x.size());
  for(int i = 0; i < x.size(); ++i) {
    abs_devs[i] = std::abs(x[i] - center);
  }
  return median_cpp(abs_devs) * 1.4826;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_rsd_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericVector out(ncol);
  // Allocate temporary vectors ONCE
  Rcpp::NumericVector temp_col(nrow);
  Rcpp::NumericVector abs_devs(nrow);
  for (int j = 0; j < ncol; ++j) {
    temp_col = x.column(j);
    // Need a copy here because the first median calculation modifies the data
    Rcpp::NumericVector col_copy_for_median = Rcpp::clone(temp_col);
    double center = rcpp_median_inplace(col_copy_for_median);
    for (int i = 0; i < nrow; ++i) {
      abs_devs[i] = std::abs(temp_col[i] - center);
    }
    out[j] = rcpp_median_inplace(abs_devs) * 1.4826;
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- MEAN ---

// [[Rcpp::export]]
double mean_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) {
    return R_NaN;
  }
  double total = std::accumulate(x.begin(), x.end(), 0.0);
  return total / n;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_mean_cpp(Rcpp::NumericMatrix x) {
  // Use R's optimized colMeans function through an Rcpp call
  Rcpp::NumericVector out = Rcpp::colMeans(x);
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- QUANTILE ---

// [[Rcpp::export]]
Rcpp::NumericVector quantile_cpp(Rcpp::NumericVector x, Rcpp::NumericVector probs) {
  int n = x.size();
  if (n == 0) {
    return Rcpp::NumericVector(probs.size(), R_NaN);
  }
  
  Rcpp::NumericVector y = Rcpp::clone(x);
  Rcpp::NumericVector out(probs.size());
  
  std::vector<int> indices;
  for (int i = 0; i < probs.size(); ++i) {
    if (!R_IsNA(probs[i])) {
      double index = (n - 1) * probs[i];
      indices.push_back(std::floor(index));
      indices.push_back(std::ceil(index));
    }
  }
  
  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  
  // Partially sort the vector by placing the needed indices in their correct positions
  for (int idx : indices) {
    if (idx < n) {
      std::nth_element(y.begin(), y.begin() + idx, y.end());
    }
  }
  
  for (int i = 0; i < probs.size(); ++i) {
    double p = probs[i];
    if (R_IsNA(p) || p < 0 || p > 1) {
      out[i] = R_NaN;
      continue;
    }
    
    double index = (n - 1) * p;
    int lo = std::floor(index);
    int hi = std::ceil(index);
    double h = index - lo;
    
    out[i] = (1.0 - h) * y[lo] + h * y[hi];
  }
  
  // Add names to the output vector
  Rcpp::CharacterVector rn(probs.size());
  for(int i = 0; i < probs.size(); ++i) {
    std::stringstream ss;
    ss << (probs[i] * 100) << "%";
    rn[i] = ss.str();
  }
  out.attr("names") = rn;
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix col_quantile_cpp(Rcpp::NumericMatrix x, Rcpp::NumericVector probs) {
  int ncol = x.ncol();
  Rcpp::NumericMatrix out(probs.size(), ncol);
  
  for (int j = 0; j < ncol; ++j) {
    out.column(j) = quantile_cpp(x.column(j), probs);
  }
  
  Rcpp::colnames(out) = Rcpp::colnames(x);
  Rcpp::CharacterVector rn(probs.size());
  for(int i = 0; i < probs.size(); ++i) {
    std::stringstream ss;
    ss << (probs[i] * 100) << "%";
    rn[i] = ss.str();
  }
  Rcpp::rownames(out) = rn;
  return out;
}

// --- COEFFICIENT OF VARIATION ---

// [[Rcpp::export]]
double cv_cpp(Rcpp::NumericVector x) {
  double m = mean_cpp(x);
  if (m == 0) return NA_REAL;
  return sd_cpp(x) / m;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_cv_cpp(Rcpp::NumericMatrix x) {
  int n_cols = x.ncol();
  Rcpp::NumericVector out(n_cols);
  for (int j = 0; j < n_cols; ++j) {
    out[j] = cv_cpp(x.column(j));
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- ROBUST COEFFICIENT OF VARIATION ---

// [[Rcpp::export]]
double rcv_cpp(Rcpp::NumericVector x) {
  double center = median_cpp(x);
  if (center == 0) return NA_REAL;
  Rcpp::NumericVector abs_devs(x.size());
  for(int i = 0; i < x.size(); ++i) {
    abs_devs[i] = std::abs(x[i] - center);
  }
  double robust_sd = median_cpp(abs_devs) * 1.4826;
  return robust_sd / center;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_rcv_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericVector out(ncol);
  // Allocate temporary vectors ONCE
  Rcpp::NumericVector temp_col(nrow);
  Rcpp::NumericVector abs_devs(nrow);
  for (int j = 0; j < ncol; ++j) {
    temp_col = x.column(j);
    Rcpp::NumericVector col_copy_for_median = Rcpp::clone(temp_col);
    double center = rcpp_median_inplace(col_copy_for_median);
    if (center == 0) {
      out[j] = R_NaN;
      continue;
    }
    for (int i = 0; i < nrow; ++i) {
      abs_devs[i] = std::abs(temp_col[i] - center);
    }
    double robust_sd = rcpp_median_inplace(abs_devs) * 1.4826;
    out[j] = robust_sd / center;
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- RANGE ---

// [[Rcpp::export]]
Rcpp::NumericVector range_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) {
    return Rcpp::NumericVector::create(R_PosInf, R_NegInf);
  }
  double min_val = x[0];
  double max_val = x[0];
  for(int i = 1; i < n; ++i) {
    if (x[i] < min_val) min_val = x[i];
    if (x[i] > max_val) max_val = x[i];
  }
  return Rcpp::NumericVector::create(min_val, max_val);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix col_range_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericMatrix out(2, ncol);
  if (nrow == 0) {
    for (int j = 0; j < ncol; ++j) {
      out(0, j) = R_PosInf;
      out(1, j) = R_NegInf;
    }
  } else {
    for (int j = 0; j < ncol; ++j) {
      Rcpp::NumericVector col = x.column(j);
      double min_val = col[0];
      double max_val = col[0];
      for (int i = 1; i < nrow; ++i) {
        if (col[i] < min_val) min_val = col[i];
        if (col[i] > max_val) max_val = col[i];
      }
      out(0, j) = min_val;
      out(1, j) = max_val;
    }
  }
  Rcpp::rownames(out) = Rcpp::CharacterVector::create("min", "max");
  Rcpp::colnames(out) = Rcpp::colnames(x);
  return out;
}

// --- MIN ---
// faster for large vectors

// [[Rcpp::export]]
double min_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) {
    return R_PosInf;
  }
  double min_val = x[0];
  for(int i = 1; i < n; ++i) {
    if (x[i] < min_val) min_val = x[i];
  }
  return min_val;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_min_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericVector out(ncol);
  
  for (int j = 0; j < ncol; ++j) {
    if (nrow == 0) {
      out[j] = R_PosInf;
    } else {
      Rcpp::NumericVector col = x.column(j);
      double min_val = col[0];
      for (int i = 1; i < nrow; ++i) {
        if (col[i] < min_val) min_val = col[i];
      }
      out[j] = min_val;
    }
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- MAX ---
// faster for large vectors

// [[Rcpp::export]]
double max_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n == 0) {
    return R_NegInf;
  }
  double max_val = x[0];
  for(int i = 1; i < n; ++i) {
    if (x[i] > max_val) max_val = x[i];
  }
  return max_val;
}

// [[Rcpp::export]]
Rcpp::NumericVector col_max_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  int nrow = x.nrow();
  Rcpp::NumericVector out(ncol);
  for (int j = 0; j < ncol; ++j) {
    if (nrow == 0) {
      out[j] = R_NegInf;
    } else {
      Rcpp::NumericVector col = x.column(j);
      double max_val = col[0];
      for (int i = 1; i < nrow; ++i) {
        if (col[i] > max_val) max_val = col[i];
      }
      out[j] = max_val;
    }
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

// --- SKEWNESS ---

// [[Rcpp::export]]
double skewness_cpp(Rcpp::NumericVector x) {
  Rcpp::NumericVector y;
  for(int i = 0; i < x.size(); ++i) {
    if(!R_IsNA(x[i])) {
      y.push_back(x[i]);
    }
  }
  
  int n = y.size();
  if (n < 3) return R_NaN;
  
  double m = mean_cpp(y);
  double m2 = 0.0;
  double m3 = 0.0;
  
  for(int i = 0; i < n; ++i) {
    double dev = y[i] - m;
    m2 += dev * dev;
    m3 += dev * dev * dev;
  }
  
  m2 /= n;
  m3 /= n;
  
  if (m2 == 0) return R_NaN;
  
  return m3 / pow(m2, 1.5);
}

// [[Rcpp::export]]
Rcpp::NumericVector col_skewness_cpp(Rcpp::NumericMatrix x) {
  int ncol = x.ncol();
  Rcpp::NumericVector out(ncol);
  for (int j = 0; j < ncol; ++j) {
    out[j] = skewness_cpp(x.column(j));
  }
  out.attr("names") = Rcpp::colnames(x);
  return out;
}

/// --- SCALE ---

// [[Rcpp::export]]
Rcpp::NumericVector scale_cpp(
    Rcpp::NumericVector x,
    Rcpp::String type = "range",
    Rcpp::Nullable<Rcpp::NumericVector> probs = R_NilValue
) {
  
  std::string type_str = type;
  std::transform(type_str.begin(), type_str.end(), type_str.begin(), ::tolower);
  
  Rcpp::NumericVector out = Rcpp::clone(x);
  int n = out.size();
  
  // Determine scaling type
  bool is_range = (type_str.rfind("r", 0) == 0);
  bool is_mean = (type_str.rfind("mea", 0) == 0);
  bool is_median = (type_str.rfind("med", 0) == 0);
  bool is_quantile = (type_str.rfind("q", 0) == 0);
  bool is_zscore = (type_str.rfind("z", 0) == 0);
  
  // Prepare probs vector for quantile scaling
  Rcpp::NumericVector probs_vec;
  if (is_quantile) {
    if (probs.isNotNull()) {
      probs_vec = Rcpp::NumericVector(probs);
      if (probs_vec.size() != 2) {
        Rcpp::stop("If provided, 'probs' must be a numeric vector of length 2.");
      }
      if (probs_vec[0] < 0 || probs_vec[0] > 1 || probs_vec[1] < 0 || probs_vec[1] > 1 || probs_vec[0] >= probs_vec[1]) {
        Rcpp::stop("'probs' must contain two values between 0 and 1, with the first being smaller than the second.");
      }
    } else {
      probs_vec = Rcpp::NumericVector::create(0.01, 0.99);
    }
  }
  
  if (is_range || is_mean || is_median) {
    Rcpp::NumericVector range_vals = range_cpp(out);
    double min_val = range_vals[0];
    double max_val = range_vals[1];
    double range_diff = max_val - min_val;
    
    if (range_diff == 0) {
      std::fill(out.begin(), out.end(), 0.0);
      return out;
    }
    
    double center = 0.0;
    if (is_range) {
      center = min_val;
    } else if (is_mean) {
      center = mean_cpp(out);
    } else if (is_median) {
      Rcpp::NumericVector col_no_na; // Create fresh vector
      for(int i = 0; i < n; ++i) {
        if(!R_IsNA(out[i])) {
          col_no_na.push_back(out[i]);
        }
      }
      center = median_cpp(col_no_na);
    }
    
    for (int i = 0; i < n; ++i) {
      out[i] = (out[i] - center) / range_diff;
    }
    
  } else if (is_quantile) {
    Rcpp::NumericVector temp_vec = Rcpp::clone(out);
    Rcpp::NumericVector quant_vals = quantile_cpp(temp_vec, probs_vec);
    
    double lower_q = quant_vals[0];
    double upper_q = quant_vals[1];
    double quant_diff = upper_q - lower_q;
    
    if (quant_diff == 0) {
      std::fill(out.begin(), out.end(), 0.5);
      return out;
    }
    
    for (int i = 0; i < n; ++i) {
      double scaled_val = (out[i] - lower_q) / quant_diff;
      if (scaled_val < 0.0) scaled_val = 0.0;
      if (scaled_val > 1.0) scaled_val = 1.0;
      out[i] = scaled_val;
    }
    
  } else if (is_zscore) {
    double sd_val = sd_cpp(out);
    if (sd_val == 0) {
      std::fill(out.begin(), out.end(), 0.0);
      return out;
    }
    double mean_val = mean_cpp(out);
    for (int i = 0; i < n; ++i) {
      out[i] = (out[i] - mean_val) / sd_val;
    }
  } else {
    Rcpp::stop("'type' must be 'range', 'mean', 'median', 'quantile' or 'zscore'!");
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix col_scale_cpp(
    Rcpp::NumericMatrix x,
    Rcpp::String type = "range",
    Rcpp::Nullable<Rcpp::NumericVector> probs = R_NilValue
) {
  
  int ncol = x.ncol();
  Rcpp::NumericMatrix out = Rcpp::clone(x);
  
  for (int j = 0; j < ncol; ++j) {
    Rcpp::NumericVector original_col = x.column(j);
    Rcpp::NumericVector scaled_col = scale_cpp(original_col, type, probs);
    out.column(j) = scaled_col;
  }
  
  Rcpp::colnames(out) = Rcpp::colnames(x);

  return out;
}

// --- RESCALE ---

Rcpp::NumericVector rescale_cpp_internal(
    Rcpp::NumericVector x,
    Rcpp::NumericVector current_scale,
    Rcpp::NumericVector current_limits
) {
  int n = x.size();
  // Create a copy to store the rescaled values
  Rcpp::NumericVector out = Rcpp::clone(x); 
  
  // Use a copy of limits to avoid modifying the input
  Rcpp::NumericVector limits_to_use = Rcpp::clone(current_limits);
  
  // Fill in NA limits from data range
  if (R_IsNA(limits_to_use[0]) || R_IsNA(limits_to_use[1])) {
    Rcpp::NumericVector data_range = range_cpp(out); // Assumes na.rm=true
    if (R_IsNA(limits_to_use[0])) limits_to_use[0] = data_range[0];
    if (R_IsNA(limits_to_use[1])) limits_to_use[1] = data_range[1];
  }
  
  double old_min = limits_to_use[0];
  double old_max = limits_to_use[1];
  double new_min = current_scale[0];
  double new_max = current_scale[1];
  
  double old_range = old_max - old_min;
  double new_range = new_max - new_min;
  
  if (old_range == 0) {
    // If original range is zero, set all values to the new min
    std::fill(out.begin(), out.end(), new_min);
    return out;
  }
  
  for (int i = 0; i < n; ++i) {
    double val = out[i];
    
    // Only process non-NA values
    if (!R_IsNA(val)) {
      // Clamp the value to the specified limits *before* scaling
      if (val < old_min) val = old_min;
      if (val > old_max) val = old_max;
      
      // Perform the rescale operation
      out[i] = new_min + ((val - old_min) / old_range) * new_range;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rescale_cpp(
    Rcpp::NumericVector x,
    Rcpp::Nullable<Rcpp::RObject> scale = R_NilValue,
    Rcpp::Nullable<Rcpp::RObject> limits = R_NilValue
) {
  
  // --- Prepare scale vector ---
  Rcpp::NumericVector current_scale(2);
  Rcpp::NumericVector default_scale = Rcpp::NumericVector::create(0, 1);
  
  if (scale.isNotNull()) {
    current_scale = Rcpp::as<Rcpp::NumericVector>(Rcpp::RObject(scale));
    if (current_scale.size() != 2) {
      Rcpp::stop("'scale' must be a numeric vector of length 2.");
    }
  } else {
    current_scale = default_scale;
  }
  
  // --- Prepare limits vector ---
  Rcpp::NumericVector current_limits(2);
  Rcpp::NumericVector default_limits = Rcpp::NumericVector::create(NA_REAL, NA_REAL);
  
  if (limits.isNotNull()) {
    current_limits = Rcpp::as<Rcpp::NumericVector>(Rcpp::RObject(limits));
    if (current_limits.size() != 2) {
      Rcpp::stop("'limits' must be a numeric vector of length 2.");
    }
  } else {
    current_limits = default_limits;
  }
  
  // --- Call the internal workhorse function ---
  return rescale_cpp_internal(x, current_scale, current_limits);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix col_rescale_cpp(
    Rcpp::NumericMatrix x,
    Rcpp::Nullable<Rcpp::RObject> scale = R_NilValue,
    Rcpp::Nullable<Rcpp::RObject> limits = R_NilValue
) {
  int ncol = x.ncol();
  
  // Create the output matrix as a clone.
  Rcpp::NumericMatrix out = Rcpp::clone(x);
  
  // --- Prepare scale and limits matrices ---
  Rcpp::NumericMatrix scale_mat(2, ncol);
  Rcpp::NumericMatrix limits_mat(2, ncol);
  
  // Handle 'scale' argument
  Rcpp::NumericVector default_scale = Rcpp::NumericVector::create(0, 1);
  if (scale.isNotNull()) {
    Rcpp::RObject scale_r_obj(scale);
    if (Rf_isMatrix(scale_r_obj)) {
      scale_mat = Rcpp::as<Rcpp::NumericMatrix>(scale_r_obj);
      if (scale_mat.nrow() != 2 || scale_mat.ncol() != ncol) {
        Rcpp::stop("If 'scale' is a matrix, it must have 2 rows and ncol(x) columns.");
      }
    } else {
      Rcpp::NumericVector scale_vec = Rcpp::as<Rcpp::NumericVector>(scale_r_obj);
      if (scale_vec.size() != 2) Rcpp::stop("'scale' vector must have length 2.");
      for(int j = 0; j < ncol; ++j) scale_mat.column(j) = scale_vec;
    }
  } else {
    for(int j = 0; j < ncol; ++j) scale_mat.column(j) = default_scale;
  }
  
  // Handle 'limits' argument
  Rcpp::NumericVector default_limits = Rcpp::NumericVector::create(NA_REAL, NA_REAL);
  if (limits.isNotNull()) {
    Rcpp::RObject limits_r_obj(limits);
    if (Rf_isMatrix(limits_r_obj)) {
      limits_mat = Rcpp::as<Rcpp::NumericMatrix>(limits_r_obj);
      if (limits_mat.nrow() != 2 || limits_mat.ncol() != ncol) {
        Rcpp::stop("If 'limits' is a matrix, it must have 2 rows and ncol(x) columns.");
      }
    } else {
      Rcpp::NumericVector limits_vec = Rcpp::as<Rcpp::NumericVector>(limits_r_obj);
      if (limits_vec.size() != 2) Rcpp::stop("'limits' vector must have length 2.");
      for(int j = 0; j < ncol; ++j) limits_mat.column(j) = limits_vec;
    }
  } else {
    for(int j = 0; j < ncol; ++j) limits_mat.column(j) = default_limits;
  }
  // --- End of parameter preparation ---
  
  
  // Loop through columns, rescale, and assign to 'out'
  for (int j = 0; j < ncol; ++j) {
    
    // Call the internal workhorse function
    out.column(j) = rescale_cpp_internal(
      x.column(j), 
      scale_mat.column(j), 
      limits_mat.column(j)
    );
  }
  
  Rcpp::colnames(out) = Rcpp::colnames(x);
  
  // Return the modified copy
  return out;
}

// --- SPECTRAL PURITY (SPARSENESS) ---

// [[Rcpp::export]]
double purity_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  if (n <= 1) {
    return NA_REAL;
  }
  // L1 NORM
  double l1_norm = Rcpp::sum(Rcpp::abs(x));
  // L2 NORM
  double l2_norm = std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));
  // L2 NORM > 0 FOR DIVISION
  if (l2_norm == 0.0) {
    return NA_REAL;
  }
  // HOYER NMF SPARSENESS
  double sparseness = (std::sqrt(static_cast<double>(n)) - (l1_norm / l2_norm)) / 
    (std::sqrt(static_cast<double>(n)) - 1.0);
  
  return sparseness;
}

// [[Rcpp::export]]
Rcpp::NumericVector row_purity_cpp(Rcpp::NumericMatrix x) {
  int n_rows = x.nrow();
  Rcpp::NumericVector out(n_rows);
  // Loop over each row and call the other C++ function
  for (int i = 0; i < n_rows; ++i) {
    out[i] = purity_cpp(x.row(i));
  }
  out.attr("names") = Rcpp::rownames(x);
  return out;
}

// --- VWQR ---

// [[Rcpp::export]]
List vwqr_cpp(DataFrame data, 
              double tau = 0.5, 
              int max_iter = 10, 
              double tol = 1e-6, 
              std::string residual_type = "both", 
              int grid_size = 100,
              double span = 0.75) {
  
  // --- 1. Setup: Load R functions and prepare data ---
  
  Function rq = Environment::namespace_env("quantreg")["rq"];
  Function loess = Environment::namespace_env("stats")["loess"];
  Function predict = Environment::namespace_env("stats")["predict"];
  
  NumericVector y = data["y"];
  NumericVector x = data["x"];
  int n_obs = data.nrows();
  NumericVector weights(n_obs, 1.0);
  
  Formula formula("y ~ x");
  List fit;
  NumericVector coef;
  NumericVector old_coef(2, 0.0);
  
  // --- 2. Iteration Loop ---
  
  for (int i = 0; i < max_iter; ++i) {
    
    // --- Step A: Fit weighted quantile regression ---
    fit = rq(formula,
             _["data"] = data,
             _["tau"] = tau,
             _["weights"] = weights,
             _["method"] = "fn");
    
    coef = as<NumericVector>(fit["coefficients"]);
    
    // --- Step B: Check for convergence ---
    if (i > 0) {
      if (sum(pow(coef - old_coef, 2)) < tol) {
        Rcout << "Converged after " << (i + 1) << " iterations." << std::endl;
        break;
      }
    }
    old_coef = clone(coef);
    
    if (i == max_iter - 1) {
      Rcout << "Reached max_iter (" << max_iter << ") without converging." << std::endl;
    }
    
    // --- Step C: Model the variance of the residuals ---
    NumericVector residuals = y - (coef[0] + coef[1] * x);
    
    // Filter residuals
    std::vector<double> filtered_res_vec;
    std::vector<double> filtered_x_vec;
    
    if (residual_type == "lower") {
      for(int j = 0; j < n_obs; ++j) {
        if(residuals[j] < 0) {
          filtered_res_vec.push_back(std::abs(residuals[j]));
          filtered_x_vec.push_back(x[j]);
        }
      }
    } else if (residual_type == "upper") {
      for(int j = 0; j < n_obs; ++j) {
        if(residuals[j] > 0) {
          filtered_res_vec.push_back(residuals[j]);
          filtered_x_vec.push_back(x[j]);
        }
      }
    } else { // "both"
      for(int j = 0; j < n_obs; ++j) {
        filtered_res_vec.push_back(std::abs(residuals[j]));
        filtered_x_vec.push_back(x[j]);
      }
    }
    
    NumericVector abs_res = wrap(filtered_res_vec);
    NumericVector filtered_x = wrap(filtered_x_vec);
    DataFrame res_data;
    
    // --- Grid-based variance modeling using MEDIAN ---
    if (grid_size > 0 && filtered_x.size() > grid_size) {
      // Define the grid
      double x_min = *std::min_element(filtered_x.begin(), filtered_x.end());
      double x_max = *std::max_element(filtered_x.begin(), filtered_x.end());
      NumericVector x_grid(grid_size);
      // Create a vector of vectors to hold residuals for each grid bin
      std::vector<std::vector<double>> grid_residuals(grid_size);
      
      for(int k = 0; k < grid_size; ++k) {
        x_grid[k] = x_min + k * (x_max - x_min) / (grid_size - 1);
      }
      
      // Assign residuals to their corresponding grid bin
      for(int k = 0; k < filtered_x.size(); ++k) {
        int idx = static_cast<int>(round((filtered_x[k] - x_min) / (x_max - x_min) * (grid_size - 1)));
        if (idx >= 0 && idx < grid_size) {
          grid_residuals[idx].push_back(abs_res[k]);
        }
      }
      
      // Calculate the median residual for each grid point
      NumericVector grid_res_median(grid_size);
      for(int k = 0; k < grid_size; ++k) {
        grid_res_median[k] = median_inplace(grid_residuals[k]);
      }
      
      // Create a data frame from the grid for the LOESS model
      res_data = DataFrame::create(_["abs_res"] = grid_res_median, _["x"] = x_grid);
      
    } else {
      // Use the original (non-gridded) data if grid_size is 0 or data is too small
      res_data = DataFrame::create(_["abs_res"] = abs_res, _["x"] = filtered_x);
    }
    
    List var_mod = loess(Formula("abs_res ~ x"),
                         _["family"] = "symmetric",
                         _["data"] = res_data, 
                         _["span"] = span);
    
    // Predict the local standard deviation for all original points
    NumericVector local_sd = predict(var_mod, _["newdata"] = data);
    
    // --- Step D: Calculate new weights ---
    const double C = 4.685;
    for(int j = 0; j < n_obs; ++j) {
      if (local_sd[j] <= 1e-6) local_sd[j] = 1e-6;
      
      double z_res = residuals[j] / local_sd[j];
      double temp_weight = 1.0 - pow(z_res / C, 2.0);
      weights[j] = (temp_weight > 0) ? pow(temp_weight, 2.0) : 0.0;
    }
  }
  
  // --- 3. Return final model ---
  
  fit.push_back(weights, "weights");
  return fit;
}

// --- 1D WASSERSTEIN DISTANCE ---

// [[Rcpp::export]]
double wasserstein1d_cpp(
    Rcpp::NumericVector a, 
    Rcpp::NumericVector b, 
    double p = 1.0, 
    Rcpp::Nullable<Rcpp::NumericVector> wa_ = R_NilValue, 
    Rcpp::Nullable<Rcpp::NumericVector> wb_ = R_NilValue
) {
  size_t m = a.size();
  size_t n = b.size();
  
  if (m == 0 || n == 0) {
    Rcpp::stop("Input vectors 'a' and 'b' must not be empty.");
  }
  
  // ========================================================================
  // Fast Path for unweighted, equal-sized inputs
  // ========================================================================
  if (wa_.isNull() && wb_.isNull() && m == n) {
    // Use Rcpp::clone to avoid modifying the original R vectors
    Rcpp::NumericVector a_sorted = Rcpp::clone(a).sort();
    Rcpp::NumericVector b_sorted = Rcpp::clone(b).sort();
    
    double total_dist = 0.0;
    for (size_t i = 0; i < m; ++i) {
      total_dist += std::pow(std::abs(a_sorted[i] - b_sorted[i]), p);
    }
    return std::pow(total_dist / m, 1.0 / p);
  }
  
  // ========================================================================
  // Handle Weights (creation, validation, and filtering)
  // ========================================================================
  Rcpp::NumericVector wa, wb;
  
  if (wa_.isNull()) {
    wa = Rcpp::NumericVector(m, 1.0);
  } else {
    wa = Rcpp::NumericVector(wa_);
    // Cast the SIGNED wa.size() to the UNSIGNED type of m
    if (static_cast<size_t>(wa.size()) != m) Rcpp::stop("Weights 'wa' must have the same size as 'a'.");  
  }
  
  if (wb_.isNull()) {
    wb = Rcpp::NumericVector(n, 1.0);
  } else {
    wb = Rcpp::NumericVector(wb_);
    // Cast the SIGNED wb.size() to the UNSIGNED type of n
    if (static_cast<size_t>(wb.size()) != n) Rcpp::stop("Weights 'wb' must have the same size as 'b'.");  
  }
  // Filter out points with zero weight
  std::vector<double> a_filtered_std, wa_filtered_std;
  for (size_t i = 0; i < m; ++i) {
    if (wa[i] > 0) {
      a_filtered_std.push_back(a[i]);
      wa_filtered_std.push_back(wa[i]);
    }
  }
  
  std::vector<double> b_filtered_std, wb_filtered_std;
  for (size_t i = 0; i < n; ++i) {
    if (wb[i] > 0) {
      b_filtered_std.push_back(b[i]);
      wb_filtered_std.push_back(wb[i]);
    }
  }
  
  // Convert std::vector back to Rcpp::NumericVector
  a = Rcpp::NumericVector(a_filtered_std.begin(), a_filtered_std.end());
  wa = Rcpp::NumericVector(wa_filtered_std.begin(), wa_filtered_std.end());
  b = Rcpp::NumericVector(b_filtered_std.begin(), b_filtered_std.end());
  wb = Rcpp::NumericVector(wb_filtered_std.begin(), wb_filtered_std.end());
  
  m = a.size();
  n = b.size();
  
  if (m == 0 || n == 0) return 0.0;
  
  // ========================================================================
  // Sort data and weights together
  // ========================================================================
  std::vector<size_t> orda(m);
  std::iota(orda.begin(), orda.end(), 0);
  std::sort(orda.begin(), orda.end(), [&](size_t i1, size_t i2) { return a[i1] < a[i2]; });
  
  Rcpp::NumericVector a_sorted(m), wa_sorted(m);
  for (size_t i = 0; i < m; ++i) {
    a_sorted[i] = a[orda[i]];
    wa_sorted[i] = wa[orda[i]];
  }
  
  std::vector<size_t> ordb(n);
  std::iota(ordb.begin(), ordb.end(), 0);
  std::sort(ordb.begin(), ordb.end(), [&](size_t i1, size_t i2) { return b[i1] < b[i2]; });
  
  Rcpp::NumericVector b_sorted(n), wb_sorted(n);
  for (size_t i = 0; i < n; ++i) {
    b_sorted[i] = b[ordb[i]];
    wb_sorted[i] = wb[ordb[i]];
  }
  
  // ========================================================================
  // Calculate ECDF jump points (cumulative weights)
  // ========================================================================
  double sum_wa = std::accumulate(wa_sorted.begin(), wa_sorted.end(), 0.0);
  double sum_wb = std::accumulate(wb_sorted.begin(), wb_sorted.end(), 0.0);
  
  std::vector<double> cua(m > 1 ? m - 1 : 0);
  double current_sum_a = 0.0;
  for (size_t i = 0; i < m - 1; ++i) {
    current_sum_a += wa_sorted[i];
    cua[i] = current_sum_a / sum_wa;
  }
  
  std::vector<double> cub(n > 1 ? n - 1 : 0);
  double current_sum_b = 0.0;
  for (size_t i = 0; i < n - 1; ++i) {
    current_sum_b += wb_sorted[i];
    cub[i] = current_sum_b / sum_wb;
  }
  
  // ========================================================================
  // Calculate the integral of the difference between quantile functions
  // ========================================================================
  double total_dist = 0.0;
  
  std::vector<double> u = cua;
  u.insert(u.end(), cub.begin(), cub.end());
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  
  double last_u = 0.0;
  size_t i_a = 0, i_b = 0;
  
  for (double current_u : u) {
    double delta_u = current_u - last_u;
    if (delta_u > 0) {
      total_dist += delta_u * std::pow(std::abs(a_sorted[i_a] - b_sorted[i_b]), p);
    }
    if (i_a < m - 1 && current_u >= cua[i_a]) i_a++;
    if (i_b < n - 1 && current_u >= cub[i_b]) i_b++;
    last_u = current_u;
  }
  
  total_dist += (1.0 - last_u) * std::pow(std::abs(a_sorted[m-1] - b_sorted[n-1]), p);
  
  return std::pow(total_dist, 1.0 / p);
}
