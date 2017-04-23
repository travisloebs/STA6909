// [[Rcpp::depends(RcppArmadillo)]]
/* This code is compiled using 'Source' in R, and one sourced it can be used in R
 * the same way a regular R function would be used. For instance, once sourced you
 * can run:
 * 
 * h <- hawkinsH(20, 0.01, 10, 200, 250000); h
 * 
 * and compare the resulting value to that in TABLE 3 of Hawkin's paper.
 * NOTE: Above code takes about 30 seconds to finish.
 */
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double hawkinsH(const int & n, 
                const double & alpha, 
                int start = 3, 
                int sample_size = 200, 
                int nsims = 100000) {
  /* Brief explanation of input parameters:
   * 'n' is the sample size parameter (see tables 2 & 3 in Hawkin's paper)
   * 'alpha' is the hazard rate
   * 'start' is the sample number to start testing for a changepoint
   * 'sample_size' is length of each observation sequence (X1, X2, ..., Xn)
   * 'nsims' is the number of runs to use to estimate h
   */
  
  
  // Think of these as scalars, I only create them as
  // matrices so I can use armadillo operators easier.
  // Names are consistent with variable names in Hawkin's paper.
  mat S(1, 1);
  mat Sn(1, 1);
  mat W(1, 1);
  mat E(1, 1);
  mat maxE(1, 1);
  mat T(1, 1);
  
  // Generate all the random values at once, where column i is the i-th sequence of
  // observationes X_1i, X_2i, ..., X_ni and the j-th row is the row vector of all the
  // j-th observations in each sequence.
  mat X(sample_size, nsims);
  X = randn<mat>(sample_size, nsims);
  
  // Create variable that will be used to track how 
  // far into the sequence we are (start at sample 'start' and 
  // end at sample 'n')
  int nt = start - 1; // minus one because C++ element access starts at zero.

  // Run this loop the required amount of times to iterate from 'start' to 'n'
  while (nt < n){

    // create storage variable for T's
    // Needs to be dynamic because our surviving sample size changes (gets smaller)
    mat T_max(1, nsims);
    T_max.zeros(); // fill with zeros
    
    // Create an xbar vector for vectorized subtraction
    mat Xbar(nt + 1, 1);
    
    // Do this loop for each column in the (surviving) sample
    for (unsigned int i = 0; i < nsims; i++){
      // Create Sn needed for calculation of E
      // (See Hawkin's paper page 7)
      Sn = sum((X.col(i)).rows(0, nt));
      
      // Initialize S at zero
      S(0, 0) = 0;
      
      // fill Xbar values with average
      Xbar.fill(Sn(0,0) / nt);
      
      // Calculate 'W' needed for calculation of T
      W = sum(pow((X.col(i)).rows(0, nt) - Xbar, 2));

      // T is maximized when E is maximized, so
      // start this value at zero
      maxE(0, 0) = 0.;
      
      // Find max E between 1 and nt
      for (unsigned int j = 0; j < nt; j++){
        // Find the maximum T by finding maximum E
        S += X(j, i);
        E = pow(((nt + 1.) * S - (j + 1.) * Sn), 2) / ((nt + 1.) * (j + 1.) * ((nt + 1.) - (j + 1.)));
        if (E(0, 0) > maxE(0, 0)){
          // Store value of max E
          maxE(0, 0) = E(0, 0);
        }
      }
      // Calculate T using formula on page 8
      T = sqrt(((nt + 1.) - 2.) * maxE(0, 0) / (W - maxE(0, 0)));
      // Store this max T
      T_max(0, i) = T(0, 0);
    }
    // Fill any NaN values with zero (just in case)
    T_max.elem(find_nonfinite(T_max)).fill(0.);
    // Sort T_max vector so we can find the value at the 1 - alpha percentile
    mat T_sort = sort(T_max, "ascend", 1);

    // If we're at the stopping point then return 'h'
    if (nt == n - 1)
      return T_sort(0, int(ceil(nsims * (1. - alpha) - 1)));
    
    // Else increment nt
    nt += 1;
    
    // Create index vector that stores the order of T_max
    uvec indices = sort_index(T_max);
    
    // Sort columns of sample matrix X by this index
    X = X.cols(indices);
    
    // Reduce nsims by size of triggered sample
    nsims = int(ceil(nsims * (1. - alpha) - 1));
    
    // Remove triggered/dead samples from X
    X = X.cols(0, nsims);
  }
  
  // If something goes wrong, return -1
  return -1;
}
