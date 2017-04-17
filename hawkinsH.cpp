// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double hawkinsH(const int & n, 
                const double & alpha, 
                int start = 3, 
                int sample_size = 200, 
                int nsims = 9999) {
  
  // Set up scalars
  mat S(1, 1);
  mat Sn(1, 1);
  mat W(1, 1);
  
  // Generate all the random values at once
  mat X(sample_size, nsims);
  X = randn<mat>(sample_size, int(ceil(nsims)));
  
  // More scalars
  mat E(1, 1);
  mat maxE(1, 1);
  mat T(1, 1);
  
  // Create variable that will be used to track how 
  // far into the sequence we are
  int nt = start - 1;

  // Run this loop the required amount of times to iterate from starting
  // 'n' to target 'n'
  while (nt < n){
    cout << X.n_cols << endl;
    // create storage variable for T's
    mat T_max(1, nsims);
    T_max.zeros();
    
    // Create an xbar vector for vectorized subtraction
    mat Xbar(nt + 1, 1);
    
    // Do this loop for each column in the surviving sample
    for (unsigned int i = 0; i < nsims; i++){
      // Create Sn needed for calculation of E
      Sn = sum((X.col(i)).rows(0, nt));
      
      // Initialize S at zero
      S(0, 0) = 0;
      
      // fill Xbar values
      Xbar.fill(Sn(0,0) / nt);
      
      // Calculate 'W' needed for calculation of T
      W = sum(pow((X.col(i)).rows(0, nt) - Xbar, 2));

      // Start this value at zero
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
      T = sqrt(((nt + 1.) - 2.) * maxE(0, 0) / (W - maxE(0, 0)));
      T_max(0, i) = T(0, 0);
    }
    T_max.elem(find_nonfinite(T_max)).fill(0.);
    mat T_sort = sort(T_max, "ascend", 1);

    if (nt == n - 1)
      return T_sort(0, int(ceil(nsims * (1. - alpha) - 1)));
    
    // Increment nt
    nt += 1;
    
    uvec indices = sort_index(T_max);
    X = X.cols(indices);
    nsims = int(ceil(nsims * (1. - alpha) - 1));
    X = X.cols(0, nsims);
  }
  return -1;
}
