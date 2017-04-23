// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "/Volumes/CompStats/STA6909/Indep Study/KFDA_Rcpp.h"
#include <cmath>
using namespace Rcpp;
using namespace arma;

mat mkI (int n);
mat randMat (const int nrows, const int ncols);
mat mkE_w (const mat & X1, const mat & X2);
double innr_prod (const mat & X1, const mat & X2);
double linear_kernel(const mat & X1, const mat & X2, double c);
typedef double (*kernelPtr) (const mat & X1, const mat & X2, double k_param);
XPtr <kernelPtr> kernelXPtr(std::string kstr);
// [[Rcpp::export]]
mat mkK (mat X, std::string kernname, double k_param);
mat mk1(const int l);
mat mkP(const int l);
mat mkNull(const int nrows, const int ncols);
mat mkN(const int n1, const int n2);
mat m(const int n1, const int n2);
double maxKFDR(const mat & X1, const mat & X2, const mat & K, const double & gamma);
double mkd(const int r, const mat & E_w, const double gamma);
double T_hat(const mat & X1, const mat & X2, const mat & K, const double gamma);
List KCpA(const mat & X, const mat & K, const double & gamma);
mat tsim(const int & nsims);
typedef mat (*distPtr) (const int & n, const int & m, const int & df);
XPtr <distPtr> distXPtr(std::string dstr);
cube randCube(const int & n, 
              const int & m, 
              const int & nsims, 
              const std::string & distname);
// [[Rcpp::export]]
mat randT(const int & n, const int & m, const int & df);
mat randMVN(const int & n, const int & m, const int & df);

// [[Rcpp::export]]
double T_max(const mat & X, const mat & K, const double gamma, const int & start_row){
  mat T_temp(X.n_rows - start_row - 1, 1);
  for (int i = 0; i < X.n_rows - start_row - 1; i++){
    T_temp(i, 0) = T_hat(X.rows(0, start_row - 1 + i), X.rows(start_row + i, X.n_rows - 1), K, gamma);
  }
  mat T_sort = sort(T_temp);
  return (T_sort(T_temp.n_rows - 1, 0));
}

// [[Rcpp::export]]
double kfdaH(const int & n, 
             const int  & m, 
             const double & alpha, 
             const int & start = 3, 
             const int & sample_size = 20, 
             const int & nsims = 10,
             const std::string & distname = "mvnorm") {
  
  // Think of these as scalars, I only create them as
  // matrices so I can use armadillo operators easier.
  // Names are consistent with variable names in Hawkin's paper.
  mat T(1, 1);
  mat maxT(1, 1);

  // Generate all the random values at once, where column i is the i-th sequence of
  // observationes X_1i, X_2i, ..., X_ni and the j-th row is the row vector of all the
  // j-th observations in each sequence.
  cube X = randCube(n, m, distname);
  
  
  // Create variable that will be used to track how 
  // far into the sequence we are (start at sample 'start' and 
  // end at sample 'n')
  int nt = start - 1; // minus one because C++ element access starts at zero.
  
  // Run this loop the required amount of times to iterate from 'start' to 'n'
  while (nt < n){
    
    // create storage variable for T's
    // Needs to be dynamic because our surviving sample size changes (gets smaller)
    mat T_max(1, nsims);
    
    // Do this loop for each column in the (surviving) sample
    for (unsigned int i = 0; i < nsims; i+=m){
      // Create storage 'scalar' for max T
      maxT(0, 0) = 0.;
      
      // Find max T between 1 and nt
      for (unsigned int j = 0; j < nt; j++){
        T = T_hat;
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
