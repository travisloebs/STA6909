// [[Rcpp::depends(RcppArmadillo)]]

//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat mkI (int n){
  mat a(n, n); 
  a.zeros();
  
  for (int i = 0; i < n; i++)
      a(i, i) = 1.;
  return a;
}

mat randMat (const int nrows, const int ncols){
  
}

// [[Rcpp::export]]
mat mkE_w (mat X1, mat X2){
  float n1 = X1.n_rows;
  float n2 = X2.n_rows;
  float n = n1 + n2;
  
  return (n1 / n * cov(X1, X1) + n2 / n * cov(X2, X2));
}

// [[Rcpp::export]]
double innr_prod (mat X1, mat X2){
  return 0;
}

// [[Rcpp::export]]
double linear_kernel(mat X1, mat X2, double k_param=0){
  mat inner_product = X1.t() * X2;

  return inner_product(0, 0);
}

typedef double (*kernelPtr) (const mat X1, const mat X2, const double k_param);

// [[Rcpp::export]]
XPtr <kernelPtr> kernelXPtr(std::string kstr) {
  // Must manually specify return for each kernel
  if (kstr == "linear")
    return (XPtr <kernelPtr> (new kernelPtr (&linear_kernel)));
  //if (kstr == "polynomial")
  //  return (XPtr <kernelPtr> (new kernelPtr (&poly_kernel)));
  else
    return (XPtr <kernelPtr> (R_NilValue));
}

// [[Rcpp::export]]
mat mkK (mat X, std::string kernname, double k_param=0){
  XPtr <kernelPtr> xpkern = kernelXPtr(kernname);
  kernelPtr func_kernel = *xpkern;
  
  int n = X.n_rows;
  mat K(n, n);
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      K(i, j) = func_kernel(X.row(i), X.row(j), k_param);
    }
  }
  return K;
}

/*

void print_test(){
  mat A(5, 5); A.fill(10.0);
  
  cout << A(2, 3) << endl;
  cout << A.row(0) << endl;
}
*/
/*** R
# X_1 will be an n1 x m matrix and X_2 will be an n2 x m matrix

*/
