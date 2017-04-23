#ifndef KFDA_RCPP
#define KFDA_RCPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
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

// [[Rcpp::export]]
mat randMat (const int nrows, const int ncols){
  return (randn<mat>(nrows, ncols));
}

// [[Rcpp::export]]
mat mkE_w (const mat & X1, const mat & X2){
  float n1 = X1.n_rows;
  float n2 = X2.n_rows;
  float n = n1 + n2;
  
  return (n1 / n * cov(X1, X1) + n2 / n * cov(X2, X2));
}

// [[Rcpp::export]]
double innr_prod (const mat & X1, const mat & X2){
  return 0;
}

// [[Rcpp::export]]
double linear_kernel(const mat & X1, const mat & X2, double c=0){
  mat inner_product = X1.t() * X2;
  
  return (inner_product(0, 0) + c);
}

typedef double (*kernelPtr) (const mat & X1, const mat & X2, double k_param);

// [[Rcpp::export]]
XPtr <kernelPtr> kernelXPtr(std::string kstr) {
  // Must manually specify return for each kernel
  if (kstr == "linear")
    return (XPtr <kernelPtr> (new kernelPtr (&linear_kernel)));
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

// [[Rcpp::export]]
mat mk1(const int l){
  mat ones(l, 1);
  ones.ones();
  
  return ones;
}

// [[Rcpp::export]]
mat mkP(const int l){
  mat one = mk1(l);
  mat P = mkI(l) - 1. / l * (one * one.t());
  
  return P;
}

// [[Rcpp::export]]
mat mkNull(const int nrows, const int ncols){
  mat null(nrows, ncols);
  null.zeros();
  
  return null;
}

// [[Rcpp::export]]
mat mkN(const int n1, const int n2){
  mat N = join_cols(join_rows(mkP(n1), mkNull(n1, n2)), join_rows(mkNull(n2, n1), mkP(n2)));
  
  return N;
}

// [[Rcpp::export]]
mat m(const int n1, const int n2){
  int n = n1 + n2;
  mat m(n, 1);
  
  
  for (int i = 0; i < n1; i++){
    m(i, 0) = -1. / n1;\
  }
  
  for (int i = 0; i < n2; i++){
    m(i + n1, 0) = 1. / n2;
  }
  
  return m;
}

// [[Rcpp::export]]
double maxKFDR(const mat & X1, const mat & X2, const mat & K, const double & gamma){
  // Dynamically compute n1, n2, and n to reduce number of function parameters
  int n1 = X1.n_rows;
  int n2 = X2.n_rows;
  int n = n1 + n2;
  mat X = join_cols(X1, X2);
  // Create the matrices we will need for final computation
  mat _m = m(n1, n2);
  mat N = mkN(n1, n2);
  // Use the formula on the last line of section 3.2 in
  // the homogeneity paper (ktest_v1-1_bible.pdf)
  mat ratio = n1 * n2 / (gamma * n) * (_m.t() * K * _m - 1. / n * _m.t() * K * N *
    inv(gamma * mkI(n) + 1. / n * N * K * N) *
    N * K * _m);
  return (ratio(0, 0));
}

// [[Rcpp::export]]
double mkd(const int r, const mat & E_w, const double gamma){
  if (r ==1)
    return trace(inv(E_w + gamma) * E_w);
  else if (r ==2)
    return (trace(inv((E_w + gamma) * (E_w + gamma)) * (E_w * E_w)));
  else
    return 0;
}

// [[Rcpp::export]]
double T_hat(const mat & X1, const mat & X2, const mat & K, const double gamma){
  
  double t_hat = (maxKFDR(X1, X2, K, gamma) - mkd(1, mkE_w(X1, X2), gamma)) / (sqrt(2.) * mkd(2, mkE_w(X1, X2), gamma));
  
  return t_hat;
}

// [[Rcpp::export]]
List KCpA(const mat & X, const mat & K, const double & gamma){
  int n = X.n_rows;
  int ncols = X.n_cols;
  
  // We will return a vector contraining all of the ratios from
  // the for loop for diagnostics
  mat ratios(n - 2, 1);
  
  double max = -(std::numeric_limits<double>::infinity());
  
  // The place variable will be used to keep track of the maximum KFDA statistic
  int place = 0;
  
  for (int k = 2; k < n-1; k++){
    // Reshape the segments based on the new k
    mat X1 = X.submat(0, 0, k - 1, ncols - 1);
    mat X2 = X.submat(k, 0, n - 1, ncols - 1);
    // Compute the new max KFDA ratio
    ratios(k - 1, 0) = T_hat(X1, X2, K, gamma);  
    // Record the new value if it is the new maximum ratio
    if (ratios(k - 1, 0) > max){
      max = ratios(k - 1, 0);
      place = k - 1;
    }
  }
  place = place + 1;
  
  return List::create(Named("ratios")=ratios, Named("max")=max, Named("place")=place);
}

// [[Rcpp::export]]
mat tsim(const int & nsims){
  mat T(nsims, 1);
  T.zeros();
  mat X(200, 1);
  mat K(200, 200);
  
  for (int i = 0; i < nsims; i++){
    X = randn<mat>(200, 1);
    K = mkK(X, "linear", 0);
    T(i, 0) = T_hat(X.rows(0, 19), X.rows(20, 199), K, 0);
  }
  
  return T;
}

mat randT(const int & n, const int & m, const int & df);
mat randMVN(const int & n, const int & m, const int & df);

typedef mat (*distPtr) (const int & n, const int & m, const int & df);

// [[Rcpp::export]]
XPtr <distPtr> distXPtr(std::string dstr) {
  // Must manually specify return for each kernel
  if (dstr == "mvnorm")
    return (XPtr <distPtr> (new distPtr (&randMVN)));
  else if (dstr == "mvt")
    return (XPtr <distPtr> (new distPtr (&randT)));
  else if (dstr == "mvcauchy")
    return (XPtr <distPtr> (new distPtr (&randT)));
  else
    return (XPtr <distPtr> (R_NilValue));
}

// [[Rcpp::export]]
cube randCube(const int & n, 
              const int & m, 
              const int & nsims, 
              const std::string & distname){
  XPtr <distPtr> xpdist = distXPtr(distname);
  distPtr mv_dist = *xpdist;
  int df;
  if (distname == "mvt")
    df = 5;
  else if (distname == "mvcauchy")
    df = 1;
  else
    df = 0;
  
  cube X(n, m, nsims);
  for (int i = 0; i < nsims; i++){
    mat temp = mv_dist(n, m, df);
    X.slice(i) = temp;
  }
  return X;
}

// [[Rcpp::export]]
mat randT(const int & n, const int & m, const int & df){
  NumericMatrix X(n, m);
  for (int i = 0; i < m; i++)
    X(_, i) = rt(n*m, df);
  mat X_arma = as<mat>(X);
  return X_arma;
}

// [[Rcpp::export]]
mat randMVN(const int & n, const int & m, const int & df = 0){
  mat mvn = randn<mat>(n, m);
  return mvn;
}

#endif











