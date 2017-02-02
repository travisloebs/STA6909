#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <iostream>
#include <math.h>
#include <vector>
/*
// Uses the lattice method to find the maximum of the function of interest
template <class T>
T lattice (T (*f)(T x), T lb, T rb, int N = 999){
  T iter = double(rb - lb) / N;
  T max = lb;
  T place = lb;
  T max_point = lb;

  for (int i = 0; i < N; i++){
  if ((*f)(place + iter) > (*f)(place)){
    max_point = place + iter;
  }
  else;
    ;
  place += iter;
  }
  return (max_point);
}
*/
 /*
template <class T>
T lattice(T (*f)(T x), const T lb, const T rb, const int N = 988){
  T iter = (rb - lb) / N;
  std::vector <int> fibos = fibonacci(N);
  T points[(fibos.size())];
  for (int i = 0; i < (fibos.size() - 1); i++){
      if (N != (fibos.size() - 1)){
        if (i % 2 == 0)
          points[(i / )]
      }
  }
  return ();
}*/

// The Golden Section Search method of univariate optimization
template <class T>
T gss (T (*f)(T x), T lb, T rb, T d){
  double golden_ratio = (1. + sqrt(5)) / 2. - 1.;
  if (lb >= rb){
    std::cout << "Left bound must be less than the right bound.\n";
    return 0.;
  }
  double left  = lb + (rb - lb) * (1. - golden_ratio);
  double right = lb + (rb - lb) * (golden_ratio);

  while ((right - left) > d){
    if ((*f)(left) > (*f)(right))
      rb = right;
    else
      lb = left;

    left  = lb + (rb - lb) * (1. - golden_ratio);
    right = lb + (rb - lb) * (golden_ratio);
  }
  return ((lb + rb) / 2.);
}


inline std::vector <int> fibonacci(const int N){
  // This requires C++11 to compile
  std::vector <int> fibos = {1, 1};

  if (N < 1){
    std::cout << "N must be at least 1.\n";
    return fibos;
  }

  while (fibos.end()[-1] < N){
    fibos.push_back(fibos.end()[-1] + fibos.end()[-2]);
  }
  return fibos;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// The Bisection Search method of univariate optimization
template <class T>
T bisection (T (*g)(T x), T lb, T rb, T d = 0.01){
  T mid_point;
  while ((rb - lb) > d){
    mid_point = (lb + rb) / 2.;
    if (sgn((*g)(mid_point)) == sgn((*g)(lb)))
      lb = mid_point;
    else if (sgn((*g)(mid_point)) == sgn((*g)(rb)))
      rb = mid_point;
    else
      return (mid_point);
  }
  return ((lb + rb) / 2.);
}

#endif
