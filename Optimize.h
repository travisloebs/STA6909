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
T gss (T (*f)(T x), T lb, T rb, T epsilon, T (*cc)(T x_0, T x_1, T epsilon)){
  double golden_ratio = (1. + sqrt(5)) / 2. - 1.;
  if (lb >= rb){
    std::cout << "Left bound must be less than the right bound.\n";
    return 0.;
  }

  T (*conv_crit)(T, T, T);
  conv_crit = &(*cc);

  double left  = lb + (rb - lb) * (1. - golden_ratio);
  double right = lb + (rb - lb) * (golden_ratio);
  std::cout << "cc is " << conv_crit(left, right, epsilon) << std::endl;
  while (conv_crit(left, right, epsilon) != 1){
    if ((*f)(left) > (*f)(right))
      rb = right;
    else
      lb = left;

    left  = lb + (rb - lb) * (1. - golden_ratio);
    right = lb + (rb - lb) * (golden_ratio);
    std::cout << "cc is " << conv_crit(left, right, epsilon) << std::endl;
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

// The Newton method for univariate optimization and root finding
template <class T>
T uNewton (T (*g)(T x), T (*g_prime)(T x), T x, T d = 0.01){
  T diff = 999999;
  while (diff > d){
    diff = x;
    x = x - (*g)(x) / (*g_prime)(x);
    diff = abs(x - diff);
  }
  return (x);
}

// The secant method for univariate optimization and root finding
template <class T>
T uSecant (T (*g)(T x), T x, T d = 0.01){
  T x_new = x * 1.01;
  T diff = 999999;
  while (diff > d){
    diff = x_new;
    x_new = x_new - ((x_new - x) / (g(x_new) - g(x))) * g(x_new);
    x = diff;
    diff = abs(x_new - x);
  }
  return (x_new);
}

// Absolute convergence criterion
template <class T>
T abs_cc (T x_0, T x_1, T epsilon){
  std::cout << abs(x_1 - x_0) << std::endl;
  return (abs(x_1 - x_0) < epsilon);
}

// Relative convergence criterion
template <class T>
bool rel_cc (T x_0, T x_1, T epsilon){
  return( abs(x_1 - x_0) / abs(x_0) < epsilon);
}

// Modified relative convergence criterion
template <class T>
bool mrel_cc (T x_0, T x_1, T epsilon){
  return( abs(x_1 - x_0) / (abs(x_0) + epsilon) < epsilon);
}

// Genetic algorithms, Nelder-Mead, simualated annealing

#endif
