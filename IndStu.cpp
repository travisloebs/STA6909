#include <iostream>
#include <vector>
#include "Optimize.h"

using namespace std;

double f(double x);
double g(double x);
double g_prime(double x);

int main(){
  double lb = -10;
  double rb = 10;
  int N = 25000;
  double x = 0;

  /* Tests the lattice function */
  //double max = lattice(&f, lb, rb, N);
  //cout << max << endl;

  /* Tests the Golden Section Search function */
  double max = gss(&f, lb, rb, 0.000000001, &(abs_cc));
  cout << max << endl;


/* Tests the bisection root finding method */
  //double root = bisection(&g, lb, rb, 0.000001);
  //cout << root << endl;

  //root = uNewton(&g, &(g_prime), 10000., .00001);
  //cout << root << endl;

  //root = uSecant(&g, 10000., .00001);
  //cout << root << endl;
  return 0;
}

double f(double x = 0){
  return (-(2.* x + 6.6)*(2. * x + 6.6) + 4.);
}

double g(double x = 0){
  return (-4.*(2. * x + 6.6));
}

double g_prime(double x = 0){
  return(-8.);
}
