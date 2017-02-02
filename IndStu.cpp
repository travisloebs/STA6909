#include <iostream>
#include <vector>
#include "Optimize.h"

using namespace std;

double f(double x);
double g(double x);

int main(){
  double lb = -2;
  double rb = 4;
  int N = 25000;
  double x = 0;

  /* Tests the lattice function */
  //double max = lattice(&f, lb, rb, N);
  //cout << max << endl;

  /* Tests the Golden Section Search function */
  //double max = gss(&f, lb, rb, 0.0000001);
  //cout << max << endl;

  /* Tests the Fibonacci sequence function
  std::vector <int> fibos = fibonacci(N);
  for (int i = 0; i < fibos.size(); i++)
    std::cout << fibos[i] << ", ";
  std::cout << std::endl;
  */

/* Tests the bisection root finding method */
  double root = bisection(&g, lb, rb, 0.0000001);
  cout << root << endl;

  return 0;
}

double f(double x = 0){
  return (-(x - 3.)*(x - 3.) + 4.);
}

double g(double x = 0){
  return ((x - 2.)*(x - 2.)*(x - 2.) + 5.);
}
