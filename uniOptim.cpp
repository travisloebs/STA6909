#include <iostream>

template <class T>
T uvLattice (T (*f)(T x), T lb, T rb, int N){
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

double function(double x);

int main(){
  double lb = -2;
  double rb = 4;
  int N = 25000;

  float max = uvLattice(&function, lb, rb, N);
  std::cout << max << std::endl;
  return 0;
}

double function(double x){
  return (-2 * x * x);
}
