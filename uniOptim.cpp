#include <iostream>
#include <vector>

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

std::vector <int> fibonacci(const int N){
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
}

double function(double x);

int main(){
  double lb = -2;
  double rb = 4;
  int N = 999;


  for (int i = 0; i < fibos.size(); i++)
    std::cout << fibos[i] << ", ";
  std::cout << std::endl;
  return 0;
}

double function(double x){
  return (-2 * x * x);
}
