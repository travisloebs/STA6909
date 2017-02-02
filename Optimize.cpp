/*******************************************************************************
*
*
*******************************************************************************/
#include "Optimize.h"

template <class T>
T lattice (T (*f)(T x), T lb, T rb, int N);

template <class T>
T gss (T (*f)(T x), T lb, T rb, T d);

std::vector <int> fibonacci(const int N);

template <class T>
T bisection (T (*g)(T x), T lb, T rb, T d);

template <typename T> int sgn(T val);
