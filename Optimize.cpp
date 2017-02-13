/*******************************************************************************
*
*
*******************************************************************************/
#include "Optimize.h"

template <class T>
T lattice (T (*f)(T x), T lb, T rb, int N);

template <class T>
T gss (T (*f)(T x), T lb, T rb, T epsilon, T (*cc)(T x_0, T x_1, T epsilon));
std::vector <int> fibonacci(const int N);

template <class T>
T bisection (T (*g)(T x), T lb, T rb, T d);

template <typename T> int sgn(T val);

template <class T>
T uNewton (T (*g)(T x), T (*g_prime)(T x), T x, T d);

template <class T>
T uSecant (T (*g)(T x), T x, T d);

template <class T>
T abs_cc (T x_0, T x_1, T epsilon);

template <class T>
bool rel_cc (T x_0, T x_1, T epsilon);

template <class T>
bool mrel_cc (T x_0, T x_1, T epsilon);
