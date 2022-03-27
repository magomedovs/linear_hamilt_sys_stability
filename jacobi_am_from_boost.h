#ifndef JACOBI_AM
#define JACOBI_AM

#include <cmath>
#include <limits>

template <class T>
T jacobi_am_iter(const T& k, const T& x, T anm1, T bnm1, unsigned N);

template <class T>
T jacobi_am(const T& k, const T& x);

// following function is the modified function "jacobi_recurse" from "boost/math/special_functions/jacobi_elliptic.hpp". Arguments "k" and "x" are swapped
template <class T>
T jacobi_am_iter(const T& k, const T& x, T anm1, T bnm1, unsigned N)
{
   ++N;
   T Tn;
   T cn = (anm1 - bnm1) / 2;
   T an = (anm1 + bnm1) / 2;
   if(cn < std::numeric_limits<T>::epsilon())
   {
      Tn = std::ldexp(T(1), (int)N) * x * an;
   }
   else
      Tn = jacobi_am_iter<T>(k, x, an, std::sqrt(anm1 * bnm1), N);

   return (Tn + std::asin((cn / an) * std::sin(Tn))) / 2;
}

template <class T>
T jacobi_am(const T& k, const T& x) {
    const T k1 = 1 - k;
    const T m = k*k;
    const T m1 = 1 - m;

    const T a = 1;
    const T b = std::sqrt(m1);
    const unsigned N = 0;
    return jacobi_am_iter(k, x, a, b, N);
}

#endif
