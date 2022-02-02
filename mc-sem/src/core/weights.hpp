#pragma once

#include <array>

#include "legendre.hpp"
#include "lgl.hpp"

namespace sem {

template <int deg>
std::array<double, deg + 1> compute_alpha() {
  std::array<double, deg + 1> alpha;
  const double k = 2. / (deg * (deg + 1));
  for (int i = 0; i < deg + 1; ++i) {
    const auto leg = sem::legendre[deg](sem::lgl[deg][i]);
    alpha[i] = k / (leg * leg);
  }
  return alpha;
}

template <int deg>
double basis_fun(int i, double x) {
  const auto xi = sem::lgl[deg][i];
  const auto lnpx = sem::legendre_prime[deg](x);
  const auto lnxi = sem::legendre[deg](xi);
  const auto nn1 = deg * (deg + 1);
  return (x * x - 1.) * lnpx / (nn1 * (x - xi) * lnxi);
}

template <int deg>
double basis_fun_prime(int i, double x) {
  const auto xi = sem::lgl[deg][i];
  const auto lnpx = sem::legendre_prime[deg](x);
  const auto lnxi = sem::legendre[deg](xi);
  const auto lnppx = sem::legendre_second[deg](x);
  const auto nn1 = deg * (deg + 1);
  const auto x2 = x * x;
  return -((x2 - 1.) * (xi - x) * lnppx + (2 * xi * x - x2 - 1.) * lnpx) /
         (nn1 * lnxi * (x - xi) * (x - xi));
}

inline double basis_fun(int i, int j) {
  return double(int(i == j));
}

// TODO: compute these values with python
template <int deg>
struct BasisFunPrimeCache {
  BasisFunPrimeCache() {
    for (int i = 0; i < deg + 1; ++i) {
      for (int j = 0; j < deg + 1; ++j) {
        if (i != j) {
          const auto xi = sem::lgl[deg][i];
          const auto xj = sem::lgl[deg][j];
          basis_fun_prime[i][j] =
              sem::legendre[deg](xj) / sem::legendre[deg](xi) / (xj - xi);
        } else if (i == 0) {
          basis_fun_prime[i][j] = -((deg + 1) * deg) / 4.;
        } else if (i == deg) {
          basis_fun_prime[i][j] = ((deg + 1) * deg) / 4.;
        } else {
          basis_fun_prime[i][j] = 0.;
        }
      }
    }
  }

  std::array<std::array<double, deg + 1>, deg + 1> basis_fun_prime;
};

}  // namespace sem