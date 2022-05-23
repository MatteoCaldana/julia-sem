#include "core/legendre.hpp"
#include "core/lgl.hpp"
#include "core/weights.hpp"

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>      

void print(std::vector<double> &v) {
  for(const auto x : v)
    std::cout << x << ", ";
  std::cout << "\n";
}

int main() {
  for(int deg = 2; deg < 20; deg++) {
    double err = 0;
    for(int i = 1; i < deg; ++i) {
      err += sem::legendre_prime[deg](sem::lgl[deg][i]) ;
    }
    std::cout << err << std::endl;
  }

  constexpr int deg = 8;
  constexpr int n = deg + 1;
  const auto alphas = sem::compute_alpha<deg>();
  const auto cache = sem::BasisFunPrimeCache<deg>();

  std::vector<double> rhs(n), mat(n*n);  

  for(int i = 0; i < deg + 1; ++i) {
    for(int j = 0; j < deg + 1; ++j) {
      mat[i + j*n] = 0.;
      for(int k = 0; k < deg + 1; ++k) {
        const auto xk = sem::lgl[deg][k];
        const auto aij = 
            cache.basis_fun_prime[i][k] 
          * cache.basis_fun_prime[j][k] 
          * alphas[k];
        mat[i + j*n] += aij; 
        
        if(j == 0) {
          rhs[i] = 0.;//M_PI * M_PI * std::cos(M_PI*xk) * sem::basis_fun<deg>(i, k) * alphas[k];
        }
      }
    }
  }

  std::cout << std::setprecision(16) << std::scientific;
  print(rhs);
  print(mat);

  return 0;
}