#include "core/legendre.hpp"
#include "core/lgl.hpp"
#include "core/weights.hpp"

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <iostream>
#include <cmath>

using namespace dealii;

int main() {
  for(int deg = 2; deg < 20; deg++) {
    double err = 0;
    for(int i = 1; i < deg; ++i) {
      err += sem::legendre_prime[deg](sem::lgl[deg][i]) ;
    }
    std::cout << err << std::endl;
  }

  constexpr int deg = 8;
  const auto alphas = sem::compute_alpha<deg>();

  Vector<double> sol(deg + 1), rhs(deg + 1);  
  FullMatrix<double> mat(deg + 1, deg + 1);

  for(int i = 0; i < deg + 1; ++i) {
    for(int j = 0; j < deg + 1; ++j) {
      mat.set(i, j, 0.);
      for(int k = 0; k < deg + 1; ++k) {
        const auto xk = sem::lgl[deg][k];
        const auto aij = 
            sem::basis_fun_prime<deg>(i, xk) 
          * sem::basis_fun_prime<deg>(j, xk) 
          * alphas[k]
        mat.add(i, j, aij); 
        std::cout << aij << std::endl;
        if(j == 0) {
          rhs[i] = M_PI * M_PI * std::cos(M_PI*xk) * sem::basis_fun<deg>(i, xk) * alphas[k];
        }
      }
    }
  }

  SolverControl            solver_control(1000, 1e-6 * rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(mat, sol, rhs, PreconditionIdentity());

  std::cout << rhs << std::endl;
  std::cout << sol << std::endl;


  return 0;
}