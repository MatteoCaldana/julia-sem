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

int main() {
  for(int deg = 2; deg < 20; deg++) {
    double err = 0;
    for(int i = 1; i < deg; ++i) {
      err += sem::legendre_prime[deg](sem::lgl[deg][i]) ;
    }
    std::cout << err << std::endl;
  }

  constexpr int deg = 2;
  const auto alphas = sem::compute_alpha<deg>();

  Vector<double> x, f;
  FullMatrix<double> A

  for(int i = 0; i < deg; ++i) {
    for(int j = 0; j < deg; ++j) {
      //A[i][j] = 0.;
      for(int k = 0; k < deg; ++k) {
        const auto xk = sem::lgl[deg][k];
        //A[i][j] += diff(xk) * sem::basis_fun_prime<deg>(i, xk) * sem::basis_fun_prime<deg>(j, xk) * alpha[k]
        //         + reac(xk) * sem::basis_fun<deg>(i, xk) * sem::basis_fun<deg>(j, xk) * alpha[k];
        if(j == 0) {
          //f[i] = forc(xk) * sem::basis_fun<deg>(i, xk);
        }
      }
    }
  }

  SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(A, x, f, PreconditionIdentity());


  return 0;
}