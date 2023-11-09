#include <iostream>
#include <quadmath.h>

#include "tests.hpp"
// #include "fespace.hpp"
#include "mesh.hpp"
#include "gmath.hpp"
#include "basis.hpp"

#include <iostream>

#include "basis.hpp"

// g++ main.cpp -I../../eigen -std=c++20 -Wall -Wextra -fext-numeric-literals -lquadmath
int main() {
  basis_test();

  // FeSpace<double, 4> fespace;
  Array<double, 3> oned = Array<double, 3>::Ones();
  Array<unsigned long, 3> onei = Array<unsigned long, 3>::Ones();
  CartesianUniformMesh<3> mesh(-1.0 * oned, oned, 7 * onei);

  return 0;
}