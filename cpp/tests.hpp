#pragma once

#include <string>
#include <iostream>

#include "basis.hpp"

void check(bool b, const std::string &msg) {
  std::cout << "TEST " << msg << " ::: " << (b ? "PASSED" : "FAILED") << std::endl;
  if(!b)
    std::exit(-1);
}

void basis_test() {
  const auto type_test = [](auto &&arg) {
    using Scalar = std::decay_t<decltype(arg)>;
    std::cout << "==============================================" << std::endl;
    std::cout << "Testing float " << sizeof(Scalar) << " bytes" << std::endl;
    for (size_t n = 1; n < 10; ++n) {
      const Scalar TOL = gmath::epsilon<Scalar>() * Scalar(8 << n);
      const auto [x, w] = xwlgl<Scalar>(n);
      const auto x_iter = ixlgl<Scalar>(n);

      const auto d = derlgl<Scalar>(x);
      const auto [V, Vinv] = vandermonde_nodes<Scalar>(x);
      const auto d_vandermonde = derlgl_vandermonde<Scalar>(V, Vinv);
      const auto z = pnleg1<Scalar>(x, n).first;

      const std::string deg = " | deg: " + std::to_string(n);
      std::cout << "Tolerance: " << TOL << std::endl;
      if (n > 1)
        check(z(Eigen::seq(1, n - 1)).abs().maxCoeff() < TOL,
              "Node residual        " + deg);
      check((x - x_iter).abs().maxCoeff() < TOL, "Node comparison      " + deg);
      check((d - d_vandermonde).abs().maxCoeff() < TOL,
            "Derivative comparison" + deg);
    }
  };

  std::apply([&](auto &&...args) { (type_test(args), ...); },
             std::tuple<float, double, long double, __float128>{});
}