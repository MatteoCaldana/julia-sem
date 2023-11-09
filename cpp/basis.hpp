#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>
#include <numbers>

#include "gmath.hpp"

template <typename T>
ArrayX<T> pnleg(const ArrayX<T> &x, size_t N) {
  const auto n = x.size();
  if (N == 0) return ArrayX<T>::Ones(n);
  ArrayX<T> p1 = ArrayX<T>::Ones(n), p2 = x, p3 = x;
  for (decltype(N) k = 1; k < N; ++k) {
    p3 = (T(2 * k + 1) * x * p2 - T(k) * p1) / T(k + 1);
    p1 = p2;
    p2 = p3;
  }
  return p3;
}

template <typename T>
std::pair<ArrayX<T>, ArrayX<T>> pnleg1(const ArrayX<T> &x, size_t N) {
  const auto n = x.size();
  if (N == 0) return {ArrayX<T>::Zero(n), ArrayX<T>::Ones(n)};
  if (N == 1) return {ArrayX<T>::Ones(n), x};
  ArrayX<T> p1 = ArrayX<T>::Ones(n), p2 = x, p3 = x;
  ArrayX<T> p11 = ArrayX<T>::Zero(n), p21 = ArrayX<T>::Ones(n), p31 = x;
  for (decltype(N) k = 1; k < N; ++k) {
    const auto k2p1 = T(2 * k + 1);
    const auto kp1 = T(1) / T(k + 1);
    p3 = (k2p1 * x * p2 - k * p1) * kp1;
    p31 = (k2p1 * (x * p21 + p2) - k * p11) * kp1;
    p11 = p21;
    p21 = p31;
    p1 = p2;
    p2 = p3;
  }
  return {p31, p3};
}

template <typename T>
std::tuple<ArrayX<T>, ArrayX<T>, ArrayX<T>> pnleg2(const ArrayX<T> &x,
                                                   size_t N) {
  const auto n = x.size();
  if (N <= 1)
    return {ArrayX<T>::Zero(n), ArrayX<T>::Zero(n), ArrayX<T>::Zero(n)};
  if (N == 2)
    return {ArrayX<T>::Zero(n), ArrayX<T>::Zero(n), ArrayX<T>::Ones(n)};
  ArrayX<T> p1 = ArrayX<T>::Ones(n), p2 = x, p3 = x;
  ArrayX<T> p11 = ArrayX<T>::Zero(n), p21 = ArrayX<T>::Ones(n), p31 = x;
  ArrayX<T> p12 = ArrayX<T>::Zero(n), p22 = ArrayX<T>::Zero(n), p32 = x;

  for (decltype(N) k = 1; k < N; ++k) {
    const auto k2p1 = T(2 * k + 1);
    const auto kp1 = T(1) / T(k + 1);
    p3 = (k2p1 * x * p2 - k * p1) * kp1;
    p31 = (k2p1 * (x * p21 + p2) - k * p11) * kp1;
    p32 = (k2p1 * (x * p22 + p21 * 2) - k * p12) * kp1;
    p12 = p22;
    p22 = p32;
    p11 = p21;
    p21 = p31;
    p1 = p2;
    p2 = p3;
  }
  return {p32, p31, p3};
}

template <typename T>
ArrayX<T> ixlgl(size_t N) {
  constexpr T TOL = T(10.0) * gmath::epsilon<T>();
  constexpr size_t MAX_ITER = 4 * sizeof(T);
  ArrayX<T> x(N - 1);
  if (N > 1) {
    for (size_t i = 0; i < N - 1; ++i)
      x(i) = gmath::cos(gmath::pi<T>() * T(i + 1) / T(N));
    for (size_t i = 0; i < MAX_ITER; ++i) {
      const auto [p2, p1, p] = pnleg2<T>(x, N);
      ArrayX<T> dx = (p1 * p2 != T(0)).select(p1 / p2, 0);
      x -= dx;
      if (dx.abs().maxCoeff() < TOL) break;
    }
  }

  ArrayX<T> x11(N + 1);
  x11(N) = T(-1);
  x11(0) = T(1);
  for (size_t i = 0; i < N - 1; ++i) x11(i + 1) = x(i);
  return x11.reverse();
}

template <typename T>
std::pair<ArrayX<T>, ArrayX<T>> jacobi_eval(const ArrayX<T> &x, size_t N,
                                            T alpha, T beta) {
  const auto n = x.size();
  const auto apb = alpha + beta;
  const auto ab2 = alpha * alpha - beta * beta;

  if (N == 0) return {ArrayX<T>::Ones(n), ArrayX<T>::Zero(n)};

  ArrayX<T> p1 = ArrayX<T>::Ones(n);
  ArrayX<T> p2 = ArrayX<T>::Zero(n);
  ArrayX<T> pd1 = ArrayX<T>::Zero(n);
  ArrayX<T> pd2 = ArrayX<T>::Zero(n);
  ArrayX<T> p = (alpha - beta + (apb + T(2)) * x) * T(0.5);
  ArrayX<T> pd = T(0.5) * (apb + T(2)) * p1;

  for (decltype(N) k = 1; k < N; ++k) {
    const T k2 = k * 2;
    const T k2ab = k2 + alpha + beta;
    const T k1 = k + 1;
    const T k2ab1 = k2ab + 1;
    const T k2ab2 = k2ab1 + 1;
    const T a1 = 2 * k1 * (k1 + apb) * k2ab;
    const T a21 = k2ab1 * ab2;
    const T a22 = k2ab2 * k2ab1 * k2ab;
    const T a3 = 2 * (k + alpha) * (k + beta) * k2ab2;
    std::swap(p2, p1);
    std::swap(p1, p);
    std::swap(pd2, pd1);
    std::swap(pd1, pd);
    p = ((a21 + a22 * x) * p1 - a3 * p2) / a1;
    pd = (a22 * p1 + (a21 + a22 * x) * pd1 - a3 * pd2) / a1;
  }

  return {p, pd};
}

template <typename T>
ArrayX<T> jacobi_roots(size_t N, T alpha, T beta) {
  constexpr T TOL = T(10.0) * gmath::epsilon<T>();
  constexpr size_t MAX_ITER = 4 * sizeof(T);
  constexpr T PI = gmath::pi<T>();
  ArrayX<T> x = ArrayX<T>::Zero(N);
  ArrayX<T> x0(1, 1);
  x0(0) = gmath::cos(PI / T(2 * N));
  T x1 = x0(0);
  for (size_t j = 0; j < N; ++j) {
    T diff = TOL * T(2);
    size_t kiter = 0;
    while (kiter <= MAX_ITER && diff >= TOL) {
      const auto [p, pd] = jacobi_eval(x0, N, alpha, beta);
      const T ss =
          j > 0 ? (T(1) / (x0(0) - x(Eigen::seq(0, j - 1)))).sum() : T(0);
      x1 = x0(0) - p(0) / (pd(0) - ss * p(0));
      diff = gmath::abs(x1 - x0(0));
      kiter++;
      x0(0) = x1;
    }
    x0(0) = (x1 + gmath::cos(T(2 * (j + 2) - 1) * PI / T(2 * N))) / T(2);
    x(j) = x1;
  }
  return x;
}

template <typename T>
std::pair<ArrayX<T>, ArrayX<T>> xwlgl(size_t N) {
  ArrayX<T> x = ArrayX<T>::Zero(N + 1), w = ArrayX<T>::Ones(N + 1);
  x(N) = T(-1);
  x(0) = T(1);
  if (N > 1) {
    x(Eigen::seq(1, N - 1)) = jacobi_roots(N - 1, T(1), T(1));
    const T coef = T(2) / T(N * (N + 1));
    w = coef / (pnleg(x, N) * pnleg(x, N));
  }
  return {x.reverse(), w};
}

template <typename T>
ArrayXX<T> derlgl(const ArrayX<T> &x) {
  const auto np = x.size();
  const auto n = np - 1;
  using idx_t = std::decay_t<decltype(np)>;
  ArrayXX<T> d = ArrayXX<T>::Zero(np, np);
  const auto lnx = pnleg(x, n);
  for (idx_t j = 0; j < np; ++j) {
    for (idx_t i = 0; i < np; ++i) {
      if (i != j) {
        d(i, j) = lnx(i) / ((x(i) - x(j)) * lnx(j));
      }
    }
  }
  d(0, 0) = T(-0.25) * T(n * np);
  d(n, n) = -d(0, 0);
  return d;
}

template <typename T>
std::pair<ArrayXX<T>, ArrayXX<T>> vandermonde_nodes(const ArrayX<T> &x) {
  const auto np = x.size();
  using idx_t = std::decay_t<decltype(np)>;
  ArrayXX<T> V = ArrayXX<T>::Ones(np, np);
  for (idx_t i = 0; i < np - 1; ++i) {
    V.col(i + 1) = V.col(i) * x.transpose();
  }
  return {V, V.matrix().inverse().array()};
}

template <typename T>
ArrayXX<T> derlgl_vandermonde(const ArrayXX<T> &V, const ArrayXX<T> &Vinv) {
  const auto np = V.rows();
  using namespace Eigen::placeholders;

  const auto derivative = Vinv(Eigen::seq(1, np - 1), all).colwise() *
                          ArrayX<T>::LinSpaced(np - 1, 1, np - 1).transpose();

  return V(all, Eigen::seq(0, np - 2)).matrix() * derivative.matrix();
}

template <typename T>
std::pair<ArrayXX<T>, ArrayXX<T>> vandermonde_lgl(const ArrayX<T> &x) {
  const auto np = x.size();
  using idx_t = std::decay_t<decltype(np)>;
  ArrayXX<T> V = ArrayXX<T>::Empty(np, np);
  ArrayXX<T> gradV = ArrayXX<T>::Empty(np, np);
  for (idx_t i = 0; i < np; ++i) {
    const auto [p, pd] = jacobi_eval(x, i, T(0.0), T(0.0));
    V.col(i) = p;
    gradV.col(i) = pd;
  }
  return {V, gradV};
}

// template <typename T>
// std::pair<ArrayXX<T>, ArrayXX<T>> pprojection(const ArrayX<T> &x) {
//   const auto np = x.size();

//     x_2p, _ = xwlgl(2 * np - 1)

//   using idx_t = std::decay_t<decltype(np)>;
//   ArrayXX<T> V = ArrayXX<T>::Empty(np, np);

//   return {Vpp, pp};
// }

// template <typename T>
// std::pair<ArrayXX<T>, ArrayXX<T>> hprojection(const ArrayX<T> &x) {
//   const auto np = x.size();
//         x_hby2 = [Float128(0.5) * (x .- Float128(1))...; Float128(0.5) *
//         (x[2:end] .+ Float128(1))...]

//   using idx_t = std::decay_t<decltype(np)>;
//   ArrayXX<T> V = ArrayXX<T>::Empty(np, np);

//   return {Vph, ph};
// }

template <typename T>
ArrayXX<T> mass(const ArrayXX<T> &m) {
  const auto iV = m.matrix().inverse();
  return (iV * iV.transpose()).array();
}
