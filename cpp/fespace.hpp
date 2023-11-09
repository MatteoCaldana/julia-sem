#pragma once

#include "basis.hpp"
#include "gmath.hpp"

#define ALIGNAS_SIZE 64

template <typename Scalar, unsigned long N>
struct alignas(ALIGNAS_SIZE) cArray {
  alignas(ALIGNAS_SIZE) Scalar data[N];
};

template <typename Scalar, unsigned long N, typename SourceScalar>
cArray<Scalar, N> to_cArray(const ArrayX<SourceScalar> &a) {
  cArray<Scalar, N> ca;
  assert(a.size() == N);
  for (unsigned long i = 0; i < N; ++i) {
    ca.data[i] = gmath::abs(a(i)) < gmath::approx_sub<Scalar>() ? 0 : a(i);
  }
  return ca;
}

template <typename Scalar, unsigned long N>
struct alignas(ALIGNAS_SIZE) cMatrix {
  alignas(ALIGNAS_SIZE) Scalar data[N][N];
};

template <typename Scalar, unsigned long N, typename SourceScalar>
cMatrix<Scalar, N> to_cMatrix(const ArrayXX<SourceScalar> &a) {
  cMatrix<Scalar, N> ca;
  assert(a.rows() == N);
  assert(a.cols() == N);
  for (unsigned long i = 0; i < N; ++i)
    for (unsigned long j = 0; j < N; ++j) {
      ca.data[i][j] =
          gmath::abs(a(i, j)) < gmath::approx_sub<Scalar>() ? 0 : a(i, j);
    }
  return ca;
}

template <typename Scalar, unsigned long Degree,
          typename SourceScalar = long double>
struct alignas(ALIGNAS_SIZE) FeSpace {
  static inline constexpr unsigned long N = Degree;
  static inline constexpr unsigned long Np = N + 1;
  // 1D reference coordinates of the interpolation nodes
  alignas(ALIGNAS_SIZE) static const cArray<Scalar, Np> x;
  // 1D weights for gll quadrature
  alignas(ALIGNAS_SIZE) static const cArray<Scalar, Np> w;
  // 1D derivative of Lagrange interpolants at the interpolation nodes the same
  // of (derlgl, derlgl_v2, transpose(V \ gradV))
  alignas(ALIGNAS_SIZE) static const cMatrix<Scalar, Np> D;
  // 1D Vandermonde matrix of Legendre polynomials
  alignas(ALIGNAS_SIZE) static const cMatrix<Scalar, Np> V;
  //    and their derivative (Nrp x Nrp)
  alignas(ALIGNAS_SIZE) static const cMatrix<Scalar, Np> gradV;
  // interpolation from this element to its 4/8 children
  alignas(ALIGNAS_SIZE) static const cMatrix<Scalar, Np> Ph;
  // interpolation from this element to its 2p version
  alignas(ALIGNAS_SIZE) static const cMatrix<Scalar, Np> Pp;
  // exact 1D Mass matrix (Nrp x Nrp) at gll
  alignas(ALIGNAS_SIZE) static const cMatrix<Scalar, Np> M;
};

template <typename Scalar, unsigned long Degree, typename SourceScalar>
const cArray<Scalar, FeSpace<Scalar, Degree, SourceScalar>::Np>
    FeSpace<Scalar, Degree, SourceScalar>::x =
        to_cArray<Scalar, Np>(xwlgl<SourceScalar>(Degree).first);

template <typename Scalar, unsigned long Degree, typename SourceScalar>
const cArray<Scalar, FeSpace<Scalar, Degree, SourceScalar>::Np>
    FeSpace<Scalar, Degree, SourceScalar>::w =
        to_cArray<Scalar, Np>(xwlgl<SourceScalar>(Degree).second);

template <typename Scalar, unsigned long Degree, typename SourceScalar>
const cMatrix<Scalar, FeSpace<Scalar, Degree, SourceScalar>::Np>
    FeSpace<Scalar, Degree, SourceScalar>::D = to_cMatrix<Scalar, Np>(
        derlgl<SourceScalar>(FeSpace<Scalar, Degree, SourceScalar>::x));

template <typename Scalar, unsigned long Degree, typename SourceScalar>
const cMatrix<Scalar, FeSpace<Scalar, Degree, SourceScalar>::Np>
    FeSpace<Scalar, Degree, SourceScalar>::V = to_cMatrix<Scalar, Np>(
        vandermonde_lgl<SourceScalar>(FeSpace<Scalar, Degree, SourceScalar>::x)
            .first);

template <typename Scalar, unsigned long Degree, typename SourceScalar>
const cMatrix<Scalar, FeSpace<Scalar, Degree, SourceScalar>::Np>
    FeSpace<Scalar, Degree, SourceScalar>::gradV = to_cMatrix<Scalar, Np>(
        vandermonde_lgl<SourceScalar>(FeSpace<Scalar, Degree, SourceScalar>::x)
            .second);

template <typename Scalar, unsigned long Degree, typename SourceScalar>
const cMatrix<Scalar, FeSpace<Scalar, Degree, SourceScalar>::Np>
    FeSpace<Scalar, Degree, SourceScalar>::M = to_cMatrix<Scalar, Np>(
        mass<SourceScalar>(FeSpace<Scalar, Degree, SourceScalar>::V));