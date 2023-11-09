#pragma once

#include <quadmath.h>

#include <numbers>
#include <ostream>
#include <sstream>
#include <unsupported/Eigen/SpecialFunctions>

namespace gmath {
template <typename T>
static inline consteval T pi() {
  if constexpr (std::is_same_v<T, __float128>)
    return M_PIq;
  else
    return std::numbers::pi_v<T>;
}

template <typename T>
static inline consteval T epsilon() {
  if constexpr (std::is_same_v<T, __float128>)
    return FLT128_EPSILON;
  else
    return std::numeric_limits<T>::epsilon();
}

// approximate subnormal
template <typename T>
static inline consteval T approx_sub() {
  return epsilon<T>() / T(2.0);
}

template <typename T>
static inline T cos(T x) {
  if constexpr (std::is_same_v<T, __float128>)
    return cosq(x);
  else
    return std::cos(x);
}

template <typename T>
static inline T abs(T x) {
  if constexpr (std::is_same_v<T, __float128>)
    return fabsq(x);
  else
    return std::abs(x);
}
}  // namespace gmath

namespace Eigen {
namespace internal {
template <>
struct erf_impl<__float128> {
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE __float128 run(__float128 x) { return ::erfq(x); }
};
}  // namespace internal
}  // namespace Eigen

template <typename T, int Size>
using Array = Eigen::Array<T, Size, 1>;

template <typename T>
using ArrayX = Eigen::Array<T, 1, -1>;

template <typename T>
using ArrayXX = Eigen::Array<T, -1, -1>;

std::ostream& operator<<(std::ostream& os, const __float128& r) {
  char buf[128];
  quadmath_snprintf(buf, sizeof(buf), "%.20Qe", r);
  os << buf;
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayX<__float128>& r) {
  for (long int i = 0; i < r.size() - 1; ++i) os << r(i) << ", ";
  os << r(r.size() - 1);
  return os;
}