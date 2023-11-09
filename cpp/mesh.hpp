#pragma once

#include "basis.hpp"

template <unsigned long dim>
class CartesianUniformMesh {
 public:
  CartesianUniformMesh(const Array<double, dim> &corner_a,
                       const Array<double, dim> &corner_b,
                       const Array<unsigned long, dim> &ns)
      : m_corner_a(corner_a), m_corner_b(corner_b), m_ns(ns) {
    build();
  }

 private:
  void build() {
    m_ncumprod = decltype(m_ncumprod)::Ones();
    m_ncumprodv = decltype(m_ncumprodv)::Ones();
    m_nvert = 1;
    for (decltype(dim) i = 0; i < dim; ++i) {
      m_ncumprod(i + 1) = m_ncumprod(i) * m_ns(i);
      m_ncumprodv(i + 1) = m_ncumprodv(i) * (m_ns(i) + 1);
    }
    m_nvert = m_ncumprodv(dim);
    m_h = (m_corner_b - m_corner_a) / (m_ns.template cast<double>());
  }
  unsigned long m_nelem;
  unsigned long m_nvert;
  Array<unsigned long, dim + 1> m_ncumprod;
  Array<unsigned long, dim + 1> m_ncumprodv;
  Array<double, dim> m_corner_a;
  Array<double, dim> m_corner_b;
  Array<double, dim> m_h;
  Array<unsigned long, dim> m_ns;
};

// todo: cartesian mesh (non-uniform) e.g. geometric refinement near wall
// todo: nd-tree cartesian mesh (non-conforming quad/oct tree refinement)

template <unsigned long dim>
static inline Array<unsigned long, dim> tuple_idx(
    const Array<unsigned long, dim + 1> &cumprod, unsigned long id) {
  Array<unsigned long, dim> tid;
  for (decltype(dim) i = 0; i < dim; ++i)
    tid(i) = (id % cumprod(i + 1)) / cumprod(i);
  return tid;
}

template <unsigned long dim>
unsigned int flat_idx(const Array<unsigned long, dim + 1> &cumprod,
                      const Array<unsigned long, dim> &id) {
  unsigned long fid = id[0];
  for (decltype(dim) i = 1; i < dim; ++i) {
    fid += id(i) * cumprod(i);
  }
  return fid;
}
