/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Algorithm:       M. Koyama, F. Mémoli, V. Robins, K. Turner, "Computation of degree-1 persistent
 *                     homology on larger point-clouds using the Reduced Vietoris-Rips filtration".
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *    - YYYY/MM Author: Description of the modification
 */

/**
 * @file Helpers.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief Shared low-level primitives for the Reduced_rips module: the flat point-cloud view, the squared
 * Euclidean metric, and index sorting.
 */

#ifndef REDUCED_RIPS_HELPERS_H_
#define REDUCED_RIPS_HELPERS_H_

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

#include <boost/version.hpp>
#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_map.hpp>
#else
#include <boost/unordered_map.hpp>
#endif

namespace Gudhi {

namespace reduced_rips {

namespace detail {

// Shared primitives for the point-cloud view, the squared-distance metric, and index sorting.

// Widens a squared search radius by a few ULP so a point lying on the search boundary is still returned
// by the kd-tree radius search. The algorithm re-tests every candidate afterwards. The slack is
// proportional to the radius, and never narrower than double's epsilon: the kd-tree searches in double, so
// a slack expressed in a more precise T (e.g. long double) would round away in the conversion and boundary
// points could be missed.
template <class T>
[[nodiscard]] constexpr T widen_radius(T squared_radius) {
  constexpr long double eps_t = std::numeric_limits<T>::epsilon();
  constexpr long double eps_d = std::numeric_limits<double>::epsilon();
  return squared_radius * (T(1) + (T(8) * static_cast<T>(eps_t > eps_d ? eps_t : eps_d)));
}

// sqrt(3) as a compile-time literal.
inline constexpr double sqrt3 = 1.7320508075688772;

// (0.5 * (2 - sqrt(3)))^2 = (7 - 4*sqrt(3))/4: squared-radius factor of the lens-inscribed ball (Lemma 3.5 of
// the reference paper). A non-empty lens-inscribed ball certifies the lune has a single connected component.
template <class T>
inline constexpr T lens_ball_factor = T(0.25) * (T(2) - T(sqrt3)) * (T(2) - T(sqrt3));

// Flat contiguous point cloud: n points of `dim` coordinates packed row-major in one buffer. point i is
// the pointer `data + i*dim`.
struct Cloud {
  const double* data;
  std::size_t dim, n;
  const double* operator[](std::size_t i) const { return data + (i * dim); }
  [[nodiscard]] std::size_t size() const { return n; }
};

// Squared Euclidean distance between two `dim`-coordinate points. Coordinates are stored as double, but the
// sum of squares is accumulated in T (the filtration type) so its precision reaches the barcode.
template <class T>
[[nodiscard]] T l2_dist_2(const double* a, const double* b, std::size_t dim) {
  T sq_norm = T(0);
  for (std::size_t i = 0; i < dim; ++i) {
    T diff = static_cast<T>(a[i]) - static_cast<T>(b[i]);
    sq_norm += diff * diff;
  }
  return sq_norm;
}

// Returns the `count` indices starting at `first` with the `k` smallest `key` values placed at the front,
// ascending by `key` with ties broken by ascending index. When k >= count the whole range is sorted; otherwise
// only the k smallest are ordered (via partial_sort) and the result is truncated to k. `key(i)` maps an index
// to its (cheap-to-read) sort key. This is the one place the module sorts indices by a distance-like key.
template <class Key>
[[nodiscard]] std::vector<std::size_t> smallest_indices_by(std::size_t first, std::size_t count, std::size_t k,
                                                           Key key) {
  std::vector<std::size_t> result(count);
  std::iota(result.begin(), result.end(), first);
  auto less = [&key](std::size_t x, std::size_t y) {
    const auto kx = key(x), ky = key(y);
    return kx != ky ? kx < ky : x < y;
  };
  if (k >= count) {
    std::sort(result.begin(), result.end(), less);
  } else {
    std::partial_sort(result.begin(), result.begin() + static_cast<std::ptrdiff_t>(k), result.end(), less);
    result.resize(k);
  }
  return result;
}

// Drop the entries <= i from `v` in place, preserving the order of the kept (above-i) indices.
inline void keep_above(std::size_t i, std::vector<std::size_t>& v) {
  v.erase(std::remove_if(v.begin(), v.end(), [i](std::size_t nb) { return nb <= i; }), v.end());
}

// Open-addressing hash map used for the packed-edge lookup tables. boost::unordered_flat_map (Boost >= 1.81)
// stores entries in a contiguous bucket array, giving markedly better cache behaviour than the node-based
// std::unordered_map for these hot, lookup-dominated tables; older Boost falls back to boost::unordered_map.
#if BOOST_VERSION >= 108100
template <class K, class V>
using Edge_map = boost::unordered_flat_map<K, V>;
#else
template <class K, class V>
using Edge_map = boost::unordered_map<K, V>;
#endif

}  // namespace detail

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_HELPERS_H_
