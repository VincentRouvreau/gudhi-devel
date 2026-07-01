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

#ifndef REDUCED_RIPS_CLOUD_H_
#define REDUCED_RIPS_CLOUD_H_

#include <algorithm>
#include <cstddef>
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

// Shared primitives for the point-cloud view, the squared-distance metric, the edge/2-simplex
// encoding, the lune-membership test, and the heap ordering.

// Tolerance added to squared search radii so a point lying on the search boundary is still
// returned by the kd-tree radius search. The algorithm re-tests every candidate afterwards.
constexpr double epsilon = 1e-14;

// (0.5 * (2 - sqrt(3)))^2 = (7 - 4*sqrt(3))/4: squared-radius factor of the lens-inscribed ball
// (Lemma 3.5 of the reference paper). sqrt(3) written as a literal so the value is a compile-time
// constant (a non-empty lens-inscribed ball certifies the lune has a single connected component).
constexpr double lens_ball_factor = 0.25 * (2.0 - 1.7320508075688772) * (2.0 - 1.7320508075688772);

// Flat contiguous point cloud: n points of `dim` coordinates packed row-major in one buffer. point i is
// the pointer `data + i*dim`.
struct Cloud {
  const double* data;
  std::size_t dim, n;
  const double* operator[](std::size_t i) const { return data + (i * dim); }
  [[nodiscard]] std::size_t size() const { return n; }
};

// Squared Euclidean distance between two `dim`-coordinate points.
inline double l2_dist_2(const double* a, const double* b, std::size_t dim) {
  double sq_norm = 0.0;
  for (std::size_t i = 0; i < dim; ++i) {
    double diff = a[i] - b[i];
    sq_norm += diff * diff;
  }
  return sq_norm;
}

// Returns the indices > ver_idx, sorted by squared distance from ver_idx (ties broken by ascending index).
inline std::vector<std::size_t> find_all_neighbors(std::size_t ver_idx, const Cloud& pm) {
  std::size_t number_of_idx = pm.n - ver_idx - 1;
  const std::size_t off = ver_idx + 1;
  std::vector<double> dist(number_of_idx);
  for (std::size_t i = off; i < pm.n; ++i) dist[i - off] = l2_dist_2(pm[ver_idx], pm[i], pm.dim);
  // Sort the absolute indices directly, looking distances up by index.
  std::vector<std::size_t> result(number_of_idx);
  std::iota(result.begin(), result.end(), off);
  std::sort(result.begin(), result.end(), [&dist, off](std::size_t x, std::size_t y) {
    return dist[x - off] != dist[y - off] ? dist[x - off] < dist[y - off] : x < y;
  });
  return result;
}

// Append the entries of `nearest` that lie above i (index > i) to `dst`, preserving their order.
inline void append_neighbors_above(std::size_t i, const std::vector<std::size_t>& nearest,
                                   std::vector<std::size_t>& dst) {
  for (std::size_t nb : nearest)
    if (nb > i) dst.push_back(nb);
}

// Drop the entries <= i from `v` in place, preserving the order of the kept (above-i) indices.
inline void keep_above(std::size_t i, std::vector<std::size_t>& v) {
  v.erase(std::remove_if(v.begin(), v.end(), [i](std::size_t nb) { return nb <= i; }), v.end());
}

// Lexicographic order of the sorted index pairs {x,y} < {p,q}. Used for the lune-boundary tie-break.
inline bool sorted_pair_less(std::size_t x, std::size_t y, std::size_t p, std::size_t q) {
  std::size_t xlo = std::min(x, y), xhi = std::max(x, y);
  std::size_t plo = std::min(p, q), phi = std::max(p, q);
  return xlo != plo ? xlo < plo : xhi < phi;
}

// True if a point k lies in the lune of edge (a,b) at threshold `thresh`: within `thresh` of both endpoints,
// where a point sitting exactly on a boundary is admitted only when the index tie-break assigns it to this
// edge. dist_ka and dist_kb are the distances from k to endpoints a and b (the paper's d(x,y) and d(x,z)) and
// thresh is the edge length (the paper's r).
inline bool in_lune(double dist_ka, double dist_kb, double thresh, std::size_t a, std::size_t b, std::size_t k) {
  if (dist_ka < thresh && dist_kb < thresh) return true;
  if (dist_ka == thresh && dist_kb < thresh) return sorted_pair_less(a, k, a, b);
  if (dist_ka < thresh && dist_kb == thresh) return sorted_pair_less(b, k, a, b);
  if (dist_ka == thresh && dist_kb == thresh) return sorted_pair_less(a, k, a, b) && sorted_pair_less(b, k, a, b);
  return false;
}

// 1-simplices are keyed by a single integer lo*n + hi (lo < hi).
inline std::size_t pack_edge(std::size_t lo, std::size_t hi, std::size_t n) { return (lo * n) + hi; }

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

#endif  // REDUCED_RIPS_CLOUD_H_
