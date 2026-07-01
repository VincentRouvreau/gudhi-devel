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

#ifndef REDUCED_RIPS_GEOMETRY_H_
#define REDUCED_RIPS_GEOMETRY_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

#include <gudhi/Reduced_rips/Cloud.h>
#include <gudhi/Reduced_rips/Euclidean_kd_tree.h>
#include <gudhi/Reduced_rips/Lune_builder.h>
#include <gudhi/Reduced_rips/Relative_neighborhood_graph.h>

namespace Gudhi {

namespace reduced_rips {

using detail::Cloud;
using detail::Edge_map;
using detail::find_all_neighbors;
using detail::keep_above;
using detail::l2_dist_2;

// ---- Geometry policies --------------------------------------------------------------------------------
// The computation core (Reduced_rips::compute_impl) is written once against a geometry policy answering purely
// metric queries: the inter-point distance on the policy's own scale (dist), a mapping from that scale back to
// a true output distance (to_distance), the k nearest points, the higher-indexed neighbors in ascending order,
// the relative-neighborhood-graph edge count, and the per-edge lune computation. The reduction relies only on
// the *ordering* of the dist values, so each policy is free to pick whichever scale it computes most cheaply
// and exactly -- dist is named for that role, not for a fixed unit. Two policies implement it.
// Euclidean_geometry (coordinates + kd-tree) keeps the paper's lens-ball / wide-angle accelerations and works
// in squared distances (dist returns d^2, to_distance is sqrt). Coordinates stay double and CGAL searches in
// double, but each candidate CGAL returns is re-checked with l2_dist_2<T>, so the barcode keeps T's precision.
// Matrix_geometry (a bare symmetric distance matrix) has no coordinates, so it gathers candidates by scanning
// a matrix row and always runs the exact union-find; it carries the supplied distances unchanged (dist returns
// d, to_distance is the identity).

template <class T>
class Euclidean_geometry {
 public:
  using Filtration_value = T;
  Euclidean_geometry(const Cloud& pm, const Euclidean_kd_tree& kd) : pm_(pm), kd_(&kd) {}
  [[nodiscard]] std::size_t size() const { return pm_.n; }
  // Ordering scale: the squared Euclidean distance.
  [[nodiscard]] T dist(std::size_t i, std::size_t j) const { return l2_dist_2<T>(pm_[i], pm_[j], pm_.dim); }
  // Maps a squared distance back to the true distance for the output barcode.
  [[nodiscard]] static T to_distance(T squared) { return std::sqrt(squared); }
  [[nodiscard]] std::vector<std::size_t> nearest(std::size_t i, std::size_t k) const {
    return kd_->nearest_neighbors(pm_[i], k);
  }
  // The k nearest points to i restricted to index > i, ascending by distance (ties by index).
  [[nodiscard]] std::vector<std::size_t> nearest_neighbors_above(std::size_t i, std::size_t k) const {
    std::vector<std::size_t> result = nearest(i, k);
    keep_above(i, result);
    return result;
  }
  [[nodiscard]] std::vector<std::size_t> neighbors_above(std::size_t i) const { return find_all_neighbors(i, pm_); }
  // Early-stop target: the RNG cycle rank (the number of finite H1 bars), from the per-dimension routine.
  [[nodiscard]] std::size_t rng_early_stop_target() const {
    if (pm_.dim == 2) return rng_cycle_rank_delaunay<T>(pm_, *kd_, delaunay_edges_2d);
    if (pm_.dim == 3) return rng_cycle_rank_delaunay<T>(pm_, *kd_, delaunay_edges_3d);
    return rng_cycle_rank_general<T>(pm_, *kd_);
  }
  [[nodiscard]] Lune_result<T> lune(const Batch_edge<T>& e, const Edge_map<std::size_t, std::size_t>& one_simp_to_idx,
                                    std::size_t n) const {
    return Lune_builder<T>(e, one_simp_to_idx, n).build_euclidean(pm_, *kd_);
  }

 private:
  Cloud pm_;                     // non-owning view; pointee outlives this
  const Euclidean_kd_tree* kd_;  // non-owning, never null; kd-tree is move-only so stored by pointer
};

// Dense symmetric matrix of distances, n by n row-major. Unlike the Euclidean policy it keeps the supplied
// distances as-is rather than squaring them: with no coordinates each value is used only for ordering and as
// the barcode birth/death, so squaring would only round-trip through sqrt at output and lose precision.
template <class T>
class Matrix_geometry {
 public:
  using Filtration_value = T;
  Matrix_geometry(std::vector<T> distances, std::size_t n) : dist_(std::move(distances)), n_(n) {}
  [[nodiscard]] std::size_t size() const { return n_; }
  // Ordering scale: the raw distance; comparisons and the barcode use it directly.
  [[nodiscard]] T dist(std::size_t i, std::size_t j) const { return dist_[(i * n_) + j]; }
  // The value already is the distance, so mapping it back to an output distance is the identity.
  [[nodiscard]] static T to_distance(T distance) { return distance; }

  // k nearest points to i (including i itself at distance 0), ascending by distance then index.
  [[nodiscard]] std::vector<std::size_t> nearest(std::size_t i, std::size_t k) const {
    const T* row = &dist_[i * n_];
    k = std::min(k, n_);
    // Sort indices directly, with the matrix row as the distance lookup.
    std::vector<std::size_t> result(n_);
    std::iota(result.begin(), result.end(), std::size_t{0});
    std::partial_sort(result.begin(), result.begin() + static_cast<std::ptrdiff_t>(k), result.end(),
                      [row](std::size_t x, std::size_t y) { return row[x] != row[y] ? row[x] < row[y] : x < y; });
    result.resize(k);
    return result;
  }

  // The k nearest points to i restricted to index > i, ascending by distance (ties by index).
  [[nodiscard]] std::vector<std::size_t> nearest_neighbors_above(std::size_t i, std::size_t k) const {
    std::vector<std::size_t> result = nearest(i, k);
    keep_above(i, result);
    return result;
  }

  // Indices > i, ascending by distance from i (ties by index).
  [[nodiscard]] std::vector<std::size_t> neighbors_above(std::size_t i) const {
    const T* row = &dist_[i * n_];
    std::vector<std::size_t> result(n_ - i - 1);
    std::iota(result.begin(), result.end(), i + 1);
    std::sort(result.begin(), result.end(),
              [row](std::size_t x, std::size_t y) { return row[x] != row[y] ? row[x] < row[y] : x < y; });
    return result;
  }

  // Early-stop target: the exact RNG cycle rank.
  [[nodiscard]] std::size_t rng_early_stop_target() const { return rng_cycle_rank_matrix(*this); }
  [[nodiscard]] Lune_result<T> lune(const Batch_edge<T>& e, const Edge_map<std::size_t, std::size_t>& one_simp_to_idx,
                                    std::size_t n) const {
    return Lune_builder<T>(e, one_simp_to_idx, n).build_matrix(*this);
  }

 private:
  std::vector<T> dist_;
  std::size_t n_;
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_GEOMETRY_H_
