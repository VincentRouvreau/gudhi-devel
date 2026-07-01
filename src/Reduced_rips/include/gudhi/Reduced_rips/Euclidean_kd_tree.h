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
 * @file Euclidean_kd_tree.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief Spatial neighbor search over a Euclidean point cloud: a CGAL kd-tree in low ambient dimension, a
 * flat brute-force scan otherwise.
 */

#ifndef REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_
#define REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <utility>
#include <vector>

#include <CGAL/Epick_d.h>

#include <gudhi/Kd_tree_search.h>
#include <gudhi/Reduced_rips/Helpers.h>

namespace Gudhi {

namespace reduced_rips {

// Spatial search over the cloud. In dimension <= 3 it uses GUDHI's Kd_tree_search. In dimension >= 4,
// it falls back to a flat brute-force scan over the Cloud. Every distance comparison uses l2_dist_2.
class Euclidean_kd_tree {
  using Cloud = detail::Cloud;

 public:
  Euclidean_kd_tree(const Cloud& pm, bool brute_force) : pm_(pm) {
    if (brute_force) return;  // leave tree_ disengaged; queries then fall back to a brute-force scan
    kd_points_.reserve(pm.n);
    for (std::size_t i = 0; i < pm.n; ++i) kd_points_.emplace_back(pm[i], pm[i] + pm.dim);
    // Kd_tree_search keeps a reference to kd_points_, which outlives the tree (both are members below).
    tree_.emplace(kd_points_);
  }

  // Indices of at least the k nearest points to query (more when distances tie at the k-th place),
  // ascending by squared distance with ties broken by index. This is the same (distance, index) order the
  // exhaustive per-point scans use, and the returned list is a closed initial segment of it: the engine's
  // heap frontier resumes positionally inside a refreshed full list, so this list must be an exact prefix
  // of it. A query that is itself a Cloud point is returned as the nearest (distance 0); callers filter
  // that out.
  std::vector<std::size_t> nearest_neighbors(const double* query, std::size_t k) const {
    if (!tree_) {
      std::vector<double> dist(pm_.n);
      for (std::size_t i = 0; i < pm_.n; ++i) dist[i] = detail::l2_dist_2<double>(pm_[i], query, pm_.dim);
      return detail::smallest_indices_by(0, pm_.n, std::min(k, pm_.n), [&dist](std::size_t x) { return dist[x]; });
    }
    // CGAL's k-nearest search orders equal distances arbitrarily, and a tie across the k-th distance even
    // makes the returned subset arbitrary. So use CGAL only to find the k-th distance, re-expressed in the
    // canonical l2_dist_2 metric, and rebuild the answer as *every* point within it (widened ball query,
    // then an exact <= cut): the full tie group is included and the order is deterministic.
    Kd_point center(query, query + pm_.dim);
    double d_k = 0.0;
    for (auto nb : tree_->k_nearest_neighbors(center, static_cast<unsigned int>(k), true))
      d_k = std::max(d_k, detail::l2_dist_2<double>(pm_[static_cast<std::size_t>(nb.first)], query, pm_.dim));
    std::vector<std::pair<std::size_t, double>> ball = points_in_squared_ball(query, detail::widen_radius(d_k));
    ball.erase(std::remove_if(ball.begin(), ball.end(),
                              [d_k](const std::pair<std::size_t, double>& pr) { return pr.second > d_k; }),
               ball.end());
    std::sort(ball.begin(), ball.end(),
              [](const std::pair<std::size_t, double>& x, const std::pair<std::size_t, double>& y) {
                return x.second != y.second ? x.second < y.second : x.first < y.first;
              });
    std::vector<std::size_t> result;
    result.reserve(ball.size());
    for (const auto& pr : ball) result.push_back(pr.first);
    return result;
  }

  // All points within the given squared radius of query, as (index, squared distance) pairs. The radius and the
  // returned squared distances are in the scalar type T; the CGAL tree searches in double, and each survivor's
  // distance is recomputed in T afterwards.
  template <class T>
  std::vector<std::pair<std::size_t, T>> points_in_squared_ball(const double* query, T squared_radius) const {
    std::vector<std::pair<std::size_t, T>> result;
    if (!tree_) {
      for (std::size_t i = 0; i < pm_.n; ++i) {
        T d = detail::l2_dist_2<T>(pm_[i], query, pm_.dim);
        if (d <= squared_radius) result.emplace_back(i, d);
      }
      return result;
    }
    Kd_point center(query, query + pm_.dim);
    std::vector<std::size_t> found;
    tree_->all_near_neighbors2(center, squared_radius, squared_radius, std::back_inserter(found));
    result.reserve(found.size());
    for (std::size_t idx : found) result.emplace_back(idx, detail::l2_dist_2<T>(pm_[idx], query, pm_.dim));
    return result;
  }

 private:
  using Kd_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Kd_point = Kd_kernel::Point_d;
  using Kd_tree = Gudhi::spatial_searching::Kd_tree_search<Kd_kernel, std::vector<Kd_point>>;

  Cloud pm_;  // non-owning view (pointer + sizes); pointee outlives this
  std::vector<Kd_point> kd_points_;
  std::optional<Kd_tree> tree_;  // disengaged in brute-force mode; engaged holds the kd-tree over kd_points_
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_
