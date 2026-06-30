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

#ifndef REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_
#define REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include <CGAL/Epick_d.h>

#include <gudhi/Kd_tree_search.h>
#include <gudhi/Reduced_rips/Cloud.h>

namespace Gudhi {

namespace reduced_rips {

using detail::Cloud;
using detail::l2_dist_2;

// Spatial search over the cloud. In dimension <= 3 it uses GUDHI's Kd_tree_search. In dimension >= 4,
// it falls back to a flat brute-force scan over the Cloud. Every distance comparison uses l2_dist_2.
class Euclidean_kd_tree {
 public:
  Euclidean_kd_tree(const Cloud& pm, bool brute_force) : pm_(pm), brute_force_(brute_force) {
    if (brute_force_) return;
    kd_points_.reserve(pm.n);
    for (std::size_t i = 0; i < pm.n; ++i) kd_points_.emplace_back(pm[i], pm[i] + pm.dim);
    // Kd_tree_search keeps a reference to kd_points_, which outlives the tree (both are members below).
    tree_ = std::make_unique<Kd_tree>(kd_points_);
  }

  // Indices of the k nearest points to query, ascending by squared distance (ties by index). A query that
  // is itself a Cloud point is returned as the nearest (distance 0); callers filter that out.
  std::vector<std::size_t> nearest_neighbors(const double* query, std::size_t k) const {
    if (brute_force_) {
      std::vector<double> dist(pm_.n);
      for (std::size_t i = 0; i < pm_.n; ++i) dist[i] = l2_dist_2(pm_[i], query, pm_.dim);
      k = std::min(k, pm_.n);
      std::vector<std::size_t> result(pm_.n);
      std::iota(result.begin(), result.end(), std::size_t{0});
      std::partial_sort(
          result.begin(), result.begin() + static_cast<std::ptrdiff_t>(k), result.end(),
          [&dist](std::size_t x, std::size_t y) { return dist[x] != dist[y] ? dist[x] < dist[y] : x < y; });
      result.resize(k);
      return result;
    }
    Kd_point center(query, query + pm_.dim);
    std::vector<std::size_t> result;
    for (auto nb : tree_->k_nearest_neighbors(center, static_cast<unsigned int>(k), true))
      result.push_back(static_cast<std::size_t>(nb.first));
    return result;
  }

  // All points within the given squared radius of query, as (index, squared distance) pairs.
  std::vector<std::pair<std::size_t, double>> points_in_squared_ball(const double* query, double squared_radius) const {
    std::vector<std::pair<std::size_t, double>> result;
    if (brute_force_) {
      for (std::size_t i = 0; i < pm_.n; ++i) {
        double d = l2_dist_2(pm_[i], query, pm_.dim);
        if (d <= squared_radius) result.emplace_back(i, d);
      }
      return result;
    }
    Kd_point center(query, query + pm_.dim);
    std::vector<std::size_t> found;
    tree_->all_near_neighbors2(center, squared_radius, squared_radius, std::back_inserter(found));
    result.reserve(found.size());
    for (std::size_t idx : found) result.emplace_back(idx, l2_dist_2(pm_[idx], query, pm_.dim));
    return result;
  }

 private:
  using Kd_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Kd_point = Kd_kernel::Point_d;
  using Kd_tree = Gudhi::spatial_searching::Kd_tree_search<Kd_kernel, std::vector<Kd_point>>;

  Cloud pm_;  // non-owning view (pointer + sizes); pointee outlives this
  bool brute_force_;
  std::vector<Kd_point> kd_points_;
  std::unique_ptr<Kd_tree> tree_;
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_
