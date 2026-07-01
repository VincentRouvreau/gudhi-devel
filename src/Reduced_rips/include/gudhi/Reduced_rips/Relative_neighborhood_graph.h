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

#ifndef REDUCED_RIPS_RELATIVE_NEIGHBORHOOD_GRAPH_H_
#define REDUCED_RIPS_RELATIVE_NEIGHBORHOOD_GRAPH_H_

#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <gudhi/Reduced_rips/Cloud.h>
#include <gudhi/Reduced_rips/Euclidean_kd_tree.h>

namespace Gudhi {

namespace reduced_rips {

using detail::Cloud;
using detail::epsilon;
using detail::l2_dist_2;

// ---- Relative neighborhood graph (RNG) cycle rank -----------------------------------------------------
// total_death = the RNG cycle rank |E| - |V| + 1 (the RNG is connected): the number of finite H1 bars, which
// drives the reduction's early-stop. Each function below returns this rank for one geometry. |V| is n, except
// the 2D/3D Delaunay can merge coincident input points into one vertex (see rng_cycle_rank_delaunay).
//
// |E| is the *true* RNG = the open lune: edge (a,b) is in the RNG iff no point lies *strictly* inside its lune
// (strict `<` on both endpoint distances); a point exactly on the boundary does NOT remove the edge. This is
// deliberately different from the lune-occupancy test inside the reduction (in_lune, in Lune_builder), whose
// lexical tie-break assigns boundary points to 2-simplices for the homology algorithm and must NOT decide RNG
// membership -- using in_lune here would over-eliminate boundary edges and undercount. Do not "unify" them.
//
// The RNG is computed with the fastest strategy for the ambient dimension (Delaunay-based in 2D/3D, direct
// otherwise).

// Deduplicated edges of the d-dimensional Delaunay triangulation, as pairs (first < second). 2D reads
// them from triangulation faces, 3D from facets. Both are an Urquhart superset of the RNG.
inline std::vector<std::pair<std::size_t, std::size_t>> delaunay_edges_2d(const Cloud& pm) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb>;
  using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>;
  using Point = K::Point_2;
  std::vector<std::pair<Point, std::size_t>> pts;
  pts.reserve(pm.size());
  for (std::size_t k = 0; k < pm.size(); ++k) pts.emplace_back(Point(pm[k][0], pm[k][1]), k);
  Delaunay t(pts.begin(), pts.end());
  std::unordered_set<std::pair<std::size_t, std::size_t>, boost::hash<std::pair<std::size_t, std::size_t>>> edges;
  for (auto it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it) {
    std::size_t i0 = it->vertex(0)->info(), i1 = it->vertex(1)->info(), i2 = it->vertex(2)->info();
    edges.insert(std::minmax(i0, i1));
    edges.insert(std::minmax(i0, i2));
    edges.insert(std::minmax(i1, i2));
  }
  return {edges.begin(), edges.end()};
}

inline std::vector<std::pair<std::size_t, std::size_t>> delaunay_edges_3d(const Cloud& pm) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_with_info_3<std::size_t, K>;
  using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
  using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
  using Point = Delaunay::Point;
  std::vector<std::pair<Point, std::size_t>> pts;
  pts.reserve(pm.size());
  for (std::size_t k = 0; k < pm.size(); ++k) pts.emplace_back(Point(pm[k][0], pm[k][1], pm[k][2]), k);
  Delaunay t(pts.begin(), pts.end());
  // The finite edge iterator visits each Delaunay edge exactly once, so we can emit
  // edges directly instead of deduplicating facet edges through a hash set.
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  edges.reserve(t.number_of_finite_edges());
  for (auto it = t.finite_edges_begin(); it != t.finite_edges_end(); ++it) {
    std::size_t a = it->first->vertex(it->second)->info();
    std::size_t b = it->first->vertex(it->third)->info();
    edges.emplace_back(std::minmax(a, b));
  }
  return edges;
}

// RNG cycle rank from the Delaunay edges (an Urquhart superset of the RNG): discard every edge whose open
// lune contains a point, then return |E| - |V| + 1. |V| is the number of *participating* vertices, counted
// from the surviving edges rather than taken as n, because CGAL's 2D/3D Delaunay merges coincident input
// points into one vertex; the n - |V| merged duplicates contribute only zero-persistence cycles, which the
// reduction drops, so counting them out keeps the rank exact. One ball query around an endpoint suffices,
// since each candidate is re-tested against the other endpoint.
template <typename DelaunayEdges>
std::size_t rng_cycle_rank_delaunay(const Cloud& pm, const Euclidean_kd_tree& kd_tree, DelaunayEdges delaunay_edges) {
  std::vector<std::pair<std::size_t, std::size_t>> possible_edges = delaunay_edges(pm);
  std::vector<char> seen(pm.n, 0);
  std::size_t kept = 0, vertices = 0;
  for (const auto& edge : possible_edges) {
    std::size_t a = edge.first, b = edge.second;
    double r = l2_dist_2(pm[a], pm[b], pm.dim);
    // epsilon widens only the ball-query radius (so the strict test below sees every candidate); occupancy is
    // the open lune: a point strictly inside both endpoint balls. Boundary points do not remove the edge.
    auto ball = kd_tree.points_in_squared_ball(pm[a], r + epsilon);
    bool lune_occupied = std::any_of(ball.begin(), ball.end(), [&](const std::pair<std::size_t, double>& pr) {
      std::size_t k = pr.first;
      if (k == a || k == b) return false;
      double dist_ka_sq = l2_dist_2(pm[a], pm[k], pm.dim);
      double dist_kb_sq = l2_dist_2(pm[b], pm[k], pm.dim);
      return dist_ka_sq < r && dist_kb_sq < r;
    });
    if (!lune_occupied) {
      ++kept;
      for (std::size_t v : {a, b})
        if (seen[v] == 0) {
          seen[v] = 1;
          ++vertices;
        }
    }
  }
  return kept - vertices + 1;  // |E| - |V| + 1 for the connected RNG; >= 0, and 1 if there are no edges
}

// A candidate edge as (lo, hi, squared length), ordered by length then index for the set operations. Shared
// by the direct (general-dimension) and distance-matrix RNG supergraph builds.
struct Rng_edge {
  std::size_t i, j;
  double length;
  Rng_edge(std::size_t a, std::size_t b, double len) : i(std::min(a, b)), j(std::max(a, b)), length(len) {}
  bool operator<(const Rng_edge& o) const {
    if (length != o.length) return length < o.length;
    if (i != o.i) return i < o.i;
    return j < o.j;
  }
  bool operator==(const Rng_edge& o) const { return i == o.i && j == o.j; }
};

// Phase 1 of the direct RNG construction: an O(n^2) supergraph of the RNG, returned sorted and deduped.
// `dist` is a squared-distance callable (i, j) -> double; both the coordinate and distance-matrix paths
// share this body and differ only in how they supply distances. e_i is a sorted vector, front() is the
// current minimum, and the survivors of each pruning pass are compacted in place, preserving sorted order
// with no per-iteration allocation. e_all is accumulated then deduped once at the end. An edge (i,j) can
// be emitted from both its endpoint passes, which the old std::set merged.
template <class Dist>
std::vector<Rng_edge> rng_supergraph(std::size_t n, Dist dist) {
  std::vector<Rng_edge> e_all;
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<Rng_edge> e_i;
    e_i.reserve(n - 1);
    for (std::size_t j = 0; j < n; ++j)
      if (i != j) e_i.emplace_back(i, j, dist(i, j));
    std::sort(e_i.begin(), e_i.end());
    while (!e_i.empty()) {
      Rng_edge min_edge = e_i.front();
      e_all.push_back(min_edge);
      std::size_t w = 0;
      for (std::size_t t = 1; t < e_i.size(); ++t) {
        const Rng_edge& edge = e_i[t];
        std::size_t other = (min_edge.i == edge.i || min_edge.i == edge.j) ? min_edge.j : min_edge.i;
        std::size_t far = (edge.i == min_edge.i || edge.i == min_edge.j) ? edge.j : edge.i;
        Rng_edge edge_1(other, far, dist(other, far));
        Rng_edge edge_2(i, far, dist(i, far));
        if (edge_2 < edge_1) e_i[w++] = edge;
      }
      e_i.erase(e_i.begin() + static_cast<std::ptrdiff_t>(w), e_i.end());
    }
  }
  std::sort(e_all.begin(), e_all.end());
  e_all.erase(std::unique(e_all.begin(), e_all.end()), e_all.end());
  return e_all;
}

// Direct RNG construction for general dimension. Phase 1 builds an O(n^2) RNG superset, phase 2 prunes
// edges whose lune is non-empty; returns the cycle rank |E| - n + 1. Edge lengths are squared distances
// throughout. No vertex merging here (unlike the Delaunay path), so |V| = n.
inline std::size_t rng_cycle_rank_general(const Cloud& pm, const Euclidean_kd_tree& kd_tree) {
  std::vector<Rng_edge> e_all =
      rng_supergraph(pm.n, [&pm](std::size_t i, std::size_t j) { return l2_dist_2(pm[i], pm[j], pm.dim); });

  // Phase 2: eliminate edges whose open lune contains a point (strictly inside both endpoint balls). epsilon
  // widens only the ball-query radius; boundary points do not remove the edge.
  std::size_t count = 0;
  for (const auto& edge : e_all) {
    std::size_t a = edge.i, b = edge.j;
    auto ball = kd_tree.points_in_squared_ball(pm[a], edge.length + epsilon);
    bool lune_occupied = std::any_of(ball.begin(), ball.end(), [&](const std::pair<std::size_t, double>& pr) {
      std::size_t k = pr.first;
      if (k == a || k == b) return false;
      double dist_ka_sq = l2_dist_2(pm[a], pm[k], pm.dim);
      double dist_kb_sq = l2_dist_2(pm[b], pm[k], pm.dim);
      return dist_ka_sq < edge.length && dist_kb_sq < edge.length;
    });
    if (!lune_occupied) ++count;
  }
  return count - pm.n + 1;
}

// Exact early-stop target for a bare distance matrix: the RNG cycle rank |E| - n + 1, i.e. the number of
// finite H1 bars. Mirrors rng_cycle_rank_general (O(n^2) supergraph build, then O(|E| n) open-lune pruning,
// reading every distance from the geometry). No vertex merging here, so |V| = n.
template <class Geom>
std::size_t rng_cycle_rank_matrix(const Geom& g) {
  std::size_t n = g.size();
  std::vector<Rng_edge> e_all = rng_supergraph(n, [&g](std::size_t i, std::size_t j) { return g.dist(i, j); });

  // Occupancy is the open lune: a point strictly inside both endpoint balls.
  std::size_t kept = 0;
  for (const auto& edge : e_all) {
    std::size_t a = edge.i, b = edge.j;
    bool occupied = false;
    for (std::size_t k = 0; k < n && !occupied; ++k) {
      if (k == a || k == b) continue;
      if (g.dist(a, k) < edge.length && g.dist(b, k) < edge.length) occupied = true;
    }
    if (!occupied) ++kept;
  }
  return kept - n + 1;
}

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_RELATIVE_NEIGHBORHOOD_GRAPH_H_
