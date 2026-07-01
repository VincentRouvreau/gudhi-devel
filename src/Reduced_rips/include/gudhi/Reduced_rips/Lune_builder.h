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

#ifndef REDUCED_RIPS_LUNE_BUILDER_H_
#define REDUCED_RIPS_LUNE_BUILDER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include <boost/pending/disjoint_sets.hpp>

#include <gudhi/Reduced_rips/Cloud.h>
#include <gudhi/Reduced_rips/Euclidean_kd_tree.h>

namespace Gudhi {

namespace reduced_rips {

using detail::Cloud;
using detail::Edge_map;
using detail::epsilon;
using detail::in_lune;
using detail::l2_dist_2;
using detail::lens_ball_factor;
using detail::pack_edge;

// ---- Per-edge "lune" computation --------------------------------------------------------------------------
// For a candidate edge (a,b) of squared length r, this determines the 2-simplices it contributes: either one
// apparent 2-simplex (whose boundary column is returned ready to file under its own pivot), or one boundary
// column per connected component of the lune (returned to be reduced).
// The same struct serves as the engine's heap entry and as the lune input. While an edge sits in the heap,
// `id` transiently holds the neighbor-list frontier position `t`. When the edge is popped into a batch that
// slot is overwritten with the assigned 1-simplex id (see Persistence_engine::pop_batch), so by the time the
// lune sees it `id` is the 1-simplex id.
struct Batch_edge {
  std::size_t a, b;  // endpoints, a < b
  double r;          // squared length of (a,b)
  std::size_t id;    // 1-simplex id (heap frontier position t while queued); also the eye-sampling RNG seed
};

// The 2-simplices a candidate edge contributes, as parallel boundary columns and their squared death values.
// The number of columns tells the engine how to file them, so no separate flags are needed:
//   - 0 columns: the lune contributes no 2-simplex (empty);
//   - 1 column:  a single apparent 2-simplex, filed directly under its own pivot (no reduction, no bar);
//   - >1 column: one boundary column per lune component, each reduced.
struct Lune_result {
  std::vector<std::vector<std::size_t>> cols;  // boundary column(s), each an ascending list of edge ids
  std::vector<double> deaths;                  // squared 2-simplex diameter for each column in `cols`

  Lune_result() = default;  // empty: no 2-simplex
  // A single 2-simplex from its ascending 3-edge boundary column and squared diameter.
  Lune_result(std::vector<std::size_t> column, double death) {
    cols.push_back(std::move(column));
    deaths.push_back(death);
  }

  // The candidate edge contributed no 2-simplex.
  [[nodiscard]] bool empty() const { return cols.empty(); }
  // Exactly one boundary column: a single apparent 2-simplex, filed under its own pivot with no recorded bar.
  [[nodiscard]] bool is_apparent() const { return cols.size() == 1; }
};

// Turns the lune points of one candidate edge into a Lune_result. The geometry-specific front-ends
// (build_euclidean / build_matrix) gather the lune points and, for the Euclidean geometry, may apply the
// lens/eye certificates of the reference paper. The work that does not depend on coordinates (the connected
// components of the lune points thresholded at r, and the boundary columns they yield) is shared by both.
class Lune_builder {
 public:
  Lune_builder(const Batch_edge& e, const Edge_map<std::size_t, std::size_t>& one_simp_to_idx, std::size_t n)
      : a_(e.a), b_(e.b), r_(e.r), id_(e.id), one_simp_to_idx_(&one_simp_to_idx), n_(n) {}

  // Euclidean front-end: gather the lune points from a midpoint ball query, apply the lens-ball fast path and
  // the wide-angle ("eye") single-component certificate, then defer the component analysis to the shared code.
  [[nodiscard]] Lune_result build_euclidean(const Cloud& pm, const Euclidean_kd_tree& kd_tree) const {
    const std::size_t a = a_, b = b_, dim = pm.dim;
    const double r = r_;

    // A zero-length edge (coincident endpoints) bounds no 2-simplex of positive diameter, so its lune is empty.
    if (r == 0.0) return Lune_result{};

    // Midpoint of (a,b): center of the candidate ball query.
    thread_local std::vector<double> mid_point;
    mid_point.resize(dim);
    for (std::size_t i = 0; i < dim; ++i) mid_point[i] = (pm[a][i] + pm[b][i]) / 2.0;
    std::vector<std::pair<std::size_t, double>> ball_mid =
        kd_tree.points_in_squared_ball(mid_point.data(), (0.75 * r) + epsilon);

    // Fast path: a candidate inside the (inscribed) lens ball certifies the lune is a single component.
    auto lens = std::find_if(ball_mid.begin(), ball_mid.end(), [&](const std::pair<std::size_t, double>& pr) {
      return pr.second <= lens_ball_factor * r && both_edges_id(pr.first);
    });
    if (lens != ball_mid.end()) return single(lens->first);

    // Slower path: prune the midpoint-ball candidates to the lune, then sort the survivors.
    std::vector<std::size_t> r_ab;
    r_ab.reserve(ball_mid.size());
    for (const auto& pr : ball_mid) {
      std::size_t k = pr.first;
      if (k == a) continue;
      double dist_ka_sq = l2_dist_2(pm[a], pm[k], dim);
      double dist_kb_sq = l2_dist_2(pm[b], pm[k], dim);
      if (in_lune(dist_ka_sq, dist_kb_sq, r, a, b, k) && both_edges_id(k)) r_ab.push_back(k);
    }
    std::sort(r_ab.begin(), r_ab.end());
    if (r_ab.empty()) return Lune_result{};

    // Heuristic (only worth sampling with more than two lune points): a wide-angle ("eye") point certifies a
    // single component without union-find. The seed is the edge id (id_), so the sampling is deterministic and
    // each parallel worker carries its own; a wrong guess only forgoes the shortcut (the exact union-find still
    // runs below) while emitting the same simplex, so the barcode does not depend on the random number generator.
    bool single_component_hint = false;
    if (r_ab.size() > 2) {
      std::mt19937 rng(id_);
      auto n_check = static_cast<std::size_t>(std::sqrt(double(r_ab.size())));
      for (std::size_t j = 0; j < n_check && !single_component_hint; ++j) {
        std::size_t temp_idx = r_ab[rng() % r_ab.size()];
        const double *pa = pm[a], *pb = pm[b], *pt = pm[temp_idx];
        double dot = 0.0, uu = 0.0, vv = 0.0;
        for (std::size_t k = 0; k < dim; ++k) {
          double u = pa[k] - pt[k], v = pb[k] - pt[k];
          dot += u * v;
          uu += u * u;
          vv += v * v;
        }
        // A single wide-angle sample (angle > 5*pi/6) certifies one component; stop sampling once one is found.
        single_component_hint = dot < 0.0 && 4.0 * dot * dot > 3.0 * uu * vv;
      }
    }

    return from_lune_points(
        r_ab, [&pm, dim](std::size_t i, std::size_t j) { return l2_dist_2(pm[i], pm[j], dim); }, single_component_hint);
  }

  // Matrix front-end: with no coordinates there is no midpoint, hence none of the Euclidean shortcuts.
  // Candidates are gathered from the row of endpoint a (the lune is contained in the closed ball of radius r
  // about a), filtered exactly by in_lune, and the components are always found by the exact union-find.
  template <class Geom>
  [[nodiscard]] Lune_result build_matrix(const Geom& g) const {
    const std::size_t a = a_, b = b_, np = g.size();
    const double r = r_;

    // A zero-length edge (coincident endpoints) bounds no 2-simplex of positive diameter, so its lune is empty.
    if (r == 0.0) return Lune_result{};

    // Lune points of (a,b): scan a's row (lune is a subset of the closed ball of radius r about a), keeping
    // those that pass the exact in_lune test against b. k increases, so r_ab is already ascending in index.
    std::vector<std::size_t> r_ab;
    for (std::size_t k = 0; k < np; ++k) {
      if (k == a || k == b) continue;
      double dist_ka = g.dist(a, k);
      if (dist_ka > r) continue;
      double dist_kb = g.dist(b, k);
      if (in_lune(dist_ka, dist_kb, r, a, b, k)) r_ab.push_back(k);
    }
    if (r_ab.empty()) return Lune_result{};

    return from_lune_points(
        r_ab, [&g](std::size_t i, std::size_t j) { return g.dist(i, j); }, /*single_component_hint=*/false);
  }

 private:
  [[nodiscard]] bool both_edges_id(std::size_t c) const {
    return one_simp_to_idx_->contains(pack_edge(std::min(a_, c), std::max(a_, c), n_)) &&
           one_simp_to_idx_->contains(pack_edge(std::min(b_, c), std::max(b_, c), n_));
  }

  // A single 2-simplex (a, b, c): its boundary column paired with the squared diameter r (c lies in the
  // closed lune, so (a,b) is the longest edge and r is the death).
  [[nodiscard]] Lune_result single(std::size_t c) const { return {column_of(c), r_}; }

  // Reduce the lune points r_ab (ascending, non-empty) to the final Lune_result. `dist` is a squared-distance
  // callable (i, j) -> double; `single_component_hint` lets the Euclidean eye certificate skip the union-find
  // when it has already certified that the lune points form a single component.
  template <class Dist>
  [[nodiscard]] Lune_result from_lune_points(const std::vector<std::size_t>& r_ab, Dist dist,
                                             bool single_component_hint) const {
    std::size_t n_rab = r_ab.size();
    if (n_rab == 1 || single_component_hint) return single(r_ab[0]);

    // One representative (a global lune-point index) per connected component of the lune points thresholded
    // at r. Two points is the cheap special case; otherwise defer to the all-pairs union-find.
    std::vector<std::size_t> reps;
    if (n_rab == 2) {
      if (dist(r_ab[0], r_ab[1]) < r_) return single(r_ab[0]);  // one component
      reps = {r_ab[0], r_ab[1]};
    } else {
      reps = component_representatives(r_ab, dist);
      if (reps.size() == 1) return single(r_ab[0]);
    }

    // Multiple components: one boundary column per component, each dying at the squared diameter r (each
    // representative lies in the closed lune of (a,b)).
    Lune_result res;
    res.cols.reserve(reps.size());
    res.deaths.reserve(reps.size());
    for (std::size_t rep : reps) {
      res.cols.push_back(column_of(rep));
      res.deaths.push_back(r_);
    }
    return res;
  }
  // All-pairs union-find over the local positions 0..|r_ab| of the lune points: merge two positions whose
  // squared distance is below r, and return each position's component root (a local position). Global indices
  // can be huge, so the disjoint-set arrays are sized by the small lune rather than by the point cloud.
  template <class Dist>
  [[nodiscard]] std::vector<std::size_t> component_roots(const std::vector<std::size_t>& r_ab, Dist dist) const {
    std::size_t n_rab = r_ab.size();
    std::vector<std::size_t> rank(n_rab), parent(n_rab);
    boost::disjoint_sets<std::size_t*, std::size_t*> ds(rank.data(), parent.data());
    for (std::size_t q = 0; q < n_rab; ++q) ds.make_set(q);
    // Count only real merges so we can stop once all points coalesce into one component (after n_rab-1 merges).
    std::size_t remaining = n_rab;
    for (std::size_t k = 0; k < n_rab && remaining > 1; ++k)
      for (std::size_t l = k + 1; l < n_rab && remaining > 1; ++l)
        if (dist(r_ab[k], r_ab[l]) < r_) {
          std::size_t rk = ds.find_set(k), rl = ds.find_set(l);
          if (rk != rl) {
            ds.link(rk, rl);
            --remaining;
          }
        }
    std::vector<std::size_t> roots(n_rab);
    for (std::size_t z = 0; z < n_rab; ++z) roots[z] = ds.find_set(z);
    return roots;
  }

  // One representative lune-point index per connected component of r_ab thresholded at r: run the union-find,
  // dedupe to one root per component, then map each local root back to its global index. The degree-1 barcode
  // is invariant to which point represents a component, so any spanning subgraph gives the same answer; for the
  // small lune sets here all-pairs is far cheaper than a Delaunay triangulation.
  template <class Dist>
  [[nodiscard]] std::vector<std::size_t> component_representatives(const std::vector<std::size_t>& r_ab,
                                                                   Dist dist) const {
    std::vector<std::size_t> reps = component_roots(r_ab, dist);
    std::sort(reps.begin(), reps.end());
    reps.erase(std::unique(reps.begin(), reps.end()), reps.end());
    for (std::size_t& rep : reps) rep = r_ab[rep];
    return reps;
  }

  // Boundary column of the 2-simplex (a, b, c): the ids of its three edges, sorted ascending so the pivot is
  // the last element. a_ < b_ by construction, so the three edges are (a_,b_) and the two (min,max) pairs with
  // c, formed without sorting the vertices. The three ids are distinct, so a fixed min/max network orders them
  // branch-free (no std::sort) -- and without assuming which edge is longest, so length ties are harmless.
  [[nodiscard]] std::vector<std::size_t> column_of(std::size_t c) const {
    const std::size_t id_ab = one_simp_to_idx_->at(pack_edge(a_, b_, n_));
    const std::size_t id_ac = one_simp_to_idx_->at(pack_edge(std::min(a_, c), std::max(a_, c), n_));
    const std::size_t id_bc = one_simp_to_idx_->at(pack_edge(std::min(b_, c), std::max(b_, c), n_));
    const std::size_t lo = std::min(id_ab, id_ac), hi = std::max(id_ab, id_ac);
    const std::size_t top = std::max(hi, id_bc), mid_hi = std::min(hi, id_bc);
    return {std::min(lo, mid_hi), std::max(lo, mid_hi), top};
  }

  std::size_t a_, b_;
  double r_;
  std::size_t id_;  // 1-simplex id of (a,b); seeds the eye-sampling RNG
  const Edge_map<std::size_t, std::size_t>* one_simp_to_idx_;  // non-owning, never null; outlives this
  std::size_t n_;
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_LUNE_BUILDER_H_
