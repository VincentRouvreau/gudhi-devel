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
 * @file Persistence_engine.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief The three-phase batched reduction that turns a geometry policy's candidate edges into the degree-1
 * persistence barcode.
 */

#ifndef REDUCED_RIPS_PERSISTENCE_ENGINE_H_
#define REDUCED_RIPS_PERSISTENCE_ENGINE_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <queue>
#include <utility>
#include <vector>

#include <gudhi/Reduced_rips/Helpers.h>
#include <gudhi/Reduced_rips/Lune_builder.h>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_for.h>
#endif

namespace Gudhi {

namespace reduced_rips {

// Computes the degree-1 persistence barcode for one geometry policy. Edges are processed in batches of three
// phases: Phase A (pop_batch) pops a batch in strict filtration order, assigning each its 1-simplex id and
// advancing the heap frontier one pop at a time; Phase B (evaluate_lunes) computes each edge's lune
// independently, reading only immutable state; Phase C (apply_result) reduces the resulting columns in order.
// After run(), the barcode (in the geometry's edge-length scale) and diagnostic counts are read back out.
template <class Geom>
class Persistence_engine {
 public:
  // The scalar type the geometry works in; births, deaths and the barcode are all in it.
  using FV = typename Geom::Filtration_value;

  Persistence_engine(Geom& geom, unsigned int num_neighbors)
      : geom_(&geom),
        n_(geom.size()),
        num_neighbors_(num_neighbors == 0 ? static_cast<std::size_t>(std::sqrt(double(n_))) : num_neighbors),
        total_death_(geom_->rng_early_stop_target()),
        neighbors_(n_) {}

  void run() {
    seed_heap();
    std::vector<Batch_edge<FV>> batch;
    std::vector<Lune_result<FV>> results;
    std::vector<std::size_t> symm_diff;  // reduction scratch, reused across columns to amortize its allocation
    while (death_counter_ < total_death_ && pop_batch(batch)) {
      results.resize(batch.size());
      evaluate_lunes(batch, results);
      for (std::size_t i = 0; i < batch.size() && death_counter_ < total_death_; ++i) {
        committed_edges_ += 1;
        apply_result(results[i], symm_diff);
      }
    }
  }

  // The (birth, death) bars in the geometry's edge-length scale (squared distances under the Euclidean policy,
  // raw distances under the matrix policy), in computed order: ascending by death (each bar dies at the
  // diameter of the 2-simplex that fills its cycle, and edges are processed shortest-first, so deaths are
  // emitted non-decreasing). The caller maps each value back to a distance via the geometry's to_distance.
  const std::vector<std::pair<FV, FV>>& barcode() const { return barcode_; }
  [[nodiscard]] std::size_t committed_edges() const { return committed_edges_; }
  [[nodiscard]] std::size_t columns_formed() const { return column_counter_; }
  [[nodiscard]] std::size_t deaths() const { return death_counter_; }

 private:
  // Min-heap ordering for the candidate edges: shortest length first, ties broken by the (a, b) index
  // pair. (While queued, Batch_edge::id carries the neighbor-list frontier position)
  struct Heap_compare {
    bool operator()(const Batch_edge<FV>& x, const Batch_edge<FV>& y) const {
      if (x.r != y.r) return x.r > y.r;
      if (x.a != y.a) return x.a > y.a;
      return x.b > y.b;
    }
  };

  // Min-heap of candidate edges, processed shortest-first.
  using Edge_heap = std::priority_queue<Batch_edge<FV>, std::vector<Batch_edge<FV>>, Heap_compare>;

  template <typename T>
  using Edge_map = detail::Edge_map<std::size_t, T>;

  // Seed the heap with, per point i, the shortest edge to a higher-indexed neighbor.
  void seed_heap() {
    for (std::size_t i = 0; i < n_; ++i) {
      neighbors_[i] = geom_->nearest_neighbors_above(i, num_neighbors_);
      // The k-nearest query may return no higher-indexed point; fall back to the full above-list.
      if (neighbors_[i].empty() && i != n_ - 1) neighbors_[i] = geom_->neighbors_above(i);
      if (!neighbors_[i].empty()) heap_.push({i, neighbors_[i][0], geom_->dist(i, neighbors_[i][0]), 0});
    }
  }

  // Phase A: pop up to batch_cap_ edges in filtration order, assign each its 1-simplex id, and advance the
  // heap frontier behind each pop. Returns false once the heap is exhausted and nothing was popped.
  bool pop_batch(std::vector<Batch_edge<FV>>& batch) {
    batch.clear();
    while (batch.size() < batch_cap_ && !heap_.empty()) {
      Batch_edge<FV> edge = heap_.top();  // edge.id currently holds the neighbor-list frontier position t
      heap_.pop();

      // Advance a's heap frontier first, while edge.id still carries the frontier position t.
      advance_frontier(edge.a, edge.id);

      // Reuse that same slot for the 1-simplex id (== the number assigned so far), then record the edge.
      edge.id = birth_by_id_.size();
      one_simp_to_idx_[detail::pack_edge(edge.a, edge.b, n_)] = edge.id;  // a < b always (stored neighbors are > a)
      birth_by_id_.push_back(edge.r);  // r == dist(a, b) by construction (a < b); cached for the barcode
      batch.push_back(edge);           // Batch_edge is trivially copyable; no move
    }
    return !batch.empty();
  }

  // Push a's next-shortest higher-indexed edge back onto the heap, refreshing a's neighbor list if exhausted.
  void advance_frontier(std::size_t a, std::size_t t) {
    if (t + 2 > neighbors_[a].size()) {
      neighbors_[a] = geom_->neighbors_above(a);
      if (neighbors_[a].size() < t + 2) {
        // Vertex a has no higher-indexed neighbor left, so no heap entry sourced at a remains and neighbors[a]
        // is never read again.
        neighbors_[a].clear();
        return;
      }
    }
    std::size_t b_next = neighbors_[a][t + 1];
    heap_.push({a, b_next, geom_->dist(a, b_next), t + 1});
  }

  // Phase B: compute each edge's lune independently (in parallel under TBB), reading only immutable state.
  void evaluate_lunes(const std::vector<Batch_edge<FV>>& batch, std::vector<Lune_result<FV>>& results) {
    auto eval_one = [&](std::size_t i) { results[i] = geom_->lune(batch[i], one_simp_to_idx_, n_); };
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(std::size_t{0}, batch.size(), eval_one);
#else
    for (std::size_t i = 0; i < batch.size(); ++i) eval_one(i);
#endif
  }

  // Phase C, one edge: file its boundary column(s) into the reduced complex and record any persistent pair.
  void apply_result(Lune_result<FV>& res, std::vector<std::size_t>& symm_diff) {
    if (res.is_apparent()) {
      // A single column is an apparent 2-simplex: file its boundary directly under its own pivot (the longest
      // edge, == the current candidate edge), with no reduction and no recorded bar (birth == death). The
      // column is the 3 ascending edge ids {v0, v1, v2} with the pivot v2 == the map key.
      column_counter_ += 1;
      const std::vector<std::size_t>& v = res.cols[0];
      root_apparent_[v[2]] = {v[0], v[1]};
      return;
    }
    // Zero columns: the lune contributed nothing. More than one: reduce one column per lune component.
    for (std::size_t c = 0; c < res.cols.size(); ++c) {
      column_counter_ += 1;
      std::vector<std::size_t>& column = res.cols[c];
      if (!reduce_column(column, symm_diff)) {
        std::size_t pivot = column.back();
        root_persistent_[pivot] = std::move(column);
        // (pivot, column) is a persistent pair: birth = pivot-edge length, death = 2-simplex diameter (both in
        // the geometry's edge-length scale). A zero-length bar is dropped.
        FV birth = birth_by_id_[pivot];
        FV death = res.deaths[c];
        if (birth != death) {
          barcode_.emplace_back(birth, death);
          death_counter_ += 1;
        }
      }
    }
  }

  // Column reduction over Z/2: while the pivot collides with a known column, add it and repick the pivot,
  // until the pivot is free or the column empties. A single hashed find() locates the colliding column
  // (apparent first, then persistent); the symmetric difference of two ascending columns is itself ascending,
  // so the new pivot is just its back element. Returns true if the column reduced to empty. `symm_diff` is
  // reusable scratch.
  bool reduce_column(std::vector<std::size_t>& column, std::vector<std::size_t>& symm_diff) {
    // Add the colliding column [first, last) into `column` over Z/2 (using symmetric difference), reusing symm_diff
    // as scratch. The symmetric difference of two ascending columns is itself ascending, so the new pivot is
    // just the back element after the swap.
    auto xor_in = [&](auto first, auto last) {
      symm_diff.clear();
      std::set_symmetric_difference(first, last, column.begin(), column.end(), std::back_inserter(symm_diff));
      column.swap(symm_diff);
    };
    // Reduce while the current pivot (column.back()) collides with a known column; stop when the column
    // empties (killed) or its pivot is free (reduced). A single hashed find() locates the colliding column,
    // apparent first then persistent.
    while (!column.empty()) {
      std::size_t pivot = column.back();
      auto it_a = root_apparent_.find(pivot);
      if (it_a != root_apparent_.end()) {
        // An apparent entry stores only its two lower edges; rebuild the full ascending 3-edge boundary
        // {e0, e1, pivot} (pivot == the key) and add it in.
        const std::array<std::size_t, 3> apparent_col{it_a->second[0], it_a->second[1], pivot};
        xor_in(apparent_col.begin(), apparent_col.end());
      } else {
        auto it_p = root_persistent_.find(pivot);
        if (it_p == root_persistent_.end()) return false;  // pivot free: the column is reduced
        xor_in(it_p->second.begin(), it_p->second.end());
      }
    }
    return true;  // column emptied: killed
  }

  Geom* geom_;  // non-owning, never null; the geometry policy outlives this engine
  std::size_t n_;
  std::size_t num_neighbors_;
  std::size_t total_death_ = 0;
  static constexpr std::size_t batch_cap_ = 4096;

  Edge_heap heap_;
  std::vector<std::vector<std::size_t>> neighbors_;     // neighbors_[i] = nearest indices > i, ascending
  Edge_map<std::size_t> one_simp_to_idx_;               // packed edge -> 1-simplex id
  std::vector<FV> birth_by_id_;                         // id -> birth length (squared under the Euclidean policy)
  Edge_map<std::array<std::size_t, 2>> root_apparent_;  // pivot -> its 2 lower edges
  Edge_map<std::vector<std::size_t>> root_persistent_;  // pivot -> reduced column
  std::vector<std::pair<FV, FV>> barcode_;              // (birth, death) bars in the geometry's scale

  std::size_t committed_edges_ = 0;  // edges actually reduced (diagnostic)
  std::size_t column_counter_ = 0;   // 2-simplex columns formed (diagnostic)
  std::size_t death_counter_ = 0;    // recorded persistent pairs
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_PERSISTENCE_ENGINE_H_
