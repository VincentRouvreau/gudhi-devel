/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief The concept Geometry describes the metric interface the Reduced_rips core
 * (@ref Gudhi::reduced_rips::Persistence_engine) runs against. Two models are included with the module:
 * @ref Gudhi::reduced_rips::Euclidean_geometry (coordinates and a kd-tree, working in squared distances) and
 * @ref Gudhi::reduced_rips::Matrix_geometry (a bare symmetric distance matrix).
 *
 * The algorithm relies only on the *ordering* of the values returned by `dist`, so a model is free to work in
 * whatever scale it computes most cheaply and exactly. `to_distance` maps that scale back to a true distance
 * for the output barcode.
 */
struct Geometry {
  /** \brief Arithmetic type of the filtration / barcode values. Must be comparable with <. */
  typedef unspecified Filtration_value;

  /** \brief Returns the number of points. */
  std::size_t size();

  /** \brief Returns the ordering scale between points i and j: any value monotone in the true distance (for
   * instance the squared distance). Only its ordering is used by the reduction. */
  Filtration_value dist(std::size_t i, std::size_t j);

  /** \brief Maps a value on the ordering scale (as returned by `dist`) back to a true distance, for the output
   * barcode. May be a static member. */
  Filtration_value to_distance(Filtration_value ordering_value);

  /** \brief Returns approximately the k nearest points to i whose index is > i, ascending by distance (ties
   * broken by ascending index). May return fewer than k points (or none), or more when distances tie; the
   * list must be a prefix of the `neighbors_above` ordering, as the reduction resumes positionally in that
   * list after a refresh. */
  std::vector<std::size_t> nearest_neighbors_above(std::size_t i, std::size_t k);

  /** \brief Returns all points with index > i, ascending by distance from i (ties by ascending index). Used as
   * the exhaustive fallback when `nearest_neighbors_above` returns nothing above i. */
  std::vector<std::size_t> neighbors_above(std::size_t i);

  /** \brief Returns the relative-neighborhood-graph cycle rank: the number of finite degree-1 bars. This serves
   * as the early-stop target. */
  std::size_t rng_early_stop_target();

  /** \brief Returns the 2-simplices the candidate edge `e` contributes, as a `Lune_result`: either a single
   * apparent 2-simplex, or one boundary column per connected component of the edge's lune. `one_simp_to_idx`
   * maps a packed edge to its 1-simplex id, and `n` is the point count. */
  Lune_result<Filtration_value> lune(const Batch_edge<Filtration_value>& e,
                                     const Edge_map<std::size_t, std::size_t>& one_simp_to_idx, std::size_t n);
};
