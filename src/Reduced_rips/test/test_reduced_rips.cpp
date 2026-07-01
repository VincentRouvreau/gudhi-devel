/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *    - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE reduced_rips
#include <boost/test/unit_test.hpp>

#include <gudhi/Reduced_rips.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/Unitary_tests_utils.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using Reduced_rips = Gudhi::reduced_rips::Reduced_rips<>;
using Bars = std::vector<std::pair<double, double>>;
using Cloud = std::vector<std::vector<double>>;

// Standard full-complex degree-1 PH pipeline, used as ground truth (with Z/2Z coefficients, matching
// Reduced_rips). A double-precision Simplex_tree keeps the comparison tight.
using Stree = Gudhi::Simplex_tree<>;
using Filtration_value = Stree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Stree, Field_Zp>;

// ---- Point-cloud generators ----------------------------------------------------------------------------

// n points evenly sampled on the unit circle in dimension `dim` (padded with zeros). Degree-1 PH of such
// a sample has exactly one prominent class (the loop), born around the sampling spacing and dying near the
// circle diameter.
static Cloud circle(unsigned n, std::size_t dim = 2) {
  Cloud pts;
  for (unsigned i = 0; i < n; ++i) {
    std::vector<double> p(dim, 0.0);
    double theta = 2.0 * M_PI * i / n;
    p[0] = std::cos(theta);
    p[1] = std::sin(theta);
    pts.push_back(std::move(p));
  }
  return pts;
}

// Two unit circles far apart on the x-axis: two independent H1 loops.
static Cloud two_circles(unsigned n) {
  Cloud a = circle(n, 2);
  Cloud b = circle(n, 2);
  for (auto& p : b) p[0] += 10.0;
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

// A torus in R^3 (nu points around the tube times nv around it): H1 of rank 2. Small enough that the full
// Rips ground truth is cheap, large enough to exercise the 3D Delaunay relative-neighborhood-graph path.
static Cloud torus(unsigned nu, unsigned nv, double R = 2.0, double r = 0.7) {
  Cloud pts;
  for (unsigned i = 0; i < nu; ++i) {
    double u = 2.0 * M_PI * i / nu;
    for (unsigned j = 0; j < nv; ++j) {
      double v = 2.0 * M_PI * j / nv;
      pts.push_back({(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)});
    }
  }
  return pts;
}

// Deterministic uniform random cloud in the unit cube of the given dimension.
static Cloud random_cloud(unsigned n, std::size_t dim, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Cloud pts(n, std::vector<double>(dim));
  for (auto& p : pts)
    for (auto& c : p) c = unit(gen);
  return pts;
}

// Full (symmetric, n x n) Euclidean distance matrix of a point cloud.
static Cloud full_distance_matrix(const Cloud& pts) {
  std::size_t n = pts.size();
  Cloud m(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (std::size_t k = 0; k < pts[i].size(); ++k) {
        double d = pts[i][k] - pts[j][k];
        s += d * d;
      }
      m[i][j] = std::sqrt(s);
    }
  return m;
}

// ---- Barcode helpers -----------------------------------------------------------------------------------

// Finite bars of positive persistence, sorted. Reduced_rips already drops zero-length bars; applying the
// same filter to both sides makes the comparison robust to near-zero numerical noise.
static Bars finite_positive(Bars b, double eps = 1e-9) {
  Bars out;
  for (const auto& x : b)
    if (std::isfinite(x.second) && x.second - x.first > eps) out.push_back(x);
  std::sort(out.begin(), out.end());
  return out;
}

// True iff two barcodes match as multisets within tolerance. Distances squared then square-rooted along the
// matrix path can differ from the coordinate path by a few ULPs, hence the tolerance rather than equality.
static bool bars_close(Bars a, Bars b, double tol = 1e-6) {
  a = finite_positive(a);
  b = finite_positive(b);
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::abs(a[i].first - b[i].first) > tol || std::abs(a[i].second - b[i].second) > tol) return false;
  return true;
}

// Number of bars more persistent than `min_persistence` (counts prominent topological features).
static std::size_t count_prominent(const Bars& b, double min_persistence) {
  std::size_t c = 0;
  for (const auto& x : b)
    if (x.second - x.first > min_persistence) ++c;
  return c;
}

// Ground-truth degree-1 barcode of the full (un-thresholded) Vietoris-Rips filtration.
static Bars full_rips_h1(const Cloud& pts) {
  Rips_complex rips(pts, std::numeric_limits<double>::infinity(), Gudhi::Euclidean_distance());
  Stree st;
  rips.create_complex(st, 2);
  Persistent_cohomology pcoh(st);
  pcoh.init_coefficients(2);
  pcoh.compute_persistent_cohomology(0.0);
  Bars bars;
  for (const auto& bd : pcoh.intervals_in_dimension(1)) bars.emplace_back(bd.first, bd.second);
  return bars;
}

// Ground-truth degree-1 barcode of the full Vietoris-Rips filtration of a (full, symmetric) distance matrix.
// 2-simplices fill at their longest edge, matching the Reduced_rips convention; no triangle inequality needed.
static Bars full_rips_h1_from_matrix(const std::vector<std::vector<double>>& full) {
  Stree st;
  for (std::size_t i = 0; i < full.size(); ++i)
    for (std::size_t j = i + 1; j < full.size(); ++j)
      st.insert_simplex_and_subfaces({static_cast<int>(i), static_cast<int>(j)}, full[i][j]);
  st.expansion(2);
  st.make_filtration_non_decreasing();
  Persistent_cohomology pcoh(st);
  pcoh.init_coefficients(2);
  pcoh.compute_persistent_cohomology(0.0);
  Bars truth;
  for (const auto& bd : pcoh.intervals_in_dimension(1)) truth.emplace_back(bd.first, bd.second);
  return truth;
}

// ---- Correctness against the full Vietoris-Rips filtration ---------------------------------------------

BOOST_AUTO_TEST_CASE(matches_full_vietoris_rips_h1) {
  // The headline guarantee of the reference paper: the reduced filtration has the same degree-1 barcode as
  // the full Vietoris-Rips filtration. Check it exactly on a spread of inputs, exercising the 2D and 3D
  // Delaunay relative-neighborhood-graph paths and the dimension-free O(n^2) path, and both backends.
  struct Case {
    const char* name;
    Cloud pts;
  };
  std::vector<Case> cases = {
      {"circle 2D", circle(40, 2)},
      {"circle 5D", circle(36, 5)},
      {"torus 3D", torus(12, 8)},
      {"two loops 2D", two_circles(30)},
      {"random 2D", random_cloud(60, 2, 1)},
      {"random 3D", random_cloud(55, 3, 2)},
      {"random 5D", random_cloud(45, 5, 3)},
  };
  for (const auto& c : cases) {
    BOOST_TEST_CONTEXT(c.name) {
      Bars truth = full_rips_h1(c.pts);
      Bars from_points = Reduced_rips::from_points(c.pts).persistence();
      Bars from_matrix = Reduced_rips::from_distance_matrix(full_distance_matrix(c.pts)).persistence();
      BOOST_CHECK(bars_close(truth, from_points));
      BOOST_CHECK(bars_close(truth, from_matrix));
    }
  }
}

// ---- Topology recovery on hand-understood inputs -------------------------------------------------------

BOOST_AUTO_TEST_CASE(circle_has_one_dominant_loop) {
  for (std::size_t dim : {std::size_t{2}, std::size_t{5}}) {  // 2D Delaunay path and >=4D general path
    BOOST_TEST_CONTEXT("dim=" << dim) {
      Bars bc = Reduced_rips::from_points(circle(60, dim)).persistence();
      // Exactly one prominent loop, and it dies past 1.5 (near the unit circle's diameter of 2).
      BOOST_CHECK_EQUAL(count_prominent(bc, 0.5), 1u);
      auto loop = *std::max_element(bc.begin(), bc.end(),
                                    [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                                      return a.second - a.first < b.second - b.first;
                                    });
      BOOST_CHECK_GT(loop.second, 1.5);
    }
  }
}

BOOST_AUTO_TEST_CASE(two_separated_circles_give_two_loops) {
  Bars bc = Reduced_rips::from_points(two_circles(30)).persistence();
  BOOST_CHECK_EQUAL(count_prominent(bc, 0.5), 2u);
}

BOOST_AUTO_TEST_CASE(unit_square_known_bar) {
  // Four points of a unit square (sides 1, diagonals sqrt 2): a single H1 loop born when the four unit edges
  // close the cycle and dying when a diagonal fills it in.
  // Use the matrix form to pin exact distances independent of trig rounding.
  const double s2 = std::sqrt(2.0);
  std::vector<std::vector<double>> lower = {
      {},             // point 0
      {1.0},          // d(1,0)
      {s2, 1.0},      // d(2,0), d(2,1)
      {1.0, s2, 1.0}  // d(3,0), d(3,1), d(3,2)
  };
  Bars bc = Reduced_rips::from_distance_matrix(lower).persistence();
  BOOST_REQUIRE_EQUAL(bc.size(), 1U);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(bc.front().first, 1.0, 1e-6);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(bc.front().second, s2, 1e-6);
}

// ---- Distance-matrix backend ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(distance_matrix_matches_point_cloud) {
  // The distance-matrix path, fed a Euclidean cloud's own distances, must reproduce the coordinate path.
  for (unsigned dim : {2U, 3U, 5U}) {
    BOOST_TEST_CONTEXT("dim=" << dim) {
      auto pts = circle(50, dim);
      auto from_points = Reduced_rips::from_points(pts);
      Reduced_rips from_matrix = Reduced_rips::from_distance_matrix(full_distance_matrix(pts));
      BOOST_CHECK(bars_close(from_points.persistence(), from_matrix.persistence()));
    }
  }
}

BOOST_AUTO_TEST_CASE(full_and_lower_triangular_agree) {
  // The same square given as a full symmetric matrix yields the same barcode as the lower-triangular form.
  const double s2 = std::sqrt(2.0);
  std::vector<std::vector<double>> lower = {{}, {1.0}, {s2, 1.0}, {1.0, s2, 1.0}};
  std::vector<std::vector<double>> full = {
      {0.0, 1.0, s2, 1.0}, {1.0, 0.0, 1.0, s2}, {s2, 1.0, 0.0, 1.0}, {1.0, s2, 1.0, 0.0}};
  auto a = Reduced_rips::from_distance_matrix(lower).persistence();
  auto b = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK(bars_close(a, b));
}

BOOST_AUTO_TEST_CASE(non_metric_dissimilarity_is_accepted) {
  // The reduction never uses the triangle inequality, so a symmetric matrix that badly violates it is still
  // a valid input and must compute the exact degree-1 barcode of that dissimilarity's Vietoris-Rips
  // filtration. Compare against the full-complex ground truth built from the same matrix.
  std::vector<std::vector<double>> full = {
      {0.0, 5.0, 1.0, 1.0}, {5.0, 0.0, 1.0, 1.0}, {1.0, 1.0, 0.0, 5.0}, {1.0, 1.0, 5.0, 0.0}};
  Bars truth = full_rips_h1_from_matrix(full);
  Bars reduced = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK(bars_close(truth, reduced));
}

// ---- Scalar (Filtration_value) type --------------------------------------------------------------------

// Re-express a Filtration_value-typed barcode as the double-valued Bars the comparison helpers consume.
template <class FV>
static Bars to_double_bars(const std::vector<std::pair<FV, FV>>& b) {
  Bars out;
  out.reserve(b.size());
  for (const auto& x : b) out.emplace_back(static_cast<double>(x.first), static_cast<double>(x.second));
  return out;
}

// The reduction is templated on Filtration_value: distances, edge lengths, lune deaths and the barcode are all
// carried in that type. A given instantiation must recover the same degree-1 barcode as the double ground
// truth, to that type's precision, on both the point-cloud and the distance-matrix paths.
template <class FV>
static void check_scalar_type(const char* tag, double tol) {
  using RR = Gudhi::reduced_rips::Reduced_rips<FV>;
  static_assert(std::is_same_v<typename RR::Filtration_value, FV>, "Reduced_rips<FV>::Filtration_value must be FV");
  static_assert(std::is_same_v<typename RR::Persistence_interval, std::pair<FV, FV>>,
                "the barcode must be pairs of FV, with no widening to double");
  BOOST_TEST_CONTEXT(tag) {
    auto pts = circle(40, 2);
    Bars truth = full_rips_h1(pts);

    std::vector<std::pair<FV, FV>> from_points = RR::from_points(pts).persistence();  // copy out of the temporary
    std::vector<std::pair<FV, FV>> from_matrix = RR::from_distance_matrix(full_distance_matrix(pts)).persistence();

    BOOST_CHECK(bars_close(truth, to_double_bars(from_points), tol));
    BOOST_CHECK(bars_close(truth, to_double_bars(from_matrix), tol));
  }
}

BOOST_AUTO_TEST_CASE(filtration_value_scalar_types) {
  // float carries ~7 significant digits (squared then square-rooted), so it needs a looser tolerance than the
  // wider types; all three must still recover the loop.
  check_scalar_type<float>("float", 1e-3);
  check_scalar_type<double>("double", 1e-6);
  check_scalar_type<long double>("long double", 1e-9);
}

BOOST_AUTO_TEST_CASE(integer_distance_matrix) {
  // The matrix geometry carries the supplied dissimilarities verbatim -- no squaring, no square root -- so an
  // integer Filtration_value keeps the whole reduction in exact integer arithmetic. A 4-cycle with edge length
  // 2 and diagonals 3 (a scaled square) has a single H1 loop: born when the four length-2 edges close the cycle,
  // dying when a diagonal triangle fills it, giving the exact bar (2, 3).
  using RR = Gudhi::reduced_rips::Reduced_rips<int>;
  static_assert(std::is_same_v<RR::Filtration_value, int>);
  std::vector<std::vector<int>> lower = {
      {},         // point 0
      {2},        // d(1,0)
      {3, 2},     // d(2,0), d(2,1)
      {2, 3, 2},  // d(3,0), d(3,1), d(3,2)
  };
  std::vector<std::pair<int, int>> bc = RR::from_distance_matrix(lower).persistence();
  BOOST_REQUIRE_EQUAL(bc.size(), 1U);
  BOOST_CHECK_EQUAL(bc.front().first, 2);   // exact integer birth
  BOOST_CHECK_EQUAL(bc.front().second, 3);  // exact integer death
}

// ---- Points exactly on lune boundaries -----------------------------------------------------------------

BOOST_AUTO_TEST_CASE(equilateral_triangle_boundary_is_trivial) {
  // Three points at mutual distance 1: each vertex lies *exactly* on the lune boundary of the opposite edge
  // (the d(a,c) == d(b,c) == d(a,b) double-boundary case). The true relative neighborhood graph uses the open
  // lune, so a boundary point does not remove an edge: all three edges are kept, the single cycle is filled by
  // the triangle at the same filtration value, and H1 is trivial. The matrix pins the distances bit-exact, so
  // the boundary comparisons are genuine equalities rather than near-ties.
  std::vector<std::vector<double>> full = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
  Bars reduced = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK_EQUAL(finite_positive(reduced).size(), 0u);
  BOOST_CHECK(bars_close(full_rips_h1_from_matrix(full), reduced));
}

BOOST_AUTO_TEST_CASE(boundary_point_on_loop_matches_full_rips) {
  // A unit-square loop (bar (1, sqrt2)) plus a fifth point sitting *exactly* on a lune boundary of one square
  // edge: distance 1 (== the edge length) from corner 0 and strictly less (0.5) from corner 1. This exercises
  // in_lune's single-boundary branch in the reduction while the true RNG (open lune) still keeps the edge, so
  // the loop must survive. Distances are pinned exactly by the matrix; the reduced barcode must match the full
  // Vietoris-Rips ground truth on the same matrix.
  const double s2 = std::sqrt(2.0);
  std::vector<std::vector<double>> full = {{0.0, 1.0, s2, 1.0, 1.0},
                                           {1.0, 0.0, 1.0, s2, 0.5},
                                           {s2, 1.0, 0.0, 1.0, 2.0},
                                           {1.0, s2, 1.0, 0.0, 2.0},
                                           {1.0, 0.5, 2.0, 2.0, 0.0}};
  Bars truth = full_rips_h1_from_matrix(full);
  Bars reduced = Reduced_rips::from_distance_matrix(full).persistence();
  BOOST_CHECK(bars_close(truth, reduced));
  BOOST_CHECK_EQUAL(count_prominent(reduced, 0.3), count_prominent(truth, 0.3));
}

// ---- Search-strategy equivalence -----------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(search_strategies_agree) {
  // kd-tree and brute-force neighbor search are an implementation choice; the barcode must not depend on it,
  // in low ambient dimension (where automatic picks kd-tree) and high (where it picks brute-force).
  for (std::size_t dim : {std::size_t{3}, std::size_t{6}}) {
    BOOST_TEST_CONTEXT("dim=" << dim) {
      auto pts = random_cloud(50, dim, 11);
      Bars various = Reduced_rips::from_points(pts, 0, Reduced_rips::Search::automatic).persistence();
      Bars kd = Reduced_rips::from_points(pts, 0, Reduced_rips::Search::kd_tree).persistence();
      Bars brute = Reduced_rips::from_points(pts, 0, Reduced_rips::Search::brute_force).persistence();
      BOOST_CHECK(bars_close(various, kd));
      BOOST_CHECK(bars_close(various, brute));
    }
  }
}

BOOST_AUTO_TEST_CASE(initial_neighbor_budget_does_not_change_result) {
  // num_neighbors only seeds the heap; the frontier grows on demand, so the barcode is independent of it.
  auto pts = random_cloud(50, 3, 21);
  Bars budget_default = Reduced_rips::from_points(pts, 0).persistence();
  Bars budget_small = Reduced_rips::from_points(pts, 3).persistence();
  Bars budget_large = Reduced_rips::from_points(pts, 40).persistence();
  BOOST_CHECK(bars_close(budget_default, budget_small));
  BOOST_CHECK(bars_close(budget_default, budget_large));
}

// ---- Robustness, caching, diagnostics ------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(persistence_is_cached_and_deterministic) {
  auto ph1 = Reduced_rips::from_points(circle(40, 2));
  auto first = ph1.persistence();  // copy
  const auto& second = ph1.persistence();
  BOOST_CHECK(first == second);

  // A fresh instance on identical input must give an identical barcode (no global RNG state leakage).
  auto other = Reduced_rips::from_points(circle(40, 2));
  BOOST_CHECK(other.persistence() == first);
}

BOOST_AUTO_TEST_CASE(coincident_points_do_not_break_the_barcode) {
  // A zero-distance pair of distinct points (here a duplicate of point 0) once crashed the point-cloud path:
  // the degenerate r == 0 lens query emitted a self-edge 2-simplex. Coincident points carry no H1, so the
  // barcode must match the duplicate-free circle, in both the 2D Delaunay and the >=4D general RNG paths.
  for (unsigned dim : {2u, 5u}) {
    BOOST_TEST_CONTEXT("dim=" << dim) {
      auto clean = circle(40, dim);
      auto dup = clean;
      dup.push_back(clean[0]);  // distance 0 between two distinct indices

      auto from_clean = Reduced_rips::from_points(clean);
      auto from_dup = Reduced_rips::from_points(dup);
      BOOST_CHECK(bars_close(from_clean.persistence(), from_dup.persistence()));

      // The distance-matrix backend must agree on the same coincident-point input.
      Reduced_rips from_matrix = Reduced_rips::from_distance_matrix(full_distance_matrix(dup));
      BOOST_CHECK(bars_close(from_dup.persistence(), from_matrix.persistence()));
    }
  }
}

BOOST_AUTO_TEST_CASE(near_tie_distances_do_not_break_the_barcode) {
  // Regression for a floating-point hazard. A certified 2-simplex (a,b,c) records the ids of its three
  // edges, and the reduction relies on the two edges shorter than the diameter edge (a,b) having already been
  // popped and id'd. The index tie-break in the lune test enforces this for exact distances, but when two
  // edges are a near-tie in length, recomputing one of them with a different rounding (e.g. FMA contraction
  // under an optimized -march=native build) could admit a third vertex whose edge had not yet been id'd,
  // throwing std::out_of_range from the boundary lookup. A dense 3D torus has many such near-ties between
  // adjacent samples and reproduces the crash on an optimized build; the reduction must instead complete and
  // recover the torus's rank-2 H1 (two prominent loops). NOTE: this only exercises the hazard under
  // optimization with FMA contraction; an -O0 build cannot contract and so cannot reproduce it.
  Bars bc = Reduced_rips::from_points(torus(60, 40)).persistence();  // must not throw
  BOOST_CHECK_EQUAL(count_prominent(bc, 0.5), 2U);
  for (const auto& bar : bc) {
    BOOST_CHECK_GE(bar.first, 0.0);
    BOOST_CHECK_GT(bar.second, bar.first);
    BOOST_CHECK(std::isfinite(bar.second));
  }
}

BOOST_AUTO_TEST_CASE(diagnostic_counters_are_consistent) {
  auto ph1 = Reduced_rips::from_points(circle(40, 2));
  const auto& bc = ph1.persistence();
  // Recorded persistent pairs are exactly the (non-trivial) bars returned.
  BOOST_CHECK_EQUAL(ph1.num_persistence_pairs(), bc.size());
  // Every persistent pair needs a reduced 2-simplex column, and the loop produces at least one bar.
  BOOST_CHECK_GE(ph1.num_two_simplices(), ph1.num_persistence_pairs());
  BOOST_CHECK_GE(ph1.num_one_simplices(), 1u);
  BOOST_CHECK_GE(bc.size(), 1u);
}

BOOST_AUTO_TEST_CASE(off_file_point_cloud_integration) {
  // End-to-end on a real OFF file (a 3D torus of 300 points): exercises file I/O, the 3D Delaunay path at a
  // larger scale, and basic barcode well-formedness.
  Gudhi::Points_off_reader<std::vector<double>> off_reader("tore3D_300.off");
  BOOST_REQUIRE(off_reader.is_valid());
  auto ph1 = Reduced_rips::from_points(off_reader.get_point_cloud());
  const auto& bc = ph1.persistence();
  BOOST_REQUIRE(!bc.empty());
  for (const auto& bar : bc) {
    BOOST_CHECK_GE(bar.first, 0.0);
    BOOST_CHECK_GT(bar.second, bar.first);
    BOOST_CHECK(std::isfinite(bar.second));
  }
}

// Fewer than two points is not an error: the barcode is simply empty (there are no edges, hence no H1).
BOOST_AUTO_TEST_CASE(too_few_points_give_an_empty_barcode) {
  std::vector<std::vector<double>> no_points;
  BOOST_CHECK(Reduced_rips::from_points(no_points).persistence().empty());

  std::vector<std::vector<double>> one_point = {{0.0, 0.0}};
  auto from_one_point = Reduced_rips::from_points(one_point);
  BOOST_CHECK(from_one_point.persistence().empty());
  BOOST_CHECK_EQUAL(from_one_point.dimension(), 2u);

  // Zero-dimensional points carry no edges either: empty barcode, no error.
  std::vector<std::vector<double>> zero_dim = {{}, {}};
  BOOST_CHECK(Reduced_rips::from_points(zero_dim).persistence().empty());

  std::vector<std::vector<double>> empty_matrix;
  BOOST_CHECK(Reduced_rips::from_distance_matrix(empty_matrix).persistence().empty());

  std::vector<std::vector<double>> one_row = {{}};
  BOOST_CHECK(Reduced_rips::from_distance_matrix(one_row).persistence().empty());
}

// The dimension-consistency check uses GUDHI_CHECK, which throws only in debug mode (GUDHI_DEBUG), so this is
// compiled only then. A genuine dimension mismatch among points is still rejected.
#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(rejects_ragged_input) {
  std::vector<std::vector<double>> ragged = {{0.0, 0.0}, {1.0}};
  BOOST_CHECK_THROW(Reduced_rips::from_points(ragged), std::invalid_argument);
}
#endif
