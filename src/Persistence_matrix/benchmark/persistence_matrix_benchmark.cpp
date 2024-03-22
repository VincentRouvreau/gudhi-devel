/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#include <gudhi/Simplex_tree.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistence_matrix.h>
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/random_graph_generators.h>
#include <gudhi/Clock.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <random>
#include <algorithm>

std::mt19937 gen(rd());

double get_random()
{
    std::uniform_real_distribution<double> dist(0., 1.);
    return dist(gen);
}

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;

using Simplex_tree = Gudhi::Simplex_tree<>;
using Persistent_cohomology_stree = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;

using Bitmap_cubical_complex_base = Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>;
using Bitmap_cubical_complex = Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base>;
using Persistent_cohomology_cub = Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp>;

const auto Column_type = Gudhi::persistence_matrix::Column_types::INTRUSIVE_SET;

struct RU_matrix_options_vine : Gudhi::persistence_matrix::Default_options<Column_type, true> {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = true;
};

struct RU_matrix_options_repr : Gudhi::persistence_matrix::Default_options<Column_type, true> {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = true;
  static const bool is_of_boundary_type = true;
};

struct Boundary_matrix_options : Gudhi::persistence_matrix::Default_options<Column_type, true> {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = true;
};

struct Chain_matrix_options : Gudhi::persistence_matrix::Default_options<Column_type, true> {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = false;
};

struct Chain_matrix_options_vine : Gudhi::persistence_matrix::Default_options<Column_type, true> {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = false;
};

using Persistence_RU_matrix_vine_stree  = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, RU_matrix_options_vine>;
using Persistence_RU_matrix_repr_stree  = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, RU_matrix_options_repr>;
using Persistence_boundary_matrix_stree = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, Boundary_matrix_options>;
using Persistence_chain_matrix_stree    = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, Chain_matrix_options>;
using Persistence_chain_matrix_vine_stree    = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, Chain_matrix_options_vine>;

using Persistence_RU_matrix_vine_cub  = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, RU_matrix_options_vine>;
using Persistence_RU_matrix_repr_cub  = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, RU_matrix_options_repr>;
using Persistence_boundary_matrix_cub = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, Boundary_matrix_options>;
using Persistence_chain_matrix_cub    = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, Chain_matrix_options>;
using Persistence_chain_matrix_vine_cub    = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, Chain_matrix_options_vine>;

template <class Persistence, class Complex>
void bench_persistence(Complex& cpx, const std::string& msg) {
  Gudhi::Clock clock_pers(msg);
  Persistence pers(cpx);
  pers.init_coefficients(2);
  pers.compute_persistent_cohomology(0.0);
  clock_pers.end();
  std::clog << cpx.num_simplices() << " simplices - " << clock_pers;

  std::string sep("_");
  std::string ext(".pers");
  std::string pers_file_name = msg + sep + std::to_string(cpx.num_simplices())
                                   + sep + std::to_string(cpx.dimension()) + ext;
  std::ofstream pers_file(pers_file_name);
  pers.output_diagram(pers_file);
}

int main() {
  Simplex_tree oneRips;
  // Insert 100% of the possible edges, with 2000 vertices -> 1-Rips
  simplex_tree_random_flag_complex(oneRips, 2000, 1.);

  Gudhi::Clock clock_if("initialize_filtration");
  oneRips.initialize_filtration();
  std::clog << clock_if;
  bench_persistence<Persistent_cohomology_stree>(oneRips,         "Persistent_cohomology_oneRips");
  //bench_persistence<Persistence_RU_matrix_vine_stree>(oneRips,    "Persistence_RU_matrix_vine_oneRips");
  //bench_persistence<Persistence_RU_matrix_repr_stree>(oneRips,    "Persistence_RU_matrix_repr_oneRips");
  bench_persistence<Persistence_boundary_matrix_stree>(oneRips,   "Persistence_boundary_matrix_oneRips");
  bench_persistence<Persistence_chain_matrix_stree>(oneRips,      "Persistence_chain_matrix_oneRips");
  bench_persistence<Persistence_chain_matrix_vine_stree>(oneRips, "Persistence_chain_matrix_vine_oneRips");

  for (int max_dim = 1; max_dim < 5; max_dim++) {
    Simplex_tree stree;
    // Insert 15% of the possible edges, with 600 vertices
    simplex_tree_random_flag_complex(stree, 600);

    Gudhi::Clock clock_expansion(std::string("Expansion dim=") + std::to_string(max_dim));
    stree.expansion(max_dim);
    std::clog << clock_expansion;
    Gudhi::Clock clock_if("initialize_filtration");
    stree.initialize_filtration();
    std::clog << clock_if;
    bench_persistence<Persistent_cohomology_stree>(stree,         "Persistent_cohomology_stree");
    //bench_persistence<Persistence_RU_matrix_vine_stree>(stree,    "Persistence_RU_matrix_vine_stree");
    //bench_persistence<Persistence_RU_matrix_repr_stree>(stree,    "Persistence_RU_matrix_repr_stree");
    bench_persistence<Persistence_boundary_matrix_stree>(stree,   "Persistence_boundary_matrix_stree");
    if (max_dim < 2) {
      bench_persistence<Persistence_chain_matrix_stree>(stree,      "Persistence_chain_matrix_stree");
      bench_persistence<Persistence_chain_matrix_vine_stree>(stree, "Persistence_chain_matrix_vine_stree");
    }
  }

  for (int max_dim = 1; max_dim < 5; max_dim++) {
    std::vector<unsigned> sizes;
    std::size_t multipliers = 1;
    const std::size_t SIZE_IN_THIS_DIM = 40;
    for (int dim = 0; dim != max_dim; ++dim) {
      sizes.push_back(SIZE_IN_THIS_DIM);
      multipliers *= SIZE_IN_THIS_DIM;
    }
    std::vector<double> data(multipliers);
    std::generate(data.begin(), data.end(), get_random);

    Bitmap_cubical_complex cub(sizes, data);
    bench_persistence<Persistent_cohomology_cub>(cub,       std::string("Persistent_cohomology_cub_")       + std::to_string(max_dim));
    //bench_persistence<Persistence_RU_matrix_vine_cub>(cub,  std::string("Persistence_RU_vine_matrix_cub_")  + std::to_string(max_dim));
    //bench_persistence<Persistence_RU_matrix_repr_cub>(cub,  std::string("Persistence_RU_repr_matrix_cub_")  + std::to_string(max_dim));
    bench_persistence<Persistence_boundary_matrix_cub>(cub, std::string("Persistence_boundary_matrix_cub_") + std::to_string(max_dim));
    if (max_dim < 4) {
      bench_persistence<Persistence_chain_matrix_cub>(cub,  std::string("Persistence_chain_matrix_cub_")  + std::to_string(max_dim));
      bench_persistence<Persistence_chain_matrix_vine_cub>(cub,  std::string("Persistence_chain_matrix_vine_cub_")  + std::to_string(max_dim));
    }
  }

  return 0;
}
