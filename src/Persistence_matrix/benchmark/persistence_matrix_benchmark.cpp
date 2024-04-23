/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistence_matrix.h>
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/random_graph_generators.h>
#include <gudhi/Clock.h>  // for Gudhi::Clock
#include <gudhi/Points_off_io.h>  // for Gudhi::Points_off_reader
#include <gudhi/pick_n_random_points.h>

#include "random_filtrations.h"

#include <CGAL/Epeck_d.h>
 
#include <boost/mpl/list.hpp>
#include <boost/mpl/for_each.hpp>

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

std::vector<std::string> column_type_str = {"INTRUSIVE_LIST", "SET", "HEAP", "VECTOR", "NAIVE_VECTOR", "UNORDERED_SET", "LIST", "INTRUSIVE_SET"};

template<typename Point_d>
std::vector<Point_d> random_activities(std::string off_file, std::size_t points_number) {
  std::vector<Point_d> results;
  Gudhi::Points_off_reader<Point_d> off_reader(off_file);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_file << "\n";
    return results;  // ----- >>
  }
  std::vector<Point_d> point_cloud = off_reader.get_point_cloud();
  Gudhi::subsampling::pick_n_random_points(point_cloud, points_number, std::back_inserter(results));
  return results;
}

template <class Persistence, class Complex>
void bench_persistence(Complex& cpx, const std::string& msg, bool& timeout) {
  if (!timeout) {
    Gudhi::Clock clock_pers(msg);
    Persistence pers(cpx, true);
    pers.init_coefficients(2);
    pers.compute_persistent_cohomology(0.0);
    clock_pers.end();
    std::clog << cpx.num_simplices() << " simplices - " << clock_pers;
  
    // If it took more than 10 seconds, sets timeout to true for not to be launched at next loop
    if (clock_pers.num_seconds() > 10.) {
      timeout = true;
    }
  
    /*std::string sep("_");
    std::string ext(".pers");
    std::string pers_file_name = msg + sep + std::to_string(cpx.num_simplices())
                                     + sep + std::to_string(cpx.dimension()) + ext;
    std::ofstream pers_file(pers_file_name);
    pers.output_diagram(pers_file);*/
  }
}

template<typename Default_option>
struct RU_matrix_options_vine : Default_option {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = true;
};
    
template<typename Default_option>
struct RU_matrix_options_repr : Default_option {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = true;
  static const bool is_of_boundary_type = true;
};
    
template<typename Default_option>
struct Boundary_matrix_options : Default_option {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = true;
};
    
template<typename Default_option>
struct Chain_matrix_options : Default_option {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = false;
};
    
template<typename Default_option>
struct Chain_matrix_options_vine : Default_option {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = false;
  static const bool is_of_boundary_type = false;
};

struct Bench_persistence {
  static inline int n {0};

  template<typename Default_option> void operator()(Default_option) {
    std::clog << " ****************************** " << column_type_str[n++] << " ******************************\n";

    using Persistence_RU_matrix_vine_stree    = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, RU_matrix_options_vine<Default_option>>;
    using Persistence_RU_matrix_repr_stree    = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, RU_matrix_options_repr<Default_option>>;
    using Persistence_boundary_matrix_stree   = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, Boundary_matrix_options<Default_option>>;
    using Persistence_chain_matrix_stree      = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, Chain_matrix_options<Default_option>>;
    using Persistence_chain_matrix_vine_stree = Gudhi::persistence_matrix::Persistence_matrix<Gudhi::Simplex_tree<>, Chain_matrix_options_vine<Default_option>>;
    
    using Persistence_RU_matrix_vine_cub    = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, RU_matrix_options_vine<Default_option>>;
    using Persistence_RU_matrix_repr_cub    = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, RU_matrix_options_repr<Default_option>>;
    using Persistence_boundary_matrix_cub   = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, Boundary_matrix_options<Default_option>>;
    using Persistence_chain_matrix_cub      = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, Chain_matrix_options<Default_option>>;
    using Persistence_chain_matrix_vine_cub = Gudhi::persistence_matrix::Persistence_matrix<Bitmap_cubical_complex, Chain_matrix_options_vine<Default_option>>;

    bool cohomology_timeout, RU_matrix_vine_timeout, RU_matrix_repr_timeout, boundary_matrix_timeout, chain_matrix_timeout, chain_matrix_vine_timeout;

    cohomology_timeout = RU_matrix_vine_timeout = RU_matrix_repr_timeout = boundary_matrix_timeout = chain_matrix_timeout = chain_matrix_vine_timeout = false;
    for (int points_number = 10000; points_number <= 30000; points_number+=20000) {
      for (int dim = 2; dim <= 4; dim++) {
        Gudhi::Simplex_tree<> stree;
        build_random(stree, points_number, points_number * dim * 2, 2);
        bench_persistence<Persistent_cohomology_stree>(stree,         std::string("Persistent_cohomology_random_filtrations________") + std::to_string(points_number) + std::string("_") + std::to_string(dim), cohomology_timeout);
        bench_persistence<Persistence_boundary_matrix_stree>(stree,   std::string("Persistence_boundary_matrix_random_filtrations__") + std::to_string(points_number) + std::string("_") + std::to_string(dim), boundary_matrix_timeout);
        bench_persistence<Persistence_chain_matrix_stree>(stree,      std::string("Persistence_chain_matrix_random_filtrations_____") + std::to_string(points_number) + std::string("_") + std::to_string(dim), chain_matrix_timeout);
        bench_persistence<Persistence_chain_matrix_vine_stree>(stree, std::string("Persistence_chain_matrix_vine_random_filtrations") + std::to_string(points_number) + std::string("_") + std::to_string(dim), chain_matrix_vine_timeout);
        bench_persistence<Persistence_RU_matrix_vine_stree>(stree,    std::string("Persistence_RU_matrix_vine_random_filtrations___") + std::to_string(points_number) + std::string("_") + std::to_string(dim), RU_matrix_vine_timeout);
        bench_persistence<Persistence_RU_matrix_repr_stree>(stree,    std::string("Persistence_RU_matrix_repr_random_filtrations___") + std::to_string(points_number) + std::string("_") + std::to_string(dim), RU_matrix_repr_timeout);
      }
    }

    // ************************* ALPHA COMPLEX *************************
    using Kernel = CGAL::Epeck_d< CGAL::Dimension_tag<3> >;
    using Point = Kernel::Point_d;
    cohomology_timeout = RU_matrix_vine_timeout = RU_matrix_repr_timeout = boundary_matrix_timeout = chain_matrix_timeout = chain_matrix_vine_timeout = false;
    for (int points_number = 10000; points_number <= 25000; points_number+=5000) {
      std::vector<Point> points = random_activities<Point>("./activities.off", points_number);

      Gudhi::alpha_complex::Alpha_complex<Kernel> alpha(points);
      Gudhi::Simplex_tree<> stree;
      alpha.create_complex(stree);
      bench_persistence<Persistent_cohomology_stree>(stree,         std::string("Persistent_cohomology_alpha________") + std::to_string(points_number), cohomology_timeout);
      bench_persistence<Persistence_boundary_matrix_stree>(stree,   std::string("Persistence_boundary_matrix_alpha__") + std::to_string(points_number), boundary_matrix_timeout);
      bench_persistence<Persistence_chain_matrix_stree>(stree,      std::string("Persistence_chain_matrix_alpha_____") + std::to_string(points_number), chain_matrix_timeout);
      bench_persistence<Persistence_chain_matrix_vine_stree>(stree, std::string("Persistence_chain_matrix_vine_alpha") + std::to_string(points_number), chain_matrix_vine_timeout);
      bench_persistence<Persistence_RU_matrix_vine_stree>(stree,    std::string("Persistence_RU_matrix_vine_alpha___") + std::to_string(points_number), RU_matrix_vine_timeout);
      bench_persistence<Persistence_RU_matrix_repr_stree>(stree,    std::string("Persistence_RU_matrix_repr_alpha___") + std::to_string(points_number), RU_matrix_repr_timeout);
    }

    // ************************* ONE RIPS *************************
    cohomology_timeout = RU_matrix_vine_timeout = RU_matrix_repr_timeout = boundary_matrix_timeout = chain_matrix_timeout = chain_matrix_vine_timeout = false;
    Simplex_tree oneRips;
    // Insert 100% of the possible edges, with 150 vertices -> 1-Rips
    simplex_tree_random_flag_complex(oneRips, 150, 1.);
    
    
    Gudhi::Clock clock_if("initialize_filtration");
    // Triangles
    oneRips.expansion(2);
    oneRips.initialize_filtration();
    std::clog << clock_if;
    bench_persistence<Persistent_cohomology_stree>(oneRips,         "Persistent_cohomology_one_rips________", cohomology_timeout       );
    bench_persistence<Persistence_boundary_matrix_stree>(oneRips,   "Persistence_boundary_matrix_one_rips__", boundary_matrix_timeout  );
    bench_persistence<Persistence_chain_matrix_stree>(oneRips,      "Persistence_chain_matrix_one_rips_____", chain_matrix_timeout     );
    bench_persistence<Persistence_chain_matrix_vine_stree>(oneRips, "Persistence_chain_matrix_vine_one_rips", chain_matrix_vine_timeout);
    /*bench_persistence<Persistence_RU_matrix_vine_stree>(oneRips,    "Persistence_RU_matrix_vine_one_rips___", RU_matrix_vine_timeout   );
    bench_persistence<Persistence_RU_matrix_repr_stree>(oneRips,    "Persistence_RU_matrix_repr_one_rips___", RU_matrix_repr_timeout   );*/

    // ************************* RANDOM FLAG COMPLEX *************************
    cohomology_timeout = RU_matrix_vine_timeout = RU_matrix_repr_timeout = boundary_matrix_timeout = chain_matrix_timeout = chain_matrix_vine_timeout = false;
    for (int max_dim = 2; max_dim < 5; max_dim++) {
      Simplex_tree stree;
      // Insert 15% of the possible edges, with 200 vertices
      simplex_tree_random_flag_complex(stree, 200);
    
      Gudhi::Clock clock_expansion(std::string("Expansion dim=") + std::to_string(max_dim));
      stree.expansion(max_dim);
      std::clog << clock_expansion;
      Gudhi::Clock clock_if("initialize_filtration");
      stree.initialize_filtration();
      std::clog << clock_if;
      bench_persistence<Persistent_cohomology_stree>(stree,         std::string("Persistent_cohomology_random_graph________") + std::to_string(max_dim), cohomology_timeout);
      bench_persistence<Persistence_boundary_matrix_stree>(stree,   std::string("Persistence_boundary_matrix_random_graph__") + std::to_string(max_dim), boundary_matrix_timeout);
      bench_persistence<Persistence_chain_matrix_stree>(stree,      std::string("Persistence_chain_matrix_random_graph_____") + std::to_string(max_dim), chain_matrix_timeout);
      bench_persistence<Persistence_chain_matrix_vine_stree>(stree, std::string("Persistence_chain_matrix_vine_random_graph") + std::to_string(max_dim), chain_matrix_vine_timeout);
      bench_persistence<Persistence_RU_matrix_vine_stree>(stree,    std::string("Persistence_RU_matrix_vine_random_graph___") + std::to_string(max_dim), RU_matrix_vine_timeout);
      bench_persistence<Persistence_RU_matrix_repr_stree>(stree,    std::string("Persistence_RU_matrix_repr_random_graph___") + std::to_string(max_dim), RU_matrix_repr_timeout);
    }
    
    // ************************* RANDOM CUBICAL *************************
    cohomology_timeout = RU_matrix_vine_timeout = RU_matrix_repr_timeout = boundary_matrix_timeout = chain_matrix_timeout = chain_matrix_vine_timeout = false;
    for (int max_dim = 1; max_dim < 4; max_dim++) {
      std::vector<unsigned> sizes;
      std::size_t multipliers = 1;
      // num_simplices = ((2 * N) + 1) ^ max_dim
      const std::size_t SIZE_IN_THIS_DIM = 200 / (max_dim * 2);
      for (int dim = 0; dim != max_dim; ++dim) {
        sizes.push_back(SIZE_IN_THIS_DIM);
        multipliers *= SIZE_IN_THIS_DIM;
      }
      std::vector<double> data(multipliers);
      std::generate(data.begin(), data.end(), get_random);
    
      Bitmap_cubical_complex cub(sizes, data);
      bench_persistence<Persistent_cohomology_cub>(cub,         std::string("Persistent_cohomology_random_cubical_________") + std::to_string(max_dim), cohomology_timeout);
      bench_persistence<Persistence_boundary_matrix_cub>(cub,   std::string("Persistence_boundary_matrix_random_cubical___") + std::to_string(max_dim), boundary_matrix_timeout);
      bench_persistence<Persistence_chain_matrix_cub>(cub,      std::string("Persistence_chain_matrix_random_cubical______") + std::to_string(max_dim), chain_matrix_timeout);
      bench_persistence<Persistence_chain_matrix_vine_cub>(cub, std::string("Persistence_chain_matrix_vine_random_cubical_") + std::to_string(max_dim), chain_matrix_vine_timeout);
      bench_persistence<Persistence_RU_matrix_vine_cub>(cub,    std::string("Persistence_RU_vine_matrix_random_cubical____") + std::to_string(max_dim), RU_matrix_vine_timeout);
      bench_persistence<Persistence_RU_matrix_repr_cub>(cub,    std::string("Persistence_RU_repr_matrix_random_cubical____") + std::to_string(max_dim), RU_matrix_repr_timeout);
    }
  }
};

int main() {
  typedef boost::mpl::list<Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::INTRUSIVE_LIST, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::SET, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::HEAP, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::VECTOR, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::UNORDERED_SET, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::LIST, true>,
                           Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::INTRUSIVE_SET, true>         
                           > Default_options;
  boost::mpl::for_each<Default_options>(Bench_persistence());

  return 0;
}

