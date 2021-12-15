/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Hasse_diagram.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Bitmap_cubical_complex.h>

// standard stuff
#include <iostream>
#include <vector>

int main() {
  // This is an example of construction of Hasse diagram from a cubical complex. In the example
  // below we first define 3 by 3 by 3 cubical complex and assign a filtration on it. Later
  // we convert it to a standard Hasse_diagram and display it. At the end, we will convert it
  // to an object of a class Hasse_diagram_persistence and compute persitence of it.
  using Cubical_complex_base = Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>;
  using Cubical_complex = Gudhi::cubical_complex::Bitmap_cubical_complex<Cubical_complex_base>;
  std::vector<unsigned> sizes = {3, 3, 3};
  // this is 3 by 3 by 3 cubical complex representing a cubical sphere embedded on two dimensional torus.
  // The lifespan of a sphere is from 1 to 10.

  std::vector<double> top_dimensional_cells_data = {
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    //
    1, 1, 1,
    1, 10, 1,
    1, 1, 1,
    //
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  };
  Cubical_complex b(sizes, top_dimensional_cells_data);

  using Hasse_diagram = Gudhi::Hasse_diagram::Hasse_diagram_persistence<>;
  Hasse_diagram* hd = Gudhi::Hasse_diagram::convert_to_hasse_diagram_persistence<Cubical_complex, Hasse_diagram>(b);

  std::cout << "Here is the Hasse diagram obtained from the cubical complex : " << *hd << std::endl;

  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diagram, Field_Zp>;

  Persistent_cohomology pcoh(*hd, true);
  int field_characteristic = 11;
  double min_persistence = 0;

  pcoh.init_coefficients(field_characteristic);
  pcoh.compute_persistent_cohomology(min_persistence);

  std::cout << "And here is the persistent homology of the Hasse diagram : " << std::endl;
  pcoh.output_diagram();

  return 0;
}
