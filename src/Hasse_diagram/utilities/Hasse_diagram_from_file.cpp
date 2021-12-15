/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "This program computes (persistent) homology of a chain complex given in a file. The format of the file "
                 "is the following:"
              << std::endl;
    std::cout << "In the first line, the number of all cells is given." << std::endl;
    std::cout << "The cells are assumed to be enumerated from 0 to the number of cells." << std::endl;
    std::cout << "Each cell is given in the following format:" << std::endl;
    std::cout << "Id of a cell followed by its dimension and optionally a filtration (in the first line)" << std::endl;
    std::cout << "Sequence of ids of boundary elements followed by the incidence coeficient between given cell and the "
                 "boundary element (all of them in the second line)"
              << std::endl
              << std::endl;
    std::cout << "The input parameters of the program are: " << std::endl;
    std::cout << "(1) Name of the file with complex, " << std::endl;
    std::cout << "(2) Optional -- name of the output file, " << std::endl;
    std::cout << "(3) Optional -- prime number p such that (persistent) homology over Zp will be computed. " << std::endl;

    return 1;
  }

  const char* filename = argv[1];
  const char* output_file = "output";
  unsigned field_characteristic = 11;
  if (argc > 2) {
    output_file = argv[2];
  }
  if (argc > 3) {
    field_characteristic = (unsigned)(atoi(argv[3]));
  }

  using Hasse_diagram = Gudhi::Hasse_diagram::Hasse_diagram_persistence<>;

  Hasse_diagram hd(filename);
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diagram, Field_Zp>;

  Persistent_cohomology pcoh(hd, true);
  double min_persistence = 0.;

  pcoh.init_coefficients(field_characteristic);
  pcoh.compute_persistent_cohomology(min_persistence);

  std::ofstream out(output_file);
  pcoh.output_diagram(out);
  out.close();

  std::cout << "Result in file: " << output_file << "\n";

  return 0;
}
