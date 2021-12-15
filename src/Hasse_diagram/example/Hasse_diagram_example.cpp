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
  using Hasse_diagram = Gudhi::Hasse_diagram::Hasse_diagram_persistence<>;
  using Cell = Hasse_diagram::Cell_type;
  // In this example we will construct a CW decomposition of two dimensional torus:
  //  _______________________
  //  |A|___B___|C|___D___|A|
  //  | |   #   | |   #   | |
  //  |M|###N###|O|###P###|M|
  //  |_|___#___|_|___#___|_|
  //  |I|___J___|K|___L___|I|
  //  | |   #   | |   #   | |
  //  |E|###F###|G|###H###|E|
  //  |_|___#___|_|___#___|_|
  //  |A|___B___|C|___D___|A|
  //
  // Here is the corresponding filtration:
  //  _______________________
  //  |0|___1___|1|___1___|0|
  //  | |   #   | |   #   | |
  //  |1|###2###|1|###2###|1|
  //  |_|___#___|_|___#___|_|
  //  |1|___1___|1|___1___|1|
  //  | |   #   | |   #   | |
  //  |1|###2###|1|###2###|1|
  //  |_|___#___|_|___#___|_|
  //  |0|___1___|1|___1___|0|

  // cells creation
  Cell* A = new Cell(0, 0.0);
  Cell* B = new Cell(1, 1.0);
  Cell* C = new Cell(0, 1.0);
  Cell* D = new Cell(1, 1.0);
  Cell* E = new Cell(1, 1.0);
  Cell* F = new Cell(2, 2.0);
  Cell* G = new Cell(1, 1.0);
  Cell* H = new Cell(2, 2.0);
  Cell* I = new Cell(0, 1.0);
  Cell* J = new Cell(1, 1.0);
  Cell* K = new Cell(0, 1.0);
  Cell* L = new Cell(1, 1.0);
  Cell* M = new Cell(1, 1.0);
  Cell* N = new Cell(2, 2.0);
  Cell* O = new Cell(1, 1.0);
  Cell* P = new Cell(2, 2.0);

  // setting up boundaries and coboundaries of cells:
  // for cell A:
  // Nothing needs to be done, coboundaries will be set up automatically.

   // for cell B:
  auto& boundary_of_B = B->boundaries();
  boundary_of_B.emplace_back(A, 1);
  boundary_of_B.emplace_back(C, 1);

  // for cell C:
  // Nothing needs to be done, coboundaries will be set up automatically.

  // for cell D
  auto& boundary_of_D = D->boundaries();
  boundary_of_D.emplace_back(A, 1);
  boundary_of_D.emplace_back(C, 1);

  // for cell E
  auto& boundary_of_E = E->boundaries();
  boundary_of_E.emplace_back(A, 1);
  boundary_of_E.emplace_back(I, 1);

  // for cell F
  auto& boundary_of_F = F->boundaries();
  boundary_of_F.emplace_back(B, 1);
  boundary_of_F.emplace_back(E, 1);
  boundary_of_F.emplace_back(G, 1);
  boundary_of_F.emplace_back(J, 1);

  // for cell G
  auto& boundary_of_G = G->boundaries();
  boundary_of_G.emplace_back(K, 1);
  boundary_of_G.emplace_back(C, 1);

  // for cell H
  auto& boundary_of_H = H->boundaries();
  boundary_of_H.emplace_back(D, 1);
  boundary_of_H.emplace_back(E, 1);
  boundary_of_H.emplace_back(L, 1);
  boundary_of_H.emplace_back(G, 1);

  // for cell I:
  // Nothing needs to be done, coboundaries will be set up automatically.

  // for cell J
  auto& boundary_of_J = J->boundaries();
  boundary_of_J.emplace_back(I, 1);
  boundary_of_J.emplace_back(K, 1);

  // for cell K:
  // Nothing needs to be done, coboundaries will be set up automatically.

  // for cell L
  auto& boundary_of_L = L->boundaries();
  boundary_of_L.emplace_back(K, 1);
  boundary_of_L.emplace_back(I, 1);

  // for cell M
  auto& boundary_of_M = M->boundaries();
  boundary_of_M.emplace_back(A, 1);
  boundary_of_M.emplace_back(I, 1);

  // for cell N
  auto& boundary_of_N = N->boundaries();
  boundary_of_N.emplace_back(J, 1);
  boundary_of_N.emplace_back(M, 1);
  boundary_of_N.emplace_back(O, 1);
  boundary_of_N.emplace_back(B, 1);

  // for cell O
  auto& boundary_of_O = O->boundaries();
  boundary_of_O.emplace_back(K, 1);
  boundary_of_O.emplace_back(C, 1);

  // for cell P
  auto& boundary_of_P = P->boundaries();
  boundary_of_P.emplace_back(L, 1);
  boundary_of_P.emplace_back(O, 1);
  boundary_of_P.emplace_back(M, 1);
  boundary_of_P.emplace_back(D, 1);

  // we can input the cells in any order into the vector of cells:
  std::vector<Cell*> vect_of_cells = {A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P};

  Hasse_diagram hd(vect_of_cells);

  // Here is a construction of a standard Hasse diagram:
  // Gudhi::Hasse_diagram::Hasse_diagram<Cell> hd( vect_of_cells );
  // std::cout << "Here is the Hasse diagam : " << std::endl << hd << std::endl;

  // Here is a construction of a Hasse_diagram_persistence and computations of
  // persistent homology of the complex above. You should see both the fundamental
  // classes of torus and the generators of the three squares.

  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diagram, Field_Zp> Persistent_cohomology;

  Persistent_cohomology pcoh_hd(hd, true);
  pcoh_hd.init_coefficients(2);
  pcoh_hd.compute_persistent_cohomology(0);

  std::stringstream output_pers;
  pcoh_hd.output_diagram(output_pers);
  std::string obtained_output = output_pers.str();

  std::cout << "Here is the obtained persistence diagram:\n" << obtained_output << std::endl;

  // Needs to delete all cells - it is convenient as they are all in a container
  for (auto cell_ptr : vect_of_cells) {
    delete cell_ptr;
  }

  return 0;
}
