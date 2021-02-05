#include <gudhi/Hasse_diagram.h>

#include <iostream>
#include <vector>

int main(int argc, char** argv) 
{ 
  using Cell = Gudhi::hasse_diagram::Hasse_diagram_cell<int, double, double>;
  Cell* vertex1 = new Cell(0, 0.0);
  Cell* vertex2 = new Cell(0, 0.0);
  Cell* edge1   = new Cell(1, 1.0);
  Cell* edge2   = new Cell(1, 1.0);

  // Nothing needs to be done for vertices, coboundaries will be set up automatically.

  edge1->get_boundary().emplace_back(vertex1,  1);
  edge1->get_boundary().emplace_back(vertex2,  1);
  edge2->get_boundary().emplace_back(vertex1, -1);
  edge2->get_boundary().emplace_back(vertex2, -1);

  // we can input the cells in any order into the vector of cells:
  std::vector< Cell* > cells_vector = {vertex1, vertex2, edge1, edge2};

  using Hasse_diagram = Gudhi::hasse_diagram::Hasse_diagram<Cell>;
  Hasse_diagram hasse_diag(cells_vector);
  std::cout << hasse_diag.full_signature_of_the_structure() << std::endl;

  for (auto cell : cells_vector)
    delete cell;
}
