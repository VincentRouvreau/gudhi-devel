/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <random>
#include <iostream>
#include <set>

template<typename Simplex_tree>
void build_random(Simplex_tree& st, unsigned int numberOfPoints, unsigned int numberOfSimplices, unsigned int maxDimension, int seed = -1) {
  std::random_device dev;
  std::mt19937 rng(dev());
  if (seed > -1) rng.seed(seed);
  std::uniform_int_distribution<int> vertexDist(0,numberOfPoints-1);
  std::uniform_int_distribution<int> dimDist(maxDimension-2,maxDimension);
  unsigned int currentFiltrationValue = 1;

  for (int vertex = 0; vertex < static_cast<int>(numberOfPoints); ++vertex) {
    st.insert_simplex({vertex});
  }

  unsigned int size = 0;
  while (size < numberOfSimplices && currentFiltrationValue < numberOfSimplices){
    unsigned int dim = dimDist(rng);
    std::set<int> simplex;

    while (simplex.size() < dim + 1){
      simplex.insert(vertexDist(rng));
    }

    st.insert_simplex_and_subfaces(simplex, currentFiltrationValue++);

    if (currentFiltrationValue % (numberOfSimplices / 20) == 0)
      size = st.num_simplices();
  }

  /*std::cout << "Number of points: " << numberOfPoints << "\n";
  std::cout << "Number of simplices: " << st.num_simplices() << "\n";
  std::cout << "Maximum dimension: " << st.dimension() << "\n";
  std::cout << "Seed: " << (seed > -1 ? std::to_string(seed) : "none") << "\n";*/
}
