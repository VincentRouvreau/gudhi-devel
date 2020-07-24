/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_EDGE_COLLAPSE_INTERFACE_H_
#define INCLUDE_EDGE_COLLAPSE_INTERFACE_H_

#include <gudhi/Flag_complex_edge_collapser.h>

#include <vector>
#include <tuple>

namespace Gudhi {

namespace collapse {

// Redefine functions with a different name in order the original name can be used in the Python version.
std::vector<std::tuple<int, int, double>> edge_collapse(const std::vector<std::tuple<int, int, double>>& edges) {
  return flag_complex_collapse_edges(edges);
}

}  // namespace collapse

}  // namespace Gudhi


#endif  // INCLUDE_EDGE_COLLAPSE_INTERFACE_H_
