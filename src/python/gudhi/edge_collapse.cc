/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Flag_complex_edge_collapser.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>
#include <tuple>

namespace py = pybind11;

std::vector<std::tuple<int, int, double>> edge_collapse(const std::vector<std::tuple<int, int, double>>& edges)
{
  return Gudhi::collapse::flag_complex_collapse_edges(edges);
}

PYBIND11_MODULE(edge_collapse, m) {
      m.attr("__license__") = "MIT";
      m.attr("__copyright__") = "Copyright (C) 2020 Inria";
      m.def("flag_complex_collapse_edges", &edge_collapse,
          py::arg("edges"),
          R"pbdoc(
    Reduces a filtration of Vietoris-Rips complex from its graph to another smaller flag filtration with same
    persistence.

    :param edges: Edges to collapse.
    :type edges: vector of tuples of shape(int, int, float)

    :returns: The remaining edges after collapse.
    :rtype: vector of tuples of shape(int, int, float)
    )pbdoc");
}
