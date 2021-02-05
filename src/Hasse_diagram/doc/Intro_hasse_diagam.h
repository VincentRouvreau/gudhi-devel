/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University, UK
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DOC_HASSE_DIAGRAM_H_
#define DOC_HASSE_DIAGRAM_H_

namespace Gudhi {

namespace hasse_diagram {

/**  \defgroup hasse_diagram Hasse diagram
 *
 * \author   Pawel Dlotko
 *
 * @{
 *
 * \section HasseDiagramIntroduction Introduction
 *
 * A Hasse diagram (implemented in the class Hasse_diagram) is a general
 * structure to store and operate on <a href="https://en.wikipedia.org/wiki/Chain_complex">chain complexes</a>.
 * Hasse_diagram is essentially a vector of elements of the class
 * Hasse_diagram_cell. They implement functionality of a cell (generator)
 * of an arbitrary chain complex. They allow to access to both boundary
 * and coboundary of cells elements in a constant time (and get the incidence coefficients).
 * In addition to that, they store information about the cell like dimension, filtration. Additional
 * information can also be stored when the template parameter Additional_information_
 * is set to a class to store the additional information.
 *
 * \section HasseDiagramBasicExample Basic example
 *
 * Please consult the picture and the corresponding code below for a simple example of a Hasse diagram representing
 * two vertices and two edges forming a loop.
 * \image html "Hasse_diag.png"
 *
 * \include Hasse_diagram/Hasse_diagram_basic_example.cpp
 *
 * \section HasseDiagramDetails Implementation details
 *
 * \subsection HasseDiagramAddDelete Cells insertions and removals
 *
 * Hasse diagram is a dynamic data structure. Addition and removal of cells
 * can be performed by using add_cell and remove_cell methods.
 * cf. <a href="_hasse_diagram_2_insertion_and_removal_operations_on__hasse_diagrams_8cpp-example.html">
 * Insertion_and_removal_operations_on_Hasse_diagrams.cpp</a>
 *
 * The Hasse_diagram class uses a 'lazy delete' philosophy,
 * i.e. cells that has been deleted are marked as deleted. They can be
 * physically removed from the data structure by invoking the
 * `Gudhi::hasse_diagram::Hasse_diagram::clean_up_the_structure()` method. Once the number of cells marked as deleted
 * reach the number defined in
 * `Gudhi::hasse_diagram::Hasse_diagram::proportion_of_removed_cells_that_triggers_reorganization_of_structure` (which
 * is set by default to 50\%) this method is called automatically. Note that in this case all the cells are permanently
 * removed from memory and any previous pointers to them became invalid. Note that the user can manipulate the value of
 * `Gudhi::hasse_diagram::Hasse_diagram::proportion_of_removed_cells_that_triggers_reorganization_of_structure`.
 * Note that all the iterators on the data structure do not automatically skip the deleted
 * elements. It is up to the user to filter the deleted cells out.
 *
 * \subsection HasseDiagramCellsFromFile Read cells from a file
 *
 * Objects of a type `Hasse_diagram_cell` can be constructed by consecutive addition
 * of new cells, by reading the data structure from a file (the format is described on
 * <a href="fileformats.html#FileFormatHasseDiagram">this page</a>).
 * One can also store the intermediate data structures in a file by using write_to_file
 * method.
 * cf. <a href="_hasse_diagram_2_hasse_diagram_from_file_8cpp-example.html">Hasse_diagram_from_file.cpp</a>.
 *
 * \subsection HasseDiagramHomology Hasse diagram homology
 *
 * In order to compute homology or persistent homology of a Hasse diagram, one should
 * use the derived class `Gudhi::hasse_diagram::Hasse_diagram_persistence`. It can be used directly with
 * `Gudhi::persistent_cohomology::Persistent_cohomology`.
 * Please consult the example
 * <a href="_hasse_diagram_2_hasse_diagram_torus_example_8cpp-example.html">Hasse_diagram_torus_example.cpp</a>.
 *
 * \subsection HasseDiagramConversion Hasse diagram conversion
 *
 * Data structures to store filtered complexes can be converted to a `Hasse_diagram` by using
 * `Gudhi::hasse_diagram::convert_to_Hasse_diagram()`
 * and `Gudhi::hasse_diagram::convert_to_Hasse_diagram_persistence()`.
 * cf. <a
 * href="_hasse_diagram_2_hasse_diagram_from_simplex_tree_8cpp-example.html">Hasse_diagram_from_simplex_tree.cpp</a> and
 * <a
 * href="_hasse_diagram_2_hasse_diagram_from_cubical_complex_8cpp-example.html">Hasse_diagram_from_cubical_complex.cpp</a>.
 */
/** @} */  // end defgroup Hasse_diagram

}  // namespace hasse_diagram
}  // namespace Gudhi

#endif  // DOC_HASSE_DIAGRAM_H_
