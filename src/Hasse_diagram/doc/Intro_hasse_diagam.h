/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#ifndef DOC_HASSE_DIAGRAM_H_
#define DOC_HASSE_DIAGRAM_H_

namespace Gudhi {

namespace Hasse_diagram {

/**  \defgroup hasse_diagram Hasse diagram
 *
 * \author   Pawel Dlotko
 *
 * @{
 *
 * A Hasse diagram (implemented in the class Hasse_diagram) is a general 
 * structure to store and operate on chain complexes. 
 * Hasse_diagram is essentially a vector of elements of the class 
 * Hasse_diagram_cell. They implement functionality of a cell (generator) 
 * of an arbitrary chain complex. They allow to access to both boundary 
 * and coboundary of cells elements in a constant time (and get the incidence coefficients). 
 * In addition to that, they store information about the cell like dimension, filtration. Additional 
 * information can also be stored when the template parameter Additional_information_ 
 * is set to a class to store the additional information. Please consult the picture
 * below for a simple example of a Hasse diagram representing two vertices
 * and two edges forming a loop.
 * \image html "Hasse_diag.png" 
 * 
 * Hasse diagram is a dynamic data structure. Addition and removal of cells 
 * can be performed by using add_cell and remove_cell methods. Please consult
 * an example folder for the details of the usage. 
 * The Hasse_diagram class uses a 'lazy delete' philosophy, 
 * i.e. cells that has been deleted are marked as deleted. They can be 
 * physically removed from the data structure by invoking the clean_up_the_structure() 
 * method. Once the number of cells marked as deleted reach the number defined 
 * in Hasse_diagram::proportion_of_removed_cells_that_triggers_reorganization_of_structure
 * (which is set by default to 50\%) this method is called automatically. Note that in this case
 * all the cells are permanently removed from memory and any previous pointers to them
 * became invalid. Note that the user can manipulate the value of 
 * Hasse_diagram::proportion_of_removed_cells_that_triggers_reorganization_of_structure
 * Note that all the iterators on the data structure do not automatically skip the deleted
 * elements. It is up to the user to filter the deleted cells out. 
 * 
 * Objects of a type Hasse_diagram_cell can be constructed by constitutive addition
 * of new cells, by reading the data strucutre from a file (the format is described below).
 * One can also store the intermediate data structures in a file by using write_to_file
 * method. 
 * 
 * The format of a file with the Hasse diagram is as follows:
 * The first line contains a positive integer being the number of cells in the complex.
 * Next lines contains the information about the cells. Information about every single
 * cells is stored in two constitutive lines. 
 * In the first line, the cell id (which is supposed to be a positive integer from the range
 * between 0 and number of cells) that uniquely identify the cell, its dimension are to be
 * provided. In addition to that one can also provide the third parameter which is the filtration
 * of that cells. The second line contains the information about boundaries of the cell. It is a 
 * vector of cell id and the incidence coefficient. 
 * Lines starting from # are ignored when reading the file. 
 * 
 * In order to compute homology or persistent homology of a Hasse diagram, one should
 * use the derived class Hasse_diagram_persistence. It can be used directly with 
 * Gudhi engine to compute (persistent) homology. Please consult the example folder
 * for further details.
 * 
 * All the data structures to store filtered complexes currently available in 
 * Gudhi can be converted to a Hasse_diagrams by using convert_to_Hasse_diagram() 
 * and convert_to_hasse_diagram_persistence(). In the example folder the conversions
 * from Simplex tree and Cubical complex are presented. 
 */
/** @} */  // end defgroup Hasse_diagram

}  // namespace hasse_diagram
}  // namespace Gudhi

#endif  // DOC_HASSE_DIAGRAM_H_
