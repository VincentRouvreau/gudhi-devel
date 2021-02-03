/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Jamin
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_COMMON_FILE_FORMAT_H_
#define DOC_COMMON_FILE_FORMAT_H_

namespace Gudhi {

/*! \page fileformats File formats

 \tableofcontents

 \section FileFormatsOFF OFF file format

 OFF files must be conform to format described here: http://www.geomview.org/docs/html/OFF.html

 OFF files are mainly used as point cloud inputs. Here is an example of 7 points in a 3-dimensional space. As edges and
 faces are not used for point set, there is no need to specify them (just set their numbers to 0):

 \include points/alphacomplexdoc.off

 For dimensions bigger than 3, the dimension can be set like here:
 \verbatim
  # Dimension is no more 3
  nOFF
  # dimension 4 - 7 vertices - 0 face - 0 edge
  4 7 0 0
  # Point set:
  1.0 1.0  0.0 0.0
  7.0 0.0  0.0 0.0
  4.0 6.0  0.0 0.0
  9.0 6.0  0.0 0.0
  0.0 14.0 0.0 0.0
  2.0 19.0 0.0 0.0
  9.0 17.0 0.0 0.0
 \endverbatim


 \section FileFormatsPers Persistence Diagram

 Such a file, whose extension is usually `.pers`, contains a list of persistence intervals.<br>
 Lines starting with `#` are ignored (comments).<br>
 Other lines might contain 2, 3 or 4 values (the number of values on each line must be the same for all lines):
 \verbatim
   [[field] dimension] birth death
 \endverbatim

 Here is a simple sample file:
 \verbatim
   # Persistence diagram example
   2 2.7 3.7
   2 9.6 14.
   # Some comments
   3 34.2 34.974
   4 3. inf
 \endverbatim

 Other sample files can be found in the `data/persistence_diagram` folder.

 Such files can be generated with `Gudhi::persistent_cohomology::Persistent_cohomology::output_diagram()` and read with
 `Gudhi::read_persistence_intervals_and_dimension()`, `Gudhi::read_persistence_intervals_grouped_by_dimension()` or
 `Gudhi::read_persistence_intervals_in_dimension()`.
 

 \section FileFormatsIsoCuboid Iso-cuboid

 Such a file describes an iso-oriented cuboid with diagonal opposite vertices (min_x, min_y, min_z,...) and (max_x, max_y, max_z, ...). The format is:<br>
 \verbatim
   min_x min_y [min_z ...]
   max_x max_y [max_z ...]
 \endverbatim

 Here is a simple sample file in the 3D case:
 \verbatim
   -1. -1. -1.
   1. 1. 1.
 \endverbatim


 \section FileFormatsPerseus Perseus

 This file format is a format inspired from the Perseus software
 (http://www.sas.upenn.edu/~vnanda/perseus/) by Vidit Nanda.
 The first line contains a number d begin the dimension of the
 bitmap (2 in the example below). Next d lines are the numbers of top dimensional cubes in each dimensions (3 and 3
 in the example below). Next, in lexicographical order, the filtration of top dimensional cubes is given (1 4 6 8
 20 4 7 6 5 in the example below).
 
 \image html "exampleBitmap.png" "Example of a input data."
 
 The input file for the following complex is:
 \verbatim
 2
 3
 3
 1
 4
 6
 8
 20
 4
 7
 6
 5
 \endverbatim

 To indicate periodic boundary conditions in a
 given direction, then number of top dimensional cells in this direction have to be multiplied by -1. For instance:

 \verbatim
 2
 -3
 3
 1
 4
 6
 8
 20
 4
 7
 6
 5
 \endverbatim

 Indicate that we have imposed periodic boundary conditions in the direction x, but not in the direction y.

 Other sample files can be found in the `data/bitmap` folder.
 
 \note Unlike in Perseus format the filtration on the maximal cubes can be any double precision number.
 Consequently one cannot mark the cubes that are not present with `-1`'s. To do that please set their filtration value
 to \f$+\infty\f$ (aka. `inf` in the file).

  
 \section FileFormatHasseDiagram Format of file based on which Hasse diagram can be created. 
  Such a Hasse diagram can be read by Hasse_diagram/utilities/Hasse_diagram_from_file.cpp
  Lines starting with `#` are ignored (comments).
  We assume that the file do not contain empty lines.
  
  The cells stored in that file are assumed to be enumerated with integers starting from zero.
  
  The first line contains a non negative integer N determining a number of cells in the Hasse diagram.  
  Subsequently the file contains N blocks. Each block represent a cell in the chain complex. Below a description 
  of a block is given.
   
  First line of a block consist of two or three numbers: number of a cell (between zero and N-1), dimension of cell (nonnegative integer). 
  The third (optional) number is a filtration of a cell.
  The next line contains sequence of ids of boundary elements of a given cell alternated by the incidence coefficient between the given cell and the boundary element.
 
  For an exampe of a file, please consult Hasse_diagram/test/cw_decomposition_of_torus.hasse    
*/
}  // namespace Gudhi

#endif  // DOC_COMMON_FILE_FORMAT_H_
