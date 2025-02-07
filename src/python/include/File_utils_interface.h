/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/02 Vincent Rouvreau: Merge Off_reader_interface.h and Reader_utils_interface.h in File_utils_interface.h
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_FILE_UTILS_INTERFACE_H_
#define INCLUDE_FILE_UTILS_INTERFACE_H_

#include <gudhi/Points_off_io.h>
#include <gudhi/reader_utils.h>

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>  // for pair<>

namespace Gudhi {

std::vector<std::vector<double>> read_points_from_OFF_file(const std::string& off_file) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
  return off_reader.get_point_cloud();
}

// Redefine functions with a different name in order the original name can be used in the Python version.
std::vector<std::vector<double>> read_matrix_from_csv_file(const std::string& filename,
                                                                const char separator = ';') {
  return read_lower_triangular_matrix_from_csv_file<double>(filename, separator);
}

inline std::map<int, std::vector<std::pair<double, double>>>
    read_pers_intervals_grouped_by_dimension(std::string const& filename) {
  return read_persistence_intervals_grouped_by_dimension(filename);
}

inline std::vector<std::pair<double, double>>
    read_pers_intervals_in_dimension(std::string const& filename, int only_this_dim = -1) {
  return read_persistence_intervals_in_dimension(filename, only_this_dim);
}

}  // namespace Gudhi

#endif  // INCLUDE_FILE_UTILS_INTERFACE_H_
