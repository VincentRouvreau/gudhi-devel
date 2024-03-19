/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_MATRIX_H
#define PERSISTENCE_MATRIX_H

#include <vector>
#include <iostream>
#include <algorithm>  // for std::sort
#include <iterator>  // for std::advance

#include <gudhi/matrix.h>
#include <gudhi/persistence_matrix_options.h>

namespace Gudhi {
namespace persistence_matrix {

template <class FilteredComplex, class Options = Gudhi::persistence_matrix::Default_options<> >
class Persistence_matrix {
 public:
  explicit Persistence_matrix(FilteredComplex& cpx)
    : cpx_(&cpx), matrix_() {}

  void init_coefficients(int charac) {
    charac_ = charac;
  }

  using Filtration_value = typename FilteredComplex::Filtration_value;
  void compute_persistent_cohomology([[maybe_unused]] Filtration_value min_interval_length = 0) {
    typename FilteredComplex::Simplex_key idx {};
    for (auto f_simplex : cpx_->filtration_simplex_range()) {
      cpx_->assign_key(f_simplex, idx++);
      std::vector<int> boundary_indexes = {};
#ifdef DEBUG_TRACES
      std::clog << "   { ";
#endif  // DEBUG_TRACES
      for (auto b_simplex : cpx_->boundary_simplex_range(f_simplex)) {
        // As the filtration is sorted, the boundary simplex key is set
        boundary_indexes.push_back(cpx_->key(b_simplex));
#ifdef DEBUG_TRACES
        std::clog << cpx_->key(b_simplex) << ", ";
#endif  // DEBUG_TRACES
      }
#ifdef DEBUG_TRACES
      std::clog << "}\n";
#endif  // DEBUG_TRACES
      std::sort(boundary_indexes.begin(), boundary_indexes.end());
      matrix_.insert_boundary(boundary_indexes);
    }
  }

  using Barcode = typename Matrix<Options>::Bar;

  struct cmp_intervals_by_length {
    explicit cmp_intervals_by_length(FilteredComplex * sc)
        : sc_(sc),
          filtration_size_(std::distance(sc_->filtration_simplex_range().begin(),
                                         sc_->filtration_simplex_range().end())) {
    }
    bool operator()(const Barcode & bc1, const Barcode & bc2) {
      //if (bc1.dim != bc2.dim) return bc1.dim > bc2.dim;
      if (bc1.death > filtration_size_) return true;
      if (bc2.death > filtration_size_) return false;
      auto filtration_begin = sc_->filtration_simplex_range().begin();
      return (sc_->filtration(*(filtration_begin + bc1.death)) - sc_->filtration(*(filtration_begin + bc1.birth))
            > sc_->filtration(*(filtration_begin + bc2.death)) - sc_->filtration(*(filtration_begin + bc2.birth)));
    }
    FilteredComplex * sc_;
    std::size_t filtration_size_;
  };

  void output_diagram(std::ostream& ostream = std::cout) {
    /*auto filtration_begin = cpx_->filtration_simplex_range().begin();
    auto filtration_size = std::distance(filtration_begin, cpx_->filtration_simplex_range().end());

    auto barcodes = matrix_.get_current_barcode();
    cmp_intervals_by_length cmp(cpx_);
    std::sort(std::begin(barcodes), std::end(barcodes), cmp);

    for (const Barcode& barcode : barcodes) {
      ostream << charac_ << "  " << barcode.dim << " "
        << cpx_->filtration(*(filtration_begin + barcode.birth)) << " ";
        if (barcode.death > filtration_size) {
          ostream << "inf" << std::endl;
        } else {
          ostream << cpx_->filtration(*(filtration_begin + barcode.death)) << std::endl;
        } 
    }*/
  }


 private:
  FilteredComplex* cpx_;
  Matrix<Options> matrix_;
  int charac_;
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PERSISTENCE_MATRIX_H
