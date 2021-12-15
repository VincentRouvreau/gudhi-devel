/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef HASSE_DIAGRAM_PERSISTENCE_H
#define HASSE_DIAGRAM_PERSISTENCE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <numeric>  // for std::iota
#include <limits>  // for std::numeric_limits<>
#include <utility>  // for std::pair<>
#include <cmath>  // for std::fabs

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <gudhi/Hasse_diagram_cell.h>
#include <gudhi/Hasse_diagram.h>

namespace Gudhi {

namespace Hasse_diagram {

template <typename HasseDiagramOptions>
class is_before_in_filtration;

/** Model of HasseDiagramOptions required to compute persistence.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Hasse_diagram_options_for_persistence : public Hasse_diagram_options_full_featured {
  using Filtration_value = double;
  using Simplex_key = std::uint32_t;
  // Mandatory to compute persistence
  static const bool store_key = true;
  static const bool store_filtration = true;
  // Not required, but can be changed if needed
  static const bool store_additional_information = false;
  using Additional_information = bool;
};

/**
 * \class Hasse_diagram_persistence
 * \brief Data structure to store Hasse diagrams and compute its persistence diagrams.
 *
 * \ingroup Hasse_diagram
 *
 * \details
 * This is a data structure derived from Hasse_diagram. The additional functionalities
 * are the ones required by Gudhi for persistent homology computations. Please
 * refer to Hasse_diagram class for further details
 *
 * Please refer to \ref Hasse_diagram for examples.
 *
 * The complex is a template class requiring the following parameters:
 * Cell_type - a parameter describing a cell of Hasse diagram. Please refer to Hasse_diagram_cell.h for further details.
 *
 */
template <typename HasseDiagramOptions = Hasse_diagram_options_for_persistence>
class Hasse_diagram_persistence : public Hasse_diagram<HasseDiagramOptions> {
 public:
  using Cell_type = typename Hasse_diagram<HasseDiagramOptions>::Cell_type;
  /**
   * Default constructor.
   **/
  Hasse_diagram_persistence() : Hasse_diagram<HasseDiagramOptions>() {}

  /**
   * Creating Hasse diagram for persistence computations from a file.
   * The file format is the following:
   * Number of cells
   * cell dimension
   * ids of cell boundary elements followed by the incidence coefficient.
   * the two lines above are repeated for each cell.
   * It is assumed that the id of a cell is its position in the file.
   **/
  Hasse_diagram_persistence(const char* filename) : Hasse_diagram<HasseDiagramOptions>(filename) { this->set_up_the_arrays(); }

  /**
   * Constructor to create a Hasse diagram for persistence computations
   * from a vector of cells. It is assumed that all the cells have
   * boundaries set up. Setting up the coboundaries will be done in the
   * constructor based on the information about boundaries.
   **/
  Hasse_diagram_persistence(const std::vector<Cell_type*>& cells_) : Hasse_diagram<HasseDiagramOptions>(cells_) {
    this->set_up_the_arrays();
  }

  friend class is_before_in_filtration<HasseDiagramOptions>;

  /**
   * A version of clean up structure for the Hasse_diagram_persistence
   * class. Note that before computations of the persistence the structure
   * should be cleaned up. This procedure have to be invoked before computations
   * of persistent homology if some removal operations has been performed on
   * the structure of Hasse diagram.
   **/
  void clean_up_the_structure() {
    Hasse_diagram<HasseDiagramOptions>::clean_up_the_structure();
    this->set_up_the_arrays();
  }

  /**
   * A procedure that need to be called before computations of persistence after some
   * addition (but not removal) operations has been performed on the structure.
   **/
  void set_up_the_arrays() {
    this->cell_associated_to_key_ = std::vector<unsigned>(this->cells.size());
    std::iota(std::begin(this->cell_associated_to_key_), std::end(this->cell_associated_to_key_), 0);
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(this->cell_associated_to_key_.begin(), this->cell_associated_to_key_.end(),
                       is_before_in_filtration<HasseDiagramOptions>(this));
#else
    std::sort(this->cell_associated_to_key_.begin(), this->cell_associated_to_key_.end(),
              is_before_in_filtration<HasseDiagramOptions>(this));
#endif
    this->key_associated_to_cell_ = std::vector<unsigned>(this->cell_associated_to_key_.size());
    for (size_t i = 0; i != this->cell_associated_to_key_.size(); ++i) {
      this->key_associated_to_cell_[this->cell_associated_to_key_[i]] = static_cast<unsigned>(i);
    }
  }

  // From here on we have implementation of methods that are required to use
  // this class with persistent homology engine.

  typedef typename HasseDiagramOptions::Filtration_value Filtration_value;
  typedef unsigned Simplex_key;
  typedef Simplex_key Simplex_handle;

  size_t num_simplices() { return this->cells.size(); }

  Simplex_key key(Simplex_handle sh) {
    if (sh != null_simplex()) {
      return this->key_associated_to_cell_[sh];
    }
    return this->null_key();
  }

  Simplex_key null_key() { return std::numeric_limits<unsigned>::max(); }

  Simplex_handle simplex(Simplex_key key) {
    if (key != null_key()) {
      return this->cell_associated_to_key_[key];
    }
    return null_simplex();
  }

  Simplex_handle null_simplex() { return this->null_key(); }

  Filtration_value filtration(Simplex_handle sh) {
    if (sh == null_simplex()) {
      return std::numeric_limits<Filtration_value>::infinity();
    }
    return this->cells[sh]->get_filtration();
  }

  int dimension(Simplex_handle sh) {
    if (sh == null_simplex()) {
      return -1;
    }
    return this->cells[sh]->dimension();
  }

  int dimension() {
    int top_dimension = 0;
    for (size_t i = 0; i != this->cells.size(); ++i) {
      int dim_of_cell = this->cells[i]->dimension();
      if (top_dimension < dim_of_cell) {
        top_dimension = dim_of_cell;
      }
    }
    return top_dimension;
  }

  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    if (sh == null_simplex()) {
      return std::pair<Simplex_handle, Simplex_handle>(null_simplex(), null_simplex());
    }
    std::vector<std::pair<Cell_type*, typename HasseDiagramOptions::Incidence_type> > boundary = this->cells[sh]->boundaries();
    return std::pair<Simplex_handle, Simplex_handle>(boundary[0].first->position(),
                                                     boundary[1].first->position());
  }

  void assign_key(Simplex_handle sh, Simplex_key key) {
    if (key == null_key()) return;
    this->key_associated_to_cell_[sh] = key;
    this->cell_associated_to_key_[key] = sh;
  }

  // ********************************************************************************************
  //                                    FILTRATION SIMPLEX ITERATOR
  // ********************************************************************************************
  class Filtration_simplex_range;
  class Filtration_simplex_iterator : std::iterator<std::input_iterator_tag, Simplex_handle> {
    // Iterator over all simplices of the complex in the order of the indexing scheme.
    // 'value_type' must be 'Simplex_handle'.
   public:
    Filtration_simplex_iterator(Hasse_diagram_persistence<HasseDiagramOptions>* hd) : hasse_diagram_(hd), position_(0) {}
    Filtration_simplex_iterator() : hasse_diagram_(NULL), position_(0) {}

    Filtration_simplex_iterator operator++() {
      ++this->position_;
      return (*this);
    }

    Filtration_simplex_iterator operator++(int) {
      Filtration_simplex_iterator result = *this;
      ++(*this);
      return result;
    }

    Filtration_simplex_iterator(const Filtration_simplex_iterator& rhs) {
      this->hasse_diagram_ = rhs.hasse_diagram_;
      this->position_ = rhs.position_;
    }

    Filtration_simplex_iterator& operator=(const Filtration_simplex_iterator& rhs) {
      this->hasse_diagram_ = rhs.hasse_diagram_;
      this->position_ = rhs.position_;
      return (*this);
    }

    bool operator==(const Filtration_simplex_iterator& rhs) const {
      return ((this->position_ == rhs.position_) && (this->hasse_diagram_ == rhs.hasse_diagram_));
    }

    bool operator!=(const Filtration_simplex_iterator& rhs) const { return !(*this == rhs); }

    Simplex_handle operator*() { return this->hasse_diagram_->cell_associated_to_key_[this->position_]; }

    friend class Filtration_simplex_range;

   private:
    Hasse_diagram_persistence<HasseDiagramOptions>* hasse_diagram_;
    unsigned position_;
  };

  /**
   * @brief Filtration_simplex_range provides the ranges for Filtration_simplex_iterator.
   **/
  class Filtration_simplex_range {
    // Range over the simplices of the complex in the order of the filtration.
    // .begin() and .end() return type Filtration_simplex_iterator.
   public:
    typedef Filtration_simplex_iterator const_iterator;
    typedef Filtration_simplex_iterator iterator;

    Filtration_simplex_range(Hasse_diagram_persistence<HasseDiagramOptions>* hd) : hasse_diagram_(hd) {}

    Filtration_simplex_iterator begin() { return Filtration_simplex_iterator(this->hasse_diagram_); }

    Filtration_simplex_iterator end() {
      Filtration_simplex_iterator it(this->hasse_diagram_);
      it.position_ = static_cast<unsigned>(this->hasse_diagram_->cell_associated_to_key_.size());
      return it;
    }

   private:
    Hasse_diagram_persistence<HasseDiagramOptions>* hasse_diagram_;
  };

  Filtration_simplex_range filtration_simplex_range() { return Filtration_simplex_range(this); }
  // ********************************************************************************************

  // ********************************************************************************************
  //                                    SKELETON SIMPLEX ITERATOR
  // ********************************************************************************************
  class Skeleton_simplex_range;
  class Skeleton_simplex_iterator : std::iterator<std::input_iterator_tag, Simplex_handle> {
   public:
    Skeleton_simplex_iterator(Hasse_diagram_persistence<HasseDiagramOptions>* hd, int d) : hasse_diagram_(hd), dimension_(d) {
      // find the position of the first cell of a dimension d
      this->position_ = 0;
      while ((this->position_ != this->hasse_diagram_->cells.size()) &&
             (this->hasse_diagram_->cells[this->position_]->dimension() != this->dimension_)) {
        ++this->position_;
      }
    }

    Skeleton_simplex_iterator() : hasse_diagram_(NULL), position_(0), dimension_(0) {}

    Skeleton_simplex_iterator(const Skeleton_simplex_iterator& rhs) {
      this->hasse_diagram_ = rhs.hasse_diagram_;
      this->position_ = rhs.position_;
      this->dimension_ = rhs.dimension_;
    }

    Skeleton_simplex_iterator operator++() {
      ++this->position_;
      while ((this->position_ != this->hasse_diagram_->cells.size()) &&
             (this->hasse_diagram_->cells[this->position_]->dimension() != this->dimension_)) {
        ++this->position_;
      }
      return (*this);
    }

    Skeleton_simplex_iterator operator++(int) {
      Skeleton_simplex_iterator result = *this;
      ++(*this);
      return result;
    }

    Skeleton_simplex_iterator& operator=(const Skeleton_simplex_iterator& rhs) {
      this->hasse_diagram_ = rhs.hasse_diagram_;
      this->position_ = rhs.position_;
      this->dimension_ = rhs.dimension_;
      return (*this);
    }

    bool operator==(const Skeleton_simplex_iterator& rhs) const { return (this->position_ == rhs.position_); }

    bool operator!=(const Skeleton_simplex_iterator& rhs) const { return !(*this == rhs); }

    Simplex_handle operator*() { return this->position_; }

    friend class Skeleton_simplex_range;

   private:
    Hasse_diagram_persistence<HasseDiagramOptions>* hasse_diagram_;
    unsigned position_;
    int dimension_;
  };

  /**
   * @brief Class needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  class Skeleton_simplex_range {
    // Range over the simplices of the complex in the order of the filtration.
    // .begin() and .end() return type Filtration_simplex_iterator.
   public:
    typedef Skeleton_simplex_iterator const_iterator;
    typedef Skeleton_simplex_iterator iterator;

    Skeleton_simplex_range(Hasse_diagram_persistence<HasseDiagramOptions>* hd, unsigned dimension)
        : hasse_diagram_(hd), dimension_(dimension) {}

    Skeleton_simplex_iterator begin() { return Skeleton_simplex_iterator(this->hasse_diagram_, this->dimension_); }

    Skeleton_simplex_iterator end() {
      Skeleton_simplex_iterator it(this->hasse_diagram_, this->dimension_);
      it.position_ = static_cast<unsigned>(this->hasse_diagram_->cells.size());
      return it;
    }

   private:
    Hasse_diagram_persistence<HasseDiagramOptions>* hasse_diagram_;
    unsigned dimension_;
  };

  /**
   * Function needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  Skeleton_simplex_range skeleton_simplex_range(unsigned dimension) { return Skeleton_simplex_range(this, dimension); }
  // ********************************************************************************************

  // ********************************************************************************************
  //                                    BOUNDARY SIMPLEX ITERATOR
  // ********************************************************************************************
  typedef typename std::vector<Simplex_handle>::iterator Boundary_simplex_iterator;
  typedef typename std::vector<Simplex_handle> Boundary_simplex_range;
  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh) {
    return this->cells[sh]->get_list_of_positions_of_boundary_elements();
  }
  //********************************************************************************************

 protected:
  std::vector<unsigned> key_associated_to_cell_;
  std::vector<unsigned> cell_associated_to_key_;
};  // Hasse_diagram

template <typename HasseDiagramOptions>
class is_before_in_filtration {
 public:
  explicit is_before_in_filtration(Hasse_diagram_persistence<HasseDiagramOptions>* hd) : hasse_diagram_(hd) {}

  bool operator()(size_t first, size_t second) const {
    typedef typename HasseDiagramOptions::Filtration_value Filtration_value;
    Filtration_value fil1 = hasse_diagram_->cells[first]->get_filtration();
    Filtration_value fil2 = hasse_diagram_->cells[second]->get_filtration();
    // Compare what can be floating point values
    if (std::fabs(fil1 - fil2) > std::numeric_limits<Filtration_value>::epsilon()) {
      return fil1 < fil2;
    }
    // in this case they are on the same filtration level, so the dimension decides.
    int dim1 = hasse_diagram_->cells[first]->dimension();
    int dim2 = hasse_diagram_->cells[second]->dimension();
    if (dim1 != dim2) {
      return dim1 < dim2;
    }
    // in this case both filtration and dimensions of the considered cells are the same. To have stable sort, we simply
    // compare their positions in the vector. At least those will have to be different.
    return first < second;
  }

 protected:
  Hasse_diagram_persistence<HasseDiagramOptions>* hasse_diagram_;
};

/**
 * This is a function to convert any representation that implements Hasse_complex interface
 * into Hasse_diagram_persistence
 **/
template <typename Complex_type, typename HasseDiagram>
HasseDiagram* convert_to_hasse_diagram_persistence(Complex_type& cmplx) {
  using Cell_type = typename HasseDiagram::Cell_type;
  return new HasseDiagram(convert_to_vector_of_Cell_type<Complex_type, Cell_type>(cmplx));
}  // convert_to_Hasse_diagram

}  // namespace Hasse_diagram
}  // namespace Gudhi

#endif  // HASSE_DIAGRAM_PERSISTENCE_H
