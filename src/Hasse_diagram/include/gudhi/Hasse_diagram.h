/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#ifndef HASSE_DIAGRAM_H
#define HASSE_DIAGRAM_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <utility>  // for std::make_pair
#include <fstream>  // for std::ifstream

#include <gudhi/Hasse_diagram_cell.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {

namespace Hasse_diagram {

template <typename Cell_type>
class is_before_in_dimension;

/** Model of HasseDiagramOptions.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Hasse_diagram_options_full_featured {
  using Filtration_value = double;
  using Simplex_key = std::uint32_t;
  // With Incidence_type = std::uint8_t, the method full_signature_of_the_structure behaves strangely as std::string(std::uint8_t)
  using Incidence_type = std::uint16_t;
  using Additional_information = bool;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool store_additional_information = true;
};

/**
 * \class Hasse_diagram
 * \brief Data structure to store Hasse diagrams.
 *
 * \ingroup Hasse_diagram
 *
 * \details
 * This is a data structure to store a Hasse diagrams. It allows to store a general
 * chain complexes in a form of Hasse diagrams. It implements insertion and
 * removal of cells. It allows to store and read Hasse diagrams from files.
 *
 * Please note that this class is not suitable to use with Gudhi engine to compute
 * persistent homology. For that purpose, please use the derived class provided in
 * Hasse_diagram_persistence.h
 *
 * Please refer to \ref Hasse_diagram for examples.
 *
 * The complex is a template class requiring the following parameters:
 * Cell_type - a parameter describing a cell of Hasse diagram. Please refer to Hasse_diagram_cell.h for further details.
 *
 */
template <typename HasseDiagramOptions = Hasse_diagram_options_full_featured>
class Hasse_diagram {
 public:
  // + Manage what is stored or not in a cell [
  using Options = HasseDiagramOptions;
  using Incidence_type = typename Options::Incidence_type;
  /** \brief Type for the value of the filtration function.
   *
   * Must be comparable with <. */
  using Filtration_value = typename Options::Filtration_value;
  /** \brief Key associated to each simplex.
   *
   * Must be an integer type. */
  using Simplex_key = typename Options::Simplex_key;
  /** \brief Type for some additionnal data.
   */
  using Additional_information = typename Options::Additional_information;

  struct Key_simplex_base_real {
    Key_simplex_base_real() : key_(-1) {}
    void assign_key(Simplex_key k) { key_ = k; }
    Simplex_key key() const { return key_; }
   private:
    Simplex_key key_;
  };
  struct Key_simplex_base_dummy {
    Key_simplex_base_dummy() {}
    // Undefined so it will not link
    void assign_key(Simplex_key);
    Simplex_key key() const;
  };
  typedef typename std::conditional<Options::store_key, Key_simplex_base_real, Key_simplex_base_dummy>::type
      Key_simplex_base;

  struct Filtration_simplex_base_real {
    Filtration_simplex_base_real() : filt_(0) {}
    void assign_filtration(Filtration_value f) { filt_ = f;}
    Filtration_value filtration() const { return filt_;  std::cout << filt_ << std::endl;}
   private:
    Filtration_value filt_;
  };
  struct Filtration_simplex_base_dummy {
    Filtration_simplex_base_dummy() {}
    void assign_filtration(Filtration_value GUDHI_CHECK_code(f)) {
      GUDHI_CHECK(f == 0, "Filtration value specified for a complex that does not store them");
    }
    Filtration_value filtration() const { return 0; }
  };
  typedef typename std::conditional<Options::store_filtration, Filtration_simplex_base_real,
    Filtration_simplex_base_dummy>::type Filtration_simplex_base;

  struct Additional_information_simplex_base_real {
    Additional_information_simplex_base_real() : ai_(0) {}
    void assign_additional_information(Additional_information f) { ai_ = f; }
    Additional_information additional_information() const { return ai_; }
   private:
    Additional_information ai_;
  };
  struct Additional_information_simplex_base_dummy {
    Additional_information_simplex_base_dummy() {}
    void assign_additional_information(Filtration_value GUDHI_CHECK_code(f)) {
      GUDHI_CHECK(f == 0, "Additional information value specified for a complex that does not store them");
    }
    Additional_information additional_information() const { return 0; }
  };
  typedef typename std::conditional<Options::store_additional_information, Additional_information_simplex_base_real,
    Additional_information_simplex_base_dummy>::type Additional_information_simplex_base;

  /* Type of Cell in the Hasse diagram. */
  using Cell_type = Hasse_diagram_cell<Hasse_diagram>;

  // - Manage what is stored or not in a cell ]
  using Cell_range = std::vector<Cell_type*>;

  /**
   * Default constructor.
   **/
  Hasse_diagram() {}

  /**
   * Creating Hasse diagram from a file. The file format is the following:
   * Number of cells
   * cell dimension
   * ids of cell boundary elements followed by the incidence coefficient.
   * the two lines above are repeated for each cell.
   * It is assumed that the id of a cell is its position in the file.
   **/
  Hasse_diagram(const char* filename);

  /**
   * Constructor to create a Hasse diagram from a vector of cells. It is assumed
   * that all the cells have boundaries set up. Setting up the coboundaries will
   * be done in the constructor based on the information about boundaries.
   **/
  Hasse_diagram(const Cell_range& cells_) : cells(cells_), number_of_deleted_cells(0) {
    this->set_up_positions();
    this->set_up_coboundaries();
  }

  /**
   * After many operation of deleting cells, this->cells vector may became
   * very fragmented. Also, the complexity of operation using all the iterators
   * depends on the actual size of a structure (where the deleted elements are still
   * stored. This procedure remove permanently all the deleted elements. Ideally,
   * it should be initialized when the proportion of deleted elements is larger than a
   * predefined constant.
   **/
  void clean_up_the_structure() {
    if (this->number_of_deleted_cells == 0) return;

#ifdef DEBUG_TRACES
    std::clog << "Calling clean_up_the_structure() procedure. \n";
#endif
    // count the number of not deleted cells:
    size_t number_of_non_deleted_cells = this->cells.size() - this->number_of_deleted_cells;

    // create a new vector to store the undeleted cells:
    std::vector<Cell_type*> new_cells;
    new_cells.reserve(number_of_non_deleted_cells);
    // fill the new vector in and adjust the new positions.
    // In the same time make sure that the boundary and coboundary vectors
    // in every cell are valid.
    size_t counter = 0;
    for (size_t i = 0; i != this->cells.size(); ++i) {
      if (!this->cells[i]->deleted()) {
        new_cells.push_back(this->cells[i]);
        this->cells[i]->position_ = static_cast<unsigned>(counter);
        this->cells[i]->remove_deleted_elements_from_boundary_and_coboundary();
        ++counter;
      } else {
        delete this->cells[i];
      }
    }
    this->cells.swap(new_cells);
    this->number_of_deleted_cells = 0;
#ifdef DEBUG_TRACES
    std::clog << "Done with the clean_up_the_structure() procedure. \n";
#endif
  }

  /**
   * Procedure that allow to add a cell into the structure. This procedure
   * automatically fill in coboundaries of boundary elements, so do not
   * duplicate it.
   **/
  void add_cell(Cell_type* cell) {
    cell->position_ = static_cast<unsigned>(this->cells.size());
    this->cells.push_back(cell);
    // we still need to check if cobounadies of boundary elements of this
    // cell are set up in the correct way:
    for (size_t bd = 0; bd != cell->boundary_.size(); ++bd) {
      cell->boundary_[bd].first->coboundary_.push_back(std::make_pair(cell, cell->boundary_[bd].second));
    }
  }

  /**
   * Procedure that allow to remove a cell into the structure.
   **/
  void remove_cell(Cell_type* cell) {
    // if the flag enable_checking_validity_of_complex is set to true,
    // we will check if the cell that is to be deleted do not have
    // a non deleted cell in the coboundary and if this is the case, we
    // will print out the warning, since this can potentially be an error.
    if (enable_checking_validity_of_complex) {
      for (size_t cbd = 0; cbd != cell->coboundary_.size(); ++cbd) {
        if (!cell->coboundary_[cbd].first->deleted()) {
          std::cout << "Warning, you are deleting cell which have non-deleted cells in the coboundary. This may lead "
                       "inconsistencies in the data structure.\n";
          break;
        }
      }
    }

    cell->delete_cell();
    this->number_of_deleted_cells++;

    // in case the structure gets too fragmented, we are calling the
    // to clean it up.
    if (this->number_of_deleted_cells / (static_cast<double>(this->cells.size())) >
        this->proportion_of_removed_cells_that_triggers_reorganization_of_structure) {
      this->clean_up_the_structure();
    }
  }

  /**
   * A procedure writng Hasse diagram to file. The Hasse diagram can be later
   * reconstructed using Hasse_diagram( const char* filename ) constructor.
   **/
  void write_to_file(const char* filename);

  /**
   * Writing to a stream operator.
   **/
  friend std::ostream& operator<<(std::ostream& out, const Hasse_diagram& c) {
    for (size_t i = 0; i != c.cells.size(); ++i) {
      // if the cell is deleted, ignore it.
      if (c.cells[i]->deleted()) continue;
      out << *(c.cells[i]);
    }
    return out;
  }

  friend class is_before_in_dimension<Cell_type>;

  /**
   * A basic iterator that iterate through all the cells in the structure. It is the
   * user's responsibility to check if the cell is deleted or not.
   **/
  typedef typename std::vector<Cell_type*>::iterator Simple_all_cells_iterator;
  typedef typename std::vector<Cell_type*> Simple_all_cells_iterator_range;
  Simple_all_cells_iterator_range simple_all_cells_iterator_range() { return this->cells; }

  /**
   * Procedure that retuns a cell in the position pos in the vector of cells.
   * Note that this cell will change after calling clean_up_the_structure()
   * procedure.
   **/
  inline Cell_type* give_me_cell_at_position(size_t pos) {
    if (pos < this->cells.size()) {
      return this->cells[pos];
    } else {
      std::cerr << "Wrong position of a cell in the give_me_cell_at_position function.\n";
      throw "Wrong position of a cell in the give_me_cell_at_position function.\n";
    }
  }

  /**
   * Function that display a string being a signature of a structure.
   * Used mainly for debugging purposes.
   **/
  std::string full_signature_of_the_structure() {
    std::string result;
    for (size_t i = 0; i != this->cells.size(); ++i) {
      result += this->cells[i]->full_signature_of_the_structure();
    }
    return result;
  }

 protected:
  Cell_range cells;

  // to check how fragmented the data structure is (as a result of removing cells).
  size_t number_of_deleted_cells;

  /**
   * This procedure assumes that the boundaries are already set up for all
   * the cells, and set up the coboundaries based on them.
   **/
  void set_up_coboundaries();

  /**
   * When cells are not constructed by the class, but given from elsewhere,
   * they may not have the positions being set up. This procedure sets them up.
   **/
  void set_up_positions();

  static double proportion_of_removed_cells_that_triggers_reorganization_of_structure;

  /**
   * This variable indicate if a warning should be given anytime a cell is
   * deleted that have nondeleted cell in the coboundary.
   **/
  static bool enable_checking_validity_of_complex;
};  // Hasse_diagram

template <typename Cell_type>
bool Hasse_diagram<Cell_type>::enable_checking_validity_of_complex = true;

template <typename Cell_type>
double Hasse_diagram<Cell_type>::proportion_of_removed_cells_that_triggers_reorganization_of_structure = 0.5;

template <typename Cell_type>
Hasse_diagram<Cell_type>::Hasse_diagram(const char* filename) {
  // We assume that the cells in the file are enumerated in increasing order.
  // The idea of the i-th cell in the file is by default i (starting from zero).
  // Moreover, the cells are provided from low dimensional to high dimensiona.
  // By doing so, we know that when constructing a given cell, all its boundary
  // has already been constructed.
  // Here is the format of a file:
  // Number of cells
  // cell dimension
  // ids of cell boundary elements followed by the incidence coefficient.
  // Note that coboundary vector will be computed based on boundary vector.
  std::string line;

  this->number_of_deleted_cells = 0;

  std::ifstream in(filename);
  if (!in.good()) {
    std::cout << "The file do not exist, program will now terminate.\n";
    throw "The file do not exist, program will now terminate.\n";
  }

  std::getline(in, line);
  while (line[0] == '#') {
    std::getline(in, line);
  }
  std::stringstream iss(line);

  unsigned number_of_cells;
  iss >> number_of_cells;
  this->cells.reserve(number_of_cells);

  // create all the cells:
  for (size_t i = 0; i != number_of_cells; ++i) {
    this->cells.push_back(new Cell_type());
  }

  std::getline(in, line);

#ifdef DEBUG_TRACES
  std::clog << "Number of cells : " << number_of_cells << std::endl;
#endif

  size_t size_of_last_boundary = 10;  // to initially reserve a vector for bounary elements.
  for (size_t i = 0; i != number_of_cells; ++i) {
    Cell_type* new_cell = this->cells[i];
    while (line[0] == '#') {
      std::getline(in, line);
    }

    iss.str("");
    iss.clear();
    iss << line;
    iss >> new_cell->position_ >> new_cell->dimension();

#ifdef DEBUG_TRACES
    std::clog << "Position and dimension of the cell : " << new_cell->position_ << " , " << new_cell->dimension()
              << std::endl;
#endif

    if (new_cell->position_ != i) {
      std::cerr << "Wrong numeration of cells in the file. Cell number : " << i
                << " is marked as : " << new_cell->position_ << " in the file." << std::endl;
      throw "Wrong numeration of cells in the file.";
    }
    if (iss.good()) {
      // in this case we still have a filtration value to be read
      // from the file.
      Filtration_value filt;
      iss >> filt;
      new_cell->set_filtration(filt);
#ifdef DEBUG_TRACES
      std::clog << "Filtration of the cell : " << new_cell->get_filtration() << std::endl;
#endif
    } else {
      new_cell->set_filtration(0.);
    }

    std::getline(in, line);
    while (line[0] == '#') {
      std::getline(in, line);
    }

    iss.str("");
    iss.clear();
    iss << line;
    unsigned cell_id;
    typename Cell_type::Incidence_type incidence_coef;
    std::vector<std::pair<unsigned, typename Cell_type::Incidence_type> > bdry;
    bdry.reserve(size_of_last_boundary);
#ifdef DEBUG_TRACES
    std::clog << "Here are the boundary elements of the cell.\n";
#endif
    while (iss.good()) {
      iss >> cell_id;
      if (!iss.good()) continue;

      if (cell_id >= i) {
        std::cerr << "Wrong format of a file. THe cell number : " << i
                  << " contain in a boundary a cell that has not been introduced yet.\n";
      }
      iss >> incidence_coef;
#ifdef DEBUG_TRACES
      std::clog << "( " << cell_id << " , " << incidence_coef << " ), ";
#endif
      bdry.push_back(std::pair<unsigned, typename Cell_type::Incidence_type>(cell_id, incidence_coef));
    }

    size_of_last_boundary = bdry.size();
    new_cell->boundary_.reserve(size_of_last_boundary);
    for (size_t bd = 0; bd != size_of_last_boundary; ++bd) {
      new_cell->boundary_.push_back(std::make_pair(this->cells[bdry[bd].first], bdry[bd].second));
    }
#ifdef DEBUG_TRACES
    std::clog << "new_cell->boundary_.size() : " << new_cell->boundary_.size() << std::endl;
    std::clog << "Done with this cell. \n";
#endif

    std::getline(in, line);
    while (line[0] == '#') {
      std::getline(in, line);
    }
  }
  // now once the boundaries are set, we are to set up the coboundaries.
  this->set_up_coboundaries();
}

template <typename Cell_type>
void Hasse_diagram<Cell_type>::set_up_coboundaries() {
  // first we check the number of coboundary elements for each cell:
  size_t number_of_cells = this->cells.size();
  std::vector<unsigned> sizes_of_coboundary(number_of_cells, 0);
  for (size_t i = 0; i != number_of_cells; ++i) {
    std::vector<std::pair<Cell_type*, typename Cell_type::Incidence_type> > bdry = this->cells[i]->boundaries();
    for (size_t bd = 0; bd != bdry.size(); ++bd) {
      sizes_of_coboundary[bdry[bd].first->position()]++;
    }
  }

  // now we reserve the space for all coboundaries
  for (size_t i = 0; i != number_of_cells; ++i) {
    this->cells[i]->coboundary_.reserve(sizes_of_coboundary[i]);
  }

  // and now we set up the coboundaries.
  for (size_t i = 0; i != number_of_cells; ++i) {
    for (size_t bd = 0; bd != this->cells[i]->boundary_.size(); ++bd) {
      this->cells[this->cells[i]->boundary_[bd].first->position_]->coboundary_.push_back(
          std::make_pair(this->cells[i], this->cells[i]->boundary_[bd].second));
    }
  }
}

template <typename Cell_type>
void Hasse_diagram<Cell_type>::set_up_positions() {
  for (unsigned pos = 0; pos != this->cells.size(); ++pos) {
    this->cells[pos]->position() = pos;
  }
}

template <typename Cell_type>
void Hasse_diagram<Cell_type>::write_to_file(const char* filename) {
  std::ofstream out(filename);
  // If there are any deleted cells, then we need to clean up the structure
  // first before writing it to a file. The reason for that is because the
  // file format assumes a continuous enumeration of cells (from zero to the
  // number od cells). This is not satisfied if we have in the structure any
  // deleted elements.
  if (this->number_of_deleted_cells != 0) {
    this->clean_up_the_structure();
  }
  // now we can write to a file the number of (non deleted) cells in the structure.
  out << this->cells.size() << std::endl;
  // and then the rest of the Hasse diagram.
  out << *this;
  out.close();
}  // template < typename Cell_type >

/**
 * This is a function that take any representation that implements Hasse_complex
 * interface and return vector of Cell_type* based on it. It is used to construct
 * objects of class Hasse_diagram and Hasse_diagram_persistence
 **/

template <typename Complex_type, typename Cell_type>
std::vector<Cell_type*> convert_to_vector_of_Cell_type(Complex_type& cmplx) {
#ifdef DEBUG_TRACES
  std::clog << "cmplx.num_simplices() : " << cmplx.num_simplices() << std::endl;
#endif

  // create vector of cells of suitable length:
  std::vector<Cell_type*> cells_of_Hasse_diag(cmplx.num_simplices());
  for (size_t i = 0; i != cmplx.num_simplices(); ++i) {
    cells_of_Hasse_diag[i] = new Cell_type();
  }

  // First we need to assign keys. It is not neccessary for cubical complexes,
  // but in simplex tree this is not done by default.
  std::vector<typename Complex_type::Simplex_key> boundary;
  typename Complex_type::Filtration_simplex_range range = cmplx.filtration_simplex_range();
  size_t pos = 0;
  for (typename Complex_type::Filtration_simplex_iterator it = range.begin(); it != range.end(); ++it) {
    cmplx.assign_key(*it, pos);
    ++pos;
  }

  size_t counter = 0;
  for (typename Complex_type::Filtration_simplex_iterator it = range.begin(); it != range.end(); ++it) {
    Cell_type* this_cell = cells_of_Hasse_diag[counter];

    this_cell->dimension() = static_cast<int>(cmplx.dimension(*it));
    this_cell->set_filtration(static_cast<typename Cell_type::Filtration_value>(cmplx.filtration(*it)));
#ifdef DEBUG_TRACES
    std::clog << "Cell [" << counter << "] of dimension " << this_cell->dimension()
              << " and filtration " << this_cell->get_filtration();
#endif

    // get the boundary:
    boundary.clear();
    boundary.reserve(10);
    typename Complex_type::Boundary_simplex_range bd_range = cmplx.boundary_simplex_range(*it);
    for (typename Complex_type::Boundary_simplex_iterator bd = bd_range.begin(); bd != bd_range.end(); ++bd) {
      boundary.push_back(cmplx.key(*bd));
    }

#ifdef DEBUG_TRACES
    std::clog << " - boundaries: ";
    for (const auto& bnd : boundary) std::clog << bnd << " ";
    std::clog << std::endl;
#endif

    // get the boundary in the Hasse diagram format:
    this_cell->boundary_.reserve(boundary.size());
    typename Cell_type::Incidence_type incidence = 1;
    for (size_t bd = 0; bd != boundary.size(); ++bd) {
      this_cell->boundary_.push_back(std::make_pair(cells_of_Hasse_diag[boundary[bd]], incidence));
      incidence *= -1;
    }
    ++counter;
  }
  return cells_of_Hasse_diag;
}  // convert_to_vector_of_Cell_type

/**
 * This is a function to convert any representation that implements Hasse_complex interface
 * into Hasse diagram
 **/
template <typename Complex_type, typename HasseDiagramOptions>
Hasse_diagram<HasseDiagramOptions>* convert_to_Hasse_diagram(Complex_type& cmplx) {
  return new Hasse_diagram<HasseDiagramOptions>(convert_to_vector_of_Cell_type(cmplx));
}  // convert_to_Hasse_diagram

}  // namespace Hasse_diagram
}  // namespace Gudhi

#endif  // HASSE_DIAGRAM_H
