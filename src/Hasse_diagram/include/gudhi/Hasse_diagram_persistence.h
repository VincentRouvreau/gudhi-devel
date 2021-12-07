/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
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
 
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <numeric>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif


#include <gudhi/Hasse_diagram_cell.h>
#include <gudhi/Hasse_diagram.h>



#ifndef HASSE_DIAGRAM_PERSISTENCE_H
#define HASSE_DIAGRAM_PERSISTENCE_H


namespace Gudhi {

namespace Hasse_diagram {
	
template <typename Cell_type> class is_before_in_filtration;


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
template < typename Cell_type >
class Hasse_diagram_persistence : public Hasse_diagram<Cell_type>
{
public:
	/**
	 * Default constructor.
	**/ 
	Hasse_diagram_persistence():Hasse_diagram<Cell_type>(){}
	
	/**
	 * Creating Hasse diagram for persistence computations from a file. 
	 * The file format is the following:
	 * Number of cells
	 * cell dimension
	 * ids of cell boundary elements followed by the incidence coefficient.
	 * the two lines above are repeated for each cell. 
	 * It is assumed that the id of a cell is its position in the file. 	 
	**/ 	
    Hasse_diagram_persistence( const char* filename ):Hasse_diagram<Cell_type>(filename)
    {		
		this->set_up_the_arrays();
	}
    
    /**
	 * Constructor to create a Hasse diagram for persistence computations 
	 * from a vector of cells. It is assumed that all the cells have 
	 * boundaries set up. Setting up the coboundaries will be done in the 
	 * constructor based on the information about boundaries. 
	**/ 
    Hasse_diagram_persistence( const std::vector< Cell_type* >& cells_ ):
    Hasse_diagram<Cell_type>(cells_)
    {		
		this->set_up_the_arrays();
	}		
		
	friend class is_before_in_filtration<Cell_type>;
	
	/**
	 * A version of clean up structure for the Hasse_diagram_persistence
	 * class. Note that before computations of the persistence the structure
	 * should be cleaned up. This procedure have to be invoked before computations
	 * of persistent homology if some removal operations has been performed on
	 * the structure of Hasse diagram.
	**/ 
	void clean_up_the_structure()
	{
		Hasse_diagram<Cell_type>::clean_up_the_structure();
		this->set_up_the_arrays();
	}
	
	/**
	 * A procedure that need to be called before computations of persistence after some
	 * addition (but not removal) operations has been performed on the structure.
	**/ 
	void set_up_the_arrays();

    
	//From here on we have implementation of methods that are required to use
	//this class with persistent homology engine.
	
	typedef typename Cell_type::Filtration_type Filtration_value;
    typedef unsigned Simplex_key;
    typedef Simplex_key Simplex_handle;
    
    size_t num_simplices() 
    {
		return this->cells.size();
    }
    
	
	Simplex_key key(Simplex_handle sh) 
	{
		if ( sh != null_simplex() )
		{
			return this->key_associated_to_cell[sh];
		}
		return this->null_key();
	}

	Simplex_key null_key() 
	{
		return std::numeric_limits<unsigned>::max();
	}

	Simplex_handle simplex(Simplex_key key) 
	{
		if (key != null_key()) 
		{
			return this->cell_associated_to_key[key];
		}
		return null_simplex();
	}

	Simplex_handle null_simplex() 
	{
		return this->null_key();
	}

	Filtration_value filtration(Simplex_handle sh) 
	{
		if (sh == null_simplex()) 
		{
		     return std::numeric_limits<Filtration_value>::infinity();
		}
		return this->cells[ sh ]->get_filtration();
	}

	int dimension(Simplex_handle sh) 
	{
		if (sh == null_simplex()) 
		{
		  return -1;
		}
		return this->cells[sh]->get_dimension();
	}

	int dimension() 
	{
		int top_dimension = 0;
		for ( size_t i = 0 ; i != this->cells.size() ; ++i )
		{
			int dim_of_cell = this->cells[i]->get_dimension();
			if ( top_dimension < dim_of_cell )
			{
				top_dimension = dim_of_cell;
			}
		}
		return top_dimension;
	}

	std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) 
	{
		if (sh == null_simplex())
		{			
			return std::pair<Simplex_handle, Simplex_handle>(null_simplex(),null_simplex());
		}
		std::vector< std::pair<Cell_type*,typename Cell_type::Incidence_type> > boundary = 
		this->cells[sh]->get_boundary();		
		return std::pair<Simplex_handle, Simplex_handle>(
							boundary[0].first->get_position(), 
							boundary[1].first->get_position()
														);																											
	}

	void assign_key(Simplex_handle sh, Simplex_key key) 
	{	
		if (key == null_key()) return;
		this->key_associated_to_cell[sh] = key;
		this->cell_associated_to_key[key] = sh;
	}

		
	
	//********************************************************************************************
	//									FILTRATION SIMPLEX ITERATOR
	//********************************************************************************************
	class Filtration_simplex_range;
	class Filtration_simplex_iterator : std::iterator<std::input_iterator_tag, Simplex_handle> 
	{
		  // Iterator over all simplices of the complex in the order of the indexing scheme.
		  // 'value_type' must be 'Simplex_handle'.
		 public:
		  Filtration_simplex_iterator(Hasse_diagram_persistence<Cell_type>* hd) : hasse_diagram_(hd), position_(0) {}
		  Filtration_simplex_iterator() : hasse_diagram_(NULL), position_(0) {}

		  Filtration_simplex_iterator operator++() 
		  {		
			++this->position_;
		   return (*this);
		  } 

		  Filtration_simplex_iterator operator++(int) 
		  {
			Filtration_simplex_iterator result = *this;
			++(*this);
			return result;
		  }
		  
		  Filtration_simplex_iterator( const Filtration_simplex_iterator& rhs )
		  {
			  this->hasse_diagram_ = rhs.hasse_diagram_;
			  this->position_ = rhs.position_;
		  }

		  Filtration_simplex_iterator& operator=(const Filtration_simplex_iterator& rhs) {
			this->hasse_diagram_ = rhs.hasse_diagram_;
			this->position_ = rhs.position_;
			return (*this);
		  }

		  bool operator==(const Filtration_simplex_iterator& rhs) const 
		  {	
			return ((this->position_ == rhs.position_)&&(this->hasse_diagram_ == rhs.hasse_diagram_));
		  }

		  bool operator!=(const Filtration_simplex_iterator& rhs) const 
		  {
			return !(*this == rhs);
		  }

		  Simplex_handle operator*() 
		  {			
			 return this->hasse_diagram_->cell_associated_to_key[this->position_];
		  }

		  friend class Filtration_simplex_range;

	private:
		Hasse_diagram_persistence<Cell_type>* hasse_diagram_;
		unsigned position_;
	};

	/**
	* @brief Filtration_simplex_range provides the ranges for Filtration_simplex_iterator.
	**/
	class Filtration_simplex_range 
	{
	// Range over the simplices of the complex in the order of the filtration.
	// .begin() and .end() return type Filtration_simplex_iterator.
	public:
		typedef Filtration_simplex_iterator const_iterator;
		typedef Filtration_simplex_iterator iterator;

		Filtration_simplex_range(Hasse_diagram_persistence<Cell_type>* hd):hasse_diagram_(hd) {}

		Filtration_simplex_iterator begin() 
		{		  
			return Filtration_simplex_iterator(this->hasse_diagram_);
		}

		Filtration_simplex_iterator end() 
		{		 
			Filtration_simplex_iterator it(this->hasse_diagram_);
			it.position_ = static_cast<unsigned>( this->hasse_diagram_->cell_associated_to_key.size() );
			return it;
		}

	private:
		Hasse_diagram_persistence<Cell_type>* hasse_diagram_;
	};
   
    Filtration_simplex_range filtration_simplex_range() 
    {    
		return Filtration_simplex_range(this);
    }
    //********************************************************************************************



	//********************************************************************************************
	//									SKELETON SIMPLEX ITERATOR
	//********************************************************************************************
	class Skeleton_simplex_range;
	class Skeleton_simplex_iterator : std::iterator<std::input_iterator_tag, Simplex_handle> 
	{	
	public:
	Skeleton_simplex_iterator(Hasse_diagram_persistence<Cell_type>* hd, int d) : hasse_diagram_(hd), dimension_(d) 
	{	 
	  //find the position of the first cell of a dimension d
	  this->position_ = 0;
	  while ((this->position_ != this->hasse_diagram_->cells.size()) && 
			 (this->hasse_diagram_->cells[this->position_]->get_dimension() != this->dimension_)) 
			 {
				++this->position_;
			  }
	}

	Skeleton_simplex_iterator() : hasse_diagram_(NULL), position_(0), dimension_(0) {}
	
	Skeleton_simplex_iterator( const Skeleton_simplex_iterator& rhs )
	{
		this->hasse_diagram_ = rhs.hasse_diagram_;
		this->position_ = rhs.position_;
	    this->dimension_ = rhs.dimension_;
	}

	Skeleton_simplex_iterator operator++() 
	{	  
	  ++this->position_;
	  while ( ( this->position_ != this->hasse_diagram_->cells.size() ) &&
			 (this->hasse_diagram_->cells[this->position_]->get_dimension() != this->dimension_)) 
	  {
		  ++this->position_;
	  }
	  return (*this);
	}

	Skeleton_simplex_iterator operator++(int) 
	{
		Skeleton_simplex_iterator result = *this;
		++(*this);
		return result;
	}

	Skeleton_simplex_iterator& operator=(const Skeleton_simplex_iterator& rhs) 
	{
	  this->hasse_diagram_ = rhs.hasse_diagram_;
	  this->position_ = rhs.position_;
	  this->dimension_ = rhs.dimension_;
	  return (*this);
	}

	bool operator==(const Skeleton_simplex_iterator& rhs) const 
	{
	  return (this->position_ == rhs.position_);
	}

	bool operator!=(const Skeleton_simplex_iterator& rhs) const 
	{	
	  return !(*this == rhs);
	}

	Simplex_handle operator*() 
	{	
	  return this->position_;
	}

	friend class Skeleton_simplex_range;

	private:
	Hasse_diagram_persistence<Cell_type>* hasse_diagram_;
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

	Skeleton_simplex_range(Hasse_diagram_persistence<Cell_type>* hd, unsigned dimension) : hasse_diagram_(hd), dimension_(dimension) {}

	Skeleton_simplex_iterator begin() 
	{	
	  return Skeleton_simplex_iterator(this->hasse_diagram_, this->dimension_);
	}

	Skeleton_simplex_iterator end() 
	{	
	  Skeleton_simplex_iterator it(this->hasse_diagram_, this->dimension_);
	  it.position_ = static_cast<unsigned>(this->hasse_diagram_->cells.size());
	  return it;
	}

	private:
	Hasse_diagram_persistence<Cell_type>* hasse_diagram_;
	unsigned dimension_;
	};

	/**
	* Function needed for compatibility with Gudhi. Not useful for other purposes.
	**/
	Skeleton_simplex_range skeleton_simplex_range(unsigned dimension) 
	{	
		return Skeleton_simplex_range(this, dimension);
	}
	//********************************************************************************************	


	//********************************************************************************************
	//									BOUNDARY SIMPLEX ITERATOR
	//********************************************************************************************
	typedef typename std::vector<Simplex_handle>::iterator Boundary_simplex_iterator;
    typedef typename std::vector<Simplex_handle> Boundary_simplex_range;
	Boundary_simplex_range boundary_simplex_range(Simplex_handle sh) 
	{ 
		return this->cells[sh]->get_list_of_positions_of_boundary_elements(); 
	}
	//********************************************************************************************	
	
protected:	
	  std::vector< unsigned > key_associated_to_cell;
      std::vector< unsigned > cell_associated_to_key;   
};//Hasse_diagram



template <typename Cell_type>
class is_before_in_filtration {
 public:
  explicit is_before_in_filtration(Hasse_diagram_persistence<Cell_type>* hd) : hasse_diagram_(hd) {}

  bool operator()(size_t first , size_t second) const 
  {    
    typedef typename Cell_type::Filtration_type Filtration_value;
    Filtration_value fil1 = hasse_diagram_->cells[first]->get_filtration();
    Filtration_value fil2 = hasse_diagram_->cells[second]->get_filtration();
    if (*(double*)&fil1 != *(double*)&fil2)//strange cast to make the compiler shut up about comparing of doubles.
     {
      return fil1 < fil2;
    }
    // in this case they are on the same filtration level, so the dimension decide.
    unsigned dim1 = hasse_diagram_->cells[first]->get_dimension();
    unsigned dim2 = hasse_diagram_->cells[second]->get_dimension();
    if (dim1 != dim2) {
      return dim1 < dim2;
    }
    // in this case both filtration and dimensions of the considered cells are the same. To have stable sort, we simply
    // compare their positions in the vector. At least those will have to be different.
    return first < second;
  }

 protected:
  Hasse_diagram_persistence<Cell_type>* hasse_diagram_;
};



template < typename Cell_type >
void Hasse_diagram_persistence<Cell_type>::set_up_the_arrays()
{		
	this->cell_associated_to_key = std::vector<unsigned>( this->cells.size() );
	std::iota (std::begin(this->cell_associated_to_key), std::end(this->cell_associated_to_key), 0); 	
	#ifdef GUDHI_USE_TBB
	  tbb::parallel_sort(this->cell_associated_to_key.begin(), this->cell_associated_to_key.end(),
						 is_before_in_filtration<Cell_type>(this));
	#else
	  std::sort(this->cell_associated_to_key.begin(), this->cell_associated_to_key.end(), is_before_in_filtration<Cell_type>(this));
	#endif	
	this->key_associated_to_cell = std::vector<unsigned>( this->cell_associated_to_key.size() );	
	for (size_t i = 0; i != this->cell_associated_to_key.size(); ++i) 
	{		
		this->key_associated_to_cell[this->cell_associated_to_key[i]] = static_cast<unsigned>(i);
    } 
    
    //Just for debugging purposes, remove later. 
    //std::cout << "cell_associated_to_key, just for degugging. : \n";
    //for ( size_t i = 0 ; i != this->cell_associated_to_key.size() ; ++i )          
    //{
	//	std::cout << this->cell_associated_to_key[i] << " ";
	//}
	//std::cout << std::endl;
}//Cell_type



/**
 * This is a function to convert any representation that implements Hasse_complex interface
 * into Hasse_diagram_persistence
**/ 
template <typename Complex_type , typename Cell_type>
Hasse_diagram_persistence<Cell_type>* convert_to_Hasse_diagram_persistence( Complex_type& cmplx )
{	
	return new Hasse_diagram_persistence<Cell_type>
	( 
	convert_to_vector_of_Cell_type< Complex_type , Cell_type >(cmplx) 
	);	
}//convert_to_Hasse_diagram


}//namespace Hasse_diagram_persistence
}//namespace Gudhi

#endif //HASSE_DIAGRAM_PERSISTENCE_H
