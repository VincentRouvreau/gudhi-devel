/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * A class to implement general cell type.
**/ 
class Cell_type
{
public:
    /**
     * Must be numerical, assignable, from a ring or field.
    **/ 
	typedef unspecified Incidence_type;
	
	/**
	 * Must provide an oder, and implement << operator.
	**/ 
	typedef unspecified Filtration_value;
	
    /**
     * Default constructor.
    **/
	Hasse_diagram_cell();

	/**
     * Procedure to obtain the boundary of a given cell. It returns a 
     * vector of pairs of pointers to boundary elements and incidence
     * coefficients.
    **/
	inline std::vector< std::pair<Cell_type*,Incidence_type> >& boundaries();
	inline std::vector< std::pair<Cell_type*,Incidence_type> > boundaries() const;

	/**
     * Cell dimension accessors.
    **/
	inline int& dimension();
	inline int dimension() const;

	/**
	 * Position of a cell in a structure accessors. Class Cell_type is required
	 * to be able to store a its position in a Hasse diagram. 
	**/
	inline unsigned& position();
	inline unsigned position() const;
	
	/**
	 * Cell filtration accessors.
	**/
	inline void set_filtration(const Filtration_value& filt);
	inline Filtration_value get_filtration() const;

	/**
	 * Implementation of Hasse diagram assume that the Cell_type allows removal
	 * and addition of cells. This procedure is used to check if the cell has 
	 * not been removed.  
	**/
	inline bool deleted();

	/**
	 * Have to be a friend to a class Hasse_diagram implemented in Hasse_diagram.h
	**/ 
	template < typename Cell_type >
	friend class Hasse_diagram;
	
	/**
	 * Have to be a friend with a class is_before_in_filtration implemented in 
	 * Hasse_diagram_persistence.h. It is required to sort the cells according
	 * filtration in the persistence algorithm.
	**/ 
	template < typename Cell_type >
	friend class is_before_in_filtration;
	
	/**
	 * Function convert_to_vector_of_Cell_type is used to cinvert any class that 
	 * implement Hasse_complex interface into Hasse diagram.
	**/ 
	template <typename Complex_type , typename Cell_type >  
	friend std::vector<Cell_type*> convert_to_vector_of_Cell_type( Complex_type& cmplx );

	/**
	 * Procedure to remove deleted boundary and coboundary elements from the
	 * vectors of boundary and coboundary elements of this cell. Should use
	 * deleted() method to check if boundary or coboudnary elements has been deleted.
	**/
	void remove_deleted_elements_from_boundary_and_coboundary();


	/**
	 * Writing to a stream operator.
	**/
	friend std::ostream& operator<<( std::ostream& out, const Cell_type& c );

	
	/**
	 * Procedure that return vector of positios of boundary elements of a given cell.
	**/  	
	inline std::vector< unsigned > get_list_of_positions_of_boundary_elements() const;

	
	/**
	 * Function that display a string being a signature of a structure. 
	 * Used mainly for debugging purposes, but required for the tests. 
	**/ 
	std::string full_signature_of_the_structure();
	
	
	/**
	 * Class Cell_type is required to store boundary elements in the following vector.  
	**/ 
	std::vector< std::pair<Cell_type*,Incidence_type> > boundary;


};//Cell_type
