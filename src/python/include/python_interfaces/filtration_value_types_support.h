/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <string>
#include <type_traits>  // for std::is_same_v

#include <boost/mpl/list.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/placeholders.hpp>

// Developper notice:
// If you need to add a new type to Filtration_value_types_supported, you will also have to add a condition to
// get_string_filtration_type() for the python class name to be suffixed with the type. The class names must be
// different (for example, _Bitmap_cubical_complex_interface and _Bitmap_cubical_complex_interface_float32)
template<class Filtration_value>
constexpr std::string get_string_filtration_type()
{
  if constexpr (std::is_same_v<Filtration_value, float>)
    return std::string("_float32");
  // for backward compatibility
  return std::string("");
}

struct Filtration_value_types_supported {
  using List = boost::mpl::list<boost::mpl::identity<float>,
                                boost::mpl::identity<double>>;
};


// Data_structure must accept Filtration_value as a template parameter (Data_structure<Filtration_value> shall compile)
// Then boost::mpl::for_each<Data_structure_filtration_supported<Data_structure>::List>(my_lambda);
// will call my_lambda for all Data_structure<Filtration_value>, where Filtration_value is one of
// Filtration_value_types_supported
// 
// Equivalent to:
//
// template<template<typename> class Data_structure>
// struct Data_structure_filtration_supported {
//   using List = boost::mpl::list<
//     boost::mpl::identity<Data_structure<float>>,
//     boost::mpl::identity<Data_structure<double>>
//   >;
// };
template<template<typename> class Data_structure>
struct Data_structure_filtration_supported {

  // metafunction: identity<T> -> identity<Data_structure<T>>
  template<typename Wrapped>
  struct wrap {
    using type = boost::mpl::identity<Data_structure<typename Wrapped::type>>;
  };

  using List = typename boost::mpl::transform<
      typename Filtration_value_types_supported::List,
      wrap<boost::mpl::_1>
    >::type;
};
