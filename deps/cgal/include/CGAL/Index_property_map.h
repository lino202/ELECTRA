// Copyright (c) 20009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-5.0.1/Property_map/include/CGAL/Index_property_map.h $
// $Id: Index_property_map.h 52164b1 2019-10-19T15:34:59+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Saboret

#ifndef CGAL_INDEX_PROPERTY_MAP_H
#define CGAL_INDEX_PROPERTY_MAP_H

#include <CGAL/property_map.h>
#include <boost/shared_ptr.hpp>

#include <iterator>
#include <map>

namespace CGAL {

/// \cond SKIP_IN_MANUAL
  
// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// Functor for operator< that compares iterators address.
template <typename Iterator>
struct Compare_iterator_address
{
  bool operator()(const Iterator& lhs, const Iterator& rhs) const
  {
    return (&*lhs < &*rhs);
  }
};


} /* namespace internal */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Template class "index" property map, which associates a 0-based index (unsigned int) 
/// to the `[first, beyond)` range of elements.
/// 
/// 2 specializations exist:
/// - if Iter is a random access iterator (typically vector and deque), 
/// get() just calls std::distance() and is very efficient;
/// - else, the property map allocates a std::map to store indices
/// and get() requires a lookup in the map.
///
/// @heading Is Model for the Concepts:
/// Model of the `ReadablePropertyMap` concept.
///
/// @heading Parameters:
/// @param Iter iterator over input elements.
/// @param iterator_tag Iter's iterator category.

// This is the default variant that creates a temporary std::map<Iter,unsigned int>.
template <class Iter,
          typename iterator_tag=typename std::iterator_traits<Iter>::iterator_category>
class Index_property_map 
  : public boost::associative_property_map< std::map<Iter,
                                                     unsigned int,
                                                     internal::Compare_iterator_address<Iter> > >
{
  // std::map to store indices
  typedef typename std::map<Iter,
                            unsigned int,
                            internal::Compare_iterator_address<Iter> >
                                            Index_map;
  // base class = property map
  typedef typename boost::associative_property_map<Index_map>
                                            Base;
public:
  Index_property_map(
    Iter first,  ///< iterator over the first element (index 0)
    Iter beyond) ///< past-the-end iterator over the elements
    
  : m_index_map(new Index_map)  // Allocate std::map
  {
    // Index elements in std::map
    Iter it;
    unsigned int index;
    for (it = first, index = 0; it != beyond; it++, index++)
      (*m_index_map)[it] = index;
      
    // Wrap std::map in property map
    (Base&)(*this) = *m_index_map;
  }

private:
  // Property maps must be lightweight classes => share std::map
  boost::shared_ptr<Index_map> m_index_map;
};

/// @cond SKIP_IN_MANUAL
// This variant is optimized for a random access container.
template <class Iter> 
class Index_property_map<Iter,
                         std::random_access_iterator_tag>
{
public:
  // Property maps' required types
  typedef boost::readable_property_map_tag  category;
  typedef std::size_t                       value_type;
  typedef value_type                        reference;
  typedef Iter                              key_type;

  Index_property_map(
    Iter first,  ///< iterator over the first element (index 0)
    Iter /*beyond*/) ///< past-the-end iterator over the elements
  : m_first(first)
  {}

  /// Free function to access the map elements.
  friend inline
  reference get(const Index_property_map& map, key_type p)
  {
    return std::distance(map.m_first, p);
  }

private:
  Iter m_first; // iterator over the first element (index 0)
};
/// @endcond


/// Free function to create an Index_property_map property map.
///
/// @tparam Iter iterator over input elements.
///
/// @return an "index" property map, which associates a 0-based index (unsigned int) 
/// to the `[first, beyond)` range of elements.
template <class Iter>
Index_property_map<Iter>
make_index_property_map(
    Iter first,  ///< iterator over the first element (index 0)
    Iter beyond) ///< past-the-end iterator over the elements
{
  return Index_property_map<Iter>(first, beyond);
}

/// \endcond

} // namespace CGAL

#endif // CGAL_INDEX_PROPERTY_MAP_H
