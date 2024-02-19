/*
 * IMP. Image and Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


/**
   \file half_facets_array.hpp
   \brief HalfFacetsArray class header file.
   \author Konstantinos A. Mountris
   \date 27/09/2018
*/

#ifndef IMP_ENGINE_ELEMENTS_HALF_FACETS_ARRAY_HPP_
#define IMP_ENGINE_ELEMENTS_HALF_FACETS_ARRAY_HPP_


#include "IMP/engine/elements/half_facet.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <vector>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <utility>
#include <algorithm>

namespace IMP {

/** \addtogroup Elements \{ */

/**
 * \class HalfFacetsArray
 * \brief Class implemmenting the half-facets array of the Array-based Half-Facet (AHF) mesh representation.
 *
 * Mapping between the vertex and its attached half facets is encapsulated in the HalfFacet class.
 *
 * \see HalfFacet
 *
 */
class HalfFacetsArray {

private:

    std::vector<HalfFacet> half_facets_;    /**< The half facets of the HalfFacetsArray */

    bool is_sorted_;                        /**< Flag for the sorting state of the HalfFacetsArray */


public:

    /**
     * \brief HalfFacetsArray constructor.
     */
    HalfFacetsArray();


    /**
     * \brief HalfFacetsArray destructor.
     */
    virtual ~HalfFacetsArray();


    /**
     * \brief Empty the HalfFacetsArray container.
     * \return [void]
     */
    inline void Clear() { this->half_facets_.clear(); }


    inline void Reserve(std::size_t n) { this->half_facets_.reserve(n); }


    /**
     * \brief Get the size of the HalfFacetsArray container.
     * \return [std::size_t] The HalfFacetsArray container's size.
     */
    inline std::size_t Size() const { return this->half_facets_.size(); }


    /**
     * \overload
     * \brief Append a new half-facet in the HalfFacetsArray.
     * \param [in] vertex_id The index of the vertex to which the appended half-facet is attached.
     * \param [in] cell_id The index of the cell to which the appended half-facet belongs.
     * \param [in] facet_id The index of the appended half-facet.
     * \return [void]
     */
    void Append(int vertex_id, int cell_id, short facet_id);


    /**
     * \overload
     * \brief Append a new half-facet in the HalfFacetsArray.
     * \param [in] half_facet The half facet to be appended in the HalfFacetsArray.
     * \return [void]
     */
    void Append(const HalfFacet &half_facet);


    /**
     * \brief Sort the half facets of the HalfFacetsArray.
     *
     * The half facets are sorted in ascending order.
     * Ascending order is specified by: ascending vertex index -> ascending cell index -> ascending half facet index
     *
     */
    void Sort();



    /**
     * \brief Get all the half-facets attached to a vertex.
     *
     * The attached half-facets, to a vertex with the given index, are specified by
     * the range [first, last) of the half-facet positions in the HalfFacetsArray.
     *
     * \param [in] vertex_id The index of the vertex to search for its attached half-facets.
     * \return [std::pair<int, int>] The range [first, last) of the first and last position
     *         of the attached half-facets in the HalfFacetsArray.
     */
    std::pair<std::vector<HalfFacet>::const_iterator,
    std::vector<HalfFacet>::const_iterator> AllAttachedToVertex(int vertex_id) const;


    /**
     * \brief Get the position in the HalfFacetsArray of the first half-facet
     *        that is attached to the vertex with the given index.
     * \param [in] vertex_id The index of the vertex to search for its first attached half-facet.
     * \return [int] The position in the HalfFacetsArray of the first attached half-facet to the given vertex.
     */
    std::vector<HalfFacet>::const_iterator FirstAttachedToVertex(int vertex_id) const;


    /**
     * \brief Get the indices of the mapped vertices by the half-facets.
     * \return [std::vector<int>] The indices of the mapped vertices by the half-facets.
     */
    std::vector<int> MappedVerticesIds() const;


    /**
     * \brief Get read-access to the half facets of the HalfFacetsArray.
     * \return [const std::vector<HalfFacet>&] The half facets of the HalfFacetsArray.
     */
    inline const std::vector<HalfFacet> &HalfFacets() const { return this->half_facets_; }


    /**
     * \brief Get write-access to the half facets of the HalfFacetsArray.
     * \return [std::vector<HalfFacet>&] The half facets of the HalfFacetsArray.
     */
    inline std::vector<HalfFacet> &EditHalfFacets() { return this->half_facets_; }


    /**
     * \brief Check if the half facets in the HalfFacetsArray are sorted
     *        in the specifieed order of the method HalfFacetsArray::Sort.
     * \see HalfFacetsArray::Sort
     * \return [bool] True if the HalfFacetsArray is sorted | False otherwise.
     */
    inline bool IsSorted() const { return this->is_sorted_; }


    /**
     * \brief Check if the HalfFacetsArray is empty.
     * \return [bool] True if the HalfFacetsArray is empty | False otherwise.
     */
    inline bool IsEmpty() const { return this->half_facets_.empty(); }

};



/** \} End of Doxygen Groups*/
} //end of namespace IMP

#endif // IMP_ENGINE_ELEMENTS_HALF_FACETS_ARRAY_HPP_
