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
   \file half_facet.hpp
   \brief HalfFacet class header file.
   \author Konstantinos A. Mountris
   \date 27/09/2018
*/

#ifndef IMP_ELEMENTS_HALF_FACET_HPP_
#define IMP_ELEMENTS_HALF_FACET_HPP_


#include "IMP/engine/utilities/logger.hpp"

#include <ostream>
#include <stdexcept>
#include <exception>


namespace IMP {

/** \addtogroup Elements \{ */

/**
 * \class HalfFacet
 * \brief Class implemmenting a half-facet of a cell in the Array-based Half-Facet (AHF) mesh representation.
 * \note Should I continue assigning cell id to the half facet and in the same time half facets in the cells?
 *
 */
class HalfFacet {

private:

    int vertex_id_;       /**< The vertex's index on which the half-facet is attached */

    int cell_id_;         /**< The cell's index on which the half-facet belongs */

    short facet_id_;      /**< The half-facet's local index in the parent cell */


public:

    /**
     * \brief HalfFacet constructor with data initialization.
     * \param [in] vertex_id The vertex index on which the half-facet is attached.
     * \param [in] cell_id The cell index to which the half-facet belong.
     * \param [in] facet_id The local index of the half-facet in the parent cell.
     */
    explicit HalfFacet(int vertex_id, int cell_id, short facet_id);


    /**
     * \brief HalfFacet default constructor.
     */
    HalfFacet();


    /**
     * \brief HalfFacet destructor.
     */
    virtual ~HalfFacet();


    /**
     * \brief Set the half-facet's data.
     * \param [in] vertex_id The vertex index on which the half-facet is attached.
     * \param [in] cell_id The cell index to which the half-facet maps.
     * \param [in] facet_id The local index of the half-facet in the mapped cell.
     * \return [void]
     * \overload
     */
    void Set(int vertex_id, int cell_id, short facet_id);


    /**
     * \brief Set the half-facet's data from another vertex-half-facet.
     * \param [in] other The other vertex-half-facet to copy its data.
     * \return [void]
     * \overload
     */
    void Set(const HalfFacet &other);


    /**
     * \brief Check if the vertex-half-facet is null.
     * \return [bool] True if the half-facet is null | False otherwise.
     */
    bool IsNull() const;


    /**
     * \brief Equality comparison operator for half-facets.
     * \param [in] hf1 The fist half-facet to check for equality.
     * \param [in] hf2 The second half-facet to check for equality.
     * \return [bool] True if the data of the two half-facets are equal | False otherwise.
     */
    friend bool operator == (const HalfFacet &hf1, const HalfFacet &hf2);


    /**
     * \brief Output stream operator for half-facet.
     * \param [out] out The output stream to print the vertex-half-facet data.
     * \param [in] hf The half-facet to be printed.
     * \return [std::ostream] The output of the half-facet's data.
     */
    friend std::ostream & operator << (std::ostream &out, const HalfFacet &vhf);


    /**
     * \brief Get the index of the vertex on which the half-facet is attached.
     * \return [int int] The index of the vertex on which the half-facet is attached.
     */
    inline int VertexId() const { return this->vertex_id_; }


    /**
     * \brief Get the index of the mapped cell of the half-facet.
     * \return [int int] The index of the mapped cell of the half-facet.
     */
    inline int CellId() const { return this->cell_id_; }


    /**
     * \brief Get the local index of the half-facet in the mapped cell.
     * \return [short int] The local index of the half-facet in the mapped cell.
     */
    inline short FacetId() const { return this->facet_id_; }


};

/** \} End of Doxygen Groups*/
} //end of namespace IMP

#endif // IMP_ENGINE_ELEMENTS_HALF_FACET_HPP_
