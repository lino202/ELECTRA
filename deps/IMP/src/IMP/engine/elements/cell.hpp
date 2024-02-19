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
   \file cell.hpp
   \brief Cell class header file.
   \author Konstantinos A. Mountris
   \date 28/09/2018
*/

#ifndef IMP_ENGINE_ELEMENTS_CELL_HPP_
#define IMP_ENGINE_ELEMENTS_CELL_HPP_


#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/elements/half_facet.hpp"
#include "IMP/engine/elements/cell_props.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <vector>
#include <algorithm>
#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <exception>


namespace IMP {

/** \addtogroup Elements \{ */


/**
 * \class Cell
 * \brief Template class implementing a generic cell of the Array-based Half-Facet (AHF) mesh representation.
 * \tparam DIM The number of the cell's dimensions.
 * \tparam CELL_VERTS The number of the cell's vertices.
 */
template <short DIM, short CELL_VERTS>
class Cell {

private:

    std::vector<HalfFacet> sib_half_facets_;    /**< The sibling half-facets to the cell's facets */

    Vec<CELL_VERTS, int> connectivity_;         /**< The global indices of the cell's vertices */

    CellShape shape_;                           /**< The type of the cell's shape */


public:

    /**
     * \brief The Cell constructor.
     */
    inline Cell();


    /**
     * \brief The Cell destructor.
     */
    inline virtual ~Cell();


    /**
     * \brief Set the type of the cell's shape.
     * \param [in] cell_shape The type of the cell's shape.
     * \return [void]
     */
    inline void SetShape(CellShape cell_shape);


    /**
     * \overload
     * \brief Set the cell's connectivity.
     * \param [in] connectivity The vertices' connectivity of the cell.
     * \return [void]
     */
    inline void SetConnectivity(const Vec<CELL_VERTS, int> &connectivity);


    /**
     * \overload
     * \brief Set the cell's connectivity.
     * \param [in] connectivity The vertices' connectivity of the cell.
     * \return [void]
     */
    inline void SetConnectivity(std::initializer_list<int> connectivity);


    /**
     * \brief Set the attached half-facets to the vertices of the cell.
     * \param [in] half_facets The cell's half-facets.
     * \return [void]
     */
    inline void SetAllSibHfacets(const std::vector<HalfFacet> &half_facets);


    /**
     * \brief Set the attached half-facet of a single node of the cell.
     * \param [in] local_node_id The local index of the cell's node.
     * \param [in] half_facet The half-facet to be assigned.
     * \return [void]
     */
    inline void SetSibHfacet(short local_node_id, const HalfFacet &half_facet);


    /**
     * \brief Get the connectivity of a half-facet in the cell.
     * \param [in] half_facet_id The local index of the half facet to extract it's connectivity.
     * \return [Vec<CELL_VERTS-1, int>] The connectivity of the cell's half facet.
     */
    inline Vec<(CELL_VERTS-1), int> FacetConnectivity(short half_facet_id) const;


    /**
     * \brief Get the indices of the adjacent nodes to a querry node in a given half-facet.
     * \param [in] node_id The index of the querry node.
     * \param [in] facet_id The index of the half-facet where this node belongs.
     */
    inline std::vector<int> AdjacentNodeIdsTo(int node_id, short hfacet_id) const;


    /**
     * \brief Get the largest global node index in the cell's half-facet.
     * \param [in] hfacet_id The local index of the half-facet for which the maximum node index will be extracted.
     * \return [int] The largest global node index in the half-facet with index hfacet_id.
     */
    inline int MaxNodeIdInHfacet(short hfacet_id) const;


    /**
     * \brief Get the local index in the cell of a given node.
     *
     * \note If the node does not belong in the cell the returned local index is -1.
     *
     * \param [in] global_vert_id The global index of the node of interest.
     * \return [short] The local index of the node in the cell. If not found in the cell, the assigned local index is -1.
     */
    inline short LocalNodeIdOf(int global_node_id) const;


    /**
     * \brief
     *
     * \param [in] local_node_id
     * \return [std::vector<short>]
     */
    inline std::vector<short> HfacetIdsOnLocalNode(short local_node_id) const;


    /**
     * \brief
     *
     * \param half_facet_id
     * \return [std::vector<short>]
     */
    inline std::vector<short> LocalNodeIdsInHfacet(short half_facet_id) const;


    /**
     * \brief
     *
     * \param
     *
     * \return [true]
     * \return [false]
     */
    inline bool IsConnectedToCell(int neigh_cell_id) const;


    /**
     * \brief Get the type of the cell's shape.
     * \return [IMP::CellShape] The type of the cell's shape.
     */
    inline CellShape Shape() const { return this->shape_; }


    /**
     * \brief Get the cell's connectivity.
     * \return [const IMP::Vec<short, int>&] The cell's connectivity.
     */
    inline const Vec<CELL_VERTS, int> & Connectivity() const { return this->connectivity_; }


    /**
     * \brief Get the global index of the node with local index id in the connectivity.
     * \param [in] id The local index of the node in the connectivity.
     * \return [int] The global index of the node.
     */
    inline int N(std::size_t id) const { return this->connectivity_[id]; }


    /**
     * \brief Get the cell's attached half facets.
     * \return [const std::vector<IMP::HalfFacet>&] The cell's attached half facets.
     */
    inline const std::vector<HalfFacet> & SibHfacets() const { return this->sib_half_facets_; }

};


/**
 * \brief Equality operator for Cell.
 * \param [in] c1 The first cell to check for equality.
 * \param [in] c2 The second cell to check for equality.
 * \return [TRUE] Cells are equal.
 * \return [FALSE] Cells are not equal.
 */
template <short DIM, short CELL_VERTS>
bool operator == (const Cell<DIM, CELL_VERTS> &c1, const Cell<DIM, CELL_VERTS> &c2);


/**
 * \brief Output stream operator for Cell.
 * \param [out] out The outstream operator.
 * \param [in] cell The cell to print in outstream.
 * \return [std::ostream] The outstream with the cell's info.
 */
template <short DIM, short CELL_VERTS>
std::ostream & operator << (std::ostream &out, const Cell<DIM, CELL_VERTS> &cell);



/** \} End of Doxygen Groups*/
} // End of namespace IMP

#endif //IMP_ENGINE_ELEMENTS_CELL_HPP_

#include "IMP/engine/elements/cell.tpp"
