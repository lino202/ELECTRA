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
   \file poly_facet.hpp
   \brief PolyFacet class header file.
   \author Konstantinos A. Mountris
   \date 07/10/2020
*/

#ifndef IMP_ELEMENTS_POLY_FACET_HPP_
#define IMP_ELEMENTS_POLY_FACET_HPP_


#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <ostream>
#include <stdexcept>
#include <exception>
#include <vector>

namespace IMP {

/** \addtogroup Elements \{ */

/**
 * \class PolyFacet
 * \brief Class implemmenting a special half-facet for PolyCell cells.
 */
class PolyFacet {

private:

    std::vector<int> connectivity_;       /**< The connectivity of the polyfacet points */

    int parent_cell_id_;                  /**< The index of the parent cell of the polyfacet */

    int neigh_cell_id_;                   /**< The index of the neighbor cell to the polyfacet */


public:

    /**
     * \brief The PolyFacet constructor with data initialization.
     * \param [in] point_ids The indices of the corner points of the polyfacet.
     * \param [in] parent_cell_id The index of the parent cell of the polyfacet.
     * \param [in] neigh_cell_id The index of the neighbor cell to the polyfacet.
     */
    explicit PolyFacet(const std::vector<int> &point_ids, int parent_cell_id, int neigh_cell_id);


    /**
     * \brief The PolyFacet default constructor.
     */
    PolyFacet();


    /**
     * \brief The PolyFacet destructor.
     */
    virtual ~PolyFacet();


    /**
     * \brief Set the polyfacet's data.
     * \param [in] point_ids The indices of the corner points of the polyfacet.
     * \param [in] parent_cell_id The index of the parent cell of the polyfacet.
     * \param [in] neigh_cell_id The index of the neighbor cell to the polyfacet.
     * \return [void]
     * \overload
     */
    void Set(const std::vector<int> &point_ids, int parent_cell_id, int neigh_cell_id);


    /**
     * \brief Set the polyfacet's data from another polyfacet.
     * \param [in] other The other polyfacet to copy its data.
     * \return [void]
     * \overload
     */
    void Set(const PolyFacet &other);


    /**
     * \brief Set the polyfacet's corner points indices.
     * \param [in] point_ids The polyfacet's corner points indices.
     * \return [void]
     * \overload
     */
    inline void SetConnectivity(const std::vector<int> &point_ids) { this->connectivity_ = point_ids; }


    /**
     * \brief Add a new index in the connectivity of the facet.
     * \param [in] id The index of the point to be added in the connectivity.
     * \return [void]
     */
    inline void AddInConnectivity(int id) { this->connectivity_.emplace_back(id); }


    /**
     * \brief Set the polyfacet's parent cell index.
     * \param [in] parent_cell_id The polyfacet's parent cell index.
     * \return [void]
     * \overload
     */
    inline void SetParentCellId(int parent_cell_id) { this->parent_cell_id_ = parent_cell_id; }


    /**
     * \brief Set the polyfacet's neighbor cell index.
     * \param [in] neigh_cell_id The polyfacet's neighbor cell index.
     * \return [void]
     * \overload
     */
    inline void SetNeighCellId(int neigh_cell_id) { this->neigh_cell_id_ = neigh_cell_id; }


    /**
     * \brief Check if the polyfacet belong to the boundary.
     * \return [True] The index of the neighbor cell to the polyfacet is equal to -1.
     * \return [False] The index of the neighbor cell to the polyfacet is greater or equal to 0.
     */
    bool IsFree() const;


    /* \brief Compute the length of the polyfacet.
     * The polyfacet is expected to be an edge.
     * \return [double] The length of the polyfacet.
     * \overload
     */
    double Length(const std::vector<Vec<1, double>> &points_coords) const;


    /**
     * \brief Compute the length of the polyfacet.
     * The polyfacet is expected to be an edge.
     * \return [double] The length of the polyfacet.
     * \overload
     */
    double Length(const std::vector<Vec<2, double>> &points_coords) const;


    /* \brief Compute the length of the polyfacet.
     * The polyfacet is expected to be an edge.
     * \return [double] The length of the polyfacet.
     * \overload
     */
    double Length(const std::vector<Vec<3, double>> &points_coords) const;


    /**
     * \brief Compute the centroid of a 1D (point) polyfacet.
     * \param point_coords 
     * \return Vec<1, double> 
     */
    Vec<1, double> Centroid(const std::vector<Vec<1, double>> &point_coords) const;


    /**
     * \brief Compute the centroid of a 2D (edge) polyfacet.
     * \param point_coords 
     * \return Vec<2, double> 
     */
    Vec<2, double> Centroid(const std::vector<Vec<2, double>> &point_coords) const;


    /**
     * \brief Compute the centroid of a 3D (polygon) polyfacet.
     * \param point_coords 
     * \return Vec<3, double> 
     */
    Vec<3, double> Centroid(const std::vector<Vec<3, double>> &point_coords) const;


    /**
     * \brief Compute the normal to the polyfacet.
     * The polyfacet is expected to be an edge.
     * \return [Vec<2, double>] The normal to the polyfacet.
     * \overload
     */
    Vec<1, double> Normal(const std::vector<Vec<1, double>> &points_coords) const;


    /**
     * \brief Compute the normal to the polyfacet.
     * The polyfacet is expected to be an edge.
     * \return [Vec<2, double>] The normal to the polyfacet.
     * \overload
     */
    Vec<2, double> Normal(const std::vector<Vec<2, double>> &points_coords) const;


    /**
     * \brief Compute the normal to the polyfacet.
     * The polyfacet is expected to be an edge.
     * \return [Vec<2, double>] The normal to the polyfacet.
     * \overload
     */
    Vec<3, double> Normal(const std::vector<Vec<3, double>> &points_coords) const;


    /**
     * \brief Get the indices of the corner points of the polyfacet.
     * \return [const std::vector<int>&] The indices of the corner points of the polyfacet.
     */
    inline const std::vector<int> & Connectivity() const { return this->connectivity_; }


    /**
     * \brief Get the index of the corner point of the polyfacet with local index id.
     * \param [in] id The local index of the corner point of the polyfacet.
     * \return [int] The index of the corner point of the polyfacet with local index id.
     */
    inline int C(std::size_t id) const { return this->connectivity_[id]; }


    /**
     * \brief Get the index of the parent cell of the polyfacet.
     * \return [int] The index of the parent cell of the polyfacet.
     */
    inline int ParentCellId() const { return this->parent_cell_id_; }


    /**
     * \brief Get the index of neighbor cell to the polyfacet.
     * \return [int] The index of neighbor cell to the polyfacet.
     */
    inline int NeighCellId() const { return this->neigh_cell_id_; }


    /**
    * \brief Less than operator.
    * \param [in] pf The polyfacet to be compared.
    * \return [TRUE] At least one point of this polyfacet has smaller index than the corresponding point in the pf polyfacet.
    * \return [FALSE] None of the points of this polyfacet has smaller index than the corresponding point in the pf polyfacet.
    */
    bool operator < (const PolyFacet &pf) const;

};


/**
 * \brief Equality operator for Polyfacet.
 * \param [in] pf1 The first polyfacet to check for equality.
 * \param [in] pf2 The second polyfacet to check for equality.
 * \return [TRUE] Polyfacets are equal.
 * \return [FALSE] Polyfacets are not equal.
 */
bool operator == (const PolyFacet &pf1, const PolyFacet &pf2);


/**
 * \brief Output stream operator for polyfacet.
 * \param [out] out The output stream to print the polyfacet data.
 * \param [in] pf The polyfacet to be printed.
 * \return [std::ostream] The output of the polyfacet's data.
 */
std::ostream & operator << (std::ostream &out, const PolyFacet &pf);


/** \} End of Doxygen Groups*/

} //end of namespace IMP

#endif // IMP_ENGINE_ELEMENTS_POLY_FACET_HPP_
