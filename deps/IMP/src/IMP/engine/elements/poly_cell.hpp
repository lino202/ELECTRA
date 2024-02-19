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
   \file polycell.hpp
   \brief PolyCell class header file.
   \author Konstantinos A. Mountris
   \date 01/10/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_POLY_CELL_HPP_
#define IMP_ENGINE_ELEMENTS_POLY_CELL_HPP_


#include "IMP/engine/elements/cell_props.hpp"
#include "IMP/engine/elements/poly_facet.hpp"
#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <exception>
#include <cmath>
#include <type_traits>


namespace IMP {

/** \addtogroup Elements \{ */


/**
 * \class PolyCell
 * \brief Template class implementing a generic polygonal/polyhedral cell.
 * \tparam DIM The number of the polygonal/polyhedral cell's dimensions.
 */
class PolyCell {

private:

    std::vector<int> connectivity_;    /**< The global indices of the polycell's points or facets */

    const CellShape shape_;            /**< The type of the polycell's shape */


public:

    /**
     * \brief The PolyCell default constructor.
     * \overload
     */
    PolyCell();


    /**
     * \brief The PolyCell copy constructor.
     * \overload
     */
    PolyCell(const PolyCell &pc);


    /**
     * \brief The PolyCell destructor.
     */
    virtual ~PolyCell();


    /**
     * \brief Set the polycell's connectivity.
     * For 1D and 2D cells refers to the connectivity of points.
     * For 3D cells refer to the connectivity of facets.
     * \param [in] ids The points' indices of the polycell.
     * \return [void]
     */
    void SetConnectivity(const std::vector<int> &ids);


    /**
     * \brief Add a new index in the connectivity of the cell.
     * \param [in] id The index of the point for 1D and 2D cells or the facet for 3D cells to be added in the connectivity.
     * \return [void]
     */
    void AddInConnectivity(int id);


    /**
     * \brief Compute the measure of the polycell.
     * Measure is a generalized concept corresponding to the polycell's length in 1D.
     * \return [double] The measure of the cell.
     * \overload
     */
    double Measure(const std::vector<Vec<1,double>> &points_coords, const std::vector<PolyFacet> &facets) const;


    /**
     * \brief Compute the measure of the polycell.
     * Measure is a generalized concept corresponding to the polycell's area in 2D.
     * \return [double] The measure of the cell.
     * \overload
     */
    double Measure(const std::vector<Vec<2,double>> &points_coords, const std::vector<PolyFacet> &facets) const;


    /**
     * \brief Compute the measure of the polycell.
     * Measure is a generalized concept corresponding to the polycell's volume in 3D.
     * \return [double] The measure of the cell.
     * \overload
     */
    double Measure(const std::vector<Vec<3,double>> &points_coords, const std::vector<PolyFacet> &facets) const;


    /**
     * \brief Compute the centroid of a 1D (line) polycell.
     * \param point_coords
     * \return Vec<1, double>
     */
    Vec<1, double> Centroid(const std::vector<Vec<1, double>> &point_coords, const std::vector<PolyFacet> &facets) const;


    /**
     * \brief Compute the centroid of a 2D (polygon) polycell.
     * \param point_coords
     * \return Vec<2, double>
     */
    Vec<2, double> Centroid(const std::vector<Vec<2, double>> &point_coords, const std::vector<PolyFacet> &facets) const;


    /**
     * \brief Compute the centroid of a 3D (polyhedral) polycell.
     * \param point_coords
     * \return Vec<3, double>
     */
    Vec<3, double> Centroid(const std::vector<Vec<3, double>> &point_coords, const std::vector<PolyFacet> &facets) const;


    /**
     * \brief Get the type of the polycell's shape.
     * \return [IMP::CellShape] The type of the polycell's shape.
     */
    inline CellShape Shape() const { return this->shape_; }


    /**
     * \brief Get the indices of the polycell points.
     * \return [const IMP::Vec<short, int>&] The indices of the polycell points.
     */
    inline const std::vector<int> & Connectivity() const { return this->connectivity_; }


    /**
     * \brief Get the index of the point (1D/2D) or facet (3D) with local index id in the polycell's connectivity.
     * \param [in] id The local index of the point (1D/2D) or facet (3D) in the polycell's connectivity.
     * \return [int] The global index of the point (1D/2D) or facet (3D) in the polycell's connectivity.
     */
    inline int C(std::size_t id) const { return this->connectivity_[id]; }


    /**
     * \brief Copy-assignment operator.
     * \param [in] pc The polycell to be copy-assigned.
     * \return [IMP::PolyCell] This polycell containing the copied data from pc polycell.
     */
    PolyCell & operator = (const PolyCell &pc);


    /**
     * \brief Move-assignment operator.
     * \param [in] pc The polycell to be move-assigned.
     * \return [IMP::PolyCell] This polycell containing the moved data from pc polycell.
     */
    PolyCell & operator = (PolyCell &&pc) noexcept;

};


/**
 * \brief Equality operator for PolyCell.
 * \param [in] pc1 The first polycell to check for equality.
 * \param [in] pc2 The second polycell to check for equality.
 * \return [TRUE] PolyCells are equal.
 * \return [FALSE] PolyCells are not equal.
 */
bool operator == (const PolyCell &pc1, const PolyCell &pc2);


/**
 * \brief No equality operator for PolyCell.
 * \param [in] pc1 The first polycell to check for no equality.
 * \param [in] pc2 The second polycell to check for no equality.
 * \return [TRUE] PolyCells are not equal.
 * \return [FALSE] PolyCells are equal.
 */
bool operator != (const PolyCell &pc1, const PolyCell &pc2);


/**
 * \brief Output stream operator for PolyCell.
 * \param [out] out The outstream operator.
 * \param [in] pc The polycell to print in outstream.
 * \return [std::ostream] The outstream with the polycell's info.
 */
std::ostream & operator << (std::ostream &out, const PolyCell &pc);



/** \} End of Doxygen Groups*/
} // End of namespace IMP

#endif //IMP_ENGINE_ELEMENTS_POLY_CELL_HPP_