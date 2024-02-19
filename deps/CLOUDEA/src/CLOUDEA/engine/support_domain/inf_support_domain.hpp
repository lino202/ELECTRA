/*
 * CLOUDEA - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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

/*!
   \file support_domain.hpp
   \brief SupportDomain class header file.
   \author Konstantinos A. Mountris
   \date 21/05/2017
*/

#ifndef CLOUDEA_SUPPORT_DOMAIN_INF_SUPPORT_DOMAIN_HPP_
#define CLOUDEA_SUPPORT_DOMAIN_INF_SUPPORT_DOMAIN_HPP_


#include "CLOUDEA/engine/elements/element_properties.hpp"
#include "CLOUDEA/engine/elements/node.hpp"
#include "CLOUDEA/engine/elements/tetrahedron.hpp"
#include "CLOUDEA/engine/vectors/vec3.hpp"
#include "CLOUDEA/engine/utilities/logger.hpp"

#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <exception>
#include <limits>


namespace CLOUDEA {

/*!
 *  \addtogroup ShapeFunctions
 *  @{
 */


/*!
 * \class InfSupportDomain
 * \brief Class implemmenting the support domain of nodes for meshless shape functions construction for MMLS 3D.
 * \deprecated 
 * \bug The calculation of min max support nodes during the closest nodes calculation can generate errors if closest nodes are generated multiple times.
 * \bug closest nodes should be part of the class probably and not returned.
 */

class InfSupportDomain
{
public:

    /*!
     * \brief SupportDomain constructor.
     */
    InfSupportDomain();


    /*!
     * \brief SupportDomain destructor.
     */
    virtual ~InfSupportDomain();


    /*!
     * \brief Set the influence nodes of the support domain.
     * \param[in] nodes The influence nodes of the support domain.
     * \return [void]
     */
    void SetInfluenceNodes(const std::vector<Node> &nodes);


    /*!
     * \brief Set the influence tetrahedra of the support domain.
     * \param[in] tetras The influence tetrahedra of the support domain.
     * \return [void]
     */
    void SetInfluenceTetrahedra(const std::vector<Tetrahedron> &tetras);


    /*!
     * \brief Compute the radiuses of influence of the influence nodes in the support domain.
     * \param[in] dilation_coeff A dilatation coefficient to increase the influence radiuses of the nodes.
     * \return [void]
     */
    void ComputeInfluenceNodesRadiuses(double dilatation_coeff = std::numeric_limits<double>::min());


    /*!
     * \brief Get the influence nodes of the support domain.
     * \return [std::vector<CLOUDEA::Node>] The influence nodes of the support domain.
     */
    inline const std::vector<Node> & InfluenceNodes() const { return this->influence_nodes_; }


    /*!
     * \brief Get the influence tetrahedra of the support domain.
     * \return [std::vector<CLOUDEA::Tetrahedron>] The influence tetrahedra of the support domain.
     */
    inline const std::vector<Tetrahedron> & InfluenceTetrahedra() const { return this->influence_tetras_; }


    /*!
     * \brief Get the influence radiuses of the influence nodes of the support domain.
     * \return [std::vector<double>] The influence radiuses of the influence nodes of the support domain.
     */
    inline const std::vector<double> & InfluenceNodesRadiuses() const { return this->influence_radiuses_; }


    /*!
     * \brief Get the indices of the closest influence nodes to each of the given evaluation nodes using exhaustive search.
     * \param[in] eval_nodes_coords The coordinates of the evaluation nodes for which the closest influence nodes indices will be returned.
     * \return [std::vector<std::vector<int> >] The indices of the influence nodes closest to each of the given evaluation node.
     */
    const std::vector<std::vector<int> > ClosestNodesIdsTo(const std::vector<Vec3<double> > &eval_nodes_coords) const;


    /*!
     * \brief Get the indices of the closest influence nodes to each of the given evaluation nodes using CGAL kNN search in sphere.
     * \param[in] eval_nodes_coords The coordinates of the evaluation nodes for which the closest influence nodes indices will be returned.
     * \return [std::vector<std::vector<int> >] The indices of the influence nodes closest to each of the given evaluation node.
     */
    const std::vector<std::vector<int> > CgalClosestNodesIdsTo(const std::vector<Vec3<double> > &eval_nodes_coords) const;


    /*!
     * \brief Get the indices of the closest influence nodes to each of the given evaluation nodes loaded from a file.
     * \param[in] closest_nodes_file The file with the neighbor nodes list for each evaluation node.
     * \return [std::vector<std::vector<int> >] The indices of the influence nodes closest to each of the given evaluation node.
     */
    const std::vector<std::vector<int> > LoadClosestNodesFrom(const std::string &closest_nodes_file) const;


    /*!
     * \brief Get the number of support nodes contained in the smallest support domain.
     * \return [int] The number of support nodes contained in the smallest support domain.
     */
    int MinSupportNodesIn(const std::vector<std::vector<int> > &neighbor_ids) const;


    /*!
     * \brief Get the number of support nodes contained in the smallest support domain.
     * \return [int] The number of support nodes contained in the smallest support domain.
     */
    int MaxSupportNodesIn(const std::vector<std::vector<int> > &neighbor_ids) const;


private:
    std::vector<Node> influence_nodes_;                /*!< The influence nodes of the support domain. */

    std::vector<Tetrahedron> influence_tetras_;        /*!< The influence tetrahedra of the support domain. */

    std::vector<double> influence_radiuses_;           /*!< The influence radiuses of the influence nodes of the support domain. */

};


/*! @} End of Doxygen Groups*/
} //end of namespace CLOUDEA

#endif //CLOUDEA_SUPPORT_DOMAIN_INF_SUPPORT_DOMAIN_HPP_
