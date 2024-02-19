/*
 * IMP. Image to Mesh Processing library.
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
   \file voronoi.hpp
   \brief Voronoi class header file.
   \author Konstantinos A. Mountris
   \date 05/10/2020
*/

#pragma once
#ifndef IMP_ENGINE_TESSELATIONS_VORONOI_HPP_
#define IMP_ENGINE_TESSELATIONS_VORONOI_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/elements/poly_cell.hpp"
#include "IMP/engine/elements/poly_facet.hpp"
#include "IMP/engine/utilities/algorithms.hpp"
#include "IMP/engine/voro_io/evt_io.hpp"
#include "IMP/engine/topology/node_set.hpp"

#include <set>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <algorithm>
#include <limits>


namespace IMP {

/** \addtogroup Tesselations \{ */


/**
 * \class Voronoi
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a Voronoi tesselation.
 * \tparam DIM the spatial dimensions of the Voronoi tesselation.
 */
template<short DIM>
class Voronoi {

private:

    std::vector<Vec<DIM, double>> nodes_;                   /**< The seed nodes of the voronoi tesselation */

    std::vector<Vec<DIM, double>> points_;                  /**< The corner points of the voronoi tesselation cells */

    std::vector<PolyCell> cells_;                           /**< The cells of the voronoi tesselation */

    std::vector<PolyFacet> facets_;                         /**< The facets of the voronoi tesselation's cells */

    std::vector<std::vector<int>> facets_on_nodes_;

    std::vector<double> cell_measures_;                     /**< The measures of the cells i.e length for 1D, area for 2D, volume for 3D, etc. */

    std::unordered_map<std::string, NodeSet> node_sets_;    /**< The node sets of the voronoi tesselation */

public:

    /**
     * \brief The Voronoi constructor.
     */
    inline Voronoi();


    /**
     * \brief The Voronoi destructor.
     */
    inline virtual ~Voronoi();


    /**
     * \brief Load a voronoi tesselation from a file.
     * The following voronoi file formats are supported: *.evt
     * \param [in] filename The filename of the voronoi file to be loaded.
     * \return [void]
     */
    inline void LoadFrom(const std::string &filename);


    /**
     * \brief Save voronoi tesselation to a file.
     * The following voronoi file formats are supported: *.evt
     * \param [in] filename The filename of the voronoi file to be saved.
     */
    inline void SaveTo(const std::string &filename);



    /**
     * @brief 
     * @param scale_value 
     */
    inline void Scale(double scale_value);


    /**
     * \brief Add a new nodeset in the voronoi tesselation.
     * \param [in] nodeset The new nodeset to be added.
     * \return [void]
     */
    inline void AddNodeSet(const NodeSet &nodeset);


    /**
     * \brief Compute the measure of each cell of the voronoi tesselation.
     * Measure equals to: length - 1D, area - 2D, volume - 3D.
     * \return [void]
     */
    inline void ComputeCellMeasures();


    /**
     * @brief
     * @param id
     * @return Vec<DIM,double>
     */
    inline Vec<DIM,double> CellCentroid(std::size_t id) const;


    /**
     * \brief Get the seed nodes of the voronoi tesselation.
     * \return [const std::vector<Vec<DIM, double>>&] The seed nodes of the voronoi tesselation.
     * \overload
     */
    inline auto & Nodes() const { return this->nodes_; }


    /**
     * \brief Get the seed node of the voronoi tesselation with index id.
     * Fast access, no range check.
     * \param [in] id The index of voronoi tesselation's seed node.
     * \return [const Vec<DIM, double>&] The seed node of the voronoi tesselation with index id.
     * \overload
     */
    inline auto & Nodes(std::size_t id) const { return this->nodes_[id]; }


    /**
     * \brief Get the seed node of the voronoi tesselation with index id.
     * Slow access, with range check.
     * \param [in] id The index of voronoi tesselation's seed node.
     * \return [const Vec<DIM, double>&] The seed node of the voronoi tesselation with index id.
     */
    inline auto & NodesAt(std::size_t id) const { return this->nodes_.at(id); }


    /**
     * \brief Get the number of nodes of the voronoi tesselation.
     * \return [int] The number of nodes of the voronoi tesselation.
     */
    inline auto NodesNum() const { return static_cast<int>(this->nodes_.size()); }


    /**
     * \brief Get the corner points of the voronoi tesselation cells.
     * \return [const std::vector<Vec<DIM, double>>&] The corner points of the voronoi tesselation cells.
     * \overload
     */
    inline auto & Points() const { return this->points_; }


    /**
     * \brief Get the corner point of the voronoi tesselation with index id.
     * Fast access, no range check.
     * \param [in] id The index of voronoi tesselation's corner point.
     * \return [const Vec<DIM, double>&] The corner point of the voronoi tesselation with index id.
     * \overload
     */
    inline auto & Points(std::size_t id) const { return this->points_[id]; }


    /**
     * \brief Get the corner point of the voronoi tesselation with index id.
     * Slow access, with range check.
     * \param [in] id The index of voronoi tesselation's corner point.
     * \return [const Vec<DIM, double>&] The corner point of the voronoi tesselation with index id.
     */
    inline auto & PointsAt(std::size_t id) const { return this->points_.at(id); }


    /**
     * \brief Get the number of points of the voronoi tesselation.
     * \return [int] The number of points of the voronoi tesselation.
     */
    inline auto PointsNum() const { return static_cast<int>(this->points_.size()); }


    /**
     * \brief Get the facets of the voronoi tesselation's cells.
     * \return [const std::vector<PolyFacet>&] The cells of the voronoi tesselation.
     */
    inline auto & Facets() const { return this->facets_; }


    /**
     * \brief Get the facet of the voronoi tesselation with index id.
     * Fast access, no range check.
     * \param [in] id The index of voronoi tesselation's facet.
     * \return [const PolyFacet&] The facet of the voronoi tesselation with index id.
     * \overload
     */
    inline auto & Facets(std::size_t id) const { return this->facets_[id]; }


    /**
     * \brief Get the facet of the voronoi tesselation with index id.
     * Slow access, with range check.
     * \param [in] id The index of voronoi tesselation's facet.
     * \return [const PolyFacet&] The facet of the voronoi tesselation with index id.
     */
    inline auto & FacetsAt(std::size_t id) const { return this->facets_.at(id); }


    /**
     * \brief Get the number of facets of the voronoi tesselation.
     * \return [int] The number of facets of the voronoi tesselation.
     */
    inline auto FacetsNum() const { return static_cast<int>(this->facets_.size()); }


    /**
     * \brief Get the cells of the voronoi tesselation.
     * \return [const std::vector<PolyCell>&] The cells of the voronoi tesselation.
     * \overload
     */
    inline auto & Cells() const { return this->cells_; }


    /**
     * \brief Get the cell of the voronoi tesselation with index id.
     * Fast access, no range check.
     * \param [in] id The index of voronoi tesselation's cell.
     * \return [const PolyCell&] The cell of the voronoi tesselation with index id.
     * \overload
     */
    inline auto & Cells(std::size_t id) const { return this->cells_[id]; }


    /**
     * \brief Get the cell of the voronoi tesselation with index id.
     * Slow access, with range check.
     * \param [in] id The index of voronoi tesselation's cell.
     * \return [const PolyCell&] The cell of the voronoi tesselation with index id.
     */
    inline auto & CellsAt(std::size_t id) const { return this->cells_.at(id); }


    /**
     * \brief Get the number of cells of the voronoi tesselation.
     * \return [int] The number of cells of the voronoi tesselation.
     */
    inline auto CellsNum() const { return static_cast<int>(this->cells_.size()); }


    inline auto & FacetsOnNodes() const { return this->facets_on_nodes_; }

    inline auto & FacetsOnNodes(std::size_t id) const { return this->facets_on_nodes_[id]; }

    inline auto & FacetsOnNodesAt(std::size_t id) const { return this->facets_on_nodes_.at(id); }


    inline auto & CellMeasures() const { return this->cell_measures_; }

    inline auto & CellMeasures(std::size_t id) const { return this->cell_measures_[id]; }

    inline auto & CellMeasuresAt(std::size_t id) const { return this->cell_measures_.at(id); }



    /**
     * \brief Get the node sets of the voronoi tesselation.
     * The node sets are groups of nodes sharing common attributes.
     * \return [const std::unordered_map<std::string, IMP::NodeSet>&] The node sets of the voronoi tesselation.
     */
    inline auto & NodeSets() const { return this->node_sets_; }


    /**
     * \brief Get the node set of the voronoi tesselation that is specified by name.
     * \param [in] name The name of the node set.
     * \return [const NodeSet&] The node set of the voronoi tesselation that is specified by name.
     */
    inline auto & NodeSets(const std::string &name) const { return this->node_sets_.at(name); }


    /**
     * \brief Get the number of node sets of the voronoi tesselation.
     * \return [int] The number of node sets of the voronoi tesselation.
     */
    inline auto NodeSetsNum() const { return static_cast<int>(this->nodes_sets_.size()); }


};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.

#endif // IMP_ENGINE_TESSELATIONS_VORONOI_HPP_

#include "IMP/engine/tesselations/voronoi.tpp"