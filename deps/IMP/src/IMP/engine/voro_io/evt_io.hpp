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
   \file abaqus_io.hpp
   \brief Abaqus input/output class header file.
   \author Konstantinos A. Mountris
   \date 07/05/2017
*/

#pragma once
#ifndef IMP_ENGINE_VORO_IO_EVT_IO_HPP_
#define IMP_ENGINE_VORO_IO_EVT_IO_HPP_

#include "IMP/engine/vectors/vec.hpp"
#include "IMP/engine/topology/node_set.hpp"
#include "IMP/engine/elements/poly_cell.hpp"
#include "IMP/engine/elements/poly_facet.hpp"
#include "IMP/engine/elements/cell_props.hpp"
#include "IMP/engine/utilities/logger.hpp"
#include "IMP/engine/utilities/algorithms.hpp"

#include <vector>
#include <unordered_map>
#include <iterator>
#include <utility>
#include <algorithm>
#include <string>
#include <filesystem>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>

namespace fs = std::filesystem;
namespace IMP {

/** \addtogroup VoroIO \{ */


/**
 * \brief Class implementing input/output for voronoi tesselation in ELECTRA Voronoi Tesselation (*.evt) format.
 * \tparam DIM The dimensions of the voronoi tesselation.
 */
template <short DIM>
class EvtIO {

private:

    std::vector<std::string> parsed_voronoi_;                           /**< The container storing the parsed voronoi tesselation for random-access to any line */

    std::unordered_map<std::size_t, std::size_t> mapped_node_ids_;      /**< The indices of the offsetted nodes from the storage order and their offsets */

    std::unordered_map<std::size_t, std::size_t> mapped_point_ids_;     /**< The indices of the offsetted points from the storage order and their offsets */

    std::unordered_map<std::size_t, std::size_t> mapped_facet_ids_;      /**< The indices of the offsetted faces from the storage order and their offsets */

    std::vector<int> nsets_startlines_;                                 /**< The indices of the starting line of the node sets in the voronoi filie */

    int nodes_startline_;                                               /**< The index of the starting line of the nodes in the parsed voronoi file */

    int points_startline_;                                              /**< The index of the starting line of the points in the parsed voronoi file */

    int facets_startline_;                                              /**< The index of the starting line of the facets in the parsed voronoi file */

    int cells_startline_;                                               /**< The index of the starting line of the cells in the parsed voronoi file */

    bool has_nsets_;                                                    /**< Conditional to check if the voronoi has node sets */


public:

    /**
     * \brief The EvtIO constructor.
     */
    inline EvtIO();


    /**
     * \brief The EvtIO destructor.
     */
    inline virtual ~EvtIO();


    /**
     * \brief Read a voronoi tesselation.
     *
     * The voronoi tessealation to be readed should be in ELECTRA Voronoi Tesselation format (.evt).
     *
     * \param [in] filename The name (full path) of the voronoi tesselation file.
     * \return [void]
     */
    inline void ReadVoronoiFrom(const std::string &filename);


    /**
     * \brief Save voronoi tesselation in ELECTRA Voronoi Tesselation format.
     * \param [in] nodes The nodes of the voronoi tessalation to be saved.
     * \param [in] points The points of the voronoi tessalation to be saved.
     * \param [in] cells The cells of the voronoi tessalation to be saved.
     * \param [in] filename The name (full path) of the file where the voronoi tesselation will be saved.
     * \return [void]
     */
    inline void SaveVoronoiTo(const std::string &filename, const std::vector<Vec<DIM, double>> &nodes,
        const std::unordered_map<std::string, NodeSet> &node_sets, const std::vector<Vec<DIM, double>> &points,
        const std::vector<PolyFacet> &facets, const std::vector<PolyCell> &cells) const;


    /**
     * \brief Load the nodes of the parsed voronoi tesselation in the given nodes container.
     * \param [out] nodes The container to load the nodes of the voronoi tesselation.
     * \return [void]
     */
    inline void LoadNodesIn(std::vector<Vec<DIM, double>> &nodes);


    /**
     * \brief Load the points of the parsed voronoi tesselation in the given points container.
     * \param [out] points The container to load the points of the voronoi tesselation.
     * \return [void]
     */
    inline void LoadPointsIn(std::vector<Vec<DIM, double>> &points);


    /**
     * \brief Load the facets of the readed mesh in the given container.
     * \param [out] facets The container to load the facets of the mesh.
     * \return [void]
     */
    inline void LoadFacetsIn(std::vector<PolyFacet> &facets);


    /**
     * \brief Load the cells of the readed mesh in the given container.
     * \param [out] cells The container to load the cells of the mesh.
     * \return [void]
     */
    inline void LoadCellsIn(std::vector<PolyCell> &cells);


    /**
     * \brief Load the voronoi tesselation node sets in the given node sets container.
     * \param [out] node_sets The container to load the node sets of the voronoi tesselation.
     * \return [void]
     */
    inline void LoadNodeSetsIn(std::unordered_map<std::string, NodeSet> &node_sets);


    /**
     * \brief Checks if node sets are available in the voronoi tesselation.
     * \return [true] Node sets are available in the voronoi tesselation.
     * \return [false] There are no node sets in the voronoi tesselation.
     */
    inline const bool & HasNodeSets() const { return this->has_nsets_; }

};


/** \} End of Doxygen Groups*/

} //end of namespace IMP

#endif // IMP_ENGINE_VORO_IO_EVT_IO_HPP_

#include "IMP/engine/voro_io/evt_io.tpp"
