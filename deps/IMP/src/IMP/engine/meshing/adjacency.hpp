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
   \file adjacency.hpp
   \brief Adjacency class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_MESHING_ADJACENCY_HPP_
#define IMP_ENGINE_MESHING_ADJACENCY_HPP_

#include "IMP/engine/elements/edge.hpp"
#include "IMP/engine/elements/triangle.hpp"
#include "IMP/engine/elements/tetrahedron.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <array>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <stdexcept>
#include <exception>

namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class Adjacency
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a Adjacency for mesh generation and processing.
 */
class Adjacency {

private:

    std::vector<std::vector<int>> faces_to_cells_;          /**< All the adjacent cells to each of the faces */

    std::vector<std::vector<int>> edges_to_faces_;          /**< All the adjacent faces to each of the edges */

    std::vector<std::vector<int>> verts_to_edges_;          /**< All the adjacent edges to each of the vertices */

public:

    /**
     * \brief The Adjacency constructor.
     */
    Adjacency();


    /**
     * \brief The Adjacency destructor.
     */
    virtual ~Adjacency();


    /**
     * @brief 
     * @param unique_faces 
     * @param cells 
     */
    void MapFacesToCells(const std::vector<Triangle> &unique_faces, const std::vector<Tetrahedron> &cells);


    /**
     * @brief 
     * @param unique_edge 
     * @param faces
     */
    void MapEdgesToFaces(const std::vector<Edge> &unique_edges, const std::vector<Triangle> &faces);


    /**
     * @brief 
     * @param unique_edge 
     * @param faces
     */
    void MapVertsToEdges(const std::vector<Vertex> &unique_verts, const std::vector<Edge> &edges);


    /**
     * @brief 
     * 
     * @return const std::vector<std::vector<int>>& 
     */
    inline const std::vector<std::vector<int>> & FacesToCells() const { return this->faces_to_cells_; }


    /**
     * @brief 
     * 
     * @return const std::vector<std::vector<int>>& 
     */
    inline const std::vector<std::vector<int>> & EdgesToFaces() const { return this->edges_to_faces_; }


    /**
     * @brief 
     * 
     * @return const std::vector<std::vector<int>>& 
     */
    inline const std::vector<std::vector<int>> & VertsToEdges() const { return this->verts_to_edges_; }


    inline bool IsReady() const { 
        return (this->faces_to_cells_.size() > 0 &&
                this->edges_to_faces_.size() > 0 &&
                this->verts_to_edges_.size() > 0);
    }
    
    
};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_MESHING_ADJACENCY_HPP_ 
