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
   \file tet_to_poly.hpp
   \brief Edge class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_MESHING_TET_TO_POLY_HPP_
#define IMP_ENGINE_MESHING_TET_TO_POLY_HPP_


#include "IMP/engine/meshing/adjacency.hpp"
#include "IMP/engine/elements/vertex.hpp"
#include "IMP/engine/elements/edge.hpp"
#include "IMP/engine/elements/polygon.hpp"
#include "IMP/engine/elements/polyhedron.hpp"
#include "IMP/engine/elements/triangle.hpp"
#include "IMP/engine/elements/tetrahedron.hpp"
#include "IMP/engine/tesselations/mesh.hpp"
#include "IMP/engine/utilities/algorithms.hpp"
#include "IMP/engine/topology/node_set.hpp"

#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <fstream>
#include <utility>
#include <iterator>


namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class TetToPoly
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting conversion of a tetrahedral mesh to polyhedral.
 */
class TetToPoly {

private:

    Adjacency primal_map_;                      /**< Map for low-dimension elements to higher-dimension elements for the primal mesh */

    Adjacency dual_map_;                        /**< Map for low-dimension elements to higher-dimension elements for the dual mesh */

    std::vector<Vertex> primal_verts_;          /**< The vertices of the primal mesh */

    std::vector<Vertex> dual_verts_;            /**< The vertices of the dual mesh */

    std::vector<Tetrahedron> primal_cells_;     /**< The cells of the primal mesh */

    std::vector<Triangle> primal_faces_;        /**< The faces of the primal mesh */

    std::vector<Polyhedron> dual_cells_;        /**< The cells of the dual mesh*/

    std::vector<Polygon> dual_faces_;           /**< The faces of the dual mesh */

    std::vector<Edge> primal_edges_;            /**< The edges of the primal mesh */

    std::vector<Edge> dual_edges_;              /**< The edges of the dual mesh */

    std::unordered_map<std::string, NodeSet> dual_node_sets;            /**< The nodesets of the dual mesh */

    std::unordered_map<int,std::vector<int>> vertex_to_dual_face_id_;

    std::unordered_map<int,int> edge_to_dual_face_id_;

    std::unordered_map<int,int> cell_to_dual_point_id_;

    std::unordered_map<int,int> face_to_dual_point_id_;

    std::unordered_map<int,int> edge_to_dual_point_id_;

    std::unordered_map<int,int> vertex_to_dual_point_id_;

    double feat_angle_;

public:

    /**
     * \brief The TetToPoly constructor.
     */
    TetToPoly();


    /**
     * \brief The TetToPoly destructor.
     */
    virtual ~TetToPoly();


    /**
     * \brief Set the angle to extract feautures in the primal mesh.
     * \param [in] angle The feature angle given in degrees. 
     * \return [void]
     */
    inline void SetFeatureAngle(double feat_angle) { this->feat_angle_ = feat_angle; }


    /**
     * @brief Clear the data of the primal mesh.
     * 
     */
    void ClearPrimalMesh();


    /**
     * @brief Load the primal mesh.
     * 
     * @param primal_mesh 
     */
    void LoadPrimalMesh(const Mesh<3,4> &primal_mesh);


    /**
     * @brief Orient the primal cells to counter-clockwise orientation.
     * 
     */
    void OrientPrimalCells();


    /**
     * @brief Extract all the unique faces in the primal cells.
     * 
     */
    void ExtractPrimalFaces();


    /**
     * @brief Extract all the unique edges in the primal faces.
     * 
     */
    void ExtractPrimalEdges();


    /**
     * @brief 
     * 
     */
    void FindEdgesOnFeature();


    /**
     * @brief 
     * 
     */
    void FindVertsOnBoundaryAndFeature();


    /**
     * @brief 
     * 
     */
    void GenerateDualVerts();


    /**
     * @brief 
     * 
     */
    void GenerateDualFaces();


    /**
     * @brief 
     * 
     */
    void GenerateDualCells();


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SavePrimalFaces(const std::string &outfile);


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SavePrimalEdges(const std::string &outfile);


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SavePrimalVerts(const std::string &outfile);


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SaveDualVerts(const std::string &outfile);


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SaveDualFaces(const std::string &outfile);


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SaveDualCells(const std::string &outfile);


    /**
     * @brief 
     * 
     * @param outfile 
     */
    void SaveDualMesh(const std::string &outfile);


    /**
     * @brief 
     * 
     * @return const Adjacency& 
     */
    inline const Adjacency & PrimalMap() const { return this->primal_map_; }
    
    
    /**
     * @brief 
     * 
     * @return const std::vector<Tetrahedron>& 
     */
    inline const std::vector<Tetrahedron> & PrimalCells() const { return this->primal_cells_; }


    /**
     * @brief 
     * 
     * @return const std::vector<Triangle>& 
     */
    inline const std::vector<Triangle> & PrimalFaces() const { return this->primal_faces_; }


    /**
     * @brief 
     * 
     * @return const std::vector<Edge>& 
     */
    inline const std::vector<Edge> & PrimalEdges() const { return this->primal_edges_; }


    inline const std::vector<Vertex> & PrimalVerts() const { return this->primal_verts_; }


    /**
     * @brief 
     * 
     * @return const std::vector<Vertex>& 
     */
    inline const std::vector<Vertex> & DualVerts() const { return this->dual_verts_; }


    /**
     * @brief 
     * 
     * @return const std::vector<Polygon>& 
     */
    inline const std::vector<Polygon> & DualFaces() const { return this->dual_faces_; }


    /**
     * @brief 
     * 
     * @return const std::vector<Polyhedron>& 
     */
    inline const std::vector<Polyhedron> & DualCells() const { return this->dual_cells_; }


};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_MESHING_TET_TO_POLY_HPP_ 
