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
   \file vertex_optimizer.hpp
   \brief VertexOptimizer class header file.
   \author Konstantinos A. Mountris
   \date 30/11/2020
*/

#ifndef IMP_ENGINE_TESSELATIONS_VERTEX_OPTIMIZER_HPP_
#define IMP_ENGINE_TESSELATIONS_VERTEX_OPTIMIZER_HPP_

#include "IMP/engine/elements/vertex.hpp"

#include <Eigen/Dense>
 
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace IMP {

/** \addtogroup Tesselations \{ */


/**
 * \class VertexOptimizer
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting an optimizer for vertex positions.
 */
class VertexOptimizer {

private:

    std::vector<std::vector<Vertex>> trival_corners_;

    Vertex vert_;

    double neigh_volume_;


protected:

    double ComputeObjectiveFunction(const Vertex &v);


    void ComputeGradient(const Vertex &v, Eigen::Vector3d &grad);


    void ComputeHessian(const Eigen::Vector3d &grad, Eigen::Matrix3d &hess);

    
    bool IsPositiveDefinite(const Eigen::Matrix3d &hess);


    void MakePositiveDefinite(Eigen::Matrix3d &hess);


    bool SatisfiesArmijo(const Vertex &v, const Eigen::Vector3d &grad, const Eigen::Vector3d &p, double a);


public:

    /**
     * \brief The VertexOptimizer constructor.
     */
    inline VertexOptimizer() : trival_corners_(), vert_(), neigh_volume_(0.) {}


    /**
     * \brief The VertexOptimizer constructor.
     */
    inline VertexOptimizer(const std::vector<std::vector<Vertex>> &trival_corners, Vertex &vert, double neigh_volume) : 
            trival_corners_(trival_corners), vert_(vert), neigh_volume_(neigh_volume) {}

    
    /**
     * \brief The VertexOptimizer constructor.
     */
    inline VertexOptimizer(const VertexOptimizer &optimizer) : 
            trival_corners_(optimizer.trival_corners_), vert_(optimizer.vert_), neigh_volume_(optimizer.neigh_volume_) {}


    /**
     * @brief Destroy the Tetra Mesh object
     * 
     */
    inline virtual ~VertexOptimizer() {}


    /**
     * @brief 
     * 
     * @param vert 
     * @param trival_corners 
     */
    void Set(const Vertex &vert, const std::vector<std::vector<Vertex>> &trival_corners, double neigh_volume);


    void Optimize();


    inline const Vertex & OptimizedVert() const { return this->vert_; }


};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_TESSELATIONS_VERTEX_OPTIMIZER_HPP_ 
