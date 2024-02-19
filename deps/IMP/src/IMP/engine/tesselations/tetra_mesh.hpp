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
   \file tetramesh.hpp
   \brief Tetrahedron class header file.
   \author Konstantinos A. Mountris
   \date 22/11/2020
*/

#ifndef IMP_ENGINE_TESSELATIONS_TETRA_MESH_HPP_
#define IMP_ENGINE_TESSELATIONS_TETRA_MESH_HPP_

#include "IMP/engine/elements/vertex.hpp"
#include "IMP/engine/elements/tetrahedron.hpp"

#include <vector>
#include <algorithm>

namespace IMP {

/** \addtogroup Tesselations \{ */


/**
 * \class TetraMesh
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting an tetrahedral mesh.
 */
class TetraMesh {

private:

    std::vector<Vertex> verts_;

    std::vector<Tetrahedron> tetras_;

public:

    /**
     * \brief The TetraMesh constructor.
     */
    inline TetraMesh() : verts_(), tetras_() {}


    /**
     * \brief The TetraMesh constructor.
     */
    inline TetraMesh(const std::vector<Vertex> &verts, std::vector<Tetrahedron> &tetras) : 
            verts_(verts), tetras_(tetras) {}

    
    /**
     * \brief The TetraMesh constructor.
     */
    inline TetraMesh(const TetraMesh &tetramesh) : 
            verts_(tetramesh.verts_), tetras_(tetramesh.tetras_) {}


    /**
     * @brief Destroy the Tetra Mesh object
     * 
     */
    inline virtual ~TetraMesh() {}


    TetraMesh & operator = (const TetraMesh &tetramesh)
    {
        this->verts_ = tetramesh.verts_;
        this->tetras_ = tetramesh.tetras_;
        return *this;
    }


    /**
     * @brief 
     * 
     * @return const std::vector<Vertex>& 
     */
    const std::vector<Vertex> & Verts() const { return this->verts_; }


    const Vertex & Verts(std::size_t id) const { return this->verts_[id]; }


    const Vertex & VertsAt(std::size_t id) const { return this->verts_.at(id); }


    /**
     * @brief 
     * 
     * @return const std::vector<Tetrahedron>& 
     */
    const std::vector<Tetrahedron> & Tetras() const { return this->tetras_; }


    const Tetrahedron & Tetras(std::size_t id) const { return this->tetras_[id]; }


    const Tetrahedron & TetrasAt(std::size_t id) const { return this->tetras_.at(id); }


};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_TESSELATIONS_TETRA_MESH_HPP_ 
