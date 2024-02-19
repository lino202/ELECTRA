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
   \file polyhedron.hpp
   \brief Polyhedron class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_POLYHEDRON_HPP_
#define IMP_ENGINE_ELEMENTS_POLYHEDRON_HPP_

#include "IMP/engine/elements/vertex.hpp"
#include "IMP/engine/elements/edge.hpp"
#include "IMP/engine/elements/polygon.hpp"
#include "IMP/engine/elements/tetrahedron.hpp"
#include "IMP/engine/tesselations/tetra_mesh.hpp"

#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <cmath>

namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class Polyhedron
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a Polyhedron element.
 */
class Polyhedron {

private:

    std::vector<int> face_ids_;

public:

    /**
     * \brief The Polyhedron constructor.
     */
    inline Polyhedron() : face_ids_() {}


    /**
     * \brief The Polyhedron constructor.
     */
    inline Polyhedron(const std::vector<int> &face_ids) : face_ids_(face_ids) {}

    
    /**
     * \brief The Polyhedron constructor.
     */
    inline Polyhedron(const Polyhedron &poly) : face_ids_(poly.face_ids_) {}


    /**
     * \brief The Polyhedron destructor.
     */
    inline virtual ~Polyhedron() {}


    /**
     * @brief Construct a new Centroid object
     * 
     * @param verts 
     */
    Vertex Centroid(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const;


    /**
     * @brief Construct a symmetric decomposition of the polyhedron. 
     * @param verts 
     * @param faces 
     * @return Mesh<3,4> 
     */
    TetraMesh SymmetricDecomposition(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const;


    /**
     * @brief 
     * @param verts 
     * @param faces 
     * @return double 
     */
    double SignedVolume(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const;


    /**
     * @brief
     * @param verts 
     * @param faces 
     * @return double 
     */
    double AbsVolume(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const;


    /**
     * @brief Compute the Frobenious condition number of the polyhedron.
     * The computation is based on the Frobenious condition numbers of the trivalent corners in the symmetric decomposition. 
     * @param verts 
     * @param faces 
     * @return double 
     */
    double ConditionNumber(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const;


    /**
     * @brief Check if a polyhedron is valid.
     * In order to be valid it should be star-shaped meaning that all the tetrahedra of its symmetric
     * decomposition have positive volume.
     * \param verts The vertices of the mesh where the polyhedron belongs. 
     * \param faces The faces of the mesh where the polyhedron belongs.
     * \return true The polyhedron is valid a.k.a star-shaped.
     * \return false The polyhedron is not valid a.k.a not star-shaped.
     */
    bool IsValid(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const;


    /**
     * @brief 
     * 
     * @param poly1 
     * @param poly2 
     * @return true 
     * @return false 
     */
    friend bool operator == (const Polyhedron &poly1, const Polyhedron &poly2);
    

    /**
     * @brief 
     * 
     * @param poly1 
     * @param poly2 
     * @return true 
     * @return false 
     */
    friend bool operator != (const Polyhedron &poly1, const Polyhedron &poly2);


     /**
     * @brief 
     * 
     * @param poly 
     * @return Triangle& 
     */
    Polyhedron & operator = (const Polyhedron &poly)
    {   
        if (*this != poly)
            this->face_ids_ = poly.face_ids_;
        return *this;
    }


    /**
     * @brief 
     * @return const std::vector<Polygon>& 
     */
    const std::vector<int> & FaceIds() const { return face_ids_; }


    /**
     * @brief 
     * @param id 
     * @return const Polygon& 
     */
    int FaceIds(std::size_t id) const { return face_ids_[id]; }


    /**
     * @brief 
     * @return std::vector<Polygon>& 
     */
    std::vector<int> & EditFaceIds() { return face_ids_; }


    /**
     * @brief 
     * @param id 
     * @return Polygon& 
     */
    int & EditFaces(std::size_t id) { return face_ids_[id]; }

};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_ELEMENTS_POLYHEDRON_HPP_ 
