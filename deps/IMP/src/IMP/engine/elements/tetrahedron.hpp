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
   \file tetrahedron.hpp
   \brief Tetrahedron class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_TETRAHEDRON_HPP_
#define IMP_ENGINE_ELEMENTS_TETRAHEDRON_HPP_

#include "IMP/engine/elements/vertex.hpp"
#include "IMP/engine/vectors/vec.hpp"

#include <vector>
#include <array>
#include <algorithm>

namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class Tetrahedron
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting an Tetrahedron for mesh generation and processing.
 */
class Tetrahedron {

private:

    int v0_;                  /**< The index of the first vertex of the triangle. */

    int v1_;                  /**< The index of the second vertex of the triangle. */

    int v2_;                  /**< The index of the third vertex of the triangle. */

    int v3_;                  /**< The index of the fourth vertex of the triangle. */

public:

    /**
     * \brief The Tetrahedron constructor.
     */
    inline Tetrahedron() : v0_(-1), v1_(-1), v2_(-1), v3_(-1) {}


    /**
     * \brief The Tetrahedron constructor.
     */
    inline Tetrahedron(int v0, int v1, int v2, int v3) : 
            v0_(v0), v1_(v1), v2_(v2), v3_(v3) {}

    
    /**
     * \brief The Tetrahedron constructor.
     */
    inline Tetrahedron(const Tetrahedron &tet) : 
            v0_(tet.v0_), v1_(tet.v1_), v2_(tet.v2_), v3_(tet.v3_) {}


    /**
     * \brief The Tetrahedron destructor.
     */
    inline virtual ~Tetrahedron() {}

    
    /**
     * \brief Set the connectivity of the edge.
     * \param [in] v0
     * \param [in] v1
     * \param [in] v2
     * \param [in] v3
     * \return [void] 
     */
    inline void SetVertIds(int v0, int v1, int v2, int v3) { this->v0_=v0; this->v1_=v1; this->v2_=v2; this->v3_=v3; }


    /**
     * @brief 
     * 
     * @param verts 
     * @return Vertex 
     */
    Vertex Centroid(const std::vector<Vertex> &verts) const;


    /**
     * @brief Compute the volume of the tetrahedron.
     * 
     * @param verts 
     * @return double 
     */
    double SignedVolume(const std::vector<Vertex> &verts) const;


    /**
     * @brief Compute the volume of the tetrahedron.
     * 
     * @param verts 
     * @return double 
     */
    double AbsVolume(const std::vector<Vertex> &verts) const;


    double Volume(const Vec<3, double> &pa, const Vec<3, double> &pb, const Vec<3, double> &pc, const Vec<3, double> &pd) const;



    /**
     * @brief Check if the tetrahedron is valid.
     * In order to be valid it should be star-shaped meaning that it must have positive volume.
     * @param verts 
     * @return true 
     * @return false 
     */
    bool IsValid(const std::vector<Vertex> &verts) const;


    /**
     * @brief 
     * 
     * @param tet 
     * @return true 
     * @return false 
     */
    bool IsSharingVertex(const Tetrahedron &tet) const;


    /**
     * @brief 
     * 
     * @param tet 
     * @return true 
     * @return false 
     */
    bool IsSharingEdge(const Tetrahedron &tet) const;


    /**
     * @brief 
     * 
     * @param tet 
     * @return true 
     * @return false 
     */
    bool IsSharingFace(const Tetrahedron &tet) const;    


    /**
     * @brief 
     * 
     * @param tet 
     * @return Tetrahedron& 
     */
    Tetrahedron & operator = (const Tetrahedron &tet)
    {
        this->v0_ = tet.v0_;
        this->v1_ = tet.v1_;
        this->v2_ = tet.v2_;
        this->v3_ = tet.v3_;
        return *this;
    }


    /**
     * @brief 
     * 
     * @return int 
     */
    inline int V0() const { return this->v0_; }


    /**
     * @brief 
     * 
     * @return double 
     */
    inline int V1() const { return this->v1_; }


    /**
     * @brief 
     * 
     * @return double 
     */
    inline int V2() const { return this->v2_; }


    /**
     * @brief 
     * 
     * @return double 
     */
    inline int V3() const { return this->v3_; }


};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_ELEMENTS_TETRAHEDRON_HPP_ 
