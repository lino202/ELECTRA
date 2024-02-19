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
   \file triangle.hpp
   \brief Triangle class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_TRIANGLE_HPP_
#define IMP_ENGINE_ELEMENTS_TRIANGLE_HPP_

#include "IMP/engine/elements/vertex.hpp"
#include "IMP/engine/elements/edge.hpp"

#include <vector>
#include <array>
#include <algorithm>
#include <cstddef>
#include <functional>

namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class Triangle
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting an Triangle for mesh generation and processing.
 */
class Triangle {

private:

    int v0_;                  /**< The index of the first vertex of the triangle. */

    int v1_;                  /**< The index of the second vertex of the triangle. */

    int v2_;                  /**< The index of the third vertex of the triangle. */

    int parent_cell_id_;      /**< The index of the parent cell of the triangle. */         

    bool on_boundary_;        /**< Boolean to determine a boundary triangle. */

public:

    /**
     * \brief The Triangle constructor.
     */
    inline Triangle() : v0_(-1), v1_(-1), v2_(-1), parent_cell_id_(-1), on_boundary_(false) {}


    /**
     * \brief The Triangle constructor.
     */
    inline Triangle(int n1, int n2, int n3, int parent_cell_id=-1, bool on_boundary=false) : 
            v0_(n1), v1_(n2), v2_(n3), parent_cell_id_(parent_cell_id), on_boundary_(on_boundary) {}

    
    /**
     * @brief Construct a new Triangle object
     * 
     * @param tri 
     */
    inline Triangle(const Triangle &tri) : 
            v0_(tri.v0_), v1_(tri.v1_), v2_(tri.v2_), parent_cell_id_(tri.parent_cell_id_), on_boundary_(tri.on_boundary_) {}


    /**
     * \brief The Triangle destructor.
     */
    inline virtual ~Triangle() {}

    
    /**
     * \brief Set the connectivity of the edge.
     * \param [in] n1
     * \param [in] n2
     * \param [in] n3
     * \return [void] 
     */
    inline void SetVertIds(int v0, int v1, int v2) { this->v0_=v0; this->v1_=v1; this->v2_=v2; }


    /**
     * \brief Set the index of the parent cell of the triangle.
     * \param [in] parent_cell_id 
     * \return [void] 
     */
    inline void SetParentCellId(int parent_cell_id) { this->parent_cell_id_ = parent_cell_id; }


    /**
     * \brief Set the On Boundary object
     * \param [in] on_boundary 
     * \return [void]
     */
    inline void SetOnBoundary(bool on_boundary) { this->on_boundary_ = on_boundary; }


    /**
     * @brief Sorts the connectivity of the triangle in ascending order.
     * 
     * @param tri 
     * @return true 
     * @return false 
     */
    inline void SortVertIds() {
        std::array<int,3> tri = {this->v0_, this->v1_, this->v2_};       
        std::sort(std::begin(tri), std::end(tri));

        this->v0_ = tri[0];
        this->v1_ = tri[1];
        this->v2_ = tri[2];
    }


    /**
     * @brief 
     * 
     * @param verts 
     * @return Vertex 
     */
    Vertex Normal(const std::vector<Vertex> &verts) const;


    /**
     * @brief 
     * 
     * @param verts 
     * @return Vertex 
     */
    Vertex Centroid(const std::vector<Vertex> &verts) const;


    /**
     * @brief 
     * 
     * @param tri 
     * @return true 
     * @return false 
     */
    bool IsSharingVertex(const Triangle &tri) const;


    /**
     * @brief 
     * 
     * @param tri 
     * @return true 
     * @return false 
     */
    bool IsSharingEdge(const Triangle &tri) const;


    /**
     * @brief 
     * 
     * @param edg 
     * @return true 
     * @return false 
     */
    bool HasEdge(const Edge &edg) const;


    /**
     * @brief 
     * 
     */
    struct HashFunc {
        size_t operator() (const Triangle& tri) const {
            size_t v0_hash = std::hash<int>()(tri.V0());
            size_t v1_hash = std::hash<int>()(tri.V1()) << 1;
            size_t v2_hash = std::hash<int>()(tri.V2()) << 1;
            return v0_hash ^ v1_hash ^ v2_hash;
        }
    };


    /**
     * @brief 
     * 
     * @param tri1 
     * @param tri2 
     * @return true 
     * @return false 
     */
    friend bool operator == (const Triangle &tri1, const Triangle &tri2);

    
    /**
     * @brief 
     * 
     * @param tri1 
     * @param tri2 
     * @return true 
     * @return false 
     */
    friend bool operator != (const Triangle &tri1, const Triangle &tri2);


    /**
     * @brief 
     * 
     * @param tri1 
     * @param tri2 
     * @return true 
     * @return false 
     */
    friend bool operator < (const Triangle &tri1, const Triangle &tri2);

    
    /**
     * @brief 
     * 
     * @param tri1 
     * @param tri2 
     * @return true 
     * @return false 
     */
    friend bool operator > (const Triangle &tri1, const Triangle &tri2);


    /**
     * @brief 
     * 
     * @param tri1 
     * @param tri2 
     * @return true 
     * @return false 
     */
    friend bool operator <= (const Triangle &tri1, const Triangle &tri2);


    /**
     * @brief 
     * 
     * @param tri1 
     * @param tri2 
     * @return true 
     * @return false 
     */
    friend bool operator >= (const Triangle &tri1, const Triangle &tri2);


    /**
     * @brief 
     * 
     * @param tri 
     * @return Triangle& 
     */
    Triangle & operator = (const Triangle &tri)
    {
        this->v0_ = tri.v0_;
        this->v1_ = tri.v1_;
        this->v2_ = tri.v2_;
        this->parent_cell_id_ = tri.parent_cell_id_;
        this->on_boundary_ = tri.on_boundary_;
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
     * @return int 
     */
    inline int ParentCellId() const { return this->parent_cell_id_; }
    

    /**
     * @brief 
     * 
     * @return true 
     * @return false 
     */
    inline bool OnBoundary() const { return this->on_boundary_; }


};



/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_ELEMENTS_TRIANGLE_HPP_ 
