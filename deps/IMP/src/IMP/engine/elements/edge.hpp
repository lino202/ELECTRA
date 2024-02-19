/*
 * IMP. Image to Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos A. Mounedgs> <konstantinos.mounedgs@gmail.com>
 *
 * This program is free software: you can redisedgbute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is disedgbuted in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


/**
   \file edge.hpp
   \brief Edge class header file.
   \author Konstantinos A. Mounedgs
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_EDGE_HPP_
#define IMP_ENGINE_ELEMENTS_EDGE_HPP_

#include "IMP/engine/elements/vertex.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace IMP {

/** \addtogroup Meshing \{ */

/**
 * \class Edge
 * \author Konstantinos A. Mounedgs
 * \brief Class implemmenting an Edge element.
 */
class Edge {

private:

    int v0_;                  /**< The index of the first vertex of the edge. */

    int v1_;                  /**< The index of the second vertex of the edge. */

    int parent_face_id_;      /**< The index of the parent face of the edge. */

    bool on_boundary_;        /**< Boolean to determine a boundary edge. */

    bool on_feature_;         /**< Boolean to determine an edge on a feuture edge. */

public:

    /**
     * \brief The Edge constructor.
     */
    inline Edge() : v0_(-1), v1_(-1), parent_face_id_(-1), on_boundary_(false), on_feature_(false) {}


    /**
     * \brief The Edge constructor.
     */
    inline Edge(int v0, int v1, int parent_face_id, bool on_boundary=false, bool on_feature=false) : 
            v0_(v0), v1_(v1), parent_face_id_(parent_face_id), on_boundary_(on_boundary), on_feature_(on_feature) {}

    
    /**
     * @brief Construct a new Edge object
     * 
     * @param edg 
     */
    inline Edge(const Edge &edg) : 
            v0_(edg.v0_), v1_(edg.v1_), parent_face_id_(edg.parent_face_id_), on_boundary_(edg.on_boundary_), on_feature_(edg.on_feature_) {}


    /**
     * \brief The Edge destructor.
     */
    inline virtual ~Edge() {}

    
    /**
     * \brief Set the connectivity of the edge.
     * \param [in] n1
     * \param [in] n2
     * \return [void] 
     */
    inline void SetVertIds(int v0, int v1) { this->v0_=v0; this->v1_=v1; }


    /**
     * \brief Set the index of the parent face of the edge.
     * \param [in] parent_cell_id 
     * \return [void] 
     */
    inline void SetParentFaceId(int parent_face_id) { this->parent_face_id_ = parent_face_id; }


    /**
     * \brief Set the On Boundary object
     * \param [in] on_boundary 
     * \return [void]
     */
    inline void SetOnBoundary(bool on_boundary) { this->on_boundary_ = on_boundary; }


    /**
     * \brief Set the On Feature edge object
     * \param [in] on_feature_edge
     * \return [void] 
     */
    inline void SetOnFeature(bool on_feature) { this->on_feature_ = on_feature; }


    /**
     * \brief Flip the connectivity of the edge. 
     * \return [void]
     */
    inline void Flip() { int temp = this->v0_; this->v0_ = this->v1_; this->v1_ = temp; }


    /**
     * @brief Sorts the connectivity of the edge in ascending order.
     */
    inline void SortVertIds() { if (this->v0_ > this->v1_) this->Flip(); }


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
    bool IsSharingVertex(const Edge &edg) const;


    /**
     * @brief 
     * 
     */
    struct HashFunc {
        size_t operator() (const Edge& edg) const {
            size_t v0_hash = std::hash<int>()(edg.V0());
            size_t v1_hash = std::hash<int>()(edg.V1()) << 1;
            return v0_hash ^ v1_hash;
        }
    };


    /**
     * @brief 
     * 
     * @param edg1 
     * @param edg2 
     * @return true 
     * @return false 
     */
    friend bool operator == (const Edge &edg1, const Edge &edg2);

    
    /**
     * @brief 
     * 
     * @param edg1 
     * @param edg2 
     * @return true 
     * @return false 
     */
    friend bool operator != (const Edge &edg1, const Edge &edg2);


    /**
     * @brief 
     * 
     * @param edg1 
     * @param edg2 
     * @return true 
     * @return false 
     */
    friend bool operator < (const Edge &edg1, const Edge &edg2);

    
    /**
     * @brief 
     * 
     * @param edg1 
     * @param edg2 
     * @return true 
     * @return false 
     */
    friend bool operator > (const Edge &edg1, const Edge &edg2);


    /**
     * @brief 
     * 
     * @param edg1 
     * @param edg2 
     * @return true 
     * @return false 
     */
    friend bool operator <= (const Edge &edg1, const Edge &edg2);


    /**
     * @brief 
     * 
     * @param edg1 
     * @param edg2 
     * @return true 
     * @return false 
     */
    friend bool operator >= (const Edge &edg1, const Edge &edg2);


    /**
     * @brief 
     * 
     * @param edg 
     * @return Edge& 
     */
    Edge & operator = (const Edge &edg)
    {
        this->v0_ = edg.v0_;
        this->v1_ = edg.v1_;
        this->parent_face_id_ = edg.parent_face_id_;
        this->on_boundary_ = edg.on_boundary_;
        this->on_feature_ = edg.on_feature_;
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
     * @return int 
     */
    inline int ParentFaceId() const { return this->parent_face_id_; }
    

    /**
     * @brief 
     * 
     * @return true 
     * @return false 
     */
    inline bool OnBoundary() const { return this->on_boundary_; }


    /**
     * @brief 
     * 
     * @return true 
     * @return false 
     */
    inline bool OnFeature() const { return this->on_feature_; }

};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_ELEMENTS_EDGE_HPP_ 
