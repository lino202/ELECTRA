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
   \file polygon.hpp
   \brief Polygon class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_POLYGON_HPP_
#define IMP_ENGINE_ELEMENTS_POLYGON_HPP_

#include "IMP/engine/elements/vertex.hpp"
#include <vector>

namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class Polygon
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a Polygon element.
 */
class Polygon {

private:

    std::vector<int> vert_ids_;

public:

    /**
     * @brief Construct a new Polygon object
     * 
     */
    inline Polygon() : vert_ids_() {}


    /**
     * @brief Construct a new Polygon object
     * 
     */
    inline Polygon(const std::vector<int> &vert_ids) : vert_ids_(vert_ids) {}


    /**
     * @brief Construct a new Polygon object
     * 
     * @param poly 
     */
    inline Polygon(const Polygon &poly) : vert_ids_(poly.vert_ids_) {}


    /**
     * @brief Destroy the Polygon object
     * 
     */
    inline virtual ~Polygon() {}


     /**
     * @brief Construct a new Centroid object
     * 
     * @param verts 
     */
    Vertex Centroid(const std::vector<Vertex> &verts) const;


    /**
     * @brief 
     * 
     * @param poly1 
     * @param poly2 
     * @return true 
     * @return false 
     */
    friend bool operator == (const Polygon &poly1, const Polygon &poly2);
    

    /**
     * @brief 
     * 
     * @param poly1 
     * @param poly2 
     * @return true 
     * @return false 
     */
    friend bool operator != (const Polygon &poly1, const Polygon &poly2);


     /**
     * @brief 
     * 
     * @param poly 
     * @return Triangle& 
     */
    Polygon & operator = (const Polygon &poly)
    {   
        if (*this != poly)
            this->vert_ids_ = poly.vert_ids_;
        return *this;
    }


    /**
     * @brief 
     * 
     * @param i 
     */
    inline void AddInVertIds(int i) { this->vert_ids_.emplace_back(i); }


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    const std::vector<int> & VertIds() const { return this->vert_ids_; }


    /**
     * @brief 
     * 
     * @param id 
     * @return int 
     */
    int V(std::size_t id) const { return this->vert_ids_[id]; }


};


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_ELEMENTS_POLYGON_HPP_ 
