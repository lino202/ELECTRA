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
   \file vertex.hpp
   \brief Vertex class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2020
*/

#ifndef IMP_ENGINE_ELEMENTS_VERTEX_HPP_
#define IMP_ENGINE_ELEMENTS_VERTEX_HPP_

#include <cmath>

namespace IMP {

/** \addtogroup Meshing \{ */


/**
 * \class Vertex
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a Vertex for mesh generation and processing.
 */
class Vertex {

private:

    double x_;              /**< The x coordinate of the vertex. */

    double y_;              /**< The y coordinate of the vertex. */

    double z_;              /**< The z coordinate of the vertex. */

    bool on_boundary_;      /**< Boolean to determine a boundary vertex. */

    bool on_feature_edge_;

    bool on_feature_corner_;    /**< Boolean to determine a vertex at the corner of a feuture edge. */

public:

    /**
     * \brief The Vertex constructor.
     */
    inline Vertex() : x_(0.), y_(0.), z_(0.), on_boundary_(false), on_feature_edge_(false), on_feature_corner_(false) {}


    /**
     * \brief The Vertex constructor.
     */
    inline Vertex(double x, double y=0., double z=0., bool on_boundary=false, bool on_feature_edge=false, bool on_feature_corner=false) : 
            x_(x), y_(y), z_(z), on_boundary_(on_boundary), on_feature_edge_(on_feature_edge), on_feature_corner_(on_feature_corner) {}


    /**
     * \brief The Vertex destructor.
     */
    inline virtual ~Vertex() {}

    
    /**
     * \brief Set the coordinates of the vertex.
     * \param [in] x 
     * \param [in] y 
     * \param [in] z
     * \return [void] 
     */
    inline void SetCoords(double x, double y=0., double z=0.) { this->x_=x; this->y_=y; this->z_=z; }


    /**
     * \brief Set the On Boundary object
     * \param [in] on_boundary 
     * \return [void]
     */
    inline void SetOnBoundary(bool on_boundary) { this->on_boundary_ = on_boundary; }


    /**
     * \brief Set the On Feature Corner object
     * \param [in] on_feature_corner
     * \return [void] 
     */
    inline void SetOnFeatureEdge(bool on_feature_edge) { this->on_feature_edge_ = on_feature_edge; }


    /**
     * \brief Set the On Feature Corner object
     * \param [in] on_feature_corner
     * \return [void] 
     */
    inline void SetOnFeatureCorner(bool on_feature_corner) { this->on_feature_corner_ = on_feature_corner; }


    /**
     * @brief 
     * 
     * @return double 
     */
    double Norm() const;


    inline Vertex & operator += (double val) {
        this->x_ += val;
        this->y_ += val;
        this->z_ += val;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator += (const Vertex &vert) {
        this->x_ += vert.x_;
        this->y_ += vert.y_;
        this->z_ += vert.z_;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator -= (double val) {
        this->x_ -= val;
        this->y_ -= val;
        this->z_ -= val;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator -= (const Vertex &vert) {
        this->x_ -= vert.x_;
        this->y_ -= vert.y_;
        this->z_ -= vert.z_;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator *= (double val) {
        this->x_ *= val;
        this->y_ *= val;
        this->z_ *= val;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator *= (const Vertex &vert) {
        this->x_ *= vert.x_;
        this->y_ *= vert.y_;
        this->z_ *= vert.z_;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator /= (double val) {
        this->x_ /= val;
        this->y_ /= val;
        this->z_ /= val;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    inline Vertex & operator /= (const Vertex &vert) {
        this->x_ /= vert.x_;
        this->y_ /= vert.y_;
        this->z_ /= vert.z_;
        this->on_boundary_ = false;
        this->on_feature_edge_ = false;
        this->on_feature_corner_ = false;
        
        return *this;
    }


    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator + (const Vertex &vert1, const Vertex &vert2);


    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator - (const Vertex &vert1, const Vertex &vert2);


    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator * (const Vertex &vert1, const Vertex &vert2);


    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator * (const Vertex &vert, double val);


    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator * (double val, const Vertex &vert);
    
    
    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator / (const Vertex &vert1, const Vertex &vert2);


    /**
     * @brief 
     * 
     * @param vert1 
     * @param vert2 
     * @return Vertex 
     */
    friend Vertex operator / (const Vertex &vert, double val);


    /**
     * @brief 
     * 
     * @return double 
     */
    inline double X() const { return this->x_; }


    /**
     * @brief 
     * 
     * @return double 
     */
    inline double Y() const { return this->y_; }


    /**
     * @brief 
     * 
     * @return double 
     */
    inline double Z() const { return this->z_; }
    

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
    inline bool OnFeatureEdge() const { return this->on_feature_edge_; }


    /**
     * @brief 
     * 
     * @return true 
     * @return false 
     */
    inline bool OnFeatureCorner() const { return this->on_feature_corner_; }

};


/**
 * \brief Computes the cross product of two vertices.
 * \param v1 
 * \param v2 
 * \return Vertex 
 */
Vertex Cross(const Vertex &v1, const Vertex &v2);


/**
 * \brief Computes the inner product of two vertices.
 * \param v1 
 * \param v2 
 * \return Vertex 
 */
double Dot(const Vertex &v1, const Vertex &v2);


/**
 * \brief Computes the determinant for two vertices.
 * \param v1 
 * \param v2 
 * \return double 
 */
double Det(const Vertex &v1, const Vertex &v2);


/**
 * \brief Computes the determinant for three vertices.
 * \param v1 
 * \param v2 
 * \param v3 
 * \return double 
 */
double Det(const Vertex &v1, const Vertex &v2, const Vertex &v3);


/** \} End of Doxygen Groups*/

} // End of namespace IMP.


#endif // IMP_ENGINE_ELEMENTS_VERTEX_HPP_ 
