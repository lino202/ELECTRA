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


#include "IMP/engine/elements/triangle.hpp"

namespace IMP {


Vertex Triangle::Normal(const std::vector<Vertex> &verts) const
{
    // Compute the normal vector of the triangle.
    Vertex normal = Cross(verts[this->v1_]-verts[this->v0_], verts[this->v2_]-verts[this->v0_]);
    double norm = normal.Norm();
    
    return normal /= norm;
}


Vertex Triangle::Centroid(const std::vector<Vertex> &verts) const
{
    return (verts[this->v0_]+verts[this->v1_]+verts[this->v2_]) / 3.;
}


bool Triangle::IsSharingVertex(const Triangle &tri) const
{   
    // Sort the connectivity of tri.
    std::array<int,3> tri_conn = {tri.V0(), tri.V1(), tri.V2()};       
    std::sort(std::begin(tri_conn), std::end(tri_conn));

    // Find number of sharing vertices.
    int sharing = 0;
    if (std::binary_search(std::begin(tri_conn), std::end(tri_conn), this->v0_)) sharing++;
    if (std::binary_search(std::begin(tri_conn), std::end(tri_conn), this->v1_)) sharing++;
    if (std::binary_search(std::begin(tri_conn), std::end(tri_conn), this->v2_)) sharing++;

    if (sharing == 1) return true;
    return false;
}


bool Triangle::IsSharingEdge(const Triangle &tri) const
{   
    // Sort the connectivity of tri.
    std::array<int,3> tri_conn = {tri.V0(), tri.V1(), tri.V2()};       
    std::sort(std::begin(tri_conn), std::end(tri_conn));

    // Find number of sharing vertices.
    int sharing = 0;
    if (std::binary_search(std::begin(tri_conn), std::end(tri_conn), this->v0_)) sharing++;
    if (std::binary_search(std::begin(tri_conn), std::end(tri_conn), this->v1_)) sharing++;
    if (std::binary_search(std::begin(tri_conn), std::end(tri_conn), this->v2_)) sharing++;

    // The triangles share an edge if they share two vertices.
    if (sharing == 2) return true;
    return false;
}


bool Triangle::HasEdge(const Edge &edg) const
{
    // Sort the connectivity of this.
    std::array<int,3> conn = {this->v0_, this->v1_, this->v2_};       
    std::sort(std::begin(conn), std::end(conn));

    // Find both vertices of the edge.
    return (std::binary_search(std::begin(conn), std::end(conn), edg.V0()) &&
            std::binary_search(std::begin(conn), std::end(conn), edg.V1()) );
}


bool operator == (const Triangle &tri1, const Triangle &tri2)
{
    return (tri1.V0() == tri2.V0() && tri1.V1() == tri2.V1() && tri1.V2() == tri2.V2() );    
}

    
bool operator != (const Triangle &tri1, const Triangle &tri2) { return !(tri1 == tri2); }


bool operator < (const Triangle &tri1, const Triangle &tri2)
{
    if (tri1.V0() != tri2.V0()) return tri1.V0() < tri2.V0();
    if (tri1.V1() != tri2.V1()) return tri1.V1() < tri2.V1();
    
    return tri1.V2() < tri2.V2();
}

    
bool operator > (const Triangle &tri1, const Triangle &tri2) { return tri2 < tri1; }


bool operator <= (const Triangle &tri1, const Triangle &tri2) { return !(tri1 > tri2); }


bool operator >= (const Triangle &tri1, const Triangle &tri2) { return !(tri1 < tri2); }


} // End of namespace IMP.