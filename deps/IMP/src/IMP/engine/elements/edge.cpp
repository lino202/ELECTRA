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


#include "IMP/engine/elements/edge.hpp"

namespace IMP {


Vertex Edge::Centroid(const std::vector<Vertex> &verts) const
{
    return 0.5*(verts[this->v0_]+verts[this->v1_]);
}


bool Edge::IsSharingVertex(const Edge &edg) const
{   
    // Sort the connectivity of tri.
    std::array<int,2> edg_conn = {edg.V0(), edg.V1()};       
    std::sort(std::begin(edg_conn), std::end(edg_conn));

    // Find number of sharing vertices.
    int sharing = 0;
    if (std::binary_search(std::begin(edg_conn), std::end(edg_conn), this->v0_)) sharing++;
    if (std::binary_search(std::begin(edg_conn), std::end(edg_conn), this->v1_)) sharing++;

    if (sharing == 1) return true;
    return false;
}


bool operator == (const Edge &edg1, const Edge &edg2)
{
    return (edg1.V0() == edg2.V0() && edg1.V1() == edg2.V1() );    
}

    
bool operator != (const Edge &edg1, const Edge &edg2) { return !(edg1 == edg2); }


bool operator < (const Edge &edg1, const Edge &edg2)
{
    if (edg1.V0() != edg2.V0()) return edg1.V0() < edg2.V0();
    
    return edg1.V1() < edg2.V1();
}

    
bool operator > (const Edge &edg1, const Edge &edg2) { return edg2 < edg1; }


bool operator <= (const Edge &edg1, const Edge &edg2) { return !(edg1 > edg2); }


bool operator >= (const Edge &edg1, const Edge &edg2) { return !(edg1 < edg2); }


} // End of namespace IMP.