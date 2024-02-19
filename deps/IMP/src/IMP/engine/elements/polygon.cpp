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


#include "IMP/engine/elements/polygon.hpp"

namespace IMP {


Vertex Polygon::Centroid(const std::vector<Vertex> &verts) const
{
    Vertex centroid_coords;
    for (const auto &v : this->vert_ids_) centroid_coords += verts[v];
    return centroid_coords / this->vert_ids_.size();
}


bool operator == (const Polygon &poly1, const Polygon &poly2)
{
    return (poly1.VertIds() == poly2.VertIds());    
}

    
bool operator != (const Polygon &poly1, const Polygon &poly2)
{ 
    return !(poly1 == poly2); 
}


} // End of namespace IMP.