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


#include "IMP/engine/elements/vertex.hpp"

namespace IMP {


double Vertex::Norm() const
{  
    return std::sqrt(this->x_*this->x_ + this->y_*this->y_ + this->z_*this->z_);
}


Vertex operator + (const Vertex &vert1, const Vertex &vert2) {
    return Vertex(vert1.x_+vert2.x_, vert1.y_+vert2.y_, vert1.z_+vert2.z_);
}


Vertex operator - (const Vertex &vert1, const Vertex &vert2) {
    return Vertex(vert1.x_-vert2.x_, vert1.y_-vert2.y_, vert1.z_-vert2.z_);
}


Vertex operator * (const Vertex &vert1, const Vertex &vert2) {
    return Vertex(vert1.x_*vert2.x_, vert1.y_*vert2.y_, vert1.z_*vert2.z_);
}


Vertex operator * (const Vertex &vert, double val) {
    return Vertex(vert.x_*val, vert.y_*val, vert.z_*val);
}


Vertex operator * (double val, const Vertex &vert) {
    return Vertex(vert.x_*val, vert.y_*val, vert.z_*val);
}


Vertex operator / (const Vertex &vert1, const Vertex &vert2) {
    return Vertex(vert1.x_/vert2.x_, vert1.y_/vert2.y_, vert1.z_/vert2.z_);
}


Vertex operator / (const Vertex &vert, double val) {
    return Vertex(vert.x_/val, vert.y_/val, vert.z_/val);
}


Vertex Cross(const Vertex &v1, const Vertex &v2)
{
    double x = v1.Y()*v2.Z() - v1.Z()*v2.Y();
    double y = v1.Z()*v2.X() - v1.X()*v2.Z();
    double z = v1.X()*v2.Y() - v1.Y()*v2.X();

    return Vertex(x,y,z);
}


double Dot(const Vertex &v1, const Vertex &v2)
{
    return v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
}


double Det(const Vertex &v1, const Vertex &v2, const Vertex &v3)
{
    return Dot(Cross(v1, v2),v3);
}



} // End of namespace IMP.
