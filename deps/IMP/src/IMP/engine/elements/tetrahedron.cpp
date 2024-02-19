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


#include "IMP/engine/elements/tetrahedron.hpp"

namespace IMP {


Vertex Tetrahedron::Centroid(const std::vector<Vertex> &verts) const
{
    return (verts[this->v0_]+verts[this->v1_]+verts[this->v2_]+verts[this->v3_]) / 4.;
}


bool Tetrahedron::IsSharingVertex(const Tetrahedron &tet) const
{
    // Sort the connectivity of tet.
    std::array<int,4> tet_conn = {tet.V0(), tet.V1(), tet.V2(), tet.V3()};
    std::sort(std::begin(tet_conn), std::end(tet_conn));

    // Find number of sharing vertices.
    int sharing = 0;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v0_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v1_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v2_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v3_)) sharing++;

    if (sharing == 1) return true;
    return false;
}


bool Tetrahedron::IsSharingEdge(const Tetrahedron &tet) const
{
    // Sort the connectivity of tet.
    std::array<int,4> tet_conn = {tet.V0(), tet.V1(), tet.V2(), tet.V3()};
    std::sort(std::begin(tet_conn), std::end(tet_conn));

    // Find number of sharing vertices.
    int sharing = 0;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v0_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v1_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v2_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v3_)) sharing++;

    // The tetrahedra share an edge if they share two vertices.
    if (sharing == 2) return true;
    return false;
}


bool Tetrahedron::IsSharingFace(const Tetrahedron &tet) const
{
    // Sort the connectivity of tet.
    std::array<int,4> tet_conn = {tet.V0(), tet.V1(), tet.V2(), tet.V3()};
    std::sort(std::begin(tet_conn), std::end(tet_conn));

    // Find number of sharing vertices.
    int sharing = 0;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v0_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v1_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v2_)) sharing++;
    if (std::binary_search(std::begin(tet_conn), std::end(tet_conn), this->v3_)) sharing++;

    // The tetrahedra share a face if they share three vertices.
    if (sharing == 3) return true;
    return false;
}


double Tetrahedron::SignedVolume(const std::vector<Vertex> &verts) const
{
    // Edges starting from firts vertex.
    Vertex e1 = verts[this->v1_] - verts[this->v0_];
    Vertex e2 = verts[this->v2_] - verts[this->v0_];
    Vertex e3 = verts[this->v3_] - verts[this->v0_];

    // Compute volume of tetrahedron.
    return (Det(e1,e2,e3)/6.);

}


double Tetrahedron::AbsVolume(const std::vector<Vertex> &verts) const
{
    return std::abs(this->SignedVolume(verts));
}


double Tetrahedron::Volume(const Vec<3, double> &pa, const Vec<3, double> &pb, const Vec<3, double> &pc, const Vec<3, double> &pd) const
{
    auto DiffOfProds = [] (double a, double b, double c, double d) {
        auto cd = c * d;
        auto err = std::fma(-c, d, cd);
        auto dop = std::fma(a, b, -cd);
        return dop + err;
    };

    // Determinant computation
    auto adx = pa[0] - pd[0];
    auto bdx = pb[0] - pd[0];
    auto cdx = pc[0] - pd[0];
    auto ady = pa[1] - pd[1];
    auto bdy = pb[1] - pd[1];
    auto cdy = pc[1] - pd[1];
    auto adz = pa[2] - pd[2];
    auto bdz = pb[2] - pd[2];
    auto cdz = pc[2] - pd[2];
    auto Det = ( adx * DiffOfProds(bdy, cdz, bdz, cdy) +
        bdx * DiffOfProds(cdy, adz, cdz, ady) +
        cdx * DiffOfProds(ady, bdz, adz, bdy) );

    return std::abs(Det) / 6.;

}


bool Tetrahedron::IsValid(const std::vector<Vertex> &verts) const
{
    if (!(this->SignedVolume(verts) > 0.)) return false;
    return true;
}


} // End of namespace IMP.
