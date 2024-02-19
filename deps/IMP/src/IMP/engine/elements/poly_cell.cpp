/*
 * IMP. Image and Mesh Processing library.
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



#include "IMP/engine/elements/poly_cell.hpp"


namespace IMP {

PolyCell::PolyCell() : connectivity_(), shape_(CellShape::poly)
{}


PolyCell::PolyCell(const PolyCell &pc) : connectivity_(pc.connectivity_), shape_(pc.shape_)
{}


PolyCell::~PolyCell()
{}


void PolyCell::SetConnectivity(const std::vector<int> &ids)
{
    // Set the points coordinates.
    this->connectivity_ = ids;
}


void PolyCell::AddInConnectivity(int id)
{
    this->connectivity_.emplace_back(id);
}


double PolyCell::Measure(const std::vector<Vec<1, double>> &points_coords, const std::vector<PolyFacet> &facets) const
{
    auto pid1 = facets[this->connectivity_[0]].C(0);
    auto pid2 = facets[this->connectivity_[1]].C(0);
    return std::sqrt(points_coords[pid1].Distance2(points_coords[pid2]));
}


double PolyCell::Measure(const std::vector<Vec<2,double>> &points_coords, const std::vector<PolyFacet> &facets) const
{
    // The polygon is degenerate.
    if (this->connectivity_.size() < 2) {
        auto err_msg = "Could not measure degenerate cell.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Get the facets of the cell.
    auto cell_facets = std::vector<PolyFacet>{};
    cell_facets.reserve(this->connectivity_.size());
    for (const auto &fid : this->connectivity_) {
        cell_facets.emplace_back(facets[fid]);
    }

    // Circular order of points in the cell.
    auto point_ids = std::vector<int>{cell_facets[0].C(0), cell_facets[0].C(1)};
    cell_facets.erase(cell_facets.begin());
    while (point_ids.size() != this->connectivity_.size()) {
        auto i = int{0}, j = int{-1};
        for (const auto &facet : cell_facets) {
            if (facet.C(0) == point_ids.back()) {
                j = 1;
                break;
            }
            if (facet.C(1) == point_ids.back()) {
                j = 0;
                break;
            }
            i++;
        }
        if (j != -1) {
            point_ids.emplace_back(cell_facets[i].C(j));
            cell_facets.erase(cell_facets.begin()+i);
        } else {
            auto err_msg = "Could not measure cell. Circular connection of cell points failed.";
            throw std::runtime_error(Logger::Error(err_msg));
        }
    }

    // Compute the area of the polygon.
    auto area = double{0.};
    auto j = std::size_t{0};
    auto p1 = std::size_t{0}, p2 = std::size_t{0};
    for (auto i = std::size_t(0); i != point_ids.size(); ++i) {
        j = (i+1) % (point_ids.size());
        p1 = point_ids[i];
        p2 = point_ids[j];
        area += (points_coords[p1][0] * points_coords[p2][1] - points_coords[p1][1] * points_coords[p2][0]);
    }

    return 0.5*std::abs(area);

    // std::size_t n = point_ids.size();
    // int id = 0;
    // std::vector<Vec<2, double>> V(n+1);
    // for (const auto &pi : point_ids) {
    //     V[id++] = points_coords[pi];
    // }
    // V[n] = points_coords[point_ids[0]];

    // // Compute the area of the polygon.
    // double area = 0.;
    // std::size_t  i, j, k;
    // for (i=1, j=2, k=0; i != n; ++i, ++j, ++k) {
    //     area += V[i][0] * (V[j][1] - V[k][1]);
    // }
    // area += V[n][0] * (V[1][1] - V[n-1][1]);  // Wrap-around term.

    // return 0.5*std::abs(area);
}


double PolyCell::Measure(const std::vector<Vec<3, double>> &points_coords, const std::vector<PolyFacet> &facets) const
{

    auto DiffOfProds = [] (double a, double b, double c, double d) {
        auto cd = c * d;
        auto err = std::fma(-c, d, cd);
        auto dop = std::fma(a, b, -cd);
        return dop + err;
    };

    // Lambda for determinant computation
    auto Det = [DiffOfProds] (const Vec<3, double> &pa, const Vec<3, double> &pb, const Vec<3, double> &pc, const Vec<3, double> &pd) {
        auto adx = pa[0] - pd[0];
        auto bdx = pb[0] - pd[0];
        auto cdx = pc[0] - pd[0];
        auto ady = pa[1] - pd[1];
        auto bdy = pb[1] - pd[1];
        auto cdy = pc[1] - pd[1];
        auto adz = pa[2] - pd[2];
        auto bdz = pb[2] - pd[2];
        auto cdz = pc[2] - pd[2];
        return ( adx * DiffOfProds(bdy, cdz, bdz, cdy) +
                 bdx * DiffOfProds(cdy, adz, cdz, ady) +
                 cdx * DiffOfProds(ady, bdz, adz, bdy) );
    };

    // Compute cell centroid.
    auto cell_ctr = Vec<3, double>{};
    auto cnt = int{0};
    for (const auto &fid : this->connectivity_) {
        for (const auto &pid : facets[fid].Connectivity()) {
            cell_ctr += points_coords[pid];
            cnt++;
        }
    }
    cell_ctr /= static_cast<double>(cnt);

    // Compute cell volume via symmetric decomposition.
    auto cell_volume = double{0.};
    auto pa = Vec<3, double>{};
    auto pb = Vec<3, double>{};
    auto face_ctr = Vec<3, double>{};
    auto nn = int{0}, n1 = int{0}, n2 = int{0};
    for (const auto &fid : this->connectivity_) {
        face_ctr = facets[fid].Centroid(points_coords);

        // Number of points in the facet.
        nn = static_cast<int>(facets[fid].Connectivity().size());

        n2 = facets[fid].C(nn-1);
        for (int i = 0; i != nn; ++i) {
            // Connectivity of current edge of the facet.
            n1 = n2;  n2 = facets[fid].C(i);

            // Add volume of partial tetrahedron.
            pa = points_coords[n1];
            pb = points_coords[n2];
            cell_volume += std::abs(Det(pa, pb, face_ctr, cell_ctr)) / 6.;
        }
    }
    return cell_volume;

}


Vec<1, double> PolyCell::Centroid(const std::vector<Vec<1, double>> &point_coords, const std::vector<PolyFacet> &facets) const
{
    // Compute mass center of the cell.
    Vec<1, double> centroid;
    auto counter = int{0};
    for (const auto &fid : this->connectivity_) {
        for (const auto &id : facets[fid].Connectivity()) {
            centroid += point_coords[id];
            counter++;
        }
    }
    return (centroid /= counter);
}


Vec<2, double> PolyCell::Centroid(const std::vector<Vec<2, double>> &point_coords, const std::vector<PolyFacet> &facets) const
{
    // Compute mass center of the cell.
    Vec<2, double> centroid;
    auto counter = int{0};
    for (const auto &fid : this->connectivity_) {
        for (const auto &id : facets[fid].Connectivity()) {
            centroid += point_coords[id];
            counter++;
        }
    }
    return (centroid /= counter);
}


Vec<3, double> PolyCell::Centroid(const std::vector<Vec<3, double>> &point_coords, const std::vector<PolyFacet> &facets) const
{
    // Compute mass center of the cell.
    Vec<3, double> centroid;
    auto counter = int{0};
    for (const auto &fid : this->connectivity_) {
        for (const auto &id : facets[fid].Connectivity()) {
            centroid += point_coords[id];
            counter++;
        }
    }
    return (centroid /= counter);
}


bool operator == (const PolyCell &pc1, const PolyCell &pc2)
{
    // Check for equality.
    return (pc1.Connectivity() == pc2.Connectivity() && pc1.Shape() == pc2.Shape());
}


bool operator != (const PolyCell &pc1, const PolyCell &pc2)
{
    // Check for equality.
    return !(pc1 == pc2);
}


std::ostream & operator << (std::ostream &out, const PolyCell &pc)
{

    for (const auto &id : pc.Connectivity()) {
        out << id << "\n";
    }
    return  out;
}


PolyCell & PolyCell::operator = (const PolyCell &pc)
{
    if (*this != pc) this->connectivity_ = pc.connectivity_;
    return *this;
}


PolyCell & PolyCell::operator = (PolyCell &&pc) noexcept
{
    if (*this != pc) this->connectivity_ = std::move(pc.connectivity_);
    return *this;
}


} // End of namespace IMP
