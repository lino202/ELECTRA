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

#include "IMP/engine/elements/poly_facet.hpp"

namespace IMP {


PolyFacet::PolyFacet(const std::vector<int> &connectivity, int parent_cell_id, int neigh_cell_id) :
    connectivity_(connectivity), parent_cell_id_(parent_cell_id), neigh_cell_id_(neigh_cell_id)
{}


PolyFacet::PolyFacet() : connectivity_(), parent_cell_id_(-1), neigh_cell_id_(-1)
{}


PolyFacet::~PolyFacet()
{}


void PolyFacet::Set(const std::vector<int> &connectivity, int parent_cell_id, int neigh_cell_id)
{
    this->connectivity_ = connectivity;
    this->parent_cell_id_ = parent_cell_id;
    this->neigh_cell_id_ = neigh_cell_id;
}


void PolyFacet::Set(const PolyFacet &other)
{
    this->connectivity_ = other.connectivity_;
    this->parent_cell_id_ = other.parent_cell_id_;
    this->neigh_cell_id_ = other.neigh_cell_id_;
}


bool PolyFacet::IsFree() const
{
    if (this->neigh_cell_id_ == -1) { return true; }
    return false;
}


double PolyFacet::Length(const std::vector<Vec<1, double>> &points_coords) const
{
    if (this->connectivity_.size() != 1) {
        std::string error_str = "Could not compute the length of the polyfacet. The polyfacet mush have one points.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Compute the length of the polyfacet.
    return 0.;
}


double PolyFacet::Length(const std::vector<Vec<2, double>> &points_coords) const
{
    if (this->connectivity_.size() != 2) {
        std::string error_str = "Could not compute the length of the polyfacet. The polyfacet mush have two points.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Extract the coordinates of the points of the polyfacet.
    Vec<2, double> v1 = points_coords[this->connectivity_[0]];
    Vec<2, double> v2 = points_coords[this->connectivity_[1]];

    // Compute the length of the polyfacet.
    return std::sqrt(v2.Distance2(v1));
}


double PolyFacet::Length(const std::vector<Vec<3, double>> &points_coords) const
{
    if (this->connectivity_.size() != 3) {
        std::string error_str = "Could not compute the length of the polyfacet. The polyfacet mush have three points.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Compute the length of the polyfacet.
    return 0.;
}


Vec<1, double> PolyFacet::Centroid(const std::vector<Vec<1, double>> &point_coords) const 
{
    // Compute mass center of the cell.
    Vec<1, double> centroid;
    for (const auto &id : this->connectivity_) centroid += point_coords[id];
    return (centroid /= this->connectivity_.size());
}


Vec<2, double> PolyFacet::Centroid(const std::vector<Vec<2, double>> &point_coords) const 
{
    // Compute mass center of the cell.
    Vec<2, double> centroid;
    for (const auto &id : this->connectivity_) centroid += point_coords[id];
    return (centroid /= this->connectivity_.size());
}


Vec<3, double> PolyFacet::Centroid(const std::vector<Vec<3, double>> &point_coords) const 
{
    // Compute mass center of the cell.
    Vec<3, double> centroid;
    for (const auto &id : this->connectivity_) centroid += point_coords[id];
    return (centroid /= this->connectivity_.size());
}


Vec<1, double> PolyFacet::Normal(const std::vector<Vec<1, double>> &points_coords) const
{
    if (this->connectivity_.size() != 1) {
        std::string error_str = "Could not compute the normal of the polyfacet. The polyfacet mush have one point.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Compute the normal of the polyfacet.
    Vec<1, double> norm;
    return norm;
}


Vec<2, double> PolyFacet::Normal(const std::vector<Vec<2, double>> &points_coords) const
{
    if (this->connectivity_.size() != 2) {
        std::string error_str = "Could not compute the normal of the polyfacet. The polyfacet mush have two points.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Extract the coordinates of the points of the polyfacet.
    Vec<2, double> v1 = points_coords[this->connectivity_[0]];
    Vec<2, double> v2 = points_coords[this->connectivity_[1]];

    // Compute the normal of the polyfacet.
    Vec<2, double> norm({v1[1]-v2[1], v2[0]-v1[0]});
    return norm / std::sqrt(norm.Length2());
}


Vec<3, double> PolyFacet::Normal(const std::vector<Vec<3, double>> &points_coords) const
{
    if (this->connectivity_.size() != 3) {
        std::string error_str = "Could not compute the normal of the polyfacet. The polyfacet mush have three points.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Compute the normal of the polyfacet.
    Vec<3, double> norm;
    return norm;
}


bool PolyFacet::operator < (const PolyFacet &pf) const
{
    return (this->connectivity_ < pf.Connectivity());
}


bool operator == (const PolyFacet &pf1, const PolyFacet &pf2)
{
    // For equality it is sufficient that the two facets have the same point indices.
    return (pf1.Connectivity() == pf2.Connectivity());
}


std::ostream & operator << (std::ostream &out, const PolyFacet &pf)
{

    for (const auto &p : pf.Connectivity()) {
        out << p << " ";
    }
    out << pf.ParentCellId() << " " << pf.NeighCellId();
    return  out;
}


} //end of namespace IMP