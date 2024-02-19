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


#include "IMP/engine/elements/half_facet.hpp"


namespace IMP {

HalfFacet::HalfFacet(int vertex_id, int cell_id, short facet_id) :
    vertex_id_(vertex_id), cell_id_(cell_id), facet_id_(facet_id)
{}


HalfFacet::HalfFacet() : HalfFacet(-1, -1, -1)
{ /* Set not initialized values to -1. */ }


HalfFacet::~HalfFacet()
{}


void HalfFacet::Set(int vertex_id, int cell_id, short facet_id)
{
    // Check if possitive vertex id was given.
    if (vertex_id < 0) {
        throw std::invalid_argument(Logger::Error("Half-facet's vertex id must be positive."));
    }

    // Check if possitive cell id was given.
    if (cell_id < 0) {
        throw std::invalid_argument(Logger::Error("Half-facet's cell id must be positive."));
    }

    // Check if possitive facet id was given.
    if (facet_id < 0) {
        throw std::invalid_argument(Logger::Error("Half-facet's facet id must be positive."));
    }

    // Set half-facet's data.
    this->vertex_id_ = vertex_id;
    this->cell_id_ = cell_id;
    this->facet_id_ = facet_id;
}


void HalfFacet::Set(const HalfFacet &other)
{
    // Copy half-facet's data from other half-facet.
    this->vertex_id_ = other.vertex_id_;
    this->cell_id_ = other.cell_id_;
    this->facet_id_ = other.facet_id_;
}


bool HalfFacet::IsNull() const
{
    // Check if the half-facet's data are not initialized.
    return (this->vertex_id_ == -1 && this->cell_id_ == -1 && this->facet_id_ == -1);
}


bool operator == (const HalfFacet &hf1, const HalfFacet &hf2)
{
    // Check half-facets for equality.
    return (hf1.vertex_id_ == hf2.vertex_id_ &&
            hf1.cell_id_ == hf2.cell_id_ &&
            hf1.facet_id_ == hf2.facet_id_);

}


std::ostream & operator << (std::ostream &out, const HalfFacet &hf)
{
    // Output the half-facet's data.
    hf.IsNull() ? out << "-1: <-,->" :
                out << hf.vertex_id_ << ": <" << hf.cell_id_ << "," << hf.facet_id_ << ">";

    return out;
}



} //end of namespace IMP
