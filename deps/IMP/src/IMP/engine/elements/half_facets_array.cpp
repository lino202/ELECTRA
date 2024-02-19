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


#include "IMP/engine/elements/half_facets_array.hpp"



namespace IMP {


HalfFacetsArray::HalfFacetsArray() : half_facets_(), is_sorted_(false)
{}


HalfFacetsArray::~HalfFacetsArray()
{}


void HalfFacetsArray::Append(int vertex_id, int cell_id, short facet_id)
{

    // Construct a half facet with the given information and append it in the array.
    this->half_facets_.emplace_back(HalfFacet(vertex_id, cell_id, facet_id));

    // The stored half facets probably are not any more sorted, so the is_sorted flag is turned to false.
    this->is_sorted_ = false;
}


void HalfFacetsArray::Append(const HalfFacet &half_facet)
{
    // Append the half facet in the array.
    this->half_facets_.emplace_back(half_facet);

    // The stored half facets probably are not any more sorted, so the is_sorted flag is turned to false.
    this->is_sorted_ = false;
}


void HalfFacetsArray::Sort()
{
    // Ascending order by vertex id -> cell id -> half facet id.
    auto ascend_order = [](const HalfFacet &hf1, const HalfFacet &hf2) {
        if (hf1.VertexId() == hf2.VertexId()) {
            if (hf1.CellId() == hf2.CellId())
                return hf1.FacetId() < hf2.FacetId();
            else
                return hf1.CellId() < hf2.CellId();
        }
        else
            return hf1.VertexId() < hf2.VertexId();
    };

    // Sort the mapped half-facets in the ascending order.
    std::sort(this->half_facets_.begin(), this->half_facets_.end(), ascend_order);

    // Set sorted flag to true.
    this->is_sorted_ = true;
}


std::pair<std::vector<HalfFacet>::const_iterator,
std::vector<HalfFacet>::const_iterator> HalfFacetsArray::AllAttachedToVertex(int vertex_id) const
{
    // Find range of half facets with mapped vertex index equal to the mapped vertex index of the temp_hf half facet.
    // The last index points out of the range.
    return std::equal_range(this->half_facets_.begin(), this->half_facets_.end(), HalfFacet(vertex_id, -1, -1),
                                      [](const HalfFacet &hf1, const HalfFacet &hf2) { return  hf1.VertexId() < hf2.VertexId(); });
}


std::vector<HalfFacet>::const_iterator HalfFacetsArray::FirstAttachedToVertex(int vertex_id) const
{
    // Get the first half facet in the range.
    return std::lower_bound(this->half_facets_.begin(), this->half_facets_.end(), HalfFacet(vertex_id, -1, -1),
                                      [](const HalfFacet &hf1, const HalfFacet &hf2) { return  hf1.VertexId() < hf2.VertexId(); });
}


std::vector<int> HalfFacetsArray::MappedVerticesIds() const
{
    // Create empty unique vertices vector.
    std::vector<int> unique_vertices;
    unique_vertices.reserve(this->half_facets_.size());

    // Fill with all the vertice indices (possibly with duplicate entries).
    for (const auto &hf : this->half_facets_) { unique_vertices.emplace_back(hf.VertexId()); }

    // Sort unique vertices vector if necessary.
    if (!this->is_sorted_) {
        std::cout << Logger::Warning("MappedVerticesIds method requested on unsorted HalfFacetsArray."
                                     "Time overhead is expected due to temporary sorting.") << std::endl;
        // Perform sorting.
        std::sort(unique_vertices.begin(), unique_vertices.end());
    }

    // Now erase the duplicate entries from unique vertices vector.
    unique_vertices.erase(std::unique(unique_vertices.begin(), unique_vertices.end()), unique_vertices.end());

    return unique_vertices;

}



} // End of namespace IMP.
