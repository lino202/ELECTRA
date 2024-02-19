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

#ifndef IMP_ENGINE_ELEMENTS_CELL_TPP_
#define IMP_ENGINE_ELEMENTS_CELL_TPP_


#include "IMP/engine/elements/cell.hpp"


namespace IMP {

template <short DIM, short CELL_VERTS>
Cell<DIM, CELL_VERTS>::Cell() : sib_half_facets_(CELL_VERTS), connectivity_(), shape_(CellShape::unknown)
{}


template <short DIM, short CELL_VERTS>
Cell<DIM, CELL_VERTS>::~Cell()
{}


template <short DIM, short CELL_VERTS>
void Cell<DIM, CELL_VERTS>::SetShape(CellShape shape)
{
    // Check if any inconsistency between the given cell type and the cell's dimensions exist.
    if ((shape == CellShape::tri || shape == CellShape::quad) && DIM < 2) {
        throw std::invalid_argument(Logger::Error("Can not set the cell's shape. The cell has lower dimensions from the given shape."));
    }
    else if ((shape == CellShape::tet || shape == CellShape::hex) && DIM < 3) {
        throw std::invalid_argument(Logger::Error("Can not set the cell's shape. The cell has lower dimensions from the given shape."));
    }

    // Set the cell's shape.
    this->shape_ = shape;
}


template <short DIM, short CELL_VERTS>
void Cell<DIM, CELL_VERTS>::SetConnectivity(const Vec<CELL_VERTS, int> &connectivity)
{
    // Set the cell's connectivity.
    this->connectivity_ = connectivity;
}


template <short DIM, short CELL_VERTS>
void Cell<DIM, CELL_VERTS>::SetConnectivity(std::initializer_list<int> connectivity)
{
    // Check size of connectivity list.
    if (connectivity.size() != CELL_VERTS) {
        throw std::invalid_argument(Logger::Error("Cannot set connectivity of AhfCell. The connectivity "
                                                  "list's size should be equal to the number of the cell's vertices."));
    }

    // Set connectivity.
    this->connectivity_.Set(connectivity);
}


template <short DIM, short CELL_VERTS>
void Cell<DIM, CELL_VERTS>::SetAllSibHfacets(const std::vector<HalfFacet> &half_facets)
{
    // Check size of half_facets container.
    if (half_facets.size() != CELL_VERTS) {
        throw std::invalid_argument(Logger::Error("Cannot set all sibling half-facets. The given half-facets container's size "
                                                  "should be equal to the number of the corresponding cell's nodes."));
    }

    // Set the cell's half-facets.
    this->sib_half_facets_ = half_facets;
}


template <short DIM, short CELL_VERTS>
void Cell<DIM, CELL_VERTS>::SetSibHfacet(short local_node_id, const HalfFacet &half_facet)
{
    // Check the validity of the local_vert_id
    if (local_node_id < 0 || local_node_id > CELL_VERTS-1) {
        throw std::invalid_argument(Logger::Error("Cannot set sibling half-facet. The given local node id must be "
                                                  "in the range [0, CELL_VERTICES_NUM-1)."));
    }

    // Set the half-facet of the cell's node with the given local id.
    this->sib_half_facets_[local_node_id].Set(half_facet);
}


template <short DIM, short CELL_VERTS>
Vec<(CELL_VERTS-1), int> Cell<DIM, CELL_VERTS>::FacetConnectivity(short half_facet_id) const
{
    // Check the validity of the half_facet_id
    if (half_facet_id < 0 || half_facet_id > CELL_VERTS-1) {
        throw std::invalid_argument(Logger::Error("Cannot get cell's half-facet connectivity. The given half-facet id must be "
                                                  "in the range [0, CELL_VERTICES_NUM-1)."));
    }

    // Container with indices of the half-facet's connectivity nodes.
    Vec<CELL_VERTS-1, int> facet_conn;
    short pos = 0;
    short local_id = 0;
    for (const auto &node_id : this->connectivity_) {

        local_id = &node_id - &this->connectivity_[0];

        if (local_id != half_facet_id) { facet_conn[pos] = node_id; pos++; }
    }

    return facet_conn;

}


template <short DIM, short CELL_VERTS>
std::vector<int> Cell<DIM, CELL_VERTS>::AdjacentNodeIdsTo(int node_id, short hfacet_id) const
{
    // The indices of the adjacent nodes to the querry node.
    std::vector<int> adjacent_node_ids;
    adjacent_node_ids.resize(CELL_VERTS-2);

    short pos = 0;
    
    // If CELL_VERTS = 2 (cell=edge) the half-facet is a node. Adjacency is not defined.
    if (CELL_VERTS > 2) {

        // Get the ids of the adjacent nodes to the querry node belonging to the same half-facet.
        // The node with id equal to hfacet_id is considered opposite to the half facet.
        for (short i = 0; i != hfacet_id; ++i) {
            if (this->connectivity_[i] != node_id) {
                adjacent_node_ids[pos] = this->connectivity_[i];
                pos++;
            }
        }
        for (short i = hfacet_id+1; i != CELL_VERTS; ++i) {
            if (this->connectivity_[i] != node_id) {
                adjacent_node_ids[pos] = this->connectivity_[i];
                pos++;
            }
        }
    }

    // Sort the adjacent nodes container.
    std::sort(adjacent_node_ids.begin(), adjacent_node_ids.end());

    return adjacent_node_ids;
}


template <short DIM, short CELL_VERTS>
int Cell<DIM, CELL_VERTS>::MaxNodeIdInHfacet(short hfacet_id) const
{

    // Initialize to the smallest possible value.
    int max_node_id = 0;

    // Iterate indices before hfacet_id.
    for (short i = 0; i != hfacet_id; ++i) {
        if (this->connectivity_[i] > max_node_id)  { max_node_id = this->connectivity_[i]; }
    }

    // Iterate indices after hfacet_id.
    for (short i = hfacet_id+1; i != CELL_VERTS; ++i) {
        if (this->connectivity_[i] > max_node_id)  { max_node_id = this->connectivity_[i]; }
    }

    // Return the largest global node index that belongs to the facet.
    return max_node_id;
}


template <short DIM, short CELL_VERTS>
inline short Cell<DIM, CELL_VERTS>::LocalNodeIdOf(int global_node_id) const
{
    // Find the local id of the node if it belongs in the cell.
    for (const auto & node_id : this->connectivity_) {
        
        short local_node_id = &node_id - &this->connectivity_[0];
        
        if (node_id == global_node_id) { return local_node_id; }
    }

    // Default local id value for a node not in the cell.
    return -1;
}


template <short DIM, short CELL_VERTS>
inline std::vector<short> Cell<DIM, CELL_VERTS>::HfacetIdsOnLocalNode(short local_node_id) const
{
    // All the half-facets are connected to the vertex except the opposite one.
    std::vector<short> hfacet_ids;  
    hfacet_ids.resize(CELL_VERTS-1);

    short pos = 0;
    // Iterate indices before local_node_id.
    for (short i = 0; i != local_node_id; ++i) {
        hfacet_ids[pos] = i;
        pos++;
    }

    // Iterate indices after local_node_id.
    for (short i = local_node_id+1; i != CELL_VERTS; ++i) {
        hfacet_ids[pos] = i;
        pos++;
    }

    return hfacet_ids;
}


template <short DIM, short CELL_VERTS>
std::vector<short> Cell<DIM, CELL_VERTS>::LocalNodeIdsInHfacet(short half_facet_id) const
{
    // All the nodes belong to the half-facet except the one opposite of it.
    std::vector<short> node_ids;  
    node_ids.resize(CELL_VERTS-1);

    short pos = 0;

    // Iterate indices before half_facet_id.
    for (short i = 0; i != half_facet_id; ++i) {
        node_ids[pos] = i;
        pos++;
    }

    // Iterate indices after half_facet_id.
    for (short i = half_facet_id+1; i != CELL_VERTS; ++i) {
        node_ids[pos] = i;
        pos++;
    }

    return node_ids;
}


template <short DIM, short CELL_VERTS>
inline bool Cell<DIM,CELL_VERTS>::IsConnectedToCell(int neigh_cell_id) const
{
    // True if any of the cell's half facets has the same id as the neighbor cell.
    for (const auto &hfacet : this->sib_half_facets_) {
        if (hfacet.CellId() == neigh_cell_id) { return true; }
    }

    // Otherwise false.
    return false;
}


template <short DIM, short CELL_VERTS>
bool operator == (const Cell<DIM, CELL_VERTS> &c1, const Cell<DIM, CELL_VERTS> &c2)
{
    // Check for equality.
    return (c1.Connectivity() == c2.Connectivity()) && (c1.SibHfacets() == c2.SibHfacets());
}


template <short DIM, short CELL_VERTS>
std::ostream & operator << (std::ostream &out, const Cell<DIM, CELL_VERTS> &cell)
{
    out << cell.Connectivity();
    return  out;
}



} // End of namespace IMP


#endif //IMP_ENGINE_ELEMENTS_CELL_TPP_
