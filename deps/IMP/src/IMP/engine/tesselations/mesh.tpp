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

#ifndef IMP_ENGINE_TESSELATIONS_MESH_TPP_
#define IMP_ENGINE_TESSELATIONS_MESH_TPP_

#include "IMP/engine/tesselations/mesh.hpp"


namespace IMP {


template <short DIM, short CELL_VERTS>
Mesh<DIM, CELL_VERTS>::Mesh() : nodes_(), cells_(), node_sets_(), nodes_to_half_facets_(), cells_shape_(CellShape::unknown),
                                are_sibling_half_facets_set_(false), are_nodes_mapped_to_half_facets_(false)
{}


template <short DIM, short CELL_VERTS>
Mesh<DIM, CELL_VERTS>::~Mesh()
{}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::LoadFrom(const std::string &mesh_filename, bool built_ahf)
{
    // Check if mesh filename is given.
    if (mesh_filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load mesh. The given mesh filename is empty."));
    }

    // Get the extension of the mesh filename.
    auto ext = mesh_filename.substr(mesh_filename.length()-4);

    // Clear the mesh containers.
    this->nodes_.clear();
    this->cells_.clear();
    this->node_sets_.clear();
    this->nodes_to_half_facets_.Clear();

    // Load mesh of given format.
    if (ext == ".inp") {
        AbaqusIO<DIM, CELL_VERTS> abaqus_io;
        abaqus_io.ReadMeshFrom(mesh_filename);
        abaqus_io.LoadNodesIn(this->nodes_);
        abaqus_io.LoadElementsIn(this->cells_);

        // Update the mesh cells type.
        this->UpdateCellsShape();

        // Load node sets if they are available.
        if (abaqus_io.NodeSetsExist()) { abaqus_io.LoadNodeSetsIn(this->node_sets_); }
    }
    else if (ext == ".vtk") {
        // Load from VTK Simple Legacy format.
        VtkIO<DIM, CELL_VERTS> vtk_io;
        vtk_io.ReadMeshFrom(mesh_filename);
        vtk_io.LoadNodesIn(this->nodes_);
        vtk_io.LoadElementsIn(this->cells_);

        // Update the mesh cells type.
        this->UpdateCellsShape();
    }
    else {
        std::string error_msg = Logger::Error("Could not load mesh of unknown format. Check given mesh filename: ") + mesh_filename;
        throw std::invalid_argument(error_msg);
    }

    // Set the Array-based HalfFacet structure of the mesh.
    if (built_ahf) {
        this->FindSiblingHalfFacets();
        this->MapNodesToHalfFacets();
    }
}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::SaveTo(const std::string &mesh_filename)
{
    // Check if mesh filename is given.
    if (mesh_filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load mesh. The given mesh filename is empty."));
    }

    // Get the extension of the mesh filename.
    auto ext = mesh_filename.substr(mesh_filename.length()-4);

    // Load mesh of given format.
    if (ext == ".inp") {
        AbaqusIO<DIM, CELL_VERTS> abaqus_io;
        abaqus_io.SaveMesh(this->nodes_, this->cells_, this->node_sets_, mesh_filename);
    } else {
        std::string error_msg = Logger::Error("Could not save mesh to unknown format. Check given mesh filename: " + mesh_filename);
        throw std::invalid_argument(error_msg);
    }
}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::Scale(double scale_value)
{
    for (auto &node : this->nodes_) {
        node *= scale_value;
    }
}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::AppendNode(std::initializer_list<double> coords)
{
    Vec<DIM, double> node;
    node.Set(coords);
    this->nodes_.emplace_back(node);
}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::AppendCell(std::initializer_list<int> connectivity)
{
    Cell<DIM, CELL_VERTS> c1;
    c1.SetConnectivity(connectivity);

    this->cells_.emplace_back(c1);
}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::AddNodeSet(const NodeSet &nodeset)
{
    this->node_sets_[nodeset.Name()] = nodeset;
}


template <short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::FindSiblingHalfFacets()
{
    // Array mapping each vertex to its incident
    // half-facets in which the vertex has the largest index.
    HalfFacetsArray v2hfs;
    v2hfs.Reserve(CELL_VERTS*this->CellsNum());

    // Iterate over the mesh cells.
    for (const auto &cell : this->cells_) {

        // The global id of the cell.
        auto cell_id = &cell - &this->cells_[0];

        // Construct the cell's half-facets.
        // The half-facets have the same id with the opposite vertex.
        for (short lf_id = 0; lf_id != CELL_VERTS; ++lf_id) {
            // Add the largest global-vertex & half-facet pair in the v2hfs container.
            v2hfs.Append(cell.MaxNodeIdInHfacet(lf_id), cell_id, lf_id);
        }

    }

    // Sort the v2hfs container.
    v2hfs.Sort();

    // Iterate over the mesh cells.
    for (auto &cell : this->cells_) {

        // The cell's index.
        auto cell_id = &cell - &this->cells_[0];

        // Iterate over the cell's half-facets.
        for (short cell_hf_id = 0; cell_hf_id != CELL_VERTS; ++cell_hf_id) {

            // Search for sibling half-facet if the cell's half-facet is null.
            if (cell.SibHfacets()[cell_hf_id].IsNull()) {

                // Get node with the largest global index in the cell's half-facet.
                auto max_vert_id = cell.MaxNodeIdInHfacet(cell_hf_id);

                // Find all the half-facets attached to the largest index vertex.
                auto attach_hfacets = v2hfs.AllAttachedToVertex(max_vert_id);

                // The number of the attached half-facets to the largest index vertex.
                auto attach_hfacets_num = std::distance(attach_hfacets.first, attach_hfacets.second);

                // Check if any attached half-facets for the max vertex of the cell's half-facet were found.
                if (attach_hfacets_num <= 0) {
                    std::string error_msg = "Could not extract sibling half-facets. No attached half-facets for vertex: " + std::to_string(max_vert_id);
                    throw std::runtime_error(Logger::Error(error_msg));
                }

                // Get adjacent vertices to the vertex with max_vert_id that belong to the cell_hf_id half-facet.
                std::vector<int> adjacents_nodes_in_cell_hf = cell.AdjacentNodeIdsTo(max_vert_id, cell_hf_id);

                // Iterate over the attached half facets to find the sibling to the cell's half-facet.
                for (auto hf_it = attach_hfacets.first; hf_it != attach_hfacets.second; ++hf_it) {

                    // Check only attached half-facets that do not belong to the cell with cell_id.
                    if ((*hf_it).CellId() != cell_id) {

                        // Get vertices in the attached half-facet that are adjacent to the vertex with max_vert_id.
                        std::vector<int> adjacent_nodes_to_attach_hf = this->cells_[(*hf_it).CellId()].AdjacentNodeIdsTo(max_vert_id, (*hf_it).FacetId());

                        // If adjacent nodes are the same for both half-facets
                        // then the sibling was found. The loop is stopped.
                        if (adjacent_nodes_to_attach_hf == adjacents_nodes_in_cell_hf) {
                            cell.SetSibHfacet(cell_hf_id, *hf_it);
                            break;
                        }
                    }
                }

            } // End Search for sibling half-facet if the cell's half-facet is null.
        } // End Iterate over the cell's half-facets.
    }  // End Iterate over the mesh cells.

    // Set status of sibling half facets set to true.
    this->are_sibling_half_facets_set_ = true;
}


template<short DIM, short CELL_VERTS>
void Mesh<DIM, CELL_VERTS>::MapNodesToHalfFacets()
{
    // Check if sibling half facets have been set.
    if (!this->are_sibling_half_facets_set_) {
        std::string error_msg = "Cannot map mesh nodes to half-facets of the mesh. The sibling half facets of the cells are not available";
        throw std::runtime_error(Logger::Error(error_msg));
    }


    // Ensure the vertex to half-facet map is empty.
    this->nodes_to_half_facets_.Clear();
    this->nodes_to_half_facets_.Reserve(this->NodesNum());

    // Create a temporary structure to mark visited half facets during vertex to half-facet mapping.
    struct MarkedCell {

        std::vector<bool> visited_hfacet;

        // Initialize the state of all visited half-facets to false.
        MarkedCell() : visited_hfacet(std::vector<bool>(CELL_VERTS, false)) {}

    } marked;

    // Create a list of marked cells for the mesh cells.
    std::vector<MarkedCell> marked_cells(this->cells_.size(), marked);

    // Create cyclic one-to-one mapping from local vertices ids to local facets ids.
    Vec<CELL_VERTS, short> local_vert_to_local_facet;
    if (CELL_VERTS == 2) {
        local_vert_to_local_facet[0] = 0;
        local_vert_to_local_facet[1] = 1;
    }
    else {
        for (short i = 0; i != CELL_VERTS; ++i) {
            // Each vertex maps to the next facet.
            local_vert_to_local_facet[i] = i+1;
        }
        // Last vertex maps to the first facet.
        local_vert_to_local_facet[CELL_VERTS-1] = 0;
    }

    // Iterate over the mesh cells to map vertices to half facets.
    for (const auto &cell : this->cells_) {

        auto cell_id = &cell - &this->cells_[0];

        // Iterate over the vertices of the current cell.
        for (const auto &local_vert : cell.Connectivity())// SmallIndType local_vi = 0; local_vi < (CELL_DIM+1); ++local_vi)
        {
            auto local_vert_id = &local_vert - &cell.Connectivity()[0];

            // Get the local id of the half-facet corresponding to the local_vert vertex.
            const short local_facet_id = local_vert_to_local_facet[local_vert_id];

            // Check if half-facet has been visited before.
            if (marked_cells[cell_id].visited_hfacet[local_facet_id] == false) {

                // Store not already visited half-facet.
                const int global_vert_id = cell.Connectivity()[local_vert_id]; // get global vertex
                this->nodes_to_half_facets_.Append(global_vert_id, cell_id, local_facet_id);

                // Mark the cell's half-facet as visited.
                marked_cells[cell_id].visited_hfacet[local_facet_id] = true;

                // Get the attached cells to the current vertex that are also sharing a half-facet.
                auto attached_cells_ids = AttachedCellIdsToNodeAndCell(global_vert_id, cell_id);

                // Mark the corresponding half-facets of the attached cells as visited.
                for (const auto & atch_cell_id : attached_cells_ids) {

                    // The local vertex index within the attached cell.
                    const short local_vert_id_attached = this->cells_[atch_cell_id].LocalNodeIdOf(global_vert_id);

                    // corresponding local facet index
                    const short local_facet_id_attached = local_vert_to_local_facet[local_vert_id_attached];
                    // mark it!
                    marked_cells[atch_cell_id].visited_hfacet[local_facet_id_attached] = true;
                }

            }

        } // End Iterate over the vertices of the current cell.

    } // End Iterate over the mesh cells to map vertices to half facets.


    // Sort the vertex to half-facet array.
    this->nodes_to_half_facets_.Sort();

    // Give priority to the boundary half-facets in the mapping.
    for (const auto &cell : this->cells_) {

        auto cell_id = &cell - &this->cells_[0];

        // Iterate over the cells half-facets.
        for (const auto &hfacet : cell.SibHfacets()) {

            auto hfacet_id = &hfacet - &cell.SibHfacets()[0];

            // Check if hfacet has no sibling (is boundary).
            if (hfacet.IsNull()) {

                // Get the vertices belonging to the half-facet.
                auto local_vert_ids = cell.LocalNodeIdsInHfacet(hfacet_id);

                // Iterate over the half-facets vertices.
                for (const auto &local_vert_id : local_vert_ids) {

                    // Get attached facets to current vertex.
                    int global_vert_id = cell.Connectivity()[local_vert_id];
                    auto attached_hfacets = this->nodes_to_half_facets_.AllAttachedToVertex(global_vert_id);

                    // Get the number of attached half facets to the current vertex.
                    std::size_t attached_num = std::distance(attached_hfacets.first, attached_hfacets.second);

                    if (attached_num > 1) {
                        for (auto att_hfacet = attached_hfacets.first; att_hfacet != attached_hfacets.second; ++att_hfacet) {

                            // Check if the half facets attached to the vertex belong also to an attached cell to the current one.
                            if (cell.IsConnectedToCell((*att_hfacet).CellId())) {

                                // Get the position of the current half facet to modify it in the vertex_to_half_facet container.
                                std::size_t mod_pos = att_hfacet - this->nodes_to_half_facets_.EditHalfFacets().begin();

                                // Modify the cell and half-facet indices of the corresponding half facet
                                this->nodes_to_half_facets_.EditHalfFacets()[mod_pos].Set((*att_hfacet).VertexId(), cell_id, hfacet_id);

                                // Break the loop since boundary facet has been assigned.
                                break;
                            }
                        }
                    }
                    else if (attached_num == 1) {
                        // Get the position of the unique half facet to modify in the vertex_to_half_facet container.
                        std::size_t mod_pos = attached_hfacets.first - this->nodes_to_half_facets_.HalfFacets().begin();

                        // Modify the cell and half-facet indices of the corresponding half facet
                        this->nodes_to_half_facets_.EditHalfFacets()[mod_pos].Set((*attached_hfacets.first).VertexId(), cell_id, hfacet_id);
                    }
                    else {
                        std::string error_msg = "Could not create Vertex to Half Facet mapping. No attached half-facets found for the vertex: " + std::to_string(global_vert_id);
                        throw std::runtime_error(Logger::Error(error_msg));
                    }
                } // End Iterate over the half-facets vertices.

            } // End Check if hfacet has no sibling (is boundary).

        } // End Iterate over the cells half-facets.

    } // End Give priority to the boundary half-facets in the mapping.

    // Set the state of the mapping of nodes to half-facets true.
    this->nodes_to_half_facets_.Sort();
    this->are_nodes_mapped_to_half_facets_ = true;
}


template<short DIM, short CELL_VERTS>
std::vector<int> Mesh<DIM, CELL_VERTS>::AttachedCellIdsToNode(int node_id) const
{
    // Get the attached half-facets to the node.
    auto attached_hfacets = this->nodes_to_half_facets_.AllAttachedToVertex(node_id);
    const auto attached_hfacets_num = std::distance(attached_hfacets.first, attached_hfacets.second);

    std::vector<int> attached_elems;
    attached_elems.reserve(5*attached_hfacets_num);

    for (auto it = attached_hfacets.first; it != attached_hfacets.second; ++it) {

        std::vector<int> temp_array = this->AttachedCellIdsToNodeAndCell(node_id, (*it).CellId());
        // store the found cells in cell_array
        const std::size_t old_size = attached_elems.size();
        attached_elems.resize(old_size + temp_array.size());
        std::copy(temp_array.begin(), temp_array.end(), attached_elems.begin()+old_size );
    }

    return attached_elems;
}


template<short DIM, short CELL_VERTS>
std::vector<int> Mesh<DIM, CELL_VERTS>::AttachedCellIdsToNodeAndCell(int node_id, int cell_id) const
{
    // The attached cells container.
    std::vector<int> attached_cells;

    // Perform a recursive search for the cells attached to the querry vertex and sharing a half-facet.
    std::function<void(int, int)> attached_search = [this, &attached_cells, &attached_search](int node_id, int cell_id) {

        // Search if the cell is already in the attached cells container.
        auto found_cell = std::find(attached_cells.begin(), attached_cells.end(), cell_id);

        // Stop if the cell was found in the container.
        if (found_cell != attached_cells.end()) { return; }

        // Search the local index of the vertex in the given cell.
        short local_node_id = this->cells_[cell_id].LocalNodeIdOf(node_id);

        // Stop if the vertex was not found in the cell.
        if (local_node_id == -1) { return;}

        // Add the cell's index in the container.
        if (attached_cells.size() == attached_cells.capacity()) {
            attached_cells.reserve(attached_cells.size() + 5*CELL_VERTS);
        }
        attached_cells.emplace_back(cell_id);

        // Get the indices of half-facets that share the vertex.
        auto hfacets_ids = this->cells_[cell_id].HfacetIdsOnLocalNode(local_node_id);

        // Search the attached cells for the vertex sharing half-facets.
        for (const auto &hfacet_id : hfacets_ids) {

            // Search the neighbor cell attached to this facet if it has been initialized.
            if (!this->cells_[cell_id].SibHfacets()[hfacet_id].IsNull())
            {
                attached_search(node_id, this->cells_[cell_id].SibHfacets()[hfacet_id].CellId());
            }
        }

    };
    attached_search(node_id, cell_id);

    attached_cells.shrink_to_fit();

    return attached_cells;
}


template<short DIM, short CELL_VERTS>
inline void Mesh<DIM, CELL_VERTS>::UpdateCellsShape()
{
    // Update the shape to the shape of the first cell.
    this->cells_shape_ = this->cells_[0].Shape();

    // Check if cells of different shape exist in the mesh.
    for (const auto &cell : this->cells_) {

        // Update to mixed shape if found a cell with different shape.
        if (cell.Shape() != this->cells_shape_) {
            this->cells_shape_ = CellShape::mixed;

            return;
        }
    }

}


template<short DIM, short CELL_VERTS>
inline std::vector<HalfFacet> Mesh<DIM, CELL_VERTS>::FreeFacets() const
{
    if (CELL_VERTS != DIM+1) {
        auto error_msg = "Can not compute free boundary. Currently available only for simplicial meshes.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize bdy container.
    auto bdy_facets = std::vector<HalfFacet>{};
    bdy_facets.reserve(CELL_VERTS*this->CellsNum());

    // Search all cells for free half facets.
    auto cid = int{0};
    auto free_hf = HalfFacet{};
    for (const auto &cell : this->cells_) {
        auto fid = int{0};
        for (const auto &cell_hf : cell.SibHfacets()) {
            if (cell_hf.IsNull()) {
                // This facet has no neighbor. It is free.
                free_hf.Set(0, cid, fid);

                // Store free facet.
                if (bdy_facets.size() >= bdy_facets.capacity())
                    bdy_facets.reserve(1.5*bdy_facets.capacity());
                bdy_facets.emplace_back(free_hf);
            }
            fid++;
        }
        cid++;
    }

    bdy_facets.shrink_to_fit();
    return bdy_facets;
}


template<short DIM, short CELL_VERTS>
inline std::vector<int> Mesh<DIM, CELL_VERTS>::FreeNodeIds() const
{
    // Get free half facets.
    auto bdy_facets = this->FreeFacets();

    // Extract node indices of free half facets.
    auto free_node_ids = std::vector<int>{};
    free_node_ids.reserve(CELL_VERTS*bdy_facets.size());
    auto facet_nids = Vec<(CELL_VERTS-1),int>{};
    for (const auto &hf : bdy_facets) {
        facet_nids = this->cells_[hf.CellId()].FacetConnectivity(hf.FacetId());
        for (const auto &id : facet_nids) {
            free_node_ids.emplace_back(id);
        }
    }

    free_node_ids.shrink_to_fit();
    return free_node_ids;
}


} // End of namespace IMP.

#endif // IMP_ENGINE_TESSELATIONS_MESH_TPP_
