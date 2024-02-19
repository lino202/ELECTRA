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

#ifndef IMP_ENGINE_MESH_IO_ABAQUS_IO_TPP_
#define IMP_ENGINE_MESH_IO_ABAQUS_IO_TPP_

#include "IMP/engine/mesh_io/abaqus_io.hpp"

namespace IMP {

template <short DIM, short CELL_NODES>
AbaqusIO<DIM, CELL_NODES>::AbaqusIO() : parsed_mesh_(), mapped_node_ids_(), mapped_elem_ids_(), parts_startlines_(), nsets_startlines_(), 
                                        nodes_startline_(0), elems_startline_(0), parsed_cells_shape_(CellShape::unknown), parts_exist_(false), nsets_exist_(false)
{}


template <short DIM, short CELL_NODES>
AbaqusIO<DIM, CELL_NODES>::~AbaqusIO()
{}


template <short DIM, short CELL_NODES>
void AbaqusIO<DIM, CELL_NODES>::ReadMeshFrom(const std::string &mesh_filename)
{

    // Clear containers.
    this->parsed_mesh_.clear();
    this->parts_startlines_.clear();
    this->nsets_startlines_.clear();

    // Check if mesh filename is given.
    if (mesh_filename.empty()) {
        std::string error_msg = Logger::Error("Cannot read mesh from empty filename.");
        throw std::invalid_argument(error_msg);
    }

    // Check if mesh is in abaqus format (.inp).
    std::string ext = mesh_filename.substr(mesh_filename.length()-4);
    if (ext != ".inp") {
        std::string error_msg = Logger::Error("Expected mesh in Abaqus format (*.inp). Check mesh filename: ") + mesh_filename;
        throw std::invalid_argument(error_msg);
    }

    // Open mesh file.
    std::ifstream mesh_file(mesh_filename, std::ios::in);

    // Check if mesh file opened successfully.
    if (!mesh_file.is_open()) {
        std::string error_msg = Logger::Error("Could not open the mesh file. Check file path: ") + mesh_filename;
        throw std::runtime_error(error_msg);
    }

    // Available mesh information conditionals to default false values.
    bool nodes_exist = false;
    bool elems_exist = false;

    //Line of the file.
    std::string line = "";
    int line_number = 0;

    //Read mesh file line by line.
    while (std::getline(mesh_file, line)) {
        // Transform mesh file line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Read mesh file lines.
        this->parsed_mesh_.emplace_back(line);

        // Check if nodes list was found in the mesh.
        if (line.find("*node") != std::string::npos) {
            this->nodes_startline_ = line_number;
            nodes_exist = true;
        }

        // Check if elements list was found in the mesh.
        if (line.find("*element") != std::string::npos) {
            this->elems_startline_ = line_number;
            elems_exist = true;
        }

        // Check if elements list was found in the mesh.
        if (line.find("dc1d2") != std::string::npos) { this->parsed_cells_shape_ = CellShape::edge; }
        else if (line.find("b31") != std::string::npos) { this->parsed_cells_shape_ = CellShape::edge; }
        else if (line.find("cpe3") != std::string::npos) { this->parsed_cells_shape_ = CellShape::tri; }
        else if (line.find("cpe4") != std::string::npos) { this->parsed_cells_shape_ = CellShape::quad; }
        else if (line.find("s3r") != std::string::npos) { this->parsed_cells_shape_ = CellShape::tri; }
        else if (line.find("s4r") != std::string::npos) { this->parsed_cells_shape_ = CellShape::quad; }
        else if (line.find("c3d4") != std::string::npos) { this->parsed_cells_shape_ = CellShape::tet; }
        else if (line.find("c3d8") != std::string::npos) { this->parsed_cells_shape_ = CellShape::hex; }

        // Find the first lines of the partitions sets.
        if (line.find("*elset") != std::string::npos) {
            this->parts_startlines_.emplace_back(line_number);
            if (!this->parts_exist_) { this->parts_exist_ = true; }
        }

        // Find the first lines of the boundary node sets.
        if (line.find("*nset") != std::string::npos) {
            this->nsets_startlines_.emplace_back(line_number);
            if (!this->nsets_exist_) { this->nsets_exist_ = true; }
        }

        // Increase line number.
        line_number++;

    }

    // Close the mesh file.
    mesh_file.close();

    // Check if nodes and elements are available in the mesh file.
    if (!nodes_exist && !elems_exist) {
        std::string error_msg = Logger::Error("Something went wrong while reading nodes and elements. Check mesh file: ") + mesh_filename;
        throw std::runtime_error(error_msg);
    }

}


template <short DIM, short CELL_NODES>
void AbaqusIO<DIM, CELL_NODES>::SaveMesh(const std::vector<Vec<DIM,double>> &nodes, const std::vector<Cell<DIM,CELL_NODES>> &cells,
                                         const std::unordered_map<std::string,NodeSet> &node_sets, std::string mesh_filename)
{
    // Check if mesh filename is given.
    if (mesh_filename.empty()) {
        std::string error_msg = "Could not save the mesh. Check output mesh filename: " + mesh_filename;
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Check if given output filename is in abaqus format (.inp).
    if (mesh_filename.substr(mesh_filename.length() - 4) != ".inp") { mesh_filename += ".inp"; }

    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(mesh_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }


    //Open mesh output file.
    std::ofstream out(mesh_filename, std::ios::out | std::ios::trunc);
    if (!out.good()) {
        std::string error_msg = "Could not save the mesh to the output mesh file. Check file path: " + mesh_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Save header information.
    out << "*Heading\n";
    out << "** Job name: imp-mesh-job Model name: Model-1\n";
    out << "** Generated by: IMP library, echo=NO, model=NO, history=NO, contact=NO\n";
    out << "**\n";
    out << "** PARTS\n";
    out << "**\n";
    out << "*Part, name=Part-1\n";
    out << "*Node\n";

    // Save mesh nodes coordinates.
    for (auto &node : nodes) {
        auto node_id = &node - &nodes[0];
        out << node_id + 1 << ", ";
        for (short i = 0; i != DIM-1; ++i)
            out << node[i] << ", ";
        out << node[DIM-1] << "\n";
    }

    // Save mesh elements type info.
    if (DIM == 2 && CELL_NODES == 3) { out << "*Element, type=CPE3\n"; }
    else if (DIM == 2 && CELL_NODES == 4) { out << "*Element, type=CPE4\n"; }
    else if (DIM == 3 && CELL_NODES == 3) { out << "*Element, type=SFM3D3\n"; }
    else if (DIM == 3 && CELL_NODES == 4) { out << "*Element, type=C3D4\n"; }
    else if (DIM == 3 && CELL_NODES == 8) { out << "*Element, type=C3D8\n"; }
    else {
        std::string error_msg = "Could not write mesh of not supported element type. "
                                "Supported: 3-node Triangle, 4-node Quadrilateral, 4-node Tetrahedron, and 8-node Hexahedron.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Save mesh elements.
    for (auto &cell : cells) {
        auto cell_id = &cell - &cells[0];
        out << cell_id+1 << ", ";
        for (short i = 0; i != CELL_NODES-1; ++i)
            out << cell.Connectivity()[i]+1 << ", ";
        out << cell.Connectivity()[CELL_NODES-1]+1 << "\n";
    }

    // Save mesh nodesets.
    for (const auto &nset : node_sets) {
        out << "*Nset, nset=" << nset.first << "\n";
        int temp = 0;
        for (const auto &nset_node : nset.second.NodeIds()) {
            out << nset_node+1;
            temp++;
            if (temp == 20 || nset_node == nset.second.NodeIds().back()) {
                out << "\n";
                temp = 0;
            } else {
                out << ", ";
            }
        }
    }

    // Save the ending mesh header.
    out << "*End Part\n" << "**\n";
    out << "** ASSEMBLY\n" << "**\n" << "*Assembly, name=Assembly\n";
    out << "**\n" << "*Instance, name=Part-1-1, part=Part-1\n";
    out << "*End Instance\n" << "**\n" << "*End Assembly";

    // Close the output file.
    out.close();

}


template <short DIM, short CELL_NODES>
void AbaqusIO<DIM, CELL_NODES>::LoadNodesIn(std::vector<Vec<DIM, double>> &nodes)
{
    // Clean the mesh nodes container.
    nodes.clear();

    // Clean the mapped node ids container.
    this->mapped_node_ids_.clear();

    // The coordinates of a node in the mesh.
    Vec<DIM, double> coords;

    // The key id of a node in the parsed mesh.
    std::size_t key_id = 0;

    // Iterate though mesh starting from the nodes set starting line.
    // Skip the first line start from the first node's coordinates.
    std::string line = "";
    for (std::size_t it = this->nodes_startline_+1; it != parsed_mesh_.size(); ++it) {

        line = this->parsed_mesh_[it];

        // Replace comma occurence in line with space.
        std::replace(line.begin(), line.end(), ',', ' ');

        // Convert line to stringstream to be passed in the coordinates Vec.
        std::stringstream ss(line);

        // Get the coordinates until reach the end of the nodes list.
        if (!(ss >> key_id >> coords)) { break; }

        // Store the current node's coordinates.
        nodes.emplace_back(coords);

        // Store the mapping from the key id of the current node to the container's id.
        this->mapped_node_ids_[key_id] = it - (this->nodes_startline_+1);

    }

}


template <short DIM, short CELL_NODES>
void AbaqusIO<DIM, CELL_NODES>::LoadElementsIn(std::vector<Cell<DIM, CELL_NODES>> &cells)
{
    // Check if mapped node ids have been set.
    if (this->mapped_node_ids_.empty()) {
        throw std::runtime_error(Logger::Error("Could not load mesh elements. Load the nodes first."));
    }

    // Clean the cells container.
    cells.clear();

    // Clean the mapped element ids container.
    this->mapped_elem_ids_.clear();

    // Initialize an empty element and connectivity.
    Cell<DIM, CELL_NODES> cell;
    Vec<CELL_NODES, int> conn;

    // Set the elements shape.
    cell.SetShape(this->parsed_cells_shape_);

    // The element's key index in the parsed mesh.
    std::size_t key_id = 0;

    // Iterate though the mesh starting from the first line of the elements list.
    // Skip a line to read the first element's connectivity.
    std::string line = "";
    for (std::size_t it = this->elems_startline_+1; it != parsed_mesh_.size(); ++it) {

        line = this->parsed_mesh_[it];

        // Replace comma occurence in line with space.
        std::replace(line.begin(), line.end(), ',', ' ');

        // Convert line to stringstream to be passed in the connectivity variables.
        std::stringstream ss(line);

        // Get the key id and the connectivity of the current element.
        if (!(ss >> key_id >> conn)) { break; }

        // Update connectivity according to the mapped node ids.
        for (auto &id : conn) { id = this->mapped_node_ids_[id]; }

        // Set updated connectivity to the cell and store it.
        cell.SetConnectivity(conn);
        cells.emplace_back(cell);

        // Store the mapping from the key id of the current element to the container's id.
        this->mapped_elem_ids_[key_id] = it - (this->elems_startline_+1);

    }
    
}


template <short DIM, short CELL_NODES>
void AbaqusIO<DIM, CELL_NODES>::LoadNodeSetsIn(std::unordered_map<std::string, NodeSet> &node_sets)
{
    // Initialize the node_sets container.
    node_sets.clear();
    node_sets.reserve(this->nsets_startlines_.size());

    // The node's index. Initialized in invalid value (-1).
    int node_id = -1;

    // Create a temporary node set.
    NodeSet nset;
    
    // Iterate though mesh starting from each node set's starting line.
    for (auto &current_nset_startline: this->nsets_startlines_) {

        // Clear temporary node set.
        nset.Clear();

        // Extract node set name from its starting line.
        std::string nset_title = this->parsed_mesh_.at(current_nset_startline);
        std::string nset_name = nset_title.substr(nset_title.find_last_of("=")+1);
        nset_name.erase(std::remove_if(nset_name.begin(), nset_name.end(), isspace), nset_name.end());

        // Assign the node set name to the temporary node set.
        nset.SetName(nset_name);
        
        // Iterate node indices in the node set.
        // Skipping the node set title line.
        for (std::vector<std::string>::size_type it = current_nset_startline+1;
                                                 it != this->parsed_mesh_.size(); ++it) {

            std::string line = this->parsed_mesh_[it];

            // Replace comma occurence in line with space.
            std::replace(line.begin(), line.end(), ',', ' ');

            // Convert line to stringstream to be passed in the connectivity variables.
            std::stringstream ss(line);

            // Stop loop if header line found.
            if (line.find('*') != std::string::npos) { break; }

            // Store the indices of the nodes belonging to the nodeset.
            while (ss >> node_id) {
                // Add the mapped index of the node id to the temporary node set.
                nset.EditNodeIds().emplace_back(this->mapped_node_ids_[node_id]);
            }
            
        } // End of Iterate node indices in the node set.

        // Assing temporary node set in the node sets container.
        node_sets[nset_name] = nset;

    } // End of Iterate though mesh starting from each node set's starting line.
}


}  //end of namespace IMP



#endif //IMP_ENGINE_MESH_IO_ABAQUS_IO_TPP_
