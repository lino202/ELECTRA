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

#ifndef IMP_ENGINE_MESH_IO_VTK_IO_TPP_
#define IMP_ENGINE_MESH_IO_VTK_IO_TPP_

#include "IMP/engine/mesh_io/vtk_io.hpp"

namespace IMP {

template <short DIM, short CELL_NODES>
VtkIO<DIM, CELL_NODES>::VtkIO() : parsed_mesh_(), parts_startlines_(), nsets_startlines_(), nodes_startline_(0), elems_startline_(0),
                                  elems_type_startline_(0), parsed_cells_shape_(CellShape::unknown), parts_exist_(false), nsets_exist_(false)
{}


template <short DIM, short CELL_NODES>
VtkIO<DIM, CELL_NODES>::~VtkIO()
{}


template <short DIM, short CELL_NODES>
void VtkIO<DIM, CELL_NODES>::ReadMeshFrom(const std::string &mesh_filename)
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

    // Check if mesh is in VTK Simple Legacy format (.vtk).
    std::string ext = mesh_filename.substr(mesh_filename.length()-4);
    if (ext != ".vtk") {
        std::string error_msg = Logger::Error("Expected mesh in VTK Simple Legacy format (*.vtk). Check mesh filename: ") + mesh_filename;
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
    bool elems_type_exist = false;

    // Line of the file.
    std::string line = "";
    int line_number = 0;

    //Read mesh file line by line.
    while (std::getline(mesh_file, line)) {

        // Transform mesh file line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Read mesh file lines.
        this->parsed_mesh_.emplace_back(line);

        // Check if nodes list was found in the mesh.
        if (line.find("points") != std::string::npos) {
            this->nodes_startline_ = line_number;
            nodes_exist = true;
        }

        // Check if elements list was found in the mesh.
        if (line.find("cells") != std::string::npos) {
            this->elems_startline_ = line_number;
            elems_exist = true;
        }

        if (line.find("cell_types") != std::string::npos) {
            this->elems_type_startline_ = line_number;
            elems_type_exist = true;
        }

        // Find the first lines of the partitions sets.
        if (line.find("*elset") != std::string::npos) {
            this->parts_startlines_.emplace_back(line_number);
            if (!this->parts_exist_) { this->parts_exist_ = true; }
        }

        // Find the first lines of the node sets.
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
    if (!nodes_exist && !elems_exist && !elems_type_exist) {
        std::string error_msg = Logger::Error("Something went wrong while reading nodes and elements. Check mesh file: ") + mesh_filename;
        throw std::runtime_error(error_msg);
    }

    // Extract the element type.
    // Iterate though mesh starting from the elements type set starting line.
    // Skip the first line start from the first element's type.
    int elem_type_id = -1; int temp_id = -1;
    std::stringstream ss;
    for (std::size_t it = this->elems_type_startline_+1; it != parsed_mesh_.size(); ++it) {

        line = this->parsed_mesh_[it];

        // Convert line to stringstream to be passed in the coordinates Vec.
        ss << line;

        // Keep the type of the previous element in the list
        temp_id = elem_type_id;

        // Get the type id until reach the end of list or encounter different element type.
        if (!(ss >> elem_type_id) || (elem_type_id != temp_id && temp_id != -1)) { break; }
        
        // Reset string stream to accept next line.
        ss.str(std::string());  ss.clear();
 
    }

    // Check uniqueness of element type.
    if (elem_type_id != temp_id) {
        std::string error_str = "Could not read VTK Simple Legacy mesh file. Meshes with more than one type of elements are not currently supported.";
        throw std::runtime_error(error_str);
    }

    // Assigned the element type according to the extracted element type index.
    if (elem_type_id == 3) { this->parsed_cells_shape_ = CellShape::edge; }
    else if (elem_type_id == 5) { this->parsed_cells_shape_ = CellShape::tri; }
    else if (elem_type_id == 9) { this->parsed_cells_shape_ = CellShape::quad; }
    else if (elem_type_id == 10) { this->parsed_cells_shape_ = CellShape::tet; }
    else if (elem_type_id == 12) { this->parsed_cells_shape_ = CellShape::hex; }

}


// template <short DIM, short CELL_NODES>
// void VtkIO<DIM, CELL_NODES>::SaveMesh(const Mesh<DIM, CELL_NODES> &mesh, std::string &mesh_filename)
// {
//     // Check if mesh filename is given.
//     if (mesh_filename.empty()) {
//         std::string error_msg = Logger::Error("Could not save the mesh. Check output mesh filename: ") + mesh_filename;
//         throw std::invalid_argument(error_msg);
//     }

//     // Check if given output filename is in abaqus format (.inp).
//     if (mesh_filename.substr(mesh_filename.length() - 4) != ".inp") {
//         // Add abaqus (.inp) extension.
//         mesh_filename += ".inp";
//     }

//     //Open mesh output file.
//     std::ofstream out(mesh_filename, std::ios::out);

//     // Check if mesh file opened successfully.
//     if (!out.is_open()) {
//         std::string error_msg = Logger::Error("Could not save the mesh to the output mesh file. Check file path: ") + mesh_filename;
//         throw std::runtime_error(error_msg);
//     }

//     // Save header information.
//     out << "*Heading\n";
//     out << "** Job name: imp-mesh-job Model name: Model-1\n";
//     out << "** Generated by: IMP library, echo=NO, model=NO, history=NO, contact=NO\n";
//     out << "**\n";
//     out << "** PARTS\n";
//     out << "**\n";
//     out << "*Part, name=Part-1\n";
//     out << "*Node\n";

//     // Save mesh nodes coordinates.
//     for (auto &node : mesh.Nodes()) {

//         auto node_id = &node - &mesh.Nodes()[0];

//         // Output node id and coordinates.
//         out << node_id + 1 << ",";
//         for (short i = 0; i != DIM-1; ++i) {
//             out << node.Coordinates()[i] << ",";
//         }
//         out << node.Coordinates()[DIM-1] << std::endl;

//     }

//     // Save mesh elements type info.
//     if (DIM == 2 && CELL_NODES == 3) { out << "*Element, type=CPE3\n"; }
//     else if (DIM == 2 && CELL_NODES == 4) { out << "*Element, type=CPE4\n"; }
//     else if (DIM == 3 && CELL_NODES == 3) { out << "*Element, type=SFM3D3\n"; }
//     else if (DIM == 3 && CELL_NODES == 4) { out << "*Element, type=C3D4\n"; }
//     else if (DIM == 3 && CELL_NODES == 8) { out << "*Element, type=C3D8\n"; }
//     else { 
//         std::string error_msg = Logger::Error("Could not write mesh of not supported element type. "
//         "Supported: 3-node Triangle, 4-node Quadrilateral, 4-node Tetrahedron, and 8-node Hexahedron.");
        
//         throw std::runtime_error(error_msg); 
//     }

//     // Save mesh elements.
//     for (auto &cell : mesh.Cells()) {

//         auto cell_id = &cell - &mesh.Cells()[0];

//         // Output element id and connectivity.
//         out << cell_id+1 << ",";
//         for (short i = 0; i != CELL_NODES-1; ++i) {
//             out << cell.Connectivity()[i] << ",";
//         }
//         out << cell.Connectivity()[CELL_NODES-1] << std::endl;

//     }

//     // Save the ending mesh header.
//     out << "*End Part\n" << "**\n";
//     out << "** ASSEMBLY\n" << "**\n" << "*Assembly, name=Assembly\n";
//     out << "**\n" << "*Instance, name=Part-1-1, part=Part-1\n";
//     out << "*End Instance\n" << "**\n" << "*End Assembly";

//     // Close the output file.
//     out.close();

// }


template <short DIM, short CELL_NODES>
void VtkIO<DIM, CELL_NODES>::LoadNodesIn(std::vector<Vec<DIM, double>> &nodes)
{
    // Get the number of nodes to load.
    std::string line = this->parsed_mesh_[this->nodes_startline_];
    int total_nodes = std::stoi(line.substr(line.find_first_of(" "), line.find_last_of("0123456789")));
    
    // Clean the mesh nodes container.
    nodes.clear();
    nodes.resize(total_nodes, Vec<DIM, double>());

    // The key id of a node in the parsed mesh.
    std::size_t node_id = 0;

    // The spatial coordinate value of a node in the mesh.
    double coord_val = 0.;

    // The index of the spatial coordinate.
    short coord_id = 0; 

    // Iterate though mesh starting from the nodes set starting line.
    // Skip the first line start from the first node's coordinates.
    std::stringstream ss;
    for (int it = this->nodes_startline_+1; it != this->elems_startline_; ++it) {

        // Read the line of the nodes coordinates list.
        line = this->parsed_mesh_[it];

        // Convert line to stringstream.
        ss << line;

        // Extract coordinates until the end of the line.
        while (ss >> coord_val) {  

            // std::cout << coord_val << " ";

            // Store the corresponding coordinate to the corresponding node.
            if (coord_id < DIM) { nodes[node_id][coord_id] = coord_val; }

            // Increase the coordinate index.
            coord_id++;

            // Move to the next node index if three coordinates have been read.
            // and reset the coordinate index. VTK files store three coordinates per node.
            if (coord_id == 3) { node_id++; coord_id = 0; }
        }

        // Clean the stringstream for the next step.
        ss.str(std::string()); 
        ss.clear();
 
    }

}


template <short DIM, short CELL_NODES>
void VtkIO<DIM, CELL_NODES>::LoadElementsIn(std::vector<Cell<DIM, CELL_NODES>> &cells)
{

    // Get the number of elements to load.
    std::string line = this->parsed_mesh_[this->elems_startline_];
    line = line.substr(line.find_first_of("0123456789"));
    int total_elems = std::stoi(line.substr(0, line.find_first_of(" ")));

    // Clean the cells container.
    cells.clear();
    cells.resize(total_elems, Cell<DIM, CELL_NODES>());

    // Initialize an empty element and connectivity.
    Cell<DIM, CELL_NODES> cell;
    Vec<CELL_NODES, int> conn;

    // Set the elements shape.
    cell.SetShape(this->parsed_cells_shape_);

    // The element's key index in the parsed mesh.
    std::size_t el_id = 0;

    // Dummy key.
    std::size_t dummy = 0;

    // Iterate though the mesh starting from the first line of the elements list.
    // Skip a line to read the first element's connectivity.
    std::stringstream ss;
    for (int it = this->elems_startline_+1; it != this->elems_startline_+total_elems+1; ++it) {

        line = this->parsed_mesh_[it];

        // Convert line to stringstream to be passed in the connectivity variables.
        ss << line;

        // Get the key id and the connectivity of the current element.
        if (!(ss >> dummy >> conn)) { break; }

        // Set updated connectivity to the cell and store it.
        cell.SetConnectivity(conn);
        cells[el_id] = cell;
        el_id++;

        // Clean the stringstream for the next step.
        ss.str(std::string()); 
        ss.clear();

    }
    
}


template <short DIM, short CELL_NODES>
void VtkIO<DIM, CELL_NODES>::LoadNodeSetsIn(std::unordered_map<std::string, NodeSet> &node_sets)
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
        
        // Remove last character of current boundary name if it is the new line character.
        // nset_name.erase(std::remove_if(nset_name.begin(), nset_name.end(), [](char c) { return !isalnum(c); }), nset_name.end());

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


#endif //IMP_ENGINE_MESH_IO_VTK_IO_TPP_
