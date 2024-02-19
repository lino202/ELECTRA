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


#include "IMP/engine/mesh_io/abaqus_io.hpp"

namespace IMP {


AbaqusIO::AbaqusIO() : parsed_mesh_(), mesh_type_(), parts_exist(false), nsets_exist(false),
                       nodes_startline_(0), elems_startline_(0), parts_startlines_(),
                       nsets_startlines_(), mapped_node_ids_(), mapped_elem_ids_()
{}


AbaqusIO::~AbaqusIO()
{}


void AbaqusIO::ReadMeshFrom(const std::string &mesh_filename)
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
        if (line.find("cpe3") != std::string::npos) { this->mesh_type_ = MeshType::Triangular; }
        else if (line.find("cpe4") != std::string::npos) { this->mesh_type_ = MeshType::Quadrilateral; }
        else if (line.find("s3r") != std::string::npos) { this->mesh_type_ = MeshType::Triangular; }
        else if (line.find("s4r") != std::string::npos) { this->mesh_type_ = MeshType::Quadrilateral; }
        else if (line.find("c3d4") != std::string::npos) { this->mesh_type_ = MeshType::Tetrahedral; }
        else if (line.find("c3d8") != std::string::npos) { this->mesh_type_ = MeshType::Hexahedral; }
        else { this->mesh_type_ = MeshType::Unknown; }

        // Find the first lines of the partitions sets.
        if (line.find("*elset") != std::string::npos) {
            this->parts_startlines_.emplace_back(line_number);
            if (!this->parts_exist) { this->parts_exist_ = true; }
        }

        // Find the first lines of the boundary node sets.
        if (line.find("*nset") != std::string::npos) {
            this->nsets_startlines_.emplace_back(line_number);
            if (!this->nsets_exist) { this->nsets_exist_ = true; }
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


// void AbaqusIO::LoadPartitionsIn(std::vector<Tetrahedron> &tetras)
// {
//     // The tetrahedron's index. Initialized in invalid value (-1).
//     int tetra_id = -1;

//     // The headerline of the current part. Used in partitions iteration.
//     std::string current_part_headerline = "";

//     // The name of the current part. Used in partitions iteration.
//     std::string current_part_name = "";

//     // The mesh file line to be processed. Used in partitions iteration.
//     std::string line = "";

//     // Iterate though mesh starting from the partition set starting line (for each partition).
//     // Skip the first line to start from the partition's elements.
//     for (auto &current_part_startline: this->parts_startlines_) {

//         auto part_id = &current_part_startline - &this->parts_startlines_[0];

//         current_part_headerline = this->input_mesh_.at(current_part_startline);

//         current_part_name = current_part_headerline.substr(current_part_headerline.find_last_of("=")+1);

//         // Remove last character of current part name if it is the new line character.
// //        if (current_part_name.back() == '\n') { current_part_name.pop_back(); }

//         current_part_name.erase(std::remove_if(current_part_name.begin(), current_part_name.end(),
//                                                [](char c) { return !isalnum(c); }), current_part_name.end());

//         for (std::vector<std::string>::size_type it = current_part_startline+1;
//              it != input_mesh_.size(); ++it) {

//             line = this->input_mesh_.at(it);

//             // Replace comma occurence in line with space.
//             std::replace(line.begin(), line.end(), ',', ' ');

//             // Convert line to stringstream to be passed in the connectivity variables.
//             std::stringstream ss(line);

//             while (ss >> tetra_id) {
//                 // Correct for c++ storage by paddling 1.
//                 tetra_id -= 1;

//                 // Correct for offsetted elements if necessary.
//                 if (this->offsetted_elems_.size() != 0) {

//                     // Check for offset at element in partition.
//                     auto offsetted_elem = std::find_if(this->offsetted_elems_.begin(), this->offsetted_elems_.end(),
//                           [&](const std::pair<int,int> &element){ return element.first == tetra_id; } );

//                     // Apply offset correction if necessary.
//                     if (offsetted_elem != this->offsetted_elems_.end()) {
//                         tetra_id -= offsetted_elem->second;
//                     }
//                 }

//                 // Set the partition.
//                 tetras.at(tetra_id).SetPartition(part_id, current_part_name);
//             }

//         }

//     }

// }





} // end of namespace IMP
