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


#ifndef IMP_ENGINE_TESSELATIONS_GRID_TPP_
#define IMP_ENGINE_TESSELATIONS_GRID_TPP_


#include "IMP/engine/tesselations/grid.hpp"



namespace IMP {


template<short DIM, short CELL_VERTS>
Grid<DIM, CELL_VERTS>::Grid() : nodes_(), ghost_cells_(), node_sets_(), node_normals_(), cells_shape_(CellShape::node)
{}


template<short DIM, short CELL_VERTS>
Grid<DIM, CELL_VERTS>::~Grid()
{}


template<short DIM, short CELL_VERTS>
void Grid<DIM, CELL_VERTS>::LoadFrom(const std::string &grid_filename)
{
    // Check if grid filename is given.
    if (grid_filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load grid. The given grid filename is empty."));
    }

    // Get the extension of the grid filename.
    auto ext = grid_filename.substr(grid_filename.length()-4);

    // Clear the grid containers.
    this->nodes_.clear();

    // Load grid of given format.
    if (ext == ".inp") {
        // Load from abaqus mesh file.
        AbaqusIO<DIM, CELL_VERTS> abaqus_io;
        abaqus_io.ReadMeshFrom(grid_filename);
        abaqus_io.LoadNodesIn(this->nodes_);
        abaqus_io.LoadElementsIn(this->ghost_cells_);

        // Update the grid ghost cells type.
        this->UpdateGhostCellsShape();

        // Load node sets if they are available.
        if (abaqus_io.NodeSetsExist()) { abaqus_io.LoadNodeSetsIn(this->node_sets_); }
    }
    else if (ext == ".vtk") {
        // Load from VTK Simple Legacy format.
        VtkIO<DIM, CELL_VERTS> vtk_io;
        vtk_io.ReadMeshFrom(grid_filename);
        vtk_io.LoadNodesIn(this->nodes_);
        vtk_io.LoadElementsIn(this->ghost_cells_);

        // Update the grid ghost cells type.
        this->UpdateGhostCellsShape();
    }
    else {
        std::string error_msg = Logger::Error("Could not load grid of unknown format. Check given grid filename: ") + grid_filename;
        throw std::invalid_argument(error_msg);
    }

}


template<short DIM, short CELL_VERTS>
void Grid<DIM, CELL_VERTS>::LoadNodeNormals(const std::string &normals_filename)
{
    // Check if normal vectors filename is given.
    if (normals_filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load normal vectors. The given normal vectors filename is empty."));
    }

    // Check if the node grids are available.
    if (this->NodesNum() <= 0) {
        throw std::runtime_error(Logger::Error("Could not load normal vectors. The grid has not any nodes."));
    }

    // Initialize node normals matrix.
    this->node_normals_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->NodesNum(), DIM);

    // Open the normals file.
    std::ifstream normals_file(normals_filename, std::ios::in);

    // Check if mesh file opened successfully.
    if (!normals_file.is_open()) {
        std::string error_msg = Logger::Error("Could not open the normal vectors file. Check file path: " + normals_filename);
        throw std::runtime_error(error_msg);
    }

     //Line of the file.
    std::string line = "";
    int line_number = 0;

    // Initialize normal vectors triplet container.
    typedef Eigen::Triplet<double> T;
    std::vector<T> normal_triplets;
    normal_triplets.reserve(DIM*this->NodesNum());

    //Read normal vectors file line by line.
    int id = 0;
    double normal_d = 0.;
    std::stringstream ss;
    while (std::getline(normals_file, line)) {
        
        // Read if not a comment line.
        if (line.find("#") == std::string::npos) {
            // Replace commas with spaces if any in the processed line.
            std::replace(line.begin(), line.end(), ',', ' ');
            
            // Convert current line to stringstream.
            ss << line;

            // Get the index of the node.
            ss >> id;
            
            // Correct for c++ storing.
            id--;

            for (short d = 0; d != DIM; ++d) {
                // Get the normal value at the coordinate of the d axis.
                ss >> normal_d;

                // Store the normal value in the normals container.
                normal_triplets.emplace_back(T(id, d, normal_d));
            }

            // Reset stringstream.
            ss.str(""); ss.clear();
        }
    }
    normal_triplets.shrink_to_fit();

    // Set the node normals matrix from the triplets.
    this->node_normals_.setFromTriplets(std::begin(normal_triplets), std::end(normal_triplets));
    
    // Close the normals file.
    normals_file.close();

}


template<short DIM, short CELL_VERTS>
void Grid<DIM, CELL_VERTS>::UpdateGhostCellsShape()
{
    // Update the shape to the shape of the first cell.
    this->cells_shape_ = this->ghost_cells_[0].Shape();

    // Check if cells of different shape exist in the mesh.
    for (const auto &cell : this->ghost_cells_) {

        // Update to mixed shape if found a cell with different shape.
        if (cell.Shape() != this->cells_shape_) { 
            this->cells_shape_ = CellShape::mixed;

            return;
        }
    }

}


} // End of namespace IMP.

#endif //IMP_ENGINE_TESSELATIONS_GRID_TPP_