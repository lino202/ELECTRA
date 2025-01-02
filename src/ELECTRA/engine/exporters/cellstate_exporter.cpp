/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019
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


#include "ELECTRA/engine/exporters/cellstate_exporter.hpp"


namespace ELECTRA
{

CellStateExporter::CellStateExporter()
{}


CellStateExporter::~CellStateExporter()
{}


void CellStateExporter::WriteCellsState(const std::vector<std::unique_ptr<EpBasic>> &cells, const std::string &filename)
{
    namespace boost_fs = boost::filesystem;

    // Get the path directory of the filename.
    boost_fs::path p(filename);
    std::string path = p.parent_path().string();

    // Create the path's directory if it doesn't exist.
    boost_fs::path dir(path);
    if (!path.empty() && !boost_fs::exists(dir)) { boost_fs::create_directories(dir); }

    // Search for the extension of the file.
    std::string ext = p.extension().string();

    // Initialize the output with proper extension.
    std::ofstream output;
    if (ext == ".elc") {
        output.open(filename, std::ios::out | std::ios::binary);
    } else {
        output.open(filename+".elc", std::ios::out | std::ios::binary);
    }

    // Check if file opened properly.
    if(!output.good()) {
        std::string error_msg = "CellStateExporter failed to write cells state in file. Check given file: " + filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the cells state header file.
    // ONLY SAVE STATES 
    // We only save vars and currs as prms and curr_blocks should be defined from the manual init file. As nodes/cells prms and block are usually common for several of them.
    // In this way we save storage space for huge meshes
    // Here cells = nodes 
    std::size_t header_size = this->header_.size();

    output.write(reinterpret_cast<char *>(&header_size), sizeof(std::size_t));
	output.write(this->header_.c_str(), header_size);
    
    // Write number of cells.
    std::size_t cells_num = cells.size();
    output.write(reinterpret_cast<char *>(&cells_num), sizeof(std::size_t));

    // Write the state of each of the cells' variables.
    int32_t cell_type_id = 0;
    double val = 0.;
    for (const auto &cell : cells) {
        cell_type_id = static_cast<int32_t>(cell->ModelType());
        output.write(reinterpret_cast<char *>(&cell_type_id), sizeof(int32_t));

        // Write variables.
        for (int32_t i = 0; i != cell->VarNum(); ++i) {
            val = cell->Var(i);
            output.write(reinterpret_cast<char *>(&val), sizeof(double));
        }

        // Write currents.
        for (int32_t i = 0; i != cell->CurrentNum(); ++i) {
            val = cell->Current(i);
            output.write(reinterpret_cast<char *>(&val), sizeof(double));
        }

    }
    
    // Close the output file.
    output.close();
}




void CellStateExporter::ReadCellsState(std::vector<std::unique_ptr<EpBasic>> &cells, const std::string &filename)
{
    namespace boost_fs = boost::filesystem;

    // Get the path directory of the filename.
    boost_fs::path p(filename);
    std::string path = p.parent_path().string();

    // Create the path's directory if it doesn't exist.
    boost_fs::path dir(path);
    if (!path.empty() && !boost_fs::exists(dir)) { boost_fs::create_directories(dir); }

    // Search for the extension of the file.
    std::string ext = p.extension().string();

    // Initialize the input with proper extension.
    std::ifstream input;
    if (ext == ".elc") {
        input.open(filename, std::ios::in | std::ios::binary);
    } else {
        input.open(filename+".elc", std::ios::in | std::ios::binary);
    }

    // Check if file opened properly.
    if(!input.good()) {
        std::string error_msg = "CellStateExporter failed to write cells state in file. Check given file: " + filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Read header
    std::size_t read_header_size;
    input.read(reinterpret_cast<char *>(&read_header_size), sizeof(std::size_t));
	std::string read_header(read_header_size, '0');
    input.read(&read_header[0], read_header_size);

    if(read_header != this->header_){
        std::string error_msg = "CellStateExporter failed to load cells state in file. Current ELECTRA version has a different header from the one read";
        throw std::runtime_error(Logger::Error(error_msg));
    }
    
    // Read number of cells.
    std::size_t read_cells_num;
    input.read(reinterpret_cast<char *>(&read_cells_num), sizeof(std::size_t));
    
    if(read_cells_num != cells.size()){
        std::string error_msg = "CellStateExporter failed to load cells state in file. Cells number differs";
        throw std::runtime_error(Logger::Error(error_msg));
    }


    // Read the state of each of the cells' variables.
    // We check if the cell type is the same and overwrite the states (vars and currs)
    int32_t cell_type_id, read_cell_type_id;
    double val;
    int32_t j=0;
    for (const auto &cell : cells) {
        cell_type_id = static_cast<int32_t>(cell->ModelType());
        input.read(reinterpret_cast<char *>(&read_cell_type_id), sizeof(int32_t));

        if(cell_type_id != read_cell_type_id){
            std::string error_msg = "CellStateExporter failed to load cells state in file. Different cell_type_id at cell/node: " + std::to_string(j);
            throw std::runtime_error(Logger::Error(error_msg));
        }

        // Write variables.
        for (int32_t i = 0; i != cell->VarNum(); ++i) {
            input.read(reinterpret_cast<char *>(&val), sizeof(double));
            cell->SetVar(i, val);
        }

        // Write currents.
        for (int32_t i = 0; i != cell->CurrentNum(); ++i) {
            input.read(reinterpret_cast<char *>(&val), sizeof(double));
            cell->SetCurrent(i, val);
        }

        j++;

    }
    
    // Close the output file.
    input.close();
}

} // End of namespace ELECTRA
