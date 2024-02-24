/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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


#ifndef ELECTRA_EXPORTERS_ENSIGHT_EXPORTER_TPP_
#define ELECTRA_EXPORTERS_ENSIGHT_EXPORTER_TPP_


#include "ELECTRA/engine/exporters/ensight_exporter.hpp"


namespace ELECTRA {

template<short DIM, short CELL_NODES>
EnsightExporter<DIM, CELL_NODES>::EnsightExporter()
{}

template<short DIM, short CELL_NODES>
EnsightExporter<DIM, CELL_NODES>::~EnsightExporter()
{}


template<short DIM, short CELL_NODES>
void EnsightExporter<DIM, CELL_NODES>::SetFiles(const std::string &geom_file, const std::string &states_file, const std::string &anim_file)
{

    //  Set the files we are going to use for saving results 
    // Get Tissue geometry file
    this->geom_file_    = geom_file;
    std::string geo_ext = std::filesystem::path(this->geom_file_).extension();
    if (geo_ext != ".geo") {this->geom_file_ += ".geo";}


    // Get Tissue states file prefix
    std::filesystem::path states_path(states_file);
    if (states_path.has_extension()) {
        std::filesystem::path e = "";
        states_path.replace_extension(e);
        this->states_file_ = states_path.string();
    }

    // Get Tissue animation file
    this->anim_file_     = anim_file;
    std::string anim_ext = std::filesystem::path(this->anim_file_).extension();
    if (anim_ext != ".case") { this->anim_file_ += ".case"; }

    // We check all results path are the same and are absolute
    namespace boost_fs = boost::filesystem;

    boost_fs::path path_geom(this->geom_file_);
    std::string parent_path_geom = path_geom.parent_path().string();

    boost_fs::path path_states(this->states_file_);
    std::string parent_path_states = path_states.parent_path().string();

    boost_fs::path path_anim(this->anim_file_);
    std::string parent_path_anim = path_anim.parent_path().string();

    if (! ( (parent_path_anim==parent_path_geom) && (parent_path_anim==parent_path_states) ) ){
        throw std::invalid_argument(Logger::Error("Animation, states and geometry need to be save in the SAME ABSOLUTE paths"));
    }

    // Create the path's directory if it doesn't exist.
    boost_fs::path dir(parent_path_anim);
    if (!boost_fs::exists(dir)) { boost_fs::create_directories(dir); }

}


template<short DIM, short CELL_NODES>
void EnsightExporter<DIM, CELL_NODES>::SaveGeo(const std::vector<IMP::Vec<DIM, double>> &nodes, const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells)
{

    // Open output geometry file.
    std::ofstream geo_out;
    geo_out.open(this->geom_file_, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

    // Check if file was properly opened.
    if (!geo_out.is_open()) {
        throw std::invalid_argument(Logger::Error("Could not create Ensight geometry file. Check file path: " + this->geom_file_));
    }

    // Write File in Binary 
    // Write header.
    WriteString_("C Binary", geo_out);
    WriteString_("Ensight binary geometry file", geo_out);
    WriteString_((std::string("Written by ELECTRA v") + ELECTRA_VERSION).c_str(), geo_out);
    WriteString_("node id off", geo_out);
    WriteString_("element id off", geo_out);
    WriteString_("part", geo_out);
    WriteInt_(1, geo_out);
    WriteString_("Main part", geo_out);  // This is a description/name of the part but we have only one for now
    WriteString_("coordinates", geo_out);

    // Write nodes total number and coordinates
    const std::size_t npts = nodes.size();
    WriteInt_(npts, geo_out);

    for (short d = 0; d != DIM; ++d) {
        for (std::size_t id = 0; id != npts; ++id){
            WriteFloat_(nodes[id][d], geo_out);    
        }
    }

    // Fill missing coordinate dimensions wit zeros when less than 3 dimensions.
    if (DIM < 3) {
        for (std::size_t id = 0; id != (3-DIM)*npts; ++id){
            WriteFloat_(0.0, geo_out);  
        }
    }

    // Cells 
    // Create a flag for each connected node of the grid.
    // NOTE: This seems unuseful, when there will be unconnected nodes?
    Eigen::SparseVector<int> conn_node_ids(npts);

    // Write cells connectivity.
    if (DIM == 1 && CELL_NODES == 2)      { WriteString_("bar2",   geo_out); }
    else if (DIM == 2 && CELL_NODES == 2) { WriteString_("bar2",   geo_out); }
    else if (DIM == 2 && CELL_NODES == 3) { WriteString_("tria3",  geo_out); }
    else if (DIM == 2 && CELL_NODES == 4) { WriteString_("quad4",  geo_out); }
    else if (DIM == 3 && CELL_NODES == 2) { WriteString_("bar2",   geo_out); }
    else if (DIM == 3 && CELL_NODES == 3) { WriteString_("tria3",  geo_out); }
    else if (DIM == 3 && CELL_NODES == 4) { WriteString_("tetra4", geo_out); }
    else if (DIM == 3 && CELL_NODES == 8) { WriteString_("hexa8",  geo_out); }
    else {
        std::string error_msg = "Could not create Ensight geometry file due to not supported elements with Dimensions: " + 
                                std::to_string(DIM) + " and Nodes: " + std::to_string(CELL_NODES);
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Cells number and node idxs that form them
    const std::size_t ncells = cells.size();
    WriteInt_(ncells, geo_out);

    for (const auto &cell : cells) {
        for (short nid = 0; nid != CELL_NODES-1; ++nid) {
            WriteInt_(cell.N(nid)+1, geo_out);
            conn_node_ids.coeffRef(cell.N(nid)) = 1;
        }
        WriteInt_(cell.N(CELL_NODES-1)+1, geo_out);
        conn_node_ids.coeffRef(cell.N(CELL_NODES-1)) = 1;
    }

    // Output 0D cells of not connected nodes. 
    // This seems unuseful, when there will be unconnected nodes?
    if (conn_node_ids.nonZeros() < static_cast<int>(npts)) {
        WriteString_("point", geo_out);
        WriteInt_(static_cast<int>(npts) - conn_node_ids.nonZeros(), geo_out);
        for (int i = 0; i != conn_node_ids.rows(); ++i) {
            if (conn_node_ids.coeff(i) == 0){
                WriteInt_(i+1, geo_out);
            }
        }
    }

    // Close the geometry output file.
    geo_out.close();

    std::cout << ELECTRA::Logger::Message("Saved Ensight geometry: " + this->geom_file_ + "\n");
}


template<short DIM, short CELL_NODES>
void EnsightExporter<DIM, CELL_NODES>::SaveStates(const Eigen::VectorXd &scalar_field, const std::string &state_number)
{

    // Open output scalar file.
    std::ofstream scalar_out;
    const std::string out_filename = this->states_file_ + state_number + ".ens";
    scalar_out.open(out_filename, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

    // Check if file was properly opened.
    if (!scalar_out.is_open()) {
        throw std::invalid_argument(Logger::Error("Could not create Ensight scalar field file. Check file path: " + out_filename));
    }
    
    // Write header file.
    WriteString_("Ensight Model Post Process", scalar_out);
    WriteString_("part", scalar_out);
    WriteInt_(1, scalar_out);
    WriteString_("coordinates", scalar_out);

    // Write scalar field.
    for (int i = 0; i != scalar_field.size(); ++i) {
        WriteFloat_(scalar_field.coeff(i), scalar_out);
    }

    // Close file.
    scalar_out.close();

}


template<short DIM, short CELL_NODES>
void EnsightExporter<DIM, CELL_NODES>::SaveAnimation(int steps_num, double time_inc)
{
    // Open output animation file. This is an ascii file
    std::ofstream case_out;
    case_out.open(this->anim_file_, std::ofstream::out | std::ofstream::trunc);


    // Isolate the name of the scalar and geom file names.
    namespace boost_fs = boost::filesystem;
    boost_fs::path scalar_p(this->states_file_);
    std::string scalar_name = scalar_p.filename().string();
    scalar_name = scalar_name.substr(0, scalar_name.find_first_of("0123456789"));

    boost_fs::path geo_p(this->geom_file_);
    std::string geom_file_name = geo_p.filename().string();

    // Set number of wildcard symbols for steps.
    //  NOTE: that steps_num comes from react_diff->states_num that have the number of states saved but the .ens names start from 0 so we 
    // have to rest 1 to states num so we have the maximum number in the naming of .ens files for generating the '**' in the animation file
    // This is specially important when we save for example 100 states, then the maximum .ens is ...99.ens and we need only two asterisks
    // not three!!
    std::size_t wildcard_num = std::to_string(steps_num-1).size();
    std::string asterisks(wildcard_num, '*');

    // Write header.
    case_out << "# Ensight output generated by: ELECTRA v." << ELECTRA_VERSION << "\n\n";
    case_out << "FORMAT\ntype: ensight gold\n\n";
    case_out << "GEOMETRY\n";
    case_out << "model: " << geom_file_name << "\n\n";
    
    case_out << "VARIABLE\n";    
    case_out << "scalar per node: Potential " << scalar_name + asterisks + ".ens" << "\n\n";
    
    case_out << "TIME\n";
    case_out << "time set: " << time_inc << "\n";
    case_out << "number of steps: " << steps_num << "\n";
    case_out << "filename start number: 0\n";
    case_out << "filename increment: 1\n";
    case_out << "time values: ";

    // Write time increments
    int skip = 0;
    for (int i = 0; i != steps_num; ++i) {
        case_out << i*time_inc << " ";
        skip++;

        if (skip == 30) {
            case_out << "\n             ";
            skip = 0;
        }
    }

    // Close file
    case_out.close();

    std::cout << ELECTRA::Logger::Message("Saved Ensight animation: " + this->anim_file_ + "\n");
}


} //end of namespace ELECTRA


#endif // ELECTRA_EXPORTERS_ENSIGHT_EXPORTER_TPP_