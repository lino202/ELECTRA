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

#ifndef ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_TPP_

#include "config_output.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigOutput<DIM, CELL_NODES>::ConfigOutput()
{}


template<short DIM, short CELL_NODES>
ConfigOutput<DIM, CELL_NODES>::~ConfigOutput()
{}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputGeneration(const Parser &parser, const std::vector<IMP::Vec<DIM,double>> &nodes,
        const std::vector<IMP::Cell<DIM,CELL_NODES>> &cells, const std::shared_ptr<ReactionDiffusion<DIM, CELL_NODES>> &react_diff, std::ostream &stream) const
{
    // Output if paraview format was requested.
    if (parser.HasAttribute("output.paraview")) {
        throw std::invalid_argument(Logger::Error("Output as paraview has been deprecated, try ensight"));
    }

    // Output if ensight format was requested.
    if (parser.HasAttribute("output.ensight")) {
        this->OutputToEnsight(nodes, cells, react_diff);
    }

    // Output if ascii format was requested.
    if (parser.HasAttribute("output.ascii")) {
        throw std::invalid_argument(Logger::Error("Output as ascii has been deprecated, try ensight"));
    }


    // Output if cell state was requested.
    if (parser.HasAttribute("output.cells state")) {
        this->OutputToCellStates(parser, react_diff, stream);
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputFibers(const Parser &parser, IMP::Mesh<DIM,CELL_NODES> mesh, IMP::Voronoi<DIM> voro, const Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const
{
    // Output partitioned mesh
    if (parser.HasAttribute("output.fibers.partitioned geometry")) {
        std::string method = parser.GetValue<std::string>("numerical approximation.method");
        std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

        std::string filename = parser.GetValue<std::string>("output.fibers.partitioned geometry");

        // Add ventricular node sets in the mesh.
        if (ldrbm.SeptumNodeIds().size() != 0) {
            // Left ventricle nodes.
            IMP::NodeSet lv_nset; lv_nset.Set("left_ventricle", ldrbm.LeftVentriNodeIds());

            // Right ventricle nodes.
            IMP::NodeSet rv_nset; rv_nset.Set("right_ventricle", ldrbm.RightVentriNodeIds());

            // Interface ventricle nodes.
            IMP::NodeSet intv_nset; intv_nset.Set("ventricle_interface", ldrbm.VentriInterfaceNodeIds());

            // Septum ventricle nodes.
            IMP::NodeSet septum_nset; septum_nset.Set("septum", ldrbm.VentriInterfaceNodeIds());

            if (method == "fem") {
                mesh.AddNodeSet(lv_nset);
                mesh.AddNodeSet(rv_nset);
                mesh.AddNodeSet(intv_nset);
                mesh.AddNodeSet(septum_nset);
            } else {
                voro.AddNodeSet(lv_nset);
                voro.AddNodeSet(rv_nset);
                voro.AddNodeSet(intv_nset);
                voro.AddNodeSet(septum_nset);
            }

        } else {
            // Add atrial node sets in the mesh.
            // Left atrium nodes.
            IMP::NodeSet la_nset; la_nset.Set("left_atrium", ldrbm.LeftAtriNodeIds());

            // Right atrium nodes.
            IMP::NodeSet ra_nset; ra_nset.Set("right_atrium", ldrbm.RightAtriNodeIds());

            if (method == "fem") {
                mesh.AddNodeSet(la_nset);
                mesh.AddNodeSet(ra_nset);
            } else {
                voro.AddNodeSet(la_nset);
                voro.AddNodeSet(ra_nset);
            }
        }

        // Save the mesh.
        // if (method == "fem") {
            mesh.SaveTo(filename);
            stream << Logger::Message("Saved partitioned mesh: " + filename + "\n");
        // } else {
        //     voro.SaveTo(filename);
        //     stream << Logger::Message("Saved partitioned voronoi tesselation: " + filename + "\n");
        // }
    }

    // Output longitudinal fibers.
    if (parser.HasAttribute("output.fibers.longitudinal fibers")) {
        std::string fibs_file = parser.GetValue<std::string>("output.fibers.longitudinal fibers");
        this->SaveFibers(ldrbm.LongFiberDirection(), fibs_file);
        stream << Logger::Message("Saved longitudinal fibers: " + fibs_file + "\n");
    }

    // Output sheet fibers.
    if (parser.HasAttribute("output.fibers.sheet fibers")) {
        std::string fibs_file = parser.GetValue<std::string>("output.fibers.sheet fibers");
        this->SaveFibers(ldrbm.SheetFiberDirection(), fibs_file);
        stream << Logger::Message("Saved sheet fibers: " + fibs_file + "\n");
    }

    // Output transversal fibers.
    if (parser.HasAttribute("output.fibers.transversal fibers")) {
        std::string fibs_file = parser.GetValue<std::string>("output.fibers.transversal fibers");

        // Compute transversal fibers as the outer product of longitudinal and sheet fibers.
        Eigen::MatrixXd trans_fibers = Eigen::MatrixXd::Zero(ldrbm.SheetFiberDirection().rows(), ldrbm.SheetFiberDirection().cols());
        Eigen::Vector3d v_l, v_s;
        for (Eigen::Index i = 0; i != ldrbm.SheetFiberDirection().rows(); ++i) {
            v_l = ldrbm.LongFiberDirection().row(i);
            v_s = ldrbm.SheetFiberDirection().row(i);
            trans_fibers.row(i) = v_l.cross(v_s);
        }

        this->SaveFibers(trans_fibers, fibs_file);
        stream << Logger::Message("Saved transversal fibers: " + fibs_file + "\n");
    }

    // Output transmural distance.
    if (parser.HasAttribute("output.fibers.transmural distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.transmural distance");
        this->SaveDistanceField(ldrbm.TransmuralDistance(), filename);
        stream << Logger::Message("Saved transmural distance: " + filename + "\n");
    }

    // Output appendage-veins distance.
    if (parser.HasAttribute("output.fibers.appendage veins distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.appendage veins distance");
        this->SaveDistanceField(ldrbm.AppendageVeinsDistance(), filename);
        stream << Logger::Message("Saved appendage veins distance: " + filename + "\n");
    }

    // Output interveins distance.
    if (parser.HasAttribute("output.fibers.interveins distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.interveins distance");
        this->SaveDistanceField(ldrbm.InterVeinsDistance(), filename);
        stream << Logger::Message("Saved interveins distance: " + filename + "\n");
    }

    // Output valve-veins distance.
    if (parser.HasAttribute("output.fibers.valve veins distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.valve veins distance");
        this->SaveDistanceField(ldrbm.ValveVeinsDistance(), filename);
        stream << Logger::Message("Saved valve veins distance: " + filename + "\n");
    }

    // Output atrium-tricuspid distance.
    if (parser.HasAttribute("output.fibers.atrium tricuspid distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.atrium tricuspid distance");
        this->SaveDistanceField(ldrbm.AtriTricuspidDistance(), filename);
        stream << Logger::Message("Saved atrium tricuspid distance: " + filename + "\n");
    }

    // Output apicobasal distance.
    if (parser.HasAttribute("output.fibers.apicobasal distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.apicobasal distance");
        this->SaveDistanceField(ldrbm.ApicobasalDistance(), filename);
        stream << Logger::Message("Saved apicobasal distance: " + filename + "\n");
    }

    // Output septal distance.
    if (parser.HasAttribute("output.fibers.septal distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.septal distance");
        this->SaveDistanceField(ldrbm.SeptalDistance(), filename);
        stream << Logger::Message("Saved septal distance: " + filename + "\n");
    }

    // Output intraventricular function.
    if (parser.HasAttribute("output.fibers.intraventricular function")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.intraventricular function");
        this->SaveDistanceField(ldrbm.IntraventricularFunction(), filename);
        stream << Logger::Message("Saved intraventricular function: " + filename + "\n");
    }

    // Output transmural direction.
    if (parser.HasAttribute("output.fibers.transmural direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.transmural direction");
        this->SaveDirectionField(ldrbm.TransmuralDirection(), filename);
        stream << Logger::Message("Saved transmural direction: " + filename + "\n");
    }

    // Output appendage-veins direction.
    if (parser.HasAttribute("output.fibers.appendage veins direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.appendage veins direction");
        this->SaveDirectionField(ldrbm.AppendageVeinsDirection(), filename);
        stream << Logger::Message("Saved appendage veins direction: " + filename + "\n");
    }

    // Output interveins direction.
    if (parser.HasAttribute("output.fibers.interveins direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.interveins direction");
        this->SaveDirectionField(ldrbm.InterVeinsDirection(), filename);
        stream << Logger::Message("Saved interveins direction: " + filename + "\n");
    }

    // Output valve-veins direction.
    if (parser.HasAttribute("output.fibers.valve veins direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.valve veins direction");
        this->SaveDirectionField(ldrbm.ValveVeinsDirection(), filename);
        stream << Logger::Message("Saved valve veins direction: " + filename + "\n");
    }

    // Output atrium-tricuspid direction.
    if (parser.HasAttribute("output.fibers.atrium tricuspid direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.atrium tricuspid direction");
        this->SaveDirectionField(ldrbm.AtriTricuspidDirection(), filename);
        stream << Logger::Message("Saved atrium tricuspid direction: " + filename + "\n");
    }

    // Output apicobasal direction.
    if (parser.HasAttribute("output.fibers.apicobasal direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.apicobasal direction");
        this->SaveDirectionField(ldrbm.ApicobasalDirection(), filename);
        stream << Logger::Message("Saved apicobasal direction: " + filename + "\n");
    }

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputToEnsight(const std::vector<IMP::Vec<int(DIM), double>> &nodes,
        const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM, CELL_NODES>> &react_diff) const
{

    if (react_diff->ens_exporter_tissue_ptr) {
        // Save GEOMETRY to file
        react_diff->ens_exporter_tissue_ptr->SaveGeo(nodes, cells);

        // Save ANIMATION to file 
        react_diff->ens_exporter_tissue_ptr->SaveAnimation(react_diff->states_num, react_diff->OutputSteps()*react_diff->Dt());
    }

    // All the same as above but for the CS-------------------------------------------------------
    if (react_diff->ens_exporter_cs_ptr) {
        
        // Save GEOMETRY to file.
        react_diff->ens_exporter_cs_ptr->SaveGeo(react_diff->ConductSystem().Nodes(), react_diff->ConductSystem().Segments());


        // Save ANIMATION to file 
        react_diff->ens_exporter_cs_ptr->SaveAnimation(react_diff->states_num, react_diff->OutputSteps()*react_diff->Dt());
        
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputToCellStates(const Parser &parser, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM, CELL_NODES>> &react_diff, std::ostream &stream) const
{
    ELECTRA::CellStateExporter cellstate_export;
    if (parser.HasAttribute("output.cells state")) {
        std::string cell_state_file = parser.GetValue<std::string>("output.cells state");
        std::string cell_state_ext = std::filesystem::path(cell_state_file).extension();
        if (cell_state_ext != ".elc") { cell_state_file += ".elc"; }

        cellstate_export.WriteCellsState(react_diff->Cells(), cell_state_file);
        stream << ELECTRA::Logger::Message("Saved cells state: " + cell_state_file + "\n");
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::SaveFibers(const Eigen::MatrixXd &fibers, const std::string &output_filename) const
{
    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(output_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }

    // Initialize the output with proper extension.
    std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
    if(!output.good()) {
        std::string error_msg = "ConfigOutput failed to write fibers in file. Check given file: " + output_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the fibers file.
    output << "##\n";
    output << "# Created by: ElectraPre\n";
    output << "# Format: Fiber X coordinate, Fiber Y coordinate, Fiber Z coordinate\n";
    output << "##\n";
    output << std::setprecision(15);
    for (Eigen::Index i = 0; i != fibers.rows(); ++i)
        output << fibers.coeff(i,0) << ", " << fibers.coeff(i,1) << ", " << fibers.coeff(i,2) << "\n";
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::SaveDistanceField(const Eigen::VectorXd &distance_field, const std::string &output_filename) const
{
    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(output_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }

    // Initialize the output with proper extension.
    std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
    if(!output.good()) {
        std::string error_msg = "ConfigOutput failed to write distance field in file. Check given file: " + output_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the distance file.
    output << "##\n";
    output << "# Created by: ElectraPre\n";
    output << "# Distance scalar field\n";
    output << "# Format: A scalar value per node\n";
    output << "##\n";
    output << std::setprecision(15);
    for (Eigen::Index i = 0; i != distance_field.rows(); ++i)
        output << distance_field.coeff(i) << "\n";
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::SaveDirectionField(const Eigen::MatrixXd &direction_field, const std::string &output_filename) const
{
    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(output_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }

    // Initialize the output with proper extension.
    std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
    if(!output.good()) {
        std::string error_msg = "ConfigOutput failed to write direction field in file. Check given file: " + output_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the fibers file.
    output << "##\n";
    output << "# Created by: ElectraPre\n";
    output << "# Direction vector field\n";
    output << "# Format: A " << DIM << "-D vector per node\n";
    output << "##\n";
    output << std::setprecision(15);
    for (Eigen::Index i = 0; i != direction_field.rows(); ++i) {
        for (Eigen::Index j = 0; j != direction_field.cols()-1; ++j)
            output << direction_field.coeff(i,j) << ", ";
        output << direction_field.coeff(i,direction_field.cols()-1) << "\n";
    }
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_TPP_