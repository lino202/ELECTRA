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

#ifndef ELECTRA_APPS_TOOLS_CONFIG_SIM_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_SIM_TPP_

#include "config_sim.hpp"

#include "ELECTRA/ELECTRA"

#include <CLOUDEA/CLOUDEA>
#include <IMP/IMP>

#include <termcolor/termcolor.hpp>

#include <boost/filesystem.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>


namespace APP_ELECTRA {


template<short DIM, short CELL_NODES>
ConfigSim<DIM, CELL_NODES>::ConfigSim()
{}


template<short DIM, short CELL_NODES>
ConfigSim<DIM, CELL_NODES>::~ConfigSim()
{}


template<short DIM, short CELL_NODES>
void ConfigSim<DIM, CELL_NODES>::CheckValid(const Parser &parser)
{
    std::string application = parser.GetValue<std::string>("application");
    std::string app_name = application.substr(0, application.find_first_of(" "));
    std::string app_version = application.substr(application.find("v")+1);

    // Check header validity.
    if (app_name != "ElectraSim") {
        throw std::invalid_argument(ELECTRA::Logger::Error("Header info in configuration file is not consistent with ElectraSim."));
    }

    // Check version validity.
    if (app_version != ELECTRA_VERSION) {
        std::string error_msg = "Could not start simulation. The configuration file version ["
                                + app_version + "] does not match the ElectraSim version [" + ELECTRA_VERSION + "]";
        throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
    }
}


template<short DIM, short CELL_NODES>
void ConfigSim<DIM, CELL_NODES>::Tissue(const Parser &parser, std::ostream &stream)
{
    // Check validity of the provided configuration file.
    this->CheckValid(parser);

    // Initialize time recorder.
    ELECTRA::Timer timer;

    // Initialize physics model.
    ConfigPhysics<DIM, CELL_NODES> config_physics;
    std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> react_diff(nullptr);
    config_physics.InitializeReactionDiffusion(parser, react_diff);

    // Get Method 
    // TODO and NOTE
    // MCM uses the monodomain attribute vout_ which is deprecated now. 
    // So for having MCM again we need to read ensight .ens binary data in the monodomain method FictitiousValuesToReal
    // If we do this we can erase this 'if' and come back to having MCM as was before version 0.6.0
    std::string method = parser.GetValue<std::string>("numerical approximation.method");
    if (method=="mcm"){
        throw std::invalid_argument(ELECTRA::Logger::Error("Invalid approximation method. MCM is temporarily deprecated. Expected: FEM | FPM "));
    }

    // Print relevant information.
    stream << termcolor::green << termcolor::bold << "[*** SIMULATION SET UP ***]\n" << termcolor::reset;
    stream << ELECTRA::Logger::Message("Name:  " + parser.GetValue<std::string>("simulation.name") + "\n");
    stream << ELECTRA::Logger::Message("Scale: " + parser.GetValue<std::string>("simulation.scale") + "\n");
    stream << ELECTRA::Logger::Message("Numerical method: " + method + "\n");

    // Set up measure units.
    ELECTRA::MeasureUnits units;
    ConfigUnits config_units;
    config_units.SetReferenceScale(parser, units, stream);
    std::cout << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Aggregate exporters to the react_diff object as compute happens there and we dynamically save the results
    ELECTRA::EnsightExporter<DIM, CELL_NODES> ens_exporter_tissue;
    ens_exporter_tissue.SetFiles(parser.GetValue<std::string>("output.ensight.tissue.geometry"),
                                parser.GetValue<std::string>("output.ensight.tissue.states"),
                                parser.GetValue<std::string>("output.ensight.tissue.animation"));  
    react_diff->SetEnsightExporterTissue(ens_exporter_tissue);
    // TODO! ens_exporter_cs must be instantiated outside the if otherwise when passing it to react_diff->SetEnsightExporterCS
    // its attributes change for some weird names that rarely change inside setFiles. This for now is not a problem as we 
    ELECTRA::EnsightExporter<DIM, 2> ens_exporter_cs;
    if (parser.HasAttribute("conduction system")) {
        ens_exporter_cs.SetFiles(parser.GetValue<std::string>("output.ensight.conduction system.geometry"),
                            parser.GetValue<std::string>("output.ensight.conduction system.states"),
                            parser.GetValue<std::string>("output.ensight.conduction system.animation"));
        react_diff->SetEnsightExporterCS(ens_exporter_cs);
    }
    
    // Set up model geometry.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** GEOMETRY SET UP ***]\n" << termcolor::reset;
    ConfigGeo<DIM, CELL_NODES> config_geo;
    IMP::Mesh<DIM, CELL_NODES> mesh;
    IMP::Grid<DIM, CELL_NODES> grid;
    IMP::Voronoi<DIM> voro;
    config_geo.SetGeometryModel(parser, mesh, grid, voro, stream);

    // Get tissue node sets and nodes number.
    std::unordered_map<std::string, IMP::NodeSet> tissue_node_sets;
    int tissue_nodes_num = 0;
    if (mesh.NodesNum() != 0) { // FEM case
        tissue_node_sets = mesh.NodeSets();
        tissue_nodes_num = mesh.NodesNum();
    } else if (voro.NodesNum()) { // FPM case
        tissue_node_sets = voro.NodeSets();
        tissue_nodes_num = voro.NodesNum();
    } else { // MCM case
        tissue_node_sets = grid.NodeSets();
        tissue_nodes_num = grid.NodesNum();
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up electric material.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** MATERIAL SET UP ***]\n" << termcolor::reset;
    ConfigElectric<DIM, CELL_NODES> config_electric;
    std::shared_ptr<ELECTRA::ElectricBasic<DIM>> tissue_mat(nullptr);
    config_electric.InitializeMaterial(parser, tissue_mat);
    tissue_mat->SetNodesNum(tissue_nodes_num);
    config_electric.SetMaterialProperties(parser, tissue_node_sets, units, tissue_mat, stream);
    config_electric.ComputeDiffusionTensors(tissue_mat, stream);
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up conduction system.
    timer.Reset();
    int conduct_sys_nodes_num = 0;
    std::vector<std::unique_ptr<ELECTRA::EpBasic>> cs_nodal_cells;
    std::vector<ELECTRA::EpVaryingParams> cs_cell_varying_params;
    if (parser.HasAttribute("conduction system")) {
        stream << termcolor::green << termcolor::bold << "[*** CONDUCTION SYSTEM SET UP ***]\n" << termcolor::reset;
        ConfigConductSys<DIM,CELL_NODES> config_cs;
        if (mesh.NodesNum() != 0) { // FEM case
            config_cs.SetConductSystem(parser, units, mesh.Nodes(), react_diff->EditConductSystem(), stream);
        } else if (voro.NodesNum()) { // FPM case
            config_cs.SetConductSystem(parser, units, voro.Nodes(), react_diff->EditConductSystem(), stream);
        } else { // MCM case
            config_cs.SetConductSystem(parser, units, grid.Nodes(), react_diff->EditConductSystem(), stream);
        }
        std::cout << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n";
        std::cout << "\n";

        // Set the conduction system nodes num paramer.
        conduct_sys_nodes_num = react_diff->ConductSystem().NodesNum();
    }

    // Set up electrophysiology.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** ELECTROPHYSIOLOGY SET UP ***]\n" << termcolor::reset;
    ConfigElectrophys config_electrophys;
    std::vector<std::unique_ptr<ELECTRA::EpBasic>> tissue_nodal_cells;
    std::vector<ELECTRA::EpVaryingParams> tissue_cell_varying_params;
    stream << termcolor::magenta << ELECTRA::Logger::Message("Tissue electrophysiology") << termcolor::reset << "\n";
    config_electrophys.SetCellElectrophysiology(parser, tissue_node_sets, "tissue", tissue_nodes_num, tissue_nodal_cells, tissue_cell_varying_params, stream);
    if (parser.HasAttribute("conduction system")) {
        stream << termcolor::magenta << ELECTRA::Logger::Message("Conduction system electrophysiology") << termcolor::reset << "\n";
        config_electrophys.SetCellElectrophysiology(parser, react_diff->ConductSystem().NodeSets(), "conduction system",
                conduct_sys_nodes_num, cs_nodal_cells, cs_cell_varying_params, stream);
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Collect all the nodal cells of the model
    std::vector<std::unique_ptr<ELECTRA::EpBasic>> nodal_cells = std::move(tissue_nodal_cells);
    std::vector<ELECTRA::EpVaryingParams> nodal_cell_varying_params = tissue_cell_varying_params;
    if (conduct_sys_nodes_num != 0) {
        nodal_cells.insert(nodal_cells.end(), std::make_move_iterator(cs_nodal_cells.begin()), std::make_move_iterator(cs_nodal_cells.end()));
        nodal_cell_varying_params.insert(nodal_cell_varying_params.end(), std::make_move_iterator(cs_cell_varying_params.begin()), std::make_move_iterator(cs_cell_varying_params.end()));
    }

    // After all nodal_cells are set up we add the load cell states if we need to start from a previous checkpoint
    if (parser.HasAttribute("simulation.load cells state")) {
        std::string cell_state_file = parser.GetValue<std::string>("simulation.load cells state");
        std::string cell_state_ext = std::filesystem::path(cell_state_file).extension();
        if (cell_state_ext != ".elc") { cell_state_file += ".elc"; }

        ELECTRA::CellStateExporter cellstate_export;
        cellstate_export.ReadCellsState(nodal_cells, cell_state_file);
        stream << ELECTRA::Logger::Message("Loaded cells state: " + cell_state_file + "\n");
    }

    // Set up stimuli.
    timer.Reset();
    ConfigStimuli config_stim;
    std::vector<ELECTRA::Stimulus> stimuli, tissue_stimuli, cs_stimuli;
    if (parser.HasAttribute("tissue.stimuli") || parser.HasAttribute("conduction system.stimuli")) {
        stream << termcolor::green << termcolor::bold << "[*** STIMULI SET UP ***]\n" << termcolor::reset;
    }
    if (parser.HasAttribute("tissue.stimuli")) {
        stream << termcolor::magenta << ELECTRA::Logger::Message("Tissue stimuli") << termcolor::reset << "\n";
        config_stim.SetStimuli(parser, "tissue", tissue_node_sets, tissue_nodes_num+conduct_sys_nodes_num, units, tissue_stimuli, stream);
        stimuli.insert(stimuli.end(), std::make_move_iterator(tissue_stimuli.begin()), std::make_move_iterator(tissue_stimuli.end()));

    }
    if (parser.HasAttribute("conduction system.stimuli")) {
        stream << termcolor::magenta << ELECTRA::Logger::Message("Conduction system stimuli") << termcolor::reset << "\n";
        config_stim.SetStimuli(parser, "conduction system", react_diff->ConductSystem().NodeSets(), tissue_nodes_num+conduct_sys_nodes_num, units, cs_stimuli, stream);
        for (auto &cs_stim : cs_stimuli) {
            cs_stim.AddStimulatedNodesPadding(tissue_nodes_num);
        }
        stimuli.insert(stimuli.end(), std::make_move_iterator(cs_stimuli.begin()), std::make_move_iterator(cs_stimuli.end()));
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up numerical approximation.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** NUMERICAL APPROXIMATION SET UP ***]\n" << termcolor::reset;
    ConfigApproximation<DIM, CELL_NODES> config_approx;
    CLOUDEA::ApproxType approx_type;
    CLOUDEA::Fpm<DIM> fpm_approx;
    std::unique_ptr<CLOUDEA::Mfree<DIM>> mcm_approx;
    IMP::NodeSet neumann_set;
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);
    if (method == "fem") {
        approx_type = CLOUDEA::ApproxType::fem;
        config_approx.SetFemApproximation(stream);
    } else if (method == "fpm") {
        approx_type = CLOUDEA::ApproxType::fpm;
        config_approx.SetFpmApproximation(parser, voro, fpm_approx, stream);
    } else if (method == "mcm") {
        approx_type = CLOUDEA::ApproxType::mcm;
        config_approx.SetMcmApproximation(parser, grid, mcm_approx, neumann_set, stream);
    } else {
        throw std::invalid_argument(ELECTRA::Logger::Error("Invalid simulation method in configuration file. Expected: FEM | FPM | MCM"));
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up physics.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** PHYSICS SOLUTION ***]\n" << termcolor::reset;
    config_physics.SetReactionDiffusion(parser, units, react_diff, stream);
    react_diff->SetNodalCells(std::move(nodal_cells));
    react_diff->SetCellVarParamGroups(nodal_cell_varying_params);
    // Assemble react_diff model matrices.
    if (approx_type == CLOUDEA::ApproxType::fem) {
        stream << ELECTRA::Logger::Message("Assembling FEM matrices... ");
        react_diff->AssembleMatrices(mesh, tissue_mat);
        stream << "OK\n";
    } else if (approx_type == CLOUDEA::ApproxType::fpm) {
        stream << ELECTRA::Logger::Message("Assembling FPM matrices... ");
        react_diff->AssembleMatrices(voro, tissue_mat, fpm_approx);
        stream << "OK\n";
    } else if (approx_type == CLOUDEA::ApproxType::mcm) {
        stream << ELECTRA::Logger::Message("Assembling MCM matrices... ");
        react_diff->AssembleMatrices(grid, tissue_mat, mcm_approx, neumann_set);
        stream << "OK\n";
    } else {
        std::string error_msg = "Could not assemble react_diff model matrices due to not supported simulation method. Supported: FEM | FPM | MCM";
        throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
    }

    // Set up reaction diffusion time stepping.
    react_diff->SetupTimeStepping();
    if (react_diff->DtCritical() < std::numeric_limits<double>::max()-2.*std::numeric_limits<double>::epsilon())
        stream << ELECTRA::Logger::Message("Critical time step: ") << react_diff->DtCritical() << " ms\n";
    stream << Logger::Message("Used time step: " + std::to_string(react_diff->Dt())) << " ms\n";
    stream << Logger::Message("Total time steps: ") << react_diff->SimulationSteps() << "\n";


    // SOLVE reaction diffusion physics
    react_diff->Compute(tissue_mat, stimuli);
    // if (approx_type == CLOUDEA::ApproxType::mcm) {
        // react_diff->FictitiousValuesToReal(mcm_approx);
        // stream << ELECTRA::Logger::Message("Applied conversion of fictitious to real nodal values.\n");
    // }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up post-process.
    // NOTE: The postprocessing is deactivated as it uses _vout which is deprecated since version 0.6.0. Also calculation of APD was not working.
    // TODO: Possible solutions:
    // 1- Make the user to set time start and end load a matrix Vout from binary saved states that can be use in postprocessing. 
    // 2- Another solution which is better, is to leave ElectraSim app to calculate so being independent from post_processing. So implement this in ElectraPost 
    // which is just a skeleton for now.
    // This decision to not have _vout and then post_process here in ElectraSim generated substancial changes and simplification of Output generation 

    // timer.Reset();
    // stream << termcolor::green << termcolor::bold << "[*** POST PROCESSING ***]\n" << termcolor::reset;
    // ELECTRA::PostProcess post_process;
    // ConfigPostProcess<DIM, CELL_NODES> config_post_process;
    // config_post_process.PerformPostProcess(parser, react_diff, tissue_node_sets, post_process, stream);
    // stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up output.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** OUTPUT ***]\n" << termcolor::reset;
    ConfigOutput<DIM, CELL_NODES> config_output;
    if (approx_type == CLOUDEA::ApproxType::fem || approx_type == CLOUDEA::ApproxType::fpm) {
        config_output.OutputGeneration(parser, mesh.Nodes(), mesh.Cells(), react_diff, stream);
    } else {
        config_output.OutputGeneration(parser, grid.Nodes(), grid.GhostCells(), react_diff, stream);
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

}


} // End of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_SIM_TPP_