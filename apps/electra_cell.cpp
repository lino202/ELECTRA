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

#include "ELECTRA/ELECTRA"
#include <IMP/IMP>

#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "tools/config_electrophys.hpp"
#include "tools/parser.hpp"
#include <termcolor/termcolor.hpp>

using namespace ELECTRA;
using namespace APP_ELECTRA;

int main(int argc, char *argv[]) {
    try {

        std::string app_name = "ElectraCell";
        std::cout << termcolor::green << "\nWelcome to " + app_name << termcolor::reset << "\n";
        std::cout << termcolor::bold << "Version:                " << termcolor::reset << ELECTRA_VERSION << "\n";

        if ((argc<10) | (argc>11)){
            std::string err_message = "Give 9 or 10 arguments: output file, model name, cell type, total duration, and stimulus start, duration, cycle length, amplitude and integration dt, for example: ";
            err_message = err_message + "/path/file_name gaur2021 ventricular 1000. 10. 0.5 1000. 80. 0.01 /path/manual_init_file.txt";
            throw std::invalid_argument(Logger::Error(err_message));
        } 

        std::string output_file = argv[1];

        // Time recorder.
        Timer timer;

        // Set measure units.
        MeasureUnits units;
        units.SetRefTimeSIValue(1.e-3);          // time:        ms
        units.SetRefLengthSIValue(1.e-3);        // length:      mm
        units.SetRefConductanceSIValue(1.e-3);   // conductance: mS
        units.SetRefCapacitanceSIValue(1.e-12);  // capacitance: pF
        units.SetRefCurrentSIValue(1.e-3);       // current amplitude: mA

        // Set time step
        double dt = atof(argv[9]);
        dt = dt*units["ms"];

        // Create the cell model calling config_electrophysiology.
        // Input should be the correct name as input in ElectraSim but the cell should be instantiate with the name of the class.
        const std::string ep_model_name = argv[2];
        const std::string cell_type = argv[3];
        std::cout << termcolor::bold << "model:                  " << termcolor::reset << ep_model_name << "\n";
        std::cout << termcolor::bold << "cell type:              " << termcolor::reset << cell_type << "\n";

        ConfigElectrophys config_electrophys;
        std::unique_ptr<ELECTRA::EpBasic> cell = ELECTRA::EpFactory::Create(config_electrophys.GetEpModelType(ep_model_name));
        cell->Initialize(config_electrophys.GetCellType(cell_type));

        if (argc==11){
            std::string manual_init_file = argv[10];
            Parser electra_cell_parser(manual_init_file, app_name);
            std::cout << termcolor::bold << "Using manual init file: " << termcolor::reset << manual_init_file << "\n";
            config_electrophys.ManualCellInitializationElectraCell(electra_cell_parser, manual_init_file, cell);
        }
        
        // Simulation time and steps.
        int total_time = atoi(argv[4]);
        total_time = total_time*units["ms"];
        int steps = static_cast<int>(std::ceil(total_time/dt));
        int steps_for_saving = static_cast <int>(round(0.01/dt));

        //  Set Stimulus
        double stim_start        = atof(argv[5]);
        double stim_dur          = atof(argv[6]);
        double stim_cycle_length = atof(argv[7]);
        double stim_amp          = atof(argv[8]);
        Stimulus stimulus;
        stimulus.SetStart(stim_start*units["ms"]);
        stimulus.SetDuration(stim_dur*units["ms"]);
        stimulus.SetCycleLengths(stim_cycle_length*units["ms"]);
        stimulus.SetAmplitude(stim_amp*units["mA"]);

        // Set the output files
        const std::string out_vm_file_path = output_file+"_Vm.txt";
        std::ofstream out_vm_file(out_vm_file_path, std::ofstream::out);
        const std::string out_manual_init_file_path = output_file+"_manual_init_file.txt";
        std::ofstream out_manual_init_file(out_manual_init_file_path, std::ofstream::out);
        // const std::string out_currs_file_path = output_file+"_currs.txt";
        // std::ofstream out_currs_file(out_currs_file_path, std::ofstream::out);
        
        // Print the remaining info
        std::cout << termcolor::bold << "timestep:               " << termcolor::reset << std::to_string(dt) << " ms\n";
        std::cout << termcolor::bold << "total time:             " << termcolor::reset << std::to_string(total_time) << " ms\n";
        
        std::cout << "\n";
        std::cout << termcolor::magenta << "Stimulus:" << termcolor::reset << "\n";
        std::cout << termcolor::bold << "start:        " << termcolor::reset << std::to_string(stim_start) << " ms\n";
        std::cout << termcolor::bold << "duration:     " << termcolor::reset << std::to_string(stim_dur) << " ms\n";
        std::cout << termcolor::bold << "cycle length: " << termcolor::reset << std::to_string(stim_cycle_length) << " ms\n";
        std::cout << termcolor::bold << "amplitude:    " << termcolor::reset << std::to_string(stim_amp) << " mA\n";

        timer.Reset();
        double stim_current = 0.;
        for (int i = 1; i <= steps; ++i) {

            if (stimulus.IsActive(i*dt, total_time)) { 
                stim_current = stimulus.Amplitude();
            } else {
                stim_current = 0.;
            }

            cell->Compute(cell->V(), dt, stim_current);
            cell->SetV(ALGORITHM::ForwardEuler(cell->V(), dt, cell->dVdt()));

            // Store new state.
            if ((i % steps_for_saving)==0){
                out_vm_file << i*dt << " " << std::setprecision(15) << cell->V() << std::endl;
                // for (int j = 0; j < cell->CurrentNum(); ++j) {
                //     out_currs_file << std::setprecision(15) << cell->Current(j) << " ";
                //     if (j==cell->CurrentNum()-1){
                //         out_currs_file<<std::endl;
                //     }
            }

            // if (i >= 79670){
            //     std::cout << cell->PrintCurrents() << "\n\n";
            // }
            // std::cout <<"step: " << i << "\n" << cell->dVdt() << "\n";
            // std::cout << cell->PrintCurrents() << "\n";
            // std::cout << cell->PrintVariables() << "\n";
            // }

            
            // In case you want a specific current or other state variable you can use this
            // out_stim << i*dt << " " << std::setprecision(15) << cell.Current(1) << std::endl;
        }

        // Save manual init files
        out_manual_init_file << "[Variables]" << std::endl;
        out_manual_init_file << "\n";
        out_manual_init_file << cell->PrintVariables() << std::endl;
        out_manual_init_file << "\n";

        out_manual_init_file << "[Parameters]" << std::endl;
        out_manual_init_file << "\n";
        out_manual_init_file << cell->PrintParameters() << std::endl;
        out_manual_init_file << "\n";

        #ifdef BLOCK_CELL_CURRS
            out_manual_init_file << "[Current Blocks]" << std::endl;
            out_manual_init_file << cell->PrintBlockCoeffs() << std::endl;
            out_manual_init_file << "\n";
        #endif
        

        // Close file.
        out_vm_file.close();
        out_manual_init_file.close();
        
        // Finish program.
        std::cout << "\n";
        std::cout << termcolor::magenta << "Saving:" << termcolor::reset << "\n";
        std::cout << termcolor::bold << "Elapsed time for cell compute: " << termcolor::reset << timer.PrintElapsedTime() << "\n";
        std::cout << termcolor::bold << "Vm saved in:                   " << termcolor::reset <<  out_vm_file_path << "\n";
        std::cout << termcolor::bold << "Manual_init_file saved in:     " << termcolor::reset <<  out_manual_init_file_path << "\n";

        std::cout << "\n";
        std::cout << termcolor::magenta << termcolor::bold;
        std::cout << Logger::Message("The simulation finished successfully. Thank you for using the " + app_name + "app.\n") << termcolor::reset;

    }
    catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    }
    catch (const std::runtime_error &e) {
        std::cerr << "Runtime error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    }
    catch (const std::out_of_range &e) {
        std::cerr << "Out of Range error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    }
    catch (const std::bad_alloc &e) {
        std::cerr << "Bad allocation error:" << termcolor::red << e.what() << termcolor::reset << std::endl;
    }
    catch (...) {
        std::cerr << termcolor::red << "Unknown exception..." << termcolor::reset << std::endl;
    }

    return 0;
}
