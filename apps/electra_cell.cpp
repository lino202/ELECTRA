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

#include "ELECTRA/ELECTRA"
#include <IMP/IMP>

#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "tools/config_electrophys.hpp"
#include <termcolor/termcolor.hpp>

using namespace ELECTRA;
using namespace APP_ELECTRA;

int main(int argc, char *argv[]) {
    try {

        std::cout << termcolor::green << "\nWelcome to ElectraCell" << termcolor::reset << "\n";

        if (argc!=9){
            std::string err_message = "Give 8 arguments: output file, model name, cell type, total duration, and stimulus start, duration, cycle length and amplitude, for example: ";
            err_message = err_message + "/path/file_name gaur2021 endo 1000. 10. 0.5 1000. 80.";
            throw std::invalid_argument(Logger::Error(err_message));
        } 

        std::string output_file = argv[1];

        // Time recorder.
        Timer timer;

        // Set measure units.
        MeasureUnits units;
        units.SetRefTimeSIValue(1.e-3);          // time:        ms
        units.SetRefLengthSIValue(1.e-2);        // length:      cm
        units.SetRefConductanceSIValue(1.e-3);   // conductance: mS
        units.SetRefCapacitanceSIValue(1.e-9);  // capacitance: nF

        // Set cell capacitance and time step
        double dt = 0.01*units["ms"];

        // Create the cell model.
        // Onput should be the correct name as input in opencarp but the cell should be instantiate with the name 
        // of the class.
        const std::string ep_model_name = argv[2];
        const std::string cell_type = argv[3];
        ConfigElectrophys config_electrophys;
        std::unique_ptr<ELECTRA::EpBasic> cell = ELECTRA::EpFactory::Create(config_electrophys.GetEpModelType(ep_model_name));
        cell->Initialize(config_electrophys.GetCellType(cell_type));
        


        // Simulation time and steps.
        int total_time = atoi(argv[4]);
        total_time = total_time*units["ms"];
        int steps = static_cast<int>(std::ceil(total_time/dt));

        //OLD left for reference 
        // std::vector<double> cycle_lengths({630*units["ms"], 630*units["ms"], 620*units["ms"], 620*units["ms"],
        //                                    610*units["ms"], 610*units["ms"], 600*units["ms"], 600*units["ms"],
        //                                    590*units["ms"], 590*units["ms"], 580*units["ms"], 580*units["ms"],
        //                                    570*units["ms"], 570*units["ms"], 560*units["ms"], 560*units["ms"],
        //                                    550*units["ms"], 550*units["ms"], 540*units["ms"], 540*units["ms"],
        //                                    500*units["ms"], 500*units["ms"], 500*units["ms"], 500*units["ms"],
        //                                    500*units["ms"], 500*units["ms"], 500*units["ms"], 500*units["ms"],
        //                                    500*units["ms"], 500*units["ms"], 540*units["ms"], 540*units["ms"],
        //                                    550*units["ms"], 550*units["ms"], 560*units["ms"], 560*units["ms"],
        //                                    570*units["ms"], 570*units["ms"], 580*units["ms"], 580*units["ms"],
        //                                    590*units["ms"], 590*units["ms"], 600*units["ms"], 600*units["ms"],
        //                                    610*units["ms"], 610*units["ms"], 620*units["ms"], 620*units["ms"],
        //                                    630*units["ms"], 630*units["ms"]});
        // Stimulus stimulus;
        // stimulus.SetStart(0.*units["ms"]);
        // stimulus.SetDuration(0.5*units["ms"]);
        // // stimulus.SetCycleLengths(cycle_lengths);
        // stimulus.SetCycleLengths(1000.);
        // stimulus.SetAmplitude(80);


        //  Set Stimulus
        double stim_start        = atof(argv[5]);
        double stim_dur          = atof(argv[6]);
        double stim_cycle_length = atof(argv[7]);
        double stim_amp          = atof(argv[8]);
        Stimulus stimulus;
        stimulus.SetStart(stim_start*units["ms"]);
        stimulus.SetDuration(stim_dur*units["ms"]);
        stimulus.SetCycleLengths(stim_cycle_length);
        stimulus.SetAmplitude(stim_amp);

        // Set the output files
        const std::string out_vm_file_path = output_file+"_Vm.txt";
        std::ofstream out_vm_file(out_vm_file_path, std::ofstream::out);
        const std::string out_manual_init_file_path = output_file+"_manual_init_file.txt";
        std::ofstream out_manual_init_file(out_manual_init_file_path, std::ofstream::out);
        
        // Print some info
        std::cout << termcolor::bold << "version:    " << termcolor::reset << ELECTRA_VERSION << "\n";
        std::cout << termcolor::bold << "timestep:   " << termcolor::reset << std::to_string(dt) << " ms\n";
        std::cout << termcolor::bold << "model:      " << termcolor::reset << ep_model_name << "\n";
        std::cout << termcolor::bold << "cell type:  " << termcolor::reset << cell_type << "\n";
        std::cout << termcolor::bold << "total time: " << termcolor::reset << std::to_string(total_time) << " ms\n";
        
        std::cout << "\n";
        std::cout << termcolor::magenta << "Stimulus:" << termcolor::reset << "\n";
        std::cout << termcolor::bold << "start:        " << termcolor::reset << std::to_string(stim_start) << " ms\n";
        std::cout << termcolor::bold << "duration:     " << termcolor::reset << std::to_string(stim_dur) << " ms\n";
        std::cout << termcolor::bold << "cycle length: " << termcolor::reset << std::to_string(stim_cycle_length) << " ms\n";
        std::cout << termcolor::bold << "amplitude:    " << termcolor::reset << std::to_string(stim_amp) << " [unit]\n";  //unit depends on the model, commonly is mA



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

            // if (i==1) {
            //     std::cout <<"step: " << i << "\n" << cell.dVdt() << "\n";
            //     std::cout << cell.PrintCurrents() << "\n";
            //     break;
            // }

            // Store new state.
            out_vm_file << i*dt << " " << std::setprecision(15) << cell->V() << std::endl;
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
        // Currents are not neccesary for setting the manual initialization
        // out_manual_init_file << "[Currents]" << std::endl;
        // out_manual_init_file << cell->PrintCurrents() << std::endl;
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
        std::cout << termcolor::bold << "Vm saved in: " << termcolor::reset <<  out_vm_file_path << "\n";
        std::cout << termcolor::bold << "Manual_init_file saved in: " << termcolor::reset <<  out_manual_init_file_path << "\n";

        std::cout << "\n";
        std::cout << termcolor::magenta << termcolor::bold;
        std::cout << Logger::Message("The simulation finished successfully. Thank you for using the ElectraCell app.\n") << termcolor::reset;

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
