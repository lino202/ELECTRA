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


#include "ELECTRA/engine/electrophysiology/passive.hpp"


namespace ELECTRA {


void Passive::SetDataMapping()
{
    using namespace PassVar;
    using namespace PassPrm;
    using namespace PassCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);

    // Set parameters mapping.
    this->mapped_data_["Cm"] = static_cast<std::size_t>(Cm);     // Membrane capacitance.
    this->mapped_data_["R"] = static_cast<std::size_t>(R);      // Membrane resistance.
    this->mapped_data_["Vrest"] = static_cast<std::size_t>(Vrest); // Resting potential.

    // Set currents mapping.
    // There is no other current than the total ionic current in the passive model.

}


Passive::Passive()
{
    this->model_type_ = EpModelType::Passive;
    this->dt_stable_ = 0.02;
    this->var_.resize(2, 0.);    // Variables: v, dvdt
    this->prm_.resize(3, 0.);
    this->cur_.resize(1, 0.);    // Currents: Iion

    // Set mapped data.
    this->SetDataMapping();
}


Passive::~Passive()
{}


void Passive::Initialize(CellType cell_type)
{
    using namespace PassVar;
    using namespace PassPrm;

    // Check required cell type.
    if (cell_type != CellType::passive) {
        std::string error_str = Logger::Error("Could not initialize Passive ap model. Expected CellType::passive.");
        throw std::invalid_argument(error_str);
    }

    //Initialize the model's data.
    this->var_.clear();          this->var_.resize(2, 0.);
    this->prm_.clear();          this->prm_.resize(3, 0.);
    this->cur_.clear();          this->cur_.resize(1, 0.);

    // Initialize the model's data, which initialized for a fibroblast as defined by Maccannell (doi: 10.1529/biophysj.106.101410)
    // Set variables.
    this->var_[v]    = -49.6;    // Membrane potential in mV.
    this->var_[dvdt] = 0.;

    // Set parameters.
    this->prm_[Cm]      = 6.3;  // Membrane capacitance in pF.
    this->prm_[R]       = 10.7;   // Membrane resistance in GOhm.
    this->prm_[Vrest]   = -49.6;  // Resting potential in mV.

}


void Passive::Compute(double v_new, double dt, double stim_current)
{
    using namespace PassVar;
    using namespace PassPrm;
    using namespace PassCur;

    // Total ionic current.
    this->cur_[PassCur::Iion] = (v_new - this->prm_[Vrest]) / this->prm_[R];

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[PassCur::Iion] - stim_current) / this->prm_[Cm]; ;
    
}


std::string Passive::PrintVariables() const
{
    using namespace PassVar;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v]    << "\n";
    oss << "dvdt = " << this->var_[dvdt];
    return oss.str();
}


std::string Passive::PrintParameters() const
{
    using namespace PassPrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Cm = " << this->prm_[Cm] << "\n";
    oss << "R = " << this->prm_[R] << "\n";
    oss << "Vrest = " << this->prm_[Vrest];
    return oss.str();
}


std::string Passive::PrintCurrents() const
{
    using namespace PassCur;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Iion = " << this->cur_[PassCur::Iion];
    return oss.str();
}


} // End of namespace ELECTRA