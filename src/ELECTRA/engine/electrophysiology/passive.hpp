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

/**
   \file passive.hpp
   \brief Passive class header file.
   \date 21/08/2025
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_PASSIVE_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_PASSIVE_HPP_

#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"
#include "ELECTRA/engine/utilities/algorithm.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <limits>


namespace ELECTRA {

/** \addtogroup Electrophysiology \{ */


/**
 * \namespace ELECTRA::PassVar
 * \brief Encaplulates enumeration to access Passive model variables.
 */
namespace PassVar {
enum { 
        v,           /**< Cell membrane potential */
        dvdt        /**< Cell membrane potential's time derivative */
    };
} // End of namespace PassVar


/**
 * \namespace ELECTRA::PassPrm
 * \brief Encaplulates enumeration to access Passive model parameters.
 */
namespace PassPrm {
enum { 
        Cm,             /**< Membrane capacitance [pF] */
        R,              /**< Membrane resistance [GOhm] */
        Vrest,          /**< Resting potential [mV] */
    };
} // End of namespace PassPrm


/**
 * \namespace ELECTRA::PassCur
 * \brief Encaplulates enumeration to access Passive model currents.
 */
namespace PassCur {
enum { 
        Iion        /**< Total ionic current */
    };
} // End of namespace PassCur


/**
 * \class Passive
 * \brief Passive model initialized for a fibroblast as defined by Maccannell (doi: 10.1529/biophysj.106.101410).
 */
class Passive : public EpBasic
{
protected:

    virtual void SetDataMapping();


public:

    /**
     * \brief The default constructor of the Passive class.
     */
    Passive();


    /**
     * \brief The default destructor of the Passive class.
     */
    virtual ~Passive();


    /**
     * \brief Initialize the variables and parameters of the Passive model. 
     * \return [void]
     */
    void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Passive model.
     * \return [void]
     */
    void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the potential.
     * \param [in] v The potential value.
     * \return [void]
     */
    inline void SetV(double v) { this->var_[PassVar::v] = v; }


    /**
     * \brief Print to std::string the model's variables and their values.
     * \return [std::string] The model's variables and their values.
    */
    virtual std::string PrintVariables() const;


    /**
     * \brief Print to std::string the model's parameters and their values.
     * \return [std::string] The model's parameters and their values.
    */
    virtual std::string PrintParameters() const;


    /**
     * \brief Print to std::string the model's currents and their values.
     * \return [std::string] The model's currents and their values.
    */
    virtual std::string PrintCurrents() const;


    /**
     * \brief Get the potential value of the model.
     * \return [double] The potential value.
     */
    inline double V() const { return this->var_[PassVar::v]; }


    /**
     * \brief Get the potential time derivative of the model.
     * \return [double] The potential time derivative.
     */
    inline double dVdt() const { return this->var_[PassVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[PassCur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_PASSIVE_HPP_