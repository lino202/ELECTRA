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
   \file config_units.hpp
   \brief ConfigUnits class header file.
   \author Konstantinos A. Mountris
   \date 13/05/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_UNITS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_UNITS_HPP_

#include "parser.hpp"

#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"

#include <termcolor/termcolor.hpp>

#include <string>
#include <iostream>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigUnits
 * \brief Class to configure the measure units of a simulation.
 */
class ConfigUnits {

public:

    /**
     * \brief ConfigUnits object constructor.
     */
    ConfigUnits();


    /**
     * \brief ConfigUnits object destructor.
     */
    virtual ~ConfigUnits();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void SetReferenceScale(const Parser &parser, ELECTRA::MeasureUnits &units, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_UNITS_HPP_