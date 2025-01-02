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
   \file config_conduct_sys.hpp
   \brief ConfigConductSys class header file.
   \author Konstantinos A. Mountris
   \date 20/10/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/ConductionSystem"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <iostream>
#include <unordered_map>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigConductSys
 * \brief Class to configure the conduction system of the heart.
 */
template <short DIM, short CELL_NODES>
class ConfigConductSys {

protected:

    void SetCsGeometry(const Parser &parser, const ELECTRA::MeasureUnits &units, const std::vector<IMP::Vec<DIM, double>> &tissue_nodes,
            ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const;


    void SetCsDiffusivity(const Parser &parser, const ELECTRA::MeasureUnits &units, ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const;


    void ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute,
            const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int values_num, std::vector<double> &values) const;

public:

    /**
     * \brief ConfigConductSys object constructor.
     */
    ConfigConductSys();


    virtual ~ConfigConductSys();


    void SetConductSystem(const Parser &parser, const ELECTRA::MeasureUnits &units,
            const std::vector<IMP::Vec<int(DIM), double>> &tissue_nodes, ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_HPP_

#include "config_conduct_sys.tpp"