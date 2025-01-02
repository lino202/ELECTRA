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
   \file config_geo.hpp
   \brief ConfigGeo class header file.
   \author Konstantinos A. Mountris
   \date 13/05/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_GEO_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_GEO_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigGeo
 * \brief Class to configure the geometry of a simulation.
 * \tparam DIM The dimensions of the geometry model of the simulation.
 * \tparam CELL_NODES The number of nodes of the geometry model's cells.
 */
template<short DIM, short CELL_NODES>
class ConfigGeo {

public:

    /**
     * \brief ConfigGeo object constructor.
     */
    ConfigGeo();


    /**
     * \brief ConfigGeo object destructor.
     */
    virtual ~ConfigGeo();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void SetGeometryModel(const Parser &parser, IMP::Mesh<DIM, CELL_NODES> &mesh, IMP::Grid<DIM, CELL_NODES> &grid,
                          IMP::Voronoi<DIM> &voro, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_GEO_HPP_

#include "config_geo.tpp"