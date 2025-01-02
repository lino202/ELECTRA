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
   \file config_fibers.hpp
   \brief ConfigFibers class header file.
   \author Konstantinos A. Mountris
   \date 02/07/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_FIBERS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_FIBERS_HPP_

#include "parser.hpp"

#include <IMP/IMP>
#include <CLOUDEA/CLOUDEA>

#include "ELECTRA/Utilities"
#include "ELECTRA/Fibers"

#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>

using namespace ELECTRA;

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigFibers
 * \brief Class to configure the preprocessing of cardiac fibers.
 */
template <short DIM, short CELL_NODES>
class ConfigFibers {

public:

    /**
     * \brief ConfigFibers object constructor.
     */
    ConfigFibers();


    virtual ~ConfigFibers();


    void SetVentriFiberRules(const Parser &parser, VentriFiberRules &ventri_fiber_rules, std::ostream &stream) const;


    void SetVentriTags(const Parser &parser, VentriTags &ventri_tags, std::ostream &stream) const;


    void SetAtriFiberRules(const Parser &parser, AtriFiberRules &atri_fiber_rules, std::ostream &stream) const;


    void SetAtriTags(const Parser &parser, AtriTags &atri_tags, std::ostream &stream) const;


    void ComputeVentriFibers(const Parser &parser, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
                             const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &tags, const VentriFiberRules &rules,
                             Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const;


    void ComputeAtriFibers(const Parser &parser, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
                           const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &tags, const AtriFiberRules &rules,
                           Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_FIBERS_HPP_

#include "config_fibers.tpp"