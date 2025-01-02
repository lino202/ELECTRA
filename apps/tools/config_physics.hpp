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
   \file config_physics.hpp
   \brief ConfigPhysics class header file.
   \author Konstantinos A. Mountris
   \date 25/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Physics"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigGeo
 * \brief Class to configure the physics of a simulation.
 */
template <short DIM, short CELL_NODES>
class ConfigPhysics {

private:

    std::unordered_map<std::string, ELECTRA::ReactionDiffusionType> react_diff_map_;        /**< Map of the reaction diffusion model types */


protected:

    void SetBidomain(const Parser &parser, const ELECTRA::MeasureUnits &units,
            std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


    void SetMonodomain(const Parser &parser, const ELECTRA::MeasureUnits &units,
            std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


public:

    /**
     * \brief ConfigPhysics object constructor.
     */
    ConfigPhysics();


    /**
     * @brief Destroy the Config Physics object
     * 
     */
    virtual ~ConfigPhysics();


    /**
     * @brief 
     * 
     * @param parser 
     * @param react_diff 
     * @param stream 
     */
    void InitializeReactionDiffusion(const Parser &parser, std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff) const;


    /**
     * @brief Set the Reaction Diffusion object
     * 
     * @param parser 
     * @param units 
     * @param react_diff 
     * @param stream 
     */
    void SetReactionDiffusion(const Parser &parser, const ELECTRA::MeasureUnits &units,
            std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_HPP_

#include "config_physics.tpp"