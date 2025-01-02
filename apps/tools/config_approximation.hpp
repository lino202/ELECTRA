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
   \file config_approximation.hpp
   \brief ConfigApproximation class header file.
   \author Konstantinos A. Mountris
   \date 24/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"

#include <IMP/IMP>
#include <CLOUDEA/CLOUDEA>

#include <termcolor/termcolor.hpp>
#include <Eigen/Dense>

#include <string>
#include <filesystem>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <functional>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigApproximation
 * \brief Class to configure the numerical approximation method of a simulation.
 */
template <short DIM, short CELL_NODES>
class ConfigApproximation {

public:

    /**
     * \brief Construct a new Config Approximation object
     */
    ConfigApproximation();


    /**
     * \brief Destroy the Config Approximation object
     */
    virtual ~ConfigApproximation();


    /**
     * @brief Set the Fem Approximation object
     * 
     * @param stream 
     */
    void SetFemApproximation(std::ostream &stream) const;


    void SetFpmApproximation(const Parser &parser, const IMP::Voronoi<DIM> &voro, CLOUDEA::Fpm<DIM> &fpm_approx, std::ostream &stream) const;


    void SetMcmApproximation(const Parser &parser, const IMP::Grid<DIM, CELL_NODES> &grid, std::unique_ptr<CLOUDEA::Mfree<DIM>> &mcm_approx,
                             IMP::NodeSet &neumann_nset, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_

#include "config_approximation.tpp"