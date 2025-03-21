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
   \file config_output.hpp
   \brief ConfigOutput class header file.
   \author Konstantinos A. Mountris
   \date 27/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Physics"
#include "ELECTRA/Fibers"
#include "ELECTRA/PostProcess"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <set>

using namespace ELECTRA;

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigOutput
 * \brief Class to configure the output of simulation results.
 */
template<short DIM, short CELL_NODES>
class ConfigOutput {

protected:


    /**
     * @brief Save results to the Ensight gold format
     * 
     * @param nodes 
     * @param cells 
     * @param react_diff 
     */
    void OutputToEnsight(const std::vector<IMP::Vec<int(DIM), double>> &nodes, const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells,
            const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff) const;


    /**
     * @brief Save final cell states for probably continuing the simulation
     * 
     * @param parser 
     * @param monodomain 
     * @param stream 
     */
    void OutputToCellStates(const Parser &parser, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff,
            std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param fibers 
     * @param output_filename 
     */
    void SaveFibers(const Eigen::MatrixXd &fibers, const std::string &output_filename) const;


    void SaveDistanceField(const Eigen::VectorXd &distance_field, const std::string &output_filename) const;

    void SaveDirectionField(const Eigen::MatrixXd &direction_field, const std::string &output_filename) const;


public:

    /**
     * \brief ConfigOutput object constructor.
     */
    ConfigOutput();


    /**
     * \brief ConfigOutput object destructor.
     */
    virtual ~ConfigOutput();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void OutputGeneration(const Parser &parser, const std::vector<IMP::Vec<DIM,double>> &nodes,
            const std::vector<IMP::Cell<DIM,CELL_NODES>> &cells, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param mesh 
     * @param voro 
     * @param ldrbm 
     * @param stream 
     */
    void OutputFibers(const Parser &parser, IMP::Mesh<DIM,CELL_NODES> mesh, IMP::Voronoi<DIM> voro,
            const Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const;
};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_HPP_

#include "config_output.tpp"