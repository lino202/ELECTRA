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
   \file config_electrophys.hpp
   \brief ConfigElectrophys class header file.
   \author Konstantinos A. Mountris
   \date 20/05/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_ELECTROPHYS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_ELECTROPHYS_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Electrophysiology"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <set>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigElectrophys
 * \brief Class to configure the electophysiology of a simulation.
 * \tparam DIM The dimensions of the geometry model of the simulation.
 * \tparam CELL_NODES The number of nodes of the geometry model's cells.
 */
class ConfigElectrophys {

private:

    std::unordered_map<std::string, ELECTRA::EpModelType> ep_model_map_;    /**< Map of the available electrophysiology models */

    std::unordered_map<std::string, ELECTRA::CellType> cell_type_map_;      /**< Map of the available cell types */


protected:

    /**
     * @brief Set the Cell Models From Node Sets object
     * 
     * @param parser 
     * @param ep_model_map 
     * @param cell_type_map 
     * @param node_sets 
     * @param body_type 
     * @param nodal_cells_num 
     * @param nodal_cells 
     * @param cell_varying_param_groups 
     * @param stream 
     */
    void SetCellModelsFromNodeSets(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &node_sets,
            const std::string &body_type, int nodal_cells_num, std::vector<std::unique_ptr<ELECTRA::EpBasic>> &nodal_cells,
            std::vector<ELECTRA::EpVaryingParams> &cell_varying_param_groups, std::ostream &stream) const;


    /**
     * @brief Set the Cell Models Individually object
     * 
     * @param parser 
     * @param ep_model_map 
     * @param cell_type_map 
     * @param body_type 
     * @param nodal_cells_num 
     * @param nodal_cells 
     * @param cell_varying_param_groups 
     * @param stream 
     */
    void SetCellModelsIndividually(const Parser &parser, const std::string &body_type, int nodal_cells_num,
            std::vector<std::unique_ptr<ELECTRA::EpBasic>> &nodal_cells,
            std::vector<ELECTRA::EpVaryingParams> &cell_varying_param_groups, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param init_file 
     * @param cell 
     */
    void ManualCellInitialization(const Parser &parser, const std::string &init_file, std::unique_ptr<ELECTRA::EpBasic> &cell) const;


    /**
     * @brief Set the Ep Varying Params object
     * 
     * @param parser 
     * @param filename 
     * @param ep_model_name 
     * @param ep_var_params 
     */
    void SetEpVaryingParams(const Parser &parser, const std::string &filename, const std::string &ep_model_name,
            ELECTRA::EpVaryingParams &ep_var_params) const;

public:

    /**
     * \brief ConfigElectrophys object constructor.
     */
    ConfigElectrophys();


    /**
     * \brief ConfigElectrophys object destructor.
     */
    virtual ~ConfigElectrophys();


    /**
     * @brief Set the Cell Electrophysiology object, called from ElectraSim
     * 
     * @param parser 
     * @param node_sets 
     * @param body_type 
     * @param nodal_cells_num 
     * @param nodal_cells 
     * @param cell_varying_param_groups 
     * @param stream 
     */
    void SetCellElectrophysiology(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &node_sets,
            const std::string &body_type, int nodal_cells_num, std::vector<std::unique_ptr<ELECTRA::EpBasic>> &nodal_cells,
            std::vector<ELECTRA::EpVaryingParams> &cell_varying_param_groups, std::ostream &stream) const;

    /**
     * @brief Manuallt init the cell states and params, called from ElectraCell
     * 
     * @param parser 
     * @param manual_init_file
     * @param cell
     */
    void ManualCellInitializationElectraCell(const Parser &parser, const std::string &manual_init_file, std::unique_ptr<ELECTRA::EpBasic> &cell){
        this->ManualCellInitialization(parser, manual_init_file, cell);
    }
    
    ELECTRA::EpModelType GetEpModelType(const std::string &ep_model_name){
        return this->ep_model_map_.at(ep_model_name);
    }

    ELECTRA::CellType GetCellType(const std::string &cell_type_name){
        return this->cell_type_map_.at(cell_type_name);
    }


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_ELECTROPHYS_HPP_