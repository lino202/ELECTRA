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
   \file cellstate_exporter.hpp
   \brief CellStateExporter class header file.
   \author Konstantinos A. Mountris
   \date 20/06/2020
*/

#ifndef ELECTRA_EXPORTERS_CELLSTATE_EXPORTER_HPP_
#define ELECTRA_EXPORTERS_CELLSTATE_EXPORTER_HPP_

#include "ELECTRA/engine/electrophysiology/ep_factory.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>


namespace ELECTRA {

/** \addtogroup Exporters \{ */

/**
 * \class CellStateExporter
 * \brief Class that exports the cell states on binary format for probably continuing the simulation.
 */
class CellStateExporter {

    private:
        std::string header_ = std::string("\n###\n# ELECTRA v") + ELECTRA_VERSION + std::string("\n# format: cells_num cell_type cell_vars cell_cur\n###\n");
    
    public:
        /**
         * \brief CellStateExporter constructor.
         */
        CellStateExporter();


        /**
         * \brief CellStateExporter destructor.
         */
        virtual ~CellStateExporter();


        void WriteCellsState(const std::vector<std::unique_ptr<EpBasic>> &cells, const std::string &filename);

        void ReadCellsState(std::vector<std::unique_ptr<EpBasic>> &cells, const std::string &filename);

};


/** \} End of Doxygen Groups*/

} // End of namespace ELECTRA.

#endif  //ELECTRA_EXPORTERS_CELLSTATE_EXPORTER_HPP_
