/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
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

};


/** \} End of Doxygen Groups*/

} // End of namespace ELECTRA.

#endif  //ELECTRA_EXPORTERS_CELLSTATE_EXPORTER_HPP_
