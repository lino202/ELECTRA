/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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
   \file ensight_exporter.hpp
   \brief EnsightExporter class header file.
   \author Konstantinos A. Mountris & Ricardo M. Rosales
   \date 24/02/2024
*/

#ifndef ELECTRA_EXPORTERS_ENSIGHT_EXPORTER_HPP_
#define ELECTRA_EXPORTERS_ENSIGHT_EXPORTER_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <boost/filesystem.hpp>
#include <boost/range.hpp>

#include <IMP/Vectors>
#include <IMP/Tesselations>

#include <Eigen/Dense>

#include <string>
#include <iostream>
#include <cstddef>
#include <vector>
#include <map>
#include <iterator>
#include <stdexcept>
#include <exception>


namespace ELECTRA {

/** \addtogroup Exporters \{ */

/**
 * \class EnsightExporter
 * \brief Class implemmenting the output of models and their solution fields for post-processing in ParaView.
 * \tparam DIM The dimensions of the domain. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the domain.
 */

template <short DIM, short CELL_NODES=1>
class EnsightExporter {

   public:

      /**
       * \brief EnsightExporter constructor.
       */
      EnsightExporter();
      
      
      /**
       * \brief EnsightExporter destructor.
       */
      virtual ~EnsightExporter();


   private:

      std::string  geom_file_;
      std::string  states_file_;
      std::string  anim_file_;

   public:

      /**
       * \brief EnsightExporter set files.
       */
      void SetFiles(const std::string &geom_file, const std::string &states_file, const std::string &anim_file);
      
      /**
       * \brief  Save a model as Ensight GEOMETRY File (.geo).
       * \param [in] mesh The model's mesh.
       * \return [void]
       */
      void SaveGeo(const std::vector<IMP::Vec<DIM, double>> &nodes, const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells);
      
      /**
       * \brief Save a Scalar Field as Ensight VARIABLE File 
       * \param scalar_field 
       * \param scalar_field_name 
       */
      void SaveStates(const Eigen::VectorXd &scalar_field, const std::string &state_number);

      /**
       * \brief Save Ensight ANIMATION File (.case) gold format.
       * \param [in] steps_num 
       * \param [in] time_inc 
       * \return [void]
       */
      void SaveAnimation(int steps_num, double time_inc);
    
   private:
      /**
       * \brief  Write chars in binary to the ensight files.
       * \param [in] val The array of chars to be written in binary.
       * \param [in] str The binary output stream to be written.
       * \return [void]
       */
      void WriteString_(const char* val, std::ofstream& str)
      {
         const int lineLength = 80;
         char buffer[lineLength] = {0};
         strncpy(buffer, val, lineLength);
         str.write(buffer, lineLength);
      }
      
      /**
       * \brief  Write ints in binary to the ensight files.
       * \param [in] val The int to be written in binary.
       * \param [in] str The binary output stream to be written.
       * \return [void]
       */
      void WriteInt_(const int32_t val, std::ofstream& str)
      {
         str.write(reinterpret_cast<const char*>(&val), sizeof(int32_t));
      }

      /**
       * \brief  Write floats in binary to the ensight files.
       * \param [in] val The float to be written in binary.
       * \param [in] str The binary output stream to be written.
       * \return [void]
       */
      void WriteFloat_(const float val, std::ofstream& str)
      {
         str.write(reinterpret_cast<const char*>(&val), sizeof(float));
      }
};


/*! \} End of Doxygen Groups*/

} // End of namespace ELECTRA.

#endif  //ELECTRA_ENSIGHT_EXPORTER_ENSIGHT_EXPORTER_HPP_

#include "ELECTRA/engine/exporters/ensight_exporter.tpp"
