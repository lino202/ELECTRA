/*
 * IMP. Image and Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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
   \file logger.hpp
   \brief Logger class header file.
   \author Konstantinos A. Mountris
   \date 16/01/2018
*/

#ifndef IMP_ENGINE_UTILITIES_LOGGER_HPP_
#define IMP_ENGINE_UTILITIES_LOGGER_HPP_

#include <string>

namespace IMP {

/** \addtogroup Utilities \{ */


/**
 * \class Logger
 * \brief Class implemmenting output messages for logging of IMP.
 */

class Logger
{
public:

    /**
     * \brief Logger default constructor.
     */
    Logger() {}


    /**
     * \brief Logger default destructor.
     */
    virtual ~Logger() {}


    /**
     * \brief Log an IMP message.
     * \param [in] msg The message to be logged.
     * \return [std::string] The logged message.
     */
    inline static std::string Message(const std::string &msg) {
        return "[IMP] " + msg;
    }


    /**
     * \brief Log an IMP error.
     * \param [in] err The error to be logged.
     * \return [std::string] The logged error.
     */
    inline static std::string Error(const std::string &err) {
        return "[IMP ERROR] " + err;
    }


    /**
     * \brief Log an IMP warning.
     * \param [in] wrng The warning to be logged.
     * \return [std::string] The logged warning.
     */
    inline static std::string Warning(const std::string &wrng) {
        return "[IMP WARNING] " + wrng;
    }


};

/** \} End of Doxygen Groups*/

} //end of namespace IMP

#endif //IMP_ENGINE_UTILITIES_LOGGER_HPP_
