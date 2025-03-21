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
 * \file parser.hpp
 * \brief Parser class header file
 * \author Konstantinos A. Mountris
 * \date 10/05/2021
 */

#ifndef ELECTRA_APPS_TOOLS_PARSER_HPP_
#define ELECTRA_APPS_TOOLS_PARSER_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <filesystem>

/** \addtogroup Application-Tools \{ */

/**
 * \namespace APP_ELECTRA
 * \brief Collects tools for configuring and executing ELECTRA applications.
 */
namespace APP_ELECTRA {

class Parser
{
private:

    nlohmann::json json_;       /**< The parser of the Electra application JSON file  */

    std::string parent_path_;     /**< The absolute path of the Electra application JSON file */


protected:

    /**
     * \brief Removes comment lines from a string.
     * Comments are expected in the style of C++ and JavaScript.
     * \param [in] str The string containing comment lines.
     * \return [std::string] A new string whith the comment lines removed.
     */
    std::string RemoveCommentLines(const std::string &str);

public:

    /**
     * \brief Construct a new Parser object, main Parser object called from ElectraSim
     * \param [in] json_filename The filename of the configuration file of the ELECTRA application including the path.
     */
    Parser(const std::string &json_filename);

    /**
     * \brief Construct a new Parser object, for being called from ElectraCell 
     * TODO This extra construct makes possible to use config_electrophys.ManualCellInitialization
     * cause this needs to have the parser as an argument. So, for not overloading the constructor I added app_name as parameter 
     * which is certaintly silly, might be improved
     * \param [in] json_filename The filename of the configuration file of the ELECTRA application including the path.
     */
    Parser(const std::string &manual_init_file, const std::string &app_name);


    /**
     * \brief Destroy the Parser object.
     */
    virtual ~Parser();


    /**
     * \brief Get the value of an attribute in the configuration file.
     * \tparam T The type of the value to be retrieved.
     * \param [in] attribute The name of the requested attribute including the complete path hierarchy.
     * \return [T] The value of an attribute in the configuration file.
     */
    template <class T>
    T GetValue(const std::string &attribute) const;


    /**
     * @brief Get the Object object
     * 
     * @param attribute 
     * @return nlohmann::json 
     */
    nlohmann::json GetObject(const std::string &attribute) const;


    /**
     * \brief Check if an attribute exists in the configuration file.
     * \param attribute The name of the requested attribute including the complete path hierarchy.
     * \return [true] The attribute exists in the configuration file.
     * \return [false] The attribute does not exist in the configuration file.
     */
    bool HasAttribute(const std::string &attribute) const;


    /**
     * @brief Check if the data of the attribute is a single value.
     * @param attribute The name of the requested attribute including the complete path hierarchy.
     * @return [true] The data of the attribute is a single value.
     * @return [false] The data of the attribute is not a single value.
     */
    bool IsSingleValue(const std::string& attribute) const;


    /**
     * \brief Check if the data of the attribute is an array.
     * \param attribute The name of the requested attribute including the complete path hierarchy.
     * \return [true] The data of the attribute is an array.
     * \return [false] The data of the attribute is not an array.
     */
    bool IsArray(const std::string& attribute) const;


    /**
     * \brief Check if the data of the attribute is an array of arrays.
     * \param attribute The name of the requested attribute including the complete path hierarchy.
     * \return [true] The data of the attribute is an array of arrays.
     * \return [false] The data of the attribute is not an array of arrays.
     */
    bool IsMultiArray(const std::string& attribute) const;


    /**
     * \brief Check and fix a filename's path if it does not match with the parsed parent path.
     * \param [in] filename The filename to check for its path.
     * \return [std::string] The filename with the resolved path to match the parsed parent path.
     */
    std::string ResolvePath(const std::string &filename) const;


    /**
     * \brief Get the parent path of the parsed simulation configuration file.
     * \return [const std::string&] The parent path of the parsed simulation configuration file.
     */
    inline const std::string & ParentPath() const { return this->parent_path_; }


};


} // end of namespace APPS
/** \} End of Doxygen Groups */

#endif //ELECTRA_APPS_TOOLS_PARSER_HPP_

#include "parser.tpp"