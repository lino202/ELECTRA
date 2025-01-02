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

#ifndef ELECTRA_APPS_TOOLS_PARSER_TPP_
#define ELECTRA_APPS_TOOLS_PARSER_TPP_

#include "parser.hpp"

namespace APP_ELECTRA
{


template <class T>
T Parser::GetValue(const std::string& attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            std::string error_message = "Required attribute is missing from JSON file: " + attribute;
            throw std::runtime_error(ELECTRA::Logger::Error(error_message));
        }
    }

    T return_value;
    if (current.is_array() && current.size()==1) return_value = current[0];
    else return_value = current;
    return return_value;
}

} // end of namespace APP_ELECTRA


#endif //ELECTRA_APPS_TOOLS_PARSER_TPP_