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


#include "IMP/engine/utilities/algorithms.hpp"

namespace IMP {
namespace ALGORITHMS {


bool StringCompCaseInsensitive(const std::string &str1, const std::string &str2)
{
    if (str1.size() != str2.size())
        return false;

    return std::equal(std::begin(str1), std::end(str1), std::begin(str1), std::end(str1),
                      [](const char &a, const char &b) { return std::tolower(a) == std::tolower(b); });
}


bool ExistExactWord(const std::string &phrase, const std::string &word)
{
    // Guess it is found.
    bool found = true;

    // First character position.
    std::size_t cp_pos = phrase.find(word[0]);
    if(cp_pos == std::string::npos) { return false; }

    // Next characters position.
    std::size_t cn_pos = std::string::npos;
    for (std::size_t i = 1; i != word.size(); ++i) {

        // Search for the charecter starting from the next char than the previously found.
        cn_pos = phrase.find(word[i], cp_pos+1);

        // If not found or not next to the previous one, return false.
        if(cn_pos == std::string::npos || cn_pos != (cp_pos+1)) { found = false; break; }

        // Set this position as previous for next iteration.
        cp_pos = cn_pos;
    }
    return found;
}


bool ExistWord(const std::string &phrase, const std::string &word)
{
    std::string low_case_phrase = phrase;
    std::transform(std::begin(low_case_phrase), std::end(low_case_phrase),
        std::begin(low_case_phrase), ::tolower);

    std::string low_case_word = word;
    std::transform(std::begin(low_case_word), std::end(low_case_word),
        std::begin(low_case_word), ::tolower);

    return ExistExactWord(low_case_phrase, low_case_word);
}


double DegToRad(double angle) noexcept
{
    return (angle*PI)/180.;
}


double InKelvin(double temp_val, const std::string &temp_unit) noexcept
{
    if (temp_unit == "C")
        return temp_val + 273.15;
    return temp_val;
}


double InCelsius(double temp_val, const std::string &temp_unit) noexcept
{
    if (temp_unit == "K")
        return temp_val - 273.15;
    return temp_val;
}


double InNewton(double force_val, const std::string &force_unit) noexcept
{
    if (force_unit == "g")
        return force_val*0.00980665;
    return force_val;
}


double InPascal(double p_val, const std::string &p_unit) noexcept
{
    if (p_unit == "g")
        return p_val*98.0665;
    return p_val;
}

} //end of namespace ALGORITHMS
} //end of namespace IMP