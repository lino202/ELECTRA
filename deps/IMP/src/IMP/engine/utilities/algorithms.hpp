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
   \file algorithms.hpp
   \brief General algorithm collection header file.
   \author Konstantinos A. Mountris
   \date 10/07/2019
*/

#ifndef IMP_ENGINE_UTILITIES_ALGORITHMS_HPP_
#define IMP_ENGINE_UTILITIES_ALGORITHMS_HPP_

#include <cmath>
#include <initializer_list>
#include <type_traits>
#include <string>
#include <algorithm>
#include <unordered_map>


namespace IMP {

/**
 * \namespace ALGORITHMS
 * \brief Collects general algorithmic functions for mathematical manipulations and more.
 */
namespace ALGORITHMS {


    const auto PI = 3.14159265359;      /*< Pi number */


/** \addtogroup Utilities \{ */


/**
 * \fn Sign(TYPE x)
 * \brief Check the sign of an arithmetic value.
 * \tparam TYPE The data type of the arithmetic value.
 * \param [in] x The arithmetic value.
 * \return [int] 1 | 0 | -1 Depending on the sign of the given arithmetic value.
 */
template <typename TYPE>
inline int Sign(TYPE x) { return (x > static_cast<TYPE>(0)) - (x < static_cast<TYPE>(0)); }


/**
 * \brief Temporal integration of differential equation using the Rush-Larsen method \cite rush1978. 
 * \tparam TYPE The data type. A standard floating point type should be preferred. 
 * \param [in] val_steady The value of the differential equation at steady state. 
 * \param [in] val_old The value of the differential equation at the previous time step.
 * \param [in] dt The time increment between time steps.
 * \param [in] time_steady The final time at steady state.
 * \return [TYPE] The new value of the differential equation at the next time step.
 */
template<typename TYPE>
inline TYPE RushLarsen(TYPE val_steady, TYPE val_old, TYPE dt, TYPE time_steady);


/**
 * \brief Temporal integration of differential equation using the Forward Euler method.
 * \tparam TYPE The data type. A standard floating point type should be preferred.
 * \param [in] func_old The value of the differential equation at the previous time step.
 * \param [in] dt The time increment between time steps.
 * \param [in] func_derivative The value of the used function to update the value of the previous time step.
 * \return [TYPE] The new value of the differential equation at the next time step.
 */
template<typename TYPE>
inline TYPE ForwardEuler(TYPE func_old, TYPE dt, TYPE func_derivative);


/**
 * \fn KahanSum(std::initializer_list<TYPE> values)
 * \brief Kahan summation algorithm for precise summation \cite kahan1965.
 * \tparam TYPE The data type. A standard floating point type should be preferred.
 * \param [in] values The list of values to be summed.
 * \return [TYPE] The sum of the given values.
 */
template<typename TYPE>
inline TYPE KahanSum(std::initializer_list<TYPE> values);


template<typename TYPE>
inline TYPE ProdsDiff(TYPE a, TYPE b, TYPE c, TYPE d);


/**
 * \brief 
 * \param str1 
 * \param str2 
 * \return true 
 * \return false 
 */
bool StringCompCaseInsensitive(const std::string &str1, const std::string &str2);


/**
 * \brief 
 * \param phrase 
 * \param word 
 * \return true 
 * \return false 
 */
bool ExistExactWord(const std::string &phrase, const std::string &word);


/**
 * \fn ExistWord(std::string phrase, std::string word)
 * \brief Checks if a word exist in a string.
 * \param [in] phrase The string to be searched.
 * \param [in] word The word that we search in the phrase.
 * \return [TRUE] The word is found in the phrase.
 * \return [FALSE] The word is not found in the phrase.
 */
bool ExistWord(const std::string &phrase, const std::string &word);


/**
 * \brief Convert an angle from degrees to radians.
 * \param [in] angle The angle in degrees.
 * \return [double] The angle in radians.
 */
double DegToRad(double angle) noexcept;


/**
 * \brief Convert temperature unit in Kelvin.
 * \param [in] temp_val The temperature value.
 * \param [in] temp_unit The temperature unit.
 * \return [double] If the unit is C, then the value is converted in Kelvin.
 */
double InKelvin(double temp_val, const std::string &temp_unit) noexcept;


/**
 * \brief Convert temperature unit in Celsius.
 * \param [in] temp_val The temperature value.
 * \param [in] temp_unit The temperature unit.
 * \return [double] If the unit is K, then the value is converted in Celsius.
 */
double InCelsius(double temp_val, const std::string &temp_unit) noexcept;


/**
 * @brief 
 * 
 * @param force_val 
 * @param force_unit 
 * @return double 
 */
double InNewton(double force_val, const std::string &force_unit) noexcept;


double InPascal(double p_val, const std::string &p_unit) noexcept;


/**
 * \brief Add a value in a 2D unordered map at the position specified by indices i, j.
 * If a value already exist the new value is summed.
 * \tparam TYPE The type of the value.
 * \param [in] i The row index of the 2D unordered map.
 * \param [in] val The value to be inserted or summed at the (i,j) entry of the 2D unordered map with index.
 * \param [in] umap The unordered map were the value will be added.
 */
template<typename TYPE>
inline void AddInUmap1D(int i, TYPE val, std::unordered_map<int, TYPE> &umap);


/**
 * \brief Add a value in a 2D unordered map at the position specified by indices i, j.
 * If a value already exist the new value is summed.
 * \tparam TYPE The type of the value.
 * \param [in] i The row index of the 2D unordered map.
 * \param [in] j The column index of the 2D unordered map.
 * \param [in] val The value to be inserted or summed at the (i,j) entry of the 2D unordered map with index.
 * \param [in] umap The unordered map were the value will be added.
 */
template<typename TYPE>
inline void AddInUmap2D(int i, int j, TYPE val, std::unordered_map<int, std::unordered_map<int,TYPE>> &umap);


/** \} End of Doxygen Groups*/

} //end of namespace ALGORITHMS

} //end of namespace IMP

#endif //IMP_ENGINE_UTILITIES_ALGORITHMS_HPP_

#include "IMP/engine/utilities/algorithms.tpp"