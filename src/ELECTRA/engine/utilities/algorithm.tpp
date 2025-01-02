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


#ifndef ELECTRA_UTILITIES_ALGORITHM_TPP_
#define ELECTRA_UTILITIES_ALGORITHM_TPP_


#include "ELECTRA/engine/utilities/algorithm.hpp"

namespace ELECTRA {

namespace ALGORITHM {


template<typename TYPE>
TYPE RushLarsen(TYPE val_steady, TYPE val_old, TYPE dt, TYPE time_steady)
{
    return val_steady - (val_steady - val_old) * std::exp(-dt/time_steady);
}


template<typename TYPE>
TYPE ForwardEuler(TYPE func_old, TYPE dt, TYPE func_derivative)
{
    return func_old + dt*func_derivative;
}


template<typename TYPE>
TYPE KahanSum(std::initializer_list<TYPE> values) {

    // Initialize sum and error.
    TYPE sum = static_cast<TYPE>(0.);
    TYPE c = static_cast<TYPE>(0.);

    for (auto val : values) {
        TYPE y = val - c;
        TYPE t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}


} //end of namespace ALGORITHM


} //end of namespace ELECTRA

#endif //ELECTRA_UTILITIES_ALGORITHM_TPP_