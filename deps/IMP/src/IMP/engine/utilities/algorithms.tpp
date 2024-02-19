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


#ifndef IMP_ENGINE_UTILITIES_ALGORITHMS_TPP_
#define IMP_ENGINE_UTILITIES_ALGORITHMS_TPP_


#include "IMP/engine/utilities/algorithms.hpp"

namespace IMP {
namespace ALGORITHMS {


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


template<typename TYPE>
TYPE ProdsDiff(TYPE a, TYPE b, TYPE c, TYPE d)
{
    TYPE cd = c * d;
    TYPE err = std::fma(-c, d, cd);
    TYPE dop = std::fma(a, b, -cd);
    return dop + err;
};


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
void AddInUmap1D(int i, TYPE val, std::unordered_map<int, TYPE> &umap)
{
    // Check if the entry exist.
    auto it = umap.find(i);
    if (it == umap.end()) {
        // Entry is missing and we insert the value.
        umap[i] = val;
    } else {
        // Entry exists and we sum the new value.
        it->second += val;
    }
}


template<typename TYPE>
void AddInUmap2D(int i, int j, TYPE val, std::unordered_map<int, std::unordered_map<int,TYPE>> &umap)
{
    // Check if the row of the entry exist.
    auto row = umap.find(i);
    if (row == umap.end()) {
        // Row is missing and we insert the value.
        umap[i][j] = val;
    } else {
        // Row exist but should check if entry exists too.
        auto it = row->second.find(j);
        if (it == row->second.end()) {
            // Entry is missing and we insert it.
            row->second[j] = val;
        } else {
            // Entry exists and we sum the new value.
            it->second += val;
        }
    }
}


} //end of namespace ALGORITHM
} //end of namespace IMP

#endif //IMP_ENGINE_UTILITIES_ALGORITHMS_TPP_