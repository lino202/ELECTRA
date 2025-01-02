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
   \file load_curve.hpp
   \brief Header file for a load curve describing time-varying loading conditions.
   \author Konstantinos A. Mountris
   \date 17/03/2020
*/

#ifndef ELECTRA_CONDITIONS_LOAD_CURVE_HPP_
#define ELECTRA_CONDITIONS_LOAD_CURVE_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <IMP/IMP>

#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>


namespace ELECTRA {

/** \addtogroup Conditions \{ */

/**
 * \class LoadCurve
 * \author Konstantinos A. Mountris
 * \brief A load curve to describe time-varying loading conditions.
 */
class LoadCurve {

private:

    std::vector<double> data_;      /**< The load curve data */

public:

    /**
     * \brief The default constructor of the LoadCurve.
     */
    LoadCurve() : data_() {}


    /**
     * \brief The destructor of the LoadCurve.
     */
    virtual ~LoadCurve() {}


    inline void SetData(const std::vector<double> &data) { this->data_ = data; }


    inline void AddData(double data_entry) { this->data_.emplace_back(data_entry); }


    inline const std::vector<double> & Data() const { return this->data_; }


    inline double Data(std::size_t id) const { return this->data_[id]; }


    inline double DataAt(std::size_t id) const { return this->data_.at(id); }

};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_CONDITIONS_LOAD_CURVE_HPP_