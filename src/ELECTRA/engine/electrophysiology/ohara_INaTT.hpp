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
   \file ohara_INaTT.hpp
   \brief Ohara_INaTT class header file.
   \date 03/07/2024
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_OHARA_INaTT_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_OHARA_INaTT_HPP_

#include "ELECTRA/engine/electrophysiology/ohara.hpp"
#include "ELECTRA/engine/utilities/algorithm.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <limits>


namespace ELECTRA {

/** \addtogroup Electrophysiology \{ */

/**
 * \class Ohara_INaTT
 * \author Konstantinos A. Mountris
 * \brief  O'hara Rudy '11, human cardiac ventricular action potential model \cite ohara2011.
 */

class Ohara_INaTT : public Ohara
{

public:

    /**
     * \brief The default constructor of the Ohara_INaTT class.
     */
    Ohara_INaTT();


    /**
     * \brief The default destructor of the Ohara_INaTT class.
     */
    virtual ~Ohara_INaTT();


    /**
     * \brief Initialize the variables and parameters of the Ohara_INaTT cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Ohara_INaTT cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);



};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_OHARA_INaTT_HPP_