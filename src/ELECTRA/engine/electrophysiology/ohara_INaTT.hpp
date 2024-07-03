/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019
 *
 * ALL RIGHTS RESERVED
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