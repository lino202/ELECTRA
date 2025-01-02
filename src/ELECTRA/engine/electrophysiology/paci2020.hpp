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
   \file paci2020.hpp
   \brief Paci2020 class header file.
   \date 20/12/2024
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_PACI2020_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_PACI2020_HPP_

#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"
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
 * \namespace ELECTRA::Paci2020Var
 * \brief Encaplulates enumeration to access Paci2020 stem cell derived  human ventricular action potential model variables.
 */
namespace Paci2020Var {
      enum { 
            v,           /**< Cell membrane action potential */
            dvdt,
            Nai,
            Cai,
            m,
            h,
            j,
            d,
            f1,
            f2,
            fCa,
            Xr1,
            Xr2,
            Xs,
            Xf,
            q,
            r,
            Ca_SR,
            m_L,
            h_L,
            RyRa,
            RyRo,
            RyRc
      };
} // End of namespace Paci2020Var


/**
 * \namespace ELECTRA::Paci2020Prm
 * \brief Encaplulates enumeration to access Paci2020 stem cell derived  human ventricular action potential model parameters. 
 */
namespace Paci2020Prm {
      enum { 
            Cm,
            R,
            T,
            F,
            Nao,
            Cao,
            Ko,
            Ki,
            E_K,
            PkNa,
            g_Na,
            Vc,
            V_SR,
            myCoefTauM,
            tauINaL,
            GNaLmax,
            Vh_hLate,
            g_f,
            fNa,
            fK,
            g_CaL,
            tau_fCa,
            g_Kr,
            L0,
            Q,
            g_Ks,
            g_K1,
            g_b_Na,
            g_b_Ca,
            Km_K,
            Km_Na,
            PNaK,
            kNaCa,
            alpha,
            local_gamma,
            Ksat,
            KmCa,
            KmNai,
            g_PCa,
            KPCa,
            g_to,
            Kup,
            Buf_C,
            Buf_SR,
            Kbuf_C,
            Kbuf_SR,
            VmaxUp,
            V_leak,
            V_half,
            g_irel_max,
            RyRa1,
            RyRa2,
            RyRahalf,
            RyRohalf,
            RyRchalf,
            RyRtauadapt
      };
} // End of namespace Paci2020Prm


/**
 * \namespace ELECTRA::Paci2020Cur
 * \brief Encaplulates enumeration to access Paci2020 stem cell derived human ventricular action potential model currents. 
 */
namespace Paci2020Cur {
      enum { 
            INa,
            INaK,
            INaCa,
            IbNa,
            ICaL,
            IK1,
            If,
            IKr,
            IKs,
            Ito,
            IpCa,
            IbCa,
            Irel,
            Iup,
            Ileak,
            INaL,
            Iion
      };
} // End of namespace Paci2020Cur


/**
 * \class Paci2020
 * \brief  Paci2020, stem cell-derived human cardiac ventricular action potential model.
 */

class Paci2020 : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the Paci2020 class.
     */
    Paci2020();


    /**
     * \brief The default destructor of the Paci2020 class.
     */
    virtual ~Paci2020();


    /**
     * \brief Initialize the variables and parameters of the Paci2020 cell model. 
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Paci2020 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[Paci2020Var::v] = v; }


    /**
     * \brief Print to std::string the cell's variables and their values.
     * \return [std::string] The cell's variables and their values.
    */
    virtual std::string PrintVariables() const;


    /**
     * \brief Print to std::string the cell's parameters and their values.
     * \return [std::string] The cell's parameters and their values.
    */
    virtual std::string PrintParameters() const;


    /**
     * \brief Print to std::string the cell's currents and their values.
     * \return [std::string] The cell's currents and their values.
    */
    virtual std::string PrintCurrents() const;


    /**
     * \brief Print to std::string the cell's currents' block coefficients and their values.
     * \return [std::string] The cell's currents' block coefficients and their values.
    */
      #ifdef BLOCK_CELL_CURRS
            virtual std::string PrintBlockCoeffs() const;
      #endif


    /**
     * \brief Get the membrane potential value of the cell model.
     * \return [double] The membrane potential value.
     */
    inline virtual double V() const { return this->var_[Paci2020Var::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[Paci2020Var::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[Paci2020Cur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_PACI2020_HPP_