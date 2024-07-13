/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file stewart.hpp
   \brief Stewart class header file.
   \author Konstantinos A. Mountris
   \date 13/03/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_STEWART_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_STEWART_HPP_

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
 * \namespace ELECTRA::StrtVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Stewart '09 human purkinje action potential model variables. 
 */
namespace StrtVar {
enum { 
      v,      
      dvdt,   
      d,      
      f2,     
      fCass,  
      f,      
      Ca_SR,  
      Ca_i,   
      Ca_ss,  
      R_prime,
      h,      
      j,      
      m,      
      y,      
      K_i,    
      Xr1,    
      Xr2,    
      Xs,     
      Na_i,   
      r,      
      s,   
      };

} // End of namespace StrtVar


/**
 * \namespace ELECTRA::StrtPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Stewart '09 human purkinje action potential model parameters. 
 */
namespace StrtPrm {
enum { 
      g_CaL,
      g_bca,
      Buf_c,
      Buf_sr,
      Buf_ss,
      Ca_o,
      EC,
      K_buf_c,
      K_buf_sr,
      K_buf_ss,
      K_up,
      V_leak,
      V_rel,
      V_sr,
      V_ss,
      V_xfer,
      Vmax_up,
      k1_prime,
      k2_prime,
      k3,
      k4,
      max_sr,
      min_sr,
      K_pCa,
      g_pCa,
      g_Na,
      g_f_K,
      g_f_Na,
      g_K1,
      Cm,
      F,
      R,
      T,
      V_c,
      K_o,
      g_pK,
      g_Kr,
      P_kna,
      g_Ks,
      g_bna,
      K_NaCa,
      K_sat,
      Km_Ca,
      Km_Nai,
      alpha,
      gamma,
      Na_o,
      K_mNa,
      K_mk,
      P_NaK,
      g_sus,
      g_to
      };
} // End of namespace StrtPrm


/**
 * \namespace ELECTRA::StrtCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Stewart '09 human purkinje action potential model currents. 
 */
namespace StrtCur {
enum { IfNa,
       IfK,
       If,
       Irel,
       Isus,
       INaK,
       INaCa,
       IpCa,
       IpK,
       Iup,
       Ileak,
       Ixfer,
       IK1,
       IKr,
       IKs,
       INa,
       IbNa,
       ICaL,
       IbCa,
       Ito,
       Iion
      };
} // End of namespace StrtCur


/**
 * \class Stewart
 * \author Konstantinos A. Mountris
 * \brief  Stewart '09, human purkinje action potential model \cite stewart2009.
 */

class Stewart : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the Stewart class.
     */
    Stewart();


    /**
     * \brief The default destructor of the Stewart class.
     */
    virtual ~Stewart();


    /**
     * \brief Initialize the variables and parameters of the cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[StrtVar::v] = v; }


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
    inline virtual double V() const { return this->var_[StrtVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[StrtVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[StrtCur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_STEWART_HPP_