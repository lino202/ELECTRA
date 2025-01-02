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


#include "ELECTRA/engine/electrophysiology/paci2020.hpp"


namespace ELECTRA {

void Paci2020::SetDataMapping()
{
    using namespace Paci2020Var;
    using namespace Paci2020Prm;
    using namespace Paci2020Cur;

    // Set variables mapping.
    this->mapped_data_["v"]     = static_cast<std::size_t>(v);
    this->mapped_data_["Nai"]   = static_cast<std::size_t>(Nai);
    this->mapped_data_["Cai"]   = static_cast<std::size_t>(Cai);
    this->mapped_data_["m"]     = static_cast<std::size_t>(m);
    this->mapped_data_["h"]     = static_cast<std::size_t>(h);
    this->mapped_data_["j"]     = static_cast<std::size_t>(j);
    this->mapped_data_["d"]     = static_cast<std::size_t>(d);
    this->mapped_data_["f1"]    = static_cast<std::size_t>(f1);
    this->mapped_data_["f2"]    = static_cast<std::size_t>(f2);
    this->mapped_data_["fCa"]   = static_cast<std::size_t>(fCa);
    this->mapped_data_["Xr1"]   = static_cast<std::size_t>(Xr1);
    this->mapped_data_["Xr2"]   = static_cast<std::size_t>(Xr2);
    this->mapped_data_["Xs"]    = static_cast<std::size_t>(Xs);
    this->mapped_data_["Xf"]    = static_cast<std::size_t>(Xf);
    this->mapped_data_["q"]     = static_cast<std::size_t>(q);
    this->mapped_data_["r"]     = static_cast<std::size_t>(r);
    this->mapped_data_["Ca_SR"] = static_cast<std::size_t>(Ca_SR);
    this->mapped_data_["m_L"]   = static_cast<std::size_t>(m_L);
    this->mapped_data_["h_L"]   = static_cast<std::size_t>(h_L);
    this->mapped_data_["RyRa"]  = static_cast<std::size_t>(RyRa);
    this->mapped_data_["RyRo"]  = static_cast<std::size_t>(RyRo);
    this->mapped_data_["RyRc"]  = static_cast<std::size_t>(RyRc);

    // Set parameters mapping.
    this->mapped_data_["Cm"]             = static_cast<std::size_t>(Cm);
    this->mapped_data_["R"]              = static_cast<std::size_t>(R);
    this->mapped_data_["T"]              = static_cast<std::size_t>(T);
    this->mapped_data_["F"]              = static_cast<std::size_t>(F);
    this->mapped_data_["Nao"]            = static_cast<std::size_t>(Nao);
    this->mapped_data_["Cao"]            = static_cast<std::size_t>(Cao);
    this->mapped_data_["Ko"]             = static_cast<std::size_t>(Ko);
    this->mapped_data_["Ki"]             = static_cast<std::size_t>(Ki);
    this->mapped_data_["E_K"]            = static_cast<std::size_t>(E_K);
    this->mapped_data_["PkNa"]           = static_cast<std::size_t>(PkNa);
    this->mapped_data_["g_Na"]           = static_cast<std::size_t>(g_Na);
    this->mapped_data_["Vc"]             = static_cast<std::size_t>(Vc);
    this->mapped_data_["V_SR"]           = static_cast<std::size_t>(V_SR);
    this->mapped_data_["myCoefTauM"]     = static_cast<std::size_t>(myCoefTauM);
    this->mapped_data_["tauINaL"]        = static_cast<std::size_t>(tauINaL);
    this->mapped_data_["GNaLmax"]        = static_cast<std::size_t>(GNaLmax);
    this->mapped_data_["Vh_hLate"]       = static_cast<std::size_t>(Vh_hLate);
    this->mapped_data_["g_f"]            = static_cast<std::size_t>(g_f);
    this->mapped_data_["fNa"]            = static_cast<std::size_t>(fNa);
    this->mapped_data_["fK"]             = static_cast<std::size_t>(fK);
    this->mapped_data_["g_CaL"]          = static_cast<std::size_t>(g_CaL);
    this->mapped_data_["tau_fCa"]        = static_cast<std::size_t>(tau_fCa);
    this->mapped_data_["g_Kr"]           = static_cast<std::size_t>(g_Kr);
    this->mapped_data_["L0"]             = static_cast<std::size_t>(L0);
    this->mapped_data_["Q"]              = static_cast<std::size_t>(Q);
    this->mapped_data_["g_Ks"]           = static_cast<std::size_t>(g_Ks);
    this->mapped_data_["g_K1"]           = static_cast<std::size_t>(g_K1);
    this->mapped_data_["g_b_Na"]         = static_cast<std::size_t>(g_b_Na);
    this->mapped_data_["g_b_Ca"]         = static_cast<std::size_t>(g_b_Ca);
    this->mapped_data_["Km_K"]           = static_cast<std::size_t>(Km_K);
    this->mapped_data_["Km_Na"]          = static_cast<std::size_t>(Km_Na);
    this->mapped_data_["PNaK"]           = static_cast<std::size_t>(PNaK);
    this->mapped_data_["kNaCa"]          = static_cast<std::size_t>(kNaCa);
    this->mapped_data_["alpha"]          = static_cast<std::size_t>(alpha);
    this->mapped_data_["local_gamma"]    = static_cast<std::size_t>(local_gamma);
    this->mapped_data_["Ksat"]           = static_cast<std::size_t>(Ksat);
    this->mapped_data_["KmCa"]           = static_cast<std::size_t>(KmCa);
    this->mapped_data_["KmNai"]          = static_cast<std::size_t>(KmNai);
    this->mapped_data_["g_PCa"]          = static_cast<std::size_t>(g_PCa);
    this->mapped_data_["KPCa"]           = static_cast<std::size_t>(KPCa);
    this->mapped_data_["g_to"]           = static_cast<std::size_t>(g_to);
    this->mapped_data_["Kup"]            = static_cast<std::size_t>(Kup);
    this->mapped_data_["Buf_C"]          = static_cast<std::size_t>(Buf_C);
    this->mapped_data_["Buf_SR"]         = static_cast<std::size_t>(Buf_SR);
    this->mapped_data_["Kbuf_C"]         = static_cast<std::size_t>(Kbuf_C);
    this->mapped_data_["Kbuf_SR"]        = static_cast<std::size_t>(Kbuf_SR);
    this->mapped_data_["VmaxUp"]         = static_cast<std::size_t>(VmaxUp);
    this->mapped_data_["V_leak"]         = static_cast<std::size_t>(V_leak);
    this->mapped_data_["V_half"]         = static_cast<std::size_t>(V_half);
    this->mapped_data_["g_irel_max"]     = static_cast<std::size_t>(g_irel_max);
    this->mapped_data_["RyRa1"]          = static_cast<std::size_t>(RyRa1);
    this->mapped_data_["RyRa2"]          = static_cast<std::size_t>(RyRa2);
    this->mapped_data_["RyRahalf"]       = static_cast<std::size_t>(RyRahalf);
    this->mapped_data_["RyRohalf"]       = static_cast<std::size_t>(RyRohalf);
    this->mapped_data_["RyRchalf"]       = static_cast<std::size_t>(RyRchalf);
    this->mapped_data_["RyRtauadapt"]    = static_cast<std::size_t>(RyRtauadapt);

    // Set currents mapping.
    this->mapped_data_["INa"]   = static_cast<std::size_t>(INa);
    this->mapped_data_["INaK"]  = static_cast<std::size_t>(INaK);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["IbNa"]  = static_cast<std::size_t>(IbNa);
    this->mapped_data_["ICaL"]  = static_cast<std::size_t>(ICaL);
    this->mapped_data_["IK1"]   = static_cast<std::size_t>(IK1);
    this->mapped_data_["If"]    = static_cast<std::size_t>(If);
    this->mapped_data_["IKr"]   = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"]   = static_cast<std::size_t>(IKs);
    this->mapped_data_["Ito"]   = static_cast<std::size_t>(Ito);
    this->mapped_data_["IpCa"]  = static_cast<std::size_t>(IpCa);
    this->mapped_data_["IbCa"]  = static_cast<std::size_t>(IbCa);
    this->mapped_data_["Irel"]  = static_cast<std::size_t>(Irel);
    this->mapped_data_["Iup"]   = static_cast<std::size_t>(Iup);
    this->mapped_data_["Ileak"] = static_cast<std::size_t>(Ileak);
    this->mapped_data_["INaL"]  = static_cast<std::size_t>(INaL);

}


Paci2020::Paci2020()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Paci2020;
    this->dt_stable_ = 0.02;
    this->var_.resize(23, 0.);
    this->prm_.resize(56, 0.);
    this->cur_.resize(17, 0.);
    #ifdef BLOCK_CELL_CURRS
        this->block_coeff_.resize(16, 0.);
    #endif

    // Set mapped data.
    this->SetDataMapping();
}


Paci2020::~Paci2020()
{}


void Paci2020::Initialize(CellType cell_type)
{
    using namespace Paci2020Var;
    using namespace Paci2020Prm;

    if (cell_type != CellType::ventricular) {
        std::string error_str = "Could not initialize Paci2020 ap model. Expected: CellType::ventricular";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(23, 0.);
    this->prm_.clear();           this->prm_.resize(56, 0.);
    this->cur_.clear();           this->cur_.resize(17, 0.);
    #ifdef BLOCK_CELL_CURRS
        this->block_coeff_.clear();   this->block_coeff_.resize(16, 0.);
    #endif

    // Set the model variables.

    // Original initial values, we need dt 0.001 ms for convergence which we use in ElectraCell
    // then you can apply the resulting prepacing states at the tissue sim and go with 0.01 dt which is the default for the tissue sims
    this->var_[v]     = -70.0;   //This is changed from the default to have mV as the other models
    this->var_[dvdt]  = 0.0;
    this->var_[Nai]   = 9.2;     
    this->var_[Cai]   = 0.0002;  
    this->var_[m]     = 0.0;      
    this->var_[h]     = 0.75;      
    this->var_[j]     = 0.75;      
    this->var_[d]     = 0.0;    
    this->var_[f1]    = 1.0;     
    this->var_[f2]    = 1.0;     
    this->var_[fCa]   = 1.0;    
    this->var_[Xr1]   = 0.0;  
    this->var_[Xr2]   = 1.0;    
    this->var_[Xs]    = 0.0;    
    this->var_[Xf]    = 0.1;     
    this->var_[q]     = 1.0;      
    this->var_[r]     = 0.0;    
    this->var_[Ca_SR] = 0.32;    
    this->var_[m_L]   = 0.0;
    this->var_[h_L]   = 0.75;
    this->var_[RyRa]  = 0.3;
    this->var_[RyRo]  = 0.9;
    this->var_[RyRc]  = 0.1;

    // Prepaced with no stimulation after 800s, after this stabilization convergence can be achieved with 0.01, 
    // ATTENTION, we do not use this as default to be in agreement with the original model and as with stim from this state the APD does not stabilize (grows linearly with time)
    // this->var_[v]     = -73.4525804324366;
    // this->var_[dvdt]  = 0.;
    // this->var_[Nai]   = 9.00417069274294;     
    // this->var_[Cai]   = 1.76731736262123e-05;  
    // this->var_[m]     = 0.0441025126145443;      
    // this->var_[h]     = 0.736080688718430;      
    // this->var_[j]     = 0.742052141479516;      
    // this->var_[d]     = 0.000101671059827320;    
    // this->var_[f1]    = 0.979415292087528;     
    // this->var_[f2]    = 0.999979398897691;     
    // this->var_[fCa]   = 0.998948743326411;    
    // this->var_[Xr1]   = 0.00894672801754690;  
    // this->var_[Xr2]   = 0.427809818209464;    
    // this->var_[Xs]    = 0.0341200913071620;    
    // this->var_[Xf]    = 0.192854555974207;     
    // this->var_[q]     = 0.829179626325527;      
    // this->var_[r]     = 0.00601271295586856;    
    // this->var_[Ca_SR] = 0.115107956531682;    
    // this->var_[m_L]   = 0.00297784296969314;
    // this->var_[h_L]   = 0.135864730293044;
    // this->var_[RyRa]  = 0.0297153296103769;
    // this->var_[RyRo]  = 0.000710450936345816;
    // this->var_[RyRc]  = 0.948119074119825;


    // Set the model parameters.            
    this->prm_[Cm]                 = 9.87109e-11;               
    this->prm_[R]                  = 8.314472;                   
    this->prm_[T]                  = 310.0;                        
    this->prm_[F]                  = 96485.3415; 
    this->prm_[Nao]                = 151.0;                      
    this->prm_[Cao]                = 1.8;  
    this->prm_[Ko]                 = 5.4; 
    this->prm_[Ki]                 = 150.0;                                      
    this->prm_[E_K]                = (this->prm_[R]*this->prm_[T]/this->prm_[F])*std::log(this->prm_[Ko]/this->prm_[Ki]);
    this->prm_[PkNa]               = 0.03;
    this->prm_[g_Na]               = 6447.1896;
    this->prm_[Vc]                 = 8800.0;                      
    this->prm_[V_SR]               = 583.73; 
    this->prm_[myCoefTauM]         = 1.0;
    this->prm_[tauINaL]            = 200.0;
    this->prm_[GNaLmax]            = 2.3*7.5;
    this->prm_[Vh_hLate]           = 87.61;
    this->prm_[g_f]                = 22.2763088;  
    this->prm_[fNa]                = 0.37;
    this->prm_[fK]                 = 1.0 - this->prm_[fNa];
    this->prm_[g_CaL]              = 8.635702e-5;            
    this->prm_[tau_fCa]            = 0.002;                
    this->prm_[g_Kr]               = 29.8667;                 
    this->prm_[L0]                 = 0.025;                     
    this->prm_[Q]                  = 2.3;                        
    this->prm_[g_Ks]               = 2.041;                   
    this->prm_[g_K1]               = 28.1492;                                  
    this->prm_[g_b_Na]             = 1.14;                   
    this->prm_[g_b_Ca]             = 0.8727264;               
    this->prm_[Km_K]               = 1.0;                       
    this->prm_[Km_Na]              = 40.0;                     
    this->prm_[PNaK]               = 2.74240;                
    this->prm_[kNaCa]              = 6514.47574;                   
    this->prm_[alpha]              = 2.16659;              
    this->prm_[local_gamma]        = 0.35;  //gamma is also a function so local_gamma                   
    this->prm_[Ksat]               = 0.1;                     
    this->prm_[KmCa]               = 1.38;                    
    this->prm_[KmNai]              = 87.5;                   
    this->prm_[g_PCa]              = 0.4125;                 
    this->prm_[KPCa]               = 0.0005;                  
    this->prm_[g_to]               = 29.9038;                                       
    this->prm_[Kup]                = 4.40435e-4;               
    this->prm_[Buf_C]              = 0.25;             
    this->prm_[Buf_SR]             = 10.0;                  
    this->prm_[Kbuf_C]             = 0.001;              
    this->prm_[Kbuf_SR]            = 0.3;             
    this->prm_[VmaxUp]             = 0.82205;               
    this->prm_[V_leak]             = 4.48209e-4;             
    this->prm_[V_half]             = 1000.0*(-this->prm_[R]*this->prm_[T]/(this->prm_[F]*this->prm_[Q])*std::log(std::pow(1.0+this->prm_[Cao]/2.6,4.0)/(this->prm_[L0]*std::pow(1.0+this->prm_[Cao]/0.58,4.0)))-0.019);                 
    this->prm_[g_irel_max]         = 55.808061;
    this->prm_[RyRa1]              = 0.05169;
    this->prm_[RyRa2]              = 0.050001;
    this->prm_[RyRahalf]           = 0.02632;
    this->prm_[RyRohalf]           = 0.00944;
    this->prm_[RyRchalf]           = 0.00167;
    this->prm_[RyRtauadapt]        = 1.0;
}


void Paci2020::Compute(double v_new, double dt, double stim_current)
{
    using namespace Paci2020Var;
    using namespace Paci2020Prm;
    using namespace Paci2020Cur;

    //IMPORTANT 
    // In the original model the stim is 5.5e-10 A (ampere) and 5 ms pulse amplitude and duration, respectively, for a single cell stim, 
    // This is valid when "reference units" field in the input json (or in electra_cell.cpp) are mA for current and ms for time and stim amp 5.5e-7 mA    
    // Convert to SI, Paci uses SI so V is in volts, Time is in seconds and Currents in Amperes
    dt           *= 1.e-3; // changed to seconds
    v_new        *= 1.e-3; // changed to volts
    stim_current *= 1.e-3; // changed to ampere

    // Nernst potential
    double E_Na = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log(this->prm_[Nao]/this->var_[Nai]);
    double E_Ca = 0.5*this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log(this->prm_[Cao]/this->var_[Cai]);
    double E_Ks = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log((this->prm_[Ko]+this->prm_[PkNa]*this->prm_[Nao])/(this->prm_[Ki]+this->prm_[PkNa]*this->var_[Nai]));

    // INa adapted from DOI:10.3389/fphys.2018.00080
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INa] =  (1.0-this->block_coeff_[INa]) * (this->prm_[g_Na]*std::pow(this->var_[m], 3.0)*this->var_[h]*this->var_[j]*(v_new - E_Na));
    #else
        this->cur_[INa] =  this->prm_[g_Na]*std::pow(this->var_[m], 3.0)*this->var_[h]*this->var_[j]*(v_new - E_Na);
    #endif

    double m_inf  = 1.0 / (1.0 + std::exp((v_new*1000.0 + 39.0)/-11.2));
    double tau_m  = (0.00001 + 0.00013*std::exp(-std::pow((v_new*1000.0 + 48.0)/15.0, 2.0)) + 0.000045 / (1.0 + std::exp((v_new*1000.0 + 42.0)/-5.0)));
    this->var_[m] = ALGORITHM::RushLarsen(m_inf, this->var_[m], dt, tau_m);

    double h_inf  = 1.0 / (1.0 + std::exp((v_new*1000.0 + 66.5)/6.8));
    double tau_h  = (0.00007 + 0.034 / (1.0 + std::exp((v_new*1000.0 + 41.0)/5.5) + std::exp(-(v_new*1000.0 + 41.0)/14.0)) + 0.0002 / (1.0 + std::exp(-(v_new*1000.0 + 79.0)/14.0)));
    this->var_[h] = ALGORITHM::RushLarsen(h_inf, this->var_[h], dt, tau_h);

    // double j_inf  = h_inf; //commented for saving memory as they are equal
    double tau_j  = 10.0*(0.0007 + 0.15 / (1.0 + std::exp((v_new*1000.0 + 41.0)/5.5) + std::exp(-(v_new*1000.0 + 41.0)/14.0)) + 0.002 / (1.0 + std::exp(-(v_new*1000.0 + 79.0)/14.0)));
    this->var_[j] = ALGORITHM::RushLarsen(h_inf, this->var_[j], dt, tau_j);

    // INaL
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INaL] = (1.0-this->block_coeff_[INaL]) * (this->prm_[GNaLmax]* std::pow(this->var_[m_L],3)*this->var_[h_L]*(v_new-E_Na));
    #else
        this->cur_[INaL] = this->prm_[GNaLmax]* std::pow(this->var_[m_L],3.0)*this->var_[h_L]*(v_new-E_Na);
    #endif

    double m_inf_L     = 1.0/(1.0+std::exp(-(v_new*1000.0+42.85)/(5.264)));
    double alpha_m_L   = 1.0/(1.0+std::exp((-60.0-v_new*1000.0)/5.0));
    double beta_m_L    = 0.1/(1.0+std::exp((v_new*1000.0+35.0)/5.0))+0.1/(1.0+std::exp((v_new*1000.0-50.0)/200.0));
    double tau_m_L     = 1.0/1000.0 * this->prm_[myCoefTauM]*alpha_m_L*beta_m_L;
    this->var_[m_L] = ALGORITHM::RushLarsen(m_inf_L, this->var_[m_L], dt, tau_m_L);

    double h_inf_L     = 1.0/(1.0+std::exp((v_new*1000.0+this->prm_[Vh_hLate])/(7.488)));
    double tau_h_L     = 1.0/1000.0 * this->prm_[tauINaL];
    this->var_[h_L] = ALGORITHM::RushLarsen(h_inf_L, this->var_[h_L], dt, tau_h_L);

    // If adapted from DOI:10.3389/fphys.2018.00080
    double i_fK  = this->prm_[fK]*this->prm_[g_f]*this->var_[Xf]*(v_new - this->prm_[E_K]);
    double i_fNa = this->prm_[fNa]*this->prm_[g_f]*this->var_[Xf]*(v_new - E_Na);
    
    #ifdef BLOCK_CELL_CURRS
        this->cur_[If] = (1.0-this->block_coeff_[If]) * (i_fK + i_fNa);
    #else
        this->cur_[If] = i_fK + i_fNa;
    #endif

    double Xf_infinity = 1.0/(1.0 + std::exp((v_new*1000.0 + 69.0)/8.0));
    double tau_Xf      = 5600.0 / (1.0 + std::exp((v_new*1000.0 + 65.0)/7.0) + std::exp(-(v_new*1000.0 + 65.0)/19.0));
    this->var_[Xf] = ALGORITHM::ForwardEuler(this->var_[Xf], dt, 1000.0*(Xf_infinity-this->var_[Xf])/tau_Xf);  //This is not rushlarsen

    // ICaL
    #ifdef BLOCK_CELL_CURRS
        this->cur_[ICaL]       = (1.0-this->block_coeff_[ICaL]) * (this->prm_[g_CaL]*4.0*v_new*std::pow(this->prm_[F],2.0)/(this->prm_[R]*this->prm_[T])*(this->var_[Cai]*std::exp(2.0*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))-0.341*this->prm_[Cao])/(std::exp(2.0*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))-1.0)*this->var_[d]*this->var_[f1]*this->var_[f2]*this->var_[fCa]);
    #else
        this->cur_[ICaL]       = this->prm_[g_CaL]*4.0*v_new*std::pow(this->prm_[F],2.0)/(this->prm_[R]*this->prm_[T])*(this->var_[Cai]*std::exp(2.0*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))-0.341*this->prm_[Cao])/(std::exp(2.0*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))-1.0)*this->var_[d]*this->var_[f1]*this->var_[f2]*this->var_[fCa];
    #endif
    

    double d_infinity  = 1.0/(1.0+std::exp(-(v_new*1000.0+9.1)/7.0));
    double alpha_d     = 0.25+1.4/(1.0+std::exp((-v_new*1000.0-35.0)/13.0));
    double beta_d      = 1.4/(1.0+std::exp((v_new*1000.0+5.0)/5.0));
    double gamma_d     = 1.0/(1.0+std::exp((-v_new*1000.0+50.0)/20.0));
    double tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
    this->var_[d] = ALGORITHM::RushLarsen(d_infinity, this->var_[d], dt, tau_d);

    double f1_inf  = 1.0/(1.0+std::exp((v_new*1000.0+26.0)/3.0));
    double constf1 = 1.0;
    if (f1_inf-this->var_[f1] > 0.0){
        constf1 = 1.0+1433.0*(this->var_[Cai]-50.0*1.0e-6);
    }
    double tau_f1      = (20.0+1102.5*std::exp(-std::pow(((v_new*1000.0+27.0)/15.0),2.0))+200.0/(1.0+std::exp((13.0-v_new*1000.0)/10.0))+180.0/(1.0+std::exp((30.0+v_new*1000.0)/10.0)))*constf1/1000.0;
    this->var_[f1] = ALGORITHM::RushLarsen(f1_inf, this->var_[f1], dt, tau_f1);

    double f2_inf      = 0.33+0.67/(1.0+std::exp((v_new*1000.0+32.0)/4.0));
    double constf2     = 1.0;
    double tau_f2      = (600.0*std::exp(-std::pow((v_new*1000.0+25.0),2.0)/170.0)+31.0/(1.0+std::exp((25.0-v_new*1000.0)/10.0))+16.0/(1.0+std::exp((30.0+v_new*1000.0)/10.0)))*constf2/1000.0;
    this->var_[f2] = ALGORITHM::RushLarsen(f2_inf, this->var_[f2], dt, tau_f2);

    double alpha_fCa   = 1.0/(1.0+std::pow(this->var_[Cai]/0.0006,8.0));
    double beta_fCa    = 0.1/(1.0+std::exp((this->var_[Cai]-0.0009)/0.0001));
    double gamma_fCa   = 0.3/(1.0+std::exp((this->var_[Cai]-0.00075)/0.0008));
    double fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
    double constfCa    = 1.0;
    if ((v_new > -0.06) && (fCa_inf > this->var_[fCa])){
        constfCa = 0.0;
    }

    this->var_[fCa] = ALGORITHM::ForwardEuler(this->var_[fCa], dt, constfCa*(fCa_inf-this->var_[fCa])/this->prm_[tau_fCa]); // This is no rushlarsen otherwise ICaL is zero for constfCa = 0

    // Ito
    #ifdef BLOCK_CELL_CURRS
        this->cur_[Ito]        = (1.0-this->block_coeff_[Ito]) * (this->prm_[g_to]*(v_new-this->prm_[E_K])*this->var_[q]*this->var_[r]);
    #else
        this->cur_[Ito]        = this->prm_[g_to]*(v_new-this->prm_[E_K])*this->var_[q]*this->var_[r];
    #endif
    
    double q_inf       = 1.0/(1.0+std::exp((v_new*1000.0+53.0)/13.0));
    double tau_q       = (6.06+39.102/(0.57*std::exp(-0.08*(v_new*1000.0+44.0))+0.065*std::exp(0.1*(v_new*1000.0+45.93))))/1000.0;
    this->var_[q] = ALGORITHM::RushLarsen(q_inf, this->var_[q], dt, tau_q);

    double r_inf       = 1.0/(1.0+std::exp(-(v_new*1000.0-22.3)/18.75));
    double tau_r       = (2.75352+14.40516/(1.037*std::exp(0.09*(v_new*1000.0+30.61))+0.369*std::exp(-0.12*(v_new*1000.0+23.84))))/1000.0;
    this->var_[r] = ALGORITHM::RushLarsen(r_inf, this->var_[r], dt, tau_r);
    
    // IKs
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IKs]        = (1.0-this->block_coeff_[IKs]) * (this->prm_[g_Ks]*(v_new-E_Ks)*std::pow(this->var_[Xs],2.0)*(1.0+0.6/(1.0+std::pow(3.8*0.00001/this->var_[Cai],1.4))));
    #else
        this->cur_[IKs]        = this->prm_[g_Ks]*(v_new-E_Ks)*std::pow(this->var_[Xs],2.0)*(1.0+0.6/(1.0+std::pow(3.8*0.00001/this->var_[Cai],1.4)));
    #endif
    
    double Xs_infinity = 1.0/(1.0+std::exp((-v_new*1000.0-20.0)/16.0));
    double alpha_Xs    = 1100.0/std::sqrt(1.0+std::exp((-10.0-v_new*1000.0)/6.0));
    double beta_Xs     = 1.0/(1.0+std::exp((-60.0+v_new*1000.0)/20.0));
    double tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
    this->var_[Xs] = ALGORITHM::RushLarsen(Xs_infinity, this->var_[Xs], dt, tau_Xs);

    // IKr
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IKr] = (1.0-this->block_coeff_[IKr]) * (this->prm_[g_Kr]*(v_new-this->prm_[E_K])*this->var_[Xr1]*this->var_[Xr2]*std::sqrt(this->prm_[Ko]/5.4));
    #else
        this->cur_[IKr] = this->prm_[g_Kr]*(v_new-this->prm_[E_K])*this->var_[Xr1]*this->var_[Xr2]*std::sqrt(this->prm_[Ko]/5.4);
    #endif
    
    double Xr1_inf      = 1.0/(1.0+std::exp((this->prm_[V_half]-v_new*1000.0)/4.9));
    double alpha_Xr1    = 450.0/(1.0+std::exp((-45.0-v_new*1000.0)/10.0));
    double beta_Xr1     = 6.0/(1.0+std::exp((30.0+v_new*1000.0)/11.5));
    double tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
    this->var_[Xr1] = ALGORITHM::RushLarsen(Xr1_inf, this->var_[Xr1], dt, tau_Xr1);

    double Xr2_infinity = 1.0/(1.0+std::exp((v_new*1000.0+88.0)/50.0));
    double alpha_Xr2    = 3.0/(1.0+std::exp((-60.0-v_new*1000.0)/20.0));
    double beta_Xr2     = 1.12/(1.0+std::exp((-60.0+v_new*1000.0)/20.0));
    double tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
    this->var_[Xr2] = ALGORITHM::RushLarsen(Xr2_infinity, this->var_[Xr2], dt, tau_Xr2);

    // IK1
    double alpha_K1    = 3.91/(1.0+std::exp(0.5942*(v_new*1000.0-this->prm_[E_K]*1000.0-200.0)));
    double beta_K1     = (-1.509*std::exp(0.0002*(v_new*1000.0-this->prm_[E_K]*1000.0+100.0))+std::exp(0.5886*(v_new*1000.0-this->prm_[E_K]*1000.0-10.0)))/(1.0+std::exp(0.4547*(v_new*1000.0-this->prm_[E_K]*1000.0)));
    double XK1_inf     = alpha_K1/(alpha_K1+beta_K1);

    // IKr
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IK1] = (1.0-this->block_coeff_[IK1]) * (this->prm_[g_K1]*XK1_inf*(v_new-this->prm_[E_K])*std::sqrt(this->prm_[Ko]/5.4));
    #else
        this->cur_[IK1] = this->prm_[g_K1]*XK1_inf*(v_new-this->prm_[E_K])*std::sqrt(this->prm_[Ko]/5.4);
    #endif

    // INaCa, INaK and IpCa
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INaCa] = (1.0-this->block_coeff_[INaCa]) * (this->prm_[kNaCa]*(std::exp(this->prm_[local_gamma]*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))*std::pow(this->var_[Nai],3.0)*this->prm_[Cao]-std::exp((this->prm_[local_gamma]-1.0)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))*std::pow(this->prm_[Nao],3.0)*this->var_[Cai]*this->prm_[alpha])/((std::pow(this->prm_[KmNai],3.0)+std::pow(this->prm_[Nao],3.0))*(this->prm_[KmCa]+this->prm_[Cao])*(1.0+this->prm_[Ksat]*std::exp((this->prm_[local_gamma]-1.0)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T])))));
        this->cur_[INaK]  = (1.0-this->block_coeff_[INaK]) * (this->prm_[PNaK]*this->prm_[Ko]/(this->prm_[Ko]+this->prm_[Km_K])*this->var_[Nai]/(this->var_[Nai]+this->prm_[Km_Na])/(1.0+0.1245*std::exp(-0.1*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))+0.0353*std::exp(-v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))));    
        this->cur_[IpCa]  = (1.0-this->block_coeff_[IpCa]) * (this->prm_[g_PCa]*this->var_[Cai]/(this->var_[Cai]+this->prm_[KPCa]));
        
    #else
        this->cur_[INaCa] = this->prm_[kNaCa]*(std::exp(this->prm_[local_gamma]*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))*std::pow(this->var_[Nai],3.0)*this->prm_[Cao]-std::exp((this->prm_[local_gamma]-1.0)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))*std::pow(this->prm_[Nao],3.0)*this->var_[Cai]*this->prm_[alpha])/((std::pow(this->prm_[KmNai],3.0)+std::pow(this->prm_[Nao],3.0))*(this->prm_[KmCa]+this->prm_[Cao])*(1.0+this->prm_[Ksat]*std::exp((this->prm_[local_gamma]-1.0)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))));
        this->cur_[INaK]  = this->prm_[PNaK]*this->prm_[Ko]/(this->prm_[Ko]+this->prm_[Km_K])*this->var_[Nai]/(this->var_[Nai]+this->prm_[Km_Na])/(1.0+0.1245*std::exp(-0.1*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))+0.0353*std::exp(-v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T])));    
        this->cur_[IpCa]  = this->prm_[g_PCa]*this->var_[Cai]/(this->var_[Cai]+this->prm_[KPCa]);
    #endif

    // Background currents and Sarcoplasmic reticulum
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IbNa]  = (1.0-this->block_coeff_[IbNa]) * (this->prm_[g_b_Na]*(v_new-E_Na));
        this->cur_[IbCa]  = (1.0-this->block_coeff_[IbCa]) * (this->prm_[g_b_Ca]*(v_new-E_Ca));
        this->cur_[Iup]   = (1.0-this->block_coeff_[Iup]) * (this->prm_[VmaxUp]/(1.0+std::pow(this->prm_[Kup],2.0)/std::pow(this->var_[Cai],2.0)));
        this->cur_[Ileak] = (1.0-this->block_coeff_[Ileak]) * ((this->var_[Ca_SR]-this->var_[Cai])*this->prm_[V_leak]);
    #else
        this->cur_[IbNa]  = this->prm_[g_b_Na]*(v_new-E_Na);
        this->cur_[IbCa]  = this->prm_[g_b_Ca]*(v_new-E_Ca);
        this->cur_[Iup]   = this->prm_[VmaxUp]/(1.0+std::pow(this->prm_[Kup],2.0)/std::pow(this->var_[Cai],2.0));
        this->cur_[Ileak] = (this->var_[Ca_SR]-this->var_[Cai])*this->prm_[V_leak];
    #endif


    // RyR
    double RyRSRCass   = (1.0 - 1.0/(1.0 +  std::exp((this->var_[Ca_SR]-0.3)/0.1)));
    #ifdef BLOCK_CELL_CURRS
        this->cur_[Irel] = (1.0-this->block_coeff_[Irel]) * (this->prm_[g_irel_max]*RyRSRCass*this->var_[RyRo]*this->var_[RyRc]*(this->var_[Ca_SR]-this->var_[Cai]));
    #else
        this->cur_[Irel] = this->prm_[g_irel_max]*RyRSRCass*this->var_[RyRo]*this->var_[RyRc]*(this->var_[Ca_SR]-this->var_[Cai]);
    #endif
    
    double RyRainfss   = this->prm_[RyRa1]-this->prm_[RyRa2]/(1.0 + std::exp((1000.0*this->var_[Cai]-(this->prm_[RyRahalf]))/0.0082));
    this->var_[RyRa] = ALGORITHM::RushLarsen(RyRainfss, this->var_[RyRa], dt, this->prm_[RyRtauadapt]);
    

    double RyRoinfss   = (1.0 - 1.0/(1.0 +  std::exp((1000.0*this->var_[Cai]-(this->var_[RyRa]+ this->prm_[RyRohalf]))/0.003)));
    double RyRtauact = 0.1*18.75e-3;
    if (RyRoinfss>= this->var_[RyRo]){
        RyRtauact = 18.75e-3;
    }
    this->var_[RyRo] = ALGORITHM::RushLarsen(RyRoinfss, this->var_[RyRo], dt, RyRtauact);

    double RyRcinfss   = (1.0/(1.0 + std::exp((1000.0*this->var_[Cai]-(this->var_[RyRa]+this->prm_[RyRchalf]))/0.001)));
    double RyRtauinact = 87.5e-3;
    if (RyRcinfss>= this->var_[RyRc]){
        RyRtauinact = 2*87.5e-3;
    }
    this->var_[RyRc] = ALGORITHM::RushLarsen(RyRcinfss, this->var_[RyRc], dt, RyRtauinact);
    
    // Ca2+ buffering
    double Cai_bufc    = 1.0/(1.0+this->prm_[Buf_C]*this->prm_[Kbuf_C]/std::pow(this->var_[Cai]+this->prm_[Kbuf_C],2.0));
    double Ca_SR_bufSR = 1.0/(1.0+this->prm_[Buf_SR]*this->prm_[Kbuf_SR]/std::pow(this->var_[Ca_SR]+this->prm_[Kbuf_SR],2.0));

    // Ionic concentrations
    double dNai = -this->prm_[Cm]*(this->cur_[INa]+this->cur_[INaL]+this->cur_[IbNa]+3.0*this->cur_[INaK]+3.0*this->cur_[INaCa]+i_fNa)/(this->prm_[F]*this->prm_[Vc]*1.0e-18);
    double dCai = Cai_bufc*(this->cur_[Ileak]-this->cur_[Iup]+this->cur_[Irel]-(this->cur_[ICaL]+this->cur_[IbCa]+this->cur_[IpCa]-2.0*this->cur_[INaCa])*this->prm_[Cm]/(2.0*this->prm_[Vc]*this->prm_[F]*1.0e-18));
    double dcaSR = Ca_SR_bufSR*this->prm_[Vc]/this->prm_[V_SR]*(this->cur_[Iup]-(this->cur_[Irel]+this->cur_[Ileak]));

    this->var_[Nai] = ALGORITHM::ForwardEuler(this->var_[Nai], dt, dNai);
    this->var_[Cai] = ALGORITHM::ForwardEuler(this->var_[Cai], dt, dCai);
    this->var_[Ca_SR] = ALGORITHM::ForwardEuler(this->var_[Ca_SR], dt, dcaSR);

    // Membrane potential    
    this->cur_[Paci2020Cur::Iion] = this->cur_[IK1] + this->cur_[Ito] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[ICaL] + 
                                  this->cur_[INaK] + this->cur_[INa] + this->cur_[INaCa] + this->cur_[IpCa] + this->cur_[If] + 
                                  this->cur_[IbNa] + this->cur_[IbCa] + this->cur_[INaL];

    this->var_[dvdt] =  - (this->cur_[Paci2020Cur::Iion] - stim_current/this->prm_[Cm]);

}


std::string Paci2020::PrintVariables() const 
{
    using namespace Paci2020Var;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << " v     = " << this->var_[v] << "\n";
    oss << " dvdt  = " << this->var_[dvdt] << "\n";
    oss << " Nai   = " << this->var_[Nai] << "\n";
    oss << " Cai   = " << this->var_[Cai] << "\n";
    oss << " m     = " << this->var_[m] << "\n";
    oss << " h     = " << this->var_[h] << "\n";
    oss << " j     = " << this->var_[j] << "\n";
    oss << " d     = " << this->var_[d] << "\n";
    oss << " f1    = " << this->var_[f1] << "\n";
    oss << " f2    = " << this->var_[f2] << "\n";
    oss << " fCa   = " << this->var_[fCa] << "\n";
    oss << " Xr1   = " << this->var_[Xr1] << "\n";
    oss << " Xr2   = " << this->var_[Xr2] << "\n";
    oss << " Xs    = " << this->var_[Xs] << "\n";
    oss << " Xf    = " << this->var_[Xf] << "\n";
    oss << " q     = " << this->var_[q] << "\n";
    oss << " r     = " << this->var_[r] << "\n";
    oss << " Ca_SR = " << this->var_[Ca_SR] << "\n";
    oss << " m_L   = " << this->var_[m_L] << "\n";
    oss << " h_L   = " << this->var_[h_L] << "\n";
    oss << " RyRa  = " << this->var_[RyRa] << "\n";
    oss << " RyRo  = " << this->var_[RyRo] << "\n";
    oss << " RyRc  = " << this->var_[RyRc];
    return oss.str();

}


std::string Paci2020::PrintParameters() const 
{
    using namespace Paci2020Prm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Cm          = " << this->prm_[Cm] << "\n";
    oss << "R           = " << this->prm_[R] << "\n";
    oss << "T           = " << this->prm_[T] << "\n";
    oss << "F           = " << this->prm_[F] << "\n";
    oss << "Nao         = " << this->prm_[Nao] << "\n";
    oss << "Cao         = " << this->prm_[Cao] << "\n";
    oss << "Ko          = " << this->prm_[Ko] << "\n";
    oss << "Ki          = " << this->prm_[Ki] << "\n";
    oss << "E_K         = " << this->prm_[E_K] << "\n";
    oss << "PkNa        = " << this->prm_[PkNa] << "\n";
    oss << "g_Na        = " << this->prm_[g_Na] << "\n";
    oss << "Vc          = " << this->prm_[Vc] << "\n";
    oss << "V_SR        = " << this->prm_[V_SR] << "\n";
    oss << "myCoefTauM  = " << this->prm_[myCoefTauM] << "\n";
    oss << "tauINaL     = " << this->prm_[tauINaL] << "\n";
    oss << "GNaLmax     = " << this->prm_[GNaLmax] << "\n";
    oss << "Vh_hLate    = " << this->prm_[Vh_hLate] << "\n";
    oss << "g_f         = " << this->prm_[g_f] << "\n";
    oss << "fNa         = " << this->prm_[fNa] << "\n";
    oss << "fK          = " << this->prm_[fK] << "\n";
    oss << "g_CaL       = " << this->prm_[g_CaL] << "\n";
    oss << "tau_fCa     = " << this->prm_[tau_fCa] << "\n";
    oss << "g_Kr        = " << this->prm_[g_Kr] << "\n";
    oss << "L0          = " << this->prm_[L0] << "\n";
    oss << "Q           = " << this->prm_[Q] << "\n";
    oss << "g_Ks        = " << this->prm_[g_Ks] << "\n";
    oss << "g_K1        = " << this->prm_[g_K1] << "\n";
    oss << "g_b_Na      = " << this->prm_[g_b_Na] << "\n";
    oss << "g_b_Ca      = " << this->prm_[g_b_Ca] << "\n";
    oss << "Km_K        = " << this->prm_[Km_K] << "\n";
    oss << "Km_Na       = " << this->prm_[Km_Na] << "\n";
    oss << "PNaK        = " << this->prm_[PNaK] << "\n";
    oss << "kNaCa       = " << this->prm_[kNaCa] << "\n";
    oss << "alpha       = " << this->prm_[alpha] << "\n";
    oss << "local_gamma = " << this->prm_[local_gamma] << "\n";
    oss << "Ksat        = " << this->prm_[Ksat] << "\n";
    oss << "KmCa        = " << this->prm_[KmCa] << "\n";
    oss << "KmNai       = " << this->prm_[KmNai] << "\n";
    oss << "g_PCa       = " << this->prm_[g_PCa] << "\n";
    oss << "KPCa        = " << this->prm_[KPCa] << "\n";
    oss << "g_to        = " << this->prm_[g_to] << "\n";
    oss << "Kup         = " << this->prm_[Kup] << "\n";
    oss << "Buf_C       = " << this->prm_[Buf_C] << "\n";
    oss << "Buf_SR      = " << this->prm_[Buf_SR] << "\n";
    oss << "Kbuf_C      = " << this->prm_[Kbuf_C] << "\n";
    oss << "Kbuf_SR     = " << this->prm_[Kbuf_SR] << "\n";
    oss << "VmaxUp      = " << this->prm_[VmaxUp] << "\n";
    oss << "V_leak      = " << this->prm_[V_leak] << "\n";
    oss << "V_half      = " << this->prm_[V_half] << "\n";
    oss << "g_irel_max  = " << this->prm_[g_irel_max] << "\n";
    oss << "RyRa1       = " << this->prm_[RyRa1] << "\n";
    oss << "RyRa2       = " << this->prm_[RyRa2] << "\n";
    oss << "RyRahalf    = " << this->prm_[RyRahalf] << "\n";
    oss << "RyRohalf    = " << this->prm_[RyRohalf] << "\n";
    oss << "RyRchalf    = " << this->prm_[RyRchalf] << "\n";
    oss << "RyRtauadapt = " << this->prm_[RyRtauadapt];
    return oss.str();

}


std::string Paci2020::PrintCurrents() const
{
    using namespace Paci2020Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa   = " << this->cur_[INa] << "\n";
    oss << "INaK  = " << this->cur_[INaK] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "IbNa  = " << this->cur_[IbNa] << "\n";
    oss << "ICaL  = " << this->cur_[ICaL] << "\n";
    oss << "IK1   = " << this->cur_[IK1] << "\n";
    oss << "If    = " << this->cur_[If] << "\n";
    oss << "IKr   = " << this->cur_[IKr] << "\n";
    oss << "IKs   = " << this->cur_[IKs] << "\n";
    oss << "Ito   = " << this->cur_[Ito] << "\n";
    oss << "IpCa  = " << this->cur_[IpCa] << "\n";
    oss << "IbCa  = " << this->cur_[IbCa] << "\n";
    oss << "Irel  = " << this->cur_[Irel] << "\n";
    oss << "Iup   = " << this->cur_[Iup] << "\n";
    oss << "Ileak = " << this->cur_[Ileak] << "\n";
    oss << "INaL  = " << this->cur_[INaL] << "\n";
    oss << "Iion  = " << this->cur_[Paci2020Cur::Iion];
    return oss.str();

}

#ifdef BLOCK_CELL_CURRS
    std::string Paci2020::PrintBlockCoeffs() const
    {
        using namespace Paci2020Cur;

        // Create output string stream to pass the currents and their values.
        std::ostringstream oss;
        oss.precision(15);
        oss << "INa   = " << this->block_coeff_[INa] << "\n";
        oss << "INaK  = " << this->block_coeff_[INaK] << "\n";
        oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
        oss << "IbNa  = " << this->block_coeff_[IbNa] << "\n";
        oss << "ICaL  = " << this->block_coeff_[ICaL] << "\n";
        oss << "IK1   = " << this->block_coeff_[IK1] << "\n";
        oss << "If    = " << this->block_coeff_[If] << "\n";
        oss << "IKr   = " << this->block_coeff_[IKr] << "\n";
        oss << "IKs   = " << this->block_coeff_[IKs] << "\n";
        oss << "Ito   = " << this->block_coeff_[Ito] << "\n";
        oss << "IpCa  = " << this->block_coeff_[IpCa] << "\n";
        oss << "IbCa  = " << this->block_coeff_[IbCa] << "\n";
        oss << "Irel  = " << this->block_coeff_[Irel] << "\n";
        oss << "Iup   = " << this->block_coeff_[Iup] << "\n";
        oss << "Ileak = " << this->block_coeff_[Ileak] << "\n";
        oss << "INaL  = " << this->block_coeff_[INaL];
        return oss.str();

    }
#endif


} // End of namespace ELECTRA
