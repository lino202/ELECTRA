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


#include "ELECTRA/engine/electrophysiology/stewart.hpp"


namespace ELECTRA {

void Stewart::SetDataMapping()
{
    using namespace StrtVar;
    using namespace StrtPrm;
    using namespace StrtCur;

    // Set variables mapping.
    this->mapped_data_["v"]        = static_cast<std::size_t>(v);     
    this->mapped_data_["d"]        = static_cast<std::size_t>(d);     
    this->mapped_data_["f2"]       = static_cast<std::size_t>(f2);    
    this->mapped_data_["fCass"]    = static_cast<std::size_t>(fCass); 
    this->mapped_data_["f"]        = static_cast<std::size_t>(f);     
    this->mapped_data_["Ca_SR"]    = static_cast<std::size_t>(Ca_SR); 
    this->mapped_data_["Ca_i"]     = static_cast<std::size_t>(Ca_i);  
    this->mapped_data_["Ca_ss"]    = static_cast<std::size_t>(Ca_ss); 
    this->mapped_data_["R_prime"]  = static_cast<std::size_t>(R_prime);
    this->mapped_data_["h"]        = static_cast<std::size_t>(h);     
    this->mapped_data_["j"]        = static_cast<std::size_t>(j);     
    this->mapped_data_["m"]        = static_cast<std::size_t>(m);     
    this->mapped_data_["y"]        = static_cast<std::size_t>(y);     
    this->mapped_data_["K_i"]      = static_cast<std::size_t>(K_i);   
    this->mapped_data_["Xr1"]      = static_cast<std::size_t>(Xr1);   
    this->mapped_data_["Xr2"]      = static_cast<std::size_t>(Xr2);   
    this->mapped_data_["Xs"]       = static_cast<std::size_t>(Xs);    
    this->mapped_data_["Na_i"]     = static_cast<std::size_t>(Na_i);  
    this->mapped_data_["r"]        = static_cast<std::size_t>(r);     
    this->mapped_data_["s"]        = static_cast<std::size_t>(s); 

    // Set parameters mapping.
    this->mapped_data_["g_CaL"]   = static_cast<std::size_t>(g_CaL);
    this->mapped_data_["g_bca"]   = static_cast<std::size_t>(g_bca);
    this->mapped_data_["Buf_c"]   = static_cast<std::size_t>(Buf_c);
    this->mapped_data_["Buf_sr"]   = static_cast<std::size_t>(Buf_sr);
    this->mapped_data_["Buf_ss"]   = static_cast<std::size_t>(Buf_ss);
    this->mapped_data_["Ca_o"]   = static_cast<std::size_t>(Ca_o);
    this->mapped_data_["EC"]   = static_cast<std::size_t>(EC);
    this->mapped_data_["K_buf_c"]   = static_cast<std::size_t>(K_buf_c);
    this->mapped_data_["K_buf_sr"]   = static_cast<std::size_t>(K_buf_sr);
    this->mapped_data_["K_buf_ss"]   = static_cast<std::size_t>(K_buf_ss);
    this->mapped_data_["K_up"]   = static_cast<std::size_t>(K_up);
    this->mapped_data_["V_leak"]   = static_cast<std::size_t>(V_leak);
    this->mapped_data_["V_rel"]   = static_cast<std::size_t>(V_rel);
    this->mapped_data_["V_sr"]   = static_cast<std::size_t>(V_sr);
    this->mapped_data_["V_ss"]   = static_cast<std::size_t>(V_ss);
    this->mapped_data_["V_xfer"]   = static_cast<std::size_t>(V_xfer);
    this->mapped_data_["Vmax_up"]   = static_cast<std::size_t>(Vmax_up);
    this->mapped_data_["k1_prime"]   = static_cast<std::size_t>(k1_prime);
    this->mapped_data_["k2_prime"]   = static_cast<std::size_t>(k2_prime);
    this->mapped_data_["k3"]   = static_cast<std::size_t>(k3);
    this->mapped_data_["k4"]   = static_cast<std::size_t>(k4);
    this->mapped_data_["max_sr"]   = static_cast<std::size_t>(max_sr);
    this->mapped_data_["min_sr"]   = static_cast<std::size_t>(min_sr);
    this->mapped_data_["K_pCa"]   = static_cast<std::size_t>(K_pCa);
    this->mapped_data_["g_pCa"]   = static_cast<std::size_t>(g_pCa);
    this->mapped_data_["g_Na"]   = static_cast<std::size_t>(g_Na);
    this->mapped_data_["g_f_K"]   = static_cast<std::size_t>(g_f_K);
    this->mapped_data_["g_f_Na"]   = static_cast<std::size_t>(g_f_Na);
    this->mapped_data_["g_K1"]   = static_cast<std::size_t>(g_K1);
    this->mapped_data_["Cm"]   = static_cast<std::size_t>(Cm);
    this->mapped_data_["F"]   = static_cast<std::size_t>(F);
    this->mapped_data_["R"]   = static_cast<std::size_t>(R);
    this->mapped_data_["T"]   = static_cast<std::size_t>(T);
    this->mapped_data_["V_c"]   = static_cast<std::size_t>(V_c);
    this->mapped_data_["K_o"]   = static_cast<std::size_t>(K_o);
    this->mapped_data_["g_pK"]   = static_cast<std::size_t>(g_pK);
    this->mapped_data_["g_Kr"]   = static_cast<std::size_t>(g_Kr);
    this->mapped_data_["P_kna"]   = static_cast<std::size_t>(P_kna);
    this->mapped_data_["g_Ks"]   = static_cast<std::size_t>(g_Ks);
    this->mapped_data_["g_bna"]   = static_cast<std::size_t>(g_bna);
    this->mapped_data_["K_NaCa"]   = static_cast<std::size_t>(K_NaCa);
    this->mapped_data_["K_sat"]   = static_cast<std::size_t>(K_sat);
    this->mapped_data_["Km_Ca"]   = static_cast<std::size_t>(Km_Ca);
    this->mapped_data_["Km_Nai"]   = static_cast<std::size_t>(Km_Nai);
    this->mapped_data_["alpha"]   = static_cast<std::size_t>(alpha);
    this->mapped_data_["gamma"]   = static_cast<std::size_t>(gamma);
    this->mapped_data_["Na_o"]   = static_cast<std::size_t>(Na_o);
    this->mapped_data_["K_mNa"]   = static_cast<std::size_t>(K_mNa);
    this->mapped_data_["K_mk"]   = static_cast<std::size_t>(K_mk);
    this->mapped_data_["P_NaK"]   = static_cast<std::size_t>(P_NaK);
    this->mapped_data_["g_sus"]   = static_cast<std::size_t>(g_sus);
    this->mapped_data_["g_to"]   = static_cast<std::size_t>(g_to);

    // Set currents mapping.
    this->mapped_data_["INa"]   = static_cast<std::size_t>(INa);
    this->mapped_data_["IfNa"]  = static_cast<std::size_t>(IfNa);
    this->mapped_data_["IfK"]   = static_cast<std::size_t>(IfK);
    this->mapped_data_["If"]    = static_cast<std::size_t>(If);
    this->mapped_data_["Irel"]  = static_cast<std::size_t>(Irel);
    this->mapped_data_["Isus"]  = static_cast<std::size_t>(Isus);
    this->mapped_data_["INaK"]  = static_cast<std::size_t>(INaK);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["IpCa"]  = static_cast<std::size_t>(IpCa);
    this->mapped_data_["IpK"]   = static_cast<std::size_t>(IpK);
    this->mapped_data_["Iup"]   = static_cast<std::size_t>(Iup);
    this->mapped_data_["Ileak"] = static_cast<std::size_t>(Ileak);
    this->mapped_data_["Ixfer"] = static_cast<std::size_t>(Ixfer);
    this->mapped_data_["IK1"]   = static_cast<std::size_t>(IK1);
    this->mapped_data_["IKr"]   = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"]   = static_cast<std::size_t>(IKs);
    this->mapped_data_["INa"]   = static_cast<std::size_t>(INa);
    this->mapped_data_["IbNa"]  = static_cast<std::size_t>(IbNa);
    this->mapped_data_["ICaL"]  = static_cast<std::size_t>(ICaL);
    this->mapped_data_["IbCa"]  = static_cast<std::size_t>(IbCa);
    this->mapped_data_["Ito"]   = static_cast<std::size_t>(Ito);

}


Stewart::Stewart()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Stewart;
    this->dt_stable_ = 0.02;
    this->var_.resize(69, 0.);
    this->prm_.resize(60, 0.);
    this->cur_.resize(21, 0.);
    #ifdef BLOCK_CELL_CURRS
        this->block_coeff_.resize(20, 0.);
    #endif

    // Set mapped data.
    this->SetDataMapping();
}


Stewart::~Stewart()
{}


void Stewart::Initialize(CellType cell_type)
{
    using namespace StrtVar;
    using namespace StrtPrm;

    if (cell_type != CellType::purkinje) {
        std::string error_str = "Could not initialize the Stewart purkinje ap model. Expected cell type: ELECTRA::CellType::purkinje";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(69, 0.);
    this->prm_.clear();           this->prm_.resize(60, 0.);
    this->cur_.clear();           this->cur_.resize(21, 0.);
    #ifdef BLOCK_CELL_CURRS
        this->block_coeff_.clear();   this->block_coeff_.resize(20, 0.);
    #endif

    // Set the model variables.
    this->var_[v]         = -69.1370441635924;
    this->var_[dvdt]      = 0.0;
    this->var_[d]         = 0.000287906256206415;
    this->var_[f2]        = 0.995474890442185;
    this->var_[fCass]     = 0.999955429598213;
    this->var_[f]         = 0.989328560287987;
    this->var_[Ca_SR]     = 3.10836886659417;
    this->var_[Ca_i]      = 0.000101878186157052;
    this->var_[Ca_ss]     = 0.000446818714055411;
    this->var_[R_prime]   = 0.991580051907845;
    this->var_[h]         = 0.190678733735145;
    this->var_[j]         = 0.238219836154029;
    this->var_[m]         = 0.0417391656294997;
    this->var_[y]         = 0.0457562667986602;
    this->var_[K_i]       = 136.781894160227;
    this->var_[Xr1]       = 0.00550281999719088;
    this->var_[Xr2]       = 0.313213286437995;
    this->var_[Xs]        = 0.00953708522974789;
    this->var_[Na_i]      = 8.80420286531673;
    this->var_[r]         = 0.00103618091196912;
    this->var_[s]         = 0.96386101799501;

 
    // Set the model parameters.
    this->prm_[g_CaL] = 3.98e-5;   
    this->prm_[g_bca] = 0.000592;  
    this->prm_[Buf_c] = 0.2;       
    this->prm_[Buf_sr] = 10.0;     
    this->prm_[Buf_ss] = 0.4;      
    this->prm_[Ca_o] = 2.0;        
    this->prm_[EC] = 1.5;          
    this->prm_[K_buf_c] = 0.001;   
    this->prm_[K_buf_sr] = 0.3;    
    this->prm_[K_buf_ss] = 0.00025;
    this->prm_[K_up] = 0.00025;    
    this->prm_[V_leak] = 0.00036;  
    this->prm_[V_rel] = 0.102;     
    this->prm_[V_sr] = 0.001094;   
    this->prm_[V_ss] = 5.468e-5;   
    this->prm_[V_xfer] = 0.0038;   
    this->prm_[Vmax_up] = 0.006375;
    this->prm_[k1_prime] = 0.15;   
    this->prm_[k2_prime] = 0.045;  
    this->prm_[k3] = 0.06;         
    this->prm_[k4] = 0.005;        
    this->prm_[max_sr] = 2.5;      
    this->prm_[min_sr] = 1.0;      
    this->prm_[K_pCa] = 0.0005;    
    this->prm_[g_pCa] = 0.1238;    
    this->prm_[g_Na] = 130.5744;   
    this->prm_[g_f_K] = 0.0234346; 
    this->prm_[g_f_Na] = 0.0145654;
    this->prm_[g_K1] = 0.065;      
    this->prm_[Cm] = 0.185;        
    this->prm_[F] = 96485.3415;    
    this->prm_[R] = 8314.472;      
    this->prm_[T] = 310.0;         
    this->prm_[V_c] = 0.016404;    
    this->prm_[K_o] = 5.4;         
    this->prm_[g_pK] = 0.0146;     
    this->prm_[g_Kr] = 0.0918;     
    this->prm_[P_kna] = 0.03;      
    this->prm_[g_Ks] = 0.2352;     
    this->prm_[g_bna] = 0.00029;   
    this->prm_[K_NaCa] = 1000.0;   
    this->prm_[K_sat] = 0.1;       
    this->prm_[Km_Ca] = 1.38;      
    this->prm_[Km_Nai] = 87.5;     
    this->prm_[alpha] = 2.5;       
    this->prm_[gamma] = 0.35;      
    this->prm_[Na_o] = 140.0;      
    this->prm_[K_mNa] = 40.0;      
    this->prm_[K_mk] = 1.0;        
    this->prm_[P_NaK] = 2.724;     
    this->prm_[g_sus] = 0.0227;    
    this->prm_[g_to] = 0.08184;  

}


void Stewart::Compute(double v_new, double dt, double stim_current)
{
    using namespace StrtVar;
    using namespace StrtPrm;
    using namespace StrtCur;

    this->cur_[ICaL] = this->prm_[g_CaL]*this->var_[d]*this->var_[f]*this->var_[f2]*this->var_[fCass]*4.0*((v_new-15.0)*std::pow(this->prm_[F],2.0)/(this->prm_[R]*this->prm_[T]))*(0.25*this->var_[Ca_ss]*std::exp(2.0*(v_new-15.0)*this->prm_[F]/(this->prm_[R]*this->prm_[T]))-this->prm_[Ca_o])/(std::exp(2.0*(v_new-15.0)*this->prm_[F]/(this->prm_[R]*this->prm_[T]))-1.0);
    
	double d_inf = 1.0/(1.0+std::exp((-8.0-v_new)/7.5));
    double alpha_d = 1.4/(1.0+std::exp((-35.0-v_new)/13.0))+0.25;
    
	double beta_d = 1.4/(1.0+std::exp((v_new+5.0)/5.0));

    double gamma_d = 1.0/(1.0+std::exp((50.0-v_new)/20.0));
    double tau_d = 1.0*alpha_d*beta_d+gamma_d;
	
    double f2_inf = 0.67/(1.0+std::exp((v_new+35.0)/7.0))+0.33;
    double tau_f2 = 562.0*std::exp(-std::pow(v_new+27.0,2.0)/240.0)+31.0/(1.0+std::exp((25.0-v_new)/10.0))+80.0/(1.0+std::exp((v_new+30.0)/10.0));
	
	double fCass_inf = 0.6/(1.0+std::pow(this->var_[Ca_ss]/0.05,2.0))+0.4;
    double tau_fCass = 80.0/(1.0+std::pow(this->var_[Ca_ss]/0.05,2.0))+2.0;
    
	double f_inf = 1.0/(1.0+std::exp((v_new+20.0)/7.0));
    double tau_f = 1102.5*std::exp(-std::pow(v_new+27.0,2.0)/225.0)+200.0/(1.0+std::exp((13.0-v_new)/10.0))+180.0/(1.0+std::exp((v_new+30.0)/10.0))+20.0;
    
    
	double E_Ca = 0.5*this->prm_[R]*this->prm_[T]/this->prm_[F]*log(this->prm_[Ca_o]/this->var_[Ca_i]);
    this->cur_[IbCa] = this->prm_[g_bca]*(v_new-E_Ca);
    
	double kcasr = this->prm_[max_sr]-(this->prm_[max_sr]-this->prm_[min_sr])/(1.0+std::pow(this->prm_[EC]/this->var_[Ca_SR],2.0));
    
	double k1 = this->prm_[k1_prime]/kcasr;
    double O = k1*std::pow(this->var_[Ca_ss],2.0)*this->var_[R_prime]/(this->prm_[k3]+k1*std::pow(this->var_[Ca_ss],2.0));
    this->cur_[Irel] = this->prm_[V_rel]*O*(this->var_[Ca_SR]-this->var_[Ca_ss]);
    this->cur_[Iup] = this->prm_[Vmax_up]/(1.0+std::pow(this->prm_[K_up],2.0)/std::pow(this->var_[Ca_i],2.0));
    this->cur_[Ileak] = this->prm_[V_leak]*(this->var_[Ca_SR]-this->var_[Ca_i]);
    this->cur_[Ixfer] = this->prm_[V_xfer]*(this->var_[Ca_ss]-this->var_[Ca_i]);
    
	double k2 = this->prm_[k2_prime]*kcasr;
    
	double Ca_i_bufc   = 1.0/(1.0+this->prm_[Buf_c]*this->prm_[K_buf_c]/std::pow(this->var_[Ca_i]+this->prm_[K_buf_c],2.0));
    double Ca_sr_bufsr = 1.0/(1.0+this->prm_[Buf_sr]*this->prm_[K_buf_sr]/std::pow(this->var_[Ca_SR]+this->prm_[K_buf_sr],2.0));
    double Ca_ss_bufss = 1.0/(1.0+this->prm_[Buf_ss]*this->prm_[K_buf_ss]/std::pow(this->var_[Ca_ss]+this->prm_[K_buf_ss],2.0));
    this->cur_[IpCa] = this->prm_[g_pCa]*this->var_[Ca_i]/(this->var_[Ca_i]+this->prm_[K_pCa]);
    this->cur_[INaCa] = this->prm_[K_NaCa]*(std::exp(this->prm_[gamma]*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))*std::pow(this->var_[Na_i],3.0)*this->prm_[Ca_o]-std::exp((this->prm_[gamma]-1.0)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))*std::pow(this->prm_[Na_o],3.0)*this->var_[Ca_i]*this->prm_[alpha])/((std::pow(this->prm_[Km_Nai],3.0)+std::pow(this->prm_[Na_o],3.0))*(this->prm_[Km_Ca]+this->prm_[Ca_o])*(1.0+this->prm_[K_sat]*std::exp((this->prm_[gamma]-1.0)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))));
    
	double E_Na = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log(this->prm_[Na_o]/this->var_[Na_i]);
    this->cur_[INa] = this->prm_[g_Na]*std::pow(this->var_[m],3.0)*this->var_[h]*this->var_[j]*(v_new-E_Na);
    
	double h_inf = 1.0/std::pow((1.0+std::exp((v_new+71.55)/7.43)),2.0);
	double alpha_h = 0.0;
	double beta_h = 0.77/(0.13*(1.0+std::exp((v_new+10.66)/-11.1)));
    if (v_new < -40.0){
       alpha_h = 0.057*std::exp(-(v_new+80.0)/6.8);
	   beta_h = 2.7*std::exp(0.079*v_new)+310000.0*std::exp(0.3485*v_new);
    }
    double tau_h = 1.0/(alpha_h+beta_h);
    
	
	double j_inf = 1.0/std::pow((1.0+std::exp((v_new+71.55)/7.43)),2.0);
	double alpha_j = 0.0;
	double beta_j = 0.6*std::exp(0.057*v_new)/(1.0+std::exp(-0.1*(v_new+32.0)));
    if (v_new < -40.0){
       alpha_j = (-25428.0*std::exp(0.2444*v_new)-6.948e-6*std::exp(-0.04391*v_new))*(v_new+37.78)/(1.0+std::exp(0.311*(v_new+79.23)));
	   beta_j = 0.02424*std::exp(-0.01052*v_new)/(1.0+std::exp(-0.1378*(v_new+40.14)));
	}
    double tau_j = 1.0/(alpha_j+beta_j);
    
    
	double m_inf = 1.0/std::pow((1.0+std::exp((-56.86-v_new)/9.03)),2.0);
    double alpha_m = 1.0/(1.0+std::exp((-60.0-v_new)/5.0));
    double beta_m = 0.1/(1.0+std::exp((v_new+35.0)/5.0))+0.1/(1.0+std::exp((v_new-50.0)/200.0));
    double tau_m = 1.0*alpha_m*beta_m;
    this->cur_[IfNa] = this->var_[y]*this->prm_[g_f_Na]*(v_new-E_Na);
    
	double E_K = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log(this->prm_[K_o]/this->var_[K_i]);
    this->cur_[IfK] = this->var_[y]*this->prm_[g_f_K]*(v_new-E_K);
    this->cur_[If] = this->cur_[IfNa]+this->cur_[IfK];
    
	double y_inf = 1.0/(1.0+std::exp((v_new+80.6)/6.8));
    double alpha_y = 1.0*std::exp(-2.9-0.04*v_new);
    double beta_y = 1.0*std::exp(3.6+0.11*v_new);
    double tau_y = 4000.0/(alpha_y+beta_y);
    
	double xK1_inf = 1.0/(1.0+std::exp(0.1*(v_new+75.44)));
    this->cur_[IK1] = this->prm_[g_K1]*xK1_inf*(v_new-8.0-E_K);
    this->cur_[Ito] = this->prm_[g_to]*this->var_[r]*this->var_[s]*(v_new-E_K);
    
	double a = 1.0/(1.0+std::exp((5.0-v_new)/17.0));
    this->cur_[Isus] = this->prm_[g_sus]*a*(v_new-E_K);
    this->cur_[IKr] = this->prm_[g_Kr]*std::pow(this->prm_[K_o]/5.4,0.5)*this->var_[Xr1]*this->var_[Xr2]*(v_new-E_K);
    
	double E_Ks = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log((this->prm_[K_o]+this->prm_[P_kna]*this->prm_[Na_o])/(this->var_[K_i]+this->prm_[P_kna]*this->var_[Na_i]));
    this->cur_[IKs] = this->prm_[g_Ks]*std::pow(this->var_[Xs],2.0)*(v_new-E_Ks);
    this->cur_[INaK] = this->prm_[P_NaK]*this->prm_[K_o]/(this->prm_[K_o]+this->prm_[K_mk])*this->var_[Na_i]/(this->var_[Na_i]+this->prm_[K_mNa])/(1.0+0.1245*std::exp(-0.1*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]))+0.0353*std::exp(-v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T])));
    this->cur_[IbNa] = this->prm_[g_bna]*(v_new-E_Na);
    this->cur_[IpK] = this->prm_[g_pK]*(v_new-E_K)/(1.0+std::exp((25.0-v_new)/5.98));
    
    this->cur_[StrtCur::Iion] = this->cur_[IK1]+this->cur_[Ito]+this->cur_[Isus]+this->cur_[IKr]+this->cur_[IKs]+this->cur_[ICaL]+this->cur_[INaK]+this->cur_[INa]+this->cur_[IbNa]+this->cur_[INaCa]+this->cur_[IbCa]+this->cur_[IpK]+this->cur_[IpCa]+this->cur_[If];
    this->var_[dvdt] =  - (this->cur_[StrtCur::Iion] - stim_current);
    
    double xr1_inf = 1.0/(1.0+std::exp((-26.0-v_new)/7.0));
    double alpha_xr1 = 450.0/(1.0+std::exp((-45.0-v_new)/10.0));
    double beta_xr1 = 6.0/(1.0+std::exp((v_new+30.0)/11.5));
    double tau_xr1 = 1.0*alpha_xr1*beta_xr1;
    
	double xr2_inf = 1.0/(1.0+std::exp((v_new+88.0)/24.0));
    double alpha_xr2 = 3.0/(1.0+std::exp((-60.0-v_new)/20.0));
    double beta_xr2 = 1.12/(1.0+std::exp((v_new-60.0)/20.0));
    double tau_xr2 = 1.0*alpha_xr2*beta_xr2;
    
	double xs_inf = 1.0/(1.0+std::exp((-5.0-v_new)/14.0));
    double alpha_xs = 1400.0/std::pow((1.0+std::exp((5.0-v_new)/6.0)),0.5);
    double beta_xs = 1.0/(1.0+std::exp((v_new-35.0)/15.0));
    double tau_xs = 1.0*alpha_xs*beta_xs+80.0;
    
	double r_inf = 1.0/(1.0+std::exp((20.0-v_new)/13.0));
    double tau_r = 10.45*std::exp(-std::pow(v_new+40.0,2.0)/1800.0)+7.3;
    
	double s_inf = 1.0/(1.0+std::exp((v_new+27.0)/13.0));
    double tau_s = 85.0*std::exp(-std::pow(v_new+25.0,2.0)/320.0)+5.0/(1.0+std::exp((v_new-40.0)/5.0))+42.0;
    
	
	// Update gating variables 
	this->var_[d]  = ALGORITHM::RushLarsen(d_inf, this->var_[d], dt, tau_d);
	this->var_[f2] = ALGORITHM::RushLarsen(f2_inf, this->var_[f2], dt, tau_f2);
    this->var_[fCass] = ALGORITHM::RushLarsen(fCass_inf, this->var_[fCass], dt, tau_fCass);
    this->var_[f] = ALGORITHM::RushLarsen(f_inf, this->var_[f], dt, tau_f);
    this->var_[h] = ALGORITHM::RushLarsen(h_inf, this->var_[h], dt, tau_h);
    this->var_[j] = ALGORITHM::RushLarsen(j_inf, this->var_[j], dt, tau_j);
    this->var_[m] = ALGORITHM::RushLarsen(m_inf, this->var_[m], dt, tau_m);
    this->var_[y] = ALGORITHM::RushLarsen(y_inf, this->var_[y], dt, tau_y);
    this->var_[Xr1] = ALGORITHM::RushLarsen(xr1_inf, this->var_[Xr1], dt, tau_xr1);
    this->var_[Xr2] = ALGORITHM::RushLarsen(xr2_inf, this->var_[Xr2], dt, tau_xr2);
    this->var_[Xs] = ALGORITHM::RushLarsen(xs_inf, this->var_[Xs], dt, tau_xs);
    this->var_[r] = ALGORITHM::RushLarsen(r_inf, this->var_[r], dt, tau_r);
    this->var_[s] = ALGORITHM::RushLarsen(s_inf, this->var_[s], dt, tau_s);
	
	//Update non-gating variables
    double dR_prime = -k2*this->var_[Ca_ss]*this->var_[R_prime]+this->prm_[k4]*(1.0-this->var_[R_prime]);
    this->var_[R_prime] = ALGORITHM::ForwardEuler(this->var_[R_prime], dt, dR_prime);

	double dCa_i = Ca_i_bufc*((this->cur_[Ileak]-this->cur_[Iup])*this->prm_[V_sr]/this->prm_[V_c]+this->cur_[Ixfer]-1.0*(this->cur_[IbCa]+this->cur_[IpCa]-2.0*this->cur_[INaCa])*this->prm_[Cm]/(2.0*1.0*this->prm_[V_c]*this->prm_[F]));
    this->var_[Ca_i] = ALGORITHM::ForwardEuler(this->var_[Ca_i], dt, dCa_i);

    double dCa_SR = Ca_sr_bufsr*(this->cur_[Iup]-(this->cur_[Irel]+this->cur_[Ileak]));
    this->var_[Ca_SR] = ALGORITHM::ForwardEuler(this->var_[Ca_SR], dt, dCa_SR);

    double dCa_ss = Ca_ss_bufss*(-1.0*this->cur_[ICaL]*this->prm_[Cm]/(2.0*1.0*this->prm_[V_ss]*this->prm_[F])+this->cur_[Irel]*this->prm_[V_sr]/this->prm_[V_ss]-this->cur_[Ixfer]*this->prm_[V_c]/this->prm_[V_ss]);
    this->var_[Ca_ss] = ALGORITHM::ForwardEuler(this->var_[Ca_ss], dt, dCa_ss);

    double dK_i = -1.0*(this->cur_[IK1]+this->cur_[Ito]+this->cur_[IfK]+this->cur_[Isus]+this->cur_[IKr]+this->cur_[IKs]+this->cur_[IpK]-2.0*this->cur_[INaK])/(1.0*this->prm_[V_c]*this->prm_[F])*this->prm_[Cm];
    this->var_[K_i] = ALGORITHM::ForwardEuler(this->var_[K_i], dt, dK_i);

    double dNa_i = -1.0*(this->cur_[INa]+this->cur_[IbNa]+this->cur_[IfNa]+3.0*this->cur_[INaK]+3.0*this->cur_[INaCa])/(1.0*this->prm_[V_c]*this->prm_[F])*this->prm_[Cm];
    this->var_[Na_i] = ALGORITHM::ForwardEuler(this->var_[Na_i], dt, dNa_i);
}


std::string Stewart::PrintVariables() const 
{
    using namespace StrtVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v       = " << this->var_[v] << "\n";
    oss << "dvdt    = " << this->var_[dvdt] << "\n";
    oss << "d       = " << this->var_[d] << "\n";
    oss << "f2      = " << this->var_[f2] << "\n";
    oss << "fCass   = " << this->var_[fCass] << "\n";
    oss << "f       = " << this->var_[f] << "\n";
    oss << "Ca_SR   = " << this->var_[Ca_SR] << "\n";
    oss << "Ca_i    = " << this->var_[Ca_i] << "\n";
    oss << "Ca_ss   = " << this->var_[Ca_ss] << "\n";
    oss << "R_prime = " << this->var_[R_prime] << "\n";
    oss << "h       = " << this->var_[h] << "\n";
    oss << "j       = " << this->var_[j] << "\n";
    oss << "m       = " << this->var_[m] << "\n";
    oss << "y       = " << this->var_[y] << "\n";
    oss << "K_i     = " << this->var_[K_i] << "\n";
    oss << "Xr1     = " << this->var_[Xr1] << "\n";
    oss << "Xr2     = " << this->var_[Xr2] << "\n";
    oss << "Xs      = " << this->var_[Xs] << "\n";
    oss << "Na_i    = " << this->var_[Na_i] << "\n";
    oss << "r       = " << this->var_[r] << "\n";
    oss << "s       = " << this->var_[s] << "\n";
    return oss.str();

}


std::string Stewart::PrintParameters() const 
{
    using namespace StrtPrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "g_CaL     = " << this->prm_[g_CaL] << "\n";
    oss << "g_bca     = " << this->prm_[g_bca] << "\n";
    oss << "Buf_c     = " << this->prm_[Buf_c] << "\n";
    oss << "Buf_sr    = " << this->prm_[Buf_sr] << "\n";
    oss << "Buf_ss    = " << this->prm_[Buf_ss] << "\n";
    oss << "Ca_o      = " << this->prm_[Ca_o] << "\n";
    oss << "EC        = " << this->prm_[EC] << "\n";
    oss << "K_buf_c   = " << this->prm_[K_buf_c] << "\n";
    oss << "K_buf_sr  = " << this->prm_[K_buf_sr] << "\n";
    oss << "K_buf_ss  = " << this->prm_[K_buf_ss] << "\n";
    oss << "K_up      = " << this->prm_[K_up] << "\n";
    oss << "V_leak    = " << this->prm_[V_leak] << "\n";
    oss << "V_rel     = " << this->prm_[V_rel] << "\n";
    oss << "V_sr      = " << this->prm_[V_sr] << "\n";
    oss << "V_ss      = " << this->prm_[V_ss] << "\n";
    oss << "V_xfer    = " << this->prm_[V_xfer] << "\n";
    oss << "Vmax_up   = " << this->prm_[Vmax_up] << "\n";
    oss << "k1_prime  = " << this->prm_[k1_prime] << "\n";
    oss << "k2_prime  = " << this->prm_[k2_prime] << "\n";
    oss << "k3        = " << this->prm_[k3] << "\n";
    oss << "k4        = " << this->prm_[k4] << "\n";
    oss << "max_sr    = " << this->prm_[max_sr] << "\n";
    oss << "min_sr    = " << this->prm_[min_sr] << "\n";
    oss << "K_pCa     = " << this->prm_[K_pCa] << "\n";
    oss << "g_pCa     = " << this->prm_[g_pCa] << "\n";
    oss << "g_Na      = " << this->prm_[g_Na] << "\n";
    oss << "g_f_K     = " << this->prm_[g_f_K] << "\n";
    oss << "g_f_Na    = " << this->prm_[g_f_Na] << "\n";
    oss << "g_K1      = " << this->prm_[g_K1] << "\n";
    oss << "Cm        = " << this->prm_[Cm] << "\n";
    oss << "F         = " << this->prm_[F] << "\n";
    oss << "R         = " << this->prm_[R] << "\n";
    oss << "T         = " << this->prm_[T] << "\n";
    oss << "V_c       = " << this->prm_[V_c] << "\n";
    oss << "K_o       = " << this->prm_[K_o] << "\n";
    oss << "g_pK      = " << this->prm_[g_pK] << "\n";
    oss << "g_Kr      = " << this->prm_[g_Kr] << "\n";
    oss << "P_kna     = " << this->prm_[P_kna] << "\n";
    oss << "g_Ks      = " << this->prm_[g_Ks] << "\n";
    oss << "g_bna     = " << this->prm_[g_bna] << "\n";
    oss << "K_NaCa    = " << this->prm_[K_NaCa] << "\n";
    oss << "K_sat     = " << this->prm_[K_sat] << "\n";
    oss << "Km_Ca     = " << this->prm_[Km_Ca] << "\n";
    oss << "Km_Nai    = " << this->prm_[Km_Nai] << "\n";
    oss << "alpha     = " << this->prm_[alpha] << "\n";
    oss << "gamma     = " << this->prm_[gamma] << "\n";
    oss << "Na_o      = " << this->prm_[Na_o] << "\n";
    oss << "K_mNa     = " << this->prm_[K_mNa] << "\n";
    oss << "K_mk      = " << this->prm_[K_mk] << "\n";
    oss << "P_NaK     = " << this->prm_[P_NaK] << "\n";
    oss << "g_sus     = " << this->prm_[g_sus] << "\n";
    oss << "g_to      = " << this->prm_[g_to] << "\n";
    return oss.str();

}


std::string Stewart::PrintCurrents() const
{
    using namespace StrtCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "IfNa = " << this->cur_[IfNa] << "\n";
    oss << "IfK = " << this->cur_[IfK] << "\n";
    oss << "If = " << this->cur_[If] << "\n";
    oss << "Irel = " << this->cur_[Irel] << "\n";
    oss << "Isus = " << this->cur_[Isus] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "IpCa = " << this->cur_[IpCa] << "\n";
    oss << "IpK = " << this->cur_[IpK] << "\n";
    oss << "Iup = " << this->cur_[Iup] << "\n";
    oss << "Ileak = " << this->cur_[Ileak] << "\n";
    oss << "Ixfer = " << this->cur_[Ixfer] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "IbNa = " << this->cur_[IbNa] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "IbCa = " << this->cur_[IbCa] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "Iion = " << this->cur_[StrtCur::Iion];
    return oss.str();

}

#ifdef BLOCK_CELL_CURRS
    std::string Stewart::PrintBlockCoeffs() const
    {
        using namespace StrtCur;

        // Create output string stream to pass the currents and their values.
        std::ostringstream oss;
        oss.precision(15);
        oss << "IfNa = " << this->block_coeff_[IfNa] << "\n";
        oss << "IfK = " << this->block_coeff_[IfK] << "\n";
        oss << "If = " << this->block_coeff_[If] << "\n";
        oss << "Irel = " << this->block_coeff_[Irel] << "\n";
        oss << "Isus = " << this->block_coeff_[Isus] << "\n";
        oss << "INaK = " << this->block_coeff_[INaK] << "\n";
        oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
        oss << "IpCa = " << this->block_coeff_[IpCa] << "\n";
        oss << "IpK = " << this->block_coeff_[IpK] << "\n";
        oss << "Iup = " << this->block_coeff_[Iup] << "\n";
        oss << "Ileak = " << this->block_coeff_[Ileak] << "\n";
        oss << "Ixfer = " << this->block_coeff_[Ixfer] << "\n";
        oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
        oss << "IKr = " << this->block_coeff_[IKr] << "\n";
        oss << "IKs = " << this->block_coeff_[IKs] << "\n";
        oss << "INa = " << this->block_coeff_[INa] << "\n";
        oss << "IbNa = " << this->block_coeff_[IbNa] << "\n";
        oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
        oss << "IbCa = " << this->block_coeff_[IbCa] << "\n";
        oss << "Ito = " << this->block_coeff_[Ito];
        return oss.str();

    }
#endif


} // End of namespace ELECTRA
