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

#include "ELECTRA/engine/electrophysiology/ohara.hpp"


namespace ELECTRA {

void Ohara::SetDataMapping()
{
    using namespace OhrVar;
    using namespace OhrPrm;
    using namespace OhrCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["nai"] = static_cast<std::size_t>(nai);
    this->mapped_data_["nass"] = static_cast<std::size_t>(nass);
    this->mapped_data_["ki"] = static_cast<std::size_t>(ki);
    this->mapped_data_["kss"] = static_cast<std::size_t>(kss);
    this->mapped_data_["cai"] = static_cast<std::size_t>(cai);
    this->mapped_data_["cass"] = static_cast<std::size_t>(cass);
    this->mapped_data_["cansr"] = static_cast<std::size_t>(cansr);
    this->mapped_data_["cajsr"] = static_cast<std::size_t>(cajsr);
    this->mapped_data_["m"] = static_cast<std::size_t>(m);
    this->mapped_data_["hf"] = static_cast<std::size_t>(hf);
    this->mapped_data_["hs"] = static_cast<std::size_t>(hs);
    this->mapped_data_["j"] = static_cast<std::size_t>(j);
    this->mapped_data_["hsp"] = static_cast<std::size_t>(hsp);
    this->mapped_data_["jp"] = static_cast<std::size_t>(jp);
    this->mapped_data_["mL"] = static_cast<std::size_t>(mL);
    this->mapped_data_["hL"] = static_cast<std::size_t>(hL);
    this->mapped_data_["hLp"] = static_cast<std::size_t>(hLp);
    this->mapped_data_["a"] = static_cast<std::size_t>(a);
    this->mapped_data_["iF"] = static_cast<std::size_t>(iF);
    this->mapped_data_["iS"] = static_cast<std::size_t>(iS);
    this->mapped_data_["ap"] = static_cast<std::size_t>(ap);
    this->mapped_data_["iFp"] = static_cast<std::size_t>(iFp);
    this->mapped_data_["iSp"] = static_cast<std::size_t>(iSp);
    this->mapped_data_["d"] = static_cast<std::size_t>(d);
    this->mapped_data_["ff"] = static_cast<std::size_t>(ff);
    this->mapped_data_["fs"] = static_cast<std::size_t>(fs);
    this->mapped_data_["fcaf"] = static_cast<std::size_t>(fcaf);
    this->mapped_data_["fcas"] = static_cast<std::size_t>(fcas);
    this->mapped_data_["jca"] = static_cast<std::size_t>(jca);
    this->mapped_data_["nca"] = static_cast<std::size_t>(nca);
    this->mapped_data_["ffp"] = static_cast<std::size_t>(ffp);
    this->mapped_data_["fcafp"] = static_cast<std::size_t>(fcafp);
    this->mapped_data_["xrf"] = static_cast<std::size_t>(xrf);
    this->mapped_data_["xrs"] = static_cast<std::size_t>(xrs);
    this->mapped_data_["xs1"] = static_cast<std::size_t>(xs1);
    this->mapped_data_["xs2"] = static_cast<std::size_t>(xs2);
    this->mapped_data_["xk1"] = static_cast<std::size_t>(xk1);
    this->mapped_data_["Jrelnp"] = static_cast<std::size_t>(Jrelnp);
    this->mapped_data_["Jrelp"] = static_cast<std::size_t>(Jrelp);
    this->mapped_data_["CaMKt"] = static_cast<std::size_t>(CaMKt);

    // Set parameters mapping.
    this->mapped_data_["gNaL"] = static_cast<std::size_t>(gNaL);
    this->mapped_data_["delta_factor"] = static_cast<std::size_t>(delta_factor);
    this->mapped_data_["gto"] = static_cast<std::size_t>(gto);
    this->mapped_data_["pCa"] = static_cast<std::size_t>(pCa);
    this->mapped_data_["gKr"] = static_cast<std::size_t>(gKr);
    this->mapped_data_["gKs"] = static_cast<std::size_t>(gKs);
    this->mapped_data_["gK1"] = static_cast<std::size_t>(gK1);
    this->mapped_data_["gncx"] = static_cast<std::size_t>(gncx);
    this->mapped_data_["pNaK"] = static_cast<std::size_t>(pNaK);
    this->mapped_data_["gKb"] = static_cast<std::size_t>(gKb);
    this->mapped_data_["Jrel_0"] = static_cast<std::size_t>(Jrel_0);
    this->mapped_data_["Jrelp_0"] = static_cast<std::size_t>(Jrelp_0);
    this->mapped_data_["Jupnp_0"] = static_cast<std::size_t>(Jupnp_0);
    this->mapped_data_["Jupp_0"] = static_cast<std::size_t>(Jupp_0);
    this->mapped_data_["nao"] = static_cast<std::size_t>(nao);
    this->mapped_data_["cao"] = static_cast<std::size_t>(cao);
    this->mapped_data_["ko"] = static_cast<std::size_t>(ko);
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["Fdy"] = static_cast<std::size_t>(Fdy);
    this->mapped_data_["l"] = static_cast<std::size_t>(l);
    this->mapped_data_["rad"] = static_cast<std::size_t>(rad);
    this->mapped_data_["vcell"] = static_cast<std::size_t>(vcell);
    this->mapped_data_["ageo"] = static_cast<std::size_t>(ageo);
    this->mapped_data_["acap"] = static_cast<std::size_t>(acap);
    this->mapped_data_["vmyo"] = static_cast<std::size_t>(vmyo);
    this->mapped_data_["vnsr"] = static_cast<std::size_t>(vnsr);
    this->mapped_data_["vjsr"] = static_cast<std::size_t>(vjsr);
    this->mapped_data_["vss"] = static_cast<std::size_t>(vss);
    this->mapped_data_["KmCaMK"] = static_cast<std::size_t>(KmCaMK);
    this->mapped_data_["aCaMK"] = static_cast<std::size_t>(aCaMK);
    this->mapped_data_["bCaMK"] = static_cast<std::size_t>(bCaMK);
    this->mapped_data_["CaMKo"] = static_cast<std::size_t>(CaMKo);
    this->mapped_data_["KmCaM"] = static_cast<std::size_t>(KmCaM);
    this->mapped_data_["BSRmax"] = static_cast<std::size_t>(BSRmax);
    this->mapped_data_["KmBSR"] = static_cast<std::size_t>(KmBSR);
    this->mapped_data_["BSLmax"] = static_cast<std::size_t>(BSLmax);
    this->mapped_data_["KmBSL"] = static_cast<std::size_t>(KmBSL);
    this->mapped_data_["cmdnmax"] = static_cast<std::size_t>(cmdnmax);
    this->mapped_data_["Kmcmdn"] = static_cast<std::size_t>(Kmcmdn);
    this->mapped_data_["trpnmax"] = static_cast<std::size_t>(trpnmax);
    this->mapped_data_["Kmtrpn"] = static_cast<std::size_t>(Kmtrpn);
    this->mapped_data_["Csqnmax"] = static_cast<std::size_t>(Csqnmax);
    this->mapped_data_["kmcsqn"] = static_cast<std::size_t>(kmcsqn);
    this->mapped_data_["gNa"] = static_cast<std::size_t>(gNa);

    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["INaL"] = static_cast<std::size_t>(INaL);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["ICaNa"] = static_cast<std::size_t>(ICaNa);
    this->mapped_data_["ICaK"] = static_cast<std::size_t>(ICaK);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["INaCa_i"] = static_cast<std::size_t>(INaCa_i);
    this->mapped_data_["INaCa_ss"] = static_cast<std::size_t>(INaCa_ss);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["IKb"] = static_cast<std::size_t>(IKb);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    this->mapped_data_["IpCa"] = static_cast<std::size_t>(IpCa);

}


Ohara::Ohara()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::OHara;
    this->dt_stable_ = 0.02;
    this->var_.resize(42, 0.);
    this->prm_.resize(45, 0.);
    this->cur_.resize(18, 0.);
    #ifdef BLOCK_CELL_CURRS
        this->block_coeff_.resize(17, 0.);
    #endif
    

    // Set mapped data.
    this->SetDataMapping();
}


Ohara::~Ohara()
{}


void Ohara::Initialize(CellType cell_type)
{
    using namespace OhrVar;
    using namespace OhrPrm;

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(42, 0.);
    this->prm_.clear();           this->prm_.resize(45, 0.);
    this->cur_.clear();           this->cur_.resize(18, 0.);
    #ifdef BLOCK_CELL_CURRS
        this->block_coeff_.clear();   this->block_coeff_.resize(17, 0.);
    #endif

    // Set data according to the cell type.
    switch (cell_type)
    {
        case CellType::endo :         
            
            // Set the model parameters for endocardium cell type.
            this->prm_[gNaL] = 0.0075;
            this->prm_[delta_factor] = 0.;
            this->prm_[gto] = 0.02;
            this->prm_[pCa] = 0.0001;
            this->prm_[gKr] = 0.046;
            this->prm_[gKs] = 0.0034;
            this->prm_[gK1] = 0.1908;
            this->prm_[gncx] = 0.0008;
            this->prm_[pNaK] = 30.;
            this->prm_[gKb]  = 0.003;
            this->prm_[Jrel_0] = 1.;
            this->prm_[Jrelp_0] = 1.;
            this->prm_[Jupnp_0] = 0.004375;
            this->prm_[Jupp_0] = 0.01203125;
            this->prm_[cmdnmax] = 0.05;

            break;

        case CellType::mid :

            // Set the model parameters for midmyocardium cell type.
            this->prm_[gNaL] = 0.0075;
            this->prm_[delta_factor] = 0.;
            this->prm_[gto] = 0.08;
            this->prm_[pCa] = 0.00025;
            this->prm_[gKr] = 0.0368;
            this->prm_[gKs] = 0.0034;
            this->prm_[gK1] = 0.24804;
            this->prm_[gncx] = 0.00112;
            this->prm_[pNaK] = 21.;
            this->prm_[gKb]  = 0.003;
            this->prm_[Jrel_0] = 1.7;
            this->prm_[Jrelp_0] = 1.7;
            this->prm_[Jupnp_0] = 0.004375;
            this->prm_[Jupp_0] = 0.01203125;
            this->prm_[cmdnmax] = 0.05;

            break;

        case CellType::epi : 

            // Set the model parameters for epicardium cell type.
            this->prm_[gNaL] = 0.0045;
            this->prm_[delta_factor] = 1.;
            this->prm_[gto] = 0.08;
            this->prm_[pCa] = 0.00012;
            this->prm_[gKr] = 0.0598;
            this->prm_[gKs] = 0.00476;
            this->prm_[gK1] = 0.22896;
            this->prm_[gncx] = 0.00088;
            this->prm_[pNaK] = 27.;
            this->prm_[gKb]  = 0.0018;
            this->prm_[Jrel_0] = 1.;
            this->prm_[Jrelp_0] = 1.;
            this->prm_[Jupnp_0] = 0.0056875;
            this->prm_[Jupp_0] =  0.015640625;
            this->prm_[cmdnmax] = 0.065;

            break;
    
        default:
            throw std::invalid_argument(Logger::Error("Could not initialize Ohara ap model. Expected: CellType::endo | CellType::mid | CellType::epi"));
            break;
    }

    // Set the cell model state variables that are independent from the cell type.
    this->var_[v]      = -87.0;
    this->var_[dvdt]   = 0.0;    
    this->var_[nai]    = 7.0;          
    this->var_[nass]   = 7.0;         
    this->var_[ki]     = 145.0;         
    this->var_[kss]    = 145.0;         
    this->var_[cai]    = 1.0e-4;
    this->var_[cass]   = 1.0e-4;
    this->var_[cansr]  = 1.2;         
    this->var_[cajsr]  = 1.2;         
    this->var_[m]      = 0.0;
    this->var_[hf]     = 1.0;          
    this->var_[hs]     = 1.0;          
    this->var_[j]      = 1.0;         
    this->var_[hsp]    = 1.0;          
    this->var_[jp]     = 1.0;          
    this->var_[mL]     = 0.0;
    this->var_[hL]     = 1.0;          
    this->var_[hLp]    = 1.0;          
    this->var_[a]      = 0.0;
    this->var_[iF]     = 1.0;          
    this->var_[iS]     = 1.0;          
    this->var_[ap]     = 0.0;
    this->var_[iFp]    = 1.0;          
    this->var_[iSp]    = 1.0;         
    this->var_[d]      = 0.0;
    this->var_[ff]     = 1.0;    
    this->var_[fs]     = 1.0;          
    this->var_[fcaf]   = 1.0;    
    this->var_[fcas]   = 1.0;          
    this->var_[jca]    = 1.0;          
    this->var_[nca]    = 0.0;
    this->var_[ffp]    = 1.0;    
    this->var_[fcafp]  = 1.0;    
    this->var_[xrf]    = 0.0;
    this->var_[xrs]    = 0.0;          
    this->var_[xs1]    = 0.0;          
    this->var_[xs2]    = 0.0;
    this->var_[xk1]    = 1.0;          
    this->var_[Jrelnp] = 0.0;
    this->var_[Jrelp]  = 0.0;
    this->var_[CaMKt]  = 0.0;  

    // Set the cell model constants that are independent from the cell type.
    this->prm_[nao] = 140.;
    this->prm_[cao] = 1.8;
    this->prm_[ko] = 5.4;
    this->prm_[R] = 8314.;
    this->prm_[T] = 310.;
    this->prm_[Fdy] = 96485.;
    this->prm_[l] =  0.01;
    this->prm_[rad] =  0.0011;
    this->prm_[vcell] = 1000. * 3.14 * this->prm_[rad]*this->prm_[rad] * this->prm_[l];
    this->prm_[ageo] = 2.*3.14*this->prm_[rad]*this->prm_[rad] + 2.*3.14*this->prm_[rad]*this->prm_[l];
    this->prm_[acap] = 2.*this->prm_[ageo];
    this->prm_[vmyo] = 0.68 * this->prm_[vcell];
    this->prm_[vnsr] = 0.0552 * this->prm_[vcell];
    this->prm_[vjsr] = 0.0048 * this->prm_[vcell];
    this->prm_[vss] = 0.02 * this->prm_[vcell];
    this->prm_[KmCaMK] = 0.15;
    this->prm_[aCaMK] = 0.05;
    this->prm_[bCaMK] = 0.00068;
    this->prm_[CaMKo] = 0.05;
    this->prm_[KmCaM] = 0.0015;
    this->prm_[BSRmax] = 0.047;
    this->prm_[KmBSR] = 0.00087;
    this->prm_[BSLmax] = 1.124;
    this->prm_[KmBSL] = 0.0087;
    this->prm_[Kmcmdn] = 0.00238;
    this->prm_[trpnmax] = 0.07;
    this->prm_[Kmtrpn] = 0.0005;
    this->prm_[Csqnmax] = 10.;
    this->prm_[kmcsqn] = 0.8;
    this->prm_[gNa] = 75.0;

}


void Ohara::Compute(double v_new, double dt, double stim_current)
{
    using namespace OhrVar;
    using namespace OhrPrm;
    using namespace OhrCur;

    // Compute the CaMKb.
    double CaMKb = this->prm_[CaMKo] * (1. - this->var_[CaMKt]) / (1. + this->prm_[KmCaM]/this->var_[cass]);

    // Compute the CaMKa.
    double CaMKa = CaMKb + this->var_[CaMKt];

    // Compute reversal potentials.
    double PKNa = 0.01833;
    double ENa = (this->prm_[R] * this->prm_[T] / this->prm_[Fdy]) * std::log(this->prm_[nao] / this->var_[nai]);
    double EK = (this->prm_[R] * this->prm_[T] / this->prm_[Fdy]) *  std::log(this->prm_[ko] / this->var_[ki]);    
    double EKs = (this->prm_[R] * this->prm_[T] / this->prm_[Fdy]) * std::log((this->prm_[ko] + PKNa*this->prm_[nao]) / (this->var_[ki] + PKNa*this->var_[nai]) );
                                                    
    // Convenient short-hand calculations.
    double vffrt = v_new * this->prm_[Fdy]*this->prm_[Fdy] / (this->prm_[R]*this->prm_[T]);
    double vfrt  = v_new * this->prm_[Fdy] / (this->prm_[R]*this->prm_[T]);

    // Compute the INa current with the original formulation form OHara-Rudy.
    double mss=1.0/(1.0+std::exp((-(v_new+39.57))/9.871));
    double tm=1.0/(6.765*std::exp((v_new+11.64)/34.77)+8.552*exp(-(v_new+77.42)/5.955));
    this->var_[m] = ALGORITHM::RushLarsen(mss, this->var_[m], dt, tm);
    
    double hss=1.0/(1+std::exp((v_new+82.90)/6.086));
    double thf=1.0/(1.432e-5*std::exp(-(v_new+1.196)/6.285)+6.149*std::exp((v_new+0.5096)/20.27));
    double ths=1.0/(0.009794*std::exp(-(v_new+17.95)/28.05)+0.3343*std::exp((v_new+5.730)/56.66));
    double Ahf=0.99;
    double Ahs=1.0-Ahf;
    this->var_[hf] = ALGORITHM::RushLarsen(hss, this->var_[hf], dt, thf);
    this->var_[hs] = ALGORITHM::RushLarsen(hss, this->var_[hs], dt, ths);
    double h=Ahf*this->var_[hf]+Ahs*this->var_[hs];
    
    double jss=hss;
    double tj=2.038+1.0/(0.02136*std::exp(-(v_new+100.6)/8.281)+0.3052*std::exp((v_new+0.9941)/38.45));
    this->var_[j] = ALGORITHM::RushLarsen(jss, this->var_[j], dt, tj);
    
    double hssp=1.0/(1+std::exp((v_new+89.1)/6.086));
    double thsp=3.0*ths;
    this->var_[hsp] = ALGORITHM::RushLarsen(hssp, this->var_[hsp], dt, thsp);
    
    double hp=Ahf*this->var_[hf]+Ahs*this->var_[hsp];
    double tjp=1.46*tj;
    this->var_[jp] = ALGORITHM::RushLarsen(jss, this->var_[jp], dt, tjp);
    
    double fINap=(1.0/(1.0+this->prm_[KmCaMK]/CaMKa));
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INa] = (1.0-this->block_coeff_[INa]) * (this->prm_[gNa]*(v_new-ENa)*std::pow(this->var_[m],3.0)*((1.0-fINap)*h*this->var_[j]+fINap*hp*this->var_[jp]));
    #else
        this->cur_[INa] = this->prm_[gNa]*(v_new-ENa)*std::pow(this->var_[m],3.0)*((1.0-fINap)*h*this->var_[j]+fINap*hp*this->var_[jp]);
    #endif
    
    // Update the mL-gate using the Rush-Larsen method.
    double mLss = 1. / (1. + std::exp((-(v_new + 42.85)) / 5.264));
    double tmL = tm;
    this->var_[mL] = ALGORITHM::RushLarsen(mLss, this->var_[mL], dt, tmL);

    // Update the hL-gate with the Rush-Larsen method.
    double hLss = 1. / (1. + std::exp((v_new+87.61)/7.488));
    double thL = 200.;
    this->var_[hL] = ALGORITHM::RushLarsen(hLss, this->var_[hL], dt, thL);
    
    // Update the hLp-gate with Rush-Larsen method.
    double hLssp = 1. / (1. + std::exp((v_new + 93.81) / 7.488));
    double thLp = 3.*thL;
    this->var_[hLp] = ALGORITHM::RushLarsen(hLssp, this->var_[hLp], dt, thLp);
    
    // Compute the INaL current.
    double fINaLp = 1. / (1. + this->prm_[KmCaMK]/CaMKa);
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INaL] = (1.-this->block_coeff_[INaL]) * (this->prm_[gNaL] * (v_new - ENa) * this->var_[mL] * ((1. - fINaLp)*this->var_[hL] + fINaLp*this->var_[hLp]));
    #else
        this->cur_[INaL] = (this->prm_[gNaL] * (v_new - ENa) * this->var_[mL] * ((1. - fINaLp)*this->var_[hL] + fINaLp*this->var_[hLp]));
    #endif

    // Update the a-gate with the Rush-Larsen method.
    double ass = 1. / (1. + std::exp((-(v_new - 14.34)) / 14.82));
    double ta = 1.0515 / (1./(1.2089*(1. + std::exp(-(v_new - 18.4099) / 29.3814))) + 3.5/(1. + std::exp((v_new + 100.) / 29.3814)));
    this->var_[a] = ALGORITHM::RushLarsen(ass, this->var_[a], dt, ta);

    double iss = 1. / (1.+std::exp((v_new + 43.94) / 5.711));
    double delta_epi = 1. - this->prm_[delta_factor]*(0.95 / (1.+std::exp((v_new + 70.) / 5.)));
    
    // Update the iF with the Rush-Larsen method.
    double tiF = 4.562 + 1./(0.3933*std::exp((-(v_new + 100.)) / 100.) + 0.08004*std::exp((v_new + 50.) / 16.59));
    tiF *= delta_epi;
    double AiF = 1. / (1. + std::exp((v_new - 213.6) / 151.2));
    this->var_[iF] = ALGORITHM::RushLarsen(iss, this->var_[iF], dt, tiF);
    
    // Update the iS with the Rush-Larsen method.
    double tiS = 23.62 + 1./(0.001416*std::exp((-(v_new + 96.52)) / 59.05) + 1.78e-8*std::exp((v_new + 114.1) / 8.079));
    tiS *= delta_epi;
    double AiS = 1. - AiF;
    this->var_[iS] = ALGORITHM::RushLarsen(iss, this->var_[iS], dt, tiS);

    double i = AiF*this->var_[iF] + AiS*this->var_[iS];
    
    // Update the ap with the Rush-Larsen method.
    double assp = 1. / (1. + std::exp((-(v_new - 24.34)) / 14.82));
    this->var_[ap] = ALGORITHM::RushLarsen(assp, this->var_[ap], dt, ta);
    
    double dti_develop = 1.354 + 1.e-4/(std::exp((v_new - 167.4) / 15.89) + std::exp(-(v_new - 12.23) / 0.2154));
    double dti_recover = 1. - 0.5 / (1. + std::exp((v_new + 70.) / 20.));
    
    // Update the iFp with the Rush-Larsen method.
    double tiFp = dti_develop * dti_recover * tiF;
    this->var_[iFp] = ALGORITHM::RushLarsen(iss, this->var_[iFp], dt, tiFp);
    
    // Update the iSp with the Rush-Larsen method.
    double tiSp = dti_develop * dti_recover * tiS;
    this->var_[iSp] = ALGORITHM::RushLarsen(iss, this->var_[iSp], dt, tiSp);
    double ip = AiF*this->var_[iFp] + AiS*this->var_[iSp];
    
    // Compute the Ito current.
    double fItop = 1. / (1. + this->prm_[KmCaMK]/CaMKa);
    #ifdef BLOCK_CELL_CURRS
        this->cur_[Ito] = (1. -this->block_coeff_[Ito]) * (this->prm_[gto] * (v_new - EK) * ((1.-fItop)*this->var_[a]*i + fItop*this->var_[ap]*ip));
    #else
        this->cur_[Ito] = (this->prm_[gto] * (v_new - EK) * ((1.-fItop)*this->var_[a]*i + fItop*this->var_[ap]*ip));
    #endif

    // Update the d-gate with the Rush-Larsen method.
    double dss = 1. / (1. + std::exp((-(v_new + 3.94)) / 4.23));
    double td = 0.6 + 1. / (std::exp(-0.05*(v_new + 6.)) + std::exp(0.09*(v_new + 14.)));
    this->var_[d] = ALGORITHM::RushLarsen(dss, this->var_[d], dt, td);
    
    // Update the ff-gate with Rush Larsen method.
    double fss = 1. / (1.+std::exp((v_new + 19.58) / 3.696));
    double tff = 7. + 1./(0.0045*std::exp(-(v_new + 20.) / 10.) + 0.0045*std::exp((v_new + 20.) / 10.));
    this->var_[ff] = ALGORITHM::RushLarsen(fss, this->var_[ff], dt, tff);

    // Update the fs-gate with the Rush-Larsen method.
    double tfs = 1000. + 1./(0.000035*std::exp(-(v_new + 5.) / 4.) + 0.000035*std::exp((v_new + 5.) / 6.));
    this->var_[fs] = ALGORITHM::RushLarsen(fss, this->var_[fs], dt, tfs);

    double Aff = 0.6;
    double Afs = 1. - Aff;
    double f = Aff*this->var_[ff] + Afs*this->var_[fs];

    // Update the fcaf-gate with the Rush-Larsen method.
    double fcass = fss;
    double tfcaf = 7. + 1./(0.04 * (std::exp(-(v_new - 4.) / 7.) + std::exp((v_new - 4.) / 7.)) );
    this->var_[fcaf] = ALGORITHM::RushLarsen(fcass, this->var_[fcaf], dt, tfcaf);
    
    // Update the fcas-gate with the Rush-Larsen method.
    double tfcas = 100. + 1./(0.00012 * (std::exp(-v_new/3.) + std::exp(v_new/7.)) );
    this->var_[fcas] = ALGORITHM::RushLarsen(fcass, this->var_[fcas], dt, tfcas);

    double Afcaf = 0.3 + 0.6/(1. + std::exp(0.1*(v_new - 10.)));
    double Afcas = 1. - Afcaf;
    double fca = Afcaf*this->var_[fcaf] + Afcas*this->var_[fcas];

    // Update the jca-gate with the Rush-Larsen method.
    double tjca = 75.;
    this->var_[jca] = ALGORITHM::RushLarsen(fcass, this->var_[jca], dt, tjca);
    
    // Update the ffp-gate with the Rush-Larsen method.
    double tffp = 2.5*tff;
    this->var_[ffp] = ALGORITHM::RushLarsen(fss, this->var_[ffp], dt, tffp);

    double fp = Aff*this->var_[ffp] + Afs*this->var_[fs];

    // Update the fcafp with the Rush-Larsen method.
    double tfcafp = 2.5*tfcaf;
    this->var_[fcafp] = ALGORITHM::RushLarsen(fcass, this->var_[fcafp], dt, tfcafp);

    double fcap = Afcaf*this->var_[fcafp] + Afcas*this->var_[fcas];
    
    // Update the nca-gate with the Rush-Larsen method.
    double Kmn = 0.002;
    double k2n = 1000.;
    double km2n = this->var_[jca];
    double anca = 1. / (k2n/km2n + std::pow(1.+Kmn/this->var_[cass],4.));
    this->var_[nca] = ALGORITHM::ForwardEuler(this->var_[nca], dt, anca*k2n-this->var_[nca]*km2n);

    double PhiCaL  = 4.*vffrt * (this->var_[cass]*std::exp(2.*vfrt) - 0.341*this->prm_[cao]) / (std::exp(2.*vfrt) - 1.);
    double PhiCaNa = vffrt * (0.75*this->var_[nass]*std::exp(vfrt) - 0.75*this->prm_[nao]) / (std::exp(vfrt) - 1.);
    double PhiCaK  = vffrt * (0.75*this->var_[kss]*std::exp(vfrt) - 0.75*this->prm_[ko]) / (std::exp(vfrt) - 1.);
    double zca = 2.;
    double PCap = 1.1 * this->prm_[pCa];
    double PCaNa = 0.00125 * this->prm_[pCa];
    double PCaK = 3.574e-4 * this->prm_[pCa];
    double PCaNap = 0.00125*PCap;
    double PCaKp = 3.574e-4*PCap;
    double fICaLp = 1. / (1. + this->prm_[KmCaMK] / CaMKa);
    
    // Calculate the ICaL current, ICaNa current and ICaK current
    #ifdef BLOCK_CELL_CURRS
        this->cur_[ICaL] =  (1. -this->block_coeff_[ICaL]) * ((1. - fICaLp)*this->prm_[pCa]*PhiCaL*this->var_[d] * (f*(1. - this->var_[nca]) + this->var_[jca]*fca*this->var_[nca]) + fICaLp*PCap*PhiCaL*this->var_[d] * (fp*(1. - this->var_[nca]) + this->var_[jca]*fcap*this->var_[nca]));
        this->cur_[ICaNa]= (1. -this->block_coeff_[ICaNa]) * ((1. - fICaLp)*PCaNa*PhiCaNa*this->var_[d] * (f*(1. - this->var_[nca]) + this->var_[jca]*fca*this->var_[nca]) + fICaLp*PCaNap*PhiCaNa*this->var_[d]*(fp*(1. - this->var_[nca]) + this->var_[jca]*fcap*this->var_[nca]));
        this->cur_[ICaK] =  (1. -this->block_coeff_[ICaK]) * ((1. - fICaLp)*PCaK*PhiCaK*this->var_[d] * (f*(1. - this->var_[nca]) + this->var_[jca]*fca*this->var_[nca]) + fICaLp*PCaKp*PhiCaK*this->var_[d]*(fp*(1. - this->var_[nca]) + this->var_[jca]*fcap*this->var_[nca]));
    #else
        this->cur_[ICaL] = ((1. - fICaLp)*this->prm_[pCa]*PhiCaL*this->var_[d] * (f*(1. - this->var_[nca]) + this->var_[jca]*fca*this->var_[nca]) + fICaLp*PCap*PhiCaL*this->var_[d] * (fp*(1. - this->var_[nca]) + this->var_[jca]*fcap*this->var_[nca]));
        this->cur_[ICaNa]= ((1. - fICaLp)*PCaNa*PhiCaNa*this->var_[d] * (f*(1. - this->var_[nca]) + this->var_[jca]*fca*this->var_[nca]) + fICaLp*PCaNap*PhiCaNa*this->var_[d]*(fp*(1. - this->var_[nca]) + this->var_[jca]*fcap*this->var_[nca]));
        this->cur_[ICaK] = ((1. - fICaLp)*PCaK*PhiCaK*this->var_[d] * (f*(1. - this->var_[nca]) + this->var_[jca]*fca*this->var_[nca]) + fICaLp*PCaKp*PhiCaK*this->var_[d]*(fp*(1. - this->var_[nca]) + this->var_[jca]*fcap*this->var_[nca]));
    #endif

    // Update the xrf-gate with the Rush-Larsen method.
    double xrss = 1. / (1. + std::exp((-(v_new+8.337))/6.789));
    double txrf = 12.98 + 1. / (0.3652*std::exp((v_new - 31.66) / 3.869) + 4.123e-5*std::exp((-(v_new - 47.78)) / 20.38));
    this->var_[xrf] = ALGORITHM::RushLarsen(xrss, this->var_[xrf], dt, txrf);
    
    // Update the xrs with the Rush-Larsen method.
    double txrs = 1.865 + 1. / (0.06629*std::exp((v_new - 34.7) / 7.355) + 1.128e-5*std::exp((-(v_new - 29.74)) / 25.94));
    this->var_[xrs] = ALGORITHM::RushLarsen(xrss, this->var_[xrs], dt, txrs);

    double Axrf = 1. / (1.+std::exp((v_new+54.81)/38.21));
    double Axrs = 1. - Axrf;
    double xr = Axrf*this->var_[xrf] + Axrs*this->var_[xrs];

    // Calculate the IKr current.
    double rkr = 1. / (1. + std::exp((v_new + 55.) / 75.)) * 1. / (1. + std::exp((v_new - 10.) / 30.));
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IKr] = (1.-this->block_coeff_[IKr]) * (this->prm_[gKr] * std::sqrt(this->prm_[ko]/5.4) * xr * rkr * (v_new - EK));
    #else
        this->cur_[IKr] = (this->prm_[gKr] * std::sqrt(this->prm_[ko]/5.4) * xr * rkr * (v_new - EK));
    #endif

    // Update the xs1-gate with the Rush-Larsen method.
    double xs1ss = 1. / (1. + std::exp((-(v_new + 11.6)) / 8.932));
    double txs1 = 817.3 + 1./(2.326e-4*std::exp((v_new + 48.28) / 17.8) + 0.001292*std::exp((-(v_new + 210.)) / 230.));
    this->var_[xs1] = ALGORITHM::RushLarsen(xs1ss, this->var_[xs1], dt, txs1);
    
    // Update the xs2-gate with the Rush-Larsen method.
    double xs2ss = xs1ss;
    double txs2 = 1. / (0.01*std::exp((v_new - 50.) / 20.) + 0.0193*std::exp((-(v_new + 66.54)) / 31.));
    this->var_[xs2] = ALGORITHM::RushLarsen(xs2ss, this->var_[xs2], dt, txs2);

    // Compute the IKs current.
    double KsCa = 1. + 0.6 / (1. + std::pow((3.8e-5 / this->var_[cai]), 1.4));
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IKs] = (1.-this->block_coeff_[IKs]) * (this->prm_[gKs] * KsCa * this->var_[xs1]* this->var_[xs2] * (v_new - EKs));
    #else
        this->cur_[IKs] = (this->prm_[gKs] * KsCa * this->var_[xs1]* this->var_[xs2] * (v_new - EKs));
    #endif

    // Update the xk1-gate with the Rush-Larsen method.
    double xk1ss = 1. / (1. + std::exp(-(v_new + 2.5538*this->prm_[ko] + 144.59) / (1.5692*this->prm_[ko] + 3.8115)));
    double txk1 = 122.2 / (std::exp((-(v_new + 127.2)) / 20.36) + std::exp((v_new + 236.8) / 69.33));
    this->var_[xk1] = ALGORITHM::RushLarsen(xk1ss, this->var_[xk1], dt, txk1);

    // Compute the IK1 current.
    double rk1 = 1. / (1. + std::exp((v_new + 105.8 - 2.6*this->prm_[ko]) / 9.493));
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IK1] = (1. -this->block_coeff_[IK1]) * (this->prm_[gK1] * std::sqrt(this->prm_[ko]) * rk1 * this->var_[xk1] * (v_new - EK));
    #else
        this->cur_[IK1] = (this->prm_[gK1] * std::sqrt(this->prm_[ko]) * rk1 * this->var_[xk1] * (v_new - EK));
    #endif

    double kna1 = 15.;     double kna2 = 5.;      double kna3 = 88.12;
    double kasymm = 12.5;  double wna = 6.e4;     double wca = 6.e4;
    double wnaca = 5.e3;   double kcaon = 1.5e6;  double kcaoff = 5.e3;
    double qna = 0.5224;   double qca = 0.1670;   
    double hca = std::exp((qca*v_new*this->prm_[Fdy]) / (this->prm_[R]*this->prm_[T]));
    double hna = std::exp((qna*v_new*this->prm_[Fdy]) / (this->prm_[R]*this->prm_[T]));
    double h1 = 1. + (this->var_[nai] / kna3*(1. + hna));
    double h2 = (this->var_[nai] * hna) / (kna3 * h1);
    double h3 = 1. / h1;
    double h4 = 1. + this->var_[nai]/kna1*(1. + this->var_[nai]/kna2);
    double h5 = this->var_[nai]*this->var_[nai] / (h4*kna1*kna2);
    double h6 = 1. / h4;
    double h7 = 1. + this->prm_[nao]/kna3*(1. + 1./hna);
    double h8 = this->prm_[nao] / (kna3*hna*h7);
    double h9 = 1. / h7;
    double h10 = kasymm + 1. + this->prm_[nao]/kna1*(1. + this->prm_[nao]/kna2);
    double h11 = this->prm_[nao]*this->prm_[nao] / (h10*kna1*kna2);
    double h12 = 1. / h10;
    double k1 = h12 * this->prm_[cao] * kcaon;
    double k2 = kcaoff;
    double k3p = h9 * wca;
    double k3pp = h8 * wnaca;
    double k3 = k3p + k3pp;
    double k4p = h3 *wca / hca;
    double k4pp = h2 * wnaca;
    double k4 = k4p + k4pp;
    double k5 = kcaoff;
    double k6 = h6 * this->var_[cai] * kcaon;
    double k7 = h5 * h2 * wna;
    double k8 = h8 * h11 * wna;
    double x1 = k2*k4*(k7+k6) + k5*k7*(k2+k3);
    double x2 = k1*k7*(k4+k5) + k4*k6*(k1+k8);
    double x3 = k1*k3*(k7+k6) + k8*k6*(k2+k3);
    double x4 = k2*k8*(k4+k5) + k3*k5*(k1+k8);
    double E1 = x1 / (x1+x2+x3+x4);
    double E2 = x2 / (x1+x2+x3+x4);
    double E3 = x3 / (x1+x2+x3+x4);
    double E4 = x4 / (x1+x2+x3+x4);
    double KmCaAct = 0.00015;
    double allo = 1. / (1. + (KmCaAct/this->var_[cai])*(KmCaAct/this->var_[cai]));
    double zna = 1.;
    double JncxNa = 3.*(E4*k7 - E1*k8) + E3*k4pp - E2*k3pp;
    double JncxCa = E2*k2 - E1*k1;

    // Compute the INaCa_i current.
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INaCa_i] = (1.-this->block_coeff_[INaCa_i]) * (0.8 * this->prm_[gncx] * allo * (zna*JncxNa + zca*JncxCa));
    #else
        this->cur_[INaCa_i] = (0.8 * this->prm_[gncx] * allo * (zna*JncxNa + zca*JncxCa));
    #endif

    // Compute the INaCa_ss current.
    h1 = 1. + this->var_[nass]/kna3*(1. + hna);
    h2 = (this->var_[nass]*hna) / (kna3*h1);
    h3 = 1. / h1;
    h4 = 1. + this->var_[nass]/kna1*(1. + this->var_[nass]/kna2);
    h5 = this->var_[nass]*this->var_[nass] / (h4*kna1*kna2);
    h6 = 1. / h4;
    h7 = 1. + this->prm_[nao]/kna3*(1. + 1./hna);
    h8 = this->prm_[nao] / (kna3*hna*h7);
    h9 = 1. / h7;
    h10 = kasymm + 1. + this->prm_[nao]/kna1*(1. + this->prm_[nao]/kna2);
    h11 = this->prm_[nao]*this->prm_[nao] / (h10*kna1*kna2);
    h12 = 1. / h10;
    k1 = h12 * this->prm_[cao] * kcaon;
    k2 = kcaoff;
    k3p = h9 * wca;
    k3pp = h8 * wnaca;
    k3 = k3p + k3pp;
    k4p = h3 * wca / hca;
    k4pp = h2 * wnaca;
    k4 = k4p + k4pp;
    k5 = kcaoff;
    k6 = h6 * this->var_[cass] * kcaon;
    k7 = h5 * h2 * wna;
    k8 = h8 * h11 * wna;
    x1 = k2*k4*(k7+k6) + k5*k7*(k2+k3);
    x2 = k1*k7*(k4+k5) + k4*k6*(k1+k8);
    x3 = k1*k3*(k7+k6) + k8*k6*(k2+k3);
    x4 = k2*k8*(k4+k5) + k3*k5*(k1+k8);
    E1 = x1 / (x1+x2+x3+x4);
    E2 = x2 / (x1+x2+x3+x4);
    E3 = x3 / (x1+x2+x3+x4);
    E4 = x4 / (x1+x2+x3+x4);
    KmCaAct = 150.e-6;
    allo = 1. / (1. + (KmCaAct/this->var_[cass])*(KmCaAct/this->var_[cass]));
    JncxNa = 3.*(E4*k7-E1*k8) + E3*k4pp - E2*k3pp;
    JncxCa = E2*k2 - E1*k1;

    // Compute the INaCa_ss current and total INaCa current.
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INaCa_ss] = (1. -this->block_coeff_[INaCa_ss]) * (0.2 * this->prm_[gncx] * allo * (zna*JncxNa + zca*JncxCa));
        this->cur_[INaCa]    = (1. -this->block_coeff_[INaCa]) * (this->cur_[INaCa_i] + this->cur_[INaCa_ss]);
    #else
        this->cur_[INaCa_ss] = (0.2 * this->prm_[gncx] * allo * (zna*JncxNa + zca*JncxCa));
        this->cur_[INaCa]    = (this->cur_[INaCa_i] + this->cur_[INaCa_ss]);
    #endif

    double k1p = 949.5;   double k1m = 182.4;  double k2p = 687.2;
    double k2m = 39.4;    
    k3p = 1899.;
    double k3m = 79300.;
    k4p = 639.;
    double k4m = 40.;    double Knai0 = 9.073;
    double Knao0 = 27.78; double delta = -0.155;

    double Knai = Knai0 * std::exp((delta*v_new*this->prm_[Fdy]) / (3.*this->prm_[R]*this->prm_[T]));
    double Knao = Knao0 * std::exp(((1.-delta)*v_new*this->prm_[Fdy]) / (3.*this->prm_[R]*this->prm_[T]));

    double Kki = 0.5;   double Kko = 0.3582;      double MgADP = 0.05;
    double MgATP = 9.8; double Kmgatp = 1.698e-7; double HH = 1.e-7;
    double eP = 4.2;    double Khp = 1.698e-7;    double Knap = 224.;
    double Kxkur = 292.;
    
    double P = eP / (1. + HH/Khp + this->var_[nai]/Knap + this->var_[ki]/Kxkur);
    double a1 = (k1p * std::pow(this->var_[nai]/Knai, 3.)) / 
                (std::pow(1. + this->var_[nai]/Knai, 3.) + (1. + this->var_[ki]/Kki)*(1. + this->var_[ki]/Kki) - 1.);

    double b1 = k1m*MgADP;
    double a2 = k2p;
    double b2 = (k2m * std::pow(this->prm_[nao]/Knao, 3.)) / 
                (std::pow(1. + this->prm_[nao]/Knao, 3.) + std::pow(1. + this->prm_[ko]/Kko, 2.) - 1.);
    
    double a3 = (k3p*(this->prm_[ko]/Kko)*(this->prm_[ko]/Kko)) / 
                (std::pow(1. + this->prm_[nao]/Knao, 3.) + (1. + this->prm_[ko]/Kko)*(1. + this->prm_[ko]/Kko) - 1.);

    double b3 = (k3m*P*HH) / (1. + MgATP/Kmgatp);
    double a4 = (k4p*MgATP/Kmgatp) / (1. + MgATP/Kmgatp);
    double b4 = (k4m*(this->var_[ki]/Kki)*(this->var_[ki]/Kki)) / 
                (std::pow(1. + this->var_[nai]/Knai, 3.) + (1. + this->var_[ki]/Kki)*(1. + this->var_[ki]/Kki) - 1.);

    x1 = a4*a1*a2 + b2*b4*b3 + a2*b4*b3 + b3*a1*a2;
    x2 = b2*b1*b4 + a1*a2*a3 + a3*b1*b4 + a2*a3*b4;
    x3 = a2*a3*a4 + b3*b2*b1 + b2*b1*a4 + a3*a4*b1;
    x4 = b4*b3*b2 + a3*a4*a1 + b2*a4*a1 + b3*b2*a1;
    E1 = x1 / (x1+x2+x3+x4);
    E2 = x2 / (x1+x2+x3+x4);
    E3 = x3 / (x1+x2+x3+x4);
    E4 = x4 / (x1+x2+x3+x4);
    double zk = 1.;
    double JnakNa = 3.*(E1*a3 - E2*b3);
    double JnakK = 2.*(E4*b1 - E3*a1);

    // Compute the INaK current.
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INaK] = (1. -this->block_coeff_[INaK]) * (this->prm_[pNaK] * (zna*JnakNa + zk*JnakK));
    #else
        this->cur_[INaK] = (this->prm_[pNaK] * (zna*JnakNa + zk*JnakK));
    #endif

    // Compute the IKb current.
    double xkb = 1. / (1. + std::exp(-(v_new-14.48)/18.34));
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IKb] = (1. -this->block_coeff_[IKb]) * (this->prm_[gKb] * xkb * (v_new - EK));
    #else
        this->cur_[IKb] = (this->prm_[gKb] * xkb * (v_new - EK));
    #endif

    // Compute the INab current.
    double PNab = 3.75e-10;
    #ifdef BLOCK_CELL_CURRS
        this->cur_[INab] = (1. -this->block_coeff_[INab]) * (PNab * vffrt * (this->var_[nai]*std::exp(vfrt) - this->prm_[nao]) / (std::exp(vfrt) - 1.));
    #else
        this->cur_[INab] = (PNab * vffrt * (this->var_[nai]*std::exp(vfrt) - this->prm_[nao]) / (std::exp(vfrt) - 1.));
    #endif

    // Compute the ICab current.
    double PCab = 2.5e-8;
    #ifdef BLOCK_CELL_CURRS
        this->cur_[ICab] = (1. -this->block_coeff_[ICab]) * (PCab * 4.*vffrt * (this->var_[cai]*std::exp(2.*vfrt) - 0.341*this->prm_[cao]) / (std::exp(2.*vfrt) - 1.));
    #else
        this->cur_[ICab] = (PCab * 4.*vffrt * (this->var_[cai]*std::exp(2.*vfrt) - 0.341*this->prm_[cao]) / (std::exp(2.*vfrt) - 1.));
    #endif

    // Compute the IpCa current.
    double GpCa = 0.0005;
    #ifdef BLOCK_CELL_CURRS
        this->cur_[IpCa] = (1. -this->block_coeff_[IpCa]) * (GpCa * this->var_[cai] / (0.0005 + this->var_[cai]));
    #else
        this->cur_[IpCa] = (GpCa * this->var_[cai] / (0.0005 + this->var_[cai]));
    #endif

    // Compute the total Iion current.
    this->cur_[OhrCur::Iion] = this->cur_[INa]  + this->cur_[INaL] +  this->cur_[Ito] + this->cur_[ICaL] + this->cur_[ICaNa] + 
                               this->cur_[ICaK] + this->cur_[IKr]  +  this->cur_[IKs] + this->cur_[IK1]  + this->cur_[INaCa] +
                               this->cur_[INaK] + this->cur_[INab] + this->cur_[IKb]  + this->cur_[IpCa] + this->cur_[ICab];

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[OhrCur::Iion] - stim_current);

    // Update the CaMKt with the Forward Euler method.
    double dCaMKt = (this->prm_[aCaMK]*CaMKb*(CaMKb+this->var_[CaMKt]) - this->prm_[bCaMK]*this->var_[CaMKt]);
    this->var_[CaMKt] = ALGORITHM::ForwardEuler(this->var_[CaMKt], dt, dCaMKt);

    // Compute the diffusion fluxes.
    double JdiffNa = 0.5 * (this->var_[nass] - this->var_[nai]);
    double JdiffK = 0.5 * (this->var_[kss] - this->var_[ki]);
    double Jdiff = (this->var_[cass] - this->var_[cai]) / 0.2;

    // Compute the ryanodione receptor calcium induced calcium release from the jsr.
    double bt = 4.75;
    double a_rel = 0.5*bt;
    double Jrel_inf = this->prm_[Jrel_0] * (a_rel * (-this->cur_[ICaL]) / (1. + std::pow(1.5/this->var_[cajsr], 8.)));

    double tau_rel = bt / (1. + 0.0123/this->var_[cajsr]);
    if (tau_rel < 0.001) { tau_rel = 0.001; }

    // Update the Jrelnp with the Rush-Larsen method.
    this->var_[Jrelnp] = ALGORITHM::RushLarsen(Jrel_inf, this->var_[Jrelnp], dt, tau_rel);

    double btp = 1.25*bt;
    double a_relp = 0.5*btp;
    double Jrel_infp = this->prm_[Jrelp_0] * (a_relp * (-this->cur_[ICaL]) / (1. + std::pow(1.5/this->var_[cajsr], 8.)));

    double tau_relp = btp / (1. + 0.0123/this->var_[cajsr]);
    if (tau_relp < 0.001) { tau_relp = 0.001; }
    
    // Update the Jrelp with the  Rush-Larsen method.
    this->var_[Jrelp] = ALGORITHM::RushLarsen(Jrel_infp, this->var_[Jrelp], dt, tau_relp);

    double fJrelp = 1. / (1. + this->prm_[KmCaMK]/CaMKa);

    double Jrel = (1. - fJrelp) * this->var_[Jrelnp] + fJrelp*this->var_[Jrelp];
    
    // Compute the serca pump, ca uptake flux.
    double Jupnp = this->prm_[Jupnp_0] * this->var_[cai] / (this->var_[cai] + 0.00092);
    double Jupp = this->prm_[Jupp_0] * this->var_[cai] / (this->var_[cai] + 0.00075);
    double fJupp = 1. / (1. + this->prm_[KmCaMK]/CaMKa);
    double Jleak = 0.0039375 * this->var_[cansr] / 15.;
    double Jup = (1. - fJupp)*Jupnp + fJupp*Jupp - Jleak;
    
    // Compute the tranlocation flux.
    double Jtr = 0.01 * (this->var_[cansr] - this->var_[cajsr]);

    // Update the Sodium intracellular concentration with Forward Euler.
    double dnai = (-(this->cur_[INa] + this->cur_[INaL] + 3.*this->cur_[INaCa_i] + 3.*this->cur_[INaK] + this->cur_[INab])) * 
                     this->prm_[acap] / (this->prm_[Fdy] * this->prm_[vmyo]) + JdiffNa*this->prm_[vss]/this->prm_[vmyo];
    
    this->var_[nai] = ALGORITHM::ForwardEuler(this->var_[nai], dt, dnai);

    // Update the Sodium subspace compartment concentration with Forward Euler.
    double dnass = (-(this->cur_[ICaNa] + 3.*this->cur_[INaCa_ss])) * this->prm_[acap] / (this->prm_[Fdy]*this->prm_[vss]) - JdiffNa;
    
    this->var_[nass] = ALGORITHM::ForwardEuler(this->var_[nass], dt, dnass);

    // Update the Potassium intracellular concentration with Forward Euler.
    double dki = (-(this->cur_[Ito] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[IK1] + this->cur_[IKb] + (-1*stim_current) - 2.*this->cur_[INaK])) *
                     this->prm_[acap] / (this->prm_[Fdy] * this->prm_[vmyo]) + JdiffK * this->prm_[vss] / this->prm_[vmyo];
    
    this->var_[ki] = ALGORITHM::ForwardEuler(this->var_[ki], dt, dki);

    // Update the Potassium subspace compartment concentration with Forward Euler.
    double dkss = (-this->cur_[ICaK]) * this->prm_[acap] / (this->prm_[Fdy]*this->prm_[vss]) - JdiffK;
    this->var_[kss] = ALGORITHM::ForwardEuler(this->var_[kss], dt, dkss);
    
    // Update the Calcium intracellular concentration with Forward Euler.                       
    double Bcai = 1. / (1. + this->prm_[cmdnmax]*this->prm_[Kmcmdn]/std::pow(this->prm_[Kmcmdn] + this->var_[cai], 2.) +
                        this->prm_[trpnmax]*this->prm_[Kmtrpn]/std::pow(this->prm_[Kmtrpn] + this->var_[cai], 2.));

    double dcai = Bcai * (-(this->cur_[IpCa] + this->cur_[ICab] - 2.*this->cur_[INaCa_i])*this->prm_[acap]/(2.*this->prm_[Fdy]*this->prm_[vmyo]) - 
                              Jup*this->prm_[vnsr]/this->prm_[vmyo] + Jdiff*this->prm_[vss]/this->prm_[vmyo]);

    this->var_[cai] = ALGORITHM::ForwardEuler(this->var_[cai], dt, dcai);

    // Update the Calcium subspace compartment concentration with Forward Euler.    
    double Bcass = 1. / (1. + this->prm_[BSRmax]*this->prm_[KmBSR] / ((this->prm_[KmBSR] + this->var_[cass]) * (this->prm_[KmBSR] + this->var_[cass]) ) + 
                         this->prm_[BSLmax]*this->prm_[KmBSL] / ((this->prm_[KmBSL] + this->var_[cass]) * (this->prm_[KmBSL] + this->var_[cass])) );

    double dcass = Bcass * (-(this->cur_[ICaL] - 2.*this->cur_[INaCa_ss])*this->prm_[acap]/(2.*this->prm_[Fdy]*this->prm_[vss]) + 
                                Jrel*this->prm_[vjsr]/this->prm_[vss] - Jdiff);

    this->var_[cass] = ALGORITHM::ForwardEuler(this->var_[cass], dt, dcass);

    // Update the Calcium network sarcoplasmic reticulum concentration with Forward Euler.
    double dcansr = Jup - Jtr * this->prm_[vjsr] / this->prm_[vnsr];
    this->var_[cansr] = ALGORITHM::ForwardEuler(this->var_[cansr], dt, dcansr);
    
    // Update the Calcium junctional sarcoplasmic reticulum concentration with Forward Euler.
    double Bcajsr = 1. / (1. + this->prm_[Csqnmax]*this->prm_[kmcsqn] / ((this->prm_[kmcsqn] + this->var_[cajsr]) * (this->prm_[kmcsqn] + this->var_[cajsr]) ));
    double dcajsr = Bcajsr * (Jtr - Jrel);
    this->var_[cajsr] = ALGORITHM::ForwardEuler(this->var_[cajsr], dt, dcajsr);
    
}


std::string Ohara::PrintVariables() const 
{
    using namespace OhrVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "nai = " << this->var_[nai] << "\n";
    oss << "nass = " << this->var_[nass] << "\n";
    oss << "ki = " << this->var_[ki] << "\n";
    oss << "kss = " << this->var_[kss] << "\n";
    oss << "cai = " << this->var_[cai] << "\n";
    oss << "cass = " << this->var_[cass] << "\n";
    oss << "cansr = " << this->var_[cansr] << "\n";
    oss << "cajsr = " << this->var_[cajsr] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "hf = " << this->var_[hf] << "\n";
    oss << "hs = " << this->var_[hs] << "\n";
    oss << "j = " << this->var_[j] << "\n";
    oss << "hsp = " << this->var_[hsp] << "\n";
    oss << "jp = " << this->var_[jp] << "\n";
    oss << "mL = " << this->var_[mL] << "\n";
    oss << "hL = " << this->var_[hL] << "\n";
    oss << "hLp = " << this->var_[hLp] << "\n";
    oss << "a = " << this->var_[a] << "\n";
    oss << "iF = " << this->var_[iF] << "\n";
    oss << "iS = " << this->var_[iS] << "\n";
    oss << "ap = " << this->var_[ap] << "\n";
    oss << "iFp = " << this->var_[iFp] << "\n";
    oss << "iSp = " << this->var_[iSp] << "\n";
    oss << "d = " << this->var_[d] << "\n";
    oss << "ff = " << this->var_[ff] << "\n";
    oss << "fs = " << this->var_[fs] << "\n";
    oss << "fcaf = " << this->var_[fcaf] << "\n";
    oss << "fcas = " << this->var_[fcas] << "\n";
    oss << "jca = " << this->var_[jca] << "\n";
    oss << "nca = " << this->var_[nca] << "\n";
    oss << "ffp = " << this->var_[ffp] << "\n";
    oss << "fcafp = " << this->var_[fcafp] << "\n";
    oss << "xrf = " << this->var_[xrf] << "\n";
    oss << "xrs = " << this->var_[xrs] << "\n";
    oss << "xs1 = " << this->var_[xs1] << "\n";
    oss << "xs2 = " << this->var_[xs2] << "\n";
    oss << "xk1 = " << this->var_[xk1] << "\n";
    oss << "Jrelnp = " << this->var_[Jrelnp] << "\n";
    oss << "Jrelp = " << this->var_[Jrelp] << "\n";
    oss << "CaMKt = " << this->var_[CaMKt];
    return oss.str();

}


std::string Ohara::PrintParameters() const 
{
    using namespace OhrPrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "gNaL = " << this->prm_[gNaL] << "\n";
    oss << "delta_factor = " << this->prm_[delta_factor] << "\n";
    oss << "gto = " << this->prm_[gto] << "\n";
    oss << "pCa = " << this->prm_[pCa] << "\n";
    oss << "gKr = " << this->prm_[gKr] << "\n";
    oss << "gKs = " << this->prm_[gKs] << "\n";
    oss << "gK1 = " << this->prm_[gK1] << "\n";
    oss << "gncx = " << this->prm_[gncx] << "\n";
    oss << "pNaK = " << this->prm_[pNaK] << "\n";
    oss << "gKb = " << this->prm_[gKb] << "\n";
    oss << "Jrel_0 = " << this->prm_[Jrel_0] << "\n";
    oss << "Jrelp_0 = " << this->prm_[Jrelp_0] << "\n";
    oss << "Jupnp_0 = " << this->prm_[Jupnp_0] << "\n";
    oss << "Jupp_0 = " << this->prm_[Jupp_0] << "\n";
    oss << "nao = " << this->prm_[nao] << "\n";
    oss << "cao = " << this->prm_[cao] << "\n";
    oss << "ko = " << this->prm_[ko] << "\n";
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "Fdy = " << this->prm_[Fdy] << "\n";
    oss << "l = " << this->prm_[l] << "\n";
    oss << "rad = " << this->prm_[rad] << "\n";
    oss << "vcell = " << this->prm_[vcell] << "\n";
    oss << "ageo = " << this->prm_[ageo] << "\n";
    oss << "acap = " << this->prm_[acap] << "\n";
    oss << "vmyo = " << this->prm_[vmyo] << "\n";
    oss << "vnsr = " << this->prm_[vnsr] << "\n";
    oss << "vjsr = " << this->prm_[vjsr] << "\n";
    oss << "vss = " << this->prm_[vss] << "\n";
    oss << "KmCaMK = " << this->prm_[KmCaMK] << "\n";
    oss << "aCaMK = " << this->prm_[aCaMK] << "\n";
    oss << "bCaMK = " << this->prm_[bCaMK] << "\n";
    oss << "CaMKo = " << this->prm_[CaMKo] << "\n";
    oss << "KmCaM = " << this->prm_[KmCaM] << "\n";
    oss << "BSRmax = " << this->prm_[BSRmax] << "\n";
    oss << "KmBSR = " << this->prm_[KmBSR] << "\n";
    oss << "BSLmax = " << this->prm_[BSLmax] << "\n";
    oss << "KmBSL = " << this->prm_[KmBSL] << "\n";
    oss << "cmdnmax = " << this->prm_[cmdnmax] << "\n";
    oss << "Kmcmdn = " << this->prm_[Kmcmdn] << "\n";
    oss << "trpnmax = " << this->prm_[trpnmax] << "\n";
    oss << "Kmtrpn = " << this->prm_[Kmtrpn] << "\n";
    oss << "Csqnmax = " << this->prm_[Csqnmax] << "\n";
    oss << "kmcsqn = " << this->prm_[kmcsqn] << "\n";
    oss << "gNa = " << this->prm_[gNa];
    return oss.str();

}


std::string Ohara::PrintCurrents() const
{
    using namespace OhrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "INaL = " << this->cur_[INaL] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "ICaNa = " << this->cur_[ICaNa] << "\n";
    oss << "ICaK = " << this->cur_[ICaK] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "INaCa_i = " << this->cur_[INaCa_i] << "\n";
    oss << "INaCa_ss = " << this->cur_[INaCa_ss] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "IKb = " << this->cur_[IKb] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "IpCa = " << this->cur_[IpCa] << "\n";
    oss << "Iion = " << this->cur_[OhrCur::Iion];
    return oss.str();

}

#ifdef BLOCK_CELL_CURRS
    std::string Ohara::PrintBlockCoeffs() const
    {
        using namespace OhrCur;

        // Create output string stream to pass the currents and their values.
        std::ostringstream oss;
        oss.precision(15);
        oss << "INa = " << this->block_coeff_[INa] << "\n";
        oss << "INaL = " << this->block_coeff_[INaL] << "\n";
        oss << "Ito = " << this->block_coeff_[Ito] << "\n";
        oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
        oss << "ICaNa = " << this->block_coeff_[ICaNa] << "\n";
        oss << "ICaK = " << this->block_coeff_[ICaK] << "\n";
        oss << "IKr = " << this->block_coeff_[IKr] << "\n";
        oss << "IKs = " << this->block_coeff_[IKs] << "\n";
        oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
        oss << "INaCa_i = " << this->block_coeff_[INaCa_i] << "\n";
        oss << "INaCa_ss = " << this->block_coeff_[INaCa_ss] << "\n";
        oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
        oss << "INaK = " << this->block_coeff_[INaK] << "\n";
        oss << "IKb = " << this->block_coeff_[IKb] << "\n";
        oss << "INab = " << this->block_coeff_[INab] << "\n";
        oss << "ICab = " << this->block_coeff_[ICab] << "\n";
        oss << "IpCa = " << this->block_coeff_[IpCa];
        return oss.str();

    }
#endif


} // End of namespace ELECTRA
