/*
 *  tc_psc_alpha.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "ht_psc_alpha.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>
#include "universal_data_logger_impl.h"

#include <iomanip>
#include <iostream>
#include <cstdio>

inline double_t vtrap(double_t x, double_t y)
{
  double_t vtrap = 0.0;
  if (std::fabs(x/y) < 1.e-6) {
    vtrap = y * (1. - x/y/2.);
  }else {
    vtrap = x / (std::exp(x/y)-1.);
  }
  return vtrap;
}


nest::RecordablesMap<nest::ht_psc_alpha> nest::ht_psc_alpha::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<ht_psc_alpha>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::V_M>);
    insert_(names::I_ex, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::I_EXC>);
    insert_(names::I_in, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::I_INH>);
    insert_(names::HT_m, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_M>);
    insert_(names::HT_h, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_H>);
    insert_(names::HT_n, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_N>);
    insert_(names::HT_a, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_A>);
    insert_(names::HT_b, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_B>);
    insert_(names::HT_c, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_C>);
    insert_(names::HT_d, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_D>);
    insert_(names::HT_e, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_E>);
    insert_(names::HT_z, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_Z>);
    insert_(names::HT_p, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_P>);
     insert_(names::HT_o2, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_O2>);
    insert_(names::HT_cai_LT, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_CAI_LT>);
    insert_(names::HT_cai_HT, 
	    &ht_psc_alpha::get_y_elem_<ht_psc_alpha::State_::HT_CAI_HT>);
  }

  extern "C"
  int ht_psc_alpha_dynamics (double, const double y[], double f[], void* pnode)
  { 
    // a shorthand
    typedef nest::ht_psc_alpha::State_ S;
  
    // get access to node so we can almost work as in a member function
    assert(pnode);
    const nest::ht_psc_alpha& node =  *(reinterpret_cast<nest::ht_psc_alpha*>(pnode));
    
    // y[] here is---and must be---the state vector supplied by the integrator,
    // not the state vector in the node, node.S_.y[]. 
    
    // The following code is verbose for the sake of clarity. We assume that a
    // good compiler will optimize the verbosity away ...
    
   // shorthand for state variables
   const double_t& V     = y[S::V_M   ];
   const double_t& m     = y[S::HT_M  ];
   const double_t& h     = y[S::HT_H  ];
   const double_t& n     = y[S::HT_N  ];
   const double_t& a     = y[S::HT_A  ];
   const double_t& b     = y[S::HT_B  ];
   const double_t& c     = y[S::HT_C  ];
   const double_t& d     = y[S::HT_D  ];
   const double_t& e     = y[S::HT_E  ];
   const double_t& z     = y[S::HT_Z  ];
   const double_t& o2    = y[S::HT_O2 ];
   const double_t& p     = y[S::HT_P ];
   const double_t& cai_LT= y[S::HT_CAI_LT];
   const double_t& cai_HT= y[S::HT_CAI_HT];
   const double_t& dI_ex = y[S::DI_EXC];
   const double_t&  I_ex = y[S::I_EXC ];
   const double_t& dI_in = y[S::DI_INH];
   const double_t&  I_in = y[S::I_INH ];

   //reversal pontential of Ca2+ for LT
   const double celsius	= 36.;
   const double R = 8.31441; //???
   const double_t FARADAY = 96485.3;
   const double cao=2.;
   const double E_Ca_LT = (1.e+3) * (R*(celsius+273.15))/(2.*FARADAY) * std::log(cao/cai_LT);

   //reversal pontential of Ca2+ for HT
   const double E_Ca_HT = (1.e+3) * (R*(celsius+273.15))/(2.*FARADAY) * std::log(cao/cai_HT);

   // Na+ channel gate variable m
   const double_t V_traub = -25.;
   const double_t v2 = V - V_traub;
   double_t temp1 = 0.32 * vtrap(13.-v2, 4.);
   double_t temp2 = 0.28 * vtrap(v2-40.,5);
   const double_t tau_m = 1. / (temp1 + temp2);
   const double_t m_inf = temp1 / (temp1 + temp2);

   // Na+ channel gate variable h
   temp1 = 0.128 * std::exp((17.-v2)/18.);
   temp2 = 4. / ( 1. + std::exp((40.-v2)/5.) );
   const double_t tau_h = 1. / (temp1 + temp2);
   const double_t h_inf = temp1 / (temp1 + temp2);

   // K+ channel gate variable n
   //temp1 = 0.032 * (15.-v2) / ( std::exp((15.-v2)/5.) - 1.);
   temp1 = 0.032 * vtrap(15. - v2, 5.);
   temp2 = 0.5 * std::exp((10.-v2)/40.);
   const double_t tau_n = 1. / (temp1 + temp2);
   const double_t n_inf = temp1 / (temp1 + temp2);

   //Low threshold T type calcium channel
   const double_t shift = 2.;
   const double_t Vm = V + shift;
   const double_t a_inf = 1.0 / ( 1. + std::exp(-(Vm+57.)/6.2) );
   const double_t b_inf = 1.0 / ( 1. + std::exp((Vm+81.)/4.0) );
   double_t tau_b = (30.8 + (211.4 + std::exp((Vm+113.2)/5.)) / (1. + std::exp((Vm+84.)/3.2)))/3.737;

   //High threshold T type calcium channel
   const double_t c_inf = 1.0 / ( 1. + std::exp(-(V+40.1)/3.5) );
   const double_t d_inf = 1.0 / ( 1. + std::exp((V+62.2)/5.5) );
   const double_t tau_d = 0.1483 * std::exp( -0.09398 * V) + 5.284 * std::exp(0.008855 * V);
   
   //H channel
   const double_t e_inf = 1. / (1. + std::exp((V + 60.)/5.5));
   const double_t tau_e = 20. + 1000./ ((std::exp((V + 56.5)/14.2)) + (std::exp(-(V+74.)/11.6)));

   //AHP channel
   const double_t z_inf = 48.* cai_HT * cai_HT / ( 48. * cai_HT * cai_HT + 0.09);
   const double_t tau_z = 1. / ( 48. * cai_HT * cai_HT + 0.09);

   //Channel currents
   const double_t I_Na =  node.P_.g_Na * m * m * m * h * (V - node.P_.E_Na);
   const double_t I_K  =  node.P_.g_K  * n * n * n * n * (V - node.P_.E_K );
   const double_t I_L  =  node.P_.g_L                  * (V - node.P_.E_L );
   const double_t I_KL =  node.P_.g_KL                 * (V - node.P_.E_K);
   const double_t I_LT =   node.P_.g_LT * a_inf * a_inf * b * (V - E_Ca_LT);
   const double_t I_HT =  node.P_.g_HT * c_inf * c_inf * d * (V - E_Ca_HT);
   const double_t I_H =   node.P_.g_H * e * (V - node.P_.E_H);
   const double_t I_AHP =  node.P_.g_AHP * z * z * (V - node.P_.E_K);

   //Calcium concentration for LT
   const double_t depth_LT = 1.;	
   const double_t taur_LT = 5.;
   const double_t cainf_LT = 2.4e-4;
   double_t drive_channel_LT =  - (10.) * I_LT / (2. * FARADAY * depth_LT);
   if (drive_channel_LT <= 0.) { drive_channel_LT = 0.; }

   //Calcium concentration for HT
   const double_t depth_HT = 1.;	
   const double_t taur_HT = 5.;
   const double_t cainf_HT = 2.4e-4;
   double_t drive_channel_HT =  - (10.) * I_HT / (2. * FARADAY * depth_HT);
   if (drive_channel_HT <= 0.) { drive_channel_HT = 0.; }

   // V dot -- synaptic input are currents, inhib current is negative
   f[S::V_M] = ( -(I_Na + I_K + I_HT + I_LT + I_H + I_KL + I_L + I_AHP) + node.B_.I_stim_ + node.P_.I_e + I_ex + I_in) / node.P_.C_m;

   //channel dynamics
   f[S::HT_M] = (m_inf - y[S::HT_M]) / tau_m; // m-variable
   f[S::HT_H] = (h_inf - y[S::HT_H]) / tau_h; // h-variable
   f[S::HT_N] = (n_inf - y[S::HT_N]) / tau_n; // n-variable
   f[S::HT_A] = 0.0;
   f[S::HT_B] = (b_inf - y[S::HT_B]) / tau_b; // n-variable
   f[S::HT_C] = 0.0;
   f[S::HT_D] = (d_inf - y[S::HT_D]) / tau_d; // n-variable
   f[S::HT_E] = (e_inf - y[S::HT_E]) / tau_e;
   f[S::HT_Z] = (z_inf - y[S::HT_Z]) / tau_z;
   f[S::HT_O2] = 0.0;//(o2_inf - y[S::HT_O2]) / tau_o2;
   f[S::HT_CAI_LT] = drive_channel_LT - (cai_LT - cainf_LT)/taur_LT; // calcium concentration
   f[S::HT_CAI_HT] = drive_channel_HT - (cai_HT - cainf_HT)/taur_HT; // calcium concentration

   // synapses: alpha functions
   f[S::DI_EXC] = -dI_ex / node.P_.tau_synE;
   f[S::I_EXC ] =  dI_ex  - (I_ex / node.P_.tau_synE);    
   f[S::DI_INH] = -dI_in / node.P_.tau_synI;
   f[S::I_INH ] =  dI_in  - (I_in / node.P_.tau_synI); 

   return GSL_SUCCESS;
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::ht_psc_alpha::Parameters_::Parameters_()
  : t_ref_  (    2.0  ),  // ms
    g_Na    (    90.  ),  // mS / cm2
    g_K     (    10.  ),  // mS / cm2
    g_L     (    0.01 ),  // mS / cm2
    g_KL    (    0.0069 ), // mS / cm2
    g_HT     (    6. ),  // mS / cm2
    g_LT     (    2. ),  // mS / cm2
    g_H     (    0.36 ),  // mS / cm2
    g_AHP   (    15.0 ),
    C_m     (    1.  ),  // uF
    E_Na    (   50.0  ),  // mV
    E_K     (  -100.0  ),  // mV
    E_H     (  -40.0  ),  // mV
    E_L     (  -70.   ),  // mV
    E_AHP     (  -100.   ),  // mV
    tau_synE(  0.2    ),  // ms
    tau_synI(  2.0    ),  // ms
    I_e     (  0.0    )   // pA

{
}

nest::ht_psc_alpha::State_::State_(const Parameters_&)
  : r_(0)
{
  y_[0] = -65.;//p.E_L;
  for ( size_t i = 1 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = 0.;
    
  //Calcium concentration
  const double_t cainf_LT = 2.4e-4;
  const double_t cainf_HT = 2.4e-4;
  
  //reversal pontential of Ca2+
  const double celsius	= 36.;
  const double R = 8.31441; //???
  const double_t FARADAY = 96485.3;
  const double cao=2.;
  const double E_Ca_LT = (1.e+3) * (R*(celsius+273.15))/(2.*FARADAY) * std::log(cao/cainf_LT);
  const double E_Ca_HT = (1.e+3) * (R*(celsius+273.15))/(2.*FARADAY) * std::log(cao/cainf_HT);
  
  // Na+ channel gate variable m
  const double_t V_traub = -25.;
  const double_t v2 = y_[0] - V_traub;
  double_t temp1 = 0.32 * vtrap(13.-v2, 4.);
  double_t temp2 = 0.28 * vtrap(v2-40.,5);
  const double_t tau_m = 1. / (temp1 + temp2);
  const double_t m_inf = temp1 / (temp1 + temp2);
  
  // Na+ channel gate variable h
  temp1 = 0.128 * std::exp((17.-v2)/18.);
  temp2 = 4. / ( 1. + std::exp((40.-v2)/5.) );
  const double_t tau_h = 1. / (temp1 + temp2);
  const double_t h_inf = temp1 / (temp1 + temp2);
  
  // K+ channel gate variable n
  temp1 = 0.032 * vtrap(15. - v2, 5.);
  temp2 = 0.5 * std::exp((10.-v2)/40.);
  const double_t tau_n = 1. / (temp1 + temp2);
  const double_t n_inf = temp1 / (temp1 + temp2);

  //Low-threshold T type calcium channel
  const double_t shift = 2.;
  const double_t Vm = y_[0]  + shift;
  const double_t a_inf = 1.0 / ( 1. + std::exp(-(Vm+57.)/6.2) );
  const double_t b_inf = 1.0 / ( 1. + std::exp((Vm+81.)/4.0) );
  double_t tau_b = 30.8 + (211.4 + std::exp((Vm+113.2)/5.)) / (1 + std::exp((Vm+84.)/3.2)) / 3.737;

  //High threshold T type calcium channel
  const double_t c_inf = 1.0 / ( 1. + std::exp(-(y_[0]+40.1)/3.5) );
  const double_t d_inf = 1.0 / ( 1. + std::exp((y_[0]+62.2)/5.5) );
  const double_t tau_d = 0.1483 * std::exp( -0.09398 * y_[0]) + 5.284 * std::exp(0.008855 * y_[0]);
  
  // H channel variable e
  const double_t e_inf = 1. / (1. + std::exp((y_[0] + 60.)/5.5));

  //AHP channel
  const double_t z_inf = 48.* cainf_HT * cainf_HT / ( 48. * cainf_HT * cainf_HT + 0.09);
   
  y_[HT_H] = h_inf;
  y_[HT_N] = n_inf;
  y_[HT_M] = m_inf;
  y_[HT_A] = a_inf;
  y_[HT_B] = b_inf;
  y_[HT_A] = c_inf;
  y_[HT_B] = d_inf;
  y_[HT_E] = e_inf;
  y_[HT_Z] = z_inf;
  y_[HT_O2] = 0.;
  y_[HT_P] = 0.;
  y_[HT_CAI_LT] = cainf_LT;
  y_[HT_CAI_HT] = cainf_HT;
}

nest::ht_psc_alpha::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
}

nest::ht_psc_alpha::State_& nest::ht_psc_alpha::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::ht_psc_alpha::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::t_ref, t_ref_);
  def<double>(d,names::g_Na, g_Na);
  def<double>(d,names::g_K, g_K);
  def<double>(d,names::g_HT, g_HT);
  def<double>(d,names::g_LT, g_LT);
  def<double>(d,names::g_H, g_H);
  def<double>(d,names::g_L, g_L);
  def<double>(d,names::g_KL, g_KL);
  def<double>(d,names::g_AHP, g_AHP);
  def<double>(d,names::E_Na, E_Na);
  def<double>(d,names::E_K, E_K);
  def<double>(d,names::E_H, E_H);
  def<double>(d,names::E_L, E_L);
  def<double>(d,names::E_AHP, E_AHP);
  def<double>(d,names::C_m, C_m);
  def<double>(d,names::tau_syn_ex, tau_synE);
  def<double>(d,names::tau_syn_in, tau_synI);
  def<double>(d,names::I_e, I_e);
}

void nest::ht_psc_alpha::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d,names::t_ref, t_ref_);
  updateValue<double>(d,names::C_m, C_m);
  updateValue<double>(d,names::g_Na,g_Na);
  updateValue<double>(d,names::E_Na,E_Na);
  updateValue<double>(d,names::g_K, g_K);
  updateValue<double>(d,names::g_HT, g_HT);
  updateValue<double>(d,names::g_LT, g_LT);
  updateValue<double>(d,names::E_K, E_K);
  updateValue<double>(d,names::g_H, g_H);
  updateValue<double>(d,names::E_H, E_H);
  updateValue<double>(d,names::g_L, g_L);
  updateValue<double>(d,names::E_L, E_L);
  updateValue<double>(d,names::g_KL, g_KL);
  updateValue<double>(d,names::g_AHP, g_AHP);
  updateValue<double>(d,names::E_AHP, E_AHP);
  updateValue<double>(d,names::tau_syn_ex,tau_synE);
  updateValue<double>(d,names::tau_syn_in,tau_synI);
  updateValue<double>(d,names::I_e, I_e);

  if ( C_m <= 0. )
    throw BadProperty("Capacitance must be strictly positive.");

  if ( t_ref_ < 0. )
    throw BadProperty("Refractory time cannot be negative.");

  if ( tau_synE <= 0. || tau_synI <= 0. )
    throw BadProperty("All time constants must be strictly positive.");
    
  if ( g_K < 0. || g_Na < 0. || g_L < 0. || g_H < 0. || g_HT < 0. || g_LT < 0. || g_KL < 0.)
    throw BadProperty("All conductances must be non-negative.");
}

void nest::ht_psc_alpha::State_::get(DictionaryDatum &d) const
{
  def<double>(d,names::V_m,  y_[V_M]);
  def<double>(d,names::HT_m, y_[HT_M]); 
  def<double>(d,names::HT_h, y_[HT_H]);
  def<double>(d,names::HT_n, y_[HT_N]);
  def<double>(d,names::HT_a, y_[HT_A]);
  def<double>(d,names::HT_b, y_[HT_B]);
  def<double>(d,names::HT_c, y_[HT_C]);
  def<double>(d,names::HT_d, y_[HT_D]);
  def<double>(d,names::HT_e, y_[HT_E]);
  def<double>(d,names::HT_z, y_[HT_Z]);
  def<double>(d,names::HT_o2, y_[HT_O2]);
  def<double>(d,names::HT_p, y_[HT_P]);
  def<double>(d,names::HT_cai_LT, y_[HT_CAI_LT]);
  def<double>(d,names::HT_cai_HT, y_[HT_CAI_HT]);
}

void nest::ht_psc_alpha::State_::set(const DictionaryDatum& d)
{
  updateValue<double>(d,names::V_m,  y_[V_M]);
  updateValue<double>(d,names::HT_m, y_[HT_M]); 
  updateValue<double>(d,names::HT_h, y_[HT_H]);
  updateValue<double>(d,names::HT_n, y_[HT_N]);
  updateValue<double>(d,names::HT_a, y_[HT_A]);
  updateValue<double>(d,names::HT_b, y_[HT_B]);
  updateValue<double>(d,names::HT_c, y_[HT_C]);
  updateValue<double>(d,names::HT_d, y_[HT_D]);
  updateValue<double>(d,names::HT_e, y_[HT_E]);
  updateValue<double>(d,names::HT_z, y_[HT_Z]);
  updateValue<double>(d,names::HT_o2, y_[HT_O2]);
  updateValue<double>(d,names::HT_p, y_[HT_P]);
  updateValue<double>(d,names::HT_cai_LT, y_[HT_CAI_LT]);
  updateValue<double>(d,names::HT_cai_HT, y_[HT_CAI_HT]);
   
  if ( y_[HT_M] < 0. || y_[HT_H] < 0. || y_[HT_N] < 0. || y_[HT_A] < 0. || y_[HT_B] < 0. || y_[HT_C] < 0. || y_[HT_D] < 0. || y_[HT_E] < 0. || y_[HT_Z] < 0. || y_[HT_CAI_LT] < 0. || y_[HT_CAI_HT] < 0. )
    throw BadProperty("All (in)activation variables must be non-negative.");
}

nest::ht_psc_alpha::Buffers_::Buffers_(ht_psc_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
    // Initialization of the remaining members is deferred to
    // init_buffers_().
}

nest::ht_psc_alpha::Buffers_::Buffers_(const Buffers_&, ht_psc_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
    // Initialization of the remaining members is deferred to
    // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::ht_psc_alpha::ht_psc_alpha()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

nest::ht_psc_alpha::ht_psc_alpha(const ht_psc_alpha& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

nest::ht_psc_alpha::~ht_psc_alpha()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::ht_psc_alpha::init_state_(const Node& proto)
{
  const ht_psc_alpha& pr = downcast<ht_psc_alpha>(proto);
  S_ = pr.S_;
}

void nest::ht_psc_alpha::init_buffers_()
{
  B_.spike_exc_.clear();       // includes resize
  B_.spike_inh_.clear();       // includes resize
  B_.currents_.clear();        // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = ht_psc_alpha_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

  B_.I_stim_ = 0.0;
}

void nest::ht_psc_alpha::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate

  V_.PSCurrInit_E_  = 1.0 * numerics::e / P_.tau_synE;
  V_.PSCurrInit_I_  = 1.0 * numerics::e / P_.tau_synI;
  V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::ht_psc_alpha::update(Time const & origin, const long_t from, const long_t to)
{
   
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    
    double_t       t = 0.0 ;
    const double_t U_old = S_.y_[State_::V_M];

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
			   &B_.sys_,             // system of ODE
			   &t,                   // from t
			    B_.step_,            // to t <= step
			   &B_.IntegrationStep_, // integration step size
			    S_.y_);              // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure(get_name(), status);
    }
    
    S_.y_[State_::DI_EXC] += B_.spike_exc_.get_value(lag) * V_.PSCurrInit_E_;
    S_.y_[State_::DI_INH] += B_.spike_inh_.get_value(lag) * V_.PSCurrInit_I_;

    // sending spikes: crossing 0 mV, pseudo-refractoriness and local maximum...
    // refractory?
    if ( S_.r_ > 0 )
	    --S_.r_;
    else
      // (    threshold    &&     maximum       )
      if ( S_.y_[State_::V_M] >= 0 && U_old > S_.y_[State_::V_M])
	    {
	      S_.r_ = V_.RefractoryCounts_;
	  
	      set_spiketime(Time::step(origin.get_steps()+lag+1));
	  
	      SpikeEvent se;
	      network()->send(*this, se, lag);
	    }
    
    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);

    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);
  
  }
}

void nest::ht_psc_alpha::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_weight() > 0.0)
    B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );  // current input, keep negative weight
}

void nest::ht_psc_alpha::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		                     w *c);
}

void nest::ht_psc_alpha::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL
