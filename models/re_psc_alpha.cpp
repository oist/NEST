/*
 *  re_psc_alpha.cpp
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


#include "re_psc_alpha.h"

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


nest::RecordablesMap<nest::re_psc_alpha> nest::re_psc_alpha::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<re_psc_alpha>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::V_M>);
    insert_(names::I_ex, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::I_EXC>);
    insert_(names::I_in, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::I_INH>);
    insert_(names::RE_m, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::RE_M>);
    insert_(names::RE_h, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::RE_H>);
    insert_(names::RE_n, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::RE_N>);
    insert_(names::RE_a, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::RE_A>);
    insert_(names::RE_b, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::RE_B>);
    insert_(names::RE_cai, 
	    &re_psc_alpha::get_y_elem_<re_psc_alpha::State_::RE_CAI>);

  }

  extern "C"
  int re_psc_alpha_dynamics (double, const double y[], double f[], void* pnode)
  { 
    // a shorthand
    typedef nest::re_psc_alpha::State_ S;
  
    // get access to node so we can almost work as in a member function
    assert(pnode);
    const nest::re_psc_alpha& node =  *(reinterpret_cast<nest::re_psc_alpha*>(pnode));
    
    // y[] here is---and must be---the state vector supplied by the integrator,
    // not the state vector in the node, node.S_.y[]. 
    
    // The following code is verbose for the sake of clarity. We assume that a
    // good compiler will optimize the verbosity away ...
    
   // shorthand for state variables
   const double_t& V     = y[S::V_M   ];
   const double_t& m     = y[S::RE_M  ];
   const double_t& h     = y[S::RE_H  ];
   const double_t& n     = y[S::RE_N  ];
   const double_t& a     = y[S::RE_A  ];
   const double_t& b     = y[S::RE_B  ];
   const double_t& cai   = y[S::RE_CAI];
   const double_t& dI_ex = y[S::DI_EXC];
   const double_t&  I_ex = y[S::I_EXC ];
   const double_t& dI_in = y[S::DI_INH];
   const double_t&  I_in = y[S::I_INH ];

   //reversal pontential of Ca2+
   const double celsius	= 36.;
   const double R = 8.31441;
   const double_t FARADAY = 96489.;
   const double cao=2.;
   const double E_Ca = (1.e+3) * (R*(celsius+273.15))/(2.*FARADAY) * std::log(cao/cai);

   // Na+ channel gate variable m
   const double_t V_traub = -63.;
   const double_t v2 = V - V_traub;
   const double_t tadj_hh = std::pow(3., (celsius-36.)/ 10.);
   double_t temp1 = 0.32 * vtrap(13.-v2, 4.);
   double_t temp2 = 0.28 * (v2-40.) / ( std::exp((v2-40.)/5.) - 1.);

   const double_t tau_m = 1. / (temp1 + temp2) / tadj_hh;
   const double_t m_inf = temp1 / (temp1 + temp2);

   // Na+ channel gate variable h
   temp1 = 0.128 * std::exp((17.-v2)/18.);
   temp2 = 4. / ( 1. + std::exp((40.-v2)/5.) );
   const double_t tau_h = 1. / (temp1 + temp2) / tadj_hh;
   const double_t h_inf = temp1 / (temp1 + temp2);

   // K+ channel gate variable n
   temp1 = 0.032 * (15.-v2) / ( std::exp((15.-v2)/5.) - 1.);
   temp2 = 0.5 * std::exp((10.-v2)/40.);
   const double_t tau_n = 1. / (temp1 + temp2) / tadj_hh;
   const double_t n_inf = temp1 / (temp1 + temp2);

   // Low threshold calcium gate variable a, b
   const double_t shift = 0.0;
   const double_t phi_a = std::pow(2.5, (celsius-24.)/10.);
   const double_t phi_b = std::pow(2.5, (celsius-24.)/10.);
   const double_t a_inf = 1.0 / ( 1. + std::exp(-(V+shift+50.)/7.4) );
   const double_t b_inf = 1.0 / ( 1. + std::exp((V+shift+78.)/5.0) );
   const double_t tau_a = ( 3. + 1.0 / ( std::exp((V+shift+25.)/10.) + exp(-(V+shift+100.)/15.) ) ) / phi_a;
   const double_t tau_b = ( 85. + 1.0 / ( std::exp((V+shift+46.)/4.) + exp(-(V+shift+405.)/50.) ) ) / phi_b;

   //Calciumu concentration
   const double_t depth = 1.;	
   const double_t taur	= 1.e+10;
   const double_t cainf = 2.4e-4;
   const double_t kt = 1.e-4;
   const double_t kd = 1.e-4;
   double drive_channel =  - (10.) * (node.P_.g_Ca * a * a * b * (V - E_Ca )) / (2. * FARADAY * depth);
   if (drive_channel <= 0.) { drive_channel = 0.; }	
   const double drive_pump = -kt * cai / (cai + kd );	

   //Currents
   const double_t I_Na =  node.P_.g_Na * m * m * m * h * (V - node.P_.E_Na);
   const double_t I_K  =  node.P_.g_K  * n * n * n * n * (V - node.P_.E_K );
   const double_t I_Ca =  node.P_.g_Ca * a * a * b     * (V - E_Ca);
   const double_t I_L  =  node.P_.g_L   * (V - node.P_.E_L );

   // V dot -- synaptic input are currents, inhib current is negative
   f[S::V_M] = ( -(I_Na + I_K + I_Ca + I_L) + node.B_.I_stim_ + node.P_.I_e + I_ex + I_in) / node.P_.C_m;

   //channel dynamics
   f[S::RE_M] = (m_inf - y[S::RE_M]) / tau_m; // m-variable
   f[S::RE_H] = (h_inf - y[S::RE_H]) / tau_h; // h-variable
   f[S::RE_N] = (n_inf - y[S::RE_N]) / tau_n; // n-variable
   f[S::RE_A] = (a_inf - y[S::RE_A]) / tau_a; // a-variable
   f[S::RE_B] = (b_inf - y[S::RE_B]) / tau_b; // b-variable
   f[S::RE_CAI] = drive_channel + drive_pump + (cainf - y[S::RE_CAI]) / taur;
  
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
    
nest::re_psc_alpha::Parameters_::Parameters_()
  : t_ref_  (    2.0  ),  // ms
    g_Na    ( 200.  ),  // mS / cm2
    g_K     ( 20.  ),  // mS / cm2
    g_K_AHP ( 0.  ),  // mS / cm2
    g_Ca    ( 3. ),  // mS / cm2
    g_CAN    ( 0. ),  // mS / cm2
    g_L     (   0.05  ),  // mS / cm2
    C_m     (  1.  ),  // uF/ cm2
    E_Na    (   50.0  ),  // mV
    E_K     (  -90.0  ),  // mV
    E_CAN    (  -20.0  ),  // mV
    E_L     (  -90.),  // mV
    tau_synE(  0.2    ),  // ms
    tau_synI(  2.0    ),  // ms
    I_e     (  0.0    )   // pA


{
}

nest::re_psc_alpha::State_::State_(const Parameters_&)
  : r_(0)
{
  y_[0] = -65;//p.E_L;

  for ( size_t i = 1 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = 0;
    
  //reversal pontential of Ca2+
  const double celsius	= 36;
  const double R = 8.31441;
  const double_t FARADAY = 96489.;
  const double cao=2.;
  const double E_Ca = (1.e+3) * (R*(celsius+273.15))/(2.*FARADAY) * std::log(cao/y_[RE_CAI]);

  const double_t V_traub = -55.;
  const double_t v2 = y_[0] - V_traub;
  double_t temp1 =  0.32 * (13.-v2) / ( std::exp((13.-v2)/4.) - 1.) ;
  double_t temp2 = 0.28 * (v2-40.) / ( std::exp((v2-40.)/5.) - 1.);
  const double_t tau_m = 1. / (temp1 + temp2);
  const double_t m_inf = temp1 / (temp1 + temp2);
  
  temp1 = 0.128 * std::exp((17.-v2)/18.);
  temp2 = 4. / ( 1. + std::exp((40.-v2)/5.) );
  const double_t tau_h = 1. / (temp1 + temp2);
  const double_t h_inf = temp1 / (temp1 + temp2);
  
  temp1 = 0.032 * (15.-v2) / ( std::exp((15.-v2)/5.) - 1.);
  temp2 = 0.5 * std::exp((10.-v2)/40.);
  const double_t tau_n = 1. / (temp1 + temp2);
  const double_t n_inf = temp1 / (temp1 + temp2);

   // Low threshold calcium gate variable a, b
  const double_t shift = 2.0;
  const double_t phi_a = std::pow(5., (celsius-24.)/10.);
  const double_t phi_b = std::pow(3., (celsius-24.)/10.);
  const double_t a_inf = 1.0 / ( 1. + std::exp(-(y_[0]+shift+50.)/7.4) );
  const double_t b_inf = 1.0 / ( 1. + std::exp((y_[0]+shift+78.)/5.0) );

  // Calcium concentration
  const double_t depth = 1.;	
  const double_t taur	= 1.e+10;
  const double_t cainf = 2.4e-4;

  y_[RE_M] = m_inf; // m-variable
  y_[RE_H] = h_inf; // h-variable
  y_[RE_N] = n_inf; // n-variable
  y_[RE_A] = a_inf; // a-variable
  y_[RE_B] = b_inf; // b-variable
  y_[RE_CAI] = cainf; // Calcium concentration
}

nest::re_psc_alpha::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
}

nest::re_psc_alpha::State_& nest::re_psc_alpha::State_::operator=(const State_& s)
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

void nest::re_psc_alpha::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::t_ref, t_ref_);
  def<double>(d,names::g_Na, g_Na);
  def<double>(d,names::g_K, g_K);
  def<double>(d,names::g_L, g_L);
  def<double>(d,names::E_Na, E_Na);
  def<double>(d,names::E_K, E_K);
  def<double>(d,names::E_L, E_L);
  def<double>(d,names::C_m, C_m);
  def<double>(d,names::tau_syn_ex, tau_synE);
  def<double>(d,names::tau_syn_in, tau_synI);
  def<double>(d,names::I_e, I_e);
}

void nest::re_psc_alpha::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d,names::t_ref, t_ref_);
  updateValue<double>(d,names::C_m, C_m);
  updateValue<double>(d,names::g_Na,g_Na);
  updateValue<double>(d,names::E_Na,E_Na);
  updateValue<double>(d,names::g_K, g_K);
  updateValue<double>(d,names::E_K, E_K);
  updateValue<double>(d,names::g_L, g_L);
  updateValue<double>(d,names::E_L, E_L);

  updateValue<double>(d,names::tau_syn_ex,tau_synE);
  updateValue<double>(d,names::tau_syn_in,tau_synI);

  updateValue<double>(d,names::I_e, I_e);

  if ( C_m <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");

  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time cannot be negative.");

  if ( tau_synE <= 0 || tau_synI <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
    
  if ( g_K < 0 || g_Na < 0 || g_L < 0 )
    throw BadProperty("All conductances must be non-negative.");
}

void nest::re_psc_alpha::State_::get(DictionaryDatum &d) const
{
  def<double>(d,names::V_m    , y_[V_M]);
  def<double>(d,names::RE_m  , y_[RE_M]); 
  def<double>(d,names::RE_h  , y_[RE_H]);
  def<double>(d,names::RE_n, y_[RE_N]);
  def<double>(d,names::RE_a  , y_[RE_A]);
  def<double>(d,names::RE_b, y_[RE_B]);
  def<double>(d,names::RE_cai, y_[RE_CAI]);
}

void nest::re_psc_alpha::State_::set(const DictionaryDatum& d)
{
  updateValue<double>(d,names::V_m    , y_[V_M]);
  updateValue<double>(d,names::RE_m  , y_[RE_M]); 
  updateValue<double>(d,names::RE_h  , y_[RE_H]);
  updateValue<double>(d,names::RE_n, y_[RE_N]);
  updateValue<double>(d,names::RE_a,   y_[RE_A]);
  updateValue<double>(d,names::RE_b, y_[RE_B]);
  updateValue<double>(d,names::RE_cai,   y_[RE_CAI]);
    
  if ( y_[RE_M] < 0 || y_[RE_H] < 0 || y_[RE_N] < 0 || y_[RE_A] < 0 || y_[RE_B] < 0 || y_[RE_CAI] < 0 )
    throw BadProperty("All (in)activation variables must be non-negative.");
}

nest::re_psc_alpha::Buffers_::Buffers_(re_psc_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
    // Initialization of the remaining members is deferred to
    // init_buffers_().
}

nest::re_psc_alpha::Buffers_::Buffers_(const Buffers_&, re_psc_alpha& n)
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

nest::re_psc_alpha::re_psc_alpha()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

nest::re_psc_alpha::re_psc_alpha(const re_psc_alpha& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

nest::re_psc_alpha::~re_psc_alpha()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::re_psc_alpha::init_state_(const Node& proto)
{
  const re_psc_alpha& pr = downcast<re_psc_alpha>(proto);
  S_ = pr.S_;
}

void nest::re_psc_alpha::init_buffers_()
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
  
  B_.sys_.function  = re_psc_alpha_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

  B_.I_stim_ = 0.0;
}

void nest::re_psc_alpha::calibrate()
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

void nest::re_psc_alpha::update(Time const & origin, const long_t from, const long_t to)
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

void nest::re_psc_alpha::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_weight() > 0.0)
    B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );  // current input, keep negative weight
}

void nest::re_psc_alpha::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		                     w *c);
}

void nest::re_psc_alpha::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


#endif //HAVE_GSL