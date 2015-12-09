/*
 *  izhikevich_cond.cpp
 *
 *  Created by Osamu Shouno	2009/02/04
 *  Copyright 2011 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *  Modified by Makoto Otsuka   2011/11/16
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */


#include "izhikevich_cond.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <iomanip>
#include <iostream>
#include <cstdio>

nest::RecordablesMap<nest::izhikevich_cond> nest::izhikevich_cond::recordablesMap_;

namespace nest {
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<izhikevich_cond>::create()
  {
    // use standard names whereever you can for consistency!
    insert_("V",
            &izhikevich_cond::get_y_elem_<izhikevich_cond::State_::V_M>);
    insert_("u",
            &izhikevich_cond::get_y_elem_<izhikevich_cond::State_::U_M>);
  }



extern "C"
inline int izhikevich_cond_dynamics (double, const double y[], double f[], void* pnode)
{ 

  // a shorthand
  typedef izhikevich_cond::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const izhikevich_cond& node =  *(reinterpret_cast<izhikevich_cond*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...


   double I_syn_exc = y[3] * (y[0] - node.P_.V_E);	// y[3] is excitatory synapse conductance
   double I_syn_inh = y[5] * (y[0] - node.P_.V_I);	// y[5] is inhibitory synapse conductance

   //V dot
   
   f[0] = ( 
			+ node.P_.k_*(y[0] - node.P_.Vr_)*(y[0] - node.P_.Vt_)
			- y[1]
			+ node.B_.I_stim_ 
			+ node.P_.I_e 
			- I_syn_exc 
			- I_syn_inh
		  ) / node.P_.Cm_;
   
   f[1] = node.P_.a_*( node.P_.b_*(y[0] - node.P_.Vr_) - y[1]);

   f[2] = -y[2] / node.P_.tau_synE;
   f[3] =  y[2] - (y[3]/node.P_.tau_synE); // Synaptic Conductance (nS)

   f[4] = -y[4] / node.P_.tau_synI;
   f[5] =  y[4] - (y[5]/node.P_.tau_synI); // Synaptic Conductance (nS)

   return GSL_SUCCESS;
 }

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::izhikevich_cond::Parameters_::Parameters_()
  : 

	Vr_			( -80.0	  ),  // mV
	Vt_			( -25.0	  ),  // mV
	Vpeak_		(  40.0	  ),  // mV

	Cm_			(  50.0	  ),  // pF/micro-meter^2

   k_			(	1.0	  ),
	a_			(	0.01  ),
	b_			( -20.0	  ),
	c_			( -55.0	  ),  // mV
	d_			( 150.0	  ),  

    tau_synE	(   1.0    ),  // ms
    tau_synI	(   1.0    ),  // ms
	V_E			(   0.0	  ),  // mV
	V_I			(-100.0	  ),  // mV
    I_e			(   0.0    )   // pA
{}

nest::izhikevich_cond::State_::State_(const Parameters_& p)
  : r_(0)
{
  y_[0] = p.Vr_;
  for ( size_t i = 1 ; i < 6 ; ++i )
    y_[i] = 0;
}

nest::izhikevich_cond::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < 6 ; ++i )
    y_[i] = s.y_[i];
}

nest::izhikevich_cond::State_& nest::izhikevich_cond::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < 6 ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::izhikevich_cond::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,"k",			  k_);
  def<double>(d,"a",			  a_);
  def<double>(d,"b",			  b_);
  def<double>(d,"c",			  c_);
  def<double>(d,"d",			  d_);
  def<double>(d,"Cm",			  Cm_);
  def<double>(d,"Vr",			  Vr_);
  def<double>(d,"Vt",			  Vt_);
  def<double>(d,"Vpeak",		  Vpeak_);
  def<double>(d,"V_E",			  V_E);
  def<double>(d,"V_I",			  V_I);
  def<double>(d,"tau_syn_ampa",  tau_synE);
  def<double>(d,"tau_syn_gaba",  tau_synI);
  def<double>(d,"I0",			  I_e);
}

void nest::izhikevich_cond::Parameters_::set(const DictionaryDatum& d)
{
  // allow setting the membrane potential
  updateValue<double>(d,"k",	  k_);
  updateValue<double>(d,"a",	  a_);
  updateValue<double>(d,"b",	  b_);
  updateValue<double>(d,"c",	  c_);
  updateValue<double>(d,"d",	  d_);
  updateValue<double>(d,"Cm",    Cm_);
  updateValue<double>(d,"Vr",    Vr_);
  updateValue<double>(d,"Vt",    Vt_);
  updateValue<double>(d,"Vpeak", Vpeak_);
  updateValue<double>(d,"V_E",   V_E);
  updateValue<double>(d,"V_I",   V_I);
  updateValue<double>(d,"tau_syn_ampa", tau_synE);
  updateValue<double>(d,"tau_syn_gaba", tau_synI);
  updateValue<double>(d,"I0",     I_e);

  if ( Vr_ >= Vt_ )
    throw nest::BadProperty("Reset potential must be smaller than threshold.");
    
  if ( Cm_ <= 0 )
    throw nest::BadProperty("Capacitance must be strictly positive.");
      
  if ( tau_synE <= 0 || tau_synI <= 0 )
    throw nest::BadProperty("All time constants must be strictly positive.");
}

void nest::izhikevich_cond::State_::get(DictionaryDatum &d) const
{
  def<double>(d, "Vm", y_[0]);	  // Membrane potential
  def<double>(d, "Um", y_[1]);	  // Membrane potential recovery variable
}

void nest::izhikevich_cond::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  updateValue<double>(d, "Vm", y_[0]);
  updateValue<double>(d, "Um", y_[1]);
}

nest::izhikevich_cond::Buffers_::Buffers_(izhikevich_cond& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
    // Initialization of the remaining members is deferred to init_buffers_().
}

nest::izhikevich_cond::Buffers_::Buffers_(const Buffers_&, izhikevich_cond& n)
    : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to init_buffers_().
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::izhikevich_cond::izhikevich_cond()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

nest::izhikevich_cond::izhikevich_cond(const izhikevich_cond& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

nest::izhikevich_cond::~izhikevich_cond()
{
  // GSL structs only allocated by init_nodes_(), so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */
/**
void nest::izhikevich_cond::init_node_(const Node& proto)
{
  const izhikevich_cond& pr = downcast<izhikevich_cond>(proto);
  P_ = pr.P_;
  S_ = pr.S_;
}
**/
void nest::izhikevich_cond::init_state_(const Node& proto)
{
  const izhikevich_cond& pr = downcast<izhikevich_cond>(proto);
  S_ = pr.S_;
}

void nest::izhikevich_cond::init_buffers_()
{
  B_.spike_exc_.clear();       // includes resize
  B_.spike_inh_.clear();       // includes resize
  
  B_.currents_.clear();        // includes resize
//  B_.potentials_.clear_data(); // includes resize
//  B_.conductances_.clear_data(); // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = nest::Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, 6);				//	6-dim system
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);		// abs, rel error
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(6);					//	6-dim system
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = izhikevich_cond_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension =  6;
  B_.sys_.params    = reinterpret_cast<void*>(this);
//  B_.sys_.params    = reinterpret_cast<void*>(&P_);

  B_.I_stim_ = 0.0;
}

void nest::izhikevich_cond::calibrate()
{
  B_.logger_.init();    // ensure initialization in case mm connected after Simulate

  V_.PSConInit_E_  = 1.0 * numerics::e / P_.tau_synE;
  V_.PSConInit_I_  = 1.0 * numerics::e / P_.tau_synI;
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::izhikevich_cond::update(nest::Time const & origin, const nest::long_t from, const nest::long_t to)
{
   
  assert(to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay());
  assert(from < to);

  for ( nest::long_t lag = from ; lag < to ; ++lag )
  {
    
    double t = 0.0;

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
			    S_.y_); 	         // neuronal state

      if ( status != GSL_SUCCESS )
        throw nest::GSLSolverFailure(get_name(), status);
    }

    S_.y_[2] += B_.spike_exc_.get_value(lag) * V_.PSConInit_E_;
    S_.y_[4] += B_.spike_inh_.get_value(lag) * V_.PSConInit_I_;
	
	if ( S_.y_[0] < P_.Vpeak_)
	{
		S_.r_ = 0;
	 } 
    else
	{
      // neuron is not absolute refractory
      if ( S_.y_[0] >= P_.Vpeak_ && S_.r_ == 0)
	  {
	      S_.y_[0] = P_.Vpeak_;
		  S_.r_	   = 1;
  
	      nest::SpikeEvent se;
	      network()->send(*this, se, lag);
	  }
	  else if (S_.y_[0] >= P_.Vpeak_ && S_.r_ == 1)
	  {
		  S_.y_[0]  = P_.c_;	// reset membrane potential after generation of a spike
		  S_.y_[1] += P_.d_;	// reset membrane potential recovery current after generation of a spike
	  }
    }

    //log state data
    B_.logger_.record_data(origin.get_steps() + lag);
	
    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);
   
    // voltage logging
//    B_.potentials_.record_data(origin.get_steps()+lag, S_.y_[0] );
//    B_.conductances_.record_data(origin.get_steps()+lag, 
//                              std::pair<double_t, double_t>(S_.y_[3], S_.y_[5]));
  }
}

void nest::izhikevich_cond::handle(nest::SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_rport() == 1)
    B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );
  
  if(e.get_rport() == 3)
    B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );  
}

void nest::izhikevich_cond::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const nest::double_t c=e.get_current();
  const nest::double_t w=e.get_weight();

  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w *c);
}

void nest::izhikevich_cond::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

} // namespace nest
#endif //HAVE_GSL
