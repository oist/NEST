/*
 *  stn_cond_neuron6.h
 *
 *
 *  Created by Osamu Shouno	2009/02/09
 *  Copyright 2011 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *  Modified by Makoto Otsuka 2012/03/15
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef STN_COND_NEURON6_H
#define STN_COND_NEURON6_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
//#include <gsl/gsl_odeiv2.h>

namespace nest
{

  /**
    * Function computing right-hand side of ODE for GSL solver.
    * @note Must be declared here so we can befriend it in class.
    * @note Must have C-linkage for passing to GSL. Internally, it is
    *       a first-class C++ function, but cannot be a member function
    *       because of the C-linkage.
    * @note No point in declaring it inline, since it is called
    *       through a function pointer.
    * @param void* Pointer to model neuron instance.
    */
  extern "C"
	int stn_cond_neuron6_dynamics(double, const double*, double*, void*);

//  extern "C"
//	int stn_cond_neuron6_nmda_dynamics(double, const double*, double*, void*);

/* BeginDocumentation

*/

  class stn_cond_neuron6: public  Archiving_Node
  {
    
  public:
    
    stn_cond_neuron6();
    stn_cond_neuron6(const stn_cond_neuron6&);
    ~stn_cond_neuron6();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */
#ifndef IS_BLUEGENE
    using  Node::check_connection;
#endif
    using  Node::connect_sender;
    using  Node::handle;

	port check_connection(Connection&,  port);
    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &);
    
	port connect_sender(SpikeEvent &,  port);
	port connect_sender(CurrentEvent &, port);
    port connect_sender(DataLoggingRequest &, port);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  private:
//    void init_node_(const Node& proto);
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void update(Time const &, const long_t, const long_t);

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------

    // make dynamics function quasi-member
    friend int stn_cond_neuron6_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<stn_cond_neuron6>;
    friend class UniversalDataLogger<stn_cond_neuron6>;

  private:

    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
      */
    struct Parameters_ {
	  double_t V_th_;       //!< threshold [mV]
	  double_t C_m;            //!<membrane capacitance [micro-F/cm2]		// 2005.7.4
	  
	  double_t V_L;   //!<leak-current reversal potential [mV]
	  double_t V_K;   //!< potassium current reversal potential (K+) [mV]  [mV]
	  double_t V_Na;   //!< sodium current reversal potential (Na+) [mV]
	  double_t V_Ca;   //!< calcium current reversal potential (Ca++) [mV] 
    
	  double_t g_L;         //!< leak conductance [mS/cm2]
	  double_t g_K;         //!< potassium conductance [mS/cm2]
	  double_t g_Na;        //!< sodium conductance [mS/cm2]
	  double_t g_IA;		//!< A-current conductance [mS/cm2]
	  double_t g_LCa;       //!< L-type calcium conductance [mS/cm2]
	  double_t g_T;         //!< T-current conductance [mS/cm2]
	  double_t g_CaK;		//!< AHP-current conductance [mS/cm2]
	
	  double_t tau_n0;   //!< time constant of K-activation [ms]
	  double_t tau_n1;   //!< time constant of K-activation [ms]
	  double_t tau_m0;   //!< time constant of Na-activation [ms]
	  double_t tau_m1;   //!< time constant of Na-activation [ms]
	  double_t tau_h0;   //!< time constant of Na-inactivation [ms]
	  double_t tau_h1;   //!< time constant of Na-inactivation [ms]
	  double_t tau_b0;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_b1;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_c0;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_c1;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_d0;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_d1;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_p0;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_p1;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_q0;   //!< time constant of CaK-inactivation [ms]
	  double_t tau_q1;   //!< time constant of CaK-inactivation [ms]

	  double_t tau_a0;   //!< time constant of K-activation [ms]
	  double_t tau_a1;   //!< time constant of K-activation [ms]

	  double_t theta_n;   //!< half activation voltage for n [mV]
	  double_t theta_m;   //!< half activation voltage for m [mV]
	  double_t theta_h;   //!< half inactivation voltage for h [mV]
	  double_t theta_a;   //!< half activation voltage for a [mV]
	  double_t theta_b;   //!< half activation voltage for b [mV]
	  double_t theta_c;   //!< half activation voltage for c [mV]
	  double_t theta_d1;   //!< half activation voltage for d1 [mV]
	  double_t theta_p;   //!< half activation voltage for p [mV]
	  double_t theta_q;   //!< half activation voltage for q [mV]
	
	  double_t sigma_n;   //!< slope factor for n
	  double_t sigma_m;   //!< slope factor for m
	  double_t sigma_h;   //!< slope factor for h
	  double_t sigma_a;   //!< slope factor for a
	  double_t sigma_b;   //!< slope factor for b
	  double_t sigma_c;   //!< slope factor for c
	  double_t sigma_d1;   //!< slope factor for d1
	  double_t sigma_p;   //!< slope factor for p
	  double_t sigma_q;   //!< slope factor for q
	  
      double_t sigma_d2;   //!< slope factor for d2

	  double_t theta_tau_m;    //!< 
	  double_t theta_tau_a;    //!< 
	  double_t theta_tau_n0;   //!< 
	  double_t theta_tau_n1;   //!< 
	  double_t theta_tau_h0;   //!< 
	  double_t theta_tau_h1;   //!< 
	  double_t theta_tau_b0;   //!< 
	  double_t theta_tau_b1;   //!<
	  double_t theta_tau_c0;   //!< 
	  double_t theta_tau_c1;   //!< 
	  double_t theta_tau_d0;   //!< 
	  double_t theta_tau_d1;   //!<
	  double_t theta_tau_p0;   //!< 
	  double_t theta_tau_p1;   //!<
	  double_t theta_tau_q0;   //!< 
	  double_t theta_tau_q1;   //!<
	
	  double_t sigma_tau_m;    //!< slope factor for the voltage dependence of the time constant tau_m
	  double_t sigma_tau_a;    //!< slope factor for the voltage dependence of the time constant tau_a
	  double_t sigma_tau_n0;   //!< slope factor for the voltage dependence of the time constant tau_n
	  double_t sigma_tau_n1;   //!< slope factor for the voltage dependence of the time constant tau_n
	  double_t sigma_tau_h0;   //!< slope factor for the voltage dependence of the time constant tau_h
	  double_t sigma_tau_h1;   //!< slope factor for the voltage dependence of the time constant tau_h
	  double_t sigma_tau_b0;   //!< slope factor for the voltage dependence of the time constant tau_b
	  double_t sigma_tau_b1;   //!< slope factor for the voltage dependence of the time constant tau_b
	  double_t sigma_tau_c0;   //!< slope factor for the voltage dependence of the time constant tau_c
	  double_t sigma_tau_c1;   //!< slope factor for the voltage dependence of the time constant tau_c
	  double_t sigma_tau_d0;   //!< slope factor for the voltage dependence of the time constant tau_d1
	  double_t sigma_tau_d1;   //!< slope factor for the voltage dependence of the time constant tau_d1
	  double_t sigma_tau_p0;   //!< slope factor for the voltage dependence of the time constant tau_p
	  double_t sigma_tau_p1;   //!< slope factor for the voltage dependence of the time constant tau_p
	  double_t sigma_tau_q0;   //!< slope factor for the voltage dependence of the time constant tau_q
	  double_t sigma_tau_q1;   //!< slope factor for the voltage dependence of the time constant tau_q
	
//	  double_t k1;		//!< dissociation constant of calcium-dependent AHP current
	  double_t k1_r;		//!< dissociation constant of calcium-dependent AHP current
	  double_t k1_d2;		//!< dissociation constant of calcium-dependent AHP current
	  double_t k_Ca;		//!< calcium pump rate constant [(coulombs-liter)/(moles-sec)]
	  double_t epsilon;	//!< combines the effects of buffers, cell volume, the molar charge of calcium [(moles-sec)/(coulombs-liter)]
	
	  double_t tau_syn_ampa;		//!<Time constant, excitatory synapse (AMPA) [mV]
	  //double_t tau_syn_gabaa;		//!<Time constant, inhibitory synapse (GABA-A) [mV]
	  double_t tau_syn_gabaa_rise;
	  double_t tau_syn_gabaa_decay;
	  double_t tau_syn_gabaa1;
	  double_t tau_syn_gabaa2;
	
	  double_t tau_syn_gabab;		//!<Time constant, inhibitory synapse (GABA-B) [mV]
	  
	  double_t ampa_conductance;
	  double_t gabaa_conductance;
	  double_t gabab_conductance;
	  
	  double_t nmda_conductance;
	  
	  double_t V_E;		//!<Reversal potential, excitatory synapse [mV]
	  double_t V_GABAA;	//!<Reversal potential, inhibitory synapse [mV]
	  double_t V_GABAB;	//!<Reversal potential, inhibitory synapse [mV]
	  
	  double_t I0;		//!<Persistent DC input current [micro-F/cm2]  ---> [pA, nS]
		
	  double_t area;    // scale parameter for external input
						//	In original model, units of currents and cunductances are micro-A/cm2 and mS/cm2.
						//	However, comparison between this model and real experiment needs to have same units such as pA or nS.
						//	According to the paper where 50 pA and 5 micro-A/cm2 are comparable, we choose 1.0 * 1e-1
						//			pA, nS x 1.0 *1e-1 ---> micro-A/cm2, mS/cm2
						//	More precisely, this value should be 1e-5 cm^2 = 1e-3 mm^2 ~= 4*PI*(0.01 mm)^2 =0.0012 
						//  However, just use this parameter just for conversion
		
		
	  double MgConc;
	  double tau_nmda_decay;
	  double tau_nmda_rise;
	  double alpha_nmda;
	  double g_nmda;
		
	  double g_ampa;

	  /** variables and parameters for dopaminergic modulation **/
	  double_t	da2_Ca_modulation;	// 0, no modulation ; 1.0, modulation
    
	  double_t	da2_Ca_mode;
	  
	  double_t	da2_k_1_;
	  double_t	da2_k_2_;
	  double_t	da2_k_3_;
	  double_t	da2_k_4_;
	  double_t	da2_k_d_;
	  double_t	da2_g_max_;
	  double_t	da2_T_;
	  double_t	da2_pulse_;
	  int_t		da2_n_;
	  
	  double_t	da2_suppress_f_;	//!< variable (< 1) which determines extent of reduction of synaptic strength
	  double_t	da2_suppress_f_const_;
	  double_t	da2_suppress_f_th_;
	  double_t	da2_suppress_f_sigma_;

	  bool	  output_all_variable_data;
	  bool	  output_conductance_data;
		
      bool    inactivate_GABAA_synapse;
  
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };
    
  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {
      /**
        * Enumeration identifying elements in state array State_::y_.
        * The state vector must be passed to GSL as a C array. This enum
        * identifies the elements of the vector. It must be public to be
        * accessible from the iteration function.
        */
      enum StateVecElems
   	  {
   	    V_M   = 0,
   	    HH_N     ,  // 1
   	    HH_M     ,  // 2
   	    HH_H     ,  // 3
   	    HH_A     ,  // 4
   	    HH_B     ,  // 5
   	    HH_C     ,  // 6
   	    HH_D1    ,  // 7
   	    HH_D2    ,  // 8
   	    HH_P     ,  // 9
   	    HH_Q     ,  // 10
	    HH_R     ,  // 11
   	    CA_C     ,  // 12
		G_AMPA_A2,
		G_AMPA_A1,
		G_NMDA_EXP,
		G_GABAA_A2,
		G_GABAA_A1,
		G_GABAB_A2,
		G_GABAB_A1,
   	    STATE_VEC_SIZE
      };


      /** soma **/
	  double_t y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
		
	  double_t I_NMDA_;
	  double_t I_AMPA_; 
	  double_t I_GABAA_;
	  double_t I_GABAB_;
		
		double_t I_K;
		double_t I_Na;
		double_t I_A;
		double_t I_LCa;
		double_t I_T;
		double_t I_CaK;
		
	  
	  /** dopaminergic response **/
	  double_t	da2_g_out_;
	  double_t	da2_r_;
	  double_t	da2_s_;

      int_t    just_after_spike_;    //!<	  variable previously called 'just_after_spike_' 
					  /*!<	  when membrane potential becomes bigger than V_th_ and 
							  a spike is generated, this variable takes value of 1. 
							  Once membrane potential decreased below V_th_, it takes 
							  value of 0. 
					  */

	  //std::vector<vector<double> > yn_;	
	  //double s_nmda_;	
		
      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

//  private:

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(stn_cond_neuron6&); //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, stn_cond_neuron6&); //!<Sets buffer pointers to 0

      //! Logger for all analog data
      UniversalDataLogger<stn_cond_neuron6> logger_;

      /** buffers and sums up incoming spikes/currents */
      RingBuffer spike_ampa_;
	  RingBuffer spike_nmda_;	
      RingBuffer spike_gabaa_;
	  RingBuffer spike_gabab_;
      RingBuffer spike_dopamine_;
      RingBuffer currents_;
		
	  /** GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system
      
      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      double_t step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL

      /**
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      double I_stim_;
    };

     // ---------------------------------------------------------------- 

     /**
      * Internal variables of the model.
      */
     struct Variables_ { 
	  	  
	  /** initial value to normalise excitatory synaptic conductance */
	  
	  /**
		* Amount added to y[5] for each arriving excitatory spike.
		* This evokes postsynaptic conductance changes with a peak
		* amplitude of 1nS.
	  **/
		double_t delta_ampa_;
   
	  /**
		* Amount added to y[7] for each arriving inhibitory spike.
		* This evokes postsynaptic conductance changes with a peak
		* amplitude of 1nS.
	  **/
		double_t delta_gabaa_;
		double_t delta_gabab_; 
	  
		double_t  da2_exp_A_00_;
		double_t  da2_exp_A_10_;
		double_t  da2_exp_A_11_;		
		double_t  da2_exp_r_add_;
		double_t  da2_r_inf_;
		 
		 
     };

    // ---------------------------------------------------------------- 
    // Access functions for UniversalDataLogger -------------------------------

    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    double get_y_elem_() const { return S_.y_[elem]; }
	double get_I_NMDA_() const { return S_.I_NMDA_;} 
	double get_I_AMPA_() const { return S_.I_AMPA_;} 
	double get_I_GABAA_() const { return S_.I_GABAA_;} 
	double get_I_GABAB_() const { return S_.I_GABAB_;} 
	  
	  double get_I_K() const { return S_.I_K;} 
	  double get_I_Na() const { return S_.I_Na;} 
	  double get_I_A() const { return S_.I_A;} 
	  double get_I_LCa() const { return S_.I_LCa;} 
	  double get_I_T() const { return S_.I_T;} 
	  double get_I_CaK() const { return S_.I_CaK;}

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to acces functions
    static RecordablesMap<stn_cond_neuron6> recordablesMap_;
  };
  
  inline
  port stn_cond_neuron6::check_connection(Connection& c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }
 
  inline
  port stn_cond_neuron6::connect_sender(CurrentEvent&, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port stn_cond_neuron6::connect_sender(DataLoggingRequest& dlr, 
                                             port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  void stn_cond_neuron6::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    (*d)[names::recordables] = recordablesMap_.get_list();

    def<double>(d, names::t_spike, get_spiketime_ms());
  }

  inline
  void stn_cond_neuron6::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);                 // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    Archiving_Node::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }

} // namespace nest


#endif //HAVE_GSL
#endif //STN_COND_NEURON6_H
