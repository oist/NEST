/*
 *  izhikevich_cond.h
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

#ifndef IZHIKEVICH_COND_H
#define IZHIKEVICH_COND_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_exp.h>


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
  int izhikevich_cond_dynamics(double, const double*, double*, void*);

 /* BeginDocumentation
     Name: izhikevich_cond - 

     Simple Model of Spiking Neurons

     Description:

     Implementation of the simple spiking neuron model introduced by Izhikevich [1].

     Remarks:	

     Parameters: 

     The following parameters can be set in the status dictionary.

	  Vm		  double - Membrane potential in mV 
	  Um		  double - Membrane potential recovery variable

	  Vr_		  double - resting potential [mV]
	  Vt_		  double - threshold potential [mV]
	  Vpeak_	
    
	  Cm		  double - membrane capacitance [nF]
	
	  k_	
	  a_		  double - describes time scale of recovery variable
	  b_		  double - sensitivity of recovery variable
	  c_		  double - after-spike reset value of V_m
	  d_		  double - after-spike reset value of U_m 
	  
	  I_e        double - Constant input current in pA. (R=1)

 
     References:
     [1] Izhikevich, Simple Model of Spiking Neurons,
     IEEE Transactions on Neural Networks (2003) 14:1569-1572 

     Sends: SpikeEvent

     Receives: SpikeEvent, CurrentEvent, PotentialRequest

  */

  /**
   * Simple Model of Spiking Neurons
   */


  class izhikevich_cond : public Archiving_Node
  {
    
  public:
    
    izhikevich_cond();
    izhikevich_cond(const izhikevich_cond&);
    ~izhikevich_cond();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */

    using Node::handles_test_event;
    using Node::handle;

    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &);

    port send_test_event(Node &, rport, synindex, bool);
    port handles_test_event(SpikeEvent &, rport);
    port handles_test_event(CurrentEvent &, rport);
    port handles_test_event(DataLoggingRequest &, rport);
    
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
    friend int izhikevich_cond_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<izhikevich_cond>;
    friend class UniversalDataLogger<izhikevich_cond>;

  private:
    
    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
     * These parameters must be passed to the iteration function that
     * is passed to the GSL ODE solvers. Since the iteration function
     * is a C++ function with C linkage, the parameters can be stored
     * in a C++ struct with member functions, as long as we just pass
     * it by void* from C++ to C++ function. The struct must be public,
     * though, since the iteration function is a function with C-linkage,
     * whence it cannot be a member function of izhikevich_cond.
     * @note One could achieve proper encapsulation by an extra level
     *       of indirection: Define the iteration function as a member
     *       function, plus an additional wrapper function with C linkage.
     *       Then pass a struct containing a pointer to the node and a
     *       pointer-to-member-function to the iteration function as void*
     *       to the wrapper function. The wrapper function can then invoke
     *       the iteration function on the node (Stroustrup, p 418). But
     *       this appears to involved, and the extra indirections cost.
     */
    struct Parameters_ {
	  double	Vr_;		//!< resting potential [mV]
	  double	Vt_;		//!< threshold potential [mV]
	  double	Vpeak_;
	  double	Cm_;        //!< membrane capacitance [nF]
	  
	  double	k_;
	  double	a_;
	  double	b_;
	  double	c_;
	  double	d_;
	  	  
      double tau_synE;    //!< Synaptic Time Constant Excitatory Synapse in ms
      double tau_synI;    //!< Synaptic Time Constant for Inhibitory Synapse in ms

	  double V_E;			//!<Reversal potential, excitatory synapse [mV]
	  double V_I;			//!<Reversal potential, inhibitory synapse [mV]

	  
      double I_e;         //!< Constant Current in pA
  
      /** 
       * External input current from CurrentEvents.
       * This is not a parameter but a variable. It is still placed here, since
       * it needs to be passed to the iteration function. We thus avoid the need
       * of an additional wrapper structure. It is not revealed or manipulateable.
       */
//      double I_stim;      //!< External Stimulus in pA
  
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
   	     U_M     ,  // 1
//   	     HH_H     ,  // 2
//   	     HH_R     ,  // 3
//   	     CA_C     ,  // 4
   	     STATE_VEC_SIZE
       };


      double	y_[6];  //!<  neuron state, must be C-array for GSL solver
	  
	  int_t		r_;		/*!<  when membrane potential becomes bigger than V_th_ and 
							  a spike is generated, this variable takes value of 1. 
							  Once membrane potential decreased below V_th_, it takes 
							  value of 0. */
									
      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

  private:

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      //Buffers_(); //!<Sets buffer pointers to 0
      Buffers_(izhikevich_cond&);                  //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, izhikevich_cond&); //!<Sets buffer pointers to 0

      //! Logger for all analog data
      UniversalDataLogger<izhikevich_cond> logger_;

      /** buffers and sums up incoming spikes/currents */
      RingBuffer spike_exc_;
      RingBuffer spike_inh_;
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
      double step_;           //!< step size in ms
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
      double PSConInit_E_; 
   
      /** initial value to normalise inhibitory synaptic conductance */
      double PSConInit_I_;    

      int_t    RefractoryCounts_;
     };

    // ---------------------------------------------------------------- 
    // Access functions for UniversalDataLogger -------------------------------

    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    double get_y_elem_() const { return S_.y_[elem]; }

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static RecordablesMap<izhikevich_cond> recordablesMap_;
  };
  
  inline
  port izhikevich_cond::send_test_event(Node& target, rport receptor_type, synindex, bool) 
  { 
    SpikeEvent e;
    e.set_sender(*this); 
    return target.handles_test_event(e, receptor_type); 
  }

  inline
  port izhikevich_cond::handles_test_event(SpikeEvent&, rport receptor_type)
  {
    if (receptor_type <= 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return receptor_type;
  }
 
  inline
  port izhikevich_cond::handles_test_event(CurrentEvent&, rport receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port izhikevich_cond::handles_test_event(DataLoggingRequest& dlr, rport receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  void izhikevich_cond::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    (*d)[names::recordables] = recordablesMap_.get_list();

    def<double>(d, names::t_spike, get_spiketime_ms());
  }

  inline
  void izhikevich_cond::set_status(const DictionaryDatum &d)
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

} // namespace bgnest


#endif //HAVE_GSL
#endif //IZHIKEVICH_COND_H
