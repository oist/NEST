/*
 *  stn_cond_neuron6.cpp
 *
 *  Created by Osamu Shouno	2009/02/09
 *  Copyright 2011 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *  Modified by Makoto Otsuka   2012/03/15
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "stn_cond_neuron6.h"

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

nest::RecordablesMap<nest::stn_cond_neuron6> nest::stn_cond_neuron6::recordablesMap_;

namespace nest 
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<stn_cond_neuron6>::create()
  {
    // use standard names whereever you can for consistency!
    insert_("V_m",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::V_M>);
    insert_("n",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_N>);
    insert_("m",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_M>);
    insert_("h",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_H>);
    insert_("a",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_A>);
    insert_("b",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_B>);
    insert_("c",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_C>);
    insert_("d1",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_D1>);
    insert_("d2",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_D2>);
    insert_("p",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_P>);
    insert_("q",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_Q>);
    insert_("r",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::HH_R>);
    insert_("Ca",
            &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::CA_C>);
	insert_("G_AMPA_A2",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_AMPA_A2>);
	insert_("G_AMPA_A1",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_AMPA_A1>);
	insert_("G_NMDA_EXP",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_NMDA_EXP>);
	insert_("G_GABAA_A2",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_GABAA_A2>);
	insert_("G_GABAA_A1",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_GABAA_A1>); 
	insert_("G_GABAB_A2",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_GABAB_A2>);
	insert_("G_GABAB_A1",
			  &stn_cond_neuron6::get_y_elem_<stn_cond_neuron6::State_::G_GABAB_A1>); 	  
	insert_("I_AMPA", &stn_cond_neuron6::get_I_AMPA_);
	insert_("I_NMDA", &stn_cond_neuron6::get_I_NMDA_);
	insert_("I_GABAA", &stn_cond_neuron6::get_I_GABAA_);
	insert_("I_GABAB", &stn_cond_neuron6::get_I_GABAB_);
	  
	insert_("I_K", &stn_cond_neuron6::get_I_K);
	insert_("I_Na", &stn_cond_neuron6::get_I_Na);
	insert_("I_A", &stn_cond_neuron6::get_I_A);
	insert_("I_LCa", &stn_cond_neuron6::get_I_LCa);
	insert_("I_T", &stn_cond_neuron6::get_I_T);
	insert_("I_CaK", &stn_cond_neuron6::get_I_CaK);	  
  }



	extern "C"
	int stn_cond_neuron6_dynamics (
			double,						// time, not needed
			const double y[],			// state
			double dydt[],				// derivative
          void* pnode)				   // parameters
	{ 
		// a shorthand
		typedef nest::stn_cond_neuron6::State_ S;

		// get access to node so we can almost work as in a member function
		assert(pnode);
		const nest::stn_cond_neuron6& node =  *(reinterpret_cast<nest::stn_cond_neuron6*>(pnode));

		// y[] here is---and must be---the state vector supplied by the integrator,
		// not the state vector in the node, node.S_.y[].

		// The following code is verbose for the sake of clarity. We assume that a
		// good compiler will optimize the verbosity away ...

		// delayed-rectified K current I_K
		double 
		I_K =	node.P_.g_K * pow(y[S::HH_N], 4.0)*(y[S::V_M] - node.P_.V_K);
  
		// transient Na current  I_Na
		double 
		I_Na =	node.P_.g_Na
					* pow( y[S::HH_M], 3.0)				// m^3
					* y[S::HH_H]						// h: Na current inactivation variable
					* (y[S::V_M] - node.P_.V_Na);
  
		// I_A
		double 
		I_A =	node.P_.g_IA
					* pow(y[S::HH_A], 2.0)				// a^2
					* y[S::HH_B]						// b
					* (y[S::V_M] - node.P_.V_K);
  
		// L-type Ca current  I_LCa
		double 
		I_LCa =	node.P_.g_LCa
					* pow(y[S::HH_C], 2.0)				// c^2
					* y[S::HH_D1]						// d1
					* y[S::HH_D2]						// d2
					* (y[S::V_M] - node.P_.V_Ca);
					
		// T current  I_T
		double 
		I_T =	node.P_.g_T*
					pow( y[S::HH_P], 2.0)				// p^2
					* y[S::HH_Q]						// q
		* (y[S::V_M] - node.P_.V_Ca) ;
					//* node.P_.da2_suppress_f_;

		// Calcium-depndent pottasium (AHP) current  I_CaK
		double 
		I_CaK =	node.P_.g_CaK
					* (y[S::V_M] - node.P_.V_K)
					* pow(y[S::HH_R], 4.0);				// r^4

		
		double I_NMDA = node.P_.g_nmda * (y[S::V_M] - node.P_.V_E) /(1.0 + node.P_.MgConc * std::exp(-0.062 * y[S::V_M])/3.57)*y[S::G_NMDA_EXP];
		double I_AMPA = y[S::G_AMPA_A2]* ( y[S::V_M] - node.P_.V_E );		
		double I_GABAA = y[S::G_GABAA_A2]* ( y[S::V_M] - node.P_.V_GABAA );
		double I_GABAB = y[S::G_GABAB_A2]* ( y[S::V_M] - node.P_.V_GABAB );
					
		// dV / dt
		// V  membrane potential
		dydt[S::V_M] =(	- node.P_.g_L * ( y[S::V_M] - node.P_.V_L )	// leak current
				- I_K								// delayed-rectified K current
				- I_Na								// transient Na current
				- I_A								// A current
				- I_T								// T current
				- I_LCa								// high threshold Ca current
				- I_CaK								// Calcium-depndent AHP current
				- I_AMPA	/ node.P_.area			// ampa synaptic input
				- I_GABAA	/ node.P_.area			// gaba synaptic input
				- I_GABAB	/ node.P_.area  
				- I_NMDA	/ node.P_.area
                + node.P_.I0 / node.P_.area			// persistent current
                + node.B_.I_stim_ / node.P_.area	// injected current
			 ) / node.P_.C_m;

		
		// dn/dt
		// n	I_K activation variable
		dydt[S::HH_N] =	(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_n)/node.P_.sigma_n)))	- y[S::HH_N] ) // (n_infinity -y[S::HH_N])
				/(	node.P_.tau_n0 + node.P_.tau_n1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_n0)/node.P_.sigma_tau_n0) 
					+ std::exp(-(y[S::V_M] - node.P_.theta_tau_n1)/node.P_.sigma_tau_n1)))
				;
		// dm/dt
		// m	I_Na activation variable
		dydt[S::HH_M] =	(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_m)/node.P_.sigma_m)))	- y[S::HH_M] ) // (m_infinity -y[S::HH_M])
					/(	node.P_.tau_m0 + node.P_.tau_m1/(1.0 + std::exp(-(y[S::V_M] - node.P_.theta_tau_m )/node.P_.sigma_tau_m) ))
					;
		// dh/dt
		// h	I_Na inactivation variable
		dydt[S::HH_H] =	(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_h)/node.P_.sigma_h)))	- y[S::HH_H] ) // (h_infinity -y[S::HH_H])
					/(	node.P_.tau_h0 + node.P_.tau_h1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_h0)/node.P_.sigma_tau_h0) 
					+ std::exp(-(y[S::V_M] - node.P_.theta_tau_h1)/node.P_.sigma_tau_h1)))
				;
		// da/dt
		// a	I_A activation variable
		dydt[S::HH_A] =		(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_a)/node.P_.sigma_a)))	- y[S::HH_A] ) // (a_infinity -y[S::HH_A])
					/(	node.P_.tau_a0 + node.P_.tau_a1/(1.0 + std::exp(-(y[S::V_M] - node.P_.theta_tau_a )/(node.P_.sigma_tau_a))))
		//				/(	1.0 + 1.0/(1.0 + std::exp(-(y[S::V_M] +40.0 )/(-0.5))))
					;
		// db/dt
		// b	I_A inactivation variable
		dydt[S::HH_B] =		(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_b)/node.P_.sigma_b)))	- y[S::HH_B] ) // (b_infinity -y[S::HH_B])
					/(	node.P_.tau_b0 + node.P_.tau_b1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_b0)/node.P_.sigma_tau_b0) 
					+ std::exp(-(y[S::V_M] - node.P_.theta_tau_b1)/node.P_.sigma_tau_b1)))
 				;
		// dc/dt
		// c	I_LCa activation variable
		dydt[S::HH_C] =		(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_c)/node.P_.sigma_c)))	- y[S::HH_C] ) // (c_infinity -y[S::HH_C])
					/(	node.P_.tau_c0 + node.P_.tau_c1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_c0)/node.P_.sigma_tau_c0) 
					+ std::exp(-(y[S::V_M] - node.P_.theta_tau_c1)/node.P_.sigma_tau_c1)))
				;
		// dd1/dt
		// d1	I_LCa inactivation variable
		dydt[S::HH_D1]	=		(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_d1 )/node.P_.sigma_d1)))	- y[S::HH_D1] ) // (d1_infinity -y[S::SS_D1])
						/(	node.P_.tau_d0 + node.P_.tau_d1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_d0)/node.P_.sigma_tau_d0) 
						+ std::exp(-(y[S::V_M] - node.P_.theta_tau_d1)/node.P_.sigma_tau_d1)))
				;
		// dd2/dt
		// d2	I_LCa inactivation variable
		dydt[S::HH_D2] =		(1.0 /(1.0 + std::exp((y[S::CA_C]- node.P_.k1_d2)/(node.P_.sigma_d2)))- y[S::HH_D2])/130.0
				;

		// dp/dt
		// p	I_T activation variable
		dydt[S::HH_P] =		(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_p)/node.P_.sigma_p)))	- y[S::HH_P] ) // (p_infinity -y[S::HH_P])
				/(	node.P_.tau_p0 + node.P_.tau_p1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_p0)/node.P_.sigma_tau_p0) 
					+ std::exp(-(y[S::V_M] - node.P_.theta_tau_p1)/node.P_.sigma_tau_p1)))
				;
		// dq/dt
		// q	I_T inactivation variable
		dydt[S::HH_Q] =		(	(1.0/(1.0 + std::exp((y[S::V_M] - node.P_.theta_q )/node.P_.sigma_q)))	- y[S::HH_Q] ) // (q_infinity -y[S::HH_Q])
				/(	node.P_.tau_q0 + node.P_.tau_q1/(std::exp(-(y[S::V_M] - node.P_.theta_tau_q0)/node.P_.sigma_tau_q0) 
					+ std::exp(-(y[S::V_M] - node.P_.theta_tau_q1)/node.P_.sigma_tau_q1)))
				;
		// dr/dt
		// r	I_CaK activation variable
		dydt[S::HH_R] =		(1.0 /(1.0 + std::exp((y[S::CA_C] - node.P_.k1_r)/(-0.08)))- y[S::HH_R])/2.0;
  
		// d[Ca]/dt
		// [Ca]
		dydt[S::CA_C] =		node.P_.epsilon*( - I_LCa - I_T) - node.P_.k_Ca*y[S::CA_C];
		
		
		dydt[S::G_AMPA_A2] =  y[S::G_AMPA_A1] - y[S::G_AMPA_A2]/node.P_.tau_syn_ampa;
		dydt[S::G_AMPA_A1] = -y[S::G_AMPA_A1]/node.P_.tau_syn_ampa;
		
		dydt[S::G_NMDA_EXP]= -y[S::G_NMDA_EXP]/node.P_.tau_nmda_decay;
		
		dydt[S::G_GABAA_A2] =  y[S::G_GABAA_A1] * (1.0/node.P_.tau_syn_gabaa2 - 1.0/node.P_.tau_syn_gabaa1) - y[S::G_GABAA_A2]/node.P_.tau_syn_gabaa1;
		dydt[S::G_GABAA_A1] = -y[S::G_GABAA_A1]/node.P_.tau_syn_gabaa2;
		
		dydt[S::G_GABAB_A2] =  y[S::G_GABAB_A1] - y[S::G_GABAB_A2]/node.P_.tau_syn_gabab;
		dydt[S::G_GABAB_A1] = -y[S::G_GABAB_A1]/node.P_.tau_syn_gabab;

		return GSL_SUCCESS;
	}
	

	/* ---------------------------------------------------------------- 
	 * Default constructors defining default parameters and state
	 * ---------------------------------------------------------------- */
    
	nest::stn_cond_neuron6::Parameters_::Parameters_()
	 :	V_th_	( -30.0   ),  // mV
		C_m		(	1.0	  ),  // pF/micro-meter^2
  
		V_L		( -60.0	  ),  // mV						 original value
		V_K		( -90.0	  ),  // mV						 original value
		V_Na	(  60.0	  ),  // mV						 original value
		V_Ca	( 140.0	  ),  // mV	missed in table 1 Otsuka et al 2004
	
		g_L		(	0.35  ),  // nS/micro-meter^2		 original value
		g_K		(  57.0	  ),  // nS/micro-meter^2	 original value
		g_Na	(  49.0	  ),  // nS/micro-meter^2	 original value
		g_IA	(	5.0	  ),  // nS/micro-meter^2	 original value
		g_LCa	(  15.0	  ),  // nS/micro-meter^2	 original value
		g_T		(	5.0	  ),  // nS/micro-meter^2	 original value
		g_CaK	(	1.0	  ),  // nS/micro-meter^2	 original value
  
		tau_n0	(	0.0	  ),  // msec	 original value
		tau_n1	(  11.0	  ),  // msec	 original value
		tau_m0	(	0.2	  ),  // msec	 original value
		tau_m1	(	3.0	  ),  // msec		  modified value 3.9	--> original value 3.0
		tau_h0	(	0.0	  ),  // msec	 original value
		tau_h1	(  24.5	  ),  // msec	 original value
		tau_b0	(	0.0	  ),  // msec	 original value
		tau_b1	( 200.0	  ),  // msec	 original value
		tau_c0	(  45.0	  ),  // msec	 original value
		tau_c1	(  10.0	  ),  // msec	 original value
		tau_d0	( 400.0	  ),  // msec	 original value
		tau_d1	( 500.0	  ),  // msec	 original value
		tau_p0	(   5.0	  ),  // msec	 original value
		tau_p1	(	0.33  ),  // msec		 original value
		tau_q0	(  10.0	  ),  // msec	  modified value 10.0 , original value 0.0
		tau_q1	( 400.0  ),  // msec		  modified value 1600.0 --> original value 400.0

		tau_a0	(   1.0 ), // msec
		tau_a1	(   1.0 ), // msec
  
		theta_n	( -41.0	  ),  // mV	 original value
		theta_m	( -40.0	  ),  // mV	 original value
		theta_h	( -45.5	  ),  // mV	 original value
		theta_a	( -45.0	  ),  // mV	 original value
		theta_b	( -90.0	  ),  // mV	 original value
		theta_c	( -30.6	  ),  // mV	 original value
		theta_d1	( -60.0	  ),  // mV	 original value
		theta_p	( -56.0	  ),  // mV	 original value
		theta_q	( -85.0	  ),  // mV	 original value
  
		sigma_n	( -14.0	  ),  // mV
		sigma_m	(  -8.0	  ),  // mV
		sigma_h	(	6.4	  ),  // mV
		sigma_a	( -14.7	  ),  // mV	 original value
		sigma_b	(	7.5	  ),  // mV	 original value
		sigma_c	(  -5.0	  ),  // mV	 original value
		sigma_d1	(	7.5	  ),  // mV	 original value
		sigma_p	(  -6.7	  ),  // mV
		sigma_q	(	5.8	  ),  // mV

		sigma_d2	(0.02),					//	 original value
  
		theta_tau_m   (	  -53.0	  ),  // mV		original value
		theta_tau_a   (	  -40.0	  ),  // mV		original value
		theta_tau_n0  (	  -40.0	  ),  // mV		original value
		theta_tau_n1  (	  -40.0	  ),  // mV		original value
		theta_tau_h0  (	  -50.0	  ),  // mV		original value
		theta_tau_h1  (	  -50.0	  ),  // mV		original value
		theta_tau_b0  (	  -60.0	  ),  // mV		original value
		theta_tau_b1  (	  -40.0	  ),  // mV		original value
		theta_tau_c0  (	  -27.0	  ),  // mV		original value
		theta_tau_c1  (	  -50.0	  ),  // mV		original value
		theta_tau_d0  (	  -40.0	  ),  // mV		original value
		theta_tau_d1  (	  -20.0	  ),  // mV		original value
		theta_tau_p0  (	  -27.0	  ),  // mV		original value
		theta_tau_p1  (	 -102.0	  ),  // mV		original value
		theta_tau_q0  (	  -50.0	  ),  // mV		original value
		theta_tau_q1  (	  -50.0	  ),  // mV		original value
  
		sigma_tau_m   (	   -0.7	  ),  // mV		original value
		sigma_tau_a   (	   -0.5	  ),  // mV		original value
		sigma_tau_n0  (	  -40.0	  ),  // mV		original value
		sigma_tau_n1  (	   50.0	  ),  // mV		original value
		sigma_tau_h0  (	  -15.0	  ),  // mV		original value
		sigma_tau_h1  (	   16.0	  ),  // mV		original value
		sigma_tau_b0  (	  -30.0	  ),  // mV		original value
		sigma_tau_b1  (	   10.0	  ),  // mV		original value
		sigma_tau_c0  (	  -20.0	  ),  // mV		original value
		sigma_tau_c1  (	   15.0	  ),  // mV		original value
		sigma_tau_d0  (	  -15.0	  ),  // mV		original value
		sigma_tau_d1  (	   20.0	  ),  // mV		original value
		sigma_tau_p0  (	  -10.0	  ),  // mV		original value
		sigma_tau_p1  (	   15.0	  ),  // mV		original value
		sigma_tau_q0  (	  -15.0	  ),  // mV		original value
		sigma_tau_q1  (	   16.0	  ),  // mV		original value
	
		k1_r		  (	  0.17	  ),			// original value
		k1_d2		  (	  0.1	  ),				// original value
		k_Ca		(	2.0	  ),
		epsilon		( 5.182 * pow(10.0, -3.0) ), // for [Ca2+] in micro-M. 5.182 * pow(10.0, -6.0) for [Ca] in mM, 
	
		tau_syn_ampa  ( 1.0	  ),  // ms
	    tau_syn_gabaa_rise  ( 0.4  ),  // ms
	    tau_syn_gabaa_decay  ( 7.7  ),  // ms
	
		tau_syn_gabab  ( 175.2  ),  // ms

		ampa_conductance (0.0 ),
		gabaa_conductance (0.0 ),
	    gabab_conductance (0.0 ),
	
	
		MgConc			(1.0),		// mM
		tau_nmda_decay	(100.0),	// msec
		tau_nmda_rise	(2.0),		// msec
		alpha_nmda		(0.5),		// msec^-1
		g_nmda			(1.0),
	
		g_ampa			(1.0),

		V_E			(	0.0	  ),  // mV
		V_GABAA		( -78.9	  ),  // mV
		V_GABAB		( -93.3	  ),  // mV
	
		I0			(	0.0	  ),
	    area(0.1),

		da2_Ca_modulation (	0.0	),

		da2_k_1_	(  90.0	  ),
		da2_k_2_	(	0.0036),
		da2_k_3_	(	0.18  ),
		da2_k_4_	(	0.034 ),
		da2_k_d_	( 100.0	  ),
		da2_g_max_	(	1.0	  ),
		da2_T_		(	1e-3  ),
		da2_pulse_	(	1.105 ),
		da2_n_		(	4	  ),
  
		da2_suppress_f_			(	1.0	  ),
		da2_suppress_f_const_	(	0.85  ),
		da2_suppress_f_th_		(	0.5	  ),
		da2_suppress_f_sigma_	(	0.1	  ),
	
	
		output_all_variable_data  (false),
		output_conductance_data	(false),
	
	    inactivate_GABAA_synapse(false)
	{
		
		tau_syn_gabaa1 = tau_syn_gabaa_decay;
		tau_syn_gabaa2 = tau_syn_gabaa1 * tau_syn_gabaa_rise /(tau_syn_gabaa1 + tau_syn_gabaa_rise);
	}

	nest::stn_cond_neuron6::State_::State_(const Parameters_& p)
	 :	da2_g_out_(0.0),
		da2_r_(0.0),
		da2_s_(0.0),
		just_after_spike_(0)
	{
		y_[State_::V_M] = p.V_L;
	}

	nest::stn_cond_neuron6::State_::State_(const State_& s)
	 :	da2_g_out_(s.da2_g_out_),
		da2_r_(s.da2_r_),
		da2_s_(s.da2_s_),
		just_after_spike_(s.just_after_spike_)
	{
		for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
			y_[i] = s.y_[i];
	}

	nest::stn_cond_neuron6::State_& nest::stn_cond_neuron6::State_::operator=(const State_& s)
	{
		assert(this != &s);  // would be bad logical error in program
  
		for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
			y_[i] = s.y_[i];
  
		da2_g_out_ =s.da2_g_out_;
		da2_r_	= s.da2_r_;
		da2_s_	= s.da2_s_;
		just_after_spike_ = s.just_after_spike_;
  
		return *this;
	}

	/* ---------------------------------------------------------------- 
	 * Paramater and state extractions and manipulation functions
	 * ---------------------------------------------------------------- */

	void nest::stn_cond_neuron6::Parameters_::get(DictionaryDatum &d) const
	{
		def<double>(d,names::V_th,	V_th_);
		def<double>(d,names::C_m,	C_m);
		def<double>(d, "V_L",		V_L);
		def<double>(d, "V_K",		V_K);
		def<double>(d, "V_Na",		V_Na);
		def<double>(d, "V_Ca",		V_Ca);
  
		def<double>(d, "g_L",		g_L);
		def<double>(d, "g_K",		g_K);
		def<double>(d, "g_Na",		g_Na);
		def<double>(d, "g_IA",		g_IA);
		def<double>(d, "g_LCa",		g_LCa);
		def<double>(d, "g_T",		g_T);
		def<double>(d, "g_CaK",		g_CaK);
  
		def<double>(d, "tau_n0",	tau_n0);
		def<double>(d, "tau_n1",	tau_n1);
		def<double>(d, "tau_m0",	tau_m0);
		def<double>(d, "tau_m1",	tau_m1);
		def<double>(d, "tau_h0",	tau_h0);
		def<double>(d, "tau_h1",	tau_h1);
		def<double>(d, "tau_b0",	tau_b0);
		def<double>(d, "tau_b1",	tau_b1);
		def<double>(d, "tau_c0",	tau_c0);
		def<double>(d, "tau_c1",	tau_c1);
		def<double>(d, "tau_d0",	tau_d0);
		def<double>(d, "tau_d1",	tau_d1);
		def<double>(d, "tau_p0",	tau_p0);
		def<double>(d, "tau_p1",	tau_p1);
		def<double>(d, "tau_q0",	tau_q0);
		def<double>(d, "tau_q1",	tau_q1);
  
		def<double>(d, "tau_a0",	tau_a0);
		def<double>(d, "tau_a1",	tau_a1);

		def<double>(d, "theta_n",	theta_n);
		def<double>(d, "theta_m",	theta_m);
		def<double>(d, "theta_h",	theta_h);
		def<double>(d, "theta_a",	theta_a);
		def<double>(d, "theta_b",	theta_b);
		def<double>(d, "theta_c",	theta_c);
		def<double>(d, "theta_d1",	theta_d1);
		def<double>(d, "theta_p",	theta_p);
		def<double>(d, "theta_q",	theta_q);
  
		def<double>(d, "sigma_n",	sigma_n);
		def<double>(d, "sigma_m",	sigma_m);
		def<double>(d, "sigma_h",	sigma_h);
		def<double>(d, "sigma_a",	sigma_a);
		def<double>(d, "sigma_b",	sigma_b);
		def<double>(d, "sigma_c",	sigma_c);
		def<double>(d, "sigma_d1",	sigma_d1); 
		def<double>(d, "sigma_d2",	sigma_d2); 
		def<double>(d, "sigma_p",	sigma_p);
		def<double>(d, "sigma_q",	sigma_q); 
  
		def<double>(d, "theta_tau_m",		 theta_tau_m);
		def<double>(d, "theta_tau_a",		 theta_tau_a);
		def<double>(d, "theta_tau_n0",	 theta_tau_n0);
		def<double>(d, "theta_tau_n1",	 theta_tau_n1);
		def<double>(d, "theta_tau_h0",	 theta_tau_h0);
		def<double>(d, "theta_tau_h1",	 theta_tau_h1);
		def<double>(d, "theta_tau_b0",	theta_tau_b0);
		def<double>(d, "theta_tau_b1",	theta_tau_b1);
		def<double>(d, "theta_tau_c0",	theta_tau_c0);
		def<double>(d, "theta_tau_c1",	theta_tau_c1);
		def<double>(d, "theta_tau_d0",	theta_tau_d0);
		def<double>(d, "theta_tau_d1",	theta_tau_d1);
		def<double>(d, "theta_tau_p0",	theta_tau_p0);
		def<double>(d, "theta_tau_p1",	theta_tau_p1);
		def<double>(d, "theta_tau_q0",	theta_tau_q0);
		def<double>(d, "theta_tau_q1",	theta_tau_q1);

		def<double>(d, "sigma_tau_m",	   sigma_tau_m);
		def<double>(d, "sigma_tau_a",	   sigma_tau_a);
		def<double>(d, "sigma_tau_n0",	sigma_tau_n0);
		def<double>(d, "sigma_tau_n1",	sigma_tau_n1);
		def<double>(d, "sigma_tau_h0",	sigma_tau_h0);
		def<double>(d, "sigma_tau_h1",	sigma_tau_h1); 
		def<double>(d, "sigma_tau_b0",	sigma_tau_b0);
		def<double>(d, "sigma_tau_b1",	sigma_tau_b1); 
		def<double>(d, "sigma_tau_c0",	sigma_tau_c0);
		def<double>(d, "sigma_tau_c1",	sigma_tau_c1); 
		def<double>(d, "sigma_tau_d0",	sigma_tau_d0);
		def<double>(d, "sigma_tau_d1",	sigma_tau_d1); 
		def<double>(d, "sigma_tau_p0",	sigma_tau_p0);
		def<double>(d, "sigma_tau_p1",	sigma_tau_p1); 
		def<double>(d, "sigma_tau_q0",	sigma_tau_q0);
		def<double>(d, "sigma_tau_q1",	sigma_tau_q1); 
  
		def<double>(d, "k1_r",			k1_r); 
		def<double>(d, "k1_d2",			k1_d2); 
		def<double>(d, "k_Ca",			k_Ca); 
		def<double>(d, "epsilon",		epsilon); 
  
		def<double>(d, "V_E",			V_E);
		def<double>(d, "V_GABAA",			V_GABAA);
		def<double>(d, "V_GABAB",			V_GABAB);
		def<double>(d, "tau_syn_ampa",	tau_syn_ampa);
		def<double>(d, "tau_syn_gabaa_rise",	tau_syn_gabaa_rise);
		def<double>(d, "tau_syn_gabaa_decay",	tau_syn_gabaa_decay);
		def<double>(d, "tau_syn_gabaa1",	tau_syn_gabaa1);
		def<double>(d, "tau_syn_gabaa2",	tau_syn_gabaa2);
		
		def<double>(d, "tau_syn_gabab",	tau_syn_gabab);
		
		def<double>(d, "MgConc",		MgConc);
		def<double>(d, "tau_nmda_decay",tau_nmda_decay);
		def<double>(d, "tau_nmda_rise",	tau_nmda_rise);
		def<double>(d, "alpha_nmda",	alpha_nmda);
		def<double>(d, "g_nmda",		g_nmda);
		
		def<double>(d, "g_ampa",		g_ampa);
		
		def<double>(d, "I0",			I0);
		def<double>(d, "area",			area);

		def<double>(d, "da2_Ca_modulation",	da2_Ca_modulation);
  
		def<double>(d, "da2_k_1",		da2_k_1_);
		def<double>(d, "da2_k_2",		da2_k_2_);
		def<double>(d, "da2_k_3",		da2_k_3_);
		def<double>(d, "da2_k_4",		da2_k_4_);
		def<double>(d, "da2_k_d",		da2_k_d_);
		def<double>(d, "da2_g_max",		da2_g_max_);
		def<double>(d, "da2_T",			da2_T_);
		def<double>(d, "da2_pulse",		da2_pulse_);
		def<double>(d, "da2_n",			da2_n_);
  
		def<double>(d, "da2_suppress_f",da2_suppress_f_);
  
		def<double>(d, "da2_suppress_f_const",  da2_suppress_f_const_);
		def<double>(d, "da2_suppress_f_th",		da2_suppress_f_th_);
		def<double>(d, "da2_suppress_f_sigma",  da2_suppress_f_sigma_);
  
		def<bool>(d, "output_all_variable_data",  output_all_variable_data);
		def<bool>(d, "output_conductance_data",  output_conductance_data);
		def<bool>(d, "inactivate_GABAA_synapse",  inactivate_GABAA_synapse);
	}

	void nest::stn_cond_neuron6::Parameters_::set(const DictionaryDatum& d)
	{
		// allow setting the membrane potential
		updateValue<double>(d,names::V_th,	V_th_);
		updateValue<double>(d,names::C_m,	C_m);
  
		updateValue<double>(d, "V_L",		V_L);
		updateValue<double>(d, "V_K",		V_K);
		updateValue<double>(d, "V_Na",		V_Na);
		updateValue<double>(d, "V_Ca",		V_Ca);
  
		updateValue<double>(d, "g_L",		g_L);
		updateValue<double>(d, "g_K",		g_K);
		updateValue<double>(d, "g_Na",		g_Na);
		updateValue<double>(d, "g_IA",		g_IA);
		updateValue<double>(d, "g_LCa",		g_LCa);
		updateValue<double>(d, "g_T",		g_T);
		updateValue<double>(d, "g_CaK",		g_CaK);
  
		updateValue<double>(d, "tau_n0",	tau_n0);
		updateValue<double>(d, "tau_n1",	tau_n1);
		updateValue<double>(d, "tau_m0",	tau_m0);
		updateValue<double>(d, "tau_m1",	tau_m1);
		updateValue<double>(d, "tau_h0",	tau_h0);
		updateValue<double>(d, "tau_h1",	tau_h1);
		updateValue<double>(d, "tau_b0",	tau_b0);
		updateValue<double>(d, "tau_b1",	tau_b1);
		updateValue<double>(d, "tau_c0",	tau_c0);
		updateValue<double>(d, "tau_c1",	tau_c1);
		updateValue<double>(d, "tau_d0",	tau_d0);
		updateValue<double>(d, "tau_d1",	tau_d1);
		updateValue<double>(d, "tau_p0",	tau_p0);
		updateValue<double>(d, "tau_p1",	tau_p1);
		updateValue<double>(d, "tau_q0",	tau_q0);
		updateValue<double>(d, "tau_q1",	tau_q1);
  
		updateValue<double>(d, "tau_a0",	tau_a0);
		updateValue<double>(d, "tau_a1",	tau_a1);

		updateValue<double>(d, "theta_n",	theta_n);
		updateValue<double>(d, "theta_m",	theta_m);
		updateValue<double>(d, "theta_h",	theta_h);
		updateValue<double>(d, "theta_a",	theta_a);
		updateValue<double>(d, "theta_b",	theta_b);
		updateValue<double>(d, "theta_c",	theta_c);
		updateValue<double>(d, "theta_d1",	theta_d1);
		updateValue<double>(d, "theta_p",	theta_p);
		updateValue<double>(d, "theta_q",	theta_q);
  
		updateValue<double>(d, "sigma_n",	sigma_n);
		updateValue<double>(d, "sigma_m",	sigma_m);
		updateValue<double>(d, "sigma_h",	sigma_h);
		updateValue<double>(d, "sigma_a",	sigma_a);
		updateValue<double>(d, "sigma_b",	sigma_b);
		updateValue<double>(d, "sigma_c",	sigma_c);
		updateValue<double>(d, "sigma_d1",	sigma_d1);
		updateValue<double>(d, "sigma_d2",	sigma_d2);
		updateValue<double>(d, "sigma_p",	sigma_p);
		updateValue<double>(d, "sigma_q",	sigma_q);
  
		updateValue<double>(d, "theta_tau_m",  theta_tau_m);
		updateValue<double>(d, "theta_tau_a",  theta_tau_a);
		updateValue<double>(d, "theta_tau_n0", theta_tau_n0);
		updateValue<double>(d, "theta_tau_n1", theta_tau_n1);
		updateValue<double>(d, "theta_tau_h0", theta_tau_h0);
		updateValue<double>(d, "theta_tau_h1", theta_tau_h1);
		updateValue<double>(d, "theta_tau_b0", theta_tau_b0);
		updateValue<double>(d, "theta_tau_b1", theta_tau_b1);
		updateValue<double>(d, "theta_tau_c0", theta_tau_c0);
		updateValue<double>(d, "theta_tau_c1", theta_tau_c1);
		updateValue<double>(d, "theta_tau_d0", theta_tau_d0);
		updateValue<double>(d, "theta_tau_d1", theta_tau_d1);
		updateValue<double>(d, "theta_tau_p0", theta_tau_p0);
		updateValue<double>(d, "theta_tau_p1", theta_tau_p1);
		updateValue<double>(d, "theta_tau_q0", theta_tau_q0);
		updateValue<double>(d, "theta_tau_q1", theta_tau_q1);

		updateValue<double>(d, "sigma_tau_m",  sigma_tau_m);
		updateValue<double>(d, "sigma_tau_a",  sigma_tau_a);
		updateValue<double>(d, "sigma_tau_n0", sigma_tau_n0);
		updateValue<double>(d, "sigma_tau_n1", sigma_tau_n1);
		updateValue<double>(d, "sigma_tau_h0", sigma_tau_h0);
		updateValue<double>(d, "sigma_tau_h1", sigma_tau_h1);
		updateValue<double>(d, "sigma_tau_b0", sigma_tau_b0);
		updateValue<double>(d, "sigma_tau_b1", sigma_tau_b1);
		updateValue<double>(d, "sigma_tau_c0", sigma_tau_c0);
		updateValue<double>(d, "sigma_tau_c1", sigma_tau_c1);
		updateValue<double>(d, "sigma_tau_d0", sigma_tau_d0);
		updateValue<double>(d, "sigma_tau_d1", sigma_tau_d1);
		updateValue<double>(d, "sigma_tau_p0", sigma_tau_p0);
		updateValue<double>(d, "sigma_tau_p1", sigma_tau_p1);
		updateValue<double>(d, "sigma_tau_q0", sigma_tau_q0);
		updateValue<double>(d, "sigma_tau_q1", sigma_tau_q1);
  
		updateValue<double>(d, "k1_r",		k1_r);
		updateValue<double>(d, "k1_d2",		k1_d2);
		updateValue<double>(d, "k_Ca",		k_Ca);
		updateValue<double>(d, "epsilon",	epsilon);
  
		updateValue<double>(d, "V_E",		V_E);
		updateValue<double>(d, "V_GABAA",		V_GABAA);
		updateValue<double>(d, "V_GABAB",		V_GABAB);
		updateValue<double>(d, "tau_syn_ampa", tau_syn_ampa);
		
		if (updateValue<double>(d, "tau_syn_gabaa_rise", tau_syn_gabaa_rise)){//
			tau_syn_gabaa1 = tau_syn_gabaa_decay;
			tau_syn_gabaa2 = tau_syn_gabaa1 * tau_syn_gabaa_rise /(tau_syn_gabaa1 + tau_syn_gabaa_rise);
		}
		
		if (updateValue<double>(d, "tau_syn_gabaa_decay", tau_syn_gabaa_decay)){
			tau_syn_gabaa1 = tau_syn_gabaa_decay;
			tau_syn_gabaa2 = tau_syn_gabaa1 * tau_syn_gabaa_rise /(tau_syn_gabaa1 + tau_syn_gabaa_rise);			
		}
		
		updateValue<double>(d, "tau_syn_gabab", tau_syn_gabab);
		
		updateValue<double>(d, "MgConc",		MgConc);
		updateValue<double>(d, "tau_nmda_decay",tau_nmda_decay);
		updateValue<double>(d, "tau_nmda_rise",	tau_nmda_rise);
		updateValue<double>(d, "alpha_nmda",	alpha_nmda);
		updateValue<double>(d, "g_nmda",		g_nmda);
		
		updateValue<double>(d, "g_ampa",		g_ampa);
  
		updateValue<double>(d, "I0",		I0);
		updateValue<double>(d, "area",		area);

		updateValue<double>(d, "da2_Ca_modulation", da2_Ca_modulation);
	
		updateValue<double>(d, "da2_k_1",	da2_k_1_);
		updateValue<double>(d, "da2_k_2",	da2_k_2_);
		updateValue<double>(d, "da2_k_3",	da2_k_3_);
		updateValue<double>(d, "da2_k_4",	da2_k_4_);
		updateValue<double>(d, "da2_k_d",	da2_k_d_);
		updateValue<double>(d, "da2_g_max",	da2_g_max_);
		updateValue<double>(d, "da2_T",		da2_T_);
		updateValue<double>(d, "da2_pulse",	da2_pulse_);
		updateValue<double>(d, "da2_n",		da2_n_);
  
		updateValue<double>(d, "da2_suppress_f",		da2_suppress_f_);
		updateValue<double>(d, "da2_suppress_f_const",  da2_suppress_f_const_);
		updateValue<double>(d, "da2_suppress_f_th",		da2_suppress_f_th_);
		updateValue<double>(d, "da2_suppress_f_sigma",  da2_suppress_f_sigma_);
  
		updateValue<bool>(d, "output_all_variable_data",output_all_variable_data);
		updateValue<bool>(d, "output_conductance_data",	output_conductance_data);
		updateValue<bool>(d, "inactivate_GABAA_synapse",inactivate_GABAA_synapse);
	}

	void nest::stn_cond_neuron6::State_::get(DictionaryDatum &d) const
	{
		def<double>(d, "V",		y_[0]); // Membrane potential
		def<double>(d, "n",		y_[1]);
		def<double>(d, "m",		y_[2]);
		def<double>(d, "h",		y_[3]);
		def<double>(d, "a",		y_[4]);
		def<double>(d, "b",		y_[5]);
		def<double>(d, "c",		y_[6]);
		def<double>(d, "d1",	y_[7]);
		def<double>(d, "d2",	y_[8]);
		def<double>(d, "p",		y_[9]);
		def<double>(d, "q",		y_[10]);
		def<double>(d, "r",		y_[11]);
		def<double>(d, "Ca",		y_[12]);
  
		def<double>(d, "da2_r",	  da2_r_);
		def<double>(d, "da2_s",	  da2_s_);
	}

	void nest::stn_cond_neuron6::State_::set(const DictionaryDatum& d, const Parameters_& p)
	{
		updateValue<double>(d, "V",  	  y_[0]);
		updateValue<double>(d, "n",		  y_[1]);
		updateValue<double>(d, "m",		  y_[2]);
		updateValue<double>(d, "h",		  y_[3]);
		updateValue<double>(d, "a",		  y_[4]);
		updateValue<double>(d, "b",		  y_[5]);
		updateValue<double>(d, "c",		  y_[6]);
		updateValue<double>(d, "d1",	  y_[7]);
		updateValue<double>(d, "d2",	  y_[8]);
		updateValue<double>(d, "p",		  y_[9]);
		updateValue<double>(d, "q",		  y_[10]);
		updateValue<double>(d, "r",		  y_[11]);
		updateValue<double>(d, "Ca",	  y_[12]);
		 
		 updateValue<double>(d, "da2_r",	  da2_r_);
		 updateValue<double>(d, "da2_s",	  da2_s_);
	}

	nest::stn_cond_neuron6::Buffers_::Buffers_(stn_cond_neuron6& n)
	 :	logger_(n),
		s_(0),
		c_(0),
		e_(0)
	{
	}

	nest::stn_cond_neuron6::Buffers_::Buffers_(const Buffers_&, stn_cond_neuron6& n)
     :	logger_(n),
		s_(0),
		c_(0),
		e_(0)
	{
	}

	/* ---------------------------------------------------------------- 
	 * Default and copy constructor for node, and destructor
	 * ---------------------------------------------------------------- */

	nest::stn_cond_neuron6::stn_cond_neuron6()
	 :	Archiving_Node(), 
		P_(), 
		S_(P_),
		B_(*this)
	{
		recordablesMap_.create();
	}

	nest::stn_cond_neuron6::stn_cond_neuron6(const stn_cond_neuron6& n)
	 :	Archiving_Node(n), 
		P_(n.P_), 
		S_(n.S_),
		B_(n.B_, *this)
	{
	}

	nest::stn_cond_neuron6::~stn_cond_neuron6()
	{
		// GSL structs only allocated by init_nodes_(), so we need to protect destruction
		if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
		if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
		if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
	}

	/* ---------------------------------------------------------------- 
	 * Node initialization functions
	 * ---------------------------------------------------------------- */

	void nest::stn_cond_neuron6::init_state_(const nest::Node& proto)
	{
		const stn_cond_neuron6& pr = downcast<stn_cond_neuron6>(proto);
		S_ = pr.S_;
	}

	void nest::stn_cond_neuron6::init_buffers_()
	{
		B_.spike_ampa_.clear();       // includes resize
		B_.spike_gabaa_.clear();      
		B_.spike_gabab_.clear();
		B_.spike_nmda_.clear(); 
		B_.spike_dopamine_.clear();

		B_.currents_.clear();
    
		Archiving_Node::clear_history();

		B_.logger_.reset();

		B_.step_ = Time::get_resolution().get_ms();
		B_.IntegrationStep_ = B_.step_;

		static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
		if ( B_.s_ == 0 )
			B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);	  // 13-dim system
		else 
			gsl_odeiv_step_reset(B_.s_);
    
		if ( B_.c_ == 0 )  
			B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);		  // abs, rel error
		else
			gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
		if ( B_.e_ == 0 )  
			B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);	  // 13-dim system
		else 
			gsl_odeiv_evolve_reset(B_.e_);
  
		B_.sys_.function  = stn_cond_neuron6_dynamics; 
		B_.sys_.jacobian  = NULL;
		B_.sys_.dimension =  State_::STATE_VEC_SIZE;					  // 13-dim system
		B_.sys_.params    = reinterpret_cast<void*>(this);

		B_.I_stim_ = 0.0;
	}

	void nest::stn_cond_neuron6::calibrate()
	{
		B_.logger_.init();    // ensure initialization in case mm connected after Simulate

		const double h = Time::get_resolution().get_ms(); 
  
		V_.delta_ampa_	= numerics::e/ P_.tau_syn_ampa;
		
		
		V_.delta_gabaa_	= 1.0/(std::pow(P_.tau_syn_gabaa2/P_.tau_syn_gabaa1, P_.tau_syn_gabaa_rise/P_.tau_syn_gabaa1) -
								std::pow(P_.tau_syn_gabaa2/P_.tau_syn_gabaa1, P_.tau_syn_gabaa_rise/P_.tau_syn_gabaa2));
		
		V_.delta_gabab_	= numerics::e/ P_.tau_syn_gabab;
  
		double_t	tau_d = 1/P_.da2_k_2_;
		double_t	tau_r = 1/P_.da2_k_4_;
	
		V_.da2_exp_A_00_	=  std::exp(-h / tau_d);
		V_.da2_exp_A_10_	= (std::exp(-h / tau_d) - std::exp(-h / tau_r)) / (1 / tau_r - 1 / tau_d);
		V_.da2_exp_A_11_	=  std::exp(-h / tau_r);
	
		double_t	da2_tau_r_add_ = 1 / (P_.da2_k_1_ * P_.da2_T_ + P_.da2_k_2_);
		V_.da2_r_inf_		= P_.da2_k_1_ * P_.da2_k_3_ * P_.da2_T_ * da2_tau_r_add_;
		V_.da2_exp_r_add_ = 1-std::exp(-P_.da2_pulse_ / da2_tau_r_add_);
	
	}

	/* ---------------------------------------------------------------- 
	 * Update and spike handling functions
	 * ---------------------------------------------------------------- */

	void nest::stn_cond_neuron6::update(Time const & origin, const long_t from, const long_t to)
	{
   
		assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
		assert(from < to);

		for ( long_t lag = from ; lag < to ; ++lag )
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
					throw GSLSolverFailure(get_name(), status);
			}
		
			/**	update synaptic conductances **/
			/** ampa **/
			S_.y_[State_::G_AMPA_A1] += V_.delta_ampa_ * B_.spike_ampa_.get_value(lag);
			S_.I_AMPA_ = S_.y_[State_::G_AMPA_A2] * ( S_.y_[State_::V_M] - P_.V_E );
	
			/** gaba **/
			if (!P_.inactivate_GABAA_synapse)
			    S_.y_[State_::G_GABAA_A1] += V_.delta_gabaa_ * B_.spike_gabaa_.get_value(lag);
			S_.y_[State_::G_GABAB_A1] += V_.delta_gabab_ * B_.spike_gabab_.get_value(lag);
			
			S_.I_GABAA_ = S_.y_[State_::G_GABAA_A2] * ( S_.y_[State_::V_M] - P_.V_GABAA );
			S_.I_GABAB_ = S_.y_[State_::G_GABAB_A2] * ( S_.y_[State_::V_M] - P_.V_GABAB );
			
			
			S_.I_K =	P_.g_K * pow(S_.y_[State_::HH_N], 4.0)*(S_.y_[State_::V_M] - P_.V_K);
			S_.I_Na =	P_.g_Na * pow( S_.y_[State_::HH_M], 3.0) * S_.y_[State_::HH_H] * (S_.y_[State_::V_M] - P_.V_Na);
			S_.I_A =	P_.g_IA * pow(S_.y_[State_::HH_A], 2.0) * S_.y_[State_::HH_B] * (S_.y_[State_::V_M] - P_.V_K);
			S_.I_LCa =		pow(S_.y_[State_::HH_C], 2.0)	* S_.y_[State_::HH_D1] * S_.y_[State_::HH_D2];
			S_.I_T =		pow(S_.y_[State_::HH_P], 2.0) * S_.y_[State_::HH_Q] ;
			S_.I_CaK =		P_.g_CaK * (S_.y_[State_::V_M] - P_.V_K) * pow(S_.y_[State_::HH_R], 4.0);

			/** dopamine **/
			S_.da2_r_	  += (V_.da2_r_inf_ - S_.da2_r_) * V_.da2_exp_r_add_ * B_.spike_dopamine_.get_value(lag); 
			S_.da2_s_	   = V_.da2_exp_A_10_ * S_.da2_r_ + V_. da2_exp_A_11_ * S_.da2_s_;
			S_.da2_r_	   = V_.da2_exp_A_00_ * S_.da2_r_;
			S_.da2_g_out_  = P_.da2_g_max_ /(P_.da2_k_d_ / std::pow(S_.da2_s_, P_.da2_n_) + 1);
			
			
			S_.y_[State_::G_NMDA_EXP] += B_.spike_nmda_.get_value(lag);
			S_.I_NMDA_ = S_.y_[State_::G_NMDA_EXP] * (S_.y_[State_::V_M] - P_.V_E) /(1.0 + P_.MgConc * std::exp(-0.062 * S_.y_[State_::V_M])/3.57);
			 
						
			/** dopaminergic modulation **/
			if ( P_.da2_Ca_modulation == 0.0)
			{
				if (S_.da2_g_out_ >= P_.da2_suppress_f_th_)	
				{
					P_.da2_suppress_f_ = P_.da2_suppress_f_const_;
				}	
			}
	
			// detect threshold crossing
			if ( S_.y_[State_::V_M] < P_.V_th_ )
			{
				S_.just_after_spike_ = 0;
			}
			else
				if ( S_.y_[State_::V_M] >= P_.V_th_ && S_.just_after_spike_ == 0)
				{
					set_spiketime(Time::step(origin.get_steps()+lag+1));
	  
					SpikeEvent se;
					network()->send(*this, se, lag);
		  
					S_.just_after_spike_ = 1 ;
				}

			//log state data
			B_.logger_.record_data(origin.get_steps() + lag);

			// set new input current
			B_.I_stim_ = B_.currents_.get_value(lag);
   
		}
	}
	
	
	port nest::stn_cond_neuron6::connect_sender(SpikeEvent& e, port receptor_type)
	{
		if (receptor_type <= 0)
			throw UnknownReceptorType(receptor_type, get_name());
			
		return receptor_type;
	}
	

	void nest::stn_cond_neuron6::handle(SpikeEvent & e)
	{
		assert(e.get_delay() > 0);
		
		//long rport = e.get_rport();

		if(e.get_rport() == 1)	/** AMPA **/
			B_.spike_ampa_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
									 e.get_weight() * e.get_multiplicity() );
			 
		if(e.get_rport() == 2)	/** NMDA **/
			B_.spike_nmda_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
									 e.get_weight() * e.get_multiplicity() );
  
		if(e.get_rport() == 3)	/** GABA-A **/
			B_.spike_gabaa_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
									 e.get_weight() * e.get_multiplicity() ); 

		if(e.get_rport() == 4)	/** GABA-B **/
			B_.spike_gabab_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
									  e.get_weight() * e.get_multiplicity() );
  
		if(e.get_rport() == 5)
			B_.spike_dopamine_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
									e.get_weight() * e.get_multiplicity() ); 
	}

	void nest::stn_cond_neuron6::handle(CurrentEvent& e)
	{
		assert(e.get_delay() > 0);

		const double_t c=e.get_current();
		const double_t w=e.get_weight();

		// add weighted current; HEP 2002-10-04
		B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w *c);
	}

	void nest::stn_cond_neuron6::handle(DataLoggingRequest& e)
	{
		B_.logger_.handle(e);
	}
	
} // namespace nest	

#endif //HAVE_GSL
 