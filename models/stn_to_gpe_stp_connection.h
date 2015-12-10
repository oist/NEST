/*
 *  stn_to_gpe_stp_connection.h
 *
 *	Created by Osamu Shouno	2013/09/25
 *	Copyright 2011 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *
 *	Permission is granted to compile and modify
 *	this file for non-commercial use.
 *	See the file LICENSE for details.
 *
 */

#ifndef STN_TO_GPE_STP_CONNECTION_H
#define STN_TO_GPE_STP_CONNECTION_H

#include "connection_het_wd.h"

/* BeginDocumentation
  Name: 

  Description:
   
   

  References:
   [1]	
*/


namespace nest 
{

	class STNtoGPeSTPConnection : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STNtoGPeSTPConnection();

  /**
   * Default Destructor.
   */
  ~STNtoGPeSTPConnection() {}

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;
  
  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, ConnectorModel &cm);

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void set_status(const DictionaryDatum & d, index p, ConnectorModel &cm);

  /**
   * Create new empty arrays for the properties of this connection in the given
   * dictionary. It is assumed that they are not existing before.
   */
  void initialize_property_arrays(DictionaryDatum & d) const;

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void append_properties(DictionaryDatum & d) const;

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp Common properties to all synapses (empty).
   */
  void send(Event& e, double_t t_lastspike, const CommonSynapseProperties &cp);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}
 
 private:
  double_t tau_fac;   //!< [ms] time constant for fascilitation
  double_t inc_fac;
  double_t inc_fac_bound;	
  double_t f0_;
  double_t f_;
  double_t f_bound;
	
  double_t tau_dep;   //!< [ms] time constant for fascilitation
  double_t inc_dep;
  double_t d0_;
  double_t d_;	

};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void STNtoGPeSTPConnection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  double_t h = e.get_stamp().get_ms() - t_lastspike;

  double_t Pf = std::exp(-h/tau_fac);
  double_t Pd = std::exp(-h/tau_dep);
  
  f_ = f0_ + (f_ - f0_)*Pf;
  inc_fac_bound = 1.0 + (inc_fac - 1.0)*(f_bound - f_)/(f_bound -1.0); 
  f_ *= inc_fac_bound;
	
  d_ = d0_ + (d_ - d0_)*Pd;
  d_ *= inc_dep;
	
  
  // send the spike to the target
  e.set_receiver(*target_);
  e.set_weight( weight_ * f_ *d_ );
  e.set_delay( delay_ );
  e.set_rport( rport_ );
  e();

}
 
} // namespace nest

#endif // STN_TO_GPE_STP_CONNECTION_H
