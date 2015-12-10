/*
 *  gpe_stn_gabaa_connection3.h
 *
 *	Created by Osamu Shouno	2013/07/01
 *	Copyright 2013 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *
 *	Permission is granted to compile and modify
 *	this file for non-commercial use.
 *	See the file LICENSE for details.
 *
 */

#ifndef GPE_STN_GABAA_CONNECTION3_H
#define GPE_STN_GABAA_CONNECTION3_H

#include "connection_het_wd.h"
#include "nestmodule.h"

/* BeginDocumentation
  
*/


namespace nest 
{

	class GPeSTN_GABAA_Connection3 : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  GPeSTN_GABAA_Connection3();

  /**
   * Default Destructor.
   */
  ~GPeSTN_GABAA_Connection3() {}

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

  double_t n_contacts_;
  double_t n_contacts_mean;
	
  double_t d_;
  double_t d0_;
  double_t d_bound;
  double_t tau_dep;
  double_t fd_;
  double_t n_bound;


  double_t myu_;
};

	
	
/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void GPeSTN_GABAA_Connection3::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
	double_t h = e.get_stamp().get_ms() - t_lastspike;
	double_t Pd_ = std::exp(-h/tau_dep);
	
	d_ =  d0_ + (d_ - d0_)*Pd_;
	d_ *= 1.0 + (fd_ - 1.0)*std::pow((d_ - d_bound)/(1.0 - d_bound), n_bound);
	
	if (d_ < d_bound)
		d_ = d_bound;
		
	myu_ = 1.0 - std::pow(1.0 - d_, 1.0/n_contacts_mean);
	
	int count = 0;
	
	for (int i =0; i < n_contacts_; i ++){
		librandom::RngPtr rng = NestModule::get_network().get_rng(target_->get_thread());
		if (rng->drand() < myu_)
			count +=1;
	}
	

	// send the spike to the target
	e.set_receiver(*target_);
	e.set_weight( weight_ * count);
	e.set_delay( delay_ );
	e.set_rport( rport_ );
	e();
}
 
} // namespace nest

#endif // GPE_STN_GABAA_CONNECTION3_H
