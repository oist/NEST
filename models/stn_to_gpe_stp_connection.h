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

#include "connection.h"

/* BeginDocumentation
Name: 

Description:



References:
[1]	
*/


namespace nest 
{
    template < typename targetidentifierT >
	class STNtoGPeSTPConnection  : public Connection< targetidentifierT >
    {
	public:
	    typedef CommonSynapseProperties CommonPropertiesType;
	    typedef Connection<targetidentifierT> ConnectionBase;
	    /**
	     * Constructor.
	     * Sets default values for all parameters. Needed by GenericConnectorModel.
	     */
	    STNtoGPeSTPConnection();

	    /**
	     * Default Destructor.
	     */
	    ~STNtoGPeSTPConnection() {}
	    
	    /**
	     * Copy constructor 
	     */
	    STNtoGPeSTPConnection( const STNtoGPeSTPConnection& );

	    using ConnectionBase::get_delay_steps;
	    using ConnectionBase::get_delay;
	    using ConnectionBase::get_rport;
	    using ConnectionBase::get_target;

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
//	    //void set_status(const DictionaryDatum & d, index p, ConnectorModel &cm);

	    /**
	     * Create new empty arrays for the properties of this connection in the given
	     * dictionary. It is assumed that they are not existing before.
	     */
//	    void initialize_property_arrays(DictionaryDatum & d) const;

	    /**
	     * Append properties of this connection to the given dictionary. If the
	     * dictionary is empty, new arrays are created first.
	     */
//	    void append_properties(DictionaryDatum & d) const;

	    /**
	     * Send an event to the receiver of this connection.
	     * \param e The event to send
	     * \param t_lastspike Point in time of last spike sent.
	     * \param cp Common properties to all synapses (empty).
	     */
	    void send(Event& e, thread t, double_t t_lastspike, const CommonSynapseProperties& cp);

	    class ConnTestDummyNode: public ConnTestDummyNodeBase
	    {
		public:
		    using ConnTestDummyNodeBase::handles_test_event;
		    port handles_test_event(SpikeEvent&, rport) { return invalid_port_; }
	    };

	    void check_connection( Node& s, Node& t, rport receptor_type, double_t, const CommonPropertiesType& )
	    {
		ConnTestDummyNode dummy_target;
		ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
	    }

	    void set_weight(double_t w) { weight_ = w; }

	    // overloaded for all supported event types
//	    using Connection::check_event;
//	    void check_event(SpikeEvent&) {}

	private:
	    double_t weight_;
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
    template < typename targetidentifierT >
    void STNtoGPeSTPConnection< targetidentifierT >::send(Event& e, thread t, double_t t_lastspike, const CommonSynapseProperties &)
    {
	Node *target = get_target(t);
	double_t delay = get_delay();	
	
	double_t h = e.get_stamp().get_ms() - t_lastspike;

	double_t Pf = std::exp(-h/tau_fac);
	double_t Pd = std::exp(-h/tau_dep);

	f_ = f0_ + (f_ - f0_)*Pf;
	inc_fac_bound = 1.0 + (inc_fac - 1.0)*(f_bound - f_)/(f_bound -1.0); 
	f_ *= inc_fac_bound;

	d_ = d0_ + (d_ - d0_)*Pd;
	d_ *= inc_dep;
    

	// send the spike to the target
	e.set_receiver(*target);
	e.set_weight( weight_ *f_ *d_ );
	e.set_delay(get_delay_steps());
	e.set_rport(get_rport());
	e();

    }

    template < typename targetidentifierT >
    STNtoGPeSTPConnection < targetidentifierT >::STNtoGPeSTPConnection() :
	ConnectionBase(),
	weight_(1.0),
	tau_fac(148.0),
	inc_fac(1.64),
	f0_(1.0),
	f_(1.0),
	f_bound(5.0),
	tau_dep(764.0),
	inc_dep(0.55),
	d0_(1.0),
	d_(1.0)
    { }

    template < typename targetidentifierT >
	STNtoGPeSTPConnection< targetidentifierT >::STNtoGPeSTPConnection( const STNtoGPeSTPConnection& rhs )
	: ConnectionBase( rhs ),
	weight_(weight_),
	tau_fac(tau_fac),
	inc_fac(inc_fac),
	f0_(f0_),
	f_(f_),
	f_bound(f_bound),
	tau_dep(tau_dep),
	inc_dep(inc_dep),
	d0_(d0_),
	d_(d_)
    { }

    template < typename targetidentifierT > void
	STNtoGPeSTPConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
	{
	    ConnectionBase::get_status(d);
	    def< double_t >( d, names::weight, weight_ );
	    def<double_t>(d, "tau_fac", tau_fac);
	    def<double_t>(d, "inc_fac", inc_fac);
	    def<double_t>(d, "f0", f0_);
	    def<double_t>(d, "f_bound", f_bound);
	    def<double_t>(d, "f", f_);
	    def<double_t>(d, "tau_dep", tau_dep);
	    def<double_t>(d, "inc_dep", d_);
	    def<double_t>(d, "d0", d0_);
	    def<double_t>(d, "d", d_);	  
	}

    template < typename targetidentifierT >
	void STNtoGPeSTPConnection< targetidentifierT >::set_status(const DictionaryDatum & d, ConnectorModel &cm)
	{
	    ConnectionBase::set_status(d, cm);
	    updateValue< double_t >( d, names::weight, weight_ );

	    updateValue<double_t>(d, "tau_fac", tau_fac);
	    updateValue<double_t>(d, "inc_fac", inc_fac);
	    updateValue<double_t>(d, "f0", f0_);
	    updateValue<double_t>(d, "f_bound", f_bound);
	    updateValue<double_t>(d, "f", f_);
	    updateValue<double_t>(d, "tau_dep", tau_dep);
	    updateValue<double_t>(d, "inc_dep", inc_dep);
	    updateValue<double_t>(d, "d0", d0_);
	    updateValue<double_t>(d, "d", d_);
	}






} // namespace nest






#endif // STN_TO_GPE_STP_CONNECTION_H
