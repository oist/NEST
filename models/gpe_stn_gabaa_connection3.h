/*
 *  gpe_stn_gabaa_connection3.h
 *
 *	Created by Osamu Shouno	2013/07/01
 *	Copyright 2013 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *	2015-12-11 Converted to NEST 2.8  Jan Moren
 *
 *	Permission is granted to compile and modify
 *	this file for non-commercial use.
 *	See the file LICENSE for details.
 *
 */

#ifndef GPE_STN_GABAA_CONNECTION3_H
#define GPE_STN_GABAA_CONNECTION3_H

#include "connection.h"

/* BeginDocumentation

*/


namespace nest 
{

    template < typename targetidentifierT >
	class GPeSTN_GABAA_Connection3  : public Connection< targetidentifierT >
    {
	public:
	    typedef CommonSynapseProperties CommonPropertiesType;
	    typedef Connection<targetidentifierT> ConnectionBase;

	    /**
	     * Default Constructor.
	     * Sets default values for all parameters. Needed by GenericConnectorModel.
	     */
	    GPeSTN_GABAA_Connection3();

	    /**
	     * Default Destructor.
	     */
	    ~GPeSTN_GABAA_Connection3() {}

	    GPeSTN_GABAA_Connection3( const GPeSTN_GABAA_Connection3& );

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
	     * Send an event to the receiver of this connection.
	     * \param e The event to send
	     * \param t thread ID
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

	private:

	    double_t weight_;
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
     * \param t thread ID
     * \param p The port under which this connection is stored in the Connector.
     * \param t_lastspike Time point of last spike emitted
     */
    template < typename targetidentifierT >
    inline void GPeSTN_GABAA_Connection3< targetidentifierT >::send(Event& e, thread t, double_t t_lastspike, const CommonSynapseProperties &)
	{
	    Node *target = get_target(t);
	    double_t delay = get_delay();	

	    double_t h = e.get_stamp().get_ms() - t_lastspike;
	    double_t Pd_ = std::exp(-h/tau_dep);

	    d_ =  d0_ + (d_ - d0_)*Pd_;
	    d_ *= 1.0 + (fd_ - 1.0)*std::pow((d_ - d_bound)/(1.0 - d_bound), n_bound);

	    if (d_ < d_bound)
		d_ = d_bound;

	    myu_ = 1.0 - std::pow(1.0 - d_, 1.0/n_contacts_mean);

	    int count = 0;

	    for (int i =0; i < n_contacts_; i ++){
		librandom::RngPtr rng = NestModule::get_network().get_rng(t);
		if (rng->drand() < myu_)
		    count +=1;
	    }

	    // send the spike to the target
	    e.set_receiver(*target);
	    e.set_weight( weight_ * count);
	    e.set_delay(get_delay_steps());
	    e.set_rport(get_rport());
	    e();
	}

    template < typename targetidentifierT >
	GPeSTN_GABAA_Connection3 < targetidentifierT >::GPeSTN_GABAA_Connection3() :
	    ConnectionBase(),
	    weight_(1.0),
	    n_contacts_(1.),
	    n_contacts_mean(15.),
	    d_(1.0),
	    d0_(1.0),
	    d_bound(0.0),
	    tau_dep(20000.0),
	    fd_(0.97),
	    n_bound(1.0),
	    myu_(1.0)
    { }
    
    template < typename targetidentifierT >
	GPeSTN_GABAA_Connection3< targetidentifierT >::GPeSTN_GABAA_Connection3( const GPeSTN_GABAA_Connection3& rhs )
	: ConnectionBase( rhs ),
	weight_(rhs.weight_),
	n_contacts_(rhs.n_contacts_),
	n_contacts_mean(rhs.n_contacts_mean),
	d_(rhs.d_),
	d0_(rhs.d0_),
	d_bound(rhs.d_bound),
	tau_dep(rhs.tau_dep),
	fd_(rhs.fd_),
	n_bound(rhs.n_bound),
	myu_(rhs.myu_)
    { }
    
    template < typename targetidentifierT > void
	GPeSTN_GABAA_Connection3< targetidentifierT >::get_status( DictionaryDatum& d ) const
	{
	    ConnectionBase::get_status(d);
	    def< double_t >( d, names::weight, weight_ );
	    def<double_t>(d, "n_contacts", n_contacts_);
	    def<double_t>(d, "n_contacts_mean", n_contacts_mean);  
	    def<double_t>(d, "d", d_);
	    def<double_t>(d, "d0", d0_);
	    def<double_t>(d, "d_bound", d_bound); 
	    def<double_t>(d, "n_bound", n_bound);
	    def<double_t>(d, "tau_dep", tau_dep);
	    def<double_t>(d, "fd", fd_);
	    def<double_t>(d, "myu", myu_); 
	}
    
    template < typename targetidentifierT >
	void GPeSTN_GABAA_Connection3< targetidentifierT >::set_status(const DictionaryDatum & d, ConnectorModel &cm)
	{
	    ConnectionBase::set_status(d, cm);
	    updateValue< double_t >( d, names::weight, weight_ );
	    updateValue<double_t>(d, "n_contacts", n_contacts_);
	    updateValue<double_t>(d, "n_contacts_mean", n_contacts_mean);  
	    updateValue<double_t>(d, "d", d_);
	    updateValue<double_t>(d, "d0", d0_);
	    updateValue<double_t>(d, "d_bound", d_bound);
	    updateValue<double_t>(d, "n_bound", n_bound);
	    updateValue<double_t>(d, "tau_dep", tau_dep);
	    updateValue<double_t>(d, "fd", fd_);

	    updateValue<double_t>(d, "myu", myu_);
	}

} // namespace nest

#endif // GPE_STN_GABAA_CONNECTION3_H
