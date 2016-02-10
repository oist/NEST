/*
 *  gap_junction.h
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


/* BeginDocumentation
Name: gap_junction - Synapse type for gap-junction connections.

Description:
 gap_junction is a connector to create gap junctions between pairs
 of neurons. Please note that gap junctions are two-way connections:
 In order to create an accurate gap-junction connection between two
 neurons i and j two connections are required:

 i j conn_spec gap_junction   Connect
 j i conn_spec gap_junction   Connect

 The value of the parameter "delay" is ignored for connections of
 type gap_junction.

Transmits: GapJEvent

References:

 Hahne, J., Helias, M., Kunkel, S., Igarashi, J.,
 Bolten, M., Frommer, A. and Diesmann, M.,
 A unified framework for spiking and gap-junction interactions
 in distributed neuronal network simulations,
 Front. Neuroinform. 9:22. (2015),
 doi: 10.3389/fninf.2015.00022

 Mancilla, J. G., Lewis, T. J., Pinto, D. J.,
 Rinzel, J., and Connors, B. W.,
 Synchronization of electrically coupled pairs
 of inhibitory interneurons in neocortex,
 J. Neurosci. 27, 2058-2073 (2007),
 doi: 10.1523/JNEUROSCI.2715-06.2007

Author: Jan Hahne, Moritz Helias, Susanne Kunkel
SeeAlso: synapsedict, hh_psc_alpha_gap
*/


#ifndef GAP_JUNCTION_H
#define GAP_JUNCTION_H

#include "connection.h"

namespace nest
{

/**
 * Class representing a gap-junction connection. A gap-junction connection
 * has the properties weight, delay and receiver port.
 */

template < typename targetidentifierT >
class GapJunction : public Connection< targetidentifierT >
{

  double_t weight_;

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;
  typedef GapJunctionEvent EventType;


  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  GapJunction()
    : ConnectionBase()
    , weight_( 1.0 )
  {
  }


  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;


  void
  check_connection( Node& s, Node& t, rport receptor_type, double_t, const CommonPropertiesType& )
  {
    EventType ge;

    s.sends_secondary_event( ge );
    ge.set_sender( s );
    Connection< targetidentifierT >::target_.set_rport( t.handles_test_event( ge, receptor_type ) );
    Connection< targetidentifierT >::target_.set_target( &t );
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted
   */
  void
  send( Event& e, thread t, double_t, const CommonSynapseProperties& )
  {
    e.set_weight( weight_ );
    e.set_delay( get_delay_steps() );
    e.set_receiver( *get_target( t ) );
    e.set_rport( get_rport() );
    e();
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  set_weight( double_t w )
  {
    weight_ = w;
  }
};

template < typename targetidentifierT >
void
GapJunction< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  // We have to include the delay here to prevent
  // errors due to internal calls of
  // this function in SLI/pyNEST
  ConnectionBase::get_status( d );
  def< double_t >( d, names::weight, weight_ );
  def< long_t >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
GapJunction< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  // If the delay is set, we throw a BadProperty
  if ( d->known( names::delay ) )
    throw BadProperty( "gap_junction connection has no delay" );

  ConnectionBase::set_status( d, cm );
  updateValue< double_t >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef GAP_JUNCTION_H */
