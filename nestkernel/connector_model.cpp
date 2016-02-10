/*
 *  connector_model.cpp
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

#include "network.h"
#include "connector_model.h"
#include "connector_base.h"

namespace nest
{

ConnectorModel::ConnectorModel( Network& net,
  const std::string name,
  bool is_primary,
  bool has_delay )
  : net_( net )
  , min_delay_( Time::pos_inf() )
  , max_delay_( Time::neg_inf() )
  , num_connections_( 0 )
  , default_delay_needs_check_( true )
  , user_set_delay_extrema_( false )
  , name_( name )
  , is_primary_( is_primary )
  , has_delay_( has_delay )
{
}

ConnectorModel::ConnectorModel( const ConnectorModel& cm, const std::string name )
  : net_( cm.net_ )
  , min_delay_( cm.min_delay_ )
  , max_delay_( cm.max_delay_ )
  , num_connections_( 0 )
  , default_delay_needs_check_( true )
  , user_set_delay_extrema_( cm.user_set_delay_extrema_ )
  , name_( name )
  , is_primary_( cm.is_primary_ )
  , has_delay_( cm.has_delay_ )
{
  min_delay_.calibrate(); // in case of change in resolution
  max_delay_.calibrate();
}

void
ConnectorModel::assert_valid_delay_ms( double_t requested_new_delay )
{
  const delay new_delay = Time::delay_ms_to_steps( requested_new_delay );
  const double new_delay_ms = Time::delay_steps_to_ms( new_delay );

  if ( new_delay < Time::get_resolution().get_steps() )
    throw BadDelay( new_delay_ms, "Delay must be greater than or equal to resolution" );

  // if already simulated, the new delay has to be checked against the
  // min_delay and the max_delay which have been used during simulation
  if ( net_.get_simulated() )
  {
    const bool bad_min_delay = new_delay < net_.get_min_delay();
    const bool bad_max_delay = new_delay > net_.get_max_delay();

    if ( bad_min_delay || bad_max_delay )
      throw BadDelay( new_delay_ms,
        "Minimum and maximum delay cannot be changed "
        "after Simulate has been called." );
  }

  const bool new_min_delay = new_delay < min_delay_.get_steps();
  const bool new_max_delay = new_delay > max_delay_.get_steps();

  if ( new_min_delay )
  {
    if ( user_set_delay_extrema_ )
    {
      throw BadDelay( new_delay_ms,
        "Delay must be greater than or equal to min_delay. "
        "You may set min_delay before creating connections." );
    }
    else
    {
      min_delay_ = Time( Time::step( new_delay ) );
    }
  }

  if ( new_max_delay )
  {
    if ( user_set_delay_extrema_ )
    {
      throw BadDelay( new_delay_ms,
        "Delay must be smaller than or equal to max_delay. "
        "You may set min_delay before creating connections." );
    }
    else
    {
      max_delay_ = Time( Time::step( new_delay ) );
    }
  }
}

void
ConnectorModel::assert_two_valid_delays_steps( delay new_delay1, delay new_delay2 )
{
  const delay ldelay = std::min( new_delay1, new_delay2 );
  const delay hdelay = std::max( new_delay1, new_delay2 );

  if ( ldelay < Time::get_resolution().get_steps() )
    throw BadDelay(
      Time::delay_steps_to_ms( ldelay ), "Delay must be greater than or equal to resolution" );

  if ( net_.get_simulated() )
  {
    const bool bad_min_delay = ldelay < net_.get_min_delay();
    const bool bad_max_delay = hdelay > net_.get_max_delay();

    if ( bad_min_delay )
      throw BadDelay( Time::delay_steps_to_ms( ldelay ),
        "Minimum delay cannot be changed after Simulate has been called." );

    if ( bad_max_delay )
      throw BadDelay( Time::delay_steps_to_ms( hdelay ),
        "Maximum delay cannot be changed after Simulate has been called." );
  }

  const bool new_min_delay = ldelay < min_delay_.get_steps();
  const bool new_max_delay = hdelay > max_delay_.get_steps();

  if ( new_min_delay )
  {
    if ( user_set_delay_extrema_ )
    {
      throw BadDelay( Time::delay_steps_to_ms( ldelay ),
        "Delay must be greater than or equal to min_delay. "
        "You may set min_delay before creating connections." );
    }
    else
    {
      min_delay_ = Time( Time::step( ldelay ) );
    }
  }

  if ( new_max_delay )
  {
    if ( user_set_delay_extrema_ )
    {
      throw BadDelay( Time::delay_steps_to_ms( hdelay ),
        "Delay must be smaller than or equal to max_delay. "
        "You may set max_delay before creating connections." );
    }
    else
    {
      max_delay_ = Time( Time::step( hdelay ) );
    }
  }
}

} // namespace nest
