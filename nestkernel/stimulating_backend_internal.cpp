/*
 *  stimulating_backend_internal.cpp
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

// C++ includes:
#include <iostream>


// Includes from nestkernel:
#include "stimulating_backend_internal.h"
#include "input_device.h"

void
nest::StimulatingBackendInternal::enroll( InputDevice& device, const DictionaryDatum& params )
{
  thread tid = device.get_thread();
  index node_id = device.get_node_id();

  auto device_it = devices_[ tid ].find( node_id );
  if ( device_it != devices_[ tid ].end() )
  {
    devices_[ tid ].erase( device_it );
  }
  devices_[ tid ].insert( std::make_pair( node_id, &device ) );
}

void
nest::StimulatingBackendInternal::disenroll( InputDevice& device )
{
  thread tid = device.get_thread();
  index node_id = device.get_node_id();

  auto device_it = devices_[ tid ].find( node_id );
  if ( device_it != devices_[ tid ].end() )
  {
    devices_[ tid ].erase( device_it );
  }
}


void
nest::StimulatingBackendInternal::initialize()
{
  device_map devices( kernel().vp_manager.get_num_threads() );
  devices_.swap( devices );
}

void
nest::StimulatingBackendInternal::prepare()
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::cleanup()
{
  // nothing to do
}


void
nest::StimulatingBackendInternal::finalize()
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::set_value_names( const InputDevice& device,
  const std::vector< Name >& double_value_names,
  const std::vector< Name >& long_value_names )
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::check_device_status( const DictionaryDatum& params ) const
{
  // nothing to do
}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */
void
nest::StimulatingBackendInternal::get_status( DictionaryDatum& d ) const
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::set_status( const DictionaryDatum& d )
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::pre_run_hook()
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::post_run_hook()
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::post_step_hook()
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::get_device_defaults( DictionaryDatum& params ) const
{
  // nothing to do
}

void
nest::StimulatingBackendInternal::get_device_status( const nest::InputDevice& device,
  DictionaryDatum& params_dictionary ) const
{
  // nothing to do
}