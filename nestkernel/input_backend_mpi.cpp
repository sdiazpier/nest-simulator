/*
 *  input_backend_mpi.cpp
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
#include "input_backend_mpi.h"
#include "input_device.h"


void
nest::InputBackendMPI::initialize()
{
  auto nthreads = kernel().vp_manager.get_num_threads();
  device_map devices( nthreads );
  devices_.swap( devices );
}

void
nest::InputBackendMPI::finalize()
{
  // clear vector of map
  for ( auto& it_device : devices_ )
  {
    it_device.clear();
  }
  devices_.clear();
  commMap_.clear();
}

void
nest::InputBackendMPI::enroll( InputDevice& device, const DictionaryDatum& params )
{
  if ( device.get_type() == InputDevice::SPIKE_GENERATOR or device.get_type() == InputDevice::STEP_CURRENT_GENERATOR )
  {
    thread tid = device.get_thread();
    index node_id = device.get_node_id();

    auto device_it = devices_[ tid ].find( node_id );
    if ( device_it != devices_[ tid ].end() )
    {
      devices_[ tid ].erase( device_it );
    }
    std::pair< MPI_Comm*, InputDevice* > pair = std::make_pair( nullptr, &device );
    devices_[ tid ].insert( std::make_pair( node_id, pair ) );
  }
  else
  {
    throw BadProperty( "Only spike generators can have input backend 'mpi'." );
  }
}

void
nest::InputBackendMPI::disenroll( InputDevice& device )
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
nest::InputBackendMPI::set_value_names( const InputDevice& device,
  const std::vector< Name >& double_value_names,
  const std::vector< Name >& long_value_names )
{
  // nothing to do
}

void
nest::InputBackendMPI::prepare()
{
  // need to be run only by the master thread : it is the case because it's not run in parallel
  thread thread_id_master = kernel().vp_manager.get_thread_id();
  // Create the connection with MPI
  // 1) take all the ports of the connections
  // get port and update the list of device only for mastecommr
  for ( auto& it_device : devices_[ thread_id_master ] ) {
    // add the link between MPI communicator and the device (devices can share the same MPI communicator)
    std::string port_name;
    get_port( it_device.second.second, &port_name );
    auto comm_it = commMap_.find( port_name );
    MPI_Comm* comm;
    std::vector<int>* vector_id_device;
    if ( comm_it != commMap_.end() ) {
      comm = comm_it->second.first;
      comm_it->second.second->push_back(it_device.second.second->get_node_id());
    } else {
      comm = new MPI_Comm;
      vector_id_device = new std::vector<int>;
      vector_id_device->push_back(it_device.second.second->get_node_id());
      std::pair< MPI_Comm *, std::vector<int>* > comm_count = std::make_pair( comm, vector_id_device );
      commMap_.insert( std::make_pair( port_name, comm_count ) );
    }
    it_device.second.first = comm;
  }

  // 2) connect the master thread to the MPI process it needs to be connected to
  for ( auto& it_comm : commMap_ ) {
    MPI_Comm_connect(it_comm.first.data(),
                     MPI_INFO_NULL,
                     0,
                     MPI_COMM_WORLD,
                     it_comm.second.first); // should use the status for handle error
    std::ostringstream msg;
    msg << "Connect to " << it_comm.first.data() << "\n";
    LOG(M_INFO, "MPI Input connect", msg.str());
  }
}

void
nest::InputBackendMPI::pre_run_hook()
{
  #pragma omp master
  {
    for ( auto& it_comm : commMap_ )
    {
      bool value [ 1 ]  = { true } ;
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 0, *it_comm.second.first );
      receive_spike_train(*it_comm.second.first,*it_comm.second.second);
    }
  }
  #pragma omp barrier
}

void
nest::InputBackendMPI::post_step_hook()
{
  // nothing to do
}

void
nest::InputBackendMPI::post_run_hook()
{
  #pragma omp master
  {
    // Send information about the end of the running part
    for ( auto& it_comm : commMap_ )
    {
      bool value [ 1 ]  = { true } ;
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 1, *it_comm.second.first );
    }
  }
  #pragma omp barrier
}

void
nest::InputBackendMPI::cleanup()
{
  // Disconnect all the MPI connection and send information about this disconnection
  // Clean all the elements in the map
  // disconnect MPI message
  #pragma omp master
  {
    for ( auto& it_comm : commMap_ ) {
      bool value[ 1 ] = { true };
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 2, *it_comm.second.first );
      MPI_Comm_disconnect( it_comm.second.first );
      delete it_comm.second.first;
    }
    // clear map of devices
    commMap_.clear();
    thread thread_id_master = kernel().vp_manager.get_thread_id();
    for ( auto& it_device : devices_[thread_id_master] ) {
      it_device.second.first = nullptr;
    }
  }
  #pragma omp barrier
}

void
nest::InputBackendMPI::check_device_status( const DictionaryDatum& params ) const
{
  // nothing to do
}

void
nest::InputBackendMPI::get_device_defaults( DictionaryDatum& params ) const
{
  // nothing to do
}

void
nest::InputBackendMPI::get_device_status( const nest::InputDevice& device, DictionaryDatum& params_dictionary ) const
{
  // nothing to do
}


void
nest::InputBackendMPI::get_status( lockPTRDatum< Dictionary, &SLIInterpreter::Dictionarytype >& ) const
{
  // nothing to do
}

void
nest::InputBackendMPI::set_status( const DictionaryDatum& d )
{
  // nothing to do
}


void
nest::InputBackendMPI::get_port( InputDevice* device, std::string* port_name )
{
  get_port( device->get_node_id(), device->get_label(), port_name );
}

void
nest::InputBackendMPI::get_port( const index index_node, const std::string& label, std::string* port_name )
{
  // path of the file : path+label+id+.txt
  // (file contains only one line with name of the port)
  std::ostringstream basename;
  const std::string& path = kernel().io_manager.get_data_path();
  if ( not path.empty() )
  {
    basename << path << '/';
  }
  basename << kernel().io_manager.get_data_prefix();

  if ( not label.empty() )
  {
    basename << label;
  }
  else
  {
    throw MPIFilePortsUnknown( index_node );
  }
  char add_path[ 150 ];
  sprintf( add_path, "/%zu.txt", index_node );
  basename << add_path;
  std::cout << basename.rdbuf() << std::endl;
  std::ifstream file( basename.str() );

  if ( file.is_open() )
  {
    getline( file, *port_name );
  }
  file.close();
}

void
nest::InputBackendMPI::receive_spike_train( const MPI_Comm& comm, std::vector<int>& devices_id )
{
  // Send size of the list id
  int size_list[1];
  size_list [0] = devices_id.size();
  MPI_Send( &size_list, 1, MPI_INT, 0, 0, comm );
  // Send the list of device ids
  MPI_Send( &devices_id[0], size_list[0], MPI_INT, 0, 0, comm );
  // Receive the size of data
  MPI_Status status_mpi;
  int* nb_size_data_per_id{ new int[ size_list[ 0 ] + 1 ]{} };
  MPI_Recv( nb_size_data_per_id, size_list[ 0 ] + 1, MPI_INT, MPI_ANY_SOURCE, devices_id[0], comm, &status_mpi );
  // Receive the data
  double* data{ new double[ nb_size_data_per_id[ 0 ] ]{} };
  MPI_Recv( data, nb_size_data_per_id[ 0 ], MPI_DOUBLE, status_mpi.MPI_SOURCE, devices_id[ 0 ], comm, &status_mpi );
  int index_data = 0; // first one is the total number spikes
  int index_id_device = 0;
  // update all the device with spikes
  for (auto &id : devices_id)
  {
    std::vector< double > data_for_device( &data[ index_data ], &data[ index_data + nb_size_data_per_id [ index_id_device+1 ] ] );

    // Update the device with the data in all the thread
    for (auto &thread_device : devices_) {
      thread_device.find( id )->second.second->update_from_backend( data_for_device );
    }
    index_id_device += 1;
    index_data += nb_size_data_per_id [ index_id_device ];
  }
  // clean the memory
  delete[] data;
  data = nullptr;
  delete[] nb_size_data_per_id;
  nb_size_data_per_id = nullptr;
}
