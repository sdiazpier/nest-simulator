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

    // for each thread, add the input device if it's not already in the map
    auto device_it = devices_[ tid ].find( node_id );
    if ( device_it != devices_[ tid ].end() )
    {
      devices_[ tid ].erase( device_it );
    }
    // the MPI communication will be initialise during the prepare function
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

  // remove the device from the map
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
  // need to be run only by the master thread : it is the case because this part is not running in parallel
  gettimeofday(&time_prepare_init[index_time],NULL);
  thread thread_id_master = kernel().vp_manager.get_thread_id();
  // Create the connection with MPI
  // 1) take all the ports of the connections
  // get port and update the list of device only for master
  for ( auto& it_device : devices_[ thread_id_master ] ) {
    // add the link between MPI communicator and the device (devices can share the same MPI communicator)
    std::string port_name;
    get_port( it_device.second.second, &port_name );
    auto comm_it = commMap_.find( port_name );
    MPI_Comm* comm;
    if ( comm_it != commMap_.end() ) {
      // it's not a new communicator
      comm = std::get<0>(comm_it->second);
      // add the id of the device if there are a connection with the device.
      if (kernel().connection_manager.get_device_connected(thread_id_master,it_device.second.second->get_local_device_id(),kernel().mpi_manager.get_rank())){
        std::get<1>(comm_it->second)->push_back(it_device.second.second->get_node_id());
        std::get<2>(comm_it->second)[thread_id_master]+=1;
      }
    } else {
      // create a new MPI communicator for the MPI process. Only the master thread is using MPI function.
      // The cause it's the management of the thread is MPI_THREAD_FUNNELED (mpi_manager 119).
      comm = new MPI_Comm;
      auto vector_id_device = new std::vector<int>; // vector of ID device for the rank
      int* vector_nb_device_th { new int[ kernel().vp_manager.get_num_threads() ] {} }; // number of device by thread
      std::fill_n(vector_nb_device_th,kernel().vp_manager.get_num_threads(),0);
      // add the id of the device if there are a connection with the device.
      if (kernel().connection_manager.get_device_connected(thread_id_master,it_device.second.second->get_local_device_id(),kernel().mpi_manager.get_rank())){
        vector_id_device->push_back(it_device.second.second->get_node_id());
        vector_nb_device_th[thread_id_master]+=1;
      }
      std::tuple< MPI_Comm *, std::vector<int>*, int* > comm_count = std::make_tuple( comm, vector_id_device, vector_nb_device_th );
      commMap_.insert( std::make_pair( port_name, comm_count ) );
    }
    it_device.second.first = comm;
  }

  // Add the id of device of the other thread in the vector_id_device and update the count of all device
  for (int id_thread = 0; id_thread<kernel().vp_manager.get_num_threads(); id_thread++)
  {
    // don't do it again for the master thread
    if ( id_thread != thread_id_master )
    {
      for ( auto& it_device : devices_[ id_thread ] )
      {
        // add the id of the device if there are a connection with the device.
        if ( kernel().connection_manager.get_device_connected(
               id_thread, it_device.second.second->get_local_device_id(), kernel().mpi_manager.get_rank() ) )
        {
          std::string port_name;
          get_port( it_device.second.second, &port_name );
          auto comm_it = commMap_.find( port_name );
          if ( comm_it != commMap_.end() )
          {
            std::get< 1 >( comm_it->second )->push_back( it_device.second.second->get_node_id() );
            std::get< 2 >( comm_it->second )[ id_thread ] += 1;
          }
          else
          {
            // should be impossible
            throw KernelException( "The MPI port was not define in the master thread" );
          }
        }
      }
    }
  }

  // 2) connect the master thread to the MPI process it needs to be connected to
  for ( auto& it_comm : commMap_ ) {
    MPI_Comm_connect(it_comm.first.data(),
                     MPI_INFO_NULL,
                     0,
                     MPI_COMM_WORLD,
                     std::get<0>(it_comm.second)); // should use the status for handle error
    std::ostringstream msg;
    msg << "Connect to " << it_comm.first.data() << "\n";
    LOG(M_INFO, "MPI Input connect", msg.str());
  }
  gettimeofday(&time_prepare_end[index_time],NULL);
}

void
nest::InputBackendMPI::pre_run_hook()
{
  // create the variable which will contains the receiving data from the communication
  auto data { new std::pair<int*,double*>[ commMap_.size() ] {} };
  int index = 0;
  #pragma omp master
  {
    gettimeofday(&time_pre_run_init[index_time],NULL);
    bool first = true;
    // receive all the information from all the MPI connection
    for ( auto& it_comm : commMap_ )
    {
      bool value [ 1 ]  = { true } ;
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 0, *std::get<0>(it_comm.second) );
      if (first){
        gettimeofday(&time_pre_run_wait[index_time],NULL);
        first = false;
      }
      data[index] = receive_spike_train(*std::get<0>(it_comm.second),*std::get<1>(it_comm.second));
      index+=1;
    }
  }
  #pragma omp barrier
  comm_map* communication_map_shared = &commMap_;
  #pragma omp parallel default(none) shared(data,communication_map_shared)
 {
   // Each thread update its own devices.
   int index_it = 0;
    for ( auto& it_comm : *communication_map_shared )
    {
      update_device(std::get<2>(it_comm.second),*std::get<1>(it_comm.second),data[index_it]);
      index_it+=1;
    }
    gettimeofday(&time_pre_run_end[index_time],NULL);
  }
  #pragma omp barrier
  #pragma omp master
  {
    // Master thread clean all the allocated memory
    clean_memory_input_data( data );
    delete[] data;
    data = nullptr;
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
    gettimeofday(&time_post_run_init[index_time],NULL);
    // Send information about the end of the running part
    for ( auto& it_comm : commMap_ )
    {
      bool value [ 1 ]  = { true } ;
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 1, *std::get<0>(it_comm.second) );
    }
    gettimeofday(&time_post_run_end[index_time],NULL);
    index_time+=1;
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
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 2, *std::get<0>(it_comm.second) );
      MPI_Comm_disconnect( std::get<0>(it_comm.second) );
      delete std::get<0>(it_comm.second);
      delete std::get<1>(it_comm.second);
      delete[] std::get<2>(it_comm.second);
      std::get<2>(it_comm.second) = nullptr;
    }
    // clear map of devices
    commMap_.clear();
    thread thread_id_master = kernel().vp_manager.get_thread_id();
    for ( auto& it_device : devices_[thread_id_master] ) {
      it_device.second.first = nullptr;
    }
    std::ofstream myfile (kernel().io_manager.get_data_path()+"/timer_input_"+std::to_string(kernel().mpi_manager.get_rank())+".txt");
    if (myfile.is_open())
    {
      myfile <<
        std::to_string((double) (time_prepare_init[0].tv_sec + time_prepare_init[0].tv_usec/1000000.0) ) <<";" <<
        std::to_string((double) (time_prepare_end[0].tv_sec + time_prepare_end[0].tv_usec/1000000.0) )<<std::endl;
      for(int count = 0; count < index_time; count ++){
          myfile <<
          std::to_string((double) (time_pre_run_init[count].tv_sec + time_pre_run_init[count].tv_usec/1000000.0) )<< ";" <<
          std::to_string((double) (time_pre_run_end[count].tv_sec + time_pre_run_end[count].tv_usec/1000000.0) )<< ";" <<
          std::to_string((double) (time_pre_run_wait[count].tv_sec + time_pre_run_wait[count].tv_usec/1000000.0) ) <<";" <<
          std::to_string((double) (time_post_run_init[count].tv_sec + time_post_run_init[count].tv_usec/1000000.0) ) <<";" <<
          std::to_string((double) (time_post_run_end[count].tv_sec + time_post_run_end[count].tv_usec/1000000.0) )<<std::endl;
      }
      myfile.close();
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
  // get the path from the kernel
  if ( not path.empty() )
  {
    basename << path << '/';
  }
  basename << kernel().io_manager.get_data_prefix();

  // add the path from the label of the device
  if ( not label.empty() )
  {
    basename << label;
  }
  else
  {
    throw MPIFilePortsUnknown( index_node );
  }
  // add the id of the device to the path
  char add_path[ 150 ];
  sprintf( add_path, "/%zu.txt", index_node );
  basename << add_path;
  std::cout << basename.rdbuf() << std::endl;

  // read the file
  std::ifstream file( basename.str() );
  if ( file.is_open() )
  {
    getline( file, *port_name );
  }
  file.close();
}

std::pair<int*,double*>
nest::InputBackendMPI::receive_spike_train( const MPI_Comm& comm, std::vector<int>& devices_id)
{
  // Send size of the list id
  int size_list[1];
  size_list [0] = devices_id.size();
  MPI_Send( &size_list, 1, MPI_INT, 0, 0, comm );
  if (size_list[0] != 0)
  {
    // Send the list of device ids
    MPI_Send( &devices_id[ 0 ], size_list[ 0 ], MPI_INT, 0, 0, comm );
    // Receive the size of data
    MPI_Status status_mpi;
    // Receive the size of the data in total and for each devices
    int* nb_size_data_per_id { new int[ size_list[ 0 ] + 1 ] {} }; // delete in the function clean_memory_input_data
    MPI_Recv( nb_size_data_per_id, size_list[ 0 ] + 1, MPI_INT, MPI_ANY_SOURCE, devices_id[ 0 ], comm, &status_mpi );
    // Receive the data
    double* data { new double[ nb_size_data_per_id[ 0 ] ] {} }; // delete in the function clean_memory_input_data
    MPI_Recv( data, nb_size_data_per_id[ 0 ], MPI_DOUBLE, status_mpi.MPI_SOURCE, devices_id[ 0 ], comm, &status_mpi );
    // return the size of the data by device and the data
    return std::make_pair(nb_size_data_per_id,data);
  }
  // if there are no data return nullptr
  return std::make_pair(nullptr ,nullptr);
}

void
nest::InputBackendMPI::update_device(int* array_index, std::vector<int>& devices_id,std::pair<int*,double*> data )
{
  if (data.first != nullptr){
    // if there are some device
    if (data.first[0] != 0)
    {
      // if there are some data
      thread thread_id = kernel().vp_manager.get_thread_id();
      int index_id_device = 0; // the index for the array of device in the data
      // get the first id of the device for the current thread
      // if the thread_id == 0, the index_id_device equals 0
      if ( thread_id != 0 )
      {
        index_id_device = std::accumulate( array_index, array_index + thread_id - 1, 0 );
      }
      // get the index of the last device by adding the number of device in the thread
      int index_id_device_end = index_id_device + array_index[ thread_id ];
      // get the index of the data receiving by summing all the size of data of the previous devices
      int index_data = std::accumulate( data.first + 1, data.first+1+index_id_device, 0 );
      if ( index_id_device != index_id_device_end )
      {
        // update all the device with spikes if there are devices
        for ( int i = index_id_device; i != index_id_device_end; i++ )
        {
          int id = devices_id[ i ];
          // create a vector with the data of the device
          std::vector< double > data_for_device(
            &data.second[ index_data ], &data.second[ index_data + data.first[ index_id_device + 1 ] ] );

          // Update the device with the data in the current thread
          devices_[ thread_id ].find( id )->second.second->update_from_backend( data_for_device );

          // update the index of the data and the index_id_device
          index_data += data.first[ index_id_device + 1 ];
          index_id_device += 1;
        }
      }
    }
  }
}

void
nest::InputBackendMPI::clean_memory_input_data(std::pair<int*,double*>* data ){
  // for all the pair of data free the memory of data and the array with teh size
  for (size_t i =0; i != commMap_.size();  i++){
    std::pair<int*,double*> pair_data = data[i];
    if (pair_data.first != nullptr)
    {
      // clean the memory allocate in the function receive_spike_train
      delete[] pair_data.first;
      pair_data.first = nullptr;
    }
    if(pair_data.second != nullptr){
      // clean the memory allocate in the function receive_spike_train
      delete[] pair_data.second;
      pair_data.second = nullptr;
    }
  }
}

