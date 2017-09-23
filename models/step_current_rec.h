/*
 *  step_current_rec.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */


/*BeginDocumentation
  Name: step_current_rec - provides a piecewise constant DC input current
  Description:
  The dc_generator provides a piecewise constant DC input to the
  connected node(s).  The amplitude of the current is changed at the
  specified times. The unit of the current is pA.

 Parameters:
     The following parameters can be set in the status dictionary:
     amplitude_times   list of doubles - Times at which current changes in ms
     amplitude_values  list of doubles - Amplitudes of step current current in
 pA

  Examples:
    The current can be altered in the following way:
    /step_current_generator Create /sc Set
    sc << /amplitude_times [0.2 0.5] /amplitude_values [2.0 4.0] >> SetStatus

    The amplitude of the DC will be 0.0 pA in the Time interval [0, 0.2),
    2.0 pA in the interval [0.2, 0.5) and 4.0 from then on.

  Sends: CurrentEvent

  Author: Jochen Martin Eppler, Jens Kremkow

  SeeAlso: ac_generator, dc_generator, step_current_generator, Device,
 StimulatingDevice
*/

#ifndef STEP_CURRENT_REC_H
#define STEP_CURRENT_REC_H

#include <vector>
#include "nest.h"
#include "event.h"
#include "node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "stimulating_device.h"
#include "universal_data_logger.h" //rec from device

namespace nest
{
  class Network; //rec from device

  class step_current_rec : public Node
  {
    
  public:        
    
    step_current_rec();
    step_current_rec(const step_current_rec&);

    using Node::handles_test_event;  //rec from device
    using Node::handle;              //rec from device

    bool has_proxies() const {return false;} 
    bool local_receiver() const { return true;  } //rec from device

    void handle(DataLoggingRequest &);                //rec from device
    port handles_test_event(DataLoggingRequest&, rport);   //rec from device
    
    port send_test_event(Node&, rport, synindex, bool);
 
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  private:
    
    void init_node_(const Node&);
    void init_state_(const Node&);
    void init_buffers_();
    void calibrate();
    
    void update(Time const &, const long, const long);
    
    struct Buffers_;
    
    /**
     * Store independent parameters of the model.
     */
    struct Parameters_ {
      std::vector<double> amp_times_;
      std::vector<double> amp_values_;
      
      Parameters_();  //!< Sets default parameter values
      Parameters_(const Parameters_&, Buffers_&);

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&, Buffers_&);  //!< Set values from dicitonary
    };

    // ------------------------------------------------------------
    
    friend class RecordablesMap<step_current_rec>;  // rec from device
    friend class UniversalDataLogger<step_current_rec>;// rec from device

    struct Buffers_ {
      Buffers_(step_current_rec&);  // rec from device
      Buffers_(const Buffers_&, step_current_rec&);  // rec from device
      UniversalDataLogger<step_current_rec> logger_;  // rec from device
      size_t   idx_;  //!< index of current amplitude
      double amp_;  //!< current amplitude
    };
   

    double get_amp_() const { return B_.amp_; } //rec from device 
    // ------------------------------------------------------------

    StimulatingDevice<CurrentEvent> device_;
    static RecordablesMap<step_current_rec> recordablesMap_; // rec from device
    Parameters_ P_;
    Buffers_    B_;
  };
  

  inline
    port step_current_rec::handles_test_event(DataLoggingRequest& dlr, 
				      rport receptor_type)
    {
      if (receptor_type != 0)
	throw UnknownReceptorType(receptor_type, get_name());
      return B_.logger_.connect_logging_device(dlr, recordablesMap_);
    }
  

  inline
  port step_current_rec::send_test_event(Node& target, rport receptor_type, 
	        			  synindex syn_id, bool dummy_target)
{
  device_.enforce_single_syn_type(syn_id);

  CurrentEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

  inline
  void step_current_rec::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    device_.get_status(d);
    (*d)[names::recordables] = recordablesMap_.get_list(); // rec from device
  }

  inline
  void step_current_rec::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d, B_);               // throws if BadProperty

    // We now know that ptmp is consistent. We do not write it back
    // to P_ before we are also sure that the properties to be set
    // in the parent class are internally consistent.
    device_.set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
  }
  
} // namespace

#endif /* #ifndef STEP_CURRENT_REC_H */
