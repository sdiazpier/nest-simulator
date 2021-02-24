/*
 *  deletion_curve.h
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

#ifndef DELETION_CURVE_H
#define DELETION_CURVE_H

/**
 * \file deletion_curve.h
 *
 */

// Includes from nestkernel:
#include "nest_types.h"
#include "exceptions.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{

/**
 * \class DeletionCurve
 * Defines the way the number of synaptic elements changes through time
 * according to the calcium concentration of the neuron.
 */
class DeletionCurve
{
public:
  virtual ~DeletionCurve()
  {
  }
  virtual void get( DictionaryDatum& d ) const = 0;
  virtual void set( const DictionaryDatum& d ) = 0;
  virtual double
  update( int z_connected, double Ca_minus, thread thrd ) const = 0;
  virtual bool
  is( Name n )
  {
    return n == name_;
  }
  Name
  get_name()
  {
    return name_;
  }

protected:
  DeletionCurve( const Name name )
    : name_( name )
  {
  }
  const Name name_;

  double deletion_probability_;
  
  // Max number of deleted synapses per update step dependent on the firing rate (precentage)
  double max_delete_z_;

  // constant deleted synaptic elements per update step
  int const_z_deletion_;
};

/** @BeginDocumentation

  Name: deletion_curve_linear - Linear version of a growth curve

  SeeAlso: SynapticElement, SPManager, SPBuilder, DeletionCurveLinear
*/
/**
 * \class DeletionCurveLinear
 */
class DeletionCurveLinear : public DeletionCurve
{
public:
  DeletionCurveLinear();
  void get( DictionaryDatum& d ) const;
  void set( const DictionaryDatum& d );
  double update( int z_connected, double Ca_minus, thread thrd ) const;

private:
  double eps_;
  double deletion_probability_;
  
  // Max number of deleted synapses per update step dependent on the firing rate (precentage)
  double max_delete_z_;

  // constant deleted synaptic elements per update step
  int const_z_deletion_;
};


} // of namespace

#endif
