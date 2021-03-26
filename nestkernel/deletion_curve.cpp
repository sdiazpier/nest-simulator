/*
 *  deletion_curve.cpp
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

/**
 * \file deletion_curve.cpp
 * Implementation of deletion_curve
 */

#include "deletion_curve.h"

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "nest_names.h"
#include "nest_time.h"
#include "kernel_manager.h"

// Includes from librandom
#include "binomial_randomdev.h"

// Includes from sli:
#include "dictutils.h"

/* ----------------------------------------------------------------
 * DeletionCurveLinear
 * ---------------------------------------------------------------- */

nest::DeletionCurveLinear::DeletionCurveLinear()
  : DeletionCurve( names::linear )
  , eps_( 0.7 )
  , max_delete_z_ ( 0.05 )
  , const_z_deletion_ ( 0.0001 )
{
}

void
nest::DeletionCurveLinear::get( DictionaryDatum& d ) const
{
  def< std::string >( d, names::deletion_curve, name_.toString() );
  def< double >( d, names::deletion_probability, deletion_probability_ );
  def< double >( d, names::max_delete_z, max_delete_z_ );
  def< double >( d, names::const_z_deletion, const_z_deletion_ );
}

void
nest::DeletionCurveLinear::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::deletion_probability, deletion_probability_ );
  updateValue< double >( d, names::max_delete_z, max_delete_z_ );
  updateValue< double >( d, names::const_z_deletion, const_z_deletion_ );
}

double
nest::DeletionCurveLinear::update( int z_connected,
  double Ca_minus,
  thread thrd ) const
{
  //std::cout << "("<< std::floor(se_it->second.get_z_connected()*std::min(const_z_deletion + 0.02 * Ca_minus_ / tau_Ca_*10000. , max_delete_z)) << "," << const_z_deletion + 0.02 * Ca_minus_ / tau_Ca_*10000. << ") "; 
  //std::cout<<Ca_minus_ * 1000.<<" ";
  librandom::RngPtr rng = kernel().rng_manager.get_rng( thrd );
  librandom::BinomialRandomDev bino_dev ;
  double pbino =  const_z_deletion_ +  max_delete_z_ /(1. + std::exp( - (Ca_minus *1000. - 30.) /1. ));
  bino_dev.set_p_n( pbino, z_connected);

  return bino_dev.ldev( rng );
}

