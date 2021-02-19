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
 * \file growth_curve.cpp
 * Implementation of growth_curve
 * \author Mikael Naveau
 * \date July 2013
 */

#include "deletion_curve.h"

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "nest_names.h"
#include "nest_time.h"

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
{
}

void
nest::DeletionCurveLinear::get( DictionaryDatum& d ) const
{
  def< std::string >( d, names::deletion_curve, name_.toString() );
  def< double >( d, names::deletion_probablity, deletion_probability_ );
}

void
nest::DeletionCurveLinear::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::deletion_probablity, deletion_probability_ );
}

double
nest::DeletionCurveLinear::update( double t,
  double t_minus,
  double Ca_minus,
  double z_connected,
  double tau_Ca,
  double deletion_rate ) const
{
  //std::cout << "("<< std::floor(se_it->second.get_z_connected()*std::min(const_z_deletion + 0.02 * Ca_minus_ / tau_Ca_*10000. , max_delete_z)) << "," << const_z_deletion + 0.02 * Ca_minus_ / tau_Ca_*10000. << ") "; 
  //std::cout<<Ca_minus_ * 1000.<<" ";
  librandom::RngPtr rng = kernel().rng_manager.get_rng( get_thread() );
  librandom::BinomialRandomDev bino_dev ;
  double pbino =  const_z_deletion +  max_delete_z /(1. + std::exp( - (Ca_minus_ *1000. - 100.) /10. ));
  int nbino = se_it->second.get_z_connected();
  bino_dev.set_p_n( pbino, nbino);

  return bino_dev.ldev( rng );
}

