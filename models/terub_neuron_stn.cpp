/* generated by template org.nest.nestml.neuron.NeuronClass*/
/*
*  terub_neuron_stn.cpp
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
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "terub_neuron_stn.h"

/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<terub_neuron_stn> terub_neuron_stn::recordablesMap_;

namespace nest {
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <> void RecordablesMap<terub_neuron_stn>::create() {
  // use standard names whereever you can for consistency!
  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("V_m", &terub_neuron_stn::get_V_m);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("g_in", &terub_neuron_stn::get_g_in);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("g_ex", &terub_neuron_stn::get_g_ex);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("gate_h", &terub_neuron_stn::get_gate_h);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("gate_n", &terub_neuron_stn::get_gate_n);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("gate_r", &terub_neuron_stn::get_gate_r);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("Ca_con", &terub_neuron_stn::get_Ca_con);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("g_ex'", &terub_neuron_stn::get___D_g_ex);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("g_in'", &terub_neuron_stn::get___D_g_in);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the PSCurrInit_E with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the PSCurrInit_I with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the refractory_counts with the domain type long

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the r with the domain type long

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the E_L with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the g_L with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the C_m with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the E_Na with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the g_Na with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the E_K with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the g_K with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the E_Ca with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the g_Ca with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the g_T with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the g_ahp with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_syn_ex with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_syn_in with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_e with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the E_gs with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the t_ref with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_stim with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_n_0 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_n_1 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_n_tau with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_n_tau with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_h_0 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_h_1 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_h_tau with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_h_tau with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_r_0 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_r_1 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_r_tau with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_r_tau with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_a with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_a with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_h with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_h with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_m with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_m with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_n with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_n with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_r with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_r with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_s with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_s with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the theta_b with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the sigma_b with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the phi_h with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the phi_n with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the phi_r with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the epsilon with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the k_Ca with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the k1 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_ex_mod with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_in_mod with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_n with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_h with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_r with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the a_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the h_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the m_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the n_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the r_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the s_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the b_inf with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_Na with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_K with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_L with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_T with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_Ca with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_ahp with the domain type double
}
}

/* ----------------------------------------------------------------
* Default constructors defining default parameters and state
* ---------------------------------------------------------------- */

terub_neuron_stn::Parameters_::Parameters_() {}

terub_neuron_stn::State_::State_() {}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

void terub_neuron_stn::Parameters_::set(const DictionaryDatum &__d) {}

void terub_neuron_stn::State_::set(const DictionaryDatum &__d,
                                   const Parameters_ &p) {}

terub_neuron_stn::Buffers_::Buffers_(terub_neuron_stn &n)
    : logger_(n), s_(0), c_(0), e_(0) {}

terub_neuron_stn::Buffers_::Buffers_(const Buffers_ &, terub_neuron_stn &n)
    : logger_(n), s_(0), c_(0), e_(0) {}

/* ----------------------------------------------------------------
* Default and copy constructor for node
* ---------------------------------------------------------------- */
// TODO inner components
terub_neuron_stn::terub_neuron_stn() : Archiving_Node(), P_(), S_(), B_(*this) {
  recordablesMap_.create();

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.E_L = ((-60));

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.g_L = 2.25;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.C_m = 1.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.E_Na = 55;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.g_Na = 37.5;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.E_K = ((-80.0));

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.g_K = 45.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.E_Ca = 140;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.g_Ca = 0.5;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.g_T = 0.5;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.g_ahp = 9;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.tau_syn_ex = 1.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.tau_syn_in = 0.3; //0.08;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.I_e = 0.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.E_gs = ((-85.0));

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.t_ref = 3;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.I_stim = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::V_m] = get_E_L();

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::g_in] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::g_ex] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::gate_h] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::gate_n] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::gate_r] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::Ca_con] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::__D_g_ex] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y[State_::__D_g_in] = 0;
}

terub_neuron_stn::terub_neuron_stn(const terub_neuron_stn &n)
    : Archiving_Node(), P_(n.P_), S_(n.S_), B_(n.B_, *this) {}

/* ----------------------------------------------------------------
* Destructors
* ---------------------------------------------------------------- */

terub_neuron_stn::~terub_neuron_stn() {
  // GSL structs may not have been allocated, so we need to protect destruction
  if (B_.s_)
    gsl_odeiv_step_free(B_.s_);
  if (B_.c_)
    gsl_odeiv_control_free(B_.c_);
  if (B_.e_)
    gsl_odeiv_evolve_free(B_.e_);
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void terub_neuron_stn::init_state_(const Node &proto) { // TODO inner components

  const terub_neuron_stn &pr = downcast<terub_neuron_stn>(proto);
  S_ = pr.S_;
}

/* generated by template org.nest.nestml.function.GSLDifferentiationFunction*/
extern "C" inline int terub_neuron_stn_dynamics(double, const double y[],
                                                double f[], void *pnode) {
  typedef terub_neuron_stn::State_ State_;
  // get access to node so we can almost work as in a member function
  assert(pnode);
  const terub_neuron_stn &node = *(reinterpret_cast<terub_neuron_stn *>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  double tau_n_0 = 1.0;
  double tau_n_1 = 100.0;
  double theta_n_tau = ((-80.0));
  double sigma_n_tau = ((-26.0));
  double tau_h_0 = 1.0;
  double tau_h_1 = 500.0;
  double theta_h_tau = ((-57.0));
  double sigma_h_tau = ((-3.0));
  double tau_r_0 = 40.0; // 7.1 or 40;
  double tau_r_1 = 17.5;
  double theta_r_tau = 68.0;
  double sigma_r_tau = ((-2.2));
  double theta_a = ((-63.0));
  double sigma_a = 7.8;
  double theta_h = ((-39.0));
  double sigma_h = ((-3.1));
  double theta_m = ((-30.0));
  double sigma_m = 15.0;
  double theta_n = ((-32.0));
  double sigma_n = 8.0;
  double theta_r = ((-67.0));
  double sigma_r = ((-2.0));
  double theta_s = ((-39.0));
  double sigma_s = 8.0;
  double theta_b = 0.4;  // 0.25 or  0.4;
  double sigma_b = -0.1; // 0.07 or -0.1;
  double phi_h = 0.75;
  double phi_n = 0.75;
  double phi_r = 0.2; // 0.5 or 0.2;
  double epsilon = 0.0000375; // 0.00005 or 0.0000375;
  double k_Ca = 22.5;
  double k1 = 15.0;
  double I_ex_mod = ((-y[State_::g_ex])) * y[State_::V_m];
  double I_in_mod = y[State_::g_in] * (y[State_::V_m] - node.get_E_gs());
  double tau_n =
      tau_n_0 +
      tau_n_1 /
          (1. + std::exp(((-(y[State_::V_m] - theta_n_tau))) / sigma_n_tau));
  double tau_h =
      tau_h_0 +
      tau_h_1 /
          (1. + std::exp(((-(y[State_::V_m] - theta_h_tau))) / sigma_h_tau));
  double tau_r =
      tau_r_0 +
      tau_r_1 /
          (1. + std::exp(((-(y[State_::V_m] - theta_r_tau))) / sigma_r_tau));
  double a_inf =
      1. / (1. + std::exp(((-(y[State_::V_m] - theta_a))) / sigma_a));
  double h_inf =
      1. / (1. + std::exp(((-(y[State_::V_m] - theta_h))) / sigma_h));
  double m_inf =
      1. / (1. + std::exp(((-(y[State_::V_m] - theta_m))) / sigma_m));
  double n_inf =
      1. / (1. + std::exp(((-(y[State_::V_m] - theta_n))) / sigma_n));
  double r_inf =
      1. / (1. + std::exp(((-(y[State_::V_m] - theta_r))) / sigma_r));
  double s_inf =
      1. / (1. + std::exp(((-(y[State_::V_m] - theta_s))) / sigma_s));
  double b_inf = 1. / (1. + std::exp((y[State_::gate_r] - theta_b) / sigma_b)) -
                 1. / (1. + std::exp(((-theta_b)) / sigma_b));
  double I_Na = node.get_g_Na() * m_inf * m_inf * m_inf * y[State_::gate_h] *
                (y[State_::V_m] - node.get_E_Na());
  double I_K = node.get_g_K() * y[State_::gate_n] * y[State_::gate_n] *
               y[State_::gate_n] * y[State_::gate_n] *
               (y[State_::V_m] - node.get_E_K());
  double I_L = node.get_g_L() * (y[State_::V_m] - node.get_E_L());
  double I_T = node.get_g_T() * a_inf * a_inf * a_inf * b_inf * b_inf *
               (y[State_::V_m] - node.get_E_Ca());
  double I_Ca =
      node.get_g_Ca() * s_inf * s_inf * (y[State_::V_m] - node.get_E_Ca());
  double I_ahp = node.get_g_ahp() *
                 (y[State_::Ca_con] / (y[State_::Ca_con] + k1)) *
                 (y[State_::V_m] - node.get_E_K());

  f[State_::V_m] = (((-(I_Na + I_K + I_L + I_T + I_Ca + I_ahp))) +
                    node.get_I_stim() + node.get_I_e() + I_ex_mod + I_in_mod) /
                   node.get_C_m();
  f[State_::g_in] =
      y[State_::__D_g_in] - (y[State_::g_in] / node.get_tau_syn_in());
  f[State_::g_ex] =
      y[State_::__D_g_ex] - (y[State_::g_ex] / node.get_tau_syn_ex());
  f[State_::gate_h] = phi_h * ((h_inf - y[State_::gate_h]) / tau_h);
  f[State_::gate_n] = phi_n * ((n_inf - y[State_::gate_n]) / tau_n);
  f[State_::gate_r] = phi_r * ((r_inf - y[State_::gate_r]) / tau_r);
  f[State_::Ca_con] = epsilon * (((-I_Ca)) - I_T - k_Ca * y[State_::Ca_con]);
  f[State_::__D_g_ex] = ((-y[State_::__D_g_ex])) / node.get_tau_syn_ex();
  f[State_::__D_g_in] = ((-y[State_::__D_g_in])) / node.get_tau_syn_in();

  return GSL_SUCCESS;
}

void terub_neuron_stn::init_buffers_() {
  get_spikeInh().clear(); // includes resize
  get_spikeExc().clear(); // includes resize
  get_currents().clear(); // includes resize
  B_.logger_.reset();     // includes resize
  Archiving_Node::clear_history();
  if (B_.s_ == 0)
    B_.s_ = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 9);
  else
    gsl_odeiv_step_reset(B_.s_);

  if (B_.c_ == 0) {
    B_.c_ = gsl_odeiv_control_y_new(1e-6, 0.0);
  } else {
    gsl_odeiv_control_init(B_.c_, 1e-6, 0.0, 1.0, 0.0);
  }

  if (B_.e_ == 0) {
    B_.e_ = gsl_odeiv_evolve_alloc(9);
  } else {
    gsl_odeiv_evolve_reset(B_.e_);
  }

  B_.sys_.function = terub_neuron_stn_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = 9;
  B_.sys_.params = reinterpret_cast<void *>(this);
}

void terub_neuron_stn::calibrate() {
  B_.logger_.init();

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.PSCurrInit_E = 1.0 * numerics::e / P_.tau_syn_ex;

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.PSCurrInit_I = 1.0 * numerics::e / P_.tau_syn_in;

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.refractory_counts =
      nest::Time(nest::Time::ms((double)P_.t_ref)).get_steps();

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.r = 0;
}

/* ----------------------------------------------------------------
* Update and spike handling functions
* ---------------------------------------------------------------- */

/*

 */
void terub_neuron_stn::update(nest::Time const &origin, const long from,
                              const long to) {

  double step_ = nest::Time::get_resolution().get_ms();
  double IntegrationStep_ = nest::Time::get_resolution().get_ms();
  double t = 0;

  for (long lag = from; lag < to; ++lag) {
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.spikeInh_last_value_ = get_spikeInh().get_value(lag);
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.spikeExc_last_value_ = get_spikeExc().get_value(lag);
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.currents_last_value_ = get_currents().get_value(lag);

    /* generated by template org.nest.spl.Block*/
    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.Declaration*/

    double U_old = S_.y[State_::V_m];

    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.FunctionCall*/
    /* generated by template org.nest.spl.GSLIntegrator*/
    t = 0;

    while (t < step_) {
      const int status =
          gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_,
                                 &B_.sys_,          // system of ODE
                                 &t,                // from t
                                 step_,             // to t <= step
                                 &IntegrationStep_, // integration step size
                                 S_.y);             // neuronal state

      if (status != GSL_SUCCESS) {
        throw nest::GSLSolverFailure(get_name(), status);
      }
    }

    /* generated by template org.nest.spl.Statement*/
    // # sending spikes: crossing 0 mV, pseudo-refractoriness and local
    // maximum...
    /* generated by template org.nest.spl.CompoundStatement*/
    /* generated by template org.nest.spl.IfStatement*/

    if (V_.r > 0) {
      /* generated by template org.nest.spl.Block*/
      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.Assignment*/
      V_.r -= 1;

    } else if (((S_.y[State_::V_m] > 0)) && ((U_old > S_.y[State_::V_m]))) {

      /* generated by template org.nest.spl.Block*/
      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.Assignment*/
      V_.r = V_.refractory_counts;

      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.FunctionCall*/
      set_spiketime(nest::Time::step(origin.get_steps() + lag + 1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
      ;

    } /* if end */

    /* generated by template org.nest.spl.Statement*/
    // # set new input current
    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.Assignment*/
    P_.I_stim = B_.currents_last_value_;

    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.Assignment*/
    S_.y[State_::__D_g_ex] += B_.spikeExc_last_value_ * V_.PSCurrInit_E;

    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.Assignment*/
    S_.y[State_::__D_g_in] += B_.spikeInh_last_value_ * V_.PSCurrInit_I;

    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void terub_neuron_stn::handle(nest::DataLoggingRequest &e) {
  B_.logger_.handle(e);
}

void terub_neuron_stn::handle(nest::SpikeEvent &e) {
  assert(e.get_delay() > 0);

  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  /* generated by template org.nest.nestml.buffer.SpikeBufferFill*/

  if (weight < 0.0) // inhibitory
  {
    get_spikeInh().add_value(
        e.get_rel_delivery_steps(
            nest::kernel().simulation_manager.get_slice_origin()),

        weight * multiplicity);
  }

  /* generated by template org.nest.nestml.buffer.SpikeBufferFill*/

  if (weight >= 0.0) // excitatory
  {
    get_spikeExc().add_value(
        e.get_rel_delivery_steps(
            nest::kernel().simulation_manager.get_slice_origin()),
        weight * multiplicity);
  }
}

void terub_neuron_stn::handle(nest::CurrentEvent &e) {
  assert(e.get_delay() > 0);

  const double current = e.get_current();
  const double weight = e.get_weight();

  // add weighted current; HEP 2002-10-04
  get_currents().add_value(
      e.get_rel_delivery_steps(
          nest::kernel().simulation_manager.get_slice_origin()),
      weight * current);
}
